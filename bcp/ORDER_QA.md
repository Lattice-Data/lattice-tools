# Order QA (`order_qa`)

QA one sequencing order in one command: verify the upload against S3, run every
check that applies to the assay, write a report, and print a machine-readable
status line. No notebook, no kernel state, no transcribing results by hand.

```bash
python -m order_qa czi-psomagen/PROJECT/AN00000001 --assay 10x
```

## Why this exists

`qa.ipynb` covers seven assays across raw and processed data, and it decides what
applies inside each cell (`if validate_raw and raw_assay == 'scale'`, twenty-odd
times). Every branch prints why it skipped, which is correct — and invisible. A
skipped cell scrolls past exactly like a passing one, and there is no point at
which the notebook can tell you that eleven of its checks never ran.

Here the checks are declared as data. The run is enumerable before it starts, and
**every declared check lands in the report with a status**, so "skipped" is a
reported outcome rather than an absence. The distinction that carries the tool:

| status | meaning | affects verdict |
| --- | --- | --- |
| `pass` | ran, found nothing | — |
| `findings` | ran, found defects | yes, exit 5 |
| `skipped` | does not apply to this order (wrong assay, no `processed/`) | no |
| `not_run` | **applies and could not be performed** | yes, exit 6 |
| `error` | raised | yes, exit 6 |

`skipped` implies nothing. `not_run` denies the order a verdict — a check that
could not run is not a check that passed.

Nothing here reimplements what a check *means*. Every one delegates to the same
`qa_checks` function the notebook calls, and gathering goes through
`qa_gather.gather_qa_data` unchanged.

## Invocation

Two equivalent forms. The module form matches every other tool in `bcp/`:

```bash
python -m order_qa czi-psomagen/PROJECT/AN00000001 --assay 10x
```

The dispatcher gives you the `bcp <command>` form. It `cd`s to `bcp/` first, so it
works from anywhere and can be symlinked onto `PATH`:

```bash
ln -s "$PWD/bcp/bcp" ~/.local/bin/bcp
bcp order-qa czi-psomagen/PROJECT/AN00000001 --assay 10x
```

Both must resolve `bcp/` on `sys.path` — the `qa_*` modules import each other flat
(`from qa_constants import ...`), which is why the module form is run from inside
`bcp/` and why the dispatcher changes directory.

## The order spec

Accepts, in order of preference:

```
czi-psomagen/PROJECT/AN00000001                     bucket/project/order
s3://czi-psomagen/PROJECT/AN00000001                same thing
czi-psomagen/PROJECT/AN00000001/raw/x.fastq.gz      a pasted key resolves to its order
AN00000001                                          bare ID, if unambiguous
```

Identity is `(bucket, project, order)`, never the order ID alone: two vendors can
issue the same-looking identifier, and the same ExperimentID appears under both a
lab bucket and the vendor bucket it was trimmed from.

Order-ID shapes (`AN########`, `NVUS...-##`) are **reported, not enforced** — a
vendor inventing a third format gets a note on a run that happened, not an error
on one that did not. What *is* enforced is the closed list of things that are not
orders: `jsonlds`, `metadata-tsv`, and the `test*` vendor smoke tests, which get
refreshed periodically and so reappear after anyone deletes them.

A bare ID is resolved by name against `WATCHED_PROJECTS` in `order_qa/spec.py`,
never by listing S3 to find it. With more than one watched project it is
ambiguous, and ambiguity is a question for the caller rather than something to
settle by taking whichever listing answered first.

## `--assay` is required

It cannot be derived. The assay determines which files an order is *expected* to
contain (`qa_constants.raw_expected` is keyed by it), and it varies between orders
from the same vendor — so there is nothing in `bucket/project/order` to read it
off, and no per-project default that would be right. A wrong default would not
fail loudly; it would check the order against the wrong expectations and report a
clean pass.

Values: `10x`, `10x_cram`, `10x_viral_ORF`, `scale`, `sci_jumbo`, `sci_plex`,
`seahub_sci`.

`--dry-run` and `--probe` need no assay, because verification is
assay-independent.

## What runs, per assay

Raw/processed scope is auto-detected from the listing (`raw/` and `processed/`
present) and can be forced with `--raw` / `--no-raw` / `--processed` /
`--no-processed`. `seahub_sci` is raw-only QA, as in the notebook.

| check | applies to | delegates to |
| --- | --- | --- |
| `raw_expected_files` | raw, all assays | `check_expected_raw_files` |
| `raw_extra_files` | raw, all assays | `check_extra_raw_files` |
| `raw_fastq_counts` | raw, all assays | `validate_fastq_counts` |
| `raw_read_metadata` | raw, all assays | `validate_read_metadata` |
| `raw_pct_q30` | raw, all assays | `validate_pct_q30` |
| `scale_cb_tag` | raw, `scale` | `validate_scale_cb_tag` |
| `processed_cellranger` | processed, all but `scale` | `validate_processed_group` |
| `scale_workflow_info` | processed, `scale` | `validate_scale_workflow_info`, `validate_scale_samples_csv` |
| `scale_processed_files` | processed, `scale` | `validate_scale_processed_files` |
| `seahub_sop_naming` | raw, `seahub_sci` | gathered `sop_violations` |
| `seahub_trimming_completeness` | raw, `seahub_sci` | **not implemented — see below** |
| `seahub_well_status` | raw, `seahub_sci` | **not implemented — see below** |

### Known gap: SeaHub cross-bucket reconciliation

The two SeaHub reconciliation checks are **declared but not built**. A
`seahub_sci` order therefore reports them as `not_run` and exits 6, rather than
passing on the other checks while silently omitting the cross-bucket
reconciliation that is the whole point of QA-ing a trimmed upload.

Declaring the gap is the difference between a known gap and an invisible one. Run
`qa.ipynb` for `seahub_sci` orders until they are wired up.

## Upload verification

Asked before QA, because QA against a live upload burns a full gather and reports
absences as defects.

**Quiescence** — the newest object must be older than `--quiescence-minutes`
(default 15). This is the entire distinction between "the vendor said the order is
done" and "the order is actually done". `--force` proceeds anyway and the report
records that the upload may have been in flight.

**Completeness** — `--manifest` compares the listing against expected
`s3://` URIs (CSV or TSV; URIs are found anywhere on a line, since column order
varies between hand-assembled manifests). Without one, the report says the
manifest check *did not run* — it never reads as a pass. There is no per-order
manifest of expected contents in the vendor buckets; what an order should contain
comes from the assay's naming convention via `raw_expected`, which is what
`raw_expected_files` checks.

A size is attributed to a key only when the line leaves no choice — one object
named, one bare integer. Taking the first integer read the *lane* out of
`sampleA,1,s3://...,1048576` and reported every file as the wrong size.

**Versioning** — the buckets are versioned, so re-uploading a key creates a new
version instead of overwriting. A silent re-upload leaves the key count and the
byte total untouched; what changes is the key's current `VersionId`. That is
recorded per run and compared on the next one, which is how a resequence shows up
at all.

## Integrity: what can actually be claimed

Only `ListObjectsV2`, `ListBucket` and `HeadObject` are confirmed on the vendor
buckets. `GetObjectAttributes` has never been exercised, and `ListObjectVersions`
needs `s3:ListBucketVersions` — a different permission from the `s3:ListBucket`
we have. The objects are also SSE-KMS encrypted with a key in a third account, so
a call can fail on the key policy rather than the bucket policy.

So capabilities are **probed against a real object in the order**, once per run,
and reported. Nothing falls back silently.

```bash
python -m order_qa czi-psomagen/PROJECT/AN00000001 --probe
```

Run this once per bucket to find out what the reports can claim. It prints the
probe result as JSON and exits without QA-ing anything.

Consequences that are stated in every report rather than assumed:

- **ETag is MD5 only for single-part uploads.** FASTQs and CRAMs are multipart,
  where the ETag is a digest of part digests. The report counts multipart objects
  so no report ever implies a checksum comparison that was not made.
- **Version identities** come from `ListObjectVersions` when permitted (one
  paginated call set, and the only way to see delete markers) and otherwise from
  `HeadObject` per key, capped by `--head-limit` (default 500). A capped sweep is
  reported as partial. A bounded check reported as complete is how a 4000-file
  order gets released on the strength of its first 500 files.
- **Delete markers** are invisible to the `HeadObject` method: a delete-marked key
  404s, indistinguishable from one that was never uploaded.

No writes to S3, ever. Reports go to EFS.

## Reports and re-runs

```
~/qa-reports/{bucket}/{project}/{order}/{timestamp}/
    report.md       the verdict, what ran, what did not, and why
    summary.json    the same thing machine-readable
    errors.txt      every finding, flat — the notebook's {order}_errors.txt
    *.csv           per-check tables (raw_missing, process_alerts, ...)
~/qa-reports/{bucket}/{project}/{order}/ledger.json
```

Run directories are timestamped and never rewritten, because a re-run after a
resequence must not destroy the report that justified the first verdict.
Same-second runs get `-2`, `-3` suffixes rather than merging.

`ledger.json` is the baseline the *next* run compares against. It is replaced only
when a run actually observed the whole settled order — a failed listing or an
in-flight upload leaves the previous baseline alone and says so, because
overwriting it would make every unseen object look "new" next time and would lose
the version identities that catch a silent re-upload.

Change `--output-dir` to write elsewhere.

## Status line and exit codes

stdout carries exactly one line of JSON. Everything human-readable goes to stderr
via logging, so this works:

```bash
python -m order_qa czi-psomagen/PROJECT/AN00000001 --assay 10x | jq .verdict
```

```json
{"order": "...", "assay": "10x", "verdict": "qa_findings", "exit_code": 5,
 "checks_passed": 5, "checks_with_findings": 1, "checks_skipped": 6,
 "checks_unanswered": 0, "findings": 2, "files": 412, "bytes": 98765432101,
 "quiet": true, "changed_since_last_run": false, "report_dir": "...",
 "timestamp": "..."}
```

| code | name | meaning |
| --- | --- | --- |
| 0 | `ok` | verified, quiescent, every applicable check passed |
| 1 | `internal_error` | the run could not be completed or recorded |
| 2 | `usage` | bad arguments |
| 3 | `still_uploading` | objects newer than the quiescence window — retry later |
| 4 | `verification_failed` | nothing there, or a manifest mismatch |
| 5 | `qa_findings` | verified, but QA found defects |
| 6 | `degraded` | an applicable check could not run, so there is no verdict |
| 7 | `resolution_failed` | the spec named nothing QA-able |

Precedence follows what to do first: 3 before 4 before 6 before 5. Findings are
not worth triaging until the gaps are closed.

Two cases that deliberately return 6 rather than 0:

- **Nothing applied.** If every check skipped — neither `raw/` nor `processed/`
  detected, or the wrong `--assay` — the order was not validated, and summing
  those skips into a pass would be the exact false clean pass this tool exists to
  prevent.
- **Incomplete coverage.** `qa_gather` warnings that mean a check did not see
  everything it should (`LISTING TRUNCATED`, `METADATA UNREADABLE`,
  `TRIM FAIL UNREADABLE`, `EXPERIMENT UNKNOWN`) degrade the verdict.

## Credentials

`boto3.client("s3")` with no arguments, exactly like every other tool here. In the
JupyterHub pod that picks up IRSA (`AWS_WEB_IDENTITY_TOKEN_FILE` + `AWS_ROLE_ARN`)
through the default chain. Do not set `AWS_PROFILE` and do not add key handling —
there is a stale `~/.aws/credentials` on those pods that anything explicit would
pick up instead.

## Tests

```bash
cd bcp && python3 -m pytest tests/test_order_qa_*.py -q
```

No network: a stub S3 client covers the listing, quiescence, capability-denial and
version-fallback paths. The load-bearing test is
`TestRegistryInvariant::test_every_declared_check_is_reported` — one result per
declaration, for every assay and scope combination. If that breaks, this CLI has
the same defect the notebook has and is worth less than the notebook.

Test identifiers must be sanitized (`AN0000000X`, `NVUS0000000000-`, no real lab
names); `tests/test_sanitized_identifiers.py` enforces it.
