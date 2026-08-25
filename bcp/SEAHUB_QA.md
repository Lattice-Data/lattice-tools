## SeaHub lab raw QA

Part of BCP tooling. Unlike the other tools documented here there is no CLI: the checks live in
`qa_seahub_sop.py`, `qa_seahub_source.py`, `qa_seahub_recon.py` and `qa_seahub_rename.py`, and are
driven from `bcp/qa.ipynb`.

Collaborator lab trimmed uploads (`czi-trapnell` / `czi-hamazaki`, `*-seahub-bcp`) are QA'd in
`bcp/qa.ipynb` with `raw_assay='seahub_sci'`; processed validation is off for this mode. The SOP is:

```
s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/{sublibrary}/{wafer}/
    {wafer}-{sublibrary}[_{well}]_{sublibrary type}-{UG}-{barcode}.trim.*
```

with `{sublibrary type}` in `GEX`, `CRI`, `hash_oligo`, `GEX_hash_oligo`. The ExperimentID is not a
filename field; it appears in a name only when it is part of that project's sublibrary name (e.g.
`REF3_P05_2`).

---

## Data source modes

Both `data_source` modes are supported, and they are the same question asked twice, so they must
return the same objects: the s3 walk lists everything under `{ExperimentID}/raw/` recursively,
matching manifest mode's `"/raw/" in key`, and uses a delimiter pass only to enumerate wafer
*folders*. Reading just the folder tree and keeping the wafer-level objects, as it once did,
dropped every object above wafer depth — so `bad_path_depth` could fire only for keys that were too
deep, never too shallow, and a `raw/` holding only loose objects was reported as missing. In
**s3** mode the ExperimentID is the last segment of the
listing prefix. Folder enumeration — the `processed/` check, the sublibrary list and the wafer list —
uses single `list_objects` calls, which return at most 1000 entries and say so only in `IsTruncated`;
that flag is now checked and reported, so a cap nobody would think to test for cannot pass silently.
Unreachable with anything observed (REF3 has 7 sublibraries, GENE7 has 9), and since objects are
paginated a truncated listing costs `discovered_wafers` entries or the `processed/` notice, never a
raw file. In **manifest** mode it comes from the `order` argument if one is given, and is
otherwise read off the manifest keys, which already contain it as a folder; a manifest mixing two
ExperimentIDs is rejected at `resolve_qa_run_context` rather than per object, because QA writes one
`{order}_*` output set per run. `ctx.order` must never be empty for a SeaHub run: it is the expected
value of the cross-experiment check on the *upload's* own paths, so an empty one turns that into one
error per object. It no longer reaches the vendor side at all — see
[Why the vendor ExperimentID is not the key](#why-the-vendor-experimentid-is-not-the-key).

---

## What is checked

Five checks are specific to this mode. All SeaHub cells form one contiguous block at the *end* of
`qa.ipynb`, tagged `seahub` and gated on `seahub_active`, so they cannot renumber the shared pipeline;
`tests/test_qa_notebook_layout.py` enforces that. Read `{order}_seahub_well_status.csv` first: the SOP
table says what is wrong, the per-well status says how many wells are affected and which need data
rather than a rename.

- **Per-well fetches share one thread pool.** A SeaHub upload carries one trim failure CSV *and* one
  `.cram-metadata.json` per well, so 336 of each on REF3 and 864 on GENE7. Both go through a
  16-worker pool; fetching the failure CSVs inside the walk instead made them one sequential
  round-trip per well (measured on REF3 at 20 ms per object: 10.6 s against 1.2 s). Only the download
  and the parse are concurrent — the parsed blocks are applied single-threaded in listing order,
  which is what keeps this lock-free and makes the output identical to the serial version rather
  than merely equivalent, since the three structures they feed are appended to. Neither worker
  raises: with one object per well, a single zero-byte or half-written file would otherwise abort the
  gather and lose every other file with it, so an unreadable one degrades — the well simply has no
  counts, which reconciles as `metadata_unavailable` — and the failures are collapsed into one
  warning so the thinner result is never silent.
- **Two suffix families.** `qa_constants.SEAHUB_TRIM_SUFFIXES` is the SOP set; `SEAHUB_BARE_SUFFIXES`
  is the same six artifacts with the `.trim` infix dropped, which real uploads have used. Completeness
  asks only whether each of the five artifact *kinds* arrived, under either spelling — naming is a
  separate axis, reported once per stem as `missing_trim_infix`. Judging by family instead let one
  optional sidecar carrying the other family's name decide the whole requirement set, and the two
  completeness paths reacted to that differently: a well with five correct `.trim.*` files plus a
  stray bare sidecar was complete to `roll_up_wells` and missing five files to
  `check_expected_raw_files`. A genuinely absent artifact is still reported, under its SOP name
  whichever spelling the upload used.
  A CRAM whose `.cram-metadata.json` sidecar is absent is reported too, since read-count QA cannot
  run without it — as one warning naming the count and two examples, broken down by delivered
  spelling, not one per well. CZI generates these sidecars for the upload as a whole, so absence is
  an upload-wide fact; the per-well form printed 336 lines on REF3 and 864 on GENE7.
- **SOP validation** (`qa_seahub_sop.py`) reports each broken rule once per distinct fact, at four
  scopes: `object` (per object), `stem` (per well), `suffix` (per distinct unrecognised
  extension) and `upload` (per distinct bucket/project fact, so a wrong bucket is one row rather
  than one per well — on a 288-well upload that is what keeps every other finding visible).
  `lab_project_mismatch` matches the project against the lab as an exact prefix, not on the first
  hyphen-separated token: a lab whose name contains a hyphen would otherwise fail against its own
  project, and being `upload` scope that is one row for the whole listing plus a banner on a
  correct upload.
  A lab that spells its artifacts something else entirely misspells *every* object, so
  `unexpected_suffix` collapses to one row per extension naming the count, an example and the SOP
  artifacts expected; if nothing at all is recognised, `no_recognized_artifacts` fires once and the
  notebook prints a FAIL banner, since no well can be identified and every table below would
  otherwise be empty and read as clean. A bare suffix is only accepted when the stem before it
  parses — `.csv` is generic enough that `<well>.trimmer_stats.csv` otherwise becomes its own well,
  which turned one real 480-well upload into 960. Rules cover wrong bucket or project, wrong
  path depth, missing `.trim` infix, a doubled leading wafer token (`duplicated_wafer_token`),
  bulk-download leftovers (`non_sequencing_artifact`), sublibrary or wafer disagreement, a bad well
  token, and a missing sublibrary type. `repeated_token` looks only inside the sublibrary/well portion of a matched stem:
  bagging the whole name put the type token beside the sublibrary name, so a sublibrary legitimately
  called `REF3_GEX_P01` under type `GEX_hash_oligo` reported a repeat on an otherwise clean name —
  and through the rename gate that flipped a whole sublibrary to `UNKNOWN`. `expected_name` carries the corrected basename repairing that rule alone.
  When a vendor delivery has been located, a missing sublibrary type is filled from it. Written to
  `{order}_raw_sop_violations.csv`; violations are non-fatal.
- **The sublibrary folder may be spelled two ways, and both are clean.** Either in full
  (`REF3/raw/REF3_P05_2/`) or with the redundant `{ExperimentID}_` prefix elided
  (`CHEM16/raw/P03/` for the sublibrary the filenames call `CHEM16_P03`). The ExperimentID is
  already an ancestor segment, so the prefix carries nothing the path does not; the filename holds
  the authoritative name; and wells match across buckets on `(wafer, UG)`, so no check depends on
  which spelling was used. Every real trimmed upload measured elides it on *every* sublibrary —
  REF3 `P04_1`…`P07_1`, GENE7 `P02`…`P10`, CHEM16 `P03`…`P07` — and not one mixes the two.
  Demanding the full form was `sublibrary_folder_truncated`, which made those folders the whole of
  GENE7's finding set: 9 SOP rows, 864 wells `RENAMEABLE` and a proposed move for all 5184 objects
  of an upload that is fine. The rule is gone, along with the `folder` scope and the
  `expected_folder` column it alone populated. `sublibrary_mismatch` still fires for a filename
  neither spelling explains, and the rename mapping now tracks the folder segment and the
  filename's sublibrary name separately — collapsing them is what proposed the move.
- **Cross-bucket trimming completeness** compares the upload against the untrimmed vendor deliveries
  found by searching the roots in the notebook's `untrimmed_search_roots` parameter — a *list* of
  projects or buckets to look in, not of orders to trust, since one experiment spans several
  Novogene orders and one order can hold several experiments. Wells match on `(wafer, UG)`,
  measured unique on a real delivery (48 CRAMs, 48 distinct UGs per wafer). The vendor layout is
  `{project}/{order}/{ExperimentID}/raw/{wafer}/`; note the segment before `raw` is the ExperimentID
  alone and carries **no** sublibrary, so the authoritative sublibrary comes from the vendor filename
  with its trailing well token stripped. `qa_seahub_source.py` builds and merges the per-well indexes;
  `qa_seahub_recon.py` classifies them (`not_trimmed`, `orphan_trimmed`, `identity_mismatch`,
  `read_count_mismatch`, `size_not_reduced`, `duplicate_source_well`, `duplicate_trimmed_well`,
  `source_prefix_empty`, `overlapping_source_prefix`, `metadata_unavailable`) into
  `{order}_trimming_completeness.csv`, with per-order coverage in
  `{order}_seahub_source_coverage.csv`. Coverage matters: wells of an order missing from the list can
  only surface as `orphan_trimmed`, which reads as a completeness failure rather than incomplete
  input, so unsourced wafers are called out explicitly. `reconcile_trimming` copies the coverage rows
  it tallies into rather than mutating the caller's, so re-running the recon cell alone — without
  re-running the indexing cell above it — reports the same numbers rather than double-counting
  `matched` and flooring `unmatched` at zero. A vendor sidecar that cannot be read — absent, truncated, or
  carrying a non-numeric `read_count` — leaves that one well's count empty and reports it as
  `metadata_unavailable`, rather than raising out of the cell and losing the several hundred
  sidecars that were fine.
- **Rename mapping** (`qa_seahub_rename.py`) composes the per-rule repairs into one corrected S3 key
  per object in `{order}_seahub_rename_mapping.csv`. Advisory only — QA never moves, deletes or
  rewrites anything. Every proposal is re-validated against the SOP and withheld unless clean, two
  objects can never be given the same destination, and an existing destination is `blocked` rather
  than overwritten. `name_source` is `vendor` when the delivery supplied a missing token, else
  `inferred`. That self-check skips `upload`-scope rules: `bad_bucket` and `lab_project_mismatch`
  hold identically before and after the move, so testing them withheld proposals that strictly
  reduce the violation set — one wrong bucket turned every rename in the upload into `unresolved`
  and every well into `UNKNOWN`. They are filtered in the rename gate only; the SOP table still
  reports them, and the rename cell prints a banner. That banner is *not* gated on there being
  anything to move: because upload-scope rules are filtered out of the per-object check, an upload
  that is otherwise SOP-clean produces no rename rows at all, so gating on them inverted the
  severity — a clean upload in the wrong bucket reported every well `COMPLIANT` and wrote nothing to
  `errors.txt`, while adding a trivial naming defect made the location defect audible.
- **Per-well status** rolls the above up to one verdict per well in
  `{order}_seahub_well_status.csv`: `COMPLIANT`, `RENAMEABLE` (complete, every defect repairable by
  renaming), `DATA_GAP` (an artifact genuinely absent) or `UNKNOWN` (unidentifiable, or no corrected
  name derivable). A CRAM that arrived and is **0 bytes** is `DATA_GAP` too: completeness asks only whether the artifact
arrived, and the trimmed-vs-vendor size check skips a zero because it cannot tell "empty" from "size
not collected" — so an interrupted copy, which is exactly what leaves a zero-byte object, used to read
`COMPLIANT`. Only the CRAM is judged this way; an empty `.stdout` or `.stderr` is ordinary. Sizes are
collected in s3 mode only, and a key absent from the mapping reads as unknown rather than empty.
  `DATA_GAP` outranks the un-nameable `UNKNOWN` so a missing CRAM stays visible even
  when its vendor order is absent from the list, and the detail names the vendor key it can be
  re-trimmed from. Wells are keyed the same way here as in `index_trimmed_upload` — on the normalized
  stem, falling back to the raw one when normalizing makes it unparseable — so the two cannot
  disagree about whether a well was uploaded. A vendor well with *nothing* uploaded is `DATA_GAP` too — the largest gap there
  is, and as `UNKNOWN` it was the only kind the notebook left out of `errors.txt`, which writes
  `DATA_GAP` rows alone. Its row lists all five required artifacts, names the vendor key, and takes
  its `sublibrary` and `well` off the vendor stem — every other row derives those from the uploaded
  objects, and this is the row that has none, so they would otherwise be blank on the one verdict an
  operator most needs to locate. Note the two sources can disagree about the same sublibrary: an
  uploaded row carries the folder the upload actually used, which in every real upload elides the
  ExperimentID prefix (`P07_1`), while a gap row carries the vendor's SOP name for it
  (`REF3_P07_1`). Both spellings are clean, so that difference is permanent rather than something a
  rename resolves — filter the CSV on wafer and UG rather than on one spelling of the sublibrary. Whether a well is renameable is decided by `RenameProposal.renameable` — the same
  predicate that decides whether an object gets a destination in the rename CSV — so the two headline
  CSVs cannot contradict each other. Testing the defect *names* against `SEAHUB_RENAMEABLE_SOP_TYPES`
  instead is how they once did: a `sublibrary_mismatch` the vendor delivery had already repaired
  produced objects the rename CSV said to move inside a well the status CSV called `UNKNOWN`. That
  set is still maintained as documentation of which defects a rename can fix, and a test reads the
  defect literals out of the module to keep it in step.

---

## Finding the untrimmed deliveries

`untrimmed_search_roots` names where to look, and the wafers in the trimmed upload
are what is looked for. A **project root** is the value to use —
`s3://czi-novogene/{proj}` — because it needs only prefix-scoped `ListBucket` and
costs roughly a tenth of the listings a bare bucket does. A bucket root is legal
for the rare cross-project case.

The search is a delimiter descent, one level at a time, stopping at a child folder
named `raw`. That is what keeps it O(folders) rather than O(objects): the object
indexer paginates *flat*, which is why `normalize_source_uris` refuses a prefix
shallower than `{project}/{order}` and why a root goes through
`normalize_search_roots` instead. Discovery only ever hands the indexer a concrete
`.../raw/`, so that guard becomes an invariant check on the descent rather than
something the new path had to route around.

Wafers are read three ways, and each rescues a case the others miss:

| reading | rescues |
| --- | --- |
| the identity index (filename) | a well whose folder is not a wafer, or disagrees with its filename |
| the object path (folder) | a wafer whose filenames do not parse — a real delivery carries 192 CRAMs with no `Z####` UG at all |
| `discovered_wafers` (folder walk) | a wafer directory holding nothing QA ingested; empty in manifest mode |

Tokens that are not wafer-shaped are reported as `wafer_seed_rejected` rather than
searched for. Where a folder and a filename disagree, **both** are searched: which
is right is not decidable there, and the wrong one costs one lookup and reports as
not found, whereas choosing silently compares against the wrong delivery.

**Sibling expansion**, `untrimmed_search_siblings`, indexes every wafer in a
delivery that a seed pointed at rather than only the seeded ones. Leave it on. A
wafer the lab never trimmed cannot be a seed, so without expansion it is never
listed, never indexed, and produces neither a `not_trimmed` row nor a vendor-only
`DATA_GAP` well — and `DATA_GAP` is the only verdict written to `errors.txt`. A
whole missed plate would read as a clean upload.

What expansion cannot reach is a delivery that *no* uploaded wafer points at. The
descent walks past those anyway, so each is reported as
`unseeded_vendor_delivery` with its ExperimentID segment: if the segment is one
you recognise, its wells were never uploaded at all. The residual blind spot after
that is a delivery in a **different project or bucket** than the roots cover, which
is the one case only a wider root can find.

**One level, deliberately.** A delivery's wafers are read from what sits directly
inside `raw/` — child folders and loose objects — while the indexer that follows
paginates recursively. A vendor layout that nested a level deeper than
`raw/{wafer}/` would therefore be located but not recognised, and reported as
`unseeded_vendor_delivery` rather than indexed. No observed delivery does this, and
the failure is loud rather than silent: the wafers it holds also report
`wafer_not_found`, which reaches `errors.txt`. If such a layout ever appears, that
pair of findings is the signal to widen the read.

Three bounds keep the walk finite, and each is reported rather than raised, because
a partial answer an operator can see beats a dead cell: `SEAHUB_SEARCH_MAX_DEPTH`
levels below a root, `SEAHUB_SEARCH_MAX_LISTINGS` calls in total, and a refused
listing recorded with its error code and stepped over so one inaccessible project
does not cost the rest. `WaferSearchPlan.complete` asks all three at once, and
`summary()` is the line the notebook prints.

---

## Why the vendor ExperimentID is not the key

The vendor index used to be filtered to the ExperimentID being QA'd, matching a
pre-`raw` path segment against `ctx.order` with `segment == id or
segment.startswith(f"{id}_")`. That is recorded here because the reasoning still
explains the shape of what replaced it, and because the rule was subtler than it
looked.

It existed for a real reason: a listed prefix was order-level, one order holds
several experiments, and without the filter the other experiments' wells were
indexed too — each then having no counterpart in the upload and reporting as
`not_trimmed` plus an `UNKNOWN` well. The `_`-prefix arm was equally deliberate,
measured against a `GENE7_reupload` delivered alongside `REF3` in one order, and
against the older `{ExperimentID}_{sublibrary}` shape (`REF5_P01`) that a bare
equality test would have orphaned entirely.

It was still the wrong key, in three ways:

- **It failed on every plausible rename except the one it was written for.**
  `REF3_reupload` matched. `REF3-reupload`, `REF3reupload`, `reupload_REF3`,
  `ref3` and a renumbered `REF2` did not, and each one orphans a whole upload —
  surfaced only as a `source_prefix_spans_experiments` note somebody had to read.
- **It silently matched things it should not.** `REF3_2`, `REF3_v2` and `REF3_redo`
  all matched `REF3` with no finding at all, and the shape is used for genuinely
  different things: the real segment `RNA3_098` matches `RNA3`.
- **It did nothing on two-thirds of real layouts.** `derive_source_experiment`
  required an exact six-segment key, which only the foldered Novogene layout has;
  the other four real vendor prefixes fell through to `unplaced_wells`, were kept,
  and produced a `source_experiment_unreadable` row instead of any filtering.

What replaces it is the **seed set**: a delivery is indexed only for the wafers
the upload actually mentions. That is a positive assertion about contents rather
than a comparison of names, so it cannot be defeated by a rename and cannot
accidentally admit a neighbour. The ExperimentID segment survives as an advisory
label — it names a delivery in the coverage table and in
`unseeded_vendor_delivery` — and never decides whether one is indexed.

`source_prefix_spans_experiments` and `source_experiment_unreadable` retired with
the filter. The check they served did not: a foreign experiment's wells staying
out of the comparison is now pinned against a delivery discovered under the same
project root.

---

## Re-running a report offline

Numbers quoted for a real upload can be re-derived without S3 access, from a saved listing of the
bucket — any CSV/TSV with a column of `s3://` URIs, such as a console export:

```bash
python -m tests.seahub_offline_report LISTING.tsv --mode both
```

It drives the same production code the notebook does, with `MockS3Client` standing in for S3, and
prints the SOP rows and rules, the per-well verdicts, the rename statuses and a parity verdict for
the two `data_source` modes. Nothing identifies the upload on the command line: the ExperimentID
comes out of the listing via `resolve_qa_run_context`, which is also a live check that manifest mode
still recovers it. Add `--vendor LISTING [LISTING ...]` to reconcile against the untrimmed
deliveries, and `--json` for a diffable snapshot — `diff before.json after.json` is how a change is
checked against real shapes before it goes near a bucket.

The three uploads this mode was built against, as of this branch:

| listing | objects | SOP rows | rules | wells |
| --- | --- | --- | --- | --- |
| REF3 | 2037 | 885 | 288 each of `duplicated_wafer_token`, `invalid_sublibrary_type`, `missing_trim_infix`; 21 `non_sequencing_artifact` | 336 — 48 `COMPLIANT`, 288 `UNKNOWN`; 48 `COMPLIANT` and 288 `RENAMEABLE` with both vendor orders |
| GENE7 | 5185 | 1 | 1 `bad_path_depth` | 864, all `COMPLIANT` |
| CHEM16 | 1450 | 14 | 10 `non_sequencing_artifact`, 3 `unexpected_suffix`, 1 `no_recognized_artifacts` | 480, all `DATA_GAP` — every artifact misspelled, so no well is *nameable* |

Accepting the elided folder prefix is what moved those numbers: REF3 shed 7 SOP rows and its 48
folder-only wells went `RENAMEABLE` to `COMPLIANT` (2016 moveable objects to 1728 with both vendor
orders), and GENE7 went from 9 rows and 864 `RENAMEABLE` wells to 1 row and none, dropping all 5184
of its proposed moves. CHEM16 did not move at all: every artifact is misspelled, so no name ever
reached the folder question. Every other field in the snapshot — including all of the vendor
reconciliation — is unchanged.

All three report identically from an S3 listing and from a manifest of the same keys. The one field
that cannot match is `discovered_wafers`: wafers are found by walking folders, which only s3 mode
does, so the wafer summary is populated there and empty under a manifest. It is excluded from
`PARITY_FIELDS` for that reason rather than silently passing.

---

## Checksums

Checksums are deliberately out of scope. S3 validates integrity on upload when the client sends a
checksum, and a multipart `ETag` is a composite of per-part digests rather than a whole-file hash, so
it cannot be compared against a lab-supplied manifest. The trimmed-vs-vendor size comparison
(`size_not_reduced`) covers the cheap part of the same question without transferring any bytes.
