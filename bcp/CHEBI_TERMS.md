## ChEBI terms

Resolve **ChEBI identifiers** to their authoritative name and synonyms, and check them for correctness, via the [ChEBI REST API](https://www.ebi.ac.uk/chebi/). Part of BCP tooling under `bcp/chebi_terms`.

This is the inverse of [CHEBI_LOOKUP.md](CHEBI_LOOKUP.md), which starts from a CAS Registry Number and has no ChEBI ID yet. Use `chebi_terms` when you **already hold a ChEBI ID** — from a submitter, a curator, or a legacy spreadsheet — and need to know whether it is real, current, and the right compound.

No local data cache and no extra dependencies: one HTTP GET per unique ID.

---

## CLI entry point

Run from the **`bcp`** directory:

```bash
cd bcp
python -m chebi_terms --help
```

Provide **either** `--input` (batch CSV) **or** `--chebi` (single lookup). Both `CHEBI:16236` and bare `16236` are accepted.

### Single-ID mode

| Flag | Description |
|------|-------------|
| `--chebi` | One ChEBI ID (no CSV) |
| `--expect-name` | Name to check against ChEBI's name and synonyms |
| `--expect-cas` | CAS number to check against ChEBI's own CAS xrefs |
| `--max-synonyms` | Cap synonyms reported (default: no cap) |
| `--format` | `json` (default) or `csv` |
| `--output`, `-o` | Output file (default: stdout) |
| `-v`, `--verbose` | Debug logging |

```bash
# name + synonyms, JSON to stdout
python -m chebi_terms --chebi CHEBI:16236

# verify an ID against a recorded name and CAS
python -m chebi_terms --chebi 18484 --expect-name "(-)-epicatechin" --expect-cas 490-46-0
```

Progress and summary messages go to **stderr** in single-ID mode so stdout stays machine-readable for JSON.

### Batch mode (CSV)

| Flag | Description |
|------|-------------|
| `--input`, `-i` | Input CSV |
| `--chebi-column` | ChEBI ID column name (default: `chebi_id`) |
| `--name-column` | Optional: column holding the recorded name to check |
| `--cas-column` | Optional: column holding the recorded CAS to check |
| `--max-synonyms` | Cap synonyms reported (default: no cap) |
| `--output`, `-o` | Output CSV (default: `<input_stem>_chebi_checked.csv`) |
| `-v`, `--verbose` | Debug logging |

```bash
python -m chebi_terms --input curated.csv \
  --chebi-column chebi_id --name-column compound_name --cas-column CAS \
  --output curated_checked.csv
```

---

## Input (batch mode)

- **Format:** CSV with `utf-8-sig` encoding.
- **Required column:** ChEBI ID (configurable via `--chebi-column`).
- All other columns are preserved in the output. An input column whose name collides with one this tool writes is overwritten, with a warning. Pointing `--chebi-column`, `--name-column`, or `--cas-column` at such a name is an error: for the ID column it would emit two columns of the same name, and for the check columns it would erase the very value being checked, leaving a verdict with no record of what it was compared against.

---

## Output

Appended columns (batch mode) or fields (single-ID mode):

- `chebi_accession` — the **primary** accession, canonicalized
- `chebi_name`
- `chebi_synonyms` (pipe-separated)
- `chebi_definition`
- `chebi_stars` — 3 means manually curated by ChEBI
- `id_status`
- `name_verdict`
- `cas_verdict`

Single-ID JSON includes a top-level `"chebi_id"` key (the ID as supplied) plus these fields.

### `id_status`

**What ChEBI says about the compound:**

| Value | Meaning |
|-------|---------|
| `ok` | Primary and released |
| `secondary` | The ID is a secondary/merged ID; `chebi_accession` holds the primary it resolves to |
| `not_released` | ChEBI has the record but has not released it |
| `not_found` | No such compound |

When more than one of these could apply, precedence is `not_found` → `not_released` → `secondary` → `ok`.

**What the input looks like** — neither costs a request:

| Value | Meaning |
|-------|---------|
| `missing` | The ID cell was empty |
| `invalid` | Not a ChEBI identifier at all |

**Whether we managed to ask:**

| Value | Meaning |
|-------|---------|
| `lookup_failed` | ChEBI could not be reached, or kept failing. **Says nothing about the compound.** |

`lookup_failed` is deliberately not `not_found`. A DNS blip or a ChEBI outage must never render as "ChEBI says this compound does not exist" — that is the one output of this tool someone is most likely to act on destructively. See [When ChEBI is unreachable](#when-chebi-is-unreachable).

### `name_verdict` and `cas_verdict`

| Column | Values |
|--------|--------|
| `name_verdict` | `match` · `synonym_match` · `mismatch` · `not_checked` |
| `cas_verdict` | `confirmed` · `not_in_chebi` · `no_cas_recorded` · `not_checked` |

Both read `not_checked` unless you supply something to check against.

`no_cas_recorded` means ChEBI holds no CAS number for that compound at all, so it is not contradicting you — it has nothing to compare against. Split out of `not_in_chebi` for the same reason `lookup_failed` is not `not_found`: silence is not a negative answer. Plenty of ChEBI entries have no CAS xref, so without the distinction a perfectly good pairing would read as a disagreement.

Name comparison is deliberately permissive: markup is stripped, case is folded, Unicode dash variants are normalized to `-`, and the value is matched against the official name *and* every synonym type — including brand names, INNs, and ChEBI's ASCII spellings. So `(-)-epicatechin` matches ChEBI's `(−)-epicatechin` (U+2212), and `alcohol etilico` matches `alcohol etílico`. Output is stricter than matching: only English `SYNONYM` and `IUPAC NAME` entries are emitted, because acetylsalicylic acid alone carries 64 brand names.

CAS comparison normalizes the same Unicode dash variants and strips stray whitespace, so a pasted `64‑17‑5` is not reported as a disagreement. Internal dashes are kept: collapsing `64-17-5` and `64175` onto one key would risk a false `confirmed`, which is the more dangerous direction to be wrong in.

`cas_verdict` checks ChEBI's **own** CAS cross-references, not PubChem's. A `confirmed` verdict means ChEBI independently agrees that the CAS and the ChEBI ID describe the same compound; `not_in_chebi` means ChEBI records CAS numbers for it and yours is not among them.

---

## When ChEBI is unreachable

A 404 is an answer — it means no such compound, and it is cached. Anything else is *not* an answer: a connection error, a timeout, a persistent 5xx, exhausted 429 backoff, a terminal 4xx such as 400/403/410 (failed immediately, since retrying will not change it), a 200 carrying a non-JSON body, or a payload missing any field the verdicts depend on (`chebi_accession`, `is_released`, `name`).

- The row gets `id_status: lookup_failed`, blank facts, and both verdicts `not_checked`.
- Nothing is cached, so a later row carrying the same ID is retried rather than served a stale failure.
- After **5 consecutive** failed lookups the run aborts with a clear error rather than burning 6s of backoff on every remaining row — plus up to three 15s connect timeouts each, if packets are being dropped rather than refused. The partial output is left in place. Only lookups that actually reached the network count: a cache hit does not reset the counter, since it is no evidence the outage is over.
- The command **exits 1** if any row ended up `lookup_failed`, so a wrapper script cannot mistake a degraded run for a clean one.

Exit codes distinguish "the tool could not do its job" from "the tool did its job and the answer was unwelcome". A `mismatch`, `not_found`, or `invalid` verdict exits **0** — that is a successful run reporting bad data.

### If the endpoint moves

A 404 on every row would otherwise be indistinguishable from a sheet of compounds that all genuinely do not exist — a clean tally and exit 0. EBI has already relocated the flat files and retired the SOAP service, so a rename here is a live possibility.

So when **every** ID that reached ChEBI came back `not_found` and none resolved, across at least 5 **distinct IDs**, the run logs an error naming the endpoint and exits **1**. Per-row `id_status` stays honest either way; only the run-level verdict changes. It also fires on a batch of genuinely bogus IDs — equally worth a human look, so the false positive is harmless.

Counted per ID rather than per row, because caching means one request can serve many rows: five rows carrying a single retired ID make one request, and "five distinct compounds all happen not to exist" is the implausible thing — not "one ID repeated five times". Those five rows are reported `not_found` and do **not** trip the guard.

### If the column is wrong

The same mistake from the cheaper direction: point `--chebi-column` at a name column or a notes column and every row comes back `invalid` without a single request. So when **nothing** reached ChEBI at all and at least 5 **distinct invalid values** were seen, the run names the column it was given and exits **1**.

Distinct values rather than rows, because 500 rows of one junk string is a single mistake rather than evidence about the column.

An entirely **blank** ID column does not trip this. A sheet whose ChEBI IDs have not been filled in yet is a plausible real input rather than an operator error, so those rows report `missing` and the run still exits 0. Only junk values — the signature of a column holding something other than ChEBI IDs — fail it.

A single `lookup_failed` disqualifies this diagnosis. Failures are never cached, so an outage also leaves both ID counts at zero and would otherwise satisfy the "nothing reached ChEBI" premise for the wrong reason — blaming the column for what the network did. The run still exits 1 either way; only the message changes.

**Exit 1** therefore means: any row `lookup_failed`, every reachable lookup missed, or no usable ID and ≥5 distinct junk values. Anything else — including a `mismatch`, a genuine `not_found`, or an unfilled ID column — exits **0**. `RunSummary.degraded` is the same predicate, for programmatic callers.

---

## Why the verdicts matter

- **Stale IDs resolve silently.** `CHEBI:18484` returns `CHEBI:90` — `(−)-epicatechin`. Without the `secondary` flag you would never notice the ID was outdated.
- **A wrong ID is usually still a valid ID.** A transposed digit lands on a real but unrelated compound, so `id_status: ok` proves nothing on its own. The name check is what catches it.
- **CAS and ChEBI can disagree.** `not_in_chebi` flags a pairing that ChEBI itself does not record.

---

## ChEBI API usage

Per **unique** ID, one REST call:

```
GET https://www.ebi.ac.uk/chebi/backend/api/public/compound/{numeric_id}/
```

Batch mode caches one payload per unique ID and evaluates each row's expectations against it, so a CSV with the same compound on many rows costs one request. Only definitive answers are cached — see [When ChEBI is unreachable](#when-chebi-is-unreachable). The client sleeps **0.25s** after each call and backs off on 429/503; ChEBI publishes no documented rate limit, so this mirrors the spacing used for PubChem in `chebi_lookup`.

Empty and malformed IDs are resolved locally and never hit the network.

---

## Package layout

```
bcp/chebi_terms/
├── __init__.py
├── __main__.py
├── cli.py              # argparse entrypoint
├── io.py               # CSV batch + single-ID output
├── client.py           # ChEBI HTTP client, parsing, verdicts
└── record_fixtures.py  # record live ChEBI JSON for tests
```

**Programmatic use:**

```python
from chebi_terms import describe_chebi_id, verify_chebi_id

facts = describe_chebi_id("CHEBI:16236")
print(facts["chebi_name"], facts["chebi_synonyms"])

check = verify_chebi_id("CHEBI:18484", expected_cas="490-46-0")
print(check["id_status"], check["chebi_accession"], check["cas_verdict"])
# secondary CHEBI:90 confirmed
```

`client.verify_payload()` is pure — hand it an already-fetched payload to evaluate expectations without making a request. Because it cannot raise, it trusts the payload's shape: anything from `fetch_compound()` has already been checked, but a payload you obtained some other way should go through `client.check_payload_shape()` first, or a changed payload will read as a verdict about the compound.

`describe_chebi_id()` and `verify_chebi_id()` never raise: an unreachable ChEBI comes back as `lookup_failed`. The lower-level `fetch_compound()` and `get_with_retry()` raise `ChebiUnavailableError` instead, so callers building their own loop can tell the two apart.

`verify_chebi_file()` returns a `RunSummary` so a caller can spot a degraded run:

```python
summary = verify_chebi_file(Path("curated.csv"), "chebi_id", Path("out.csv"))
summary.status_counts["lookup_failed"]  # Counter, per-status row tally
summary.missed_ids  # distinct IDs that reached ChEBI and 404'd
summary.resolved_ids  # distinct IDs that came back
summary.invalid_values  # distinct unparseable values seen
summary.missing_rows  # rows with an empty ID cell
summary.name_mismatches  # rows whose recorded name disagreed
summary.cas_disagreements  # rows whose CAS is not among the ones ChEBI records
summary.cas_not_recorded  # rows ChEBI holds no CAS for, so unverifiable
summary.suspect_endpoint  # every reachable lookup missed
summary.suspect_column  # nothing was usable as a ChEBI ID
summary.degraded  # any of the above — mirrors exit 1
```

---

## Recording test fixtures

To refresh committed ChEBI response snapshots (e.g. after API shape changes):

```bash
cd bcp
python -m chebi_terms.record_fixtures
```

Records ethanol, a secondary ID, and a markup-heavy name by default; pass `--chebi` (repeatable) for others. Writes JSON under `bcp/tests/fixtures/chebi_terms/chebi_live/{numeric_id}.json`, keyed by the **requested** ID so secondary-ID behaviour stays reproducible offline. Commit the updated files.

---

## Testing

**Default (offline, mocked HTTP):**

```bash
cd bcp
pytest tests/test_chebi_terms.py tests/test_chebi_terms_cli.py -v
```

**Live ChEBI (opt-in, requires network):**

```bash
pytest -m chebi tests/test_chebi_terms_chebi_live.py -v
```

Live tests are excluded from default `pytest` runs (see `pytest.ini`: `not e2e and not pubchem and not chebi`).
