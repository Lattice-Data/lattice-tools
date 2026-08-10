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
- All other columns are preserved in the output. An input column whose name collides with one this tool writes is overwritten, with a warning; naming the *ID* column one of them is an error.

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

| Value | Meaning |
|-------|---------|
| `ok` | Primary and released |
| `secondary` | The ID is a secondary/merged ID; `chebi_accession` holds the primary it resolves to |
| `not_released` | ChEBI has the record but has not released it |
| `not_found` | No such compound |
| `invalid` | Not a ChEBI identifier at all — no request is made |

When more than one could apply, precedence is `not_found` → `not_released` → `secondary` → `ok`.

### `name_verdict` and `cas_verdict`

| Column | Values |
|--------|--------|
| `name_verdict` | `match` · `synonym_match` · `mismatch` · `not_checked` |
| `cas_verdict` | `confirmed` · `not_in_chebi` · `not_checked` |

Both read `not_checked` unless you supply something to check against.

Name comparison is deliberately permissive: markup is stripped, case is folded, Unicode dash variants are normalized to `-`, and the value is matched against the official name *and* every synonym type — including brand names, INNs, and ChEBI's ASCII spellings. So `(-)-epicatechin` matches ChEBI's `(−)-epicatechin` (U+2212), and `alcohol etilico` matches `alcohol etílico`. Output is stricter than matching: only English `SYNONYM` and `IUPAC NAME` entries are emitted, because acetylsalicylic acid alone carries 64 brand names.

`cas_verdict` checks ChEBI's **own** CAS cross-references, not PubChem's. A `confirmed` verdict means ChEBI independently agrees that the CAS and the ChEBI ID describe the same compound.

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

Batch mode caches one payload per unique ID and evaluates each row's expectations against it, so a CSV with the same compound on many rows costs one request. The client sleeps **0.25s** after each call and backs off on 429/503; ChEBI publishes no documented rate limit, so this mirrors the spacing used for PubChem in `chebi_lookup`.

Malformed IDs are rejected locally and never hit the network.

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

`client.verify_payload()` is pure — hand it an already-fetched payload to evaluate expectations without making a request.

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
