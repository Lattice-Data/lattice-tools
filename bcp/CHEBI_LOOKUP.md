## ChEBI lookup

Map **CAS Registry Numbers** to ChEBI identifiers and PubChem compound properties via the [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest). Part of BCP tooling under `bcp/chebi_lookup`.

ChEBI IDs are extracted from PubChem cross-references when present. A compound may resolve to a PubChem CID without a ChEBI xref.

To go the other way — from a ChEBI ID you already hold to its authoritative name, synonyms, and a correctness verdict — see [CHEBI_TERMS.md](CHEBI_TERMS.md).

When the question is whether a whole sheet's identifiers agree with each other **by structure** rather than by name, see [STRUCTURE_CHECK.md](STRUCTURE_CHECK.md): it resolves the CAS, the name, and the ChEBI ID independently to InChIKeys and compares those.

When a whole batch has **no trustworthy identifier at all** — only a compound name and a CAS number per row — see [CHEBI_IDENTIFICATION.md](CHEBI_IDENTIFICATION.md), the phased method that establishes one from scratch and builds on this tool.

### Request behaviour

`get_with_retry` keeps its signature and return semantics — a `Response`, or `None` for both "no such thing" and "could not ask". `get_with_retry_status` returns the same response alongside an outcome (`ok` / `not_found` / `unreachable`) for callers that need to tell those apart; `structure_check` does, because reporting an outage as "no such compound" would turn a network problem into a finding about the chemistry.

Three things changed for **all** callers when that was added:

- **Malformed-request 4xx are no longer retried.** `400`, `414` and `422` are the statuses a bad *cell value* can provoke, so asking twice more gets the same reply. Every other 4xx is deliberately excluded, because no cell can cause it: `401`/`403`/`407` is a blocked client, `405` the endpoint refusing the method, `410` the endpoint retired, `413` a proxy. Treating any of those as "no such compound" would misreport a whole run of good CAS numbers and blame the column for it.
- **No backoff sleep after the final attempt.** It delayed the next try, and there isn't one; this was ~8s per exhausted request.
- **`cas_to_cid` percent-encodes the CAS.** These values come from untrusted spreadsheets and land in the URL path, where a stray `/` or `?` would rewrite or truncate the request. Quoting happens inside the function so every caller is covered rather than each one remembering.

---

## CLI entry point

Run from the **`bcp`** directory:

```bash
cd bcp
python -m chebi_lookup --help
```

Provide **either** `--input` (batch CSV) **or** `--cas` (single lookup).

### Batch mode (CSV)

| Flag | Description |
|------|-------------|
| `--input`, `-i` | Input CSV |
| `--cas-column`, `-c` | CAS column name (default: `CAS`) |
| `--output`, `-o` | Output CSV (default: `<input_stem>_chebi_mapped.csv`) |
| `-v`, `--verbose` | Debug logging |

```bash
python -m chebi_lookup --input chemicals.csv --cas-column "CAS Number" --output mapped.csv
```

### Single-CAS mode

| Flag | Description |
|------|-------------|
| `--cas` | One CAS Registry Number (no CSV) |
| `--format` | `json` (default) or `csv` |
| `--output`, `-o` | Output file (default: stdout) |
| `-v`, `--verbose` | Debug logging |

```bash
# JSON to stdout
python -m chebi_lookup --cas 64-17-5

# JSON to file
python -m chebi_lookup --cas 64-17-5 -o result.json

# One-row CSV
python -m chebi_lookup --cas 64-17-5 --format csv -o result.csv
```

Progress and summary messages go to **stderr** in single-CAS mode so stdout stays machine-readable for JSON.

---

## Input (batch mode)

- **Format:** CSV with `utf-8-sig` encoding.
- **Required column:** CAS Registry Number (configurable via `--cas-column`).
- All other columns are preserved in the output.

---

## Output

Appended columns (batch mode) or fields (single-CAS mode):

- `pubchem_cid`
- `chebi_id`
- `preferred_name`
- `iupac_name`
- `molecular_formula`
- `molecular_weight`
- `isomeric_smiles` — PubChem's `SMILES` property (stereo-bearing)
- `canonical_smiles` — PubChem's `ConnectivitySMILES` property
- `inchi`
- `inchikey`
- `xlogp`
- `tpsa`
- `synonyms` (pipe-separated, capped at 20)

Single-CAS JSON includes a top-level `"CAS"` key plus these fields.

**The two SMILES column names predate PubChem's rename and are kept on purpose.** PubChem
retired `IsomericSMILES` in favour of `SMILES` and `CanonicalSMILES` in favour of
`ConnectivitySMILES`. The request asks by the new names; the output columns keep the old
ones, because downstream consumers read them. The mapping is not the one the names
suggest: `SMILES` is the **stereo-bearing** property and belongs in `isomeric_smiles`, so
putting it in `canonical_smiles` would quietly swap the two columns.

**Two tracked outputs still carry the empty columns the rename produced.**
`bcp/cas_batch_input_chebi_mapped.csv` and `bcp/cas_batch_input_set2_mapped.csv` were
generated while the code asked by the retired names *and* read the response by them, so
every lookup returned `""` for both SMILES fields with no error, no non-200 and no log
line. The empty `isomeric_smiles` and `canonical_smiles` cells in those files are that
bug, **not** PubChem's answer for those compounds. Regenerating them needs live PubChem —
about 4 calls for each of 224 rows — and has not been done; every other column in them is
unaffected.

---

## PubChem API usage

Per CAS, the tool makes up to **3 REST calls**:

1. CAS → CID (`/compound/name/{CAS}/cids`)
2. CID → properties + ChEBI xref (`/property/...` and `/xrefs/RegistryID`)
3. CID → synonyms

PubChem rate limits: ~5 req/s, 400 req/min. The client sleeps **0.25s** after each call (~1.3 req/s for a full row).

---

## Package layout

```
bcp/chebi_lookup/
├── __init__.py
├── __main__.py
├── cli.py              # argparse entrypoint
├── io.py               # CSV batch + single-CAS output
├── client.py           # PubChem HTTP client
└── record_fixtures.py  # record live PubChem JSON for tests
```

**Programmatic use:**

```python
from chebi_lookup import lookup_cas

result = lookup_cas("64-17-5")
print(result["chebi_id"], result["preferred_name"])
```

---

## Recording test fixtures

To refresh committed PubChem response snapshots (e.g. after API shape changes):

```bash
cd bcp
python -m chebi_lookup.record_fixtures --cas 64-17-5
```

Writes JSON under `bcp/tests/fixtures/chebi_lookup/pubchem_live/{CAS}/` (`cids.json`, `properties.json`, `registry_ids.json`, `synonyms.json`). Commit the updated files. Default offline tests prefer these recordings when URLs match.

---

## Testing

**Default (offline, mocked HTTP):**

```bash
cd bcp
pytest tests/test_chebi_lookup.py tests/test_chebi_lookup_cli.py -v
```

**Live PubChem (opt-in, requires network):**

```bash
pytest -m pubchem tests/test_chebi_lookup_pubchem_live.py -v
```

Live tests are excluded from default `pytest` runs (see `pytest.ini`: `not e2e and not pubchem`).
