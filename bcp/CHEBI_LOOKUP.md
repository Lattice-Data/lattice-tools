## ChEBI lookup

Map **CAS Registry Numbers** to ChEBI identifiers and PubChem compound properties via the [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest). Part of BCP tooling under `bcp/chebi_lookup`.

ChEBI IDs are extracted from PubChem cross-references when present. A compound may resolve to a PubChem CID without a ChEBI xref.

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
- `isomeric_smiles`
- `canonical_smiles`
- `inchi`
- `inchikey`
- `xlogp`
- `tpsa`
- `synonyms` (pipe-separated, capped at 20)

Single-CAS JSON includes a top-level `"CAS"` key plus these fields.

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
