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

**Arguments:**

| Flag | Description |
|------|-------------|
| `--input`, `-i` | Input CSV (required) |
| `--cas-column`, `-c` | CAS column name (default: `CAS`) |
| `--output`, `-o` | Output CSV (default: `<input_stem>_chebi_mapped.csv`) |
| `-v`, `--verbose` | Debug logging |

**Example:**

```bash
python -m chebi_lookup --input chemicals.csv --cas-column "CAS Number" --output mapped.csv
```

---

## Input

- **Format:** CSV with `utf-8-sig` encoding.
- **Required column:** CAS Registry Number (configurable via `--cas-column`).
- All other columns are preserved in the output.

---

## Output

Appended columns (after original CSV columns):

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

---

## PubChem API usage

Per CAS row, the tool makes up to **3 REST calls**:

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
├── cli.py      # argparse entrypoint
├── io.py       # CSV batch read/write
└── client.py   # PubChem HTTP client
```

**Programmatic use:**

```python
from chebi_lookup import lookup_cas

result = lookup_cas("50-00-0")
print(result["chebi_id"], result["preferred_name"])
```

---

## Testing

```bash
cd bcp
pytest tests/test_chebi_lookup.py tests/test_chebi_lookup_cli.py -v
```

Tests mock HTTP; no live PubChem calls in the default suite.
