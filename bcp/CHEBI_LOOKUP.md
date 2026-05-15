## ChEBI lookup

Batch lookup of ChEBI records from a file of identifiers. Part of the BCP tooling under `bcp/chebi_lookup`.

Implementation is in progress; this document will be expanded once the lookup logic is in place.

---

## CLI entry point

Run from the **`bcp`** directory (so Python resolves the package):

```bash
cd bcp
python -m chebi_lookup --help
```

**Planned arguments:**

- `--input`, `-i` — path to file with identifiers (one per line or CSV/TSV)
- `--output`, `-o` — optional output path (default: stdout)
- `--id-column` — column name for tabular input

---

## Package layout

```
bcp/chebi_lookup/
├── __init__.py
├── __main__.py
├── cli.py      # argparse entrypoint
├── io.py       # read identifier files
└── client.py   # web/API queries
```

Tests: `bcp/tests/test_chebi_lookup.py`, `test_chebi_lookup_cli.py`  
Fixtures: `bcp/tests/fixtures/chebi_lookup/`

---

## Testing

```bash
cd bcp
pytest tests/test_chebi_lookup.py tests/test_chebi_lookup_cli.py -v
```
