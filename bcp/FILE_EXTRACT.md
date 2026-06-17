## File extract

Extract S3 metadata for **FASTQ.gz** deliverables and **Cell Ranger h5** matrices. Part of BCP tooling under `bcp/file_extract`.

Requires AWS credentials with read access to the target bucket (standard `boto3` credential chain).

---

## CLI entry point

Run from the **`bcp`** directory:

```bash
cd bcp
python -m file_extract --help
python -m file_extract fastq --help
python -m file_extract h5 --help
```

### FASTQ mode

Walk an S3 order prefix and emit one row per deliverable `.fastq.gz` with CRC64NVME and companion `-metadata.json` `read_count`.

```bash
python -m file_extract fastq s3://example-bucket/path/to/order
python -m file_extract fastq s3://example-bucket/path/to/order -o order_fastq_info.tsv
python -m file_extract fastq s3://example-bucket/path/to/order --no-require-raw
```

| Flag | Description |
|------|-------------|
| `-o`, `--output` | Output TSV (default: `<order>_fastq_info.tsv`) |
| `--no-require-raw` | Don't require `/raw/` in the S3 key |
| `--workers` | Process pool size (default: min(64, n_files)) |
| `--retries` | Max attempts per transient S3 error (default: 5) |
| `--strict` | Exit 1 if any per-file CRC or metadata fetch fails |
| `-v`, `--verbose` | Debug logging |
| `-q`, `--quiet` | Disable progress bars |

**Output columns:** `filename`, `s3_uri`, `read`, `lane`, `size_bytes`, `crc64nvme_base64`, `read_count`, `crc_error`, `metadata_error`

**Guardrails:** warns when R1 and R2 file counts differ; suggests `--no-require-raw` when zero files match.

### H5 mode

Point at a `per_sample_outs` prefix. By default matches `sample_filtered_feature_bc_matrix.h5`, fetches CRC64NVME, and introspects matrix shape (cell count, feature types).

```bash
python -m file_extract h5 s3://example-bucket/.../outs/per_sample_outs
python -m file_extract h5 s3://.../per_sample_outs --no-introspect
python -m file_extract h5 s3://.../per_sample_outs --genome --metrics
```

| Flag | Description |
|------|-------------|
| `-o`, `--output` | Output TSV (default: `<run-or-dir>_h5_info.tsv`) |
| `--target-filename` | h5 basename to match (default: `sample_filtered_feature_bc_matrix.h5`) |
| `--no-introspect` | Checksums and listing only |
| `--genome` | Add `gene_counts_by_genome` JSON column |
| `--metrics` | Cross-check against sibling `metrics_summary.csv` |
| `--workers` | Thread count (default: 16 with introspection, 64 without) |
| `--retries` | Max attempts per transient S3 error (default: 5) |
| `--strict` | Exit 1 if any per-file enrichment fails |
| `-v`, `--verbose` | Debug logging |
| `-q`, `--quiet` | Disable progress bars |

**Optional introspection dependencies** (only needed without `--no-introspect`):

```bash
pip install h5py fsspec s3fs
```

---

## Testing

```bash
cd bcp
pytest tests/test_file_extract_*.py -v
pytest --cov=file_extract tests/test_file_extract_*.py
```

All tests use mocked S3; no AWS credentials required for the default suite.

---

## Migration from prototypes

The standalone prototypes in `bcp/docs/file_extractor.py` and `bcp/docs/extract_h5.py` are superseded by this package. Use `python -m file_extract fastq` and `python -m file_extract h5` instead.
