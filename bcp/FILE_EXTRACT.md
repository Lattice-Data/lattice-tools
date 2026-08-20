## File extract

Extract S3 metadata for **FASTQ.gz**, **CRAM**, **Cell Ranger h5**, and **Scale h5ad** matrices. Part of BCP tooling under `bcp/file_extract`.

Requires AWS credentials with read access to the target bucket (standard `boto3` credential chain).

---

## CLI entry point

Run from the **`bcp`** directory:

```bash
cd bcp
python -m file_extract --help
python -m file_extract fastq --help
python -m file_extract cram --help
python -m file_extract h5 --help
python -m file_extract scale_h5ad --help
```

### FASTQ mode

Walk an S3 order prefix and emit one row per deliverable `.fastq.gz` with CRC64NVME and companion `-metadata.json` `read_count`, shaped into the Lattice submission sheets.

```bash
python -m file_extract fastq s3://example-bucket/path/to/ORDER \
  --lab example-lab \
  --cro-order ORDER \
  --is-pilot-order false \
  --sequencing-platform "Ultima Genomics UG 100"
```

| Flag | Description |
|------|-------------|
| `--lab` | **Required.** `example-lab` or `/labs/example-lab/` |
| `--cro-order` | **Required.** Order identifier written to `CRO_order` |
| `--is-pilot-order` | **Required.** `true` or `false` (case-insensitive) |
| `--sequencing-platform` | **Required.** e.g. `"Ultima Genomics UG 100"` |
| `--alias-namespace` | Alias namespace (default: the lab name) |
| `--out-dir` | Directory for the sheet TSVs (default: alongside `--output`) |
| `-o`, `--output` | Diagnostic TSV (default: `<order>_fastq_info.tsv`) |
| `--no-require-raw` | Don't require `/raw/` in the S3 key |
| `--workers` | Process pool size (default: min(64, n_files)) |
| `--retries` | Max attempts per transient S3 error (default: 5) |
| `--strict` | Exit 1 if any per-file CRC or metadata fetch fails |
| `-v`, `--verbose` | Debug logging |
| `-q`, `--quiet` | Disable progress bars |

**Three files per run:**

| File | Purpose |
|------|---------|
| `<order>_SequenceFile.tsv` | One row per file, in submission-sheet column order |
| `<order>_SequenceFileSet.tsv` | One row per set of files sequenced together |
| `<order>_fastq_info.tsv` | Diagnostics: per-file errors and lookup helpers |

**Diagnostic columns:** `filename`, `s3_uri`, `read`, `lane`, `size_bytes`, `crc64nvme_base64`, `read_count`, `crc_error`, `metadata_error`, `sample_dir`, `set_alias`

**Guardrails:** warns when R1 and R2 file counts differ, when R2 has no R1, when a set's R1 and R2 report different `read_count` (a truncation signal), when a file carries no read designator, and when `--cro-order` does not appear in the S3 prefix. Exits 1 on an alias collision, before any per-file request is spent. Suggests `--no-require-raw` when zero files match.

### CRAM mode

Walk an S3 order prefix and emit one row per deliverable `.cram` under `/raw/` with CRC64NVME and companion `-metadata.json` `read_count`. CRAMs are single-file, so each one becomes its own set with `run_cardinality` `single-end`.

```bash
python -m file_extract cram s3://example-bucket/path/to/ORDER \
  --lab example-lab \
  --cro-order ORDER \
  --is-pilot-order false \
  --sequencing-platform "Ultima Genomics UG 100" \
  --cram-slot trimmed
```

Takes every flag from FASTQ mode, plus:

| Flag | Description |
|------|-------------|
| `--cram-slot` | **Required.** `trimmed` or `untrimmed` — which SequenceFileSet slot the deliverable CRAM fills. Deliberately has no default: the S3 layout does not say whether a delivered CRAM is trimmer output. |

**Diagnostic columns:** `filename`, `s3_uri`, `size_bytes`, `crc64nvme_base64`, `read_count`, `crc_error`, `metadata_error`, `sample_dir`, `set_alias`

**Guardrails:** the FASTQ guardrails that apply, plus warnings when only unmatched CRAMs are found (excluded from output) and when `.ucram` files are present.

---

## Submission sheets

Both sheets carry the full column list of their Lattice tab, including the server-written `#response`, `#response_time` and `uuid` columns as empty leading cells, so a run's rows **paste into the workbook at A2** without reordering. Rows are emitted in S3-key order with R1 next to its R2, and two runs over the same order produce identical files.

Alias arrays are written unquoted (`["example-lab:file.fastq.gz"]`). An unquoted array survives both a spreadsheet import and a raw-text paste; the quoted form only survives the former.

### SequenceFile

`#response` · `#response_time` · `uuid` · `aliases` · `lab` · `file_format` · `s3_uri` · `crc64nvme_base64` · `no_file_available` · `derived_from` · `status` · `submitter_comment` · `description` · `read_count` · `file_size`

| Column | Value |
|--------|-------|
| `aliases` | `["<namespace>:<filename>"]`, extension included |
| `lab` | `/labs/<lab>/` |
| `file_format` | `fastq` or `cram` |
| `no_file_available` / `status` | `FALSE` / `current` |
| `read_count` | From the file's companion `-metadata.json` |
| `file_size` | From the S3 listing |

A failed checksum or `read_count` fetch leaves the cell empty and records the reason in the diagnostic TSV, so a partial run still pastes cleanly.

### SequenceFileSet

`#response` · `#response_time` · `uuid` · `aliases` · `lab` · `library` · `run_cardinality` · `status` · `submitter_comment` · `description` · `CRO_order` · `is_pilot_order` · `read1` · `read2` · `read3` · `index1` · `index2` · `untrimmed_cram` · `trimmed_cram` · `sequencing_platform`

| Column | Value |
|--------|-------|
| `aliases` | `["<namespace>:<stem>"]`, the basename minus `_R<n>_<chunk>.fastq.gz` (or minus `.cram`) |
| `run_cardinality` | `paired-end` when R2 is present, `single-end` otherwise. R3 and index files occupy their own slots and do not change the term. |
| `CRO_order` / `is_pilot_order` / `sequencing_platform` | From the required flags |
| `read1` … `index2` | Member file aliases, **bare** rather than JSON arrays |
| `untrimmed_cram` / `trimmed_cram` | The CRAM's alias, per `--cram-slot` |

**One row schema covers both delivery formats; a single order never uses both.** The read slots and the CRAM slots share this sheet because one `SequenceFileSet` shape has to describe either delivery, not because one set mixes them. Delivery format is a per-library choice declared on the sequencing information form — `UG - FASTQ` or `UG - CRAM` — so an order arrives as one or the other, and only the matching subcommand finds anything. Running the other one against the same prefix matches zero files and exits before writing, leaving any existing sheets untouched.

One set per (sample directory, read stem), so separate lanes become separate sets.

**Chunked reads are rejected.** A delivery holding `_R1_001` alongside `_R1_002` for one lane is refused during planning, before any per-file request: those are two pieces of one read, a SequenceFileSet holds a single file per slot, and splitting them across two sets would claim two sequencing runs happened where there was one. Neither Novogene nor Psomagen is expected to ship split FASTQs, so this means something upstream changed. The error names each affected slot and its chunks; concatenate to one file per read and re-run.

### Columns left blank on purpose

- **`library`** links to a Library object whose alias is not derivable from the S3 layout: a sample directory `…/GROUP_003_003/raw/` may correspond to library alias `…:library-GROUP-003`, a mapping only the submitter holds. Fill it by lookup — the diagnostic TSV carries `sample_dir` (the directory naming the CRO group) and `set_alias` (the key into the SequenceFileSet sheet) for exactly this.
- **`derived_from`**, **`description`**, **`submitter_comment`**, **`uuid`** are left for the submitter or the server.

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

### Scale h5ad mode

Point at a ScaleRna **rundate** directory (one level past `processed/`). Pairs rundate `samples.csv` barcodes to the Google Sheet `sample template` `RT_index` column, warns about control samples that have no pairing, and writes one TSV row per non-control `{rundate}/samples/*.h5ad`.

```bash
python -m file_extract scale_h5ad \
  s3://czi-cro/project/order/processed/run_date/ \
  --metadata-gid <google-sheet-uuid> \
  --cro-order NVUS0000000000-04 NVUS0000000000-05 \
  --wafers 426971 441969
```

| Flag | Description |
|------|-------------|
| `s3_uri` | **Required.** Rundate prefix, e.g. `s3://czi-cro/project/order/processed/run_date/` |
| `--metadata-gid` | **Required.** Google Sheet UUID (spreadsheet id in the URL). Reads tab `sample template` |
| `--cro-order` | **Required.** One or more CRO order identifiers |
| `--wafers` | **Required.** One or more wafer / RunIDs |
| `-o`, `--output` | Output TSV (default: `<run_date>_scale_h5ad_info.tsv`) |
| `--workers` | Thread count (default: min(64, n_files)) |
| `--retries` | Max attempts per transient S3 error (default: 5) |
| `--strict` | Exit 1 if any per-file CRC fetch fails |
| `-v`, `--verbose` | Debug logging |
| `-q`, `--quiet` | Disable progress bars |

**Pairing:** `RT_index` values such as `SCALEQUANT-A11` are stripped and flipped to `11A`. `samples.csv` `barcodes` such as `1A-2C` expand in column-wise 96-well order. A `samples.csv` row with no matching sheet well is a control: it is printed as a warning and its h5ad files are omitted.

**TSV columns:** `filename` · `s3_uri` · `crc64nvme_base64` · `sample` (first filename segment split on `.`)

`scale_cram` is not implemented yet; it can reuse the same `--metadata-gid`, `--cro-order`, and `--wafers` flags later.

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

The standalone prototypes in `bcp/docs/file_extractor.py` and `bcp/docs/extract_h5.py` are superseded by this package. Use `python -m file_extract fastq`, `python -m file_extract cram`, `python -m file_extract h5`, and `python -m file_extract scale_h5ad` instead.
