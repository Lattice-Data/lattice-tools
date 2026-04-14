## Mapping validation guide

This document explains how to run the `mapping_validation` checks from `bcp` and what to expect in the output, for providers **with** and **without** SIF files.

The validator is implemented in `mapping_validation.cli.main` and is exercised end‑to‑end in `test_mapping_validation_e2e.py` and `test_mapping_validation_cli.py`.

---

## Inputs and basic concepts

- **Mapping file**
  - **Format**: text file with **two columns per non‑empty line**:
    - column 1: S3 path
    - column 2: local filesystem path
  - **Provider-specific column order**:
    - Novogene exports are typically `S3 path,local path`.
    - Psomagen exports are sometimes `Local path,S3 path` and may include a header row
      like `Local Path,S3 Path`.
    - The parser uses `--provider` as a hint, and falls back to heuristics based on
      which column starts with `s3://`.
  - **Ignored lines**:
    - empty lines
    - lines starting with `#` (comments)
    - header/meta rows starting with `@` (e.g. `@NVUS...mapping_processed.csv`)
    - header rows like `Local Path,S3 Path` / `S3 Path,Local Path`
  - **Separators supported**:
    - comma: `s3://...,/ORPROJ1/...`
    - tab: `s3://...\t/ORPROJ1/...`
    - or `s3://...   /absolute/local/path` (space‑separated)
  - Empty lines and lines that cannot be split into two non‑empty fields are ignored.

- **SIF file (optional but strongly recommended for many modes)**
  - Provider “Sample Information Form” / intake sheet.
  - Used to derive:
    - GroupID → assay‑type expectations.
    - Library name → assay‑type expectations.
    - Library names when GroupIDs are not present (Ultima‑style).
  - Supported formats:
    - Excel: `.xlsx`, `.xlsm`, `.xls` (preferred and most realistic).
    - CSV: header‑based parsing with `DictReader`.
    - Ultima intake‑style CSV (Scale‑specific heuristic).

- **Providers and assay families**
  - `--provider`: `novogene` or `psomagen`.
  - `--assay`:
    - `10x` (10x Genomics)
    - `sci` (sci‑plex Ultima)
    - `scale` (Scale / Quantum Ultima)
  - `--data`:
    - `raw`: raw sequencing data.
    - `processed`: processed 10x outputs.

---

## CLI entry point

The package lives under `bcp/mapping_validation`. Run the module from the **`bcp` directory** (so Python resolves the package), or set `PYTHONPATH` to include `bcp`:

```bash
cd bcp
python -m mapping_validation --help
```

If your environment installs a console entry point named `mapping_validation`, you can use that instead:

```bash
mapping_validation --help
```

(The root `lattice-tools` `pyproject.toml` does not define that script; local runs typically use `python -m mapping_validation` from `bcp`.)

**Required arguments:**

- `--mapping PATH`  
  Path to the mapping CSV/TSV file (two columns: S3, local).

- `--provider {novogene,psomagen}`  
  Sequencing provider; controls bucket name, order patterns, and allowed assays.

- `--data {raw,processed}`  
  Kind of data the mapping represents.

- `--assay {10x,sci,scale}`  
  High‑level assay family:
  - `10x`: 10x Genomics
  - `sci`: Novogene sci Ultima
  - `scale`: Novogene Scale / Quantum Ultima

**Optional arguments:**

- `--sif PATH`  
  Path to SIF file (CSV or Excel). Enables SIF‑based completeness and consistency checks where supported.

---

## High‑level validation flow

Regardless of mode, the CLI runs the following phases in order:

1. **Uniqueness checks**
   - Implemented via `validate_uniqueness`.
   - Validates:
     - No duplicated S3 paths.
     - No duplicated local paths.
   - Output:
     - Summary line:  
       `Uniqueness: {total} mappings, {N} duplicate locals, {M} duplicate S3 paths`
   - **If any duplicates exist**, the exit code is forced to non‑zero and a failure reason is recorded.

2. **Mode‑specific validation**
   - Selected based on the `(provider, data, assay)` combination.
   - Currently implemented modes:
     - `provider in {novogene, psomagen}, data=raw, assay=10x`
     - `provider=novogene, data=raw, assay in {scale,sci}`  
       (**Scale and sci raw are Novogene-only**; `psomagen` with `scale` or `sci` is not implemented.)
     - `provider in {novogene, psomagen}, data=processed, assay=10x`
   - For any other combination, the tool prints:
     - `Mode (provider=..., data=..., assay=...) is not implemented yet.`
     - And fails with exit code `1`.

3. **Final verdict**
   - If all checks pass:
     - `VERDICT: PASS — mapping validates successfully against SOP rules.`
   - If any check fails:
     - `VERDICT: FAIL — mapping has issues that need attention:`
       - Followed by a bullet list of failure reasons (e.g. duplicate mappings, SIF completeness failures, S3/local inconsistencies).

Exit code semantics:

- `0`: all checks passed for the chosen mode.
- `1`: at least one check failed or mode unsupported.

---

## Supported modes and what they validate

In the command shapes below, run from **`bcp`** and use `python -m mapping_validation` (or a `mapping_validation` executable if installed).

### 10x raw data (Novogene and Psomagen)

Command shape:

```bash
python -m mapping_validation \
  --mapping PATH_TO_MAPPING.csv \
  --provider {novogene|psomagen} \
  --data raw \
  --assay 10x \
  [--sif PATH_TO_SIF]
```

Checks run:

- **S3 SOP validation (`validate_s3_10x_raw`)**
  - Ensures S3 paths follow provider‑specific 10x raw layout:
    - `s3://czi-{provider}/{project}/{order}/{GroupID}/raw/{RunID}-{GroupID}_{Assay}-{UG}-{Barcode}{suffix}`
  - Validates:
    - Bucket / provider (`czi-novogene` or `czi-psomagen`).
    - Order pattern (e.g. `NVUS...` for Novogene, `AN...` for Psomagen).
    - Assay name is allowed for that provider.
    - Barcodes consist only of `A/C/G/T`.
    - GroupID in the directory matches the filename prefix.
    - Project name casing / underscores (emits **warnings**, not hard errors).
  - Output:
    - Summary:  
      `10x raw SOP: matched {matched} S3 paths, {E} errors, {W} warnings[, {meta} run-metadata files]`
    - GroupID → assays and RunIDs summary.
    - Sampled error and warning lines.
  - Behavior:
    - If **no** S3 paths match the expected pattern:
      - Warns and records a failure reason.
    - If any SOP errors exist:
      - Marks the run as failed.

- **SIF‑based checks (only when `--sif` is provided)**

  When `--sif` is supplied, additional checks are run:

  - **GroupID vs assay completeness (`compare_groupid_assays`)**
    - Uses SIF to build expected GroupID → set of assays.
    - Normalizes SIF GroupIDs (e.g. `A + AF` → `A_AF`) and assay names to lower case.
    - Compares with actual S3 GroupID → assay mapping.
    - **Fails the run** when:
      - GroupIDs are present in SIF but missing from S3.
      - GroupIDs are present in S3 but absent from SIF.
      - A GroupID is in both, but S3 is **missing** assay types that the SIF expects.
    - **Reports only (does not fail)** when S3 has **extra** assay types for a GroupID that the SIF does not list (still printed under “unexpected assay types”).

  - **Library‑assay consistency (`validate_library_assay_consistency`)**
    - Runs only when the SIF can supply both **Library name** and **Assay type** columns (Excel or CSV); otherwise this step is skipped with no summary line.
    - Loads `Library name → assay type` from SIF.
    - Infers library names from local paths.
    - Validates:
      - Library name appears in the S3 GroupID.
      - S3 assay for that GroupID matches the assay expected in the SIF.

  - **Per‑path SIF coverage (`find_unmatched_sif_paths_10x`)**
    - Ensures that S3 paths map to SIF GroupIDs when possible.
    - Flags:
      - S3 paths whose GroupID is not present in SIF.
      - S3 paths that cannot be parsed (non‑metadata).
    - Reports:
      - Counts by GroupID.
      - Example paths per anomalous GroupID.

- **Without SIF (`--sif` omitted)**
  - Only uniqueness and S3 SOP checks are run.
  - The CLI prints:
    - `10x mode: no --sif provided, skipping SIF completeness checks.`
  - This is useful for a fast, structure‑only sanity check where SIF is not yet available.

---

### Novogene Scale raw data

Command shape:

```bash
python -m mapping_validation \
  --mapping PATH_TO_MAPPING.csv \
  --provider novogene \
  --data raw \
  --assay scale \
  [--sif PATH_TO_SIF]
```

Checks run:

- **S3 SOP validation (`validate_s3_seahub_raw` with family=`scale`)**
  - Supports both SOP and index forms of Scale S3 paths.
  - Validates:
    - Project name casing and formatting.
    - RunID consistency within the path.
    - GroupID structure.
    - Assay name is allowed and uses canonical spelling (e.g. `GEX`, `hash_oligo`, `GEX_hash_oligo`).
    - Common typos (e.g. `hash_oliga`) are treated as errors.
    - Family‑specific UG/Scaleplex logic (e.g. `hash_oligo` must have `SCALEPLEX` present in UG_RT; `GEX` must not).

- **Local path sanity (`validate_local_paths_scale_raw`)**
  - Validates Novogene Scale local layouts (two accepted versions).
  - Ensures internal consistency of:
    - Wafer IDs across directory and filename.
    - QSR numbers.
    - `SCALEPLEX` flags.
    - Index sequences.

- **SIF completeness (if `--sif` provided) (`validate_sif_completeness_seahub`)**
  - Loads expected GroupID → assays from Scale SIF.
  - Compares against GroupID → assays derived from S3 paths.
  - **Fails** on missing/extra GroupIDs and on **missing** assays per GroupID; **extra** assays in S3 vs SIF are printed but do not fail the run (same rules as 10x raw SIF comparison above).

- **S3/local fuzzy consistency (`validate_s3_local_consistency_scale`)**
  - Compares high‑level attributes extracted from S3 and local paths:
    - RunID vs wafer.
    - QSR numbers.
    - `SCALEPLEX` presence.
  - Emits errors when clearly inconsistent; warnings for partial mismatches.

- **Without SIF (`--sif` omitted)**
  - S3 SOP, local path checks, and S3/local consistency still run.
  - SIF completeness is **skipped**, with a message:
    - `scale mode: no --sif provided, skipping SIF completeness checks.`

---

### Novogene sci raw data

Command shape:

```bash
python -m mapping_validation \
  --mapping PATH_TO_MAPPING.csv \
  --provider novogene \
  --data raw \
  --assay sci \
  [--sif PATH_TO_SIF]
```

Checks run:

- **S3 SOP validation (`validate_s3_seahub_raw` with family=`sci`)**
  - Validates S3 path layout for sci raw data.
  - Ensures:
    - S3 prefix (`czi-novogene`) and order pattern.
    - GroupID and assay extraction.
    - UG and barcode formats (barcode limited to `A/C/G/T`).
    - Project naming conventions.

- **Local path sanity (`validate_local_paths_sci_raw`)**
  - Validates internal consistency of sci local directory structures:
    - RunID agreement across directory segments and filename.
    - GroupID consistency.
    - UG and barcode consistency across directory and filename.

- **SIF completeness (if `--sif` provided) (`validate_sif_completeness_seahub`)**
  - Same pattern as Scale:
    - SIF GroupID → assays vs S3 GroupID → assays.
    - **Fails** on missing/extra GroupIDs and **missing** assays; **extra** S3 assays vs SIF are informational only.

- **S3/local consistency (`validate_s3_local_consistency_sci`)**
  - Cross‑checks S3 vs local for:
    - RunID.
    - GroupID.
    - UG.
    - Barcode.
  - Reports and counts mismatches.

- **Without SIF (`--sif` omitted)**
  - S3 SOP, local path, and S3/local consistency still run.
  - SIF completeness is **skipped**, with:
    - `sci mode: no --sif provided, skipping SIF completeness checks.`

---

### 10x processed data (Novogene and Psomagen)

Command shape:

```bash
python -m mapping_validation \
  --mapping PATH_TO_MAPPING.csv \
  --provider {novogene|psomagen} \
  --data processed \
  --assay 10x \
  [--sif PATH_TO_SIF]
```

Checks run:

- **S3 SOP validation (`validate_s3_10x_processed`)**
  - Expected layout:
    - `s3://czi-{provider}/{project}/{order}/{GroupID}/processed/{pipeline}/{Run_YYYY-MM-DD}/outs/{file_path}`
  - Validates:
    - Pipeline name (currently `cellranger` only).
    - Run‑date format (`Run_YYYY-MM-DD`).
    - Project naming conventions.
  - Summarizes:
    - Unique GroupIDs.
    - Pipelines.
    - Run dates.

- **SIF completeness (if `--sif` provided) (`validate_sif_completeness_10x_processed`)**
  - Uses SIF to derive expected identifiers:
    - Prefer rows keyed by **Group Identifier** (plus **Assay type**) when those columns exist.
    - If that structure is not present, fall back to unique **Library name** values (Ultima-style SIFs).
  - Compares SIF identifiers vs GroupIDs present in processed S3.
  - **Fails** when identifiers are in the SIF but missing from the mapping’s processed S3 paths.
  - **Prints** GroupIDs that appear in S3 but not in the SIF (`Extra in S3`); that output is **informational** and does not by itself set a failing exit code.

- **S3/local consistency (`validate_s3_local_consistency_10x_processed`)**
  - Ensures:
    - GroupID from the S3 path appears in the local path.
    - Relative path after `/outs/` matches between S3 and local.
  - Emits errors when mismatches are found; warnings when `/outs/` is not present in the local path.
  - When the GroupID string differs only by `-` vs `_`, the CLI may still match and emit a **NOTE** counting warnings matched after normalization.

- **Without SIF (`--sif` omitted)**
  - S3 SOP and S3/local consistency checks still run.
  - SIF completeness is skipped, with:
    - `10x processed mode: no --sif provided, skipping SIF completeness checks.`

---

## Interpreting output

The CLI prints:

- **Summary lines** for each validator:
  - Uniqueness.
  - S3 SOP (per mode).
  - Local path checks.
  - SIF completeness.
  - S3/local consistency.

- **Grouped counts of errors and warnings**:
  - Each labeled section shows:
    - `{label} errors:` / `{label} warnings:`
    - `{type}: {count}` by error/warning type.

- **Sampled examples**
  - For each section with issues, up to a fixed number of examples are printed.
  - Each example includes:
    - Line number in the mapping.
    - S3 path and/or local path.
    - Short `detail` description.

- **Final verdict**
  - Always at the end, making it easy to grep for `VERDICT`.

**Informational-only lines:** Some SIF-related messages (for example **extra** assay types in S3 vs SIF for a GroupID, and **extra** processed GroupIDs in S3 vs SIF) are printed for review but **do not** add a failure reason or change `VERDICT` / exit code. Treat `VERDICT: FAIL` and the bulleted failure reasons as the source of truth for whether the run failed.

Typical workflow:

1. From `bcp`, run `python -m mapping_validation` (or your installed `mapping_validation` command) with appropriate `--provider`, `--data`, `--assay`, and `--sif` whenever a SIF exists.
2. If it fails:
   - Start with the **summary counts** to see which category is failing (SOP vs SIF completeness vs S3/local vs uniqueness).
   - Inspect the **sample issues** for concrete examples.
3. Fix mapping/S3/local path issues or SIF contents as appropriate.
4. Re‑run until `VERDICT: PASS` and exit code `0`.

---

## Example commands

Run these from the **`bcp`** directory (paths are illustrative). Replace `python -m mapping_validation` with `mapping_validation` if your environment provides that executable.

- **Novogene 10x raw with SIF**

```bash
python -m mapping_validation \
  --mapping novogene_10x_raw.csv \
  --sif novogene_10x_sif.xlsx \
  --provider novogene \
  --data raw \
  --assay 10x
```

- **Novogene 10x raw without SIF (structure‑only check)**

```bash
python -m mapping_validation \
  --mapping novogene_10x_raw.csv \
  --provider novogene \
  --data raw \
  --assay 10x
```

- **Psomagen 10x raw with SIF**

```bash
python -m mapping_validation \
  --mapping psomagen_10x_raw.csv \
  --sif psomagen_10x_sif.csv \
  --provider psomagen \
  --data raw \
  --assay 10x
```

- **Novogene sci raw with SIF**

```bash
python -m mapping_validation \
  --mapping novogene_sci_raw.csv \
  --sif novogene_sci_sif.xlsx \
  --provider novogene \
  --data raw \
  --assay sci
```

- **Novogene Scale raw with SIF**

```bash
python -m mapping_validation \
  --mapping novogene_scale_raw.csv \
  --sif novogene_scale_sif.xlsx \
  --provider novogene \
  --data raw \
  --assay scale
```

- **10x processed with SIF (Novogene)**

```bash
python -m mapping_validation \
  --mapping novogene_10x_processed.csv \
  --sif novogene_10x_processed_sif.xlsx \
  --provider novogene \
  --data processed \
  --assay 10x
```

- **10x processed with SIF (Psomagen)**

```bash
python -m mapping_validation \
  --mapping psomagen_10x_processed.csv \
  --sif psomagen_10x_processed_sif.csv \
  --provider psomagen \
  --data processed \
  --assay 10x
```

Use the e2e fixtures under `bcp/tests/fixtures/mapping_validation` as concrete examples of mappings and SIFs that are expected to **pass** and **fail** for each mode.

