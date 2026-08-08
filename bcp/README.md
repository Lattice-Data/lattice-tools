# CRISPR Guide Processing Pipeline

A Python pipeline for processing CRISPR guide sequences using guidescan2 to find genomic matches and calculate coordinates.

## Overview

This pipeline takes CRISPR guide sequences as input, uses guidescan2 to find genomic matches, retains all exact matches for each guide, and outputs formatted genomic coordinates.

**Key Features:**
- Validates PAM uniformity across all guides
- Supports variable PAM lengths (NGG, NNGG, etc.)
- Finds exact genomic matches using guidescan2
- **Retains ALL exact matches** - if a guide has multiple exact genomic matches, all are retained
- Only exact matches (distance=0) are reported; non-exact matches result in NA
- Calculates genomic coordinates for both + and - strands
- Handles guides with no genomic matches
- Comprehensive error handling and logging

## Requirements

### Software Dependencies

1. **Python 3.8+** with conda environment
2. **guidescan2** - Must be installed and available in PATH
   - Installation: See [guidescan2 documentation](https://github.com/pritykinlab/guidescan-cli)
3. **Python packages**:
   - pandas
   - pathlib (standard library)

### Installation

```bash
# Activate your conda environment
conda activate lattice  # or your environment name

# Verify guidescan is installed
guidescan --help

# Navigate to the bcp directory
cd bcp/
```

## Input Format

The pipeline expects a CSV file with the following columns:

```csv
id,sequence,pam,chromosome,position,sense
guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,
guide2,GAGTTCGCTGCGCGCTGTT,NGG,,,
guide3,GCCCGCTCCCCGCGATCCC,NGG,,,
```

**Column descriptions:**
- `id`: Unique identifier for each guide
- `sequence`: Protospacer sequence (without PAM)
- `pam`: PAM sequence (e.g., NGG, NNGG)
- `chromosome`, `position`, `sense`: Initially empty, filled by pipeline

**Important:** All guides in the file must have identical PAM sequences.

## Usage

### Basic Usage

```bash
python guidescan_pipeline.py \
  -i /path/to/genome.fa.index \
  -g /path/to/guides.txt \
  -o /path/to/output_dir/
```

### Advanced Options

```bash
# Keep intermediate files for inspection
python guidescan_pipeline.py \
  -i genome.fa.index \
  -g guides.txt \
  -o results/ \
  --keep-intermediate

# Enable verbose logging
python guidescan_pipeline.py \
  -i genome.fa.index \
  -g guides.txt \
  -o results/ \
  --verbose

# View help
python guidescan_pipeline.py --help
```

### Command-Line Arguments

**Required:**
- `-i, --genome-index`: Path to genome index file (`.fa.index`)
- `-g, --guides`: Path to input guides file (CSV format)
- `-o, --output-dir`: Directory for output files

**Optional:**
- `--keep-intermediate`: Keep intermediate files (default: delete after completion)
- `-v, --verbose`: Enable verbose output for debugging
- `--version`: Show version information

## Output Files

### Final Output: `exact_matches_formatted.csv`

Contains genomic coordinates for exact matches. **If a guide has multiple exact matches, multiple rows are output:**

```csv
id,sequence,pam,chromosome,start,end,sense
guide1,TGCCTCGCGCAGCTCGCGG,NGG,1,1000,1018,+
guide1,TGCCTCGCGCAGCTCGCGG,NGG,5,2340,2358,+
guide2,GAGTTCGCTGCGCGCTGTT,NGG,3,1982,2000,-
guide3,GCCCGCTCCCCGCGATCCC,NGG,NA,NA,NA,NA
```

**Column descriptions:**
- `id`: Guide identifier (can appear multiple times if multiple exact matches exist)
- `sequence`: Protospacer sequence
- `pam`: PAM sequence
- `chromosome`: Chromosome where exact match found (NA if no exact match)
- `start`, `end`: Genomic coordinates (1-based, inclusive)
- `sense`: Strand orientation (+ or -, NA if no exact match)

**Note:** Only exact matches (distance=0) are included. Guides with no exact matches will have a single row with NA values.

### Intermediate Files (if `--keep-intermediate` is used)

1. **`all_matches.csv`**: All genomic matches from guidescan (includes exact and non-exact matches)
2. **`best_matches.csv`**: All exact matches (distance=0) for each guide; non-exact matches result in NA

## Example Workflow

### 1. Prepare Input File

Create `guides.txt`:
```csv
id,sequence,pam,chromosome,position,sense
BRCA1-guide1,GCACTGATCTAGCTAGCTAG,NGG,,,
BRCA1-guide2,ATCGATCGATCGATCGATCG,NGG,,,
TP53-guide1,GGCCAATTCCGGAATTCCGG,NGG,,,
```

### 2. Run Pipeline

```bash
python guidescan_pipeline.py \
  -i /data/genomes/hg38.fa.index \
  -g guides.txt \
  -o results/ \
  --verbose
```

### 3. Check Output

```bash
# View results
cat results/exact_matches_formatted.csv

# Count exact matches
grep -v "NA" results/exact_matches_formatted.csv | wc -l
```

### 4. Example Output

```
2025-01-15 10:30:00 - INFO - Starting CRISPR guide processing pipeline...
2025-01-15 10:30:00 - INFO - Validating PAM sequences from input guides...
2025-01-15 10:30:00 - INFO - PAM sequence validated: NGG (length: 3)
2025-01-15 10:30:00 - INFO - Running guidescan enumerate...
2025-01-15 10:30:15 - INFO - guidescan completed. Results saved to results/all_matches.csv
2025-01-15 10:30:15 - INFO - Selecting matches for each guide (retaining all exact matches only)...
2025-01-15 10:30:15 - INFO - Matches saved to results/best_matches.csv
2025-01-15 10:30:15 - INFO - Processed 3 guides: 5 exact matches, 0 no matches
2025-01-15 10:30:15 - INFO - Formatting exact matches...
2025-01-15 10:30:15 - INFO - Processing complete. Output saved to: results/exact_matches_formatted.csv
2025-01-15 10:30:15 - INFO - Pipeline completed successfully!
2025-01-15 10:30:15 - INFO - Summary: 3 guides processed, 5 exact matches found
```

## Using as a Python Module

```python
from guidescan_pipeline import GuidescanPipeline

# Initialize pipeline
pipeline = GuidescanPipeline(
    genome_index="genome.fa.index",
    guides_file="guides.txt",
    output_dir="results/",
    keep_intermediate=True,
    verbose=True,
)

# Run complete pipeline
results_df = pipeline.run()

# Access results
print(f"Total output rows: {len(results_df)}")
print(f"Rows with exact matches: {len(results_df[results_df['chromosome'] != 'NA'])}")
print(f"Unique guides processed: {results_df['id'].nunique()}")

# Results are also saved to results/exact_matches_formatted.csv
```

## Testing

The pipeline includes comprehensive tests covering:
- PAM validation and length calculation
- Best match selection logic
- Coordinate calculation for both strands
- Integration tests with mocked guidescan
- **End-to-end tests with real guidescan** (opt-in)

### Running Tests

```bash
# Run regular tests (default, fast)
cd bcp/
conda activate lattice
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_pam_length.py -v

# Run E2E tests with real guidescan (opt-in)
pytest -m e2e -v

# Run ALL tests (regular + E2E)
pytest -m "" -v

# Run with coverage (if pytest-cov installed)
pytest --cov=. --cov-report=html
```

### Test Data Sanitization

Test fixtures and inline test literals in `bcp/tests` must use anonymized identifiers only.

- Use synthetic bucket/provider names such as `vendor-novogene` and `vendor-psomagen`.
- Use synthetic project/lab slugs such as `project-*` or `lab-*`; avoid real surnames or institution names.
- Use format-preserving sanitized order IDs:
  - `NVUS0000000000-<suffix>`
  - `AN0000000<digit>`
- Use sanitized local paths like `/local/user_001/...`; do not use real usernames or internal mount roots.
- Keep path structure and file naming shape intact so parser/validator behavior remains realistic.

An anti-regression test (`bcp/tests/test_sanitized_identifiers.py`) enforces these rules and fails when blocked patterns are reintroduced.

**Test Statistics**: 59 total tests (51 regular + 8 E2E)

### Test Organization

```
bcp/tests/
├── conftest.py                 # Shared fixtures
├── fixtures/                   # Test data files
│   ├── test_guides.txt        # Sample guides (5 guides)
│   ├── toy.fa.index.*         # Toy genome index for E2E tests
│   ├── genes.gtf              # Gene annotations for testing
│   ├── all_matches_sample.csv
│   ├── best_matches_sample.csv
│   └── expected_all_matches.csv
├── test_pam_length.py         # PAM validation tests (7 tests)
├── test_best_match.py         # Match selection tests (8 tests)
├── test_formatting.py         # Coordinate calculation tests (10 tests)
├── test_integration.py        # Integration tests (11 tests)
├── test_guide_to_gene.py      # Gene annotation tests (15 tests)
└── test_e2e_real_guidescan.py # E2E tests with real guidescan (8 tests)
```

For detailed testing documentation, see:
- **[TESTING.md](TESTING.md)** - Complete testing guide
- **[E2E_TESTS.md](E2E_TESTS.md)** - End-to-end testing documentation

## Troubleshooting

### Common Errors

**Error: "guidescan not found in PATH"**
```
Solution: Install guidescan2 and ensure it's in your PATH
```

**Error: "All PAM sequences must be identical"**
```
Solution: Check your input file - all guides must use the same PAM (e.g., all NGG)
```

**Error: "No PAM sequences found in input file"**
```
Solution: Ensure your CSV has a 'pam' column and it's not empty
```

**Error: "Genome index not found"**
```
Solution: Check the path to your genome index file exists
```

### Getting Help

```bash
# View command-line help
python guidescan_pipeline.py --help

# Check guidescan version
guidescan --version

# Run tests to verify installation
pytest bcp/tests/
```

## Pipeline Logic

### Workflow Steps

1. **Validation**
   - Check input files exist
   - Validate all guides have identical PAM
   - Calculate PAM length

2. **Guidescan Execution**
   - Run guidescan enumerate with parameters:
     - Max mismatches: 1
     - Max off-targets: 1
     - Mode: complete

3. **Match Selection**
   - Group matches by guide ID
   - **Keep ALL exact matches** (distance = 0) for each guide
   - Non-exact matches (distance > 0) result in NA values
   - Treats guides without exact matches the same as guides with no genomic matches

4. **Coordinate Formatting**
   - Process all exact matches from selection step
   - Split sequence into protospacer + PAM
   - Calculate genomic coordinates:
     - **+ strand**: start = position, end = position + length - 1
     - **- strand**: end = position, start = position - length + 1
   - Output multiple rows if guide has multiple exact matches

5. **Output Generation**
   - Save formatted results to CSV
   - Optionally clean up intermediate files

### Coordinate Calculation

For a 19bp protospacer at position 1000:

**+ strand:**
```
Position: 1000
Start:    1000
End:      1018  (1000 + 19 - 1)
```

**- strand:**
```
Position: 1000
End:      1000
Start:    982   (1000 - 19 + 1)
```

## Gene Annotation with guide_to_gene.py

After running the main pipeline, you can annotate guides with gene identifiers using the `guide_to_gene.py` script. This tool maps guide coordinates to genes using a GTF annotation file.

### Usage

```bash
python guide_to_gene.py \
  results/exact_matches_formatted.csv \
  /path/to/genes.gtf \
  results/guides_with_genes.csv
```

### Input Requirements

1. **Guide file**: Output from guidescan_pipeline.py (`exact_matches_formatted.csv`)
2. **GTF file**: Gene annotation file (e.g., GENCODE, Ensembl)
3. **Output file**: Path for annotated results

### Output Format

The script creates **one row per guide-gene overlap**. If a guide overlaps multiple genes, multiple rows are created:

```csv
id,sequence,pam,chromosome,start,end,sense,gene_id,gene_name
guide1,TGCCTCGCGCAGCTCGCGG,NGG,1,1000,1018,+,ENSG00000123456,BRCA1
guide1,TGCCTCGCGCAGCTCGCGG,NGG,5,2340,2358,+,ENSG00000234567,TP53
guide2,GAGTTCGCTGCGCGCTGTT,NGG,3,1982,2000,-,ENSG00000345678,EGFR
guide3,GCCCGCTCCCCGCGATCCC,NGG,10,5000,5018,+,intergenic,intergenic
guide4,ATCGATCGATCGATCGATC,NGG,NA,NA,NA,NA,NA,NA
```

### Annotation Categories

- **Gene overlap**: Includes `gene_id` (Ensembl ID) and `gene_name` (symbol)
- **Intergenic**: Guides between genes get `intergenic` for both fields
- **Unmapped**: Guides with NA coordinates (no genomic match) get `NA` for both fields
- **Invalid coordinates**: Malformed coordinates get `INVALID_COORDS`

### Complete Workflow Example

```bash
# Step 1: Find genomic matches for guides
python guidescan_pipeline.py \
  -i /data/genomes/hg38.fa.index \
  -g my_guides.csv \
  -o results/ \
  --verbose

# Step 2: Annotate with gene information
python guide_to_gene.py \
  results/exact_matches_formatted.csv \
  /data/annotations/gencode.v44.annotation.gtf \
  results/guides_annotated.csv

# Step 3: View results
head results/guides_annotated.csv
```

### Statistics Output

The script provides summary statistics:

```
Parsing GTF file: gencode.v44.annotation.gtf
Found genes on 25 chromosomes

Annotating guides from: results/exact_matches_formatted.csv

Annotation complete!
  Total guides processed: 100
  Guides with gene overlaps: 85
  Guides in intergenic regions: 10
  Guides with unmapped/invalid coordinates: 5
  Total gene overlaps found: 92
  Average overlaps per guide (for guides with overlaps): 1.08

Output written to: results/guides_annotated.csv
```

### Features

- **Chromosome normalization**: Automatically handles chromosome naming (e.g., "1" vs "chr1")
- **Multiple overlaps**: If a guide overlaps 3 genes, creates 3 output rows
- **Edge case handling**: Properly handles NA coordinates, invalid values, and intergenic regions
- **GTF compatibility**: Works with GENCODE, Ensembl, and RefSeq GTF files

### Command-Line Help

```bash
python guide_to_gene.py --help
```

## Gene Annotation with Bioframe (Alternative Implementation)

An alternative implementation using the `bioframe` library is available in `guide_to_gene_bioframe.py`. This provides the same functionality but uses bioframe's optimized overlap detection algorithms.

### Overview

The bioframe implementation leverages the `bioframe` library for genomic interval operations, providing:
- Standardized genomic data structures (bedframes)
- Optimized interval tree-based overlap detection
- Built-in GTF/BED file parsing
- Extensibility for additional genomic operations

**Status**: ✅ Production-ready with full test coverage (25/25 tests passing)

### Installation

```bash
# Install bioframe using conda (recommended)
conda install -c bioconda bioframe

# Or using pip
pip install bioframe>=0.4.0

# Install all dependencies including bioframe
cd bcp/
pip install -r requirements.txt
```

### Usage

**Basic command-line usage:**

```bash
python guide_to_gene_bioframe.py \
  results/exact_matches_formatted.csv \
  /path/to/genes.gtf \
  results/guides_with_genes.csv
```

**Python module usage:**

```python
from guide_to_gene_bioframe import annotate_guides_bioframe

annotate_guides_bioframe(
    guide_file="results/exact_matches_formatted.csv",
    gtf_file="/path/to/genes.gtf",
    output_file="results/guides_annotated.csv",
)
```

### Key Features

- ✅ **Same output format**: Produces identical output to `guide_to_gene.py`
- ✅ **Same statistics**: Provides identical statistics and error handling
- ✅ **Optimized performance**: Uses bioframe's interval tree algorithms
- ✅ **Chromosome normalization**: Automatically handles "1" vs "chr1" formats
- ✅ **Edge case handling**: NA coordinates, invalid values, intergenic regions
- ✅ **Multiple overlaps**: Creates one row per guide-gene overlap
- ✅ **Fully compatible**: Drop-in replacement for `guide_to_gene.py`

### Implementation Details

**Core functions:**
- `csv_to_bed_dataframe()`: Converts guide CSV to BED format with chromosome normalization
- `load_genes_gtf()`: Loads GTF and extracts gene_id/gene_name from attributes
- `find_overlaps_bioframe()`: Uses `bioframe.overlap()` for efficient overlap detection
- `format_output()`: Formats results to match current implementation
- `annotate_guides_bioframe()`: Main function with same interface as `guide_to_gene.py`

**Coordinate handling:**
- Converts 1-based inclusive coordinates to 0-based start, 1-based end (BED format)
- Ensures int64 data types for bioframe compatibility
- Maintains accuracy across coordinate systems

**GTF attribute parsing:**
- Extracts `gene_id` and `gene_name` using regex from attributes column
- Handles missing attributes gracefully (defaults to 'NA')
- Compatible with GENCODE, Ensembl, and RefSeq GTF formats

### Testing

The bioframe implementation includes comprehensive test coverage:

**Test suite: 25 tests across 2 test files**

```bash
# Run bioframe implementation tests (17 tests)
pytest tests/test_guide_to_gene_bioframe.py -v

# Run compatibility tests comparing implementations (8 tests)
pytest tests/test_compatibility_bioframe_vs_current.py -v

# Run all bioframe tests
pytest tests/test_guide_to_gene_bioframe.py tests/test_compatibility_bioframe_vs_current.py -v
```

**Test coverage:**
- CSV to BED format conversion (3 tests)
- GTF loading and parsing (1 test)
- Bioframe overlap detection (3 tests)
- Output formatting (2 tests)
- End-to-end annotation (4 tests)
- Edge cases: NA coordinates, invalid values (2 tests)
- Compatibility with current implementation (8 tests)
- Multiple overlaps handling (2 tests)

### Compatibility Verification

The bioframe implementation has been verified to produce **identical results** to `guide_to_gene.py`:

| Aspect | Verified |
|--------|----------|
| Output column names and order | ✅ |
| Number of output rows | ✅ |
| Guide coverage | ✅ |
| Gene annotations (gene_id, gene_name) | ✅ |
| Statistics output format | ✅ |
| Edge case handling (NA, invalid, intergenic) | ✅ |
| Multiple overlap behavior | ✅ |

### When to Use

**Use `guide_to_gene.py` (current implementation) if:**
- You prefer minimal dependencies
- You want a simple, straightforward implementation
- You don't need additional genomic operations

**Use `guide_to_gene_bioframe.py` (bioframe implementation) if:**
- You already use bioframe in your workflow
- You prefer standardized genomic interval libraries
- You want optimized overlap algorithms for large datasets
- You plan to extend with additional bioframe operations

Both implementations produce identical results and handle edge cases identically.

### Example Output

Both implementations produce the same output format:

```csv
id,sequence,pam,chromosome,start,end,sense,gene_id,gene_name
FBXW2-1,GGCTGCGGACCGGGAGCAG,NGG,9,120793348,120793366,-,ENSG00000119402,FBXW2
FBXW2-1,GGCTGCGGACCGGGAGCAG,NGG,9,120793348,120793366,-,ENSG00000214654,B3GALT9
ARMC5-1,TGCCTCGCGCAGCTCGCGG,NGG,16,31459567,31459585,+,ENSG00000260267,ENSG00000260267
ARMC5-1,TGCCTCGCGCAGCTCGCGG,NGG,16,31459567,31459585,+,ENSG00000140691,ARMC5
```

### Future Enhancements

Potential improvements using bioframe:
1. **Strand-specific overlaps**: Filter overlaps by strand orientation
2. **Partial overlap metrics**: Calculate overlap percentages
3. **Distance to nearest gene**: For intergenic guides
4. **Multiple GTF sources**: Merge annotations from multiple files
5. **Additional bedtools-style operations**: Leverage bioframe's full feature set

## Protospacer set signature (`guidesig`)

`guidesig` computes a deterministic string for the **set** of protospacer sequences in a
guide-template TSV or CSV. Two files with the same protospacer set (ignoring row order and
all non-protospacer columns) produce the same signature; otherwise they differ.

### What it guarantees

- Same set of protospacer strings → identical `gsig1:set:n=…:…` signature
- Different set → different signature
- For structurally valid rows, layout, column order, case, padding, duplicate rows, and
  empty protospacer cells do not affect the result
- The input format does not affect the result: the same sequences stored as TSV and as CSV
  produce the same signature

### Input format is mandatory

`--format tsv` or `--format csv` is **required**, and the library equivalent
`file_format=` is a required keyword argument. The format is never inferred from the file
extension and never sniffed from the contents, because content sniffing misclassifies the
sparse single-column rows these templates contain. Naming the wrong format fails loudly
(the protospacer column will be reported as absent from the header) rather than
misparsing.

`--compare` applies one `--format` to both files, so a mixed CSV/TSV pair is rejected,
naming the offending file and the format it looks like. Signatures themselves are
format-independent, so comparing across formats means signaturing each file under its own
`--format` and reading the two lines.

### Input validation

Input must be a UTF-8 (optionally BOM-prefixed) TSV or CSV with a header containing
`guide_protospacer`, or the name supplied through `--column`. In CSV input, quoted cells
containing commas are parsed correctly and do not shift the protospacer column. The tool:

- ignores fully blank rows and rows whose first cell starts with `#` after leading
  whitespace is removed; the comment marker is always read from the first column,
  independent of `--column`, so a `#` in any other cell is data
- strips and uppercases protospacers, skipping empty values
- rejects non-empty rows shorter than the header, including rows where trailing tabs were
  dropped during export
- rejects rows carrying values in fields beyond the header, which indicate a truncated
  header or a shifted export; a surplus field that is empty or whitespace-only is
  tolerated, because a trailing delimiter is a routine export artifact
- rejects a missing signature column, non-ACGT protospacer values, invalid UTF-8, and
  files with no usable protospacers

Row numbers in rejection messages are **physical file line numbers**, counted the way a
text editor counts them. A quoted cell containing an embedded newline spans several lines,
and such a record is reported at the line it starts on, so the number always points at
what the user sees when they open the file.

### What it deliberately does **not** guarantee

This is a **string-identity** check on one column, not a biological equivalence check.

- Sequence length is not normalized
- A leading/trailing `G` is not stripped
- Reverse-complement is not applied
- Strand / orientation is not canonicalized
- Sequences are not mapped to genomic coordinates

Two files that differ only by a trimmed base or a strand flip **must** produce different
signatures. A false match is the failure mode this tool exists to prevent.

### Version prefix

Signatures are prefixed with `gsig1` so a future change to the rules cannot silently make
old and new signatures comparable. The `n=` field reports unique sequence count for
mismatch diagnosis without re-parsing the files.

### Digest width

The digest is SHA-256 truncated to 32 hex characters, or 128 bits.

The risk being managed is **accidental** divergence between two curated libraries, not an
adversarial collision. This is not a tamper-evidence mechanism: anyone who can edit a
library can recompute its signature, so a matching signature attests that two files hold
the same protospacer set, not that either is authentic. Under that threat model 128 bits
is far more than enough — the birthday bound sits at roughly 2^64 libraries, so the
probability of two distinct guide sets colliding stays negligible at any library count
this project will ever hold.

Truncation buys legibility, which matters because these signatures get pasted into
tickets, spreadsheet cells, and commit messages: 32 hex characters stay one token a person
can compare by eye, where 64 invite a truncated copy-paste that silently drops the half
that differs. If the threat model ever changes — signatures crossing a trust boundary, say
— the `gsig1` prefix is the escape hatch, and `gsig2` can carry the full digest without
either version being mistaken for the other.

### Usage

Stdlib only (`csv`, `hashlib`, `argparse`, `pathlib`). From the `bcp/` directory:

```bash
# One signature line per file: <signature>\t<path>
python -m guidesig --format tsv path/to/lib.tsv

# CSV input
python -m guidesig --format csv path/to/lib.csv

# Compare two libraries (both must be in the given format)
python -m guidesig --format tsv --compare lib_a.tsv lib_b.tsv

# Alternate column
python -m guidesig --format tsv lib.tsv --column guide_rc_sequence
```

Output is always `<signature>\t<path>`, tab-separated regardless of the input format.

Library API:

```python
from guidesig import signature, protospacer_set

sig = signature("lib.tsv", file_format="tsv")
seqs = protospacer_set("lib.csv", file_format="csv")
```

Regression fixtures and golden digests live in `tests/fixtures/lib_*.tsv` and
`tests/fixtures/lib_two_layout_a.csv`, and are covered by `tests/test_guidesig.py` and
`tests/test_guidesig_csv.py`. The CSV fixture is a conversion of the layout A TSV and is
pinned to the same digest, which is what proves the two readers agree.

## Contributing

When modifying the pipeline:

1. Write tests for new functionality
2. Run the test suite: `pytest`
3. Update this README if adding features
4. Follow existing code style (use Ruff for linting)

### QA validation note (S3 gather)

For `bcp/qa.ipynb` and `qa_gather` raw-file validation, 10x group comparisons now normalize group IDs before checking equality:

- `-` and `_` are treated as equivalent (for example, `CUIMC-001` equals `CUIMC_001`)
- comparison is case-insensitive
- true mismatches still produce `WRONG GROUP` errors

This avoids false-positive group errors when providers use different separators between folder names and file names.

### SeaHub lab raw QA (`raw_assay='seahub_sci'`)

Collaborator lab trimmed uploads (`czi-trapnell` / `czi-hamazaki`, `*-seahub-bcp`) are QA'd in
`bcp/qa.ipynb` with `raw_assay='seahub_sci'`; processed validation is off for this mode. The SOP is:

```
s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/{sublibrary}/{wafer}/
    {wafer}-{sublibrary}[_{well}]_{sublibrary type}-{UG}-{barcode}.trim.*
```

with `{sublibrary type}` in `GEX`, `CRI`, `hash_oligo`, `GEX_hash_oligo`. The ExperimentID is not a
filename field; it appears in a name only when it is part of that project's sublibrary name (e.g.
`REF3_P05_2`).

Both `data_source` modes are supported, and they are the same question asked twice, so they must
return the same objects: the s3 walk lists everything under `{ExperimentID}/raw/` recursively,
matching manifest mode's `"/raw/" in key`, and uses a delimiter pass only to enumerate wafer
*folders*. Reading just the folder tree and keeping the wafer-level objects, as it once did,
dropped every object above wafer depth — so `bad_path_depth` could fire only for keys that were too
deep, never too shallow, and a `raw/` holding only loose objects was reported as missing. In
**s3** mode the ExperimentID is the last segment of the
listing prefix. In **manifest** mode it comes from the `order` argument if one is given, and is
otherwise read off the manifest keys, which already contain it as a folder; a manifest mixing two
ExperimentIDs is rejected at `resolve_qa_run_context` rather than per object, because QA writes one
`{order}_*` output set and filters the vendor index to a single ExperimentID. `ctx.order` must never
be empty for a SeaHub run: it is both the expected value of the cross-experiment check and the guard
that enables the vendor index's experiment filter, so an empty one silently turns the first into one
error per object and the second off entirely.

Five checks are specific to this mode. All SeaHub cells form one contiguous block at the *end* of
`qa.ipynb`, tagged `seahub` and gated on `seahub_active`, so they cannot renumber the shared pipeline;
`tests/test_qa_notebook_layout.py` enforces that. Read `{order}_seahub_well_status.csv` first: the SOP
table says what is wrong, the per-well status says how many wells are affected and which need data
rather than a rename.

- **Two suffix families.** `qa_constants.SEAHUB_TRIM_SUFFIXES` is the SOP set; `SEAHUB_BARE_SUFFIXES`
  is the same six artifacts with the `.trim` infix dropped, which real uploads have used. Completeness
  asks only whether each of the five artifact *kinds* arrived, under either spelling — naming is a
  separate axis, reported once per stem as `missing_trim_infix`. Judging by family instead let one
  optional sidecar carrying the other family's name decide the whole requirement set, and the two
  completeness paths reacted to that differently: a well with five correct `.trim.*` files plus a
  stray bare sidecar was complete to `roll_up_wells` and missing five files to
  `check_expected_raw_files`. A genuinely absent artifact is still reported, under its SOP name
  whichever spelling the upload used.
- **SOP validation** (`qa_seahub_sop.py`) reports each broken rule once per distinct fact, at five
  scopes: `object` (per object), `stem` (per well), `folder` (per sublibrary directory, so one
  truncated folder is one row rather than one per well), `suffix` (per distinct unrecognised
  extension) and `upload` (per distinct bucket/project fact, so a wrong bucket is one row rather
  than one per well — on a 288-well upload that is what keeps every other finding visible).
  A lab that spells its artifacts something else entirely misspells *every* object, so
  `unexpected_suffix` collapses to one row per extension naming the count, an example and the SOP
  artifacts expected; if nothing at all is recognised, `no_recognized_artifacts` fires once and the
  notebook prints a FAIL banner, since no well can be identified and every table below would
  otherwise be empty and read as clean. A bare suffix is only accepted when the stem before it
  parses — `.csv` is generic enough that `<well>.trimmer_stats.csv` otherwise becomes its own well,
  which turned one real 480-well upload into 960. Rules cover wrong bucket or project, wrong
  path depth, missing `.trim` infix, a doubled leading wafer token (`duplicated_wafer_token`), a
  truncated sublibrary folder (`sublibrary_folder_truncated`), bulk-download leftovers
  (`non_sequencing_artifact`), sublibrary or wafer disagreement, a bad well token, and a missing
  sublibrary type. `expected_name` carries the corrected basename repairing that rule alone, and
  `expected_folder` the corrected path segment. A folder/filename disagreement is blamed on the
  **folder** when the filename says `{ExperimentID}_{folder}`, because the vendor delivery is the
  source of truth for the sublibrary name. Given `untrimmed_s3_paths`, a missing sublibrary type is
  filled from the vendor delivery. Written to `{order}_raw_sop_violations.csv`; violations are
  non-fatal.
- **Cross-bucket trimming completeness** compares the upload against the untrimmed vendor deliveries
  listed in the notebook's `untrimmed_s3_paths` parameter — a *list*, because one experiment spans
  several Novogene orders and one order can hold several experiments. Wells match on `(wafer, UG)`,
  measured unique on a real delivery (48 CRAMs, 48 distinct UGs per wafer). The vendor layout is
  `{project}/{order}/{ExperimentID}/raw/{wafer}/`; note the segment before `raw` is the ExperimentID
  alone and carries **no** sublibrary, so the authoritative sublibrary comes from the vendor filename
  with its trailing well token stripped. Since the listed prefix is order-level, the index is filtered
  to the ExperimentID being QA'd (`ctx.order`, not the notebook's `order`, which `run_label` overrides);
  without that filter the other experiments' wells are indexed too and each reports as `not_trimmed`.
  Both segment shapes count as a match — the ExperimentID alone (`REF3`) and the older
  `{ExperimentID}_{sublibrary}` (`REF5_P01`) — because a bare equality test would orphan every well of
  an old-shape delivery. `qa_seahub_source.py` builds and merges the per-well indexes;
  `qa_seahub_recon.py` classifies them (`not_trimmed`, `orphan_trimmed`, `identity_mismatch`,
  `read_count_mismatch`, `size_not_reduced`, `duplicate_source_well`, `duplicate_trimmed_well`,
  `source_prefix_empty`, `overlapping_source_prefix`, `source_prefix_spans_experiments`,
  `source_experiment_unreadable`, `metadata_unavailable`) into
  `{order}_trimming_completeness.csv`, with per-order coverage in
  `{order}_seahub_source_coverage.csv`. Coverage matters: wells of an order missing from the list can
  only surface as `orphan_trimmed`, which reads as a completeness failure rather than incomplete
  input, so unsourced wafers are called out explicitly.
- **Rename mapping** (`qa_seahub_rename.py`) composes the per-rule repairs into one corrected S3 key
  per object in `{order}_seahub_rename_mapping.csv`. Advisory only — QA never moves, deletes or
  rewrites anything. Every proposal is re-validated against the SOP and withheld unless clean, two
  objects can never be given the same destination, and an existing destination is `blocked` rather
  than overwritten. `name_source` is `vendor` when the delivery supplied a missing token, else
  `inferred`.
- **Per-well status** rolls the above up to one verdict per well in
  `{order}_seahub_well_status.csv`: `COMPLIANT`, `RENAMEABLE` (complete, every defect repairable by
  renaming), `DATA_GAP` (an artifact genuinely absent) or `UNKNOWN` (unidentifiable, or no corrected
  name derivable). `DATA_GAP` outranks the un-nameable `UNKNOWN` so a missing CRAM stays visible even
  when its vendor order is absent from the list, and the detail names the vendor key it can be
  re-trimmed from. Whether a well is renameable is decided by `RenameProposal.renameable` — the same
  predicate that decides whether an object gets a destination in the rename CSV — so the two headline
  CSVs cannot contradict each other. Testing the defect *names* against `SEAHUB_RENAMEABLE_SOP_TYPES`
  instead is how they once did: a `sublibrary_mismatch` the vendor delivery had already repaired
  produced objects the rename CSV said to move inside a well the status CSV called `UNKNOWN`. That
  set is still maintained as documentation of which defects a rename can fix, and a test reads the
  defect literals out of the module to keep it in step.

Checksums are deliberately out of scope. S3 validates integrity on upload when the client sends a
checksum, and a multipart `ETag` is a composite of per-part digests rather than a whole-file hash, so
it cannot be compared against a lab-supplied manifest. The trimmed-vs-vendor size comparison
(`size_not_reduced`) covers the cheap part of the same question without transferring any bytes.

## Version

Current version: 1.0.0

## License

See LICENSE file in the repository root.

