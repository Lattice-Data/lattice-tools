# CRISPR Guide Processing Pipeline

A Python pipeline for processing CRISPR guide sequences using guidescan2 to find genomic matches and calculate coordinates.

## Overview

This pipeline takes CRISPR guide sequences as input, uses guidescan2 to find genomic matches, selects the best match for each guide, and outputs formatted genomic coordinates.

**Key Features:**
- Validates PAM uniformity across all guides
- Supports variable PAM lengths (NGG, NNGG, etc.)
- Finds exact genomic matches using guidescan2
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

Contains genomic coordinates for exact matches:

```csv
id,sequence,pam,chromosome,start,end,sense
guide1,TGCCTCGCGCAGCTCGCGG,NGG,1,1000,1018,+
guide2,GAGTTCGCTGCGCGCTGTT,NGG,3,1982,2000,-
guide3,GCCCGCTCCCCGCGATCCC,NGG,NA,NA,NA,NA
```

**Column descriptions:**
- `id`: Guide identifier
- `sequence`: Protospacer sequence
- `pam`: PAM sequence
- `chromosome`: Chromosome where exact match found (NA if no match)
- `start`, `end`: Genomic coordinates (1-based, inclusive)
- `sense`: Strand orientation (+ or -, NA if no match)

### Intermediate Files (if `--keep-intermediate` is used)

1. **`all_matches.csv`**: All genomic matches from guidescan
2. **`best_matches.csv`**: Best match selected for each guide

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
2025-01-15 10:30:15 - INFO - Selecting best matches for each guide...
2025-01-15 10:30:15 - INFO - Best matches saved to results/best_matches.csv
2025-01-15 10:30:15 - INFO - Formatting exact matches...
2025-01-15 10:30:15 - INFO - Processing complete. Output saved to: results/exact_matches_formatted.csv
2025-01-15 10:30:15 - INFO - Pipeline completed successfully!
2025-01-15 10:30:15 - INFO - Summary: 3 guides processed, 2 exact matches found
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
    verbose=True
)

# Run complete pipeline
results_df = pipeline.run()

# Access results
print(f"Total guides: {len(results_df)}")
print(f"Exact matches: {len(results_df[results_df['chromosome'] != 'NA'])}")

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

**Test Statistics**: 43 total tests (35 regular + 8 E2E)

### Test Organization

```
bcp/tests/
├── conftest.py                 # Shared fixtures
├── fixtures/                   # Test data files
│   ├── test_guides.txt        # Sample guides (5 guides)
│   ├── toy.fa.index.*         # Toy genome index for E2E tests
│   ├── all_matches_sample.csv
│   ├── best_matches_sample.csv
│   └── expected_all_matches.csv
├── test_pam_length.py         # PAM validation tests (7 tests)
├── test_best_match.py         # Match selection tests (7 tests)
├── test_formatting.py         # Coordinate calculation tests (10 tests)
├── test_integration.py        # Integration tests (11 tests)
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

3. **Best Match Selection**
   - Group matches by guide ID
   - Select match with lowest edit distance
   - Prefer actual genomic matches over NA

4. **Coordinate Formatting**
   - Filter for exact matches (distance = 0)
   - Split sequence into protospacer + PAM
   - Calculate genomic coordinates:
     - **+ strand**: start = position, end = position + length - 1
     - **- strand**: end = position, start = position - length + 1

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

## Contributing

When modifying the pipeline:

1. Write tests for new functionality
2. Run the test suite: `pytest`
3. Update this README if adding features
4. Follow existing code style (use Ruff for linting)

## Version

Current version: 1.0.0

## License

See LICENSE file in the repository root.

