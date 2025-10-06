# Testing Guide for Guidescan Pipeline

## Quick Start

```bash
# Run all regular tests (default, E2E tests skipped)
cd bcp/
conda activate lattice
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_pam_length.py -v

# Run integration tests only
pytest tests/test_integration.py -v

# Run E2E tests only (opt-in)
pytest -m e2e -v

# Run ALL tests (regular + E2E)
pytest -m "" -v
```

## Test Organization

```
bcp/tests/
├── conftest.py                 # Shared pytest fixtures and configuration
├── pytest.ini                  # Pytest settings (E2E opt-in configuration)
├── fixtures/                   # Test data files
│   ├── test_guides.txt        # Sample guides with uniform NGG PAM (5 guides)
│   ├── test_guides_mixed_pam.txt  # Guides with mixed PAMs (for error testing)
│   ├── test_guides_empty.txt  # Empty file (edge case)
│   ├── all_matches_sample.csv # Mock guidescan output
│   ├── best_matches_sample.csv # Mock best matches
│   ├── expected_all_matches.csv # Ground truth for E2E tests
│   ├── empty_genome.fa.index  # Dummy genome index for unit tests
│   ├── toy.fa.index.forward   # Real toy genome index (519KB)
│   ├── toy.fa.index.reverse   # Real toy genome index (519KB)
│   └── toy.fa.index.gs        # Real toy genome index (91B)
├── test_pam_length.py         # PAM validation tests (7 tests)
├── test_best_match.py         # Match selection tests (7 tests)
├── test_formatting.py         # Coordinate calculation tests (10 tests)
├── test_integration.py        # Integration tests with mocked guidescan (11 tests)
└── test_e2e_real_guidescan.py # E2E tests with real guidescan (8 tests, opt-in)
```

## Test Coverage

### Total: 51 tests (regular) + 8 E2E tests, all passing ✅

| Test File | Tests | Focus Area | Default |
|-----------|-------|------------|----------|
| `test_pam_length.py` | 7 | PAM validation, length calculation | ✅ Run |
| `test_best_match.py` | 8 | Match selection algorithm (retain all exact matches) | ✅ Run |
| `test_formatting.py` | 10 | Sequence splitting, coordinate calculation | ✅ Run |
| `test_integration.py` | 11 | Pipeline orchestration, file management | ✅ Run |
| `test_guide_to_gene.py` | 15 | Gene annotation from GTF files | ✅ Run |
| `test_e2e_real_guidescan.py` | 8 | Real guidescan execution, E2E workflow | ❌ Opt-in |

**Note**: E2E tests are skipped by default. Run with `pytest -m e2e` to include them.

## What We Test

### ✅ Tested Components

1. **PAM Validation** (`test_pam_length.py`)
   - Uniform PAM sequences (NGG, NNGG, NNNGG)
   - Mixed PAM error detection
   - Empty file handling
   - Missing column error detection
   - Length calculation from PAM string

2. **Match Selection** (`test_best_match.py`)
   - Single exact match selection
   - **Multiple exact matches - retains ALL** (not just one)
   - **Non-exact matches result in NA** (treated as no match)
   - Prefers exact matches over any non-exact matches
   - Handles all NA results
   - Invalid distance value handling

3. **Coordinate Formatting** (`test_formatting.py`)
   - Sequence splitting (protospacer + PAM)
   - Variable PAM length support (3bp, 4bp, 5bp)
   - + strand coordinate calculation
   - - strand coordinate calculation
   - Exact match filtering (distance=0 only)
   - No match handling (NA values)
   - Output file structure

4. **Integration** (`test_integration.py`)
   - Pipeline initialization
   - Output directory creation
   - Input validation (missing files)
   - PAM length setting during validation
   - Cleanup functionality
   - End-to-end workflow
   - Mocked guidescan execution
   - Error propagation

5. **End-to-End (E2E)** (`test_e2e_real_guidescan.py`) - **Opt-in with `-m e2e`**
   - Real guidescan binary execution
   - Direct command-line execution
   - Full pipeline with toy genome
   - Output property validation
   - Cleanup with real files
   - Mixed match/no-match handling
   - Genomic coordinate accuracy

### ❌ NOT Tested in Regular Tests

- Actual guidescan binary execution (see E2E tests above)
- Real genome alignment (see E2E tests above)
- Biological accuracy of matches
- Large-scale genome performance

## Testing Strategy

### Unit Tests
Test individual functions in isolation using:
- Synthetic data
- pandas DataFrames
- Temporary directories

### Integration Tests
Test workflow orchestration using:
- Mocked guidescan output
- Fixture CSV files
- Temporary file system

### Mocking Approach
```python
# Mock guidescan execution
with patch.object(GuidescanPipeline, 'run_guidescan', mock_function):
    result = pipeline.run()

# Mock subprocess calls
with patch('subprocess.run'):
    pipeline.validate_inputs()
```

## Shared Fixtures (conftest.py)

Available in all test files:

```python
def test_example(temp_dir, dummy_genome_index, sample_guides_file_ngg):
    """Use shared fixtures across tests."""
    # temp_dir: Temporary directory (auto-cleanup)
    # dummy_genome_index: Fake genome.fa.index file
    # sample_guides_file_ngg: Sample guides with NGG PAM
    pass
```

## Running Specific Test Groups

```bash
# Run only PAM validation tests
pytest tests/test_pam_length.py -v

# Run only integration tests
pytest tests/test_integration.py -v

# Run only E2E tests (requires guidescan installed)
pytest -m e2e -v

# Run ALL tests including E2E
pytest -m "" -v

# Run tests matching a pattern
pytest -k "pam" -v

# Run tests with coverage (regular tests only)
pytest --cov=. --cov-report=html

# Run with detailed output on failure
pytest -vv --tb=long

# Stop after first failure
pytest -x
```

## Continuous Integration

### Fast CI (Skip E2E)
```yaml
# Example GitHub Actions - Fast unit/integration tests only
- name: Run Tests
  run: |
    conda activate lattice
    cd bcp/
    pytest --tb=short  # E2E tests skipped by default
```

### Comprehensive CI (Include E2E)
```yaml
# Example GitHub Actions - All tests including E2E
- name: Install Guidescan
  run: |
    # Install guidescan2 binary
    
- name: Run All Tests
  run: |
    conda activate lattice
    cd bcp/
    pytest -m "" --tb=short  # Include E2E tests
```

## Adding New Tests

### 1. Create Test File

```python
# tests/test_new_feature.py
import pytest
from guidescan_pipeline import GuidescanPipeline

class TestNewFeature:
    def test_new_functionality(self, temp_dir, dummy_genome_index):
        """Test description."""
        # Arrange
        pipeline = GuidescanPipeline(...)
        
        # Act
        result = pipeline.new_method()
        
        # Assert
        assert result == expected
```

### 2. Run Tests

```bash
pytest tests/test_new_feature.py -v
```

### 3. Verify Coverage

```bash
pytest --cov=. --cov-report=term-missing
```

## Test Data Formats

### Input Guides Format
```csv
id,sequence,pam,chromosome,position,sense
guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,
```

### Guidescan Output Format (all_matches.csv)
```csv
id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
```

### Final Output Format (exact_matches_formatted.csv)
```csv
id,sequence,pam,chromosome,start,end,sense
guide1,TGCCTCGCGCAGCTCGCGG,NGG,1,1000,1018,+
```

## Troubleshooting Tests

### Import Errors
```bash
# Ensure you're in the right directory
cd bcp/
pytest

# Or specify Python path
PYTHONPATH=. pytest tests/
```

### Fixture Not Found
```bash
# Check conftest.py is present
ls tests/conftest.py

# Clear pytest cache
rm -rf .pytest_cache
pytest --cache-clear
```

### Tests Hang
```bash
# Run with timeout
pytest --timeout=30

# Stop after first failure
pytest -x
```

## Test Maintenance

- **After code changes**: Run full test suite
- **Before committing**: Ensure all tests pass
- **After adding features**: Add corresponding tests
- **Monthly**: Review test coverage and add missing tests

## Performance

Current test suite performance:

| Test Suite | Tests | Time | Default |
|------------|-------|------|----------|
| **Regular (Unit + Integration)** | 51 | ~0.6s | ✅ Run |
| **E2E (Real Guidescan)** | 8 | ~0.5s | ❌ Skip |
| **Total** | 59 | ~0.9s | 51 run, 8 skip |

- All tests use temporary directories (auto-cleanup)
- Regular tests use mocked guidescan (fast)
- E2E tests use real guidescan with toy genome (small, fast)
- Suitable for CI/CD pipelines

## E2E Tests (Opt-in)

For comprehensive end-to-end testing with the actual guidescan binary, see:
- **[E2E_TESTS.md](E2E_TESTS.md)** - Complete E2E testing documentation
- Tests execute real guidescan with toy genome index
- Validates actual coordinate calculation
- Opt-in design: run with `pytest -m e2e`
- Auto-skips if guidescan not installed

## Questions?

- Check test documentation in each test file
- Review `conftest.py` for available fixtures
- See `README.md` for pipeline usage examples
- See `E2E_TESTS.md` for E2E testing details

