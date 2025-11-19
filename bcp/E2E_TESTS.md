# End-to-End (E2E) Tests with Real Guidescan

## Overview

E2E tests execute the actual guidescan binary with a toy genome index to verify the complete pipeline works end-to-end. These tests are **opt-in** and skipped by default.

## Running E2E Tests

### Run Only E2E Tests
```bash
cd bcp/
conda activate lattice
pytest -m e2e -v
```

### Run ALL Tests (Regular + E2E)
```bash
pytest -m ""
```

### Skip E2E Tests (Default)
```bash
pytest  # E2E tests automatically skipped
```

## Test Coverage

**8 E2E Tests** covering:

### 1. Guidescan Availability (2 tests)
- ✅ Guidescan executable exists
- ✅ Guidescan enumerate command is available

### 2. Direct Guidescan Execution (2 tests)
- ✅ Run guidescan command directly
- ✅ Verify output format has expected columns

### 3. Full Pipeline with Real Guidescan (3 tests)
- ✅ Complete pipeline execution
- ✅ Output properties validation (retains all exact matches)
- ✅ Cleanup functionality with real files

### 4. Edge Cases (1 test)
- ✅ Mixed results (guides with and without exact matches)

## Test Data

### Toy Genome Index Files
Located in `tests/fixtures/`:
- `toy.fa.index.forward` (519KB, 2305 lines)
- `toy.fa.index.reverse` (519KB, 2193 lines)
- `toy.fa.index.gs` (91B, 21 lines)

### Test Guides
File: `tests/fixtures/test_guides.txt`

Contains 5 guides:
- **guide1-4**: No matches in toy genome (returns NA)
- **guide5**: Has exact match on chromosome 9, position 93966

### Expected Output
File: `tests/fixtures/expected_all_matches.csv`

Pre-generated guidescan output for comparison and validation.

## What E2E Tests Verify

### ✅ Tested
1. **Real Guidescan Execution**
   - Subprocess call succeeds
   - Command-line arguments correct
   - Output files created

2. **Output Format**
   - All expected columns present
   - CSV structure valid
   - Data types correct

3. **Pipeline Integration**
   - Guidescan → Match Selection → Format workflow
   - **All exact matches retained** (not just one per guide)
   - Non-exact matches result in NA
   - Intermediate files created correctly
   - Final output has correct structure

4. **Coordinate Calculation**
   - Start/end positions valid
   - Strand orientation correct
   - Protospacer/PAM splitting works
   - Multiple rows output for guides with multiple exact matches

5. **Mixed Results**
   - Handles guides with exact matches
   - Handles guides without exact matches (NA)
   - Both types in same output

6. **File Management**
   - Cleanup removes intermediate files
   - Final output preserved
   - Output directory creation

### ❌ NOT Tested (Per Requirements)
- Biological accuracy of matches
- Guidescan algorithm correctness
- Large-scale genome performance

## Requirements

1. **guidescan must be installed**
   ```bash
   guidescan --help  # Verify installation
   ```

2. **Toy genome index files** must exist in `tests/fixtures/`
   - Automatically skipped if guidescan not available
   - Automatically skipped if index files missing

## Test Execution Time

- **E2E tests**: ~0.5 seconds total
- **Regular tests**: ~0.6 seconds total
- **All tests combined**: ~0.9 seconds total

Very fast because toy genome is small (~500KB) and tests are well-optimized.

## CI/CD Integration

### Skip E2E in CI (Fast)
```yaml
- name: Run Unit Tests
  run: |
    pytest  # E2E skipped by default
```

### Run E2E in CI (Comprehensive)
```yaml
- name: Install Guidescan
  run: |
    # Install guidescan2
    
- name: Run All Tests
  run: |
    pytest -m ""  # Include E2E tests
```

## Test Results Example

```
tests/test_e2e_real_guidescan.py::TestGuidescanAvailability::test_guidescan_executable_exists PASSED [ 12%]
tests/test_e2e_real_guidescan.py::TestGuidescanAvailability::test_guidescan_enumerate_command_exists PASSED [ 25%]
tests/test_e2e_real_guidescan.py::TestDirectGuidescanExecution::test_direct_guidescan_command PASSED [ 37%]
tests/test_e2e_real_guidescan.py::TestDirectGuidescanExecution::test_guidescan_output_format PASSED [ 50%]
tests/test_e2e_real_guidescan.py::TestFullPipelineWithRealGuidescan::test_full_pipeline_execution PASSED [ 62%]
tests/test_e2e_real_guidescan.py::TestFullPipelineWithRealGuidescan::test_pipeline_output_properties PASSED [ 75%]
tests/test_e2e_real_guidescan.py::TestFullPipelineWithRealGuidescan::test_pipeline_with_cleanup PASSED [ 87%]
tests/test_e2e_real_guidescan.py::TestE2EEdgeCases::test_guides_with_matches_and_no_matches PASSED [100%]

================= 8 passed, 51 deselected, 4 warnings in 0.48s =================
```

## Troubleshooting

### Tests Skipped
```
tests/test_e2e_real_guidescan.py::test_xxx SKIPPED (guidescan not available in PATH)
```
**Solution**: Install guidescan2 and ensure it's in PATH

### Toy Genome Not Found
```
AssertionError: assert False  # Index files don't exist
```
**Solution**: Ensure toy genome index files are in `tests/fixtures/`

### All Guides Return NA
**Expected**: Guides 1-4 should return NA (no matches in toy genome)
**Unexpected**: Guide5 should have a match on chr9

If guide5 returns NA, check that test_guides.txt includes the 5th guide.

## Implementation Notes

### Key Features

1. **Opt-in Design**
   - E2E tests marked with `@pytest.mark.e2e`
   - pytest.ini configured to skip by default: `-m "not e2e"`
   - Must explicitly request: `pytest -m e2e`

2. **Auto-Skip Logic**
   ```python
   pytestmark = [
       pytest.mark.e2e,
       pytest.mark.skipif(
           not is_guidescan_available(),
           reason="guidescan not available in PATH"
       )
   ]
   ```

3. **Backward Compatible Validation**
   - Accepts real genome indices (`.forward`, `.reverse`, `.gs` files)
   - Also accepts dummy indices (single file) for unit tests
   - Ensures E2E tests don't break regular tests

4. **Float Parsing Fix**
   - Handles positions as floats or ints: `int(float(position_str))`
   - Prevents errors when pandas reads positions as floats

## Complete Test Suite Summary

| Test Type | Count | Speed | Default |
|-----------|-------|-------|---------|
| Unit Tests | 40 | 0.4s | ✅ Run |
| Integration Tests (Mocked) | 11 | 0.2s | ✅ Run |
| **E2E Tests (Real Guidescan)** | **8** | **0.5s** | **❌ Skip** |
| **Total** | **59** | **0.9s** | **51 run, 8 skip** |

## Future Enhancements

Potential additions:
1. E2E tests with larger genomes (performance testing)
2. E2E tests with different PAM types
3. E2E tests with multiple mismatches
4. E2E tests with off-target analysis
5. Comparison tests (output vs. known good results)

## Questions?

- See `README.md` for general pipeline documentation
- See `TESTING.md` for regular testing guide
- Check pytest output for detailed error messages

