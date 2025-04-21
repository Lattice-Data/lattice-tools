"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1303
PR:
"""
import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]
NON_HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" not in file]

HUMAN_VALID_VALUES = [
    "HANCESTRO:0005",
    "unknown",
    "HANCESTRO:0005 || HANCESTRO:0014",
    "HANCESTRO:0005 || HANCESTRO:0014 || HANCESTRO:0015"
]

HUMAN_INVALID_VALUES = [
    "HANCESTRO:0014 || HANCESTRO:00005",
    "HANCESTRO:00005,HANCESTRO:0014",
    "HANCESTRO:00005||HANCESTRO:0014 ",
    "HANCESTRO:0005 || unknown",
    "HANCESTRO:0005 || HANCESTRO:0005",
    "HANCESTRO:0005 ",
    "HANCESTRO:0560",
    "HANCESTRO:0031",
    "na"
]

ERROR_MESSAGES = (
    #*** NEED TO FILL IN ***
)

# Testing human h5ads
@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
@pytest.mark.parametrize("valid_term", HUMAN_VALID_VALUES)
def test_human_ethnicity_vaild_values(validator_with_adatas,valid_term):
    validator = validator_with_adatas
    validator.adata.obs['self_reported_ethnicity_term_id'] = valid_term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []

@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
@pytest.mark.parametrize("invalid_term", HUMAN_INVALID_VALUES)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_human_ethnicity_invaild_values(validator_with_adatas,invalid_term,error):
    validator = validator_with_adatas
    validator.adata.obs['self_reported_ethnicity_term_id'] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors

# Testing non-human h5ads
@pytest.mark.parametrize("test_h5ads", NON_HUMAN_H5ADS)
@pytest.mark.parametrize("valid_term", ('na'))
def test_ethnicity_vaild_values(validator_with_adatas,valid_term):
    validator = validator_with_adatas
    validator.adata.obs['self_reported_ethnicity_term_id'] = valid_term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []

@pytest.mark.parametrize("test_h5ads", NON_HUMAN_H5ADS)
@pytest.mark.parametrize("invalid_term", HUMAN_VALID_VALUES)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_ethnicity_invaild_values(validator_with_adatas,invalid_term,error):
    validator = validator_with_adatas
    validator.adata.obs['self_reported_ethnicity_term_id'] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors