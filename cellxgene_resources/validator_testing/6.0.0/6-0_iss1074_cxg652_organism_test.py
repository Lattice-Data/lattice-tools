"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
"""
import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]


ERROR_MESSAGES = (
    #*** NEED TO FILL IN ***
)

@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_passes(validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

# Testing human term in uns is valid
@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
def test_human_term_in_uns(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert "NCBITaxon:9606" == validator.adata.uns["organism_ontology_term_id"]
    assert validator.is_valid
    assert validator.errors == []

# Test organism in obs and in uns is invalid
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_organism_invalid_in_obs(validator_with_adatas, error):
    validator = validator_with_adatas
    validator.adata.obs["organism_ontology_term_id"] = validator.adata.uns["organism_ontology_term_id"]
    validator.validate_adata()
    assert validator.adata.obs["organism_ontology_term_id"] == validator.adata.uns["organism_ontology_term_id"]
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id"]}" {error}'
        ) in validator.errors

# Test organism in obs is invalid
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_organism_invalid_in_obs(validator_with_adatas, error):
    validator = validator_with_adatas
    validator.adata.obs["organism_ontology_term_id"] = validator.adata.uns["organism_ontology_term_id"]
    assert validator.adata.obs["organism_ontology_term_id"] == validator.adata.uns["organism_ontology_term_id"]
    del validator.adata.uns["organism_ontology_term_id"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id"]}" {error}'
        ) in validator.errors