"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue:

Testing conditions:
Should not pass
- absent uns.organism_ontology_term_id X
- present obs.organism X
- present obs.organism_ontology_term_id X
- present uns.organism X
- present uns.organism_ontology_term_id_colors
- present uns.organism_ontology_colors
- uns.organism_ontology_term_id - not on accepted organisms (esp COVID)

Should pass
- uns.organism_ontology_term_id - from list of accepted organisms
"""

import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]

ACCEPTED_ORGANISMS = {"NCBITaxon:6293","NCBITaxon:9483","NCBITaxon:7955","NCBITaxon:7227","NCBITaxon:9595",
    "NCBITaxon:9606","NCBITaxon:9541","NCBITaxon:9544","NCBITaxon:30608","NCBITaxon:10090","NCBITaxon:9986",
    "NCBITaxon:9598","NCBITaxon:10116","NCBITaxon:9823","NCBITaxon:1747270","NCBITaxon:1654737",
    "NCBITaxon:947985","NCBITaxon:230741", "NCBITaxon:756884","NCBITaxon:947987","NCBITaxon:1170810"}

EXEMPT_ORGANISMS = {"NCBITaxon:2697049","NCBITaxon:2697049"}  #covid, spike-ins

ERROR_MESSAGES = (
    #*** NEED TO FILL IN ***
)

# All fixtures should be valid
@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_passes(validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

# Test human term in uns is valid
@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_human_term_in_uns(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert "NCBITaxon:9606" == validator.adata.uns["organism_ontology_term_id"]
    assert validator.is_valid
    assert validator.errors == []

# Test organism not in uns is invalid
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_organism_term_not_in_uns(validator_with_adatas, error):
    validator = validator_with_adatas
    del validator.adata.uns["organism_ontology_term_id"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id"]}" {error}'
        ) in validator.errors

# Test human-readable term in uns and obs is invalid  -> should this check more than just human?
@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
def test_human_in_uns(validator_with_adatas,error):
    validator = validator_with_adatas
    validator.adata.uns["organism"] = "Homo sapiens"
    validator.adata.obs["organism"] = "Homo sapiens"
    assert "Homo sapiens" == validator.adata.uns["organism"]
    assert "Homo sapiens" == validator.adata.obs["organism"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.uns["organism"]}" {error}'
            f'ERROR: "{validator.adata.obs["organism"]}" {error}'
        ) in validator.errors

# Test colors in uns is invalid
@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_colors_in_uns(validator_with_adatas,error):
    validator = validator_with_adatas
    validator.adata.uns["organism_ontology_colors"] = ['aqua']
    validator.adata.uns["organism_ontology_term_id_colors"] = ['aqua']
    assert ['aqua'] == validator.adata.uns["organism_ontology_colors"]
    assert ['aqua'] == validator.adata.uns["organism_ontology_term_id_colors"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.uns["organism_ontology_colors"]}" {error}'
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id_colors"]}" {error}'
        ) in validator.errors

# Test organism in obs and uns is invalid
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_organism_in_obs_and_uns(validator_with_adatas, error):
    validator = validator_with_adatas
    validator.adata.obs["organism_ontology_term_id"] = validator.adata.uns["organism_ontology_term_id"]
    assert validator.adata.obs["organism_ontology_term_id"] == validator.adata.uns["organism_ontology_term_id"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id"]}" {error}'
        ) in validator.errors

# Test organism just in obs is invalid
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_organism_in_obs(validator_with_adatas, error):
    validator = validator_with_adatas
    validator.adata.obs["organism_ontology_term_id"] = validator.adata.uns["organism_ontology_term_id"]
    assert validator.adata.obs["organism_ontology_term_id"] == validator.adata.uns["organism_ontology_term_id"]
    del validator.adata.uns["organism_ontology_term_id"]
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.obs["organism_ontology_term_id"]}" {error}'
        ) in validator.errors

# Test other organism terms in uns are valid
@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
@pytest.mark.parametrize("organismss", ACCEPTED_ORGANISMS)
def test_organism_term_in_uns(validator_with_adatas,organisms):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.adata.uns["organism_ontology_term_id"] in organisms  # will this cover the descendants in accepted_organisms set?
    assert validator.is_valid
    assert validator.errors == []

# Test exempt organism in uns is invalid
@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
@pytest.mark.parametrize("organism", EXEMPT_ORGANISMS)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_e_organism_term_in_uns(validator_with_adatas,organisms,error):
    validator = validator_with_adatas
    validator.adata.uns['organism_ontology_term_id'] = organisms
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f'ERROR: "{validator.adata.uns["organism_ontology_term_id"]}" {error}'
        ) in validator.errors