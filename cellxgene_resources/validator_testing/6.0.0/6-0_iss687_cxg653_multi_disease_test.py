"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/687
PR:
"""
import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)


VALID_VALUES = [
    "PATO:0000461",
    "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008",  # in lexical order
    "MONDO:0004604",
]

INVALID_VALUES = [
    "PATO:0000461 ",
    "MONDO:0043004 || MONDO:0004604 || MONDO:1030008 || MONDO:0800349",  # out of lexical order
    "MONDO:0004604 , MONDO:0043004 , MONDO:0800349 , MONDO:1030008",  # different delimiter
    "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008 || MONDO:1030008",  # last term duplicated
    "MONDO:0004604||MONDO:0043004||MONDO:0800349||MONDO:1030008",  # no spaces
    "PATO:0000461 || MONDO:0004604"  # normal term w/ valid MONDO term
    "MONDO:0100355 || MONDO:0043004"  # invalid MONDO term + valid MONDO term

]

ERROR_MESSAGES = (
    #*** NEED TO FILL IN ***
)

def test_passes(validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

@pytest.mark.parametrize("valid_term", VALID_VALUES)
def test_disease_vaild_values(validator_with_adatas,valid_term):
    validator = validator_with_adatas
    validator.adata.obs['disease_ontology_term_id'] = valid_term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []

@pytest.mark.parametrize("invalid_term", INVALID_VALUES)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_disease_invaild_values(validator_with_adatas,invalid_term,error):
    validator = validator_with_adatas
    validator.adata.obs['disease_ontology_term_id'] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (
            f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors

# add-labels check: obs.disease order matches disease_ontology_term_id order