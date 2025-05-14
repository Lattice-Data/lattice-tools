"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/687
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1356

Should pass:
(Y) various lists of valid MONDO terms
    (y) just normal term
    (y) just 1 valid MONDO term
(Y) add labels check: obs.disease order matches disease_ontology_term_id order

Should not pass:
(Y) valid MONDO terms out of lexical order
(Y) valid MONDO terms with different delimiters - e.g. || (without spaces)
(Y) normal term w/ valid MONDO term
(Q) invalid MONDO term + valid MONDO term
"""


import pytest
import anndata as ad
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

VALID_VALUES = [
    "PATO:0000461",
    "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008",  # in lexical order
    "MONDO:0004604"
    ]

INVALID_VALUES = [
    "PATO:0000461 ",
    "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008 || MONDO:1030008",  # last term duplicated
    "MONDO:0043004 || MONDO:0004604 || MONDO:1030008 || MONDO:0800349",  # out of lexical order
    "MONDO:0004604 , MONDO:0043004 , MONDO:0800349 , MONDO:1030008",  # different delimiter
    "MONDO:0004604||MONDO:0043004||MONDO:0800349||MONDO:1030008",  # no spaces
    "PATO:0000461 || MONDO:0004604",  # normal term w/ valid MONDO term
    "MONDO:0100355 || MONDO:0043004"  # invalid MONDO term + valid MONDO term
    ]

error_message_suffix = " Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of " \
    + "'MONDO:0000001' (disease) are allowed. Multiple MONDO terms are supported if in ascending lexical order with the delimiter ` || `."


ERROR_MESSAGE = [
    "in 'disease_ontology_term_id' contains duplicates." + error_message_suffix,
    "in 'disease_ontology_term_id' is not in ascending lexical order." + error_message_suffix,
    "in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'." + error_message_suffix
    ]


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
def test_passes(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("valid_term", VALID_VALUES)
def test_disease_valid_values(validator_with_adatas,valid_term):

    # various lists of valid MONDO terms => valid

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = valid_term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("invalid_term", INVALID_VALUES[1:2])
@pytest.mark.parametrize("error", ERROR_MESSAGE[0:1])
def test_disease_invalid_duplicates(validator_with_adatas,invalid_term,error):

    # valid MONDO terms with duplicates => invalid

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("invalid_term", INVALID_VALUES[2:3])
@pytest.mark.parametrize("error", ERROR_MESSAGE[1:2])
def test_disease_invalid_order(validator_with_adatas,invalid_term,error):

    # valid MONDO terms out of lexical order => invalid

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("invalid_term", INVALID_VALUES[3:5])
@pytest.mark.parametrize("error", ERROR_MESSAGE[2:3])
def test_disease_invalid_delimiter(validator_with_adatas,invalid_term,error):

    # valid MONDO terms with different delimiters - e.g. || (without spaces and commas) => invalid

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("invalid_term", INVALID_VALUES[5:6])
@pytest.mark.parametrize("error", ERROR_MESSAGE[2:3])
def test_disease_valid_term_with_normal_term(validator_with_adatas,invalid_term,error):

     # normal term w/ valid MONDO term => invalid

    ### Validation error is confusing here and should probably be fixed.

    # Correct me if I'm wrong, but the error message is slicing out PATO:0000461 as an invalid term and
    # then also erroring on the lexical order of the terms:

        # First error:
        # "ERROR: 'PATO:0000461' in 'disease_ontology_term_id' is not a valid ontology term id
        # of 'MONDO'. Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof,
        # or descendant terms of 'MONDO:0000001' (disease) are allowed. Multiple MONDO terms are supported if in
        # ascending lexical order with the delimiter ` || `.",

        # Second error:
        # "ERROR: 'PATO:0000461 || MONDO:0004604'
        # in 'disease_ontology_term_id' is not in ascending lexical order. Only 'PATO:0000461' (normal),
        # 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease)
        # are allowed. Multiple MONDO terms are supported if in ascending lexical order with the delimiter ` || `."]"

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = invalid_term
    validator.validate_adata()
    assert not validator.is_valid
    assert (f"ERROR: '{invalid_term}' {error}"
        ) in validator.errors

@pytest.mark.parametrize("test_h5ads", ALL_H5ADS[0:1])
@pytest.mark.parametrize("test_term", VALID_VALUES[1:2])
def test_disease_add_label_order_valid_values(validator_with_adatas,test_term):

    # add_labels check: obs.disease order matches disease_ontology_term_id order

    validator = validator_with_adatas
    validator.adata.obs["disease_ontology_term_id"] = validator.adata.obs["disease_ontology_term_id"].cat.add_categories(test_term)
    validator.adata.obs.loc[validator.adata.obs.index[0], "disease_ontology_term_id"] = test_term
    validator.validate_adata()
    labeler = AnnDataLabelAppender(validator.adata)
    labeler._add_labels()
    map = { "Hodgkin's lymphoma, lymphocytic-histiocytic predominance":"MONDO:0004604",
            "Weil's disease":"MONDO:0043004",
            "atrial fibrillation, familial, 16":"MONDO:0800349",
            "mitral valve insufficiency":"MONDO:1030008"}
    assert labeler.adata.obs["disease_ontology_term_id"][0].split(" || ") == [map[t] for t in labeler.adata.obs["disease"][0].split(" || ")]
