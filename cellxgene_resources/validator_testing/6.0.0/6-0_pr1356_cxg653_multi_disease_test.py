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
(Y) valid MONDO terms with different delimiters - e.g. || without spaces
(Y) valid MONDO terms with different delimiters - e.g. commas
(Y) normal term w/ valid MONDO term
(Y) invalid MONDO term + valid MONDO term
(Y) valid MONDO terms with duplicates
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

VALID_VALUES = {
    "normal": "PATO:0000461",
    "in lexical order": "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008",  # in lexical order
    "valid MONDO": "MONDO:0004604"
    }

INVALID_VALUES = {
    "normal term with space": "PATO:0000461 ",
    "last term duplicated": "MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008 || MONDO:1030008",
    "out of lexical order": "MONDO:0043004 || MONDO:0004604 || MONDO:1030008 || MONDO:0800349",
    "different delimiter": "MONDO:0004604 , MONDO:0043004 , MONDO:0800349 , MONDO:1030008",
    "no spaces": "MONDO:0004604||MONDO:0043004||MONDO:0800349||MONDO:1030008",
    "invalid MONDO term": "MONDO:0100355"
    }

error_message_suffix = " Individual terms 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of " \
    + "'MONDO:0000001' (disease) are allowed. Multiple terms are supported if in ascending lexical order with the delimiter ` || ` if all terms are valid MONDO terms."


ERROR_MESSAGE = {
    "contains duplicates": "in 'disease_ontology_term_id' contains duplicates." + error_message_suffix,
    "not in lexical order": "in 'disease_ontology_term_id' is not in ascending lexical order." + error_message_suffix,
    "not a valid term": "in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'." + error_message_suffix,
    "not an allowed term": "in 'disease_ontology_term_id' is not an allowed term id." + error_message_suffix,
    "not a valid MONDO term": "in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO'." + error_message_suffix
    }


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestDiseaseOntologyValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_passes(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("valid_term", VALID_VALUES.values())
    def test_disease_valid_values(self, valid_term):

        # various lists of valid MONDO terms => valid

        self.validator.adata.obs["disease_ontology_term_id"] = valid_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []

    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["normal term with space"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_disease_valid_term_with_space(self, invalid_term, error):

        # valid MONDO terms with space => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["last term duplicated"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["contains duplicates"]])
    def test_disease_invalid_duplicates(self, invalid_term, error):

        # valid MONDO terms with duplicates => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["out of lexical order"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not in lexical order"]])
    def test_disease_invalid_order(self, invalid_term, error):

        # valid MONDO terms out of lexical order => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["different delimiter"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_disease_invalid_delimiter(self, invalid_term, error):

        # valid MONDO terms with different delimiters - e.g. || (commas) => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["no spaces"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_disease_invalid_no_spaces(self, invalid_term, error):

        # valid MONDO terms with different delimiters - e.g. || (without spaces) => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("valid_term", [VALID_VALUES["valid MONDO"]])
    @pytest.mark.parametrize("normal_term", [VALID_VALUES["normal"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid MONDO term"]])
    def test_disease_valid_term_with_normal_term(self, valid_term,normal_term, error):

        # normal term w/ valid MONDO term => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = valid_term + " || " + normal_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{normal_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [INVALID_VALUES["invalid MONDO term"]])
    @pytest.mark.parametrize("valid_term", [VALID_VALUES["valid MONDO"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not an allowed term"]])
    def test_disease_invalid_valid_terms(self, invalid_term,valid_term, error):

        # invalid MONDO term with valid MONDO term => invalid

        self.validator.adata.obs["disease_ontology_term_id"] = valid_term + " || " + invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("test_term", [VALID_VALUES["in lexical order"]])
    def test_disease_add_label_order_valid_values(self, test_term):

        # add_labels check: obs.disease order matches disease_ontology_term_id order

        self.validator.adata.obs["disease_ontology_term_id"] = self.validator.adata.obs["disease_ontology_term_id"].cat.add_categories(test_term)
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "disease_ontology_term_id"] = test_term
        self.validator.validate_adata()
        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        map = { "Hodgkin's lymphoma, lymphocytic-histiocytic predominance":"MONDO:0004604",
                "Weil's disease":"MONDO:0043004",
                "atrial fibrillation, familial, 16":"MONDO:0800349",
                "mitral valve insufficiency":"MONDO:1030008"}
        assert labeler.adata.obs["disease_ontology_term_id"][0].split(" || ") == [map[t] for t in labeler.adata.obs["disease"][0].split(" || ")]

    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not an allowed term"]])
    def test_disease_invalid_term(self, error):

        # invalid MONDO term -> fail

        self.validator.adata.obs["disease_ontology_term_id"] = "MONDO:0012153"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{self.validator.adata.obs['disease_ontology_term_id'].unique()[0]}' {error}"
        ) in self.validator.errors