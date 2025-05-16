"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1303
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1356

Should pass:
() various lists of valid HANCESTRO terms
    () just normal term
    () just 1 valid HANCESTRO term
() add labels check: obs.disease order matches  self_reported_ethnicity_ontology_term_id order

Should not pass:
() valid HANCESTRO terms out of lexical order
() valid HANCESTRO terms with different delimiters - e.g. || (without spaces)
() normal term w/ valid HANCESTRO term
() invalid HANCESTRO term + valid HANCESTRO term

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

HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]
NON_HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" not in file]

HUMAN_VALID_VALUES = [
    "HANCESTRO:0005",  # single valid term
    "unknown",  # valid human unknown
    "HANCESTRO:0005 || HANCESTRO:0014",  # in lexical order
    "HANCESTRO:0005 || HANCESTRO:0014 || HANCESTRO:0015"  # multiple valid terms
]

HUMAN_INVALID_VALUES = [
    "HANCESTRO:0014 || HANCESTRO:00005",  # out of lexical order
    "HANCESTRO:00005,HANCESTRO:0014",  # different delimiter
    "HANCESTRO:00005||HANCESTRO:0014",  # no spaces
    "HANCESTRO:0005 || unknown",  # valid term with unknown
    "HANCESTRO:0005 || HANCESTRO:0005",  # duplicate terms
    "HANCESTRO:0005 ",  # valid term with whitespace typo
    "HANCESTRO:0031",  # forbidden term
    "na"  # invalid for human
]

NON_HUMAN_VALID_VALUES = ["na"] # nonhuman valid term

NON_HUMAN_INVALID_VALUES = ["HANCESTRO:0005"] # nonhuman invalid term

error_message_suffix = " When 'organism_ontology_term_id' is 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be formatted as one or more HANCESTRO terms in ascending lexical order with the delimiter ` || `, or 'unknown' if unavailable. Cannot match any forbidden HANCESTRO terms listed in schema definition."

ERROR_MESSAGE = [
    "in 'self_reported_ethnicity_ontology_term_id' contains duplicates." + error_message_suffix,  # contains duplicates
    "in 'self_reported_ethnicity_ontology_term_id' is not in ascending lexical order." + error_message_suffix,  # not in lexical order
    "in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id of 'HANCESTRO'." + error_message_suffix  # not a valid term
]

NON_HUMAN_ERROR_MESSAGE = ["in 'self_reported_ethnicity_ontology_term_id' is not a valid value of 'self_reported_ethnicity_ontology_term_id'. When 'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'."]

@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
class TestHumanEthnicityOntologyValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_passes(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("valid_term", HUMAN_VALID_VALUES)
    def test_ethnicity_valid_values(self, valid_term):

        # various lists of valid HANCESTRO terms => valid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = valid_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES[4]])  # duplicate terms
    @pytest.mark.parametrize("error", [ERROR_MESSAGE[0]])  # contains duplicates error
    def test_ethnicity_invalid_duplicates(self, invalid_term, error):

        # valid HANCESTRO terms with duplicates => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        print(self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"])
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES[0]])  # out of lexical order
    @pytest.mark.parametrize("error", [ERROR_MESSAGE[1]])  # not in lexical order error
    def test_ethnicity_invalid_order(self, invalid_term, error):

        # valid HANCESTRO terms out of lexical order => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES[1]])  # different delimiter
    @pytest.mark.parametrize("error", [ERROR_MESSAGE[2]])  # not a valid term error
    def test_ethnicity_invalid_delimiter(self, invalid_term, error):

        # valid HANCESTRO terms with different delimiters - e.g. || (without spaces and commas) => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES[2]])  # no spaces
    @pytest.mark.parametrize("error", [ERROR_MESSAGE[2]])  # not a valid term error
    def test_ethnicity_invalid_space(self, invalid_term, error):

        # valid HANCESTRO terms with different delimiters - e.g. || (without spaces and commas) => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("unknown_valid_term", [HUMAN_VALID_VALUES[1]])  # unknown term
    @pytest.mark.parametrize("valid_term", [HUMAN_VALID_VALUES[0]])  # valid term
    @pytest.mark.parametrize("error", [ERROR_MESSAGE[2]])  # not a valid term error
    def test_ethnicity_valid_term_with_unknown(self, unknown_valid_term,valid_term, error):

        # unknown w/ valid HANCESTRO term => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = unknown_valid_term + " || " + valid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{unknown_valid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("test_term", [HUMAN_VALID_VALUES[3]])
    def test_ethnicity_add_label_order_valid_values(self, test_term):

        # add_labels check: obs.ethnicity order matches self_reported_ethnicity_ontology_term_id order

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"].cat.add_categories(test_term)
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "self_reported_ethnicity_ontology_term_id"] = test_term
        self.validator.validate_adata()
        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        map = {"European":"HANCESTRO:0005",
                "Hispanic or Latin American":"HANCESTRO:0014",
                "Greater Middle Eastern  (Middle Eastern or North African or Persian)":"HANCESTRO:0015"}
        assert labeler.adata.obs["self_reported_ethnicity_ontology_term_id"][0].split(" || ") == [map[t] for t in labeler.adata.obs["self_reported_ethnicity"][0].split(" || ")]


@pytest.mark.parametrize("test_h5ads", NON_HUMAN_H5ADS)
class TestNonHumanEthnicityOntologyValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas

    def test_nonhuman_passes(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("valid_term", [NON_HUMAN_VALID_VALUES[0]])
    def test_nonhuman_ethnicity_valid_values(self, valid_term):

        # various lists of valid HANCESTRO terms => valid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = valid_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []

    @pytest.mark.parametrize("invalid_term",[NON_HUMAN_INVALID_VALUES[0]])
    @pytest.mark.parametrize("error", [NON_HUMAN_ERROR_MESSAGE[0]])
    def test_nonhuman_ethnicity_invalid_human(self, invalid_term, error):

        # human HANCESTRO term in non-human dataset => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors
