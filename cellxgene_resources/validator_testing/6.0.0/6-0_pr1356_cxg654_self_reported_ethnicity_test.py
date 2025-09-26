"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1303
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1356

Should pass:
(y) various lists of valid HANCESTRO terms
    (y) just 1 valid HANCESTRO term
    (y) just unknown term
    (y) multiple valid terms in lexical order
(y) add labels check: obs.self_reported_ethnicity order matches  self_reported_ethnicity_ontology_term_id order

Should not pass:
(y) valid HANCESTRO terms out of lexical order
(y) valid HANCESTRO terms with different delimiters - e.g. commas
(y) valid HANCESTRO terms with different delimiters - e.g. || (without spaces)
(y) valid terms duplicated
(y) unknown term w/ valid HANCESTRO term
(y) invalid HANCESTRO term + valid HANCESTRO term
(y) forbidden term
(y) "na" for human

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

HUMAN_H5ADS = [f for f in ALL_H5ADS if any(s in f for s in ["human","integrated"])]
NON_HUMAN_H5ADS = [f for f in ALL_H5ADS if not any(s in f for s in ["human","integrated"])]

HUMAN_VALID_VALUES = {
    "single valid term":"HANCESTRO:0861",
    "unknown":"unknown",
    "multiple valid terms":"HANCESTRO:0847 || HANCESTRO:0861 || HANCESTRO:0862"
}

HUMAN_INVALID_VALUES = {
    "out of lexical order":"HANCESTRO:0861 || HANCESTRO:0847",
    "different delimiter":"HANCESTRO:0847,HANCESTRO:0861",
    "no spaces":"HANCESTRO:0847||HANCESTRO:0861",
    "duplicate terms":"HANCESTRO:0847 || HANCESTRO:0847",
    "valid term with whitespace typo":"HANCESTRO:0847 ",
    "forbidden term":"HANCESTRO:0031",
    "na":"na"
}

error_message_suffix = " When 'organism_ontology_term_id' is 'NCBITaxon:9606' (Homo sapiens) and 'tissue_type' is not 'cell line', self_reported_ethnicity_ontology_term_id MUST be formatted as one or more AfPO or HANCESTRO terms that are descendants of 'HANCESTRO:0601' for ethnicity category or 'HANCESTRO:0602' for geography-based population category, in ascending lexical order with the delimiter ` || `, or 'unknown' if unavailable."

ERROR_MESSAGE = {
    "contains duplicates":"in 'self_reported_ethnicity_ontology_term_id' contains duplicates." + error_message_suffix,
    "not in lexical order":"in 'self_reported_ethnicity_ontology_term_id' is not in ascending lexical order." + error_message_suffix,
    "not a valid term":"in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id of 'HANCESTRO, AFPO'." + error_message_suffix,
    "forbidden":"in 'self_reported_ethnicity_ontology_term_id' is a deprecated term id of 'HANCESTRO'." + error_message_suffix
}

NON_HUMAN_VALID_VALUES = {
    "na":"na"
    }

NON_HUMAN_INVALID_VALUES = {
    "nonhuman invalid term":"HANCESTRO:0005"
    }

NON_HUMAN_ERROR_MESSAGE = {
    "not a valid term":"in 'self_reported_ethnicity_ontology_term_id' is not a valid value of 'self_reported_ethnicity_ontology_term_id'. When 'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'."
    }


@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
class TestHumanEthnicityOntologyValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_passes(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("valid_term", HUMAN_VALID_VALUES.values())
    def test_ethnicity_valid_values(self, valid_term):

        # various lists of valid HANCESTRO terms => valid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = valid_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["out of lexical order"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not in lexical order"]])
    def test_ethnicity_invalid_order(self, invalid_term, error):

        # valid HANCESTRO terms out of lexical order => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["different delimiter"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_invalid_delimiter(self, invalid_term, error):

        # valid HANCESTRO terms with different delimiters - e.g. commas => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["no spaces"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_invalid_space(self, invalid_term, error):

        # valid HANCESTRO terms with different delimiters - e.g. || (without spaces) => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["duplicate terms"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["contains duplicates"]])
    def test_ethnicity_invalid_duplicates(self, invalid_term, error):

        # valid HANCESTRO terms with duplicates => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        print(self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"])
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors


    @pytest.mark.parametrize("unknown_valid_term", [HUMAN_VALID_VALUES["unknown"]])
    @pytest.mark.parametrize("valid_term", [HUMAN_VALID_VALUES["single valid term"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_valid_term_with_unknown(self, unknown_valid_term,valid_term, error):

        # unknown w/ valid HANCESTRO term => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = unknown_valid_term + " || " + valid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{unknown_valid_term}' {error}"
            ) in self.validator.errors

    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["na"]])
    @pytest.mark.parametrize("valid_term", [HUMAN_VALID_VALUES["single valid term"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_invalid_term_with_valid(self, invalid_term,valid_term, error):

        # invalid w/ valid HANCESTRO term => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term + " || " + valid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["valid term with whitespace typo"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_whitespace_for_human(self, invalid_term, error):

        # whitespace in self_reported_ethnicity_ontology_term_id => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["forbidden term"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["forbidden"]])
    def test_ethnicity_forbidden_term(self, invalid_term, error):

        # forbidden in self_reported_ethnicity_ontology_term_id => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("invalid_term", [HUMAN_INVALID_VALUES["na"]])
    @pytest.mark.parametrize("error", [ERROR_MESSAGE["not a valid term"]])
    def test_ethnicity_na_for_human(self, invalid_term, error):

        # na in self_reported_ethnicity_ontology_term_id => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: '{invalid_term}' {error}"
            ) in self.validator.errors


    @pytest.mark.parametrize("test_term", [HUMAN_VALID_VALUES["multiple valid terms"]])
    def test_ethnicity_add_label_order_valid_values(self, test_term):

        # add_labels check: obs.ethnicity order matches self_reported_ethnicity_ontology_term_id order

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"].cat.add_categories(test_term)
        self.validator.adata.obs.loc[self.validator.adata.obs.index[0], "self_reported_ethnicity_ontology_term_id"] = test_term
        self.validator.validate_adata()
        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        map = {"Asian":"HANCESTRO:0847",
                "Black British":"HANCESTRO:0861",
                "Asian-American":"HANCESTRO:0862"}
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


    @pytest.mark.parametrize("valid_term", NON_HUMAN_VALID_VALUES.values())
    def test_nonhuman_ethnicity_valid_values(self, valid_term):

        # 'na' in self_reported_ethnicity_ontology_term_id => valid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = valid_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []

    @pytest.mark.parametrize("invalid_term",[NON_HUMAN_INVALID_VALUES["nonhuman invalid term"]])
    @pytest.mark.parametrize("error", [NON_HUMAN_ERROR_MESSAGE["not a valid term"]])
    def test_nonhuman_ethnicity_invalid_human(self, invalid_term, error):

        # human HANCESTRO term in non-human dataset => invalid

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = invalid_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert f"ERROR: '{invalid_term}' {error}" in self.validator.errors
