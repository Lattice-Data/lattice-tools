"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1444

Testing conditions:
Should pass
() - tissue_type == "cell line" + sex term == “na”
() - tissue_type == "cell line" + sex term == “na” +  sex label == “na”
() - all tissue_types + normal sex terms

Should not pass
() - tissue_type == "cell line" + normal sex terms
() - tissue_type != "cell line" + sex term == “na”
() - tissue_type == "cell line" + sex term == mix of “na” and “unknown”
() - all tissue_types + sex terms == mix of one random 'na' and valid sex terms

"""

import numpy as np
import pytest
from cellxgene_schema.write_labels import AnnDataLabelAppender
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
)


NON_CELL_LINE_TISSUES = [
    "organoid",
    "primary cell culture",
    "tissue",
]
ALL_TISSUE_TYPES = [
    "organoid",
    "primary cell culture",
    "tissue",
    "cell line",
]


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestDevStageValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_tissue_cell_line_sex_na_valid(self):

        # tissue_type == "cell line" + sex term == “na”

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"  # need dev stage set to na to be valid
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_label_na_for_cell_line_sex_na_valid(self):

        # tissue_type == "cell line" with na for sex term has na for sex label

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"  # need dev stage set to na to be valid
        self.validator.validate_adata()
        assert self.validator.is_valid

        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        assert labeler.adata.obs["sex"].unique()[0] == "na"


    def test_tissue_all_tissue_type_normal_sex_terms_valid(self,ALL_TISSUE_TYPES):

        # all tissue_type + normal sex terms

        self.validator.adata.obs["tissue_type"] = ALL_TISSUE_TYPES
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_normal_sex_term_invalid(self):

        # tissue_type == "cell line" + normal sex terms

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"  # testing sex term, need dev stage set to na
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each sex term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'sex_ontology_term_id' MUST be 'na'.")



    def test_tissue_cell_line_with_unknown_invalid(self):

        # tissue_type == "cell line" and sex term unknown

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        self.validator.adata.obs["sex_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'sex_ontology_term_id' MUST be 'na'.")


    def test_tissue_cell_line_with_unknown_and_na_invalid(self):

        # tissue_type == "cell line" and dev stage mix of na and unknown
        # UBERON dev stage organisms now passing

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "sex_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'sex_ontology_term_id' MUST be 'na'.")


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_sex_term_na_invalid(self, tissue_type):

        # tissue_type != "cell line" with na for dev term
        # primary cell culture has further rules for tissue terms, need to take account of those

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'sex_ontology_term_id' is not allowed. "
            "When 'tissue_type' is not 'cell line', 'sex_ontology_term_id' cannot be 'na'."
        ) in self.validator.errors


    @pytest.mark.parametrize("tissue_type", ALL_TISSUE_TYPES)
    def test_tissue_not_cell_line_dev_term_one_na_fails(self, tissue_type):

        # all tissue_types with one random 'na' mixed in with valid dev terms
        # primary cell culture has further rules for tissue terms, need to take account of those

        self.validator.adata.obs["tissue_type"] = tissue_type
        if tissue_type == "cell line":
            self.validator.adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["sex_ontology_term_id"] = self.validator.adata.obs["sex_ontology_term_id"].cat.add_categories("na")
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "sex_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid

        # this works, but probably should be 2 tests to have less complicated logic for asserts
        if tissue_type != "cell line":
            assert (
                "ERROR: 'na' in 'sex_ontology_term_id' is not allowed. "
                "When 'tissue_type' is not 'cell line', 'sex_ontology_term_id' cannot be 'na'."
            ) in self.validator.errors
        else:
            for error in self.validator.errors:
                assert error.endswith("When 'tissue_type' is 'cell line', 'sex_ontology_term_id' MUST be 'na'.")

