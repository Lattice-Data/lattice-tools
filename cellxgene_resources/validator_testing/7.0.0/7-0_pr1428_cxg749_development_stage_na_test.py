"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1428

Testing conditions:
Should not pass
(Y) - tissue_type == "cell line" with normal dev terms
(Y) - tissue_type != "cell line" with na for dev term
(Y) - tissue_type == "cell line" and dev stage unknown
(Y) - tissue_type == "cell line" and dev stage mix of na and unknown
(Y) - all tissue_types with one random 'na' mixed in with valid dev terms

Should pass
(Y) - tissue_type == "cell line" with na for dev term
(Y) - tissue_type == "cell line" with na for dev term has na for dev label
"""

import numpy as np
import pytest
from cellxgene_schema.write_labels import AnnDataLabelAppender
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from testing_internals.utils import (
    make_valid_cell_line_fixture,
    NA_COLUMNS
)
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


    def test_tissue_cell_line_dev_na_passes(self):

        # tissue_type == "cell line" with na for dev term

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_valid_dev_term_fails(self):

        # tissue_type == "cell line" with normal dev terms

        # helper function sets dev terms to na, reset to original values in fixture
        original_dev_terms = self.validator.adata.obs["development_stage_ontology_term_id"]
        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.adata.obs["development_stage_ontology_term_id"] = original_dev_terms

        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'development_stage_ontology_term_id' MUST be 'na'.")


    def test_tissue_cell_line_with_unknown_fails(self):

        # tissue_type == "cell line" and dev stage unknown

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.adata.obs["development_stage_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'development_stage_ontology_term_id' MUST be 'na'.")


    def test_tissue_cell_line_with_unknown_and_na_fails(self):

        # tissue_type == "cell line" and dev stage mix of na and unknown
        # UBERON dev stage organisms now passing

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "development_stage_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith("When 'tissue_type' is 'cell line', 'development_stage_ontology_term_id' MUST be 'na'.")


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_dev_term_na_fails(self, tissue_type):

        # tissue_type != "cell line" with na for dev term

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'development_stage_ontology_term_id' is not allowed. "
            "When 'tissue_type' is not 'cell line', 'development_stage_ontology_term_id' cannot be 'na'." 
        ) in self.validator.errors


    @pytest.mark.parametrize("tissue_type", ALL_TISSUE_TYPES)
    def test_tissue_not_cell_line_dev_term_one_na_fails(self, tissue_type):

        # all tissue_types with one random 'na' mixed in with valid dev terms

        self.validator.adata.obs["tissue_type"] = tissue_type
        if tissue_type == "cell line":
            self.validator.adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'

        # errors for each unique value in ontologies, generate list of tuples
        # for error string matching
        ontology_values = []
        for ontology in NA_COLUMNS:
            for value in self.validator.adata.obs[ontology].unique():
                ontology_values.append((ontology, value))

        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        for col in NA_COLUMNS:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].cat.add_categories("na")
            self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], col] = "na"

        self.validator.validate_adata()
        assert not self.validator.is_valid

        # this works, but probably should be 2 tests to have less complicated logic for asserts
        if tissue_type != "cell line":
            assert (
                "ERROR: 'na' in 'development_stage_ontology_term_id' is not allowed. "
                "When 'tissue_type' is not 'cell line', 'development_stage_ontology_term_id' cannot be 'na'." 
            ) in self.validator.errors
        else:
            for ontology, value in ontology_values:
                error_string = (
                    f"ERROR: '{value}' in '{ontology}' is not a valid value of '{ontology}'. "
                    f"When 'tissue_type' is 'cell line', '{ontology}' MUST be 'na'."
                )
                assert error_string in self.validator.errors


    def test_label_na_for_cell_line_dev_stage_na(self):

        # tissue_type == "cell line" with na for dev term has na for dev label
        # all organisms now pass

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid

        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        assert labeler.adata.obs["development_stage"].unique()[0] == "na"
