"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1438

Testing conditions:
Should not pass
() - tissue_type != "cell line" with na for cell term
() - tissue_type == "cell line" and cell term mix of na and unknown
(N) - all tissue_types with one random 'na' mixed in with valid cell terms

Should pass
(Y) - tissue_type == "cell line" with na for all cell terms
(Y) - tissue_type == "cell line" with normal cell terms
(Y) - tissue_type == "cell line" and cell term all unknown - SHOULD THIS BE ALLOWED?
(Y) - tissue_type == "cell line" with na for cell term has na for cell label
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
class TestCellTermValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_tissue_cell_line_cell_term_na_passes(self):

        # tissue_type == "cell line" with na for cell term
        # visium is_single == True will fail since cell_type should be unknown for in_tissue == 0
        # i think this is ok, visium slices will almost always be tissue/organoid sections
        # not sure if we should have tissue_type restrictions for visium is_single

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_valid_cell_terms_passes(self):

        # tissue_type == "cell line" with normal cell terms

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_all_unknown_passes(self):

        # tissue_type == "cell line" and cell term all unknown
        # do we want this behavior or instead us na?
        # some visium fixtures fail on mismatched uns color column after setting cell term to unknown

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_cell_term_with_unknown_and_na_fails(self):

        # tissue_type == "cell line" and cell term mix of na and unknown
        # currently pass for non-visium datasets
        # do we need to test split datasets with mix of "cell line" and other tissue_types?

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "cell_type_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_cell_term_na_fails(self, tissue_type):

        # tissue_type != "cell line" with na for cell term
        # should error message be more specific?

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'cell_type_ontology_term_id' is not a valid ontology term id of 'CL, ZFA, FBbt, WBbt'."
        ) in self.validator.errors


    @pytest.mark.parametrize("tissue_type", ALL_TISSUE_TYPES)
    def test_tissue_not_cell_line_cell_term_one_na_fails(self, tissue_type):

        # all tissue_types with one random 'na' mixed in with valid cell terms
        # tissue_type == "cell line" fails for missing dev term == na, not due to just one na value

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")

        self.validator.adata.obs["cell_type_ontology_term_id"] = self.validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories("na")
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'cell_type_ontology_term_id' is not a valid ontology term id of 'CL, ZFA, FBbt, WBbt'."
        ) in self.validator.errors


    def test_label_na_for_cell_line_cell_type_na(self):

        # tissue_type == "cell line" with na for cell term has na for cell label

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid

        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        assert labeler.adata.obs["cell_type"].unique()[0] == "na"
