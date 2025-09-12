"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1438

Testing conditions:
Should not pass
(Y) - tissue_type != "cell line" with na for cell term
(Y) - tissue_type == "cell line" and cell term mix of na and unknown
(N) - all tissue_types with one random 'na' mixed in with valid cell terms

Should pass
(Y) - tissue_type == "cell line" with na for all cell terms
(Y) - tissue_type != "cell line" 
(Y) - tissue_type == "cell line" with normal cell terms
(Y) - tissue_type == "cell line" and cell term all unknown - SHOULD THIS BE ALLOWED?
(Y) - tissue_type == "cell line" with na for cell term has na for cell label
(Y) - tissue_type == "cell line" with all na for cell term, other tissue_type with normal CL
"""

import numpy as np
import pytest
from cellxgene_schema.write_labels import AnnDataLabelAppender
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from fixtures.valid_adatas import (
    NON_SPATIAL_H5ADS,
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


# use non-spatial to reduce errors for visium in_tissue == 0
@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestCellTermValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_tissue_cell_line_cell_term_na_passes(self):

        # tissue_type == "cell line" with na for cell term

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_passes(self, tissue_type):

        # tissue_type != "cell line" passes

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")

        # set tissue term to random CL term to pass validation
        if tissue_type == "primary cell culture":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CL:0000617"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_valid_cell_terms_passes(self):

        # tissue_type == "cell line" with normal cell terms

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_with_all_unknown_passes(self):

        # tissue_type == "cell line" and cell term all unknown

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "unknown"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_cell_term_with_unknown_and_na_fails(self):

        # tissue_type == "cell line" and cell term mix of na and unknown

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "cell_type_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_cell_term_na_fails(self, tissue_type):

        # tissue_type != "cell line" with na for cell term

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        if tissue_type == "cell line":
            assert (
                "ERROR: When tissue_type is 'cell line', 'na' is allowed for 'cell_type_ontology_term_id' "
                "but then all observations where tissue_type is 'cell line' MUST be 'na'."
            ) in self.validator.errors
        else:
            assert (
                "ERROR: 'na' in 'cell_type_ontology_term_id' is not a valid ontology term id of "
                "'CL, ZFA, FBbt, WBbt'. cell_type_ontology_term_id must be a valid ontology term of "
                "CL or 'unknown'. WBbt, ZFA, and FBbt terms can be allowed depending on the organism. "
                "'na' is allowed if tissue_type is 'cell line'."
            ) in self.validator.errors


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_mix_of_tissue_types_cell_line_na_passes(self, tissue_type):

        # tissue_type == "cell line" with all na for cell term, other tissue_type with normal CL

        # clear categories for modified cols, re-categorize before validating
        cols = [
            "tissue_type",
            "cell_type_ontology_term_id",
            "development_stage_ontology_term_id",
            "sex_ontology_term_id",
            "tissue_ontology_term_id",
        ]

        half_point = self.validator.adata.shape[0] // 2
        half_index = self.validator.adata.obs.iloc[half_point].name
        for col in cols:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].astype("object")

        # set first half to not cell line
        self.validator.adata.obs.loc[:half_index, "tissue_type"] = tissue_type

        # set tissue term to random CL term to pass validation
        if tissue_type == "primary cell culture":
            self.validator.adata.obs.loc[:half_index, "tissue_ontology_term_id"] = "CL:0000617"

        # set second half to cell line and cell type na
        self.validator.adata.obs.loc[half_index:, "tissue_type"] = "cell line"
        self.validator.adata.obs.loc[half_index:, "tissue_ontology_term_id"] = "CVCL_2830"
        self.validator.adata.obs.loc[half_index:, "cell_type_ontology_term_id"] = "na"
        self.validator.adata.obs.loc[half_index:, "development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs.loc[half_index:, "sex_ontology_term_id"] = "na"

        for col in cols:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].astype("category")

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("tissue_type", ALL_TISSUE_TYPES)
    def test_tissue_not_cell_line_cell_term_one_na_fails(self, tissue_type):

        # all tissue_types with one random 'na' mixed in with valid cell terms
        # currently one unexpected test result for tissue_type == "cell line" in valid_mouse.h5ad

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        if tissue_type == "cell line":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
            self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
            self.validator.adata.obs["sex_ontology_term_id"] = "na"

        self.validator.adata.obs["cell_type_ontology_term_id"] = self.validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories("na")
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "cell_type_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        if tissue_type == "cell line":
            assert (
                "ERROR: When tissue_type is 'cell line', 'na' is allowed for 'cell_type_ontology_term_id' "
                "but then all observations where tissue_type is 'cell line' MUST be 'na'."
            ) in self.validator.errors
        else:
            assert (
                "ERROR: 'na' in 'cell_type_ontology_term_id' is not a valid ontology term id of "
                "'CL, ZFA, FBbt, WBbt'. cell_type_ontology_term_id must be a valid ontology term of "
                "CL or 'unknown'. WBbt, ZFA, and FBbt terms can be allowed depending on the organism. "
                "'na' is allowed if tissue_type is 'cell line'."
            ) in self.validator.errors


    def test_label_na_for_cell_line_cell_type_na(self):

        # tissue_type == "cell line" with na for cell term has na for cell label

        self.validator.adata.obs["tissue_type"] = "cell line"
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["tissue_ontology_term_id"] = "CVCL_2830"
        # need dev stage set to na to be valid
        self.validator.adata.obs["development_stage_ontology_term_id"] = "na"
        self.validator.adata.obs["cell_type_ontology_term_id"] = "na"
        self.validator.adata.obs["sex_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid

        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        assert labeler.adata.obs["cell_type"].unique()[0] == "na"
