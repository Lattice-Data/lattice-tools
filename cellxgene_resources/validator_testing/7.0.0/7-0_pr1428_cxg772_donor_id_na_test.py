"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1455

Testing conditions:

Should pass
(Y) - tissue_type != "cell line" with normal donor_ids
(Y) - tissue_type == "cell line" with na for donor_id
(Y) - mixed tissue_types with equal counts of donor_id != "na" and == "na"


Should not pass
(Y) - tissue_type == "cell line" with all normal donor_ids
(Y) - tissue_type == "cell line" with donor_id mix of one random 'na' and normal ids
(Y) - tissue_type != "cell line" with all na for donor_ids
(Y) - tissue_type != "cell line"  with donor_id mix of one random 'na' and normal ids
(Y) - mixed tissue_types with diiferent counts of donor_id != and == "na"
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
class TestDonorIDValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_non_cell_line_valid(self, tissue_type):

        # tissue_type != "cell line" with normal donor_ids

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")

        if tissue_type == "primary cell culture":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CL:0000617"

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_cell_line_donor_na_valid(self):

        # tissue_type == "cell line" with na for donor_id

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_tissue_type_mixed_valid(self):

        # mixed tissue_types with equal counts of donor_id != "na" and == "na"

        cell_line_values = {
            "tissue_type":"cell line",
            "tissue_ontology_term_id":"CVCL_2830",
            "donor_id":"na",
            "development_stage_ontology_term_id":"na",
            "sex_ontology_term_id":"na",
            "self_reported_ethnicity_ontology_term_id":"na"
        }

        original_donor_ids = self.validator.adata.obs["donor_id"]
        idx = self.validator.adata.obs.sample(frac=0.5, random_state=0).index

        # set tissue rows
        self.validator.adata.obs.loc[idx, "tissue_type"] = "tissue"
        self.validator.adata.obs.loc[idx, "donor_id"] = original_donor_ids[idx]

        # set cell line rows
        for k,v in cell_line_values.items():
            if v not in self.validator.adata.obs[k].unique():
                self.validator.adata.obs[k] = self.validator.adata.obs[k].cat.add_categories(v)

            self.validator.adata.obs.loc[self.validator.adata.obs.index.difference(idx), k] = v
            self.validator.adata.obs[k] = self.validator.adata.obs[k].astype("category")

        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_tissue_cell_line_with_normal_ids_invalid(self):

        # tissue_type == "cell line" with normal donor_ids

        # helper function sets donor_ids to na, reset to original values in fixture
        original_donor_ids = self.validator.adata.obs["donor_id"]
        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.adata.obs["donor_id"] = original_donor_ids
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        assert (
            "ERROR: Column 'donor_id' in dataframe 'obs' contains invalid values "
            f"['{self.validator.adata.obs['donor_id'].unique().tolist()}']. "
            "Values must be one of ['na'] when 'tissue_type' is 'cell line'.")


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cell_line_id_na_invalid(self, tissue_type):

        # tissue_type != "cell line" with donor_id == "na"

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")
        self.validator.adata.obs["donor_id"] = "na"
        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].astype("category")
        if tissue_type == "primary cell culture":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CL:0000617"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: Column 'donor_id' in dataframe 'obs' contains forbidden values '['na']'. Values must not be one of ['na']"
        ) in self.validator.errors


    def test_tissue_cell_line_with_normal_ids_and_na_invalid(self):

        # tissue_type == "cell line" and donor_id with one random na + normal ids

        original_donor_ids = self.validator.adata.obs["donor_id"]
        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.adata.obs["donor_id"] = original_donor_ids
        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].cat.add_categories(["na"])
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "donor_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        # trying assert for common error string ending instead of matching on full error
        # string for each dev term in test fixture
        for error in self.validator.errors:
            assert error.endswith(
                "Values must be one of ['na'] when 'tissue_type' is 'cell line'."
                )


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cellline_donor_id_one_na_invalid(self, tissue_type):

        # tissue_type != cell line and donor_id with one random 'na' + normal ids

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")

        if tissue_type == "primary cell culture":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CL:0000617"

        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].cat.add_categories(["na"])
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "donor_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: Column 'donor_id' in dataframe 'obs' contains forbidden values '['na']'. Values must not be one of ['na']"
        ) in self.validator.errors


    def test_tissue_type_mixed_donor_unequal_invalid(self):

        # mixed tissue_types with diiferent counts of donor_id != and == "na"
        # this test is not using utils.py to side step conflicts with other tests
        # I may have overcomplicated this test

        cell_line_values = {
            "tissue_type":"cell line",
            "tissue_ontology_term_id":"CVCL_2830",
            "development_stage_ontology_term_id":"na",
            "sex_ontology_term_id":"na",
            "self_reported_ethnicity_ontology_term_id":"na"
        }

        original_donor_ids = self.validator.adata.obs["donor_id"]
        idx_tissue = self.validator.adata.obs.sample(frac=0.5, random_state=0).index
        idx_tissue_donor_id = self.validator.adata.obs.sample(frac=0.5, random_state=1).index

        # set tissue rows
        self.validator.adata.obs.loc[idx_tissue, "tissue_type"] = "tissue"
        self.validator.adata.obs.loc[idx_tissue_donor_id, "donor_id"] = original_donor_ids[idx_tissue_donor_id]

        # set cell line rows
        idx_cl = self.validator.adata.obs.sample(frac=0.5, random_state=2).index
        idx_cl_donor_id = self.validator.adata.obs.sample(frac=0.5, random_state=3).index

        for k,v in cell_line_values.items():
            if v not in self.validator.adata.obs[k].unique():
                self.validator.adata.obs[k] = self.validator.adata.obs[k].cat.add_categories(v)

            self.validator.adata.obs.loc[idx_cl, k] = v
            self.validator.adata.obs[k] = self.validator.adata.obs[k].astype("category")

        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].cat.add_categories("na")
        self.validator.adata.obs.loc[idx_cl_donor_id, "donor_id"] = "na"
        self.validator.adata.obs["donor_id"] = self.validator.adata.obs["donor_id"].astype("category")
        unique_donors_for_error = [x for x in self.validator.adata.obs.loc[self.validator.adata.obs["tissue_type"] == "cell line"]["donor_id"].unique().tolist() if x != "na"]
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Column 'donor_id' in dataframe 'obs' contains invalid values '{unique_donors_for_error}'. "
            "Values must be one of ['na'] when 'tissue_type' is 'cell line'."
        ) in self.validator.errors
        assert (
            "ERROR: Column 'donor_id' in dataframe 'obs' contains forbidden values '['na']'. Values must not be one of ['na']"
        ) in self.validator.errors