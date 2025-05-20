'''
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1355

Testing conditions for is_single = False (non-spatial and integrated visium datasets):
- valid_human.h5ad converted to integrated visium (is_single = False)
- wild integrated visium datasets: https://cellxgene.cziscience.com/e/22ec5e4c-a2db-4769-8c4a-9530b5f6962d.cxg/ and https://cellxgene.cziscience.com/e/83ec9e14-87a4-4d41-aa3b-f7c51af70a64.cxg/

##################################
-> ISSUE with submission of integrated visium datasets: are all dots included in integrated dataset?
If so, edge dots with the low duplicate counts or all zero counts will be included, and then fail the duplicate check,
or are they filtering for only meaningful dots and the duplicate check should pass?

    Currently:
    - valid_human.h5ad doesn't include duplicates
    - wild integrated visium datasets don't include duplicates
##################################

Should pass:
    - is_single: False (both bool and npbool), is_primary_data: False

Should not pass:
    - is_single: False (both bool and npbool), is_primary_data: True

'''

import pytest
import anndata as ad
import numpy as np
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

NON_SPATIAL_HUMAN = [file for file in ALL_H5ADS if "human" in file and "visium" not in file and "slide_seq" not in file]
INTEGRATED_VISIUM_H5ADS = [file for file in ALL_H5ADS if "integrated" in file]


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_HUMAN)
class TestNonVisiumData:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_valid(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_nonV_all_false(self):

        # is_single: False (bool), is_primary_data: False -> pass

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0022860"
        self.validator.adata.uns["spatial"] = {'is_single': False}
        self.validator.adata.obs["is_primary_data"] = False
        self.validator.adata.obs["suspension_type"] = "na"
        self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_nonV_all_false_npbool(self):

        # is_single: False (npbool), is_primary_data: False -> pass

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0022860"
        self.validator.adata.uns["spatial"] = {'is_single': np.bool_(False)}
        self.validator.adata.obs["is_primary_data"] = False
        self.validator.adata.obs["suspension_type"] = "na"
        self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_nonV_ipd_true(self):

        # is_single: False (bool), is_primary_data: True -> fail

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0022860"
        self.validator.adata.uns["spatial"] = {'is_single': False}
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.adata.obs["suspension_type"] = "na"
        self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: When uns['spatial']['is_single'] is {self.validator.adata.uns['spatial']['is_single']}, obs['is_primary_data'] must be False for all rows."
        ) in self.validator.errors


    def test_nonV_ipd_true_npbool(self):

        # is_single: False (npbool), is_primary_data: True -> fail

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0022860"
        self.validator.adata.uns["spatial"] = {'is_single': np.bool_(False)}
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.adata.obs["suspension_type"] = "na"
        self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: When uns['spatial']['is_single'] is {self.validator.adata.uns['spatial']['is_single']}, obs['is_primary_data'] must be False for all rows."
        ) in self.validator.errors


@pytest.mark.parametrize("test_h5ads", INTEGRATED_VISIUM_H5ADS)
class TestVisiumData:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_V_valid(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_valid_V_all_false(self):

         # is_single: False (bool), is_primary_data: False -> pass

        self.validator.adata.uns["spatial"]["is_single"] = False
        self.validator.adata.obs["is_primary_data"] = False
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_valid_V_all_false_npbool(self):

        # is_single: False (npbool), is_primary_data: False -> pass

        self.validator.adata.uns["spatial"]["is_single"] = np.bool_(False)
        self.validator.adata.obs["is_primary_data"] = False
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_invalid_V_ipd_true(self):

        # is_single: False (bool), is_primary_data: True -> fail

        self.validator.adata.uns["spatial"]["is_single"] = False
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: When uns['spatial']['is_single'] is {self.validator.adata.uns['spatial']['is_single']}, obs['is_primary_data'] must be False for all rows."
        ) in self.validator.errors


    def test_invalid_V_ipd_true_npbool(self):

        # is_single: False (npbool), is_primary_data: True -> fail

        self.validator.adata.uns["spatial"]["is_single"] = np.bool_(False)
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: When uns['spatial']['is_single'] is {self.validator.adata.uns['spatial']['is_single']}, obs['is_primary_data'] must be False for all rows."
    ) in self.validator.errors
