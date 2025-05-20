'''
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1355

Testing conditions:

- valid_integrated_visium.h5ad from: https://cellxgene.cziscience.com/e/83ec9e14-87a4-4d41-aa3b-f7c51af70a64.cxg/ with cxg-labels removed

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

INTEGRATED_VISIUM_H5ADS = [file for file in ALL_H5ADS if "integrated" in file]

@pytest.mark.parametrize("test_h5ads", INTEGRATED_VISIUM_H5ADS)
class TestVisiumData:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_V_valid(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("is_single_value", [False, np.bool_(False)])
    @pytest.mark.parametrize("is_primary_data_value", [False, np.bool_(False)])
    def test_valid_V_all_false(self,is_single_value, is_primary_data_value):

         # is_single: False, is_primary_data: False -> pass

        self.validator.adata.uns["spatial"]["is_single"] = is_single_value
        self.validator.adata.obs["is_primary_data"] = is_primary_data_value
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("is_single_value", [False, np.bool_(False)])
    @pytest.mark.parametrize("is_primary_data_value", [True, np.bool_(True)])
    def test_invalid_V_ipd_true(self, is_single_value, is_primary_data_value):

        # is_single: False, is_primary_data: True -> fail

        self.validator.adata.uns["spatial"]["is_single"] = is_single_value
        self.validator.adata.obs["is_primary_data"] = is_primary_data_value
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (f"ERROR: When uns['spatial']['is_single'] is {self.validator.adata.uns['spatial']['is_single']}, obs['is_primary_data'] must be False for all rows."
        ) in self.validator.errors
