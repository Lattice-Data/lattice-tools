"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/828
https://github.com/chanzuckerberg/single-cell-curation/pull/876
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from dask.array import from_array
from scipy import sparse
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]




@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestVisiumZeroValues:
    def test_visiums_pass(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # only raw layer is checked
    def test_visiums_with_normalized_all_zero_row(self, validator_with_adatas):
        validator = validator_with_adatas
        normalized = np.log1p(validator.adata.X, dtype=np.float32)
        raw_adata = ad.AnnData(validator.adata.X, dtype=np.float32, var=validator.adata.var, obs=validator.adata.obs)
        index_name = validator.adata.obs.index[100]
        validator.adata.obs.loc[index_name, "in_tissue"] = 0
        normalized[100] = 0.0
        validator.adata.X = normalized
        validator.adata.raw = raw_adata
        del validator.adata.raw.var["feature_is_filtered"]
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # only raw layer is checked
    def test_visiums_with_normalized_all_zero_row_in_tissue_1(self, validator_with_adatas):
        validator = validator_with_adatas
        normalized = np.log1p(validator.adata.X)
        raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
        normalized[100] = 0.0     # this row in_tissue == 0
        validator.adata.X = normalized
        validator.adata.raw = raw_adata
        validator.adata.obs.loc[validator.adata.obs.index[100], "in_tissue"] = 1
        del validator.adata.raw.var["feature_is_filtered"]
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_visiums_with_normalized_and_raw_all_zero_row(self, validator_with_adatas):
        validator = validator_with_adatas
        normalized = np.log1p(validator.adata.X, dtype=np.float32)
        raw_adata = ad.AnnData(validator.adata.X, dtype=np.float32, var=validator.adata.var, obs=validator.adata.obs)
        validator.adata.obs.loc[validator.adata.obs.index[100], "in_tissue"] = 0
        normalized[100] = 0.0     # this row in_tissue == 0
        raw_adata.X[100] = 0.0
        validator.adata.X = normalized
        validator.adata.raw = raw_adata
        del validator.adata.raw.var["feature_is_filtered"]
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_visiums_with_normalized_and_raw_all_zero_row_not_valid(self, validator_with_adatas):
        validator = validator_with_adatas
        normalized = np.log1p(validator.adata.X)
        raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
        normalized[100] = 0.0     # this row in_tissue == 0
        raw_adata.X[100] = 0.0
        validator.adata.obs.loc[validator.adata.obs.index[100], "in_tissue"] = 1
        validator.adata.X = normalized
        validator.adata.raw = raw_adata
        del validator.adata.raw.var["feature_is_filtered"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: Each observation with obs['in_tissue'] == 1 must have at least one non-zero value in its row in the raw matrix.", 
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements."
        ]

    def test_in_tissue_one(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[2895], "in_tissue"] = 1
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: Each observation with obs['in_tissue'] == 1 must have at least one non-zero value in its row in the raw matrix.", 
            "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
        ]

    def test_make_all_zero_row(self, validator_with_adatas):
        validator = validator_with_adatas
        index_name = validator.adata.obs.index[100]
        validator.adata.obs.loc[index_name, "in_tissue"] = 0
        matrix = validator.adata.X.compute()
        matrix[100] = 0.0
        validator.adata.X = from_array(sparse.csr_matrix(matrix, dtype=np.float32))
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_is_single_false(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        validator.adata.obs["is_primary_data"] = False
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]
        del validator.adata.obs["array_col"]
        del validator.adata.obs["array_row"]
        del validator.adata.obs["in_tissue"]
        validator.validate_adata()
        assert validator.is_valid is False
        # expected_errors = [
        #     "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
        #     "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
        # ]
        # for error in expected_errors:
        #     assert error in validator.errors
        assert len(validator.errors) == 2


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestNonVisiumMatrixAllZero:
    def test_make_all_zero_row_raw_non_visium(self, validator_with_adatas):
        validator = validator_with_adatas
        matrix = validator.adata.raw.X.compute()
        matrix = matrix.toarray()
        matrix[0] = 0
        raw_adata = ad.AnnData(
            X=from_array(sparse.csr_matrix(matrix, dtype=np.float32)),
            var=validator.adata.raw.var, 
            obs=validator.adata.obs
        )
        validator.adata.raw = raw_adata
        validator.validate_adata()
        assert validator.is_valid is False
        expected_errors = [
            "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
            "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements."
        ]
        for error in expected_errors:
            assert error in validator.errors

    def test_make_all_zero_row_non_visium(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.X = validator.adata.raw.X
        del validator.adata.raw
        validator.adata.var["feature_is_filtered"] = False
        matrix = validator.adata.X.compute()
        matrix = matrix.toarray()
        matrix[0] = 0
        validator.adata.X = from_array(sparse.csr_matrix(matrix, dtype=np.float32))
        validator.validate_adata()
        assert validator.is_valid is False
        expected_errors = [
            "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
            "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
        ]
        for error in expected_errors:
            assert error in validator.errors
