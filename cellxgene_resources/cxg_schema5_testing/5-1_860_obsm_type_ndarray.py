"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/860
https://github.com/chanzuckerberg/single-cell-curation/issues/853
https://github.com/chanzuckerberg/single-cell-curation/issues/855
https://github.com/chanzuckerberg/single-cell-curation/pull/859
"""

import pytest
import numpy as np
import pandas as pd
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


class TestAllObsmTypeNdarray:
    def test_obsm_no_changes(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_obsm_obs_to_obsm(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["test"] = validator.adata.obs
        obsm = validator.adata.obsm["test"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['test']' is {type(obsm)}')."
        ]


# anndata uses shape attribute when setting obsm values. as far as I can tell, only numpy ndarrays and pandas dfs have
# shape attribute
@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestSpatialObsmTypeNdarray:
    def test_obsm_value_as_df(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["spatial"] = pd.DataFrame(validator.adata.obsm["spatial"], index=validator.adata.obs.index)
        obsm = validator.adata.obsm["spatial"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['spatial']' is {type(obsm)}')."
        ]

    # obsm X_ and spatial need at least 2 columns
    def test_obsm_shape_spatial(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:,:1]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['spatial']' has columns='1'."
        ]


# obsm X_ and spatial need at least 2 columns
@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestNonSpatialObsmTypeNdArray:
    def test_obsm_shape_non_spatial(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["X_umap"] = validator.adata.obsm["X_umap"][:,:1]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['X_umap']' has columns='1'."
        ]

    # other obsm can have only 1 column
    def test_obsm_shape_one_column(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["test"] = np.zeros((validator.adata.n_obs, 1))
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_obsm_shape_zero_columns(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["test"] = np.zeros((validator.adata.n_obs, 0))
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: The size of the ndarray stored for a 'adata.obsm['test']' MUST NOT be zero.", 
            "ERROR: All unspecified embeddings must have at least one column. 'adata.obsm['test']' has columns='0'."
        ]
