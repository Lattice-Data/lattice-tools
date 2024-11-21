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
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"


def test_obsm_no_changes(validator_with_all_adatas):
    validator = validator_with_all_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_obsm_obs_to_obsm(validator_with_all_adatas):
    validator = validator_with_all_adatas
    validator.adata.obsm["test"] = validator.adata.obs
    obsm = validator.adata.obsm["test"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['test']' is {type(obsm)}')."
    ]


# anndata uses shape attribute when setting obsm values. as far as I can tell, only numpy ndarrays and pandas dfs have
# shape attribute
def test_obsm_value_as_df(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obsm["spatial"] = pd.DataFrame(validator.adata.obsm["spatial"], index=validator.adata.obs.index)
    obsm = validator.adata.obsm["spatial"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: All embeddings have to be of 'numpy.ndarray' type, 'adata.obsm['spatial']' is {type(obsm)}')."
    ]


# obsm X_ and spatial need at least 2 columns
def test_obsm_shape_spatial(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:,:1]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['spatial']' has columns='1'."
    ]


# obsm X_ and spatial need at least 2 columns
def test_obsm_shape_non_spatial(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obsm["X_umap"] = validator.adata.obsm["X_umap"][:,:1]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. 'adata.obsm['X_umap']' has columns='1'."
    ]


# other obsm can have only 1 column
def test_obsm_shape_one_column(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obsm["test"] = np.zeros((validator.adata.n_obs, 1))
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_obsm_shape_zero_columns(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obsm["test"] = np.zeros((validator.adata.n_obs, 0))
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: The size of the ndarray stored for a 'adata.obsm['test']' MUST NOT be zero.", 
        "ERROR: All unspecified embeddings must have at least one column. 'adata.obsm['test']' has columns='0'."
    ]
