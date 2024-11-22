"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1105
https://github.com/chanzuckerberg/single-cell-curation/pull/1119/
"""

import numpy as np
import pytest
from fixtures.valid_adatas import (
    validator_with_all_visiums,
    validator_with_spatial_adatas,
    validator_with_non_spatial_adata
)

def test_spatial_embeddings_nan_fails(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"].astype("float32")
    validator.adata.obsm["spatial"][0][0] = np.nan
    validator.validate_adata()
    assert not validator.is_valid
    assert "ERROR: adata.obs['spatial] contains at least one NaN value." in validator.errors


def test_non_visium_nan_embedding_passes(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    obsm_key = [k for k in validator.adata.obsm.keys()][0]
    validator.adata.obsm[obsm_key] = validator.adata.obsm[obsm_key].astype("float32")
    validator.adata.obsm[obsm_key][0][0] = np.nan
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_non_visium_all_nan_embedding_fails(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    obsm_key = [k for k in validator.adata.obsm.keys()][0]
    shape = validator.adata.obsm[obsm_key].shape
    empty_array = np.empty(shape=shape, dtype="float32")
    empty_array[:] = np.nan
    validator.adata.obsm[obsm_key] = empty_array
    validator.validate_adata()
    assert not validator.is_valid
    assert f"ERROR: adata.obsm['{obsm_key}'] contains all NaN values." in validator.errors


def test_spatial_all_nan_embedding_fails(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    shape = validator.adata.obsm["spatial"].shape
    empty_array = np.empty(shape=shape, dtype="float32")
    empty_array[:] = np.nan
    validator.adata.obsm["spatial"] = empty_array
    validator.validate_adata()
    assert not validator.is_valid
    assert f"ERROR: adata.obs['spatial] contains at least one NaN value." in validator.errors
