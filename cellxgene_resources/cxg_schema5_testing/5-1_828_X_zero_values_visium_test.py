"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/828
https://github.com/chanzuckerberg/single-cell-curation/pull/876
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_non_visium_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_visium_some,
    validator_with_all_visiums,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"  # all spots library_id
SOME_LIB_ID = "Visium_11_CK289"     # some spots library_id


def test_visiums_pass(validator_with_all_visiums):
    validator = validator_with_all_visiums
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_visiums_pass_np_array(validator_with_all_visiums):
    validator = validator_with_all_visiums
    validator.adata.X = validator.adata.X.toarray()
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# only raw layer is checked
def test_visiums_with_normalized_all_zero_row(validator_with_all_visiums):
    validator = validator_with_all_visiums
    normalized = np.log1p(validator.adata.X)
    raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
    normalized[100] = 0     # this row in_tissue == 0
    validator.adata.X = normalized
    validator.adata.raw = raw_adata
    del validator.adata.raw.var["feature_is_filtered"]
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# only raw layer is checked
def test_visiums_with_normalized_all_zero_row_in_tissue_1(validator_with_all_visiums):
    validator = validator_with_all_visiums
    normalized = np.log1p(validator.adata.X)
    raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
    normalized[100] = 0     # this row in_tissue == 0
    validator.adata.X = normalized
    validator.adata.raw = raw_adata
    validator.adata.obs.loc[validator.adata.obs.index[100], "in_tissue"] = 1
    del validator.adata.raw.var["feature_is_filtered"]
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_visiums_with_normalized_and_raw_all_zero_row(validator_with_all_visiums):
    validator = validator_with_all_visiums
    normalized = np.log1p(validator.adata.X)
    raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
    normalized[100] = 0     # this row in_tissue == 0
    raw_adata.X[100] = 0
    validator.adata.X = normalized
    validator.adata.raw = raw_adata
    del validator.adata.raw.var["feature_is_filtered"]
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_visiums_with_normalized_and_raw_all_zero_row_not_valid(validator_with_all_visiums):
    validator = validator_with_all_visiums
    normalized = np.log1p(validator.adata.X)
    raw_adata = ad.AnnData(validator.adata.X, var=validator.adata.var, obs=validator.adata.obs)
    normalized[100] = 0     # this row in_tissue == 0
    raw_adata.X[100] = 0
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


def test_in_tissue_one(validator_with_visium_some):
    validator = validator_with_visium_some
    validator.adata.obs.loc[validator.adata.obs.index[2895], "in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Each observation with obs['in_tissue'] == 1 must have at least one non-zero value in its row in the raw matrix.", 
        "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
    ]


def test_make_all_zero_row(validator_with_visium):
    validator = validator_with_visium
    validator.adata.X = validator.adata.X.toarray()
    validator.adata.X[100] = 0  # this row has in_tissue == 0
    validator.adata.X = sparse.csr_matrix(validator.adata.X)
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_make_all_zero_row_raw_non_visium(validator_with_non_visium_adatas):
    validator = validator_with_non_visium_adatas
    raw_adata = validator.adata.raw.X.toarray()
    raw_adata[100] = 0
    raw_adata = ad.AnnData(sparse.csr_matrix(raw_adata), var=validator.adata.raw.var, obs=validator.adata.obs)
    validator.adata.raw = raw_adata
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
        "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements."
    ]


def test_make_all_zero_row_non_visium(validator_with_non_visium_adatas):
    validator = validator_with_non_visium_adatas
    validator.adata.X = validator.adata.raw.X
    del validator.adata.raw
    validator.adata.var["feature_is_filtered"] = False
    validator.adata.X = validator.adata.X.toarray()
    validator.adata.X[100] = 0
    validator.adata.X = sparse.csr_matrix(validator.adata.X)
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
        "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
    ]


def test_is_single_false(validator_with_visium_some):
    validator = validator_with_visium_some
    validator.adata.uns["spatial"]["is_single"] = False
    validator.adata.obs["is_primary_data"] = False
    del validator.adata.uns["spatial"][SOME_LIB_ID]
    del validator.adata.obs["array_col"]
    del validator.adata.obs["array_row"]
    del validator.adata.obs["in_tissue"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
        "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
    ]


def test_is_single_np_array(validator_with_visium_some):
    validator = validator_with_visium_some
    validator.adata.uns["spatial"]["is_single"] = False
    validator.adata.obs["is_primary_data"] = False
    del validator.adata.uns["spatial"][SOME_LIB_ID]
    del validator.adata.obs["array_col"]
    del validator.adata.obs["array_row"]
    del validator.adata.obs["in_tissue"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.", 
        "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
    ]
