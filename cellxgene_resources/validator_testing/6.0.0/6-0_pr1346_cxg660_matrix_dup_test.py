"""
QA testing for this issue: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-curation/1317
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1346

should not pass:
(Y) duplicate a cell’s raw counts once (raw is in .X / raw is in raw.X)
(Y) duplicate many cell's raw counts (raw is in .X / raw is in raw.X)
(Y) duplicate a cell's raw counts many times (raw is in .X / raw is in raw.X)

should pass:
(Y) raw counts in raw.X contain no duplicates + normalized in .X contains duplicates
"""

import pytest
import numpy as np
import anndata as ad
from testing_internals.create_duplications import create_duplications
from fixtures.valid_adatas import (
    ALL_H5ADS,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
    test_h5ads,
    validator_with_adatas
)


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_a_cell_one_dup_raw(validator_with_adatas):

    # duplicate a cell’s raw counts once (raw is in raw.X) -> fail

    validator = validator_with_adatas
    barcode = [validator.adata.obs.index[0]]
    if "in_tissue" in validator.adata.obs.columns:
        barcode = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[0]]

    dup_adata = create_duplications(validator.adata,barcode,n=1)
    validator.adata.raw = dup_adata
    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs adata.raw.X."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_a_cell_one_dup_raw_in_X(validator_with_adatas):

    # duplicate a cell’s raw counts once (raw is in .X) -> fail

    validator = validator_with_adatas
    barcode = [validator.adata.obs.index[0]]
    if "in_tissue" in validator.adata.obs.columns:
        barcode = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[0]]

    dup_adata = create_duplications(validator.adata,barcode,n=1)
    if validator.adata.raw:
        del validator.adata.raw
        validator.adata = dup_adata

    else:
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs adata.X."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_many_cells_dup_raw(validator_with_adatas):

    # duplicate many cell's raw counts (raw is in raw.X) -> fail

    validator = validator_with_adatas
    raw_loc = []
    n_barcodes = 5
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    if "in_tissue" in validator.adata.obs.columns:
        barcodes = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[i] for i in range(n_barcodes)]

    dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
    validator.adata.raw = dup_adata
    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs adata.raw.X."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_many_cells_dup_raw_in_X(validator_with_adatas):

    # duplicate many cell's raw counts (raw is in .X) -> fail

    validator = validator_with_adatas
    n_barcodes = 5
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    if "in_tissue" in validator.adata.obs.columns:
        barcodes = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[i] for i in range(n_barcodes)]

    dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
    if validator.adata.raw:
        del validator.adata.raw
        validator.adata = dup_adata

    else:
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs adata.X."
    ) in validator.errors

@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_one_cell_many_dups_raw(validator_with_adatas):

    # duplicate a cell's raw counts many times (raw is in raw.X) -> fail

    validator = validator_with_adatas
    raw_loc = []
    n_barcodes = 1
    n_duplications = 8
    total_duplicates = n_barcodes * (n_duplications + 1)
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    if "in_tissue" in validator.adata.obs.columns:
        barcodes = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[i] for i in range(n_barcodes)]

    dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
    validator.adata.raw = dup_adata
    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs adata.raw.X."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_one_cell_many_dups_raw_in_X(validator_with_adatas):

    # duplicate a cell's raw counts many times (raw is in .X) -> fail

    validator = validator_with_adatas
    n_barcodes = 1
    n_duplications = 8
    total_duplicates = n_barcodes * (n_duplications + 1)
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    if "in_tissue" in validator.adata.obs.columns:
        barcodes = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[i] for i in range(n_barcodes)]

    dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
    if validator.adata.raw:
        del validator.adata.raw
        validator.adata = dup_adata

    else:
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs adata.X."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_dup_in_norm(validator_with_adatas):

    # raw counts in raw.X contain no duplicates + normalized in .X contains duplicates -> pass

    validator = validator_with_adatas
    n_barcodes = 3
    n_duplications = 2
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    if "in_tissue" in validator.adata.obs.columns:
        barcodes = [validator.adata.obs[validator.adata.obs["in_tissue"]==1].index[i] for i in range(n_barcodes)]

    if validator.adata.raw:
        norm_adata = ad.AnnData(X=validator.adata.X.copy(),
                                obs=validator.adata.obs.copy(),
                                var=validator.adata.var.copy(),
                                uns=validator.adata.uns.copy(),
                                obsm=validator.adata.obsm.copy())

        dup_adata = create_duplications(norm_adata,barcodes,n=n_duplications)
        dup_adata.raw = ad.AnnData(X=validator.adata.raw.X.copy(),
                                   var=validator.adata.raw.var.copy())
        validator.adata = dup_adata
        zero_cols = np.where(np.array(validator.adata.X.compute().sum(axis=0)).flatten() == 0)[0]
        validator.adata.var.loc[validator.adata.var.index[zero_cols],"feature_is_filtered"] = True

    else:
        print("No normalized data found in AnnData object.")
        pass

    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
