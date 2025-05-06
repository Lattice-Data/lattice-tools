"""
QA testing for this issue: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-curation/1317
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1346

should not pass:
 - duplicate a cell’s raw counts once (raw is in .X / raw is in raw.X)
 - duplicate many cell's raw counts (raw is in .X / raw is in raw.X)
 - duplicate a cell's raw counts many times (raw is in .X / raw is in raw.X)

should pass:
 - raw counts in raw.X contain no duplicates + normalized in .X contains duplicates
"""

import pytest
from testing_internals.create_duplications import create_duplications
import anndata as ad
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


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_a_cell_one_dup_non_visium(validator_with_adatas):
    validator = validator_with_adatas
    barcode = [validator.adata.obs.index[0]]
    dup_adata = create_duplications(validator.adata,barcode,n=1)
    raw_loc = []

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        validator.adata.raw = dup_adata
    else:
        raw_loc.append("adata.X")
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_a_cell_one_dup_non_visium_raw_in_X(validator_with_adatas):
    validator = validator_with_adatas
    barcode = [validator.adata.obs.index[0]]
    dup_adata = create_duplications(validator.adata,barcode,n=1)
    raw_loc = []

    if validator.adata.raw:
        del validator.adata.raw
        validator.adata = dup_adata
        raw_loc.append("adata.X")
    else:
        raw_loc.append("adata.X")
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_a_cell_one_dup_visium(validator_with_adatas):
    validator = validator_with_adatas

    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"]!= 0].index
        validator.adata = validator.adata[obs_to_keep, :].copy()

    barcode = [validator.adata.obs.index[0]]
    dup_adata = create_duplications(validator.adata,barcode,n=1)
    raw_loc = []

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        validator.adata.raw = dup_adata
    else:
        raw_loc.append("adata.X")
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_a_cell_one_dup_visium_raw_in_X(validator_with_adatas):
    validator = validator_with_adatas

    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"]!= 0].index
        validator.adata = validator.adata[obs_to_keep, :].copy()

    barcode = [validator.adata.obs.index[0]]
    dup_adata = create_duplications(validator.adata,barcode,n=1)
    raw_loc = []

    if validator.adata.raw:
        del validator.adata.raw
        validator.adata = dup_adata
        raw_loc.append("adata.X")
    else:
        raw_loc.append("adata.X")
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found 2 duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_many_cells_dup_non_visium(validator_with_adatas):
    validator = validator_with_adatas
    raw_loc = []
    n_barcodes = 5
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        raw_adata = ad.AnnData(validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.var, uns=validator.adata.uns, obsm=validator.adata.obsm)
        dup_adata = create_duplications(raw_adata,barcodes,n=n_duplications)
        validator.adata.raw = dup_adata

    else:
        raw_loc.append("adata.X")
        dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_many_cells_dup_visium(validator_with_adatas):
    validator = validator_with_adatas

    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"]!= 0].index
        validator.adata = validator.adata[obs_to_keep, :].copy()

    raw_loc = []
    n_barcodes = 5
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        raw_adata = ad.AnnData(validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.var, uns=validator.adata.uns, obsm=validator.adata.obsm)
        dup_adata = create_duplications(raw_adata,barcodes,n=n_duplications)
        validator.adata.raw = dup_adata

    else:
        raw_loc.append("adata.X")
        dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_one_cell_many_dups_non_visium(validator_with_adatas):
    validator = validator_with_adatas
    raw_loc = []
    n_barcodes = 1
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 8
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        raw_adata = ad.AnnData(validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.var, uns=validator.adata.uns, obsm=validator.adata.obsm)
        dup_adata = create_duplications(raw_adata,barcodes,n=n_duplications)
        validator.adata.raw = dup_adata

    else:
        raw_loc.append("adata.X")
        dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_one_cell_many_dups_visium(validator_with_adatas):
    validator = validator_with_adatas

    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"]!= 0].index
        validator.adata = validator.adata[obs_to_keep, :].copy()

    raw_loc = []
    n_barcodes = 1
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 8
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw:
        raw_loc.append("adata.raw.X")
        raw_adata = ad.AnnData(validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.var, uns=validator.adata.uns, obsm=validator.adata.obsm)
        dup_adata = create_duplications(raw_adata,barcodes,n=n_duplications)
        validator.adata.raw = dup_adata

    else:
        raw_loc.append("adata.X")
        dup_adata = create_duplications(validator.adata,barcodes,n=n_duplications)
        validator.adata = dup_adata

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Found {total_duplicates} duplicated raw counts in obs {raw_loc[0]}."
    ) in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_dup_in_norm_non_visium(validator_with_adatas):
    validator = validator_with_adatas
    n_barcodes = 3
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw is not None:
        norm_adata = ad.AnnData(X=validator.adata.X.copy(),
                                obs=validator.adata.obs.copy(),
                                var=validator.adata.var.copy(),
                                uns=validator.adata.uns.copy(),
                                obsm=validator.adata.obsm.copy())
        dup_adata = create_duplications(norm_adata,barcodes,n=n_duplications)
        dup_adata.raw = ad.AnnData(X=validator.adata.raw.X.copy(),
                                   var=validator.adata.raw.var.copy())

        validator.adata = dup_adata

    else:
        print("No normalized data found in AnnData object.")
        pass

    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_dup_in_norm_visium(validator_with_adatas):
    validator = validator_with_adatas

    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"]!= 0].index
        validator.adata = validator.adata[obs_to_keep, :].copy()

    n_barcodes = 3
    barcodes = [validator.adata.obs.index[i] for i in range(n_barcodes)]
    n_duplications = 2
    total_duplicates = n_barcodes * (n_duplications + 1)

    if validator.adata.raw is not None:
        norm_adata = ad.AnnData(X=validator.adata.X.copy(),
                                obs=validator.adata.obs.copy(),
                                var=validator.adata.var.copy(),
                                uns=validator.adata.uns.copy(),
                                obsm=validator.adata.obsm.copy())
        dup_adata = create_duplications(norm_adata,barcodes,n=n_duplications)
        dup_adata.raw = ad.AnnData(X=validator.adata.raw.X.copy(),
                                   var=validator.adata.raw.var.copy())

        validator.adata = dup_adata

    else:
        print("No normalized data found in AnnData object.")
        pass

    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []