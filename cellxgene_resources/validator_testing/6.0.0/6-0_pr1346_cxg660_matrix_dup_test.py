"""
QA testing for this issue: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-curation/1317
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1346

should not pass:
 - duplicate a cell’s raw counts once (raw is in .X / raw is in raw.X)
 - duplicate many cell's raw counts (raw is in .X / raw is in raw.X)
 - duplicate a cell's raw counts many times (raw is in .X / raw is in raw.X)

"""

import pytest
from testing_internals.create_duplications import create_duplications
from fixtures.valid_adatas import (
    ALL_H5ADS,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
    test_h5ads,
    validator_with_adatas

)

ERROR_MESSAGES = (
    #*** NEED TO FILL IN ***
)


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_one_duplication_non_visium(validator_with_adatas,error):
    validator = validator_with_adatas
    barcode = validator.adata.obs.index[0]
    print(barcode)
    dup_adata = create_duplications(validator.adata, i=barcode,n=1)
    if validator.adata.raw:
        validator.adata.raw = dup_adata
    else:
        validator.adata = dup_adata
    validator.validate_adata()
    assert not validator.is_valid


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
@pytest.mark.parametrize("error", ERROR_MESSAGES)
def test_one_duplication_visium(validator_with_adatas,error):
    validator = validator_with_adatas
    if "in_tissue" in validator.adata.obs.columns:
        obs_to_keep = validator.adata.obs[validator.adata.obs["in_tissue"] != 0].index

    validator.adata = validator.adata[obs_to_keep, :]
    barcode = validator.adata.obs.index[0]
    print(barcode)
    dup_adata = create_duplications(validator.adata, i=barcode,n=1)
    if validator.adata.raw:
        validator.adata.raw = dup_adata
    else:
        validator.adata = dup_adata
    validator.validate_adata()
    assert not validator.is_valid