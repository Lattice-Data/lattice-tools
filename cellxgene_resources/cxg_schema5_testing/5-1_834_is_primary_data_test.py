"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/834
https://github.com/chanzuckerberg/single-cell-curation/pull/865
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


class TestIsPrimaryDataAll:
    def test_is_primary_data_true(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = True
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # spatial fixtures start at is_primary_data = True, non-spatial has some False values
    def test_is_primary_data_mixed(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[:100], "is_primary_data"] = False
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # already present checks, don't think we have a specific test for this yet
    @pytest.mark.parametrize(
        "type", (None, 1, 1.0, "True", "False")
    )
    def test_is_primary_data_types(self, validator_with_adatas, type):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = type
        col = validator.adata.obs["is_primary_data"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: Column 'is_primary_data' in dataframe 'obs' must be boolean, not '{col.dtype}'."
        ]


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
class TestSlideSeqIsPrimaryData:
    def test_is_primary_data_false_slide_seq(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
        ]

    def test_one_random_is_primary_false_slide_seq(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
        validator.adata.obs.loc[validator.adata.obs.index[random_index], "is_primary_data"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
        ]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestVisiumIsPrimaryData:
    def test_is_primary_data_false(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors[4] == "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."

    def test_one_random_is_primary_false_visium(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
        validator.adata.obs.loc[validator.adata.obs.index[random_index], "is_primary_data"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors[4] == "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_non_spatial_is_primary_false(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# could be other instances of data in corpus even though uns["spatial"]["is_single"] is True
@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_spatial_is_primary_false(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
