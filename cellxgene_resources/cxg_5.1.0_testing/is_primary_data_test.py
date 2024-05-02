"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/834
https://github.com/chanzuckerberg/single-cell-curation/pull/865
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"


def test_is_primary_data_false_slide_seq(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
    ]


def test_is_primary_data_false(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[4] == "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."


def test_non_spatial_is_primary_false(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# could be other instances of data in corpus even though uns["spatial"]["is_single"] is True
def test_spatial_is_primary_false(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obs["is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_one_random_is_primary_false_visium(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["is_single"] = False
    random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
    validator.adata.obs.loc[validator.adata.obs.index[random_index], "is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[4] == "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."


def test_one_random_is_primary_false_slide_seq(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    validator.adata.uns["spatial"]["is_single"] = False
    random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
    validator.adata.obs.loc[validator.adata.obs.index[random_index], "is_primary_data"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: When uns['spatial']['is_single'] is False, obs['is_primary_data'] must be False for all rows."
    ]
