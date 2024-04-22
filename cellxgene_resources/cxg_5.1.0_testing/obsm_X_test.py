"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/pull/856
"""

import numpy as np
import pytest
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_spatial_adatas,
)


def test_obsm_all_is_valid(validator_with_all_adatas):
    validator = validator_with_all_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_spatial_is_valid(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_X_spatial_fails(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obsm["X_spatial"] = validator.adata.obsm["spatial"]
    validator.adata.uns["default_embedding"] = "X_spatial"
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: Embedding key in 'adata.obsm' X_spatial cannot be used.",
    ]


def test_non_spatial_assay(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obs["assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.obs["suspension_type"] = "cell"
    validator.adata.obs.loc[:, ["suspension_type"]] = validator.adata.obs.astype("category")
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: At least one embedding in 'obsm' has to have a key with an 'X_' prefix.",
    ]
