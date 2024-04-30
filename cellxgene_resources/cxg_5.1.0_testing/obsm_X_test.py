"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/pull/856
"""

import pytest
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_spatial_adatas,
    validator_with_non_spatial_adata,
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


# X_ reserved, other key names should be allowed
@pytest.mark.parametrize(
    "obsm_key,expected", 
    (
        ("X_spatial", False),
        ("x_spatial", True),
        ("X_Spatial", False),
        ("X__Spatial", False),
        ("x__Spatial", True),
        ("x__spatial", True),
    )
)
def test_X_spatial_variations(validator_with_spatial_adatas, obsm_key, expected):
    validator = validator_with_spatial_adatas
    validator.adata.obsm[obsm_key] = validator.adata.obsm["spatial"]
    validator.adata.uns["default_embedding"] = obsm_key
    validator.validate_adata()
    assert validator.is_valid is expected


# will need further tie-in with uns['spatial']['is_single']
def test_non_spatial_assay(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obs["assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.obs["suspension_type"] = "cell"
    validator.adata.obs.loc[:, ["suspension_type"]] = validator.adata.obs.astype("category")
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] is only allowed for obs['assay_ontology_term_id'] values 'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2).",
    ]


@pytest.mark.parametrize("obsm_key", ("spatial", "X_spatial", "random_name"))
def test_non_spatial_adata(validator_with_non_spatial_adata, obsm_key):
    validator = validator_with_non_spatial_adata
    validator.adata.obsm[obsm_key] = validator.adata.obsm["X_umap"]
    del validator.adata.obsm["X_umap"]
    del validator.adata.obsm["X_pca"]
    del validator.adata.obsm["X_varimax"]
    validator.validate_adata()
    assert validator.is_valid is False
