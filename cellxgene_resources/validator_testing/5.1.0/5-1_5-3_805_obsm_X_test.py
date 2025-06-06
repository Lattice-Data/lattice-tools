"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/805
https://github.com/chanzuckerberg/single-cell-curation/pull/856
"""

import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)

SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


def test_obsm_all_is_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_spatial_is_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
def test_slide_seq_no_X(validator_with_adatas):
    validator = validator_with_adatas
    del validator.adata.obsm['X_pca']
    del validator.adata.obsm['X_umap']
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# X_ reserved, other key names should be allowed
@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestSpatialAdatasObsm:
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
    def test_X_spatial_variations(self, validator_with_adatas, obsm_key, expected):
        validator = validator_with_adatas
        validator.adata.obsm[obsm_key] = validator.adata.obsm["spatial"]
        validator.adata.uns["default_embedding"] = obsm_key
        validator.validate_adata()
        assert validator.is_valid is expected

    # will need further tie-in with uns['spatial']['is_single']
    def test_non_spatial_assay(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0009922"
        validator.adata.obs["suspension_type"] = "cell"
        validator.adata.obs.loc[:, ["suspension_type"]] = validator.adata.obs.astype("category")
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'] is only allowed when obs['assay_ontology_term_id'] is either a "
            "descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)"
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
@pytest.mark.parametrize("obsm_key", ("spatial", "X_spatial", "random_name"))
def test_non_spatial_adata(validator_with_adatas, obsm_key):
    validator = validator_with_adatas
    validator.adata.obsm[obsm_key] = validator.adata.obsm["X_umap"]
    for key in [
        "X_umap",
        "X_pca",
        "X_harmony",    # random embedding in valid_mouse
        "X_varimax"
    ]:
        if key in validator.adata.obsm.keys():
            del validator.adata.obsm[key]
    validator.validate_adata()
    assert validator.is_valid is False
