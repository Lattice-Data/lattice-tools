"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/806
https://github.com/chanzuckerberg/single-cell-curation/pull/863
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    get_library_id,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)

VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestSpatialObsmSpatial:
    def test_obsm_spatial_is_valid(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    #changing X shape raises AnnData exception, this should be alright instead of validator error
    def test_less_rows_X(self, validator_with_adatas):
        validator = validator_with_adatas
        with pytest.raises(ValueError):
            validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:1000, :]
            validator.validate_adata()
        assert validator.is_valid is False

    def test_less_than_two_cols(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:, :1]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. "
            f"'adata.obsm['spatial']' has columns='{validator.adata.obsm['spatial'].shape[1]}'."
        ]

    def test_obsm_np_nan_spatial(self, validator_with_adatas):
        validator = validator_with_adatas
        # coerce to float, visium test datasets seem to have int for obsm['spatial']
        validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"].astype("float32")
        validator.adata.obsm["spatial"][0,0] = np.nan
        validator.validate_adata()
        assert validator.is_valid is False
        assert "ERROR: adata.obsm['spatial'] contains at least one NaN value." in validator.errors


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
class TestObsmSpatialSlideSeq:
    @pytest.mark.parametrize(
        "obsm_value", (np.inf, np.NINF)
    )
    def test_obsm_inf_slide_seq(self, validator_with_adatas, obsm_value):
        validator = validator_with_adatas
        validator.adata.obsm["spatial"][0, 0] = obsm_value
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == ["ERROR: adata.obsm['spatial'] contains positive infinity or negative infinity values."]

    def test_obsm_np_nan_slide_seq(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["spatial"][:, :] = np.nan
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: adata.obsm['spatial'] contains at least one NaN value."
        ]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestObsmSpatialVisium:
    def test_obsm_warnings(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obsm["new_spatial"] = validator.adata.obsm["spatial"]
        del validator.adata.obsm["spatial"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "WARNING: Embedding key in 'adata.obsm' new_spatial is not 'spatial' nor does it start with 'X_'. "
            "Thus, it will not be available in Explorer"
        ) in validator.warnings

    def test_obsm_visium_is_single_false(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]
        validator.validate_adata()
        assert validator.is_valid is False
        errors = [
            (
                "ERROR: obs['array_col'] is only allowed for obs['assay_ontology_term_id'] "
                "is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
            ),
            (
                "ERROR: obs['array_row'] is only allowed for obs['assay_ontology_term_id'] "
                "is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
            ),
            (
                "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] "
                "is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
            )
        ]
        for error in errors:
            assert error in validator.errors

    # OverflowError when trying to set to inf
    @pytest.mark.parametrize(
        "obsm_value", (np.inf, np.NINF)
    )
    def test_obsm_inf_visium(self, validator_with_adatas, obsm_value):
        validator = validator_with_adatas
        with pytest.raises(OverflowError):
            validator.adata.obsm["spatial"][0, 0] = obsm_value
        # validator.validate_adata()
        # assert validator.is_valid is False
        # assert validator.errors == ["ERROR: adata.obsm['spatial'] contains positive infinity or negative infinity values."]

    # visium obsm raises ValueError when set to np.nan
    def test_obsm_np_nan_visium(self, validator_with_adatas):
        validator = validator_with_adatas
        with pytest.raises(ValueError):
            validator.adata.obsm["spatial"][:, :] = np.nan
        # validator.validate_adata()
        # assert validator.is_valid is False
        # assert validator.errors == [
        #     "ERROR: adata.obsm['spatial'] contains all NaN values."
        # ]
