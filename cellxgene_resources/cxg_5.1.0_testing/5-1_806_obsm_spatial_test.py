"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/806
https://github.com/chanzuckerberg/single-cell-curation/pull/863
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_all_visiums,
    validator_with_non_spatial_adata,
)

def get_library_id(adata):
    return [key for key in adata.uns['spatial'].keys() if 'is_single' not in key][0]


def test_obsm_spatial_is_valid(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_obsm_warnings(validator_with_all_visiums):
    validator = validator_with_all_visiums
    validator.adata.obsm["new_spatial"] = validator.adata.obsm["spatial"]
    del validator.adata.obsm["spatial"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        "WARNING: Embedding key in 'adata.obsm' new_spatial is not 'spatial' nor does it start with 'X_'. "
        "Thus, it will not be available in Explorer"
    ) in validator.warnings


def test_obsm_visium_is_single_false(validator_with_visium):
    validator = validator_with_visium
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


#changing X shape raises AnnData exception, this should be alright instead of validator error
def test_less_rows_X(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    with pytest.raises(ValueError):
        validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:1000, :]
        validator.validate_adata()
    assert validator.is_valid is False


def test_less_than_two_cols(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"][:, :1]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: All 'X_' and 'spatial' embeddings must have at least two columns. "
        f"'adata.obsm['spatial']' has columns='{validator.adata.obsm['spatial'].shape[1]}'."
    ]


# OverflowError when trying to set to inf
@pytest.mark.parametrize(
    "obsm_value", (np.inf, np.NINF)
)
def test_obsm_inf_visium(validator_with_visium, obsm_value):
    validator = validator_with_visium
    with pytest.raises(OverflowError):
        validator.adata.obsm["spatial"][0, 0] = obsm_value
    # validator.validate_adata()
    # assert validator.is_valid is False
    # assert validator.errors == ["ERROR: adata.obsm['spatial'] contains positive infinity or negative infinity values."]


@pytest.mark.parametrize(
    "obsm_value", (np.inf, np.NINF)
)
def test_obsm_inf_slide_seq(validator_with_slide_seq_adatas, obsm_value):
    validator = validator_with_slide_seq_adatas
    validator.adata.obsm["spatial"][0, 0] = obsm_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == ["ERROR: adata.obsm['spatial'] contains positive infinity or negative infinity values."]


def test_obsm_np_nan_slide_seq(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    validator.adata.obsm["spatial"][:, :] = np.nan
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: adata.obsm['spatial'] contains at least one NaN value."
    ]


def test_obsm_np_nan_spatial(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    # coerce to float, visium test datasets seem to have int for obsm['spatial']
    validator.adata.obsm["spatial"] = validator.adata.obsm["spatial"].astype("float32")
    validator.adata.obsm["spatial"][0,0] = np.nan
    validator.validate_adata()
    assert validator.is_valid is False
    assert "ERROR: adata.obsm['spatial'] contains at least one NaN value." in validator.errors


# visium obsm raises ValueError when set to np.nan
def test_obsm_np_nan_visium(validator_with_all_visiums):
    validator = validator_with_all_visiums
    with pytest.raises(ValueError):
        validator.adata.obsm["spatial"][:, :] = np.nan
    # validator.validate_adata()
    # assert validator.is_valid is False
    # assert validator.errors == [
    #     "ERROR: adata.obsm['spatial'] contains all NaN values."
    # ]
