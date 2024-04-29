"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/827
https://github.com/chanzuckerberg/single-cell-curation/pull/858
Further bug fixes here:
https://github.com/chanzuckerberg/single-cell-curation/issues/862
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"


def test_all_passes(validator_with_all_adatas):
    validator = validator_with_all_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_uns_spatial_fails(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.uns["spatial"] = True
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] is only allowed for obs['assay_ontology_term_id'] "
        "values 'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' "
        "(Slide-seqV2)."
    ]


def test_delete_uns_spatial(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    del validator.adata.uns["spatial"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] is required for obs['assay_ontology_term_id'] values "
        "'EFO:0010961' (Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)."
    ]


@pytest.mark.parametrize(
    "value", (None, 1, 1.0, "string value", "", [], True, False, np.array([]), pd.DataFrame([]))
)
def test_uns_spatial_is_not_dict(validator_with_visium, value):
    validator = validator_with_visium
    validator.adata.uns["spatial"] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial']['is_single'] must be of boolean type, it is {type(value)}."
    ]


def test_delete_uns_is_single_visium(validator_with_visium):
    validator = validator_with_visium
    del validator.adata.uns["spatial"]["is_single"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] must contain the key 'is_single'."
    ]


def test_delete_uns_is_single_slide_seq(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    del validator.adata.uns["spatial"]["is_single"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] must contain the key 'is_single'.",
        "ERROR: uns['spatial'] cannot be an empty value.",
    ]


@pytest.mark.parametrize(
    "value", (None, 1, 1.0, "string value", "", {}, [], np.bool_, np.array([]))
)
def test_uns_is_single_non_boolean(validator_with_spatial_adatas, value):
    validator = validator_with_spatial_adatas
    validator.adata.uns["spatial"]["is_single"] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial']['is_single'] must be of boolean type, it is {type(value)}."
    ]


def test_uns_spatial_library_id(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id] is only allowed for obs['assay_ontology_term_id'] "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


@pytest.mark.parametrize(
    "value", (None, 1, 1.0, {}, [], np.bool_, np.array([]))
)
def test_uns_library_id_type(validator_with_visium, value):
    validator = validator_with_visium
    del validator.adata.uns["spatial"][LIBRARY_ID]
    validator.adata.uns["spatial"][value] = {}
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id] must contain the key 'images'.",
        "ERROR: uns['spatial'][library_id] must contain the key 'scalefactors'."
    ]


@pytest.mark.parametrize(
    "value", (None, 1, 1.0, "string", [], np.bool_, np.array([]), pd.DataFrame([]))
)
def test_uns_library_id_value_type(validator_with_visium, value):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id] must contain the key 'images'.",
        "ERROR: uns['spatial'][library_id] must contain the key 'scalefactors'."
    ]


def test_multiple_library_ids(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["second_library_id"] = {}
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] must contain only two top-level keys: 'is_single' and a library_id. "
        "More than two top-level keys detected: ['spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna', "
        "'second_library_id']."
    ]


def test_add_library_id_keys(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["metadata"] = validator.adata.uns["spatial"][LIBRARY_ID]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == ["ERROR: uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'.Detected keys: ['images', 'scalefactors', 'metadata']."]


def test_add_library_id_keys_no_fullres(validator_with_visium):
    validator = validator_with_visium
    del validator.adata.uns["spatial"][LIBRARY_ID]["images"]["fullres"]
    validator.adata.uns["spatial"][LIBRARY_ID]["metadata"] = validator.adata.uns["spatial"][LIBRARY_ID]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == ["ERROR: uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'.Detected keys: ['images', 'scalefactors', 'metadata']."]


@pytest.mark.parametrize(
    "image,expected",
    (
        pytest.param("fullres", True, id="Delete fullres, valid"),
        pytest.param("hires", False, id="Delete hires, not valid"),
    ),
)
def test_delete_images(validator_with_visium, image, expected):
    validator = validator_with_visium
    del validator.adata.uns["spatial"][LIBRARY_ID]["images"][image]
    validator.validate_adata()
    assert validator.is_valid is expected
    if image == "hires":
        assert validator.errors == ["ERROR: uns['spatial'][library_id]['images'] must contain the key 'hires'."]


def test_add_images(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["images"]["lowres"] = validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id]['images'] can only contain the keys "
        "'fullres' and 'hires'.Detected keys: ['fullres', 'hires', 'lowres']."
    ]


def test_images_with_non_visium_assay(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.uns["spatial"] = {LIBRARY_ID: {"images": {"hires": np.ndarray([])}}}
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'] is only allowed for obs['assay_ontology_term_id'] values 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and 'EFO:0030062' (Slide-seqV2)."
    ]


@pytest.mark.parametrize(
    "image,value",
    (
        pytest.param("fullres", {}, id="fullres dict"),
        pytest.param("fullres", [], id="fullres list"),
        pytest.param("fullres", "string", id="fullres string"),
        pytest.param("fullres", pd.DataFrame([]), id="fullres dataframe"),
        pytest.param("fullres", True, id="fullres True"),
        pytest.param("fullres", False, id="fullres False"),
        pytest.param("fullres", None, id="fullres None"),
        pytest.param("fullres", 1, id="fullres int"),
        pytest.param("fullres", 1.0, id="fullres float"),
        pytest.param("hires", {}, id="hires dict"),
        pytest.param("hires", [], id="hires list"),
        pytest.param("hires", "string", id="hires string"),
        pytest.param("hires", pd.DataFrame([]), id="hires dataframe"),
        pytest.param("hires", True, id="hires True"),
        pytest.param("hires", False, id="hires False"),
        pytest.param("hires", None, id="hires None"),
        pytest.param("hires", 1, id="hires int"),
        pytest.param("hires", 1.0, id="hires float"),
    ),
)
def test_non_array_images(validator_with_visium, image, value):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["images"][image] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['images']['{image}'] must be of numpy.ndarray type, it is {type(value)}."
    ]


@pytest.mark.parametrize("image", ("fullres", "hires"))
def test_image_array_shape(validator_with_visium, image):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["images"][image] = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image][:,:,2]
    new_image_shape = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image].shape
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have shape (,,3), "
        f"it has shape {new_image_shape}."
    ]


def test_hires_dimension_too_large(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"] = validator.adata.uns["spatial"][LIBRARY_ID]["images"]["fullres"]
    new_image_shape = validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"].shape
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels, "
        f"it has a largest dimension of {max(new_image_shape)} pixels."
    ]


def test_hires_dimension_too_small(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"] = validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"][:1999,:,:]
    new_image_shape = validator.adata.uns["spatial"][LIBRARY_ID]["images"]["hires"].shape
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels, "
        f"it has a largest dimension of {max(new_image_shape)} pixels."
    ]


@pytest.mark.parametrize(
    "value", (None, 1, 1.0, "string value", "", [], True, False, np.array([]), pd.DataFrame([]))
)
def test_scalefactor_type(validator_with_visium, value):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["scalefactors"] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial']['is_single'] must be of boolean type, it is {type(value)}."
    ]


@pytest.mark.parametrize(
    "scalefactor,value",
    (
        pytest.param("spot_diameter_fullres", {}, id="spot_diameter_fullres dict"),
        pytest.param("spot_diameter_fullres", [], id="spot_diameter_fullres list"),
        pytest.param("spot_diameter_fullres", "string", id="spot_diameter_fullres string"),
        pytest.param("spot_diameter_fullres", pd.DataFrame([]), id="spot_diameter_fullres dataframe"),
        pytest.param("spot_diameter_fullres", True, id="spot_diameter_fullres True"),
        pytest.param("spot_diameter_fullres", False, id="spot_diameter_fullres False"),
        pytest.param("spot_diameter_fullres", None, id="spot_diameter_fullres None"),
        pytest.param("spot_diameter_fullres", 1, id="spot_diameter_fullres int"),
        pytest.param("tissue_hires_scalef", {}, id="tissue_hires_scalef dict"),
        pytest.param("tissue_hires_scalef", [], id="tissue_hires_scalef list"),
        pytest.param("tissue_hires_scalef", "string", id="tissue_hires_scalef string"),
        pytest.param("tissue_hires_scalef", pd.DataFrame([]), id="tissue_hires_scalef dataframe"),
        pytest.param("tissue_hires_scalef", True, id="tissue_hires_scalef True"),
        pytest.param("tissue_hires_scalef", False, id="tissue_hires_scalef False"),
        pytest.param("tissue_hires_scalef", None, id="tissue_hires_scalef None"),
        pytest.param("tissue_hires_scalef", 1, id="tissue_hires_scalef int"),
    ),
)
def test_scale_factor_types(validator_with_visium, scalefactor, value):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["scalefactors"][scalefactor] = value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['scalefactors']['{scalefactor}'] must be of type float, it is {type(value)}. "
        f"This must be the value of the {scalefactor} field from scalefactors_json.json"
    ]


def test_add_scalefactor(validator_with_visium):
    validator = validator_with_visium
    validator.adata.uns["spatial"][LIBRARY_ID]["scalefactors"]["tissue_lowres_scalef"] = 1.234
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id]['scalefactors'] can only contain the keys "
        "'spot_diameter_fullres' and 'tissue_hires_scalef'.Detected keys: ['spot_diameter_fullres', "
        "'tissue_hires_scalef', 'tissue_lowres_scalef']."
    ]
