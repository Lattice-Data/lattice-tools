"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/827
https://github.com/chanzuckerberg/single-cell-curation/pull/858
Further bug fixes here:
https://github.com/chanzuckerberg/single-cell-curation/issues/862
https://github.com/chanzuckerberg/single-cell-curation/pull/875

Lots of overlap with updating to visium CytAssist and 11mm, so relevant links:
https://github.com/chanzuckerberg/single-cell-curation/pull/1129
https://github.com/chanzuckerberg/single-cell-curation/issues/1107
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"
VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]

def get_library_id(adata):
    return [key for key in adata.uns['spatial'].keys() if 'is_single' not in key][0]


def test_all_passes(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestNonSpatialUns:
    def test_uns_spatial_fails(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"] = True
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'] is only allowed when obs['assay_ontology_term_id']"
            " is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression)"
            " or 'EFO:0030062' (Slide-seqV2)"
        ) in validator.errors

    def test_images_with_non_visium_assay(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"] = {LIBRARY_ID: {"images": {"hires": np.ndarray([])}}}
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: uns['spatial'] is only allowed when obs['assay_ontology_term_id'] is either"
            " a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)"
        ]


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestSpatialAdatasUns:
    @pytest.mark.parametrize(
        "value", (None, 1, 1.0, "string value", "", {}, [], np.bool_, np.array([]))
    )
    def test_uns_is_single_non_boolean(self, validator_with_adatas, value):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial']['is_single'] must be of boolean type, it is {type(value)}."
        ) in validator.errors

    def test_delete_uns_spatial(self, validator_with_adatas):
        validator = validator_with_adatas
        del validator.adata.uns["spatial"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: A dict in uns['spatial'] is required when obs['assay_ontology_term_id']"
            " is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)."
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
def test_delete_uns_is_single_slide_seq(validator_with_adatas):
    validator = validator_with_adatas
    del validator.adata.uns["spatial"]["is_single"]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[:2] == [
        "ERROR: uns['spatial'] must contain the key 'is_single'.",
        "ERROR: uns['spatial'] cannot be an empty value.",
    ]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestVisiumUns:
    @pytest.mark.parametrize(
        "value", (None, 1, 1.0, "string value", "", [], True, False, np.array([]), pd.DataFrame([1, 2, 3], index=["a", "b", "c"]))
    )
    def test_uns_spatial_is_not_dict(self, validator_with_adatas, value):
        validator = validator_with_adatas
        validator.adata.uns["spatial"] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: A dict in uns['spatial'] is required when obs['assay_ontology_term_id'] "
            "is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2)."
        ) in validator.errors

    def test_delete_uns_is_single_visium(self, validator_with_adatas):
        validator = validator_with_adatas
        del validator.adata.uns["spatial"]["is_single"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'] must contain the key 'is_single'."
        ) in validator.errors

    def test_uns_spatial_library_id(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors[0] == "ERROR: uns['spatial'][library_id] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."

    @pytest.mark.parametrize(
        "value", (True, False, None, 1, 1.0, np.bool_)
    )
    def test_uns_library_id_type(self, validator_with_adatas, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]
        validator.adata.uns["spatial"][value] = {}
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: uns['spatial'][library_id] must contain the key 'images'.",
            "ERROR: uns['spatial'][library_id] must contain the key 'scalefactors'."
        ]

    # mutable types raise TypeError: unhashable type
    @pytest.mark.parametrize(
        "value", ({}, [], np.array([]))
    )
    def test_uns_library_id_unhashable_type(self, validator_with_adatas, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]
        with pytest.raises(TypeError):
            validator.adata.uns["spatial"][value] = {}
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: uns['spatial'] must contain at least one key representing the library_id when "
            "obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) "
            "and uns['spatial']['is_single'] is True."
        ]

    @pytest.mark.parametrize(
        "value", (None, True, False, 1, 1.0, "string", [], np.bool_, np.array([]), pd.DataFrame([]))
    )
    def test_uns_library_id_value_type(self, validator_with_adatas, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: uns['spatial'][library_id] must be a dictionary.",
        ]

    def test_multiple_library_ids(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"]["second_library_id"] = {}
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial'] must contain only two top-level keys: 'is_single' and a library_id. "
            f"More than two top-level keys detected: ['{library_id}', "
            f"'second_library_id']."
        ) in validator.errors

    def test_add_library_id_keys(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["metadata"] = validator.adata.uns["spatial"][library_id]
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'."
            "Detected keys: ['images', 'scalefactors', 'metadata']." 
        ) in validator.errors

    def test_add_library_id_keys_no_fullres(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]["images"]["fullres"]
        validator.adata.uns["spatial"][library_id]["metadata"] = validator.adata.uns["spatial"][library_id]
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'][library_id] can only contain the keys 'images' and 'scalefactors'."
            "Detected keys: ['images', 'scalefactors', 'metadata']." 
        ) in validator.errors

    @pytest.mark.parametrize(
        "image,expected",
        (
            pytest.param("fullres", True, id="Delete fullres, valid"),
            pytest.param("hires", False, id="Delete hires, not valid"),
        ),
    )
    def test_delete_images(self, validator_with_adatas, image, expected):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        del validator.adata.uns["spatial"][library_id]["images"][image]
        validator.validate_adata()
        assert validator.is_valid is expected
        if image == "hires":
            assert validator.errors == ["ERROR: uns['spatial'][library_id]['images'] must contain the key 'hires'."]

    def test_add_images(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"]["lowres"] = validator.adata.uns["spatial"][library_id]["images"]["hires"]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: uns['spatial'][library_id]['images'] can only contain the keys "
            "'fullres' and 'hires'.Detected keys: ['fullres', 'hires', 'lowres']."
        ]

    @pytest.mark.parametrize("image", ["fullres", "hires"])
    @pytest.mark.parametrize(
        "value",
        ([{}, [], "string", pd.DataFrame([]), True, False, None, 1, 1.0])
    )
    def test_non_array_images(self, validator_with_adatas, image, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"][image] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] "
            f"must be of numpy.ndarray type, it is {type(value)}." 
        ) in validator.errors

    @pytest.mark.parametrize("image", ("fullres", "hires"))
    def test_image_array_shape(self, validator_with_adatas, image):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"][image] = validator.adata.uns["spatial"][library_id]["images"][image][:,:,2]
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"][image].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and "
            "either 3 (RGB color model for example) or 4 (RGBA color model for example) for its "
            f"last dimension, it has shape {new_image_shape}."
        ) in validator.errors
        
    def test_hires_dimension_too_large(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"]["hires"] = validator.adata.uns["spatial"][library_id]["images"]["fullres"]
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"]["hires"].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels, "
            f"it has a largest dimension of {max(new_image_shape)} pixels."
        ) in validator.errors

    def test_hires_dimension_too_small(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"]["hires"] = validator.adata.uns["spatial"][library_id]["images"]["hires"][:1999,:,:]
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"]["hires"].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: The largest dimension of uns['spatial'][library_id]['images']['hires'] must be 2000 pixels, "
            f"it has a largest dimension of {max(new_image_shape)} pixels."
        ) in validator.errors

    @pytest.mark.parametrize(
        "value", (None, 1, 1.0, "string value", "", [], True, False, np.array([]), pd.DataFrame([]))
    )
    def test_scalefactor_type(self, validator_with_adatas, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["scalefactors"] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'][library_id]['scalefactors'] must be a dictionary."
        ) in validator.errors

    @pytest.mark.parametrize("scalefactor", ["spot_diameter_fullres", "tissue_hires_scalef"])
    @pytest.mark.parametrize(
        "value",
        (
            pytest.param({}, id= "dict"),
            pytest.param([], id="list"),
            pytest.param("string", id="string"),
            pytest.param(pd.DataFrame([]), id="dataframe"),
            pytest.param(True, id="True"),
            pytest.param(False, id="False"),
            pytest.param(None, id="None"),
            pytest.param(1, id="int"),
        ),
    )
    def test_scale_factor_types(self, validator_with_adatas, scalefactor, value):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["scalefactors"][scalefactor] = value
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial'][library_id]['scalefactors']['{scalefactor}'] must be of type float, it is {type(value)}. "
            f"This must be the value of the {scalefactor} field from scalefactors_json.json"
        ) in validator.errors

    def test_add_scalefactor(self, validator_with_adatas):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["scalefactors"]["tissue_lowres_scalef"] = 1.234
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            "ERROR: uns['spatial'][library_id]['scalefactors'] can only contain the keys "
            "'spot_diameter_fullres' and 'tissue_hires_scalef'.Detected keys: ['spot_diameter_fullres', "
            "'tissue_hires_scalef', 'tissue_lowres_scalef']."
        ) in validator.errors
