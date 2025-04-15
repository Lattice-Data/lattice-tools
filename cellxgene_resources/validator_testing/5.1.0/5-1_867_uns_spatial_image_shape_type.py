"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/867
https://github.com/chanzuckerberg/single-cell-curation/pull/882
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    get_library_id,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestVisiumUnsSpatialImageShape:
    @pytest.mark.parametrize("image", ("fullres", "hires"))
    def test_image_array_smaller_shape(self, validator_with_adatas, image):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        if image not in validator.adata.uns["spatial"][library_id]["images"]:
            validator.adata.uns["spatial"][library_id]["images"][image] = validator.adata.uns["spatial"][library_id]["images"]["hires"]
        validator.adata.uns["spatial"][library_id]["images"][image] = validator.adata.uns["spatial"][library_id]["images"][image][:,:,2]
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"][image].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert (
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and either 3 "
            "(RGB color model for example) or 4 (RGBA color model for example) for its last dimension, it has shape"
            f" {new_image_shape}."
        ) in validator.errors

    @pytest.mark.parametrize("image", ["fullres", "hires"])
    @pytest.mark.parametrize(
        "dtype", [np.uint, np.int8, np.int16, np.uint16, np.float32, int, float]
    )
    def test_image_array_type(self, validator_with_adatas, image, dtype):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        if image in validator.adata.uns["spatial"][library_id]["images"]:
            uns_image = validator.adata.uns["spatial"][library_id]["images"][image]
        else:
            uns_image = validator.adata.uns["spatial"][library_id]["images"]["hires"]
        uns_image = uns_image.astype(dtype)
        validator.adata.uns["spatial"][library_id]["images"][image] = uns_image
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] must be of type numpy.uint8, it is {uns_image.dtype}."
        ]


@pytest.mark.parametrize("image", ["hires", "fullres"])
class TestVisiumImageSize:
    @pytest.mark.parametrize(
        "shape", 
            (
                (2000, 1000, 5),
                (2000, 1000, 5),
                (2000, 1000, 5, 3),
                (2000, 1000, 5, 3),
            )
    )
    @pytest.mark.parametrize("test_h5ads", ["visium_human_all_spots.h5ad", "visium_human_some_spots.h5ad"])
    def test_image_array_larger_shape_v1(self, validator_with_adatas, image, shape):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"][image] = np.zeros(shape=shape, dtype=np.uint8)
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"][image].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and either 3 "
            "(RGB color model for example) or 4 (RGBA color model for example) for its last dimension, it has shape"
            f" {new_image_shape}."
        ]

    @pytest.mark.parametrize(
        "shape", 
            (
                (4000, 1000, 5),
                (4000, 1000, 5),
                (4000, 1000, 5, 3),
                (4000, 1000, 5, 3),
            )
    )
    @pytest.mark.parametrize("test_h5ads", ["visium_v2_11mm_human.h5ad"])
    def test_image_array_larger_shape_v2(self, validator_with_adatas, image, shape):
        validator = validator_with_adatas
        library_id = get_library_id(validator.adata)
        validator.adata.uns["spatial"][library_id]["images"][image] = np.zeros(shape=shape, dtype=np.uint8)
        new_image_shape = validator.adata.uns["spatial"][library_id]["images"][image].shape
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and either 3 "
            "(RGB color model for example) or 4 (RGBA color model for example) for its last dimension, it has shape"
            f" {new_image_shape}."
        ]
