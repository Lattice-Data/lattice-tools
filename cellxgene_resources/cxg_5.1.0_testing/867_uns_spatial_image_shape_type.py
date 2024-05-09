"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/867
https://github.com/chanzuckerberg/single-cell-curation/pull/882
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_all_visiums,
    validator_with_visium_some,
    validator_with_non_spatial_adata,
)


def test_image_passes(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("image", ("fullres", "hires"))
def test_image_array_smaller_shape(validator_with_visium, image):
    validator = validator_with_visium
    LIBRARY_ID = [k for k in validator.adata.uns['spatial'].keys() if k not in 'is_single'][0]
    validator.adata.uns["spatial"][LIBRARY_ID]["images"][image] = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image][:,:,2]
    new_image_shape = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image].shape
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and either 3 "
        "(RGB color model for example) or 4 (RGBA color model for example) for its last dimension, it has shape"
        f" {new_image_shape}."
    ]


@pytest.mark.parametrize(
    "image,shape", 
        (
            ("hires", (2000, 1000, 5)),
            ("fullres", (2000, 1000, 5)),
            ("hires", (2000, 1000, 5, 3)),
            ("fullres", (2000, 1000, 5, 3)),
        )
)
def test_image_array_larger_shape(validator_with_all_visiums, image, shape):
    validator = validator_with_all_visiums
    LIBRARY_ID = [k for k in validator.adata.uns['spatial'].keys() if k not in 'is_single'][0]
    validator.adata.uns["spatial"][LIBRARY_ID]["images"][image] = np.zeros(shape=shape, dtype=np.uint8)
    new_image_shape = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image].shape
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['images']['{image}'] must have a length of 3 and either 3 "
        "(RGB color model for example) or 4 (RGBA color model for example) for its last dimension, it has shape"
        f" {new_image_shape}."
    ]


@pytest.mark.parametrize(
    "image,dtype", 
        (
            ("fullres", np.uint),
            ("fullres", np.int8),
            ("fullres", np.int16),
            ("fullres", np.uint16),
            ("fullres", np.float32),
            ("fullres", int),
            ("fullres", float),
            ("hires", np.uint),
            ("hires", np.int8),
            ("hires", np.int16),
            ("hires", np.uint16),
            ("hires", np.float32),
            ("hires", int),
            ("hires", float),
        )
)
def test_image_array_type(validator_with_visium, image, dtype):
    validator = validator_with_visium
    LIBRARY_ID = [k for k in validator.adata.uns['spatial'].keys() if k not in 'is_single'][0]
    uns_image = validator.adata.uns["spatial"][LIBRARY_ID]["images"][image]
    uns_image = uns_image.astype(dtype)
    validator.adata.uns["spatial"][LIBRARY_ID]["images"][image] = uns_image
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: uns['spatial'][library_id]['images']['{image}'] must be of type numpy.uint8, it is {uns_image.dtype}."
    ]
