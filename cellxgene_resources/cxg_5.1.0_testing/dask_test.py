"""
Tests for basics with dask arrays from CXG validator.

If tests do not touch the matrices, then nothing extra to be done

Tests that modify a matrix need to wrap matrix assignment back to AnnData object with the
from_array() function as follows (adding in new matrix):
    ones = np.ones(shape=validator.adata.shape, dtype=np.float32)
    validator.adata.X = from_array(sparse.csr_matrix(ones))

Also possible to copy matrix, modify, then reassign to X or raw.X:
    x = validator.adata.raw.X.copy()
    x[100] = 1.0
    del validator.adata.raw
    validator.adata.X = x

Can also call compute() method, modify, then reassign with from_array():
    x = validator.adata.raw.X.compute()
    x[100] = 1.0
    del validator.adata.raw
    validator.adata.X = from_array(x)

Copy or compute might fully load the matrix into memory but this should not cause issue
as most of the test fixtures are small in size.

This seems to then allow validation to occur as expected:
    The matrix is modified per the test
    Matrix is then dask array as expected for the validator

Some further examples from this PR, look to updated tests:
https://github.com/chanzuckerberg/single-cell-curation/pull/1152/files
"""

import anndata as ad
import numpy as np
from numpy.core.multiarray import dtype
import pandas as pd
import pytest
from dask.array import from_array
from scipy import sparse
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_non_visium_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_visium_some,
    validator_with_all_visiums,
    validator_with_non_spatial_adata,
)


def test_np_zeros_from_array(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    var = validator.adata.raw.var
    raw_adata = ad.AnnData(
        X=from_array(np.zeros(shape=validator.adata.shape, dtype=np.float32)),
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    nnz = validator.count_matrix_nonzero(validator.adata.raw.X)
    sparsity = 1 - nnz / np.prod(validator.adata.raw.X.shape)

    validator.validate_adata()
    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'raw.X' is {sparsity} which is greater than 0.5, and it is "
        "not a 'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_csr_matrix_from_array(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    var = validator.adata.raw.var
    zeros = np.zeros(shape=validator.adata.shape, dtype=np.float32)
    x = from_array(sparse.csr_matrix(zeros))
    raw_adata = ad.AnnData(
        X=x,
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    validator.validate_adata()
    assert not validator.is_valid
    errors = [
        "ERROR: Each cell must have at least one non-zero value in its row in the raw matrix.",
        "ERROR: Raw data may be missing: data in 'raw.X' does not meet schema requirements."
    ]
    for error in errors:
        assert error in validator.errors


def test_only_x_from_array(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    zeros = np.zeros(shape=validator.adata.shape, dtype=np.float32)
    del validator.adata.raw
    validator.adata.X = from_array(zeros)
    validator.validate_adata()
    
    nnz = validator.count_matrix_nonzero(validator.adata.X)
    sparsity = 1 - nnz / np.prod(validator.adata.X.shape)

    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'X' is {sparsity} which is greater than 0.5, and it is not a "
        "'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_change_one_matrix_row(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    zeros = np.zeros(shape=validator.adata.shape, dtype=np.float32)
    zeros[0] = 1
    del validator.adata.raw
    validator.adata.X = from_array(zeros)
    validator.validate_adata()
    
    nnz = validator.count_matrix_nonzero(validator.adata.X)
    sparsity = 1 - nnz / np.prod(validator.adata.X.shape)

    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'X' is {sparsity} which is greater than 0.5, and it is not a "
        "'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_ones_np_matrix(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    ones = np.ones(shape=validator.adata.shape, dtype=np.float32)
    del validator.adata.raw
    validator.adata.X = from_array(ones)
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# looks like copy() works to grab matrix and then modify it
def test_change_original_csr_matrix(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    x = validator.adata.raw.X.copy()
    x[100] = 1.0
    del validator.adata.raw
    validator.adata.X = x
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# can also call compute(), change matrix, then use from_array() to reassign
def test_change_original_csr_matrix_passes(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    x = validator.adata.raw.X.compute()
    x[100] = 1.0
    del validator.adata.raw
    validator.adata.X = from_array(x)
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
