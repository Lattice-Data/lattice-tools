"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1158
https://github.com/chanzuckerberg/single-cell-curation/pull/1169/
"""

import anndata as ad
import numpy as np
import pytest
from dask.array import from_array
from scipy import sparse
from fixtures.valid_adatas import validator_with_all_adatas


# might be memory intensive, probably best not to run in parallel
def test_np_in_raw_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas
    var = validator.adata.raw.var if validator.adata.raw else validator.adata.var

    # .compute() changes from dask array to csr matrix within the validator
    if validator.adata.raw:
        matrix = validator.adata.raw.X.compute()
    else:
        matrix = validator.adata.X.compute()

    # then can convert to dense array
    matrix = matrix.toarray()
    assert isinstance(matrix, np.ndarray)

    # reset raw layer to have dense array
    raw_adata = ad.AnnData(
        X=from_array(matrix),
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    # validator check for sparsity, use here to match error string
    nnz = validator.count_matrix_nonzero(validator.adata.raw.X)
    sparsity = 1 - nnz / np.prod(matrix.shape)

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'raw.X' is {sparsity} which is greater than 0.5, and it is "
        "not a 'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_np_in_x_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas

    # .compute() changes from dask array to csr matrix within the validator
    matrix = validator.adata.X.compute()

    # then can convert to dense array
    matrix = matrix.toarray()
    assert isinstance(matrix, np.ndarray)

    # reset X to the dense matrix and as dask array
    validator.adata.X = from_array(matrix)

    # validator check for sparsity, use here to match error string
    nnz = validator.count_matrix_nonzero(validator.adata.X)
    sparsity = 1 - nnz / np.prod(matrix.shape)

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'X' is {sparsity} which is greater than 0.5, and it is "
        "not a 'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_np_in_layer_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas

    # .compute() changes from dask array to csr matrix within the validator
    matrix = validator.adata.X.compute()

    # then can convert to dense array
    matrix = matrix.toarray()
    assert isinstance(matrix, np.ndarray)

    # set dense X to layer
    layer = "my_layer"
    validator.adata.layers[layer] = from_array(matrix)

    # validator check for sparsity, use here to match error string
    nnz = validator.count_matrix_nonzero(validator.adata.layers[layer])
    sparsity = 1 - nnz / np.prod(matrix.shape)

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        f"ERROR: Sparsity of 'layers['{layer}']' is {sparsity} which is greater than 0.5, and it is "
        "not a 'scipy.sparse.csr_matrix'. The matrix MUST use this type of matrix for the given sparsity."
    ) in validator.errors


def test_csc_sparse_in_raw_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas
    var = validator.adata.raw.var if validator.adata.raw else validator.adata.var

    # .compute() changes from dask array to csr matrix within the validator
    if validator.adata.raw:
        matrix = validator.adata.raw.X.compute()
    else:
        matrix = validator.adata.X.compute()

    # then can convert to csc matrix
    matrix = sparse.csc_matrix(matrix)
    assert isinstance(matrix, sparse.csc_matrix)

    # reset raw layer to have csc matrix
    raw_adata = ad.AnnData(
        X=from_array(matrix),
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        "ERROR: Invalid sparse encoding for raw.X with encoding csc. Only csr sparse encodings are supported."
    ) in validator.errors


def test_csc_in_x_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas

    # .compute() changes from dask array to csr matrix within the validator
    matrix = validator.adata.X.compute()

    # then can convert to csc matrix
    matrix = sparse.csc_matrix(matrix)
    assert isinstance(matrix, sparse.csc_matrix)

    # reset X to the csc matrix and as dask array
    validator.adata.X = from_array(matrix)

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        "ERROR: Invalid sparse encoding for X with encoding csc. Only csr sparse encodings are supported."
    ) in validator.errors


def test_csc_in_layer_fails(validator_with_all_adatas):
    validator = validator_with_all_adatas

    # .compute() changes from dask array to csr matrix within the validator
    matrix = validator.adata.X.compute()

    # then can convert to csc matrix
    matrix = sparse.csc_matrix(matrix)
    assert isinstance(matrix, sparse.csc_matrix)

    # set csc X to layer
    layer = "my_layer"
    validator.adata.layers[layer] = from_array(matrix)

    validator.validate_adata()

    assert not validator.is_valid
    assert (
        f"ERROR: Invalid sparse encoding for layers['{layer}'] with encoding csc. "
        "Only csr sparse encodings are supported."
    ) in validator.errors


def test_total_dense_in_raw_passes(validator_with_all_adatas):
    validator = validator_with_all_adatas
    var = validator.adata.raw.var if validator.adata.raw else validator.adata.var

    matrix = np.ones(shape=validator.adata.shape, dtype=np.float32)
    assert isinstance(matrix, np.ndarray)

    # reset raw layer to have csc matrix
    raw_adata = ad.AnnData(
        X=from_array(matrix),
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    if "feature_is_filtered" in validator.adata.raw.var.columns:
        validator.adata.raw.var.drop(columns="feature_is_filtered", inplace=True)

    validator.validate_adata()

    assert validator.is_valid
    assert validator.errors == []


def test_total_dense_x_passes(validator_with_all_adatas):
    validator = validator_with_all_adatas

    matrix = np.ones(shape=validator.adata.shape, dtype=np.float32)
    assert isinstance(matrix, np.ndarray)

    validator.adata.X = from_array(matrix)

    # set all feature_is_filtered to False for total density
    validator.adata.var["feature_is_filtered"] = False

    validator.validate_adata()

    assert validator.is_valid
    assert validator.errors == []


def test_total_dense_layer_passes(validator_with_all_adatas):
    validator = validator_with_all_adatas

    matrix = np.ones(shape=validator.adata.shape, dtype=np.float32)
    assert isinstance(matrix, np.ndarray)

    validator.adata.layers["my_layer"] = from_array(matrix)

    validator.validate_adata()

    assert validator.is_valid
    assert validator.errors == []


def test_half_dense_x(validator_with_all_adatas):
    validator = validator_with_all_adatas

    matrix = np.ones(shape=validator.adata.shape, dtype=np.float32)

    # should be sparsity of 0.5, every other element set to 0
    # adatas with non-even shape will be slightly off 0.5 sparsity
    matrix[:,::2] = 0
    assert isinstance(matrix, np.ndarray)

    validator.adata.X = from_array(matrix)

    # validator check for sparsity, use here to match error string
    nnz = validator.count_matrix_nonzero(validator.adata.X)
    sparsity = 1 - nnz / np.prod(matrix.shape)

    validator.validate_adata()

    # hacky way to get around sparsity of 0.500002342 or something
    if sparsity > 0.5:
        assert not validator.is_valid
    else:
        assert validator.is_valid

    assert validator.errors == []


def test_half_dense_raw(validator_with_all_adatas):
    validator = validator_with_all_adatas
    var = validator.adata.raw.var if validator.adata.raw else validator.adata.var

    matrix = np.ones(shape=validator.adata.shape, dtype=np.float32)
    # slice and take every row, set every other value to 0
    # with even-shaped adata, end up with 0.5 sparsity
    matrix[:,::2] = 0
    assert isinstance(matrix, np.ndarray)

    # reset raw layer to have csc matrix
    raw_adata = ad.AnnData(
        X=from_array(matrix),
        obs=validator.adata.obs,
        var=var,
    )
    validator.adata.raw = raw_adata

    if "feature_is_filtered" in validator.adata.raw.var.columns:
        validator.adata.raw.var.drop(columns="feature_is_filtered", inplace=True)

    # validator check for sparsity, use here to match error string
    nnz = validator.count_matrix_nonzero(validator.adata.raw.X)
    sparsity = 1 - nnz / np.prod(matrix.shape)

    validator.validate_adata()

    # hacky way to get around sparsity of 0.500002342 or something
    if sparsity > 0.5:
        assert not validator.is_valid
    else:
        assert validator.is_valid

    assert validator.errors == []
