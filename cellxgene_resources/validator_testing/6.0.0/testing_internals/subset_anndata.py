import pytest
from IPython.display import display
import anndata as ad
import numpy as np
import dask.array as da
from scipy import sparse
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    FIXTURES_ROOT
)
from cellxgene_schema.validate import Validator


def back_to_dask(adata:ad.AnnData) -> Validator:
    adata.X = da.from_array(sparse.csr_matrix(adata.X))
    if adata.raw:
        raw_adata = ad.AnnData(
            X = da.from_array(sparse.csr_matrix(adata.raw.X)),
            obs=adata.obs,
            var=adata.raw.var,
        )
        adata.raw = raw_adata
    validator=Validator()
    validator.adata = adata
    return validator


@pytest.fixture
def subset_adata(fixture_file:str):
    '''
    Subsets anndata at both raw.X and .X and replaces matrices with valid matrices according to schema 6.0

    :rtype Anndata n_obs = 5 and n_var = 5.
    '''

    adata = read_h5ad(FIXTURES_ROOT / fixture_file)
    new_adata = adata[:5 , :5].copy()  # subsets all attributes, except adata.raw
    raw_matrix = np.array([
        [1, 0, 0, 0, 1],
        [1, 2, 0, 0, 0],
        [1, 0, 3, 0, 0],
        [1, 0, 0, 4, 0],
        [1, 0, 0, 0, 5],
    ],dtype=np.float32)
    norm_matrix = np.array([
        [1.1, 0, 0, 0, 1.1],
        [1.1, 2.2, 0, 0, 0],
        [1.1, 0, 3.3, 0, 0],
        [1.1, 0, 0, 4.4, 0],
        [1.1, 0, 0, 0, 5.5]
    ],dtype=np.float32)

    if new_adata.raw:
        raw_var = new_adata.raw.var

    else:
        raw_var = new_adata.var
        raw_var.drop(columns=["feature_is_filtered"],inplace=True)

    new_adata.X = norm_matrix
    new_raw_adata = ad.AnnData(
        X = raw_matrix,
        obs = new_adata.obs[:5],
        var =  raw_var[:5]
    )
    new_adata.raw = new_raw_adata
    new_adata.var['feature_is_filtered'] = False
    return new_adata
