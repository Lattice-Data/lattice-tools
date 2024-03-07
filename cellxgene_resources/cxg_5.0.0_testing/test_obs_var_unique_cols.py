"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/729

Ways pandas DF can allow duplicate column names:
    concat same DF onto itself with second argument containing name of col that already exists
    two or more cols are renamed to the same string
"""

import anndata as ad
import pandas as pd
import pytest
from fixtures.valid_adatas import validator_with_adata, Validator, FIXTURES_ROOT


@pytest.mark.parametrize(
    "duplicate_col",
    (
        pytest.param('suspension_type', id="duplicated obs col 'suspension_type'"),
        pytest.param('donor_id', id="duplicated obs col 'donor_id'"),
        pytest.param('author_test', id="duplicated obs author col 'author_test'"),
    )
)
def test_obs_concat_unique_fails(validator_with_adata, duplicate_col):
    validator = validator_with_adata
    validator.adata.obs['author_test'] = 'test'
    validator.adata.obs = pd.concat([validator.adata.obs, validator.adata.obs[duplicate_col]], axis=1)

    with pytest.raises(( ValueError, AssertionError )):
        validator.validate_adata()


@pytest.mark.parametrize(
    "duplicate_col",
    (
        pytest.param('feature_is_filtered', id="duplicated var col 'feature_is_filtered'"),
    )
)
def test_var_concat_unique_fails(validator_with_adata, duplicate_col):
    validator = validator_with_adata
    validator.adata.var = pd.concat([validator.adata.var, validator.adata.var[duplicate_col]], axis=1)

    with pytest.raises(( ValueError, AssertionError )):
        validator.validate_adata()


@pytest.mark.parametrize(
    "duplicate_col,col_rename_dict",
    (
        pytest.param('cluster', {'cluster_color': 'cluster', 'size': 'cluster'}, id="rename author obs col"),
        pytest.param('suspension_type', {'cluster_color': 'suspension_type'}, id="rename to required obs col"),
    )
)
def test_obs_rename_unique_fails(validator_with_adata, duplicate_col, col_rename_dict):
    validator = validator_with_adata

    # no clear common columns between mouse and human h5ads, using this to match human obs to mouse obs
    for col in ['cluster_color', 'size']:
        try:
            validator.adata.obs[col]
        except KeyError:
            validator.adata.obs[col] = True

    validator.adata.obs.rename(columns=col_rename_dict, inplace=True)

    with pytest.raises(( ValueError, AssertionError )):
        validator.validate_adata()


@pytest.mark.parametrize(
    "duplicate_col,col_rename_dict",
    (
        pytest.param('feature_is_filtered', {'test': 'feature_is_filtered'}, id="rename var col to required name"),
        pytest.param('gene_symbol', {'test': 'gene_symbol'}, id="rename var col to author name"),
    )
)
def test_var_rename_unique_fails(duplicate_col, validator_with_adata, col_rename_dict):
    validator = validator_with_adata
    validator.adata.var['test'] = 'test'
    validator.adata.var['gene_symbol'] = 'test'
    validator.adata.var.rename(columns=col_rename_dict, inplace=True)

    with pytest.raises(( ValueError, AssertionError )):
        validator.validate_adata()


# testing for raw.var duplicates, more complicated pathway we take to merge in raw
# matrix and raw.var
@pytest.mark.parametrize(
    "h5ad,duplicate_col",
    (
        pytest.param('valid_human.h5ad', 'gene_symbol', id="duplicated raw var 'gene_symbol'"),
        pytest.param('valid_human.h5ad', 'ensembl_version', id="duplicated raw var 'ensembl_version'"),
        pytest.param('valid_mouse.h5ad', 'gene_symbol', id="duplicated raw var 'gene_symbol'"),
        pytest.param('valid_mouse.h5ad', 'ensembl_version', id="duplicated raw var 'ensembl_version'"),
    )
)
def test_raw_var_concat_fails(h5ad, duplicate_col, fixtures_root=FIXTURES_ROOT):
    adata = ad.read_h5ad(f'{fixtures_root}/{h5ad}')
    var = adata.raw.var
    var[duplicate_col] = 'test'
    new_var = pd.concat([var, var[duplicate_col]], axis=1)
    raw_adata = ad.AnnData(adata.raw.X, obs=adata.obs, var=new_var)
    adata.raw = raw_adata
    
    # quick check to make sure only raw.var contains duplicates
    print(adata.var)
    print(adata.raw.var)

    validator = Validator()
    validator.adata = adata

    with pytest.raises(
        ValueError, 
        match=f"Duplicate column name '{duplicate_col}' detected in 'adata.raw.var' DataFrame. "
        "All DataFrame column names must be unique."
    ):
        validator.validate_adata()


@pytest.mark.parametrize(
    "h5ad,duplicate_col,rename_dict",
    (
        pytest.param('valid_human.h5ad', 'feature_is_filtered', {'test': 'feature_is_filtered'}, id="rename raw.var col to required name"),
        pytest.param('valid_human.h5ad', 'gene_symbol', {'test': 'gene_symbol'}, id="rename raw.var col to author name"),
        pytest.param('valid_mouse.h5ad', 'feature_is_filtered', {'test': 'feature_is_filtered'}, id="rename raw.var col to required name"),
        pytest.param('valid_mouse.h5ad', 'gene_symbol', {'test': 'gene_symbol'}, id="rename raw.var col to author name"),
    )
)
def test_raw_var_rename_fails(h5ad, duplicate_col, rename_dict, fixtures_root=FIXTURES_ROOT):
    adata = ad.read_h5ad(f'{fixtures_root}/{h5ad}')
    var = adata.raw.var
    var['test'] = 'test'
    var[duplicate_col] = 'test'
    var.rename(columns=rename_dict, inplace=True)
    raw_adata = ad.AnnData(adata.raw.X, obs=adata.obs, var=var)
    adata.raw = raw_adata
    
    # quick check to make sure only raw.var contains duplicates
    print(adata.var)
    print(adata.raw.var)

    validator = Validator()
    validator.adata = adata

    with pytest.raises(
        ValueError, 
        match=f"Duplicate column name '{duplicate_col}' detected in 'adata.raw.var' DataFrame. "
        "All DataFrame column names must be unique."
    ):
        validator.validate_adata()


def test_del_raw_var(validator_with_adata):
    validator = validator_with_adata
    del validator.adata.raw
    validator.validate_adata()
    assert validator.is_valid is False
