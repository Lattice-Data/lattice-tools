"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1298
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1349


should pass:
(Y) only .X + feature_is_filtered = False for all
(Y) raw.X & .X - feature_is_filtered = False for all, all .X genes have at least 1 non-0 value in .X
(Y) raw.X & .X - feature_is_filtered = False for all, 1 gene have all 0 values in .X & all 0 values in raw.X
(Y) raw.X & .X - feature_is_filtered = True for 1 gene, which has all 0 values in .X & at least 1 non-0 value in raw.X

should not pass:
only .X + feature_is_filtered = True
only .X + feature_is_filtered = True for just a single var that sums to 0 for all cells
raw.X & .X - a gene ID in raw.var.index that is not present in var.index, and vice versa
raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ALL .X counts != 0
raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + SOME .X counts != 0
raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ONE .X count != 0
raw.X & .X + feature present in raw.X & .X + feature_is_filtered = False + ALL .X counts == 0 + SUM raw.X counts != 0

"""

import pytest
import anndata as ad
from testing_internals.subset_anndata import subset_adata,back_to_dask
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    MULTISPECIES_H5ADS,
    FIXTURES_ROOT,
    test_h5ads,
    validator_with_adatas
)


@pytest.mark.parametrize("fixture_file", MULTISPECIES_H5ADS)
class TestSubset:
    def test_valid(self,fixture_file):
        adata = read_h5ad(FIXTURES_ROOT / fixture_file)
        validator = Validator()
        validator.adata = adata
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_onlyX_all_false(self,subset_adata):
        '''
        only .X + feature_is_filtered = False for all
        '''
        adata = subset_adata
        del adata.X
        adata.X = adata.raw.X
        del adata.raw
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_X_raw_all_false(self,subset_adata):
        '''
        raw.X & .X - feature_is_filtered = False for all, all .X genes have at least 1 non-0 value in .X --> see subset_adata()
        '''
        adata = subset_adata
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_1var_all0_x_raw(self,subset_adata):
        '''
        raw.X & .X - feature_is_filtered = False for all, 1 gene have all 0 values in .X & all 0 values in raw.X
        '''
        adata = subset_adata
        gene_index = 0
        adata.X[:, gene_index] = 0
        raw_matrix = adata.raw.X
        raw_matrix[:, gene_index] = 0
        adata.raw = ad.AnnData(X=raw_matrix,obs = adata.obs[:],var = adata.raw.var[:])
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    def test_pass_0s_in_X_true(self,subset_adata):
        '''
        raw.X & .X - feature_is_filtered = True for 1 gene, which has all 0 values in .X & at least 1 non-0 value in raw.X
        '''
        adata = subset_adata
        gene_index = 0
        adata.X[:, gene_index] = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_fail_X_true(self,subset_adata):
        '''
        only .X + feature_is_filtered = True
        '''
        adata = subset_adata
        adata.X = adata.raw.X
        del adata.raw
        adata.var['feature_is_filtered'] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (f"ERROR: 'feature_is_filtered' must be False for all features if 'adata.raw' is not present."
        ) in validator.errors


    def test_fail_X_true_single_var(self,subset_adata):

        #only .X + feature_is_filtered = True for just a single var that sums to 0 for all cells

        adata = subset_adata
        adata.X = adata.raw.X
        del adata.raw
        gene_index = 0
        adata.X[:, gene_index] = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (f"ERROR: 'feature_is_filtered' must be False for all features if 'adata.raw' is not present."
        ) in validator.errors

'''
    def test_fail_(self,subset_adata):

        #raw.X & .X - a gene ID in raw.var.index that is not present in var.index, and vice versa

        adata = subset_adata
        adata.X[] = 0 #modify dense arrays
        raw_matrix = adata.raw.X
        raw_matrix[] = 0 #modify = ad.AnnData(modified )
        adata.raw = ad.AnnData(X=raw_matrix, ... )
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []



    '''def test_raw_(validator_with_adata):
        validator = validator_with_adata

        if validator.adata.raw:
            raw_adata = ad.AnnData(validator.adata.raw.X, obs=validator.adata.obs, var=validator.adata.var, uns=validator.adata.uns, obsm=validator.adata.obsm)
            raw_adata.var['feature_is_filtered'] = True
            validator.adata = raw_adata
            del validator.adata.raw

        else:
            validator.adata.var['feature_is_filtered'] = True

        n_var = validator.adata.shape[1]
        validator.validate_adata()
        assert not validator.is_valid
        assert validator.errors == [
            f"ERROR: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', but there are {n_var} non-zero values in the corresponding columns of the matrix 'X'. All values for these features must be 0."
        ]'''
