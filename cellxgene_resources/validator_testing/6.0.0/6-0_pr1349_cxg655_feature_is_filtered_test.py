"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1298
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1349


should pass:
(Y) only .X + feature_is_filtered = False for all
(Y) raw.X & .X - feature_is_filtered = False for all, all .X genes have at least 1 non-0 value in .X
(Y) raw.X & .X - feature_is_filtered = False for all, 1 gene have all 0 values in .X & all 0 values in raw.X
(Y) raw.X & .X - feature_is_filtered = True for 1 gene, which has all 0 values in .X & at least 1 non-0 value in raw.X

Warnings:
(Errors addressed by the following tests were changed to warnings to bypass Visium "background spots" conflict (introduction of normalized data counts=0)
(Y) raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ALL .X counts != 0
(Y) raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + SOME .X counts != 0
(Y) raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ONE .X count != 0
(Y) raw.X & .X + feature present in raw.X & .X + feature_is_filtered = False + ALL .X counts == 0 + SUM raw.X counts != 0

should not pass:
(Y) only .X + feature_is_filtered = True
(Y) only .X + feature_is_filtered = True for just a single var that sums to 0 for all cells
(Y) raw.X & .X - a gene ID in raw.var.index that is not present in var.index, and vice versa

"""

import pytest
import anndata as ad
from testing_internals.subset_anndata import subset_adata,back_to_dask
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    ALL_H5ADS,
    FIXTURES_ROOT,
    test_h5ads,
    validator_with_adatas
)


@pytest.mark.parametrize("fixture_file", ALL_H5ADS)
class TestSubset:
    def test_valid(self,fixture_file):
        adata = read_h5ad(FIXTURES_ROOT / fixture_file)
        validator = Validator()
        validator.adata = adata
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_onlyX_all_false(self,subset_adata):

        # only .X + feature_is_filtered = False for all

        adata = subset_adata
        del adata.X
        adata.X = adata.raw.X
        del adata.raw
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_X_raw_all_false(self,subset_adata):

        # raw.X & .X - feature_is_filtered = False for all, all .X genes have at least 1 non-0 value in .X --> see subset_adata()

        adata = subset_adata
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_pass_1var_all0_x_raw(self,subset_adata):

        # raw.X & .X - feature_is_filtered = False for all, 1 gene have all 0 values in .X & all 0 values in raw.X

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

        # raw.X & .X - feature_is_filtered = True for 1 gene, which has all 0 values in .X & at least 1 non-0 value in raw.X

        adata = subset_adata
        gene_index = 0
        adata.X[:, gene_index] = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


    def test_fail_X_true(self,subset_adata):

        # only .X + feature_is_filtered = True

        adata = subset_adata
        adata.X = adata.raw.X
        del adata.raw
        adata.var['feature_is_filtered'] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: 'feature_is_filtered' must be False for all features if 'adata.raw' is not present."
            ) in validator.errors


    def test_fail_X_true_1_var(self,subset_adata):

        # only .X + feature_is_filtered = True for just a single var that sums to 0 for all cells

        adata = subset_adata
        adata.X = adata.raw.X
        del adata.raw
        gene_index = 0
        adata.X[:, gene_index] = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: 'feature_is_filtered' must be False for all features if 'adata.raw' is not present."
            ) in validator.errors


    def test_fail_gene_in_raw(self,subset_adata):

        # raw.X & .X - a gene ID in raw.var.index that is not present in var.index should not pass

        adata = subset_adata
        gene_index = 0
        gene_name = adata.var_names[gene_index]
        var_to_keep = [v for v in adata.var_names if v != gene_name]
        new_adata = adata[:,var_to_keep]
        raw_adata = ad.AnnData(X=adata.raw.X,obs=adata.obs,var=adata.raw.var)
        new_adata.raw = raw_adata
        validator = back_to_dask(new_adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: Number of genes in X ({len(new_adata.var)}) is different than raw.X ({len(new_adata.raw.var)})."
            ) in validator.errors


    def test_fail_gene_in_X(self,subset_adata):

        # raw.X & .X - a gene ID in var.index that is not present in raw.var.index should not pass

        adata = subset_adata
        gene_index = 0
        gene_name = adata.raw.var_names[gene_index]
        var_to_keep = [v for v in adata.raw.var_names if v != gene_name]
        raw_adata = ad.AnnData(
            X = adata.raw.X,
            obs = adata.obs,
            var = adata.raw.var
        )
        new_raw_adata = raw_adata[:,var_to_keep]
        adata.raw = new_raw_adata
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: Could not complete full validation of feature_is_filtered because of size differences between var and raw.var."
            ) in validator.errors


    def test_warn_true_allnot0(self,subset_adata):

        # raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ALL .X counts != 0

        adata = subset_adata
        gene_index = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert (
            f"WARNING: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', but there are {adata.X.shape[0]} non-zero "
                "values in the corresponding columns of the matrix 'X'. All values for these features must be 0."
                ) in validator.warnings


    @pytest.mark.parametrize("range_end", [1,2])
    def test_warn_true_some_1_not0(self,subset_adata,range_end):

        # raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + SOME .X counts != 0
        # raw.X & .X + feature present in raw.X & .X + feature_is_filtered = True + ONE .X counts != 0

        adata = subset_adata
        gene_index = 0
        r = range(0,range_end)
        adata.X[r, gene_index] = 0
        adata.var.iloc[gene_index, adata.var.columns.get_loc('feature_is_filtered')] = True
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert (
            f"WARNING: Some features are 'True' in 'feature_is_filtered' of dataframe 'var', but there are {adata.X.shape[0] - len(r)} non-zero "
                "values in the corresponding columns of the matrix 'X'. All values for these features must be 0."
                ) in validator.warnings


    def test_warn_false_X0_rawnot0(self,subset_adata):

        # raw.X & .X + feature present in raw.X & .X + feature_is_filtered = False + ALL .X counts == 0 + SUM raw.X counts != 0

        adata = subset_adata
        gene_indices = [0,1]
        gene_names = []
        for gene_index in gene_indices:
            gene_name = adata.var_names[gene_index]
            gene_names.append(gene_name)
            adata.X[:, gene_index] = 0
            assert adata.raw.X[:, gene_index].any() != 0
        validator = back_to_dask(adata)
        validator.validate_adata()
        assert validator.is_valid
        assert (
            f"WARNING: Genes '{', '.join(gene_names)}' have all-zero values in adata.X. Either feature_is_filtered should be set to True "
                "or adata.raw.X should be set to all-zero values."
                ) in validator.warnings