'''
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1355

Testing conditions for visium/spatial datasets:
Datasets that are is_single:FALSE but is_primary_data:TRUE were allowed, contrary to the schema

Should pass:
- is_single:False, is_primary_data:False
- is_single:True, is_primary_data:True

Shoul not pass:
- is_single:False (bool), is_primary_data:True
- is_single:False (numpy._bool), is_primary_data:True

'''

import pytest
import anndata as ad
import numpy as np
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    SPATIAL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestSpatialData:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_valid(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []

    def test_valid_both_true(self):

        # is_single:True, is_primary_data:True -> pass

        self.validator.adata.uns["spatial"]["is_single"] = True
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_valid_npbool_both_true(self):

        # is_single:True, is_primary_data:True -> pass

        self.validator.adata.uns["spatial"]["is_single"] = np.bool_(True)
        self.validator.adata.obs["is_primary_data"] = True
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_valid_both_false(self):

        # is_single:False, is_primary_data:False -> pass

        self.validator.adata.uns["spatial"]["is_single"] = False
        self.validator.adata.obs["is_primary_data"] = False

        visium_assays = ["EFO:0010961","EFO:0022860","EFO:0022859","EFO:0022857"]
        if self.validator.adata.obs["assay_ontology_term_id"].isin(visium_assays).all():
            del self.validator.adata.obs["in_tissue"]
            del self.validator.adata.obs["array_col"]
            del self.validator.adata.obs["array_row"]
            try:
                library_id = [k for k in self.validator.adata.uns["spatial"].keys() if k != "is_single"][0]
                del self.validator.adata.uns["spatial"][library_id]
            except:
                pass

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_valid_npbool_both_false(self):

        # is_single:False, is_primary_data:False -> pass

        self.validator.adata.uns["spatial"]["is_single"] = np.bool_(False)
        self.validator.adata.obs["is_primary_data"] = False

        visium_assays = ["EFO:0010961","EFO:0022860","EFO:0022859","EFO:0022857"]
        if self.validator.adata.obs["assay_ontology_term_id"].isin(visium_assays).all():
            del self.validator.adata.obs["in_tissue"]
            del self.validator.adata.obs["array_col"]
            del self.validator.adata.obs["array_row"]
            try:
                library_id = [k for k in self.validator.adata.uns["spatial"].keys() if k != "is_single"][0]
                del self.validator.adata.uns["spatial"][library_id]
            except:
                pass

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_invalid_single_false_primary_true(self):

        # is_single:False, is_primary_data:True -> fail

        self.validator.adata.uns["spatial"]["is_single"] = False
        self.validator.adata.obs["is_primary_data"] = True

        visium_assays = ["EFO:0010961","EFO:0022860","EFO:0022859","EFO:0022857"]
        if self.validator.adata.obs["assay_ontology_term_id"].isin(visium_assays).all():
            del self.validator.adata.obs["in_tissue"]
            del self.validator.adata.obs["array_col"]
            del self.validator.adata.obs["array_row"]
            try:
                library_id = [k for k in self.validator.adata.uns["spatial"].keys() if k != "is_single"][0]
                del self.validator.adata.uns["spatial"][library_id]
            except:
                pass

        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert self.validator.errors == []


    def test_invalid_npbool_single_false_primary_true(self):

            # is_single:False, is_primary_data:True -> fail

            self.validator.adata.uns["spatial"]["is_single"] = np.bool_(False)
            self.validator.adata.obs["is_primary_data"] = True

            visium_assays = ["EFO:0010961","EFO:0022860","EFO:0022859","EFO:0022857"]
            if self.validator.adata.obs["assay_ontology_term_id"].isin(visium_assays).all():
                del self.validator.adata.obs["in_tissue"]
                del self.validator.adata.obs["array_col"]
                del self.validator.adata.obs["array_row"]
                try:
                    library_id = [k for k in self.validator.adata.uns["spatial"].keys() if k != "is_single"][0]
                    del self.validator.adata.uns["spatial"][library_id]
                except:
                    pass

            self.validator.validate_adata()
            assert not self.validator.is_valid
            #assert self.validator.errors == []
