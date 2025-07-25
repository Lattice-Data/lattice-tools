"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1254
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1359


From tickets:
- https://lattice.atlassian.net/browse/CXG-656



Testing conditions:

- check add-labels
    - get genes from a variety of gtfs that do not have gene_name specified. feature_name should be redundant with the index value
    - ensure that feature_name is not redundant for handful of cases where gene_name is specified

Should pass:
- feature_name is redundant with index value for genes that do not have gene_name specified

Shouldn't pass:
- feature_name is redundant with index value for genes that do have gene_name specified



"""

import pytest
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import dask.array as da
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)


ORGANISM_GENE_VALUES = {"NCBITaxon:9606":"ENSG00000290826",
                        "NCBITaxon:10090":"ENSMUSG00000121346",
                        "NCBITaxon:6239": "WBGene00001706",
                        "NCBITaxon:9483": "ENSCJAG00000060785",
                        "NCBITaxon:7955": "ENSDARG00000031930",
                        "NCBITaxon:7227": "FBtr0304423_df_nrg",
                        "NCBITaxon:9595": "ENSGGOG00000008017",
                        "NCBITaxon:9541": "ENSMFAG00000026972",
                        "NCBITaxon:9544": "ENSMMUG00000017368",
                        "NCBITaxon:30608": "ENSMICG00000032327",
                        "NCBITaxon:9986": "ENSOCUG00000032260",
                        "NCBITaxon:9598": "ENSPTRG00000045894",
                        "NCBITaxon:10116": "ENSRNOG00000070886",
                        "NCBITaxon:9823": "ENSSSCG00000057774",
                        "NCBITaxon:1747270": "ENSSSCG00000057774",
                        "NCBITaxon:1654737": "ENSMMUG00000053841",
                        "NCBITaxon:9825": "ENSSSCG00000059238",
                        }


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestFeatureNameValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def fixture_pass(self):
        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_feature_name_add_label(self):

        # add_labels check: var.feature_name should be redundant with the var index value -> pass


        organism = self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: ORGANISM_GENE_VALUES[organism]})

        if self.validator.adata.raw:
            raw_adata = ad.AnnData(X = da.from_array(self.validator.adata.raw.X.compute(),chunks=self.validator.adata.raw.X.chunks),
                                    var = self.validator.adata.raw.var.copy(),
                                    obs = self.validator.adata.obs.copy()
                                    )
            raw_adata.var = raw_adata.var.rename(index={raw_adata.var.index[0]:ORGANISM_GENE_VALUES[organism]})
            self.validator.adata.raw = raw_adata

        self.validator.validate_adata()
        assert self.validator.is_valid
        labeler = AnnDataLabelAppender(self.validator.adata)
        labeler._add_labels()
        assert labeler.adata.var.loc[ORGANISM_GENE_VALUES[organism], 'feature_name'] == labeler.adata.var.index[0]


