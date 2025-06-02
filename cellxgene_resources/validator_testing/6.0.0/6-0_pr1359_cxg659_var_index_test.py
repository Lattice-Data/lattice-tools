"""
QA testing for this issue:
PR for this issue:


Should pass:
(Y) ENSG00000290826 in var index
(Y) all genes from organism - same as fixture_pass
(Y) all genes from organism + covid
(Y) all genes from organism + spike-ins
(Y) all genes from organism + covid + spike-ins
(Y) any gene names that contain "." (edge case)

* Check add-labels: feature_name for ENSG00000290826 is ENSG00000290826 -> found in other script: 6-0_pr1359_cxg656_feature_name_test.py: currently gene_version
isn't being split from gene_id for slide-seq, visium and valid_human.h5ad

Shouldn't pass:
(Y) ENSG00000290826.1 in var index
(Y) mix-and-match genes
(Q) all genes from one organism, different uns.organism -> Need a clearer error message.
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


ORGANISM_GENE_VALUES = {'NCBITaxon:9606':'ENSG00000290826',
                        'NCBITaxon:10090':'ENSMUSG00000121346',
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
                        "NCBITaxon:9825":"ENSSSCG00000057774"
                        }

EXEMPT_ORGANISMS = {
    "NCBITaxon:2697049":"ENSSASG00005000004",  # Severe acute respiratory syndrome coronavirus 2
    "NCBITaxon:32630":"ERCC-00003"  # synthetic construct
}

EDGE_CASES_GENE_NAMES = {
    "NCBITaxon:9606":["ENSG00000245857", "GS1-24F4.2"],
    "NCBITaxon:10090":["ENSMUSG00000039337", "Tex19.2"],
    "NCBITaxon:9986":["ENSOCUG00000017166", "KAP6.1.1"],
    "NCBITaxon:1747270":["ENSSSCG00000029160", "HSP70.2"],
    "NCBITaxon:6239":["WBGene00007102", "B0024.13"]
}


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestVarIndexValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_fixture_pass(self):

        # same as all genes from organism -> pass

        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_all_genes_from_organism_covid(self):

        # all genes from organism + covid -> pass

        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:2697049"]})
        if self.validator.adata.raw:
            new_var = self.validator.adata.raw.var.rename(index={self.validator.adata.raw.var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:2697049"]})
            raw_adata = ad.AnnData(self.validator.adata.raw.X, var=new_var, obs=self.validator.adata.obs)
            self.validator.adata.raw = raw_adata

        else:
            pass

        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_all_genes_from_organism_spike_ins(self):

        # all genes from organism + spike-ins -> pass

        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:32630"]})
        if self.validator.adata.raw:
            new_var = self.validator.adata.raw.var.rename(index={self.validator.adata.raw.var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:32630"]})
            raw_adata = ad.AnnData(self.validator.adata.raw.X, var=new_var, obs=self.validator.adata.obs)
            self.validator.adata.raw = raw_adata

        else:
            pass

        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_all_genes_from_organism_both_exempt(self):

        # all genes from organism + covid + spike-ins -> pass

        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:2697049"]})
        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[1]: EXEMPT_ORGANISMS["NCBITaxon:32630"]})
        if self.validator.adata.raw:
            new_var = self.validator.adata.raw.var.copy()
            new_var = new_var.rename(index={new_var.index[0]: EXEMPT_ORGANISMS["NCBITaxon:2697049"]})
            new_var_2 = new_var.rename(index={new_var.index[1]: EXEMPT_ORGANISMS["NCBITaxon:32630"]})
            raw_adata = ad.AnnData(self.validator.adata.raw.X, var=new_var_2, obs=self.validator.adata.obs)
            self.validator.adata.raw = raw_adata

        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_var_index_ensembl(self):

        # ensembl gene ID in var index -> pass

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: ORGANISM_GENE_VALUES[organism]})

        if self.validator.adata.raw:
            new_var = self.validator.adata.raw.var.rename(index={self.validator.adata.raw.var.index[0]: ORGANISM_GENE_VALUES[organism]})
            raw_adata = ad.AnnData(self.validator.adata.raw.X, var=new_var, obs=self.validator.adata.obs)
            self.validator.adata.raw = raw_adata

        else:
            pass

        assert self.validator.adata.var.index[0] == ORGANISM_GENE_VALUES[organism]
        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_human_var_index_ensembl_decimal(self):

        # ENSG00000290826.1 in var index -> fail

        if self.validator.adata.uns["organism_ontology_term_id"] == "NCBITaxon:9606":
            decimal_id = "ENSG00000290826.1"
            self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: decimal_id})

            if self.validator.adata.raw:
                new_var = self.validator.adata.raw.var.rename(index={self.validator.adata.raw.var.index[0]: decimal_id})
                raw_adata = ad.AnnData(self.validator.adata.raw.X, var=new_var, obs=self.validator.adata.obs)
                self.validator.adata.raw = raw_adata

            self.validator.validate_adata()
            assert self.validator.adata.var.index[0] == decimal_id
            assert (
                f"ERROR: Could not infer organism from feature ID '{decimal_id}' in 'var', make sure it is a valid ID."
                ) in self.validator.errors
            assert not self.validator.is_valid

        else:
            self.validator.validate_adata()
            assert self.validator.is_valid


    def test_mix_match_genes(self):

        # mix-and-match genes -> fail

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        if organism == "NCBITaxon:9606":
            self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: ORGANISM_GENE_VALUES["NCBITaxon:10090"]})
            assert self.validator.adata.var.index[0] == ORGANISM_GENE_VALUES["NCBITaxon:10090"]
            assert not self.validator.adata.var.index[0].startswith("ENSG0")
            assert self.validator.adata.var.index[1].startswith("ENSG0")

            if self.validator.adata.raw:
                raw_adata = ad.AnnData(
                    X = da.from_array(self.validator.adata.raw.X.compute(), chunks=self.validator.adata.raw.X.chunks),
                    var = self.validator.adata.raw.var,
                    obs = self.validator.adata.obs
                    )
                raw_adata.var = raw_adata.var.rename(index={raw_adata.var.index[0]: ORGANISM_GENE_VALUES["NCBITaxon:10090"]})
                self.validator.adata.raw = raw_adata
                assert self.validator.adata.raw.var.index[0] == ORGANISM_GENE_VALUES["NCBITaxon:10090"]
                assert not self.validator.adata.raw.var.index[0].startswith("ENSG0")
                assert self.validator.adata.raw.var.index[1].startswith("ENSG0")
                self.validator.validate_adata()
                assert not self.validator.is_valid
                assert (
                    f"ERROR: uns['organism_ontology_term_id'] is '{organism}' but feature_ids are from [<SupportedOrganisms.MUS_MUSCULUS: 'NCBITaxon:10090'>]."
                    ) in self.validator.errors

        else:
            self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: ORGANISM_GENE_VALUES["NCBITaxon:9606"]})
            assert self.validator.adata.var.index[0] == ORGANISM_GENE_VALUES["NCBITaxon:9606"]
            assert self.validator.adata.var.index[0].startswith("ENSG0")
            assert not self.validator.adata.var.index[1].startswith("ENSG0")

            if self.validator.adata.raw:
                raw_adata = ad.AnnData(
                    X = da.from_array(self.validator.adata.raw.X.compute(), chunks=self.validator.adata.raw.X.chunks),
                    var = self.validator.adata.raw.var,
                    obs = self.validator.adata.obs
                    )
                raw_adata.var = raw_adata.var.rename(index={raw_adata.var.index[0]: ORGANISM_GENE_VALUES["NCBITaxon:9606"]})
                self.validator.adata.raw = raw_adata
                assert self.validator.adata.raw.var.index[0] == ORGANISM_GENE_VALUES["NCBITaxon:9606"]
                assert self.validator.adata.raw.var.index[0].startswith("ENSG0")
                assert not self.validator.adata.raw.var.index[1].startswith("ENSG0")
                self.validator.validate_adata()
                assert not self.validator.is_valid
                assert (
                    f"ERROR: uns['organism_ontology_term_id'] is '{organism}' but feature_ids are from [<SupportedOrganisms.HOMO_SAPIENS: 'NCBITaxon:9606'>]."
                    ) in self.validator.errors


    def test_different_uns_organism(self):

        # all genes from organism, different uns.organism -> fail   ###COMMENT: Passing but the error messages are dependency errors. Need clearer error message.

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        if organism == "NCBITaxon:9606":
            self.validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
            assert self.validator.adata.uns["organism_ontology_term_id"] == "NCBITaxon:10090"
        else:
            self.validator.adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
            assert self.validator.adata.uns["organism_ontology_term_id"] == "NCBITaxon:9606"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR:")


    @pytest.mark.parametrize("edge_case", [EDGE_CASES_GENE_NAMES])
    def test_edge_case(self, edge_case):

        # any gene names that contain "." -> pass

        for organism, gene_list in edge_case.items():
            if organism == self.validator.adata.uns["organism_ontology_term_id"]:
                if gene_list[0] not in self.validator.adata.var.index:
                    self.validator.adata.var = self.validator.adata.var.rename(index={self.validator.adata.var.index[0]: gene_list[0]})

                if self.validator.adata.raw:
                    if gene_list[0] not in self.validator.adata.raw.var.index:
                        raw_adata = ad.AnnData(
                            X = da.from_array(self.validator.adata.raw.X.compute(), chunks=self.validator.adata.raw.X.chunks),
                            var = self.validator.adata.raw.var,
                            obs = self.validator.adata.obs
                            )
                        raw_adata.var = raw_adata.var.rename(index={raw_adata.var.index[0]: gene_list[0]})
                        self.validator.adata.raw = raw_adata

                self.validator.validate_adata()
                assert self.validator.is_valid
                labeler = AnnDataLabelAppender(self.validator.adata)
                labeler._add_labels()
                assert labeler.adata.var.loc[gene_list[0], "feature_name"] == gene_list[1]
