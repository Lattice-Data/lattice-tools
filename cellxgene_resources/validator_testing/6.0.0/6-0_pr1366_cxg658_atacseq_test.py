"""
QA testing for this issue:
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1366

Should pass:
(Y) human ATAC fragment w/ human h5ad
(Y) mouse ATAC fragment w/ mouse h5ad

Shouldn't pass:
(Y) ATAC fragment associated with a dataset with uns.organism_ontology_term_id on approved list but not human or mouse
(Y) mouse ATAC fragment w/ human h5ad
(Y) human ATAC fragment w/ mouse h5ad
"""

import anndata as ad
import gc
import pandas as pd
import numpy as np
import pytest
import pyarrow
from cellxgene_schema.atac_seq import (
    process_fragment,
)
from pandas._libs.parsers import STR_NA_VALUES
from pathlib import Path
from fixtures.create_fixtures import Organism
from fixtures.valid_adatas import (
    ATAC_H5ADS,
    FIXTURES_ROOT,
    NON_SPATIAL_H5ADS,
    read_h5ad,
    yield_atac_fixture_data,
    yield_atac_h5ads,
    _to_anndata_file,
    to_temp_files
)

HUMAN_H5ADS = [file for file in NON_SPATIAL_H5ADS if "human" in file]
MOUSE_H5ADS = [file for file in NON_SPATIAL_H5ADS if "mouse" in file]
NON_HUMAN_MOUSE_H5ADS = [file for file in NON_SPATIAL_H5ADS if "human" not in file and "mouse" not in file]

ASSAYS = ["EFO:0030007", "EFO:0008925", "EFO:0008904", "EFO:0022045", "EFO:0030059"]

@pytest.mark.parametrize("yield_atac_h5ads", HUMAN_H5ADS)
class TestHumanATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, yield_atac_fixture_data):
        self.atac = yield_atac_fixture_data


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_human_h5ad(self, assay, tmpdir):

        # human ATAC fragment w/ human h5ad -> pass

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:9606" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_human_h5ad(self, assay, tmpdir):

        # mouse ATAC fragment w/ human h5ad -> fail

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:9606" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        frag_df = pd.read_csv(str(FIXTURES_ROOT / f"valid_mouse_fragments.tsv.gz"), sep='\t',header=None)
        self.atac.fragment_df = frag_df
        self.atac.fragment_file_name = "valid_mouse_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results != []
        assert (
            "Errors found in Fragment and/or Anndata file"
        ) in results
        assert (
            "Barcodes don't match anndata.obs.index"
        ) in results
        assert (
            f"Chromosomes in the fragment do not match the organism({organism})."
        ) in [r for r in results if "Chromosomes" in r][0].split("\n")



@pytest.mark.parametrize("yield_atac_h5ads", MOUSE_H5ADS)
class TestMouseATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, yield_atac_fixture_data):
        self.atac = yield_atac_fixture_data


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_mouse_h5ad(self, assay, tmpdir):

        # mouse ATAC fragment w/ mouse h5ad -> pass

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_mouse_h5ad(self, assay, tmpdir):

        # human ATAC fragment w/ mouse h5ad -> fail

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        frag_df = pd.read_csv(str(FIXTURES_ROOT / f"valid_human_fragments.tsv.gz"), sep='\t',header=None)
        self.atac.fragment_df = frag_df
        self.atac.fragment_file_name = "valid_human_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results != []
        assert (
            "Errors found in Fragment and/or Anndata file"
        ) in results
        assert (
            "Barcodes don't match anndata.obs.index"
        ) in results
        assert (
            f"Chromosomes in the fragment do not match the organism({organism})."
        ) in [r for r in results if "Chromosomes" in r][0].split("\n")


@pytest.mark.parametrize("yield_atac_h5ads", NON_HUMAN_MOUSE_H5ADS)
class TestOtherATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, yield_atac_fixture_data):
        self.atac = yield_atac_fixture_data


    @pytest.mark.parametrize("assay", ASSAYS)
    @pytest.mark.parametrize("organism_h_m", ["human","mouse"])
    def test_human_mouse_atac_fragment_w_other_h5ad_unmatched_barcodes(self, assay, tmpdir, organism_h_m):

       # human/mouse ATAC fragment w/ other organism h5ads without match barcodes -> fail

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" != organism
        assert "NCBITaxon:9606" != organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        frag_df = pd.read_csv(str(FIXTURES_ROOT / f"valid_{organism_h_m}_fragments.tsv.gz"), sep='\t',header=None)
        self.atac.fragment_df = frag_df
        self.atac.fragment_file_name = f"valid_{organism_h_m}_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results != []
        assert (
            "Errors found in Anndata file. Skipping fragment validation."
        ) in results
        assert (
            f"Anndata.obs.organism_ontology_term_id must be one of ['NCBITaxon:9606', 'NCBITaxon:10090']. Got {organism}."
        ) in results



    @pytest.mark.parametrize("assay",ASSAYS)
    @pytest.mark.parametrize("organism_h_m", ["human","mouse"])
    def test_human_mouse_atac_fragment_w_other_h5ad_matching_barcodes(self, assay, tmpdir, organism_h_m):

        # human/mouse ATAC fragment w/ other organism h5ads with matching barcodes -> fail

        organism = self.atac.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" != organism
        assert "NCBITaxon:9606" != organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.atac.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        frag_df = pd.read_csv(str(FIXTURES_ROOT / f"valid_{organism_h_m}_fragments.tsv.gz"), sep='\t',header=None)
        self.atac.fragment_df = frag_df
        self.atac.adata.obs.reset_index(inplace=True)
        self.atac.fragment_df[3] = self.atac.adata.obs["obs_index"].sample(n=len(self.atac.fragment_df),replace=True).values
        self.atac.adata.obs.set_index("obs_index",inplace=True)
        assert self.atac.fragment_df[3][0] in self.atac.adata.obs.index
        self.atac.fragment_file_name = f"valid_{organism_h_m}_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results != []
        assert (
            "Errors found in Anndata file. Skipping fragment validation."
        ) in results
        assert (
            f"Anndata.obs.organism_ontology_term_id must be one of ['NCBITaxon:9606', 'NCBITaxon:10090']. Got {organism}."
        ) in results
