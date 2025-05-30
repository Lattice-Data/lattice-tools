"""
QA testing for this issue:
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1366

Should pass:
(Y) human ATAC fragment w/ human h5ad
(Y) mouse ATAC fragment w/ mouse h5ad

Shouldn't pass:
(N) ATAC fragment associated with a dataset with uns.organism_ontology_term_id on approved list but not human or mouse
(N) mouse ATAC fragment w/ human h5ad -> no errors
(N) human ATAC fragment w/ mouse h5ad -> no errors
"""

import anndata as ad
import gc
import pandas as pd
import numpy as np
import pytest
import pyarrow
from cellxgene_schema.atac_seq import (
    check_anndata_requires_fragment,
    process_fragment,
    validate_anndata_with_fragment
)
from pandas._libs.parsers import STR_NA_VALUES
from pathlib import Path
from fixtures.create_fixtures import Organism
from cellxgene_schema.validate import Validator
from fixtures.valid_adatas import (
    ATAC_H5ADS,
    FIXTURES_ROOT,
    NON_SPATIAL_H5ADS,
    read_h5ad,
    validator_with_adatas,
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
        self.validator = Validator()
        self.validator.adata = read_h5ad(FIXTURES_ROOT / self.atac.h5ad_file_name)


    def test_valid_human(self):
        print(self.atac)
        print(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_human_h5ad(self, assay, tmpdir):

        # human ATAC fragment w/ human h5ad -> pass

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:9606" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_human_h5ad(self, assay, tmpdir):

        # mouse ATAC fragment w/ human h5ad -> fail  ### doesn't error

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:9606" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        self.atac.fragment_file_name = "valid_mouse_fragments.tsv.gz"
        print(self.atac.fragment_df)
        temp_files = to_temp_files(self.atac, tmpdir)
        print(temp_files)
        results = process_fragment(**temp_files)
        print(results)
        assert results != []


@pytest.mark.parametrize("yield_atac_h5ads", MOUSE_H5ADS)
class TestMouseATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, yield_atac_fixture_data):
        self.atac = yield_atac_fixture_data
        self.validator = Validator()
        self.validator.adata = read_h5ad(FIXTURES_ROOT / self.atac.h5ad_file_name)


    def test_valid_mouse(self):
        print(self.atac)
        print(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_mouse_h5ad(self, assay, tmpdir):

        # mouse ATAC fragment w/ mouse h5ad -> pass

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        temp_files = to_temp_files(self.atac, tmpdir)
        results = process_fragment(**temp_files)
        assert results == []


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_mouse_h5ad(self, assay, tmpdir):

        # human ATAC fragment w/ mouse h5ad -> fail ### doesnt errror

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" == organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        self.atac.fragment_file_name = "valid_human_fragments.tsv.gz"
        print(self.atac.fragment_df)
        temp_files = to_temp_files(self.atac, tmpdir)
        print(temp_files)
        results = process_fragment(**temp_files)
        print(results)
        assert results != []


@pytest.mark.parametrize("yield_atac_h5ads", NON_HUMAN_MOUSE_H5ADS)
class TestOtherATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, yield_atac_fixture_data):
        self.atac = yield_atac_fixture_data
        self.validator = Validator()
        self.validator.adata = read_h5ad(FIXTURES_ROOT / self.atac.h5ad_file_name)


    def test_valid_other(self):
        print(self.atac)
        print(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("assay", ASSAYS)
    @pytest.mark.parametrize("organism_h_m", ["human","mouse"])
    def test_human_mouse_atac_fragment_w_other_h5ad(self, assay, tmpdir, organism_h_m):

       # ATAC fragment associated with a dataset with uns.organism_ontology_term_id on approved list but not human or mouse -> fail

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" != organism
        assert "NCBITaxon:9606" != organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        print(self.atac.fragment_df)
        self.atac.fragment_file_name = f"valid_{organism_h_m}_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        #print(temp_files)
        results = process_fragment(**temp_files)
        #print(results)
        assert results != []


    @pytest.mark.parametrize("assay",ASSAYS)
    @pytest.mark.parametrize("organism_h_m", ["human","mouse"])
    def test_human_mouse_atac_fragment_w_other_h5ad_matching_barcodes(self, assay, tmpdir, organism_h_m):

        # human/mouse ATAC fragment w/ other organism h5ads with matching barcodes -> fail

        organism = self.validator.adata.uns["organism_ontology_term_id"]
        assert "NCBITaxon:10090" != organism
        assert "NCBITaxon:9606" != organism
        self.atac.adata.obs["assay_ontology_term_id"] = assay
        self.atac.adata.obs["suspension_type"] = "nucleus"
        self.atac.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
        self.atac.adata.obs["is_primary_data"] = True
        print(self.atac.fragment_df)

        self.atac.adata.obs.reset_index(inplace=True)
        self.atac.fragment_df[3] = self.atac.adata.obs["obs_index"].sample(n=len(self.atac.fragment_df),replace=True).values

        self.atac.fragment_file_name = f"valid_{organism_h_m}_fragments.tsv.gz"
        temp_files = to_temp_files(self.atac, tmpdir)
        #print(temp_files)
        results = process_fragment(**temp_files)
        #print(results)
        assert results != []
