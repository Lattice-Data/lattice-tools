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

import pytest
import anndata as ad
import pandas as pd
from cellxgene_schema import atac_seq
import tempfile
import shutil

from fixtures.valid_adatas import (
    FIXTURES_ROOT,
    NON_SPATIAL_H5ADS,
    validator_with_adatas
)


HUMAN_H5ADs = [file for file in NON_SPATIAL_H5ADS if "human" in file]
MOUSE_H5ADs = [file for file in NON_SPATIAL_H5ADS if "mouse" in file]
NON_HUMAN_MOUSE_H5ADS = [file for file in NON_SPATIAL_H5ADS if "human" not in file and "mouse" not in file]

ASSAYS = ["EFO:0030007", "EFO:0008925", "EFO:0008904", "EFO:0022045", "EFO:0030059"]


@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADs)
class TestHumanATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas, test_h5ads):
        self.validator = validator_with_adatas
        self.h5ad_file_name = str(FIXTURES_ROOT / test_h5ads)
        self.fragement_file_name = str(FIXTURES_ROOT / test_h5ads).replace(".h5ad","_fragments.tsv.gz")


    def test_human_fixture_pass(self):

        self.validator.validate_adata()
        assert self.validator.is_valid


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_human_h5ad(self,assay):

        # human ATAC fragment w/ human h5ad -> pass

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:9606" == organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            parquet_file = atac_seq.convert_to_parquet(fragment_file = self.fragement_file_name,tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id=organism)
            assert result == []

        finally:
            shutil.rmtree(temp_dir)


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_human_h5ad(self,assay):

        # mouse ATAC fragment w/ human h5ad -> fail

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:9606" == organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            self.fragement_file_name = self.fragement_file_name.replace("human","mouse")
            assert self.fragement_file_name == str(FIXTURES_ROOT / "valid_mouse_fragments.tsv.gz")
            parquet_file = atac_seq.convert_to_parquet(fragment_file=self.fragement_file_name, tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id=organism)
            assert "Errors found in Fragment and/or Anndata file" in result
            assert "Barcodes don't match anndata.obs.index" in result
            assert any(f"Chromosomes in the fragment do not match the organism({organism})." in e for e in result)

        finally:
            shutil.rmtree(temp_dir)


@pytest.mark.parametrize("test_h5ads", MOUSE_H5ADs)
class TestMouseATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas, test_h5ads):
        self.validator = validator_with_adatas
        self.h5ad_file_name = str(FIXTURES_ROOT / test_h5ads)
        self.fragement_file_name = str(FIXTURES_ROOT / test_h5ads).replace(".h5ad","_fragments.tsv.gz")


    def test_mouse_fixture_pass(self):

        self.validator.validate_adata()
        assert self.validator.is_valid


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_mouse_atac_fragment_w_mouse_h5ad(self,assay):

        # mouse ATAC fragment w/ mouse h5ad -> pass

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:10090" == organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            parquet_file = atac_seq.convert_to_parquet(fragment_file = self.fragement_file_name,tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id=organism)
            assert result == []

        finally:
            shutil.rmtree(temp_dir)


    @pytest.mark.parametrize("assay", ASSAYS)
    def test_human_atac_fragment_w_mouse_h5ad(self,assay):

        # human ATAC fragment w/ mouse h5ad -> fail

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:10090" == organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            self.fragement_file_name = self.fragement_file_name.replace("mouse","human")
            assert self.fragement_file_name == str(FIXTURES_ROOT / "valid_human_fragments.tsv.gz")
            parquet_file = atac_seq.convert_to_parquet(fragment_file=self.fragement_file_name, tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id=organism)
            assert "Errors found in Fragment and/or Anndata file" in result
            assert "Barcodes don't match anndata.obs.index" in result
            assert any(f"Chromosomes in the fragment do not match the organism({organism})." in e for e in result)

        finally:
            shutil.rmtree(temp_dir)


@pytest.mark.parametrize("test_h5ads", NON_HUMAN_MOUSE_H5ADS)
class TestNonHumanMouseATACValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas, test_h5ads):
        self.validator = validator_with_adatas
        self.h5ad_file_name = str(FIXTURES_ROOT / test_h5ads)
        self.fragement_file_name = str(FIXTURES_ROOT / test_h5ads).replace(".h5ad","_fragments.tsv.gz")


    def test_other_fixture_pass(self):

        self.validator.validate_adata()
        assert self.validator.is_valid


    @pytest.mark.parametrize("assay",ASSAYS)
    def test_human_atac_fragment_w_other_h5ad(self,assay):

        # human ATAC fragment w/ other organism h5ads -> fail

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:9606" != organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            self.fragement_file_name = str(FIXTURES_ROOT / "valid_human_fragments.tsv.gz")
            parquet_file = atac_seq.convert_to_parquet(fragment_file=self.fragement_file_name, tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id="NCBITaxon:9606")
            assert "Errors found in Fragment and/or Anndata file" in result
            assert "Barcodes don't match anndata.obs.index" in result

        finally:
            shutil.rmtree(temp_dir)


    @pytest.mark.parametrize("assay",ASSAYS)
    def test_mouse_atac_fragment_w_other_h5ad(self,assay):

        # mouse ATAC fragment w/ other organism h5ads -> fail

        temp_dir = tempfile.mkdtemp()
        try:
            organism = self.validator.adata.uns["organism_ontology_term_id"]
            assert "NCBITaxon:10090" != organism
            self.validator.adata.obs["assay_ontology_term_id"] = assay
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
            self.validator.adata.obs["is_primary_data"] = True
            self.fragement_file_name = str(FIXTURES_ROOT / "valid_mouse_fragments.tsv.gz")
            parquet_file = atac_seq.convert_to_parquet(fragment_file=self.fragement_file_name, tempdir=temp_dir)
            result = atac_seq.validate_anndata_with_fragment(parquet_file, self.h5ad_file_name, organism_ontology_term_id="NCBITaxon:10090")
            assert "Errors found in Fragment and/or Anndata file" in result
            assert "Barcodes don't match anndata.obs.index" in result

        finally:
            shutil.rmtree(temp_dir)