"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1245
https://github.com/chanzuckerberg/single-cell-curation/pull/1250
"""

import pytest
import anndata as ad
from fixtures.create_fixtures import VAR_META_DF
from fixtures.valid_adatas import (
    validator_with_multispecies_adatas,
    validator_with_lemur_adata
)


def test_multispecies_valid(validator_with_multispecies_adatas):
    validator = validator_with_multispecies_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_var_meta_df():
    print(VAR_META_DF)
    assert VAR_META_DF.shape[1] == 2
    assert "ensembl_id" in VAR_META_DF.columns
    assert "organism" in VAR_META_DF.columns


@pytest.mark.parametrize(
    "ensembl_prefix_term,error_variable",
    (
        pytest.param("LRTM2_ENSMICG00000042507","make sure it is a valid ID",id="Cxg unique identifier doesn't start with feature id"),
        pytest.param("NSMICG00000042507","make sure it is a valid ID",id="Not a valid prefix for feature id"),
        pytest.param(" ENSMICG00000042507","make sure it is a valid ID",id="Not a valid prefix for feature id"),
        pytest.param("LRTM2","is not a valid feature ID",id="Invalid feature id"),
        pytest.param("ENSMICT00000068998.1","is not a valid feature ID",id="Invalid feature id"),
        pytest.param("ENSMICG00000042507 ","is not a valid feature ID",id="Invalid feature id"),
        pytest.param("ENSMICG-00000042507","is not a valid feature ID",id="Invalid feature id"),
        pytest.param("ENSMICG_00000042507","is not a valid feature ID",id="Invalid feature id"),
    )
)
def test_feature_ids(validator_with_lemur_adata,ensembl_prefix_term,error_variable):
    validator = validator_with_lemur_adata
    raw_adata = ad.AnnData(validator.adata.raw.X, dtype='float32', var=validator.adata.raw.var)
    raw_adata.var.reset_index(inplace=True)
    raw_adata.var.replace({'ensembl_id':{'ENSMICG00000042507':ensembl_prefix_term}},inplace=True)
    raw_adata.var.set_index('ensembl_id',inplace=True)
    validator.adata.raw = raw_adata
    validator.adata.var.reset_index(inplace=True)
    validator.adata.var.replace({'ensembl_id':{'ENSMICG00000042507':ensembl_prefix_term}},inplace=True)
    validator.adata.var.set_index('ensembl_id',inplace=True)
    validator.validate_adata()
    print(validator.adata.var)
    assert not validator.is_valid
    assert (
        f"ERROR: 'Could not infer organism from feature ID '{ensembl_prefix_term}' in 'var',"
        f"{error_variable}."
        )


@pytest.mark.parametrize(
    "organism_term,error_variable",
    (
        pytest.param("NCBITaxon:6239","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for worm"),
        pytest.param("NCBITaxon:9483","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for marmoset"),
        pytest.param("NCBITaxon:7955","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for zebrafish"),
        pytest.param("NCBITaxon:7227","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for fly"),
        pytest.param("NCBITaxon:9595","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for gorilla"),
        pytest.param("NCBITaxon:9606","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for human"),
        pytest.param("NCBITaxon:10090","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for mouse"),
        pytest.param("NCBITaxon:179238","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for mouse descendant"),
        pytest.param("NCBITaxon:9986","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rabbit"),
        pytest.param("NCBITaxon:230741","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rabbit descendant"),
        pytest.param("NCBITaxon:10116","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rat"),
        pytest.param("NCBITaxon:947987","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rat descendant"),
        pytest.param("NCBITaxon:9598","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for chimp"),
        pytest.param("NCBITaxon:756884","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for chimp descendant"),
        pytest.param("NCBITaxon:9823","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for domestic pig"),
        pytest.param("NCBITaxon:9825","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for domestic pig descendant"),
        pytest.param("NCBITaxon:9544","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rhesus"),
        pytest.param("NCBITaxon:1654737","is not a valid ontology term id of 'NCBITaxon'", id="NCBITaxon term for rhesus descendant")
    )
)
def test_ncbi_term(validator_with_lemur_adata,organism_term,error_variable):
    validator = validator_with_lemur_adata
    validator.adata.obs["organism_ontology_term_id"] = organism_term
    validator.validate_adata()
    assert not validator.is_valid
    assert(
        f"ERROR: {organism_term} in 'organism_ontology_term_id' {error_variable}."
        "Only explicitly enumerated species are allowed. See Schema"
    )
