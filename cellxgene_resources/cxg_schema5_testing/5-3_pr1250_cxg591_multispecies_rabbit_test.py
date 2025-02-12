"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1246
"""

import pytest
from fixtures.valid_adatas import validator_with_rabbit_adata


# Pass test that rabbit fixture is valid
def test_valid_adata(validator_with_rabbit_adata):
	validator = validator_with_rabbit_adata
	validator.validate_adata()
	assert validator.is_valid
	assert validator.errors == []


# Pass test if there is are genes from other organisms
@pytest.mark.parametrize(
	"other_genes",
	(
		pytest.param("ENSG00000187730", id="Human Ensembl is valid"),
		pytest.param("ENSDARG00000095324", id="Zebrafish is valid")
	)
)
def test_other_species_ensembl(validator_with_rabbit_adata, other_genes):
	validator = validator_with_rabbit_adata
	validator.adata.var.rename(index={validator.adata.var.index[0]:other_genes}, inplace=True)
	validator.adata.raw.var.rename(index={validator.adata.raw.var.index[0]:other_genes}, inplace=True)
	validator.validate_adata()
	assert validator.is_valid


# Pass test for descendent terms of Oryctolagus cuniculus
@pytest.mark.parametrize(
	"organism_descendent",
	(
		pytest.param("NCBITaxon:230741", id="Oryctolagus cuniculus algirus NCBITaxon"),
		pytest.param("NCBITaxon:568996", id="Oryctolagus cuniculus cuniculus NCBITaxon")
	)
)
def test_descendent_organism(validator_with_rabbit_adata, organism_descendent):
	validator = validator_with_rabbit_adata
	validator.adata.obs['organism_ontology_term_id'] = organism_descendent
	validator.validate_adata()
	assert validator.is_valid


# Fail test for invalid organism ontology
@pytest.mark.parametrize(
	"invalid_organism",
	(
		pytest.param("NCBITaxon:9984", id="Oryctolagus NCBITaxon"),
		pytest.param("NCBITaxon:2641529", id="unclassified Oryctolagus NCBITaxon")
	)
)
def test_descendent_organism(validator_with_rabbit_adata, invalid_organism):
	validator = validator_with_rabbit_adata
	validator.adata.obs['organism_ontology_term_id'] = invalid_organism
	validator.validate_adata()
	assert not validator.is_valid

# Dev stages
# Organism ontology