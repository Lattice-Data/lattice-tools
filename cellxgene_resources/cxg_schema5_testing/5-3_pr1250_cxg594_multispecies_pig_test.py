"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1222
https://github.com/chanzuckerberg/single-cell-curation/pull/1250
"""


import pytest
from fixtures.valid_adatas import validator_with_pig_adata


# Pass test that pig fixture is valid
def test_valid_adata(validator_with_pig_adata):
	validator = validator_with_pig_adata
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
def test_other_species_ensembl(validator_with_pig_adata, other_genes):
	validator = validator_with_pig_adata
	validator.adata.var.rename(index={validator.adata.var.index[0]:other_genes}, inplace=True)
	validator.adata.raw.var.rename(index={validator.adata.raw.var.index[0]:other_genes}, inplace=True)
	validator.validate_adata()
	assert validator.is_valid


# Pass test for descendent terms of Sus scrofa
@pytest.mark.parametrize(
	"organism_descendent",
	(
		pytest.param("NCBITaxon:1170810", id="Sus scrofa affinis NCBITaxon"),
		pytest.param("NCBITaxon:310260", id="Sus scrofa andamensis NCBITaxon"),
		pytest.param("NCBITaxon:291050", id="Sus scrofa coreanus NCBITaxon"),
		pytest.param("NCBITaxon:309913", id="Sus scrofa cristatus NCBITaxon"),
		pytest.param("NCBITaxon:9825", id="Sus scrofa domesticus NCBITaxon"),
		pytest.param("NCBITaxon:375578", id="Sus scrofa leucomystax NCBITaxon"),
		pytest.param("NCBITaxon:1611879", id="Sus scrofa lybicus NCBITaxon"),
		pytest.param("NCBITaxon:1611880", id="Sus scrofa majori NCBITaxon"),
		pytest.param("NCBITaxon:1611878", id="Sus scrofa meridionalis NCBITaxon"),
		pytest.param("NCBITaxon:309914", id="Sus scrofa papuensis NCBITaxon"),
		pytest.param("NCBITaxon:375579", id="Sus scrofa riukiuanus NCBITaxon"),
		pytest.param("NCBITaxon:415978", id="Sus scrofa scrofa NCBITaxon"),
		pytest.param("NCBITaxon:310261", id="Sus scrofa taivanus NCBITaxon"),
		pytest.param("NCBITaxon:490583", id="Sus scrofa ussuricus NCBITaxon"),
		pytest.param("NCBITaxon:2485929", id="Sus scrofa vittatus NCBITaxon")
	)
)
def test_descendent_organism_pass(validator_with_pig_adata, organism_descendent):
	validator = validator_with_pig_adata
	validator.adata.obs['organism_ontology_term_id'] = organism_descendent
	validator.validate_adata()
	assert validator.is_valid


# Fail test for invalid organism ontology
@pytest.mark.parametrize(
	"invalid_organism",
	(
		pytest.param("NCBITaxon:9822", id="Sus NCBITaxon"),
		pytest.param("NCBITaxon:35497", id="Suina NCBITaxon")
	)
)
def test_descendent_organism_fail(validator_with_pig_adata, invalid_organism):
	validator = validator_with_pig_adata
	validator.adata.obs['organism_ontology_term_id'] = invalid_organism
	validator.validate_adata()
	assert not validator.is_valid
	assert (
		f"ERROR: '{invalid_organism}' in 'organism_ontology_term_id' is not an allowed "
		"term id. Only explicitly enumerated species are allowed. See Schema"
	) in validator.errors


# Fail test for model organism dev stage ontology terms
@pytest.mark.parametrize(
	"invalid_dev_stage",
	(
		pytest.param("WBls:0000125", id="C. elegans dev stage not valid"),
		pytest.param("FBdv:00005336", id="Drosophila dev stage not valid"),
		pytest.param("ZFS:0000022", id="Zebrafish dev stage not valid")
	)
)
def test_invalid_dev(validator_with_pig_adata, invalid_dev_stage):
	validator = validator_with_pig_adata
	validator.adata.obs['development_stage_ontology_term_id'] = invalid_dev_stage
	validator.validate_adata()
	assert not validator.is_valid
	assert (
		f"ERROR: '{invalid_dev_stage}' in 'development_stage_ontology_term_id' is not a valid "
		"ontology term id of 'UBERON'. When 'organism_ontology_term_id'-specific requirements "
		"are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be "
		"a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
	) in validator.errors


# Fail test for invalid dev stage ontology terms
@pytest.mark.parametrize(
	"invalid_uberon_stage, error_wording",
	(
		pytest.param("UBERON:0000105", "an allowed term id", id="'life cycle stage' UBERON term is not valid"),
		pytest.param("UBERON:0000071", "allowed", id="'death stage' UBERON term is not valid")
	)
)
def test_invalid_dev_uberon(validator_with_pig_adata, invalid_uberon_stage, error_wording):
	validator = validator_with_pig_adata
	validator.adata.obs['development_stage_ontology_term_id'] = invalid_uberon_stage
	validator.validate_adata()
	assert not validator.is_valid
	assert (
		f"ERROR: '{invalid_uberon_stage}' in 'development_stage_ontology_term_id' is not {error_wording}. "
		"When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, "
		"'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding "
		"'UBERON:0000071', or unknown."
	) in validator.errors


# Fail test for model organism cell type ontology terms
@pytest.mark.parametrize(
	"invalid_celltype",
	(
		pytest.param("WBbt:0003672", id="'epithelial cell' for C. elegans is not valid"),
		pytest.param("ZFA:0009277", id="'acinar cell' for Zebrafish is not valid"),
		pytest.param("FBbt:00007436", id = "'endocrine neuron' for Drosophila")
	)
)
def test_invalid_celltype(validator_with_pig_adata, invalid_celltype):
	validator = validator_with_pig_adata
	validator.adata.obs['cell_type_ontology_term_id'] = invalid_celltype
	validator.validate_adata()
	assert not validator.is_valid
	assert (
		f"ERROR: cell_type_ontology_term_id must be a valid CL term. If organism is NCBITaxon:6239, it can "
		"be a valid CL term or a valid WBbt term. If organism is NCBITaxon:7955, it can be a valid CL term or "
		"a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid CL term or a valid FBbt term."
	) in validator.errors


# Fail test for model organism tissue ontology terms
@pytest.mark.parametrize(
	"invalid_tissue",
	(
		pytest.param("WBbt:0005742", id="'body wall' for C. elegans is not valid"),
		pytest.param("ZFA:0001284", id="'optic fissure' for Zebrafish is not valid"),
		pytest.param("FBbt:00004895", id="'basal stalk' for Drosophila is not valid")
	)
)
def test_invalid_tissue(validator_with_pig_adata, invalid_tissue):
	validator = validator_with_pig_adata
	validator.adata.obs['tissue_ontology_term_id'] = invalid_tissue
	validator.validate_adata()
	assert not validator.is_valid
	assert (
		f"ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. If "
		"organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. If organism is NCBITaxon:7955, "
		"it can be a valid UBERON term or a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid UBERON term "
		"or a valid FBbt term. When tissue_type is cell culture, tissue_ontology_term_id must follow the validation rules "
		"for cell_type_ontology_term_id."
	) in validator.errors
