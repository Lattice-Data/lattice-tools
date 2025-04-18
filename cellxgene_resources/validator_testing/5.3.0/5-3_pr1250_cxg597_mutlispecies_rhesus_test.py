"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1245
https://github.com/chanzuckerberg/single-cell-curation/pull/1250
"""

import pytest
import anndata as ad
from fixtures.create_fixtures import VAR_META_DF
from fixtures.valid_adatas import (
    MULTISPECIES_H5ADS,
    test_h5ads,
    validator_with_adatas
)


@pytest.mark.parametrize("test_h5ads", MULTISPECIES_H5ADS)
def test_multispecies_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


# Test invalid variation in feature ids
@pytest.mark.parametrize("test_h5ads", ["valid_rhesus.h5ad"])
class TestRhesusObsAndVar:
    @pytest.mark.parametrize(
        "ensembl_prefix_term",
        (
            pytest.param("PPP2R5A_ENSMMUG00000019771",id="Cxg unique identifier doesn't start with feature id"),
            pytest.param("NSMMUG00000019771",id="Not a valid prefix for feature id"),
            pytest.param(" ENSMMUG00000019771",id="Not a valid prefix for feature id"),
            pytest.param("PPP2R5A",id="Invalid feature id"),
            pytest.param("ENSMMUG00000019771.1",id="Invalid feature id"),
            pytest.param("ENSMMUG00000019771 ",id="Invalid feature id"),
            pytest.param("ENSMMUG-00000019771",id="Invalid feature id"),
            pytest.param("ENSMMUG_00000019771",id="Invalid feature id"),
        )
    )
    def test_feature_ids(self, validator_with_adatas,ensembl_prefix_term):
        validator = validator_with_adatas
        raw_adata = ad.AnnData(validator.adata.raw.X, dtype='float32', var=validator.adata.raw.var)
        raw_adata.var.reset_index(inplace=True)
        raw_adata.var.replace({'ensembl_id':{'ENSMMUG00000019771':ensembl_prefix_term}},inplace=True)
        raw_adata.var.set_index('ensembl_id',inplace=True)
        validator.adata.raw = raw_adata
        validator.adata.var.reset_index(inplace=True)
        validator.adata.var.replace({'ensembl_id':{'ENSMMUG00000019771':ensembl_prefix_term}},inplace=True)
        validator.adata.var.set_index('ensembl_id',inplace=True)
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: Could not infer organism from feature ID '{ensembl_prefix_term}' in 'var', "
            f"make sure it is a valid ID."
            ) in validator.errors
        assert (
            f"ERROR: Could not infer organism from feature ID '{ensembl_prefix_term}' in 'raw.var', "
            f"make sure it is a valid ID."
            ) in validator.errors

    @pytest.mark.parametrize(
        "organism_term",
        (
            pytest.param("NCBITaxon:6239", id="NCBITaxon term for worm not valid with rhesus obs"),
            pytest.param("NCBITaxon:7955", id="NCBITaxon term for zebrafish not valid with rhesus obs"),
            pytest.param("NCBITaxon:7227", id="NCBITaxon term for fly not valid with rhesus obs"),
            pytest.param("NCBITaxon:9606", id="NCBITaxon term for human not valid with rhesus obs"),
            pytest.param("NCBITaxon:10090", id="NCBITaxon term for mouse not valid with rhesus obs"),
            pytest.param("NCBITaxon:179238", id="NCBITaxon term for mouse descendant not valid with rhesus obs"),
        )
    )
    def test_ncbi_term(self, validator_with_adatas,organism_term):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = organism_term
        validator.validate_adata()
        assert not validator.is_valid

    # Test lack of dependencies between organism_ontology_term_id and gene ids - this shouldn't pass in the next schema bump.
    @pytest.mark.parametrize(
        "organism_term",
        (
            pytest.param("NCBITaxon:9483",id="valid NCBITaxon term for marmoset"),
            pytest.param("NCBITaxon:9595",id="valid NCBITaxon term for gorilla"),
            pytest.param("NCBITaxon:9986",id="valid NCBITaxon term for rabbit"),
            pytest.param("NCBITaxon:230741",id="valid NCBITaxon term for rabbit descendant"),
            pytest.param("NCBITaxon:10116",id="valid NCBITaxon term for rat"),
            pytest.param("NCBITaxon:947987",id="valid NCBITaxon term for rat descendant"),
            pytest.param("NCBITaxon:9598",id="valid NCBITaxon term for chimp"),
            pytest.param("NCBITaxon:756884",id="valid NCBITaxon term for chimp descendant"),
            pytest.param("NCBITaxon:9823",id="valid NCBITaxon term for domestic pig"),
            pytest.param("NCBITaxon:9825",id="valid NCBITaxon term for domestic pig descendant"),
            pytest.param("NCBITaxon:30608",id="valid NCBITaxon term for lemur"),
        )
    )
    def test_orthologs(self, validator_with_adatas,organism_term):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = organism_term
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # Test rhesus descendant
    @pytest.mark.parametrize(
        "organism_descendant",
        (
            pytest.param("NCBITaxon:1654737", id="Macaca mulatta brevicaudus NCBITaxon"),
            pytest.param("NCBITaxon:1449913", id="Macaca mulatta lasiotus NCBITaxon"),
            pytest.param("NCBITaxon:2008792", id="Macaca mulatta vestita NCBITaxon"),
        )
    )
    def test_descendant_organism_pass(self, validator_with_adatas, organism_descendant):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = organism_descendant
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []


# Test errors for ethnicity, tissue, development stage dependencies in human, mouse, and other model organisms with rhesus taxon term
@pytest.mark.parametrize("test_h5ads", ["valid_human.h5ad"])
class TestRhesusFixtureHumanTerms:
    @pytest.mark.parametrize(
        "rhesus_term,error_variable",
        (
            pytest.param("NCBITaxon:9544","unknown", id="self_reported_ethnicity_ontology_term_id invalid for rhesus"),
            pytest.param("NCBITaxon:9544","HANCESTRO:0005", id="self_reported_ethnicity_ontology_term_id invalid for rhesus"),
        )
    )
    def test_rhesus_ncbi_term_in_human_self_reported_ethnicity(self, validator_with_adatas,rhesus_term,error_variable):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = rhesus_term
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            f"ERROR: '{error_variable}' in 'self_reported_ethnicity_ontology_term_id' is not a valid value of "
            "'self_reported_ethnicity_ontology_term_id'. When 'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens),"
            " self_reported_ethnicity_ontology_term_id MUST be 'na'."
        ) in validator.errors

    @pytest.mark.parametrize(
        "rhesus_term,error_variable",
        (
            pytest.param("NCBITaxon:9544","HsapDv:0000258", id="development_stage_ontology_term_id invalid for rhesus"),
        )
    )
    def test_rhesus_ncbi_term_in_human_dev_stage(self, validator_with_adatas,rhesus_term,error_variable):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = rhesus_term
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            f"ERROR: '{error_variable}' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
            "When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant "
            "term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ["valid_mouse.h5ad"])
def test_rhesus_ncbi_term_in_mouse_dev_stage(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
    validator.validate_adata()
    assert not validator.is_valid
    assert(
        "ERROR: 'MmusDv:0000136' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
        "When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant "
        "term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ["valid_worm.h5ad"])
class TestWormFixtureRhesusTerms:
    def test_rhesus_ncbi_term_in_worm_tissue(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            "ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. If organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid UBERON term or a valid ZFA term. If organism is NCBITaxon:7227, "
            "it can be a valid UBERON term or a valid FBbt term. When tissue_type is cell culture, "
            "tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id."
        ) in validator.errors


    def test_rhesus_ncbi_term_in_worm_cl(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: cell_type_ontology_term_id must be a valid CL term. If organism is NCBITaxon:6239, it can be a valid CL term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid CL term or a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid CL term or a valid FBbt term."
        ) in validator.errors


    def test_rhesus_ncbi_term_in_worm_dev(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            "ERROR: 'WBls:0000544' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
            "When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, "
            "'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ["valid_fly.h5ad"])
class TestFlyFixtureRhesusTerms:
    def test_rhesus_ncbi_term_in_fly_tissue(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            "ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. If organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid UBERON term or a valid ZFA term. If organism is NCBITaxon:7227, "
            "it can be a valid UBERON term or a valid FBbt term. When tissue_type is cell culture, "
            "tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id."
        ) in validator.errors

    def test_rhesus_ncbi_term_in_fly_cl(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: cell_type_ontology_term_id must be a valid CL term. If organism is NCBITaxon:6239, it can be a valid CL term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid CL term or a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid CL term or a valid FBbt term."
        ) in validator.errors

    def test_rhesus_ncbi_term_in_fly_dev(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            f"ERROR: 'FBdv:00005370' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
            "When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, "
            "'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ) in validator.errors


@pytest.mark.parametrize("test_h5ads", ["valid_zebrafish.h5ad"])
class TestZebrafishFixtureRhesusTerms:
    def test_rhesus_ncbi_term_in_fish_tissue(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            "ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. If organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid UBERON term or a valid ZFA term. If organism is NCBITaxon:7227, "
            "it can be a valid UBERON term or a valid FBbt term. When tissue_type is cell culture, "
            "tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id."
        ) in validator.errors

    def test_rhesus_ncbi_term_in_fish_cl(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: cell_type_ontology_term_id must be a valid CL term. If organism is NCBITaxon:6239, it can be a valid CL term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid CL term or a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid CL term or a valid FBbt term."
        ) in validator.errors

    def test_rhesus_ncbi_term_in_fish_dev(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["organism_ontology_term_id"] = "NCBITaxon:9544"
        validator.validate_adata()
        assert not validator.is_valid
        assert(
            "ERROR: 'ZFS:0000033' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
            "When 'organism_ontology_term_id'-specific requirements are not defined in the schema definition, "
            "'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
        ) in validator.errors
