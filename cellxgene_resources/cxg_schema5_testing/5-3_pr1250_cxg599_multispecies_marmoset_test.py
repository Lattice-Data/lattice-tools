"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1218
https://github.com/chanzuckerberg/single-cell-curation/pull/1250
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_marmoset_adata
)


def test_passes(validator_with_marmoset_adata):
    validator = validator_with_marmoset_adata
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize(
    "genes", ('ENSDARG00000076014','WBGene00000003', 'FBgn0038542', 
        'ENSGGOG00000033645', 'ENSRNOG00000029315', 'ENSPTRG00000026694',
        'ENSMFAG00000033997', 'ENSMMUG00000052613', 'ENSSSCG00000005267',
        'ENSOCUG00000022999', 'ENSMICG00000010131')
)

# Should pass, test that ensembl_ids from other species are allowed
def test_different_species_ensembl(validator_with_marmoset_adata,genes):
    validator = validator_with_marmoset_adata
    validator.adata.var.rename(index={validator.adata.var.index[0]:genes}, inplace=True)
    validator.adata.raw.var.rename(index={validator.adata.raw.var.index[0]:genes}, inplace=True)
    validator.validate_adata()
    assert validator.is_valid

@pytest.mark.parametrize(
    "non_descendents", ('NCBITaxon:1965096','NCBITaxon:1476500','NCBITaxon:57378')
)

# Should fail, test if non descendent terms of Callithrix jacchus
def test_not_descendent(validator_with_marmoset_adata, non_descendents):
    validator = validator_with_marmoset_adata
    validator.adata.obs['organism_ontology_term_id'] = \
    validator.adata.obs['organism_ontology_term_id'].cat.add_categories([non_descendents])
    validator.adata.obs.loc[validator.adata.obs.index[0], "organism_ontology_term_id"] = non_descendents
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: '{non_descendents}' in 'organism_ontology_term_id' is not an allowed "
        "term id. Only explicitly enumerated species are allowed. See Schema"
    ) in validator.errors

# Should fail, test if death stage allowed for development_stage_ontology_term_id
def test_death_stage(validator_with_marmoset_adata):
    validator = validator_with_marmoset_adata
    validator.adata.obs['development_stage_ontology_term_id'] = \
    validator.adata.obs['development_stage_ontology_term_id'].cat.add_categories(['UBERON:0000071'])
    validator.adata.obs.loc[validator.adata.obs.index[0], "development_stage_ontology_term_id"] = 'UBERON:0000071'
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: 'UBERON:0000071' in 'development_stage_ontology_term_id' is not "
        "allowed. When 'organism_ontology_term_id'-specific requirements are not "
        "defined in the schema definition, 'development_stage_ontology_term_id' "
        "MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
    ) in validator.errors

# Should fail, test if life cycle stage allowed for development_stage_ontology_term_id
def test_life_cycle_stage(validator_with_marmoset_adata):
    validator = validator_with_marmoset_adata
    validator.adata.obs['development_stage_ontology_term_id'] = \
    validator.adata.obs['development_stage_ontology_term_id'].cat.add_categories(['UBERON:0000105'])
    validator.adata.obs.loc[validator.adata.obs.index[0], "development_stage_ontology_term_id"] = 'UBERON:0000105'
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: 'UBERON:0000105' in 'development_stage_ontology_term_id' is not an "
        "allowed term id. When 'organism_ontology_term_id'-specific requirements are not "
        "defined in the schema definition, 'development_stage_ontology_term_id' "
        "MUST be a descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or unknown."
    ) in validator.errors

@pytest.mark.parametrize(
    "stages", ('WBls:0000811','ZFS:0000008','FBdv:00007088')
)

# Should fail, test if invalid (fly/worm/zebrafish) stages allowed for development_stage_ontology_term_id
def test_invalid_stage(validator_with_marmoset_adata, stages):
    validator = validator_with_marmoset_adata
    validator.adata.obs['development_stage_ontology_term_id'] = \
    validator.adata.obs['development_stage_ontology_term_id'].cat.add_categories([stages])
    validator.adata.obs.loc[validator.adata.obs.index[0], "development_stage_ontology_term_id"] = stages
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: '{stages}' in 'development_stage_ontology_term_id' is not a "
         "valid ontology term id of 'UBERON'. When "
         "'organism_ontology_term_id'-specific requirements are not defined in the "
         "schema definition, 'development_stage_ontology_term_id' MUST be a "
         "descendant term id of 'UBERON:0000105' excluding 'UBERON:0000071', or "
         "unknown."
    ) in validator.errors

@pytest.mark.parametrize(
    "ctypes", ('WBbt:0003672','ZFA:0009277','FBbt:00007436')
)
# Should fail, test if invalid (fly/worm/zebrafish) celltypes allowed for cell_type_ontology_term_id
def test_invalid_celltype(validator_with_marmoset_adata, ctypes):
    validator = validator_with_marmoset_adata
    validator.adata.obs['cell_type_ontology_term_id'] = \
    validator.adata.obs['cell_type_ontology_term_id'].cat.add_categories([ctypes])
    validator.adata.obs.loc[validator.adata.obs.index[0], "cell_type_ontology_term_id"] = ctypes
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: cell_type_ontology_term_id must be a valid CL term. "
        "If organism is NCBITaxon:6239, it can be a valid CL term or "
        "a valid WBbt term. If organism is NCBITaxon:7955, it can be "
        "a valid CL term or a valid ZFA term. If organism is "
        "NCBITaxon:7227, it can be a valid CL term or a valid FBbt term."
    ) in validator.errors



@pytest.mark.parametrize(
    "tissues", ('WBbt:0005742','ZFA:0001284','FBbt:00004895')
)

# Should fail, test if invalid (fly/worm/zebrafish) tissue types allowed for tissue_ontology_term_id
def test_invalid_tissue(validator_with_marmoset_adata, tissues):
    validator = validator_with_marmoset_adata
    validator.adata.obs['tissue_ontology_term_id'] = \
    validator.adata.obs['tissue_ontology_term_id'].cat.add_categories([tissues])
    validator.adata.obs.loc[validator.adata.obs.index[0], "tissue_ontology_term_id"] = tissues
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id "
        "must be a valid UBERON term. If organism is NCBITaxon:6239, it can be a "
        "valid UBERON term or a valid WBbt term. If organism is NCBITaxon:7955, it "
        "can be a valid UBERON term or a valid ZFA term. If organism is NCBITaxon:7227, "
        "it can be a valid UBERON term or a valid FBbt term. When tissue_type is cell "
        "culture, tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id."
    ) in validator.errors