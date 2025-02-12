"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1247
https://github.com/chanzuckerberg/single-cell-curation/pull/1250
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_rat_adata
)


def test_passes(validator_with_rat_adata):
    validator = validator_with_rat_adata
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize(
    "genes", ('ENSDARG00000076014','WBGene00000003', 'FBgn0038542', 
        'ENSCJAG00000077381', 'ENSGGOG00000011263', 'ENSPTRG00000026694',
        'ENSMFAG00000033997', 'ENSMMUG00000052613', 'ENSSSCG00000005267',
        'ENSOCUG00000022999', 'ENSMICG00000010131')
)

# Should pass, test that ensembl_ids from other species are allowed
def test_different_species_ensembl(validator_with_rat_adata,genes):
    validator = validator_with_rat_adata
    validator.adata.var.rename(index={validator.adata.var.index[0]:genes}, inplace=True)
    validator.adata.raw.var.rename(index={validator.adata.raw.var.index[0]:genes}, inplace=True)
    validator.validate_adata()
    assert validator.is_valid

@pytest.mark.parametrize(
    "non_descendents", ('NCBITaxon:10114','NCBITaxon:83752', 'NCBITaxon:39107')
)

# Should fail, test if non descendent terms of Rattus norvegicus
def test_not_descendent(validator_with_rat_adata, non_descendents):
    validator = validator_with_rat_adata
    validator.adata.obs['organism_ontology_term_id'] = \
    validator.adata.obs['organism_ontology_term_id'].cat.add_categories([non_descendents])
    validator.adata.obs.loc[validator.adata.obs.index[0], "organism_ontology_term_id"] = non_descendents
    validator.validate_adata()
    assert validator.is_valid is False
    assert (
        f"ERROR: '{non_descendents}' in 'organism_ontology_term_id' is not an allowed "
        "term id. Only explicitly enumerated species are allowed. See Schema"
    ) in validator.errors


# Should pass, test if descendent terms of Rattus norvegicus
def test_descendent(validator_with_rat_adata):
    validator = validator_with_rat_adata
    validator.adata.obs['organism_ontology_term_id'] = \
    validator.adata.obs['organism_ontology_term_id'].cat.add_categories(['NCBITaxon:947987'])
    validator.adata.obs.loc[validator.adata.obs.index[0], "organism_ontology_term_id"] = 'NCBITaxon:947987'
    validator.validate_adata()
    assert validator.is_valid is True

# Should fail, test if death stage allowed for development_stage_ontology_term_id
def test_death_stage(validator_with_rat_adata):
    validator = validator_with_rat_adata
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
def test_life_cycle_stage(validator_with_rat_adata):
    validator = validator_with_rat_adata
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

# Should fail, test if invalid stages allowed for development_stage_ontology_term_id
def test_invalid_stage(validator_with_rat_adata, stages):
    validator = validator_with_rat_adata
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

