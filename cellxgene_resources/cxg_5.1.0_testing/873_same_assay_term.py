"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/873
https://github.com/chanzuckerberg/single-cell-curation/pull/885
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_all_visiums,
    validator_with_visium_some,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"


def test_mixed_assays_spatial(validator_with_spatial_adatas):
    validator = validator_with_spatial_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert validator.errors[0] == "ERROR: When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations must contain the same value."


def test_mixed_assays_slide_seq_is_single_false(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert validator.errors[0] == "ERROR: When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations must contain the same value."


def test_mixed_assays_visium_is_single_false(validator_with_all_visiums):
    validator = validator_with_all_visiums
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert validator.errors[1] == "ERROR: When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations must contain the same value."


def test_mixed_assays_non_spatial(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0010961")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0010961"
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() >= 2
    assert validator.is_valid is False
    assert validator.errors[1] == "ERROR: When obs['assay_ontology_term_id'] is either 'EFO:0010961' (Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations must contain the same value."
    
