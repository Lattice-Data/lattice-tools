"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/830
https://github.com/chanzuckerberg/single-cell-curation/pull/864
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"


def test_change_unknown_to_cell_type(validator_with_visium):
    validator = validator_with_visium
    validator.adata.obs.replace({"cell_type_ontology_term_id": {"unknown": "CL:0002601"}}, inplace=True)
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
    ]


def test_change_all_in_tissue_to_zero(validator_with_visium):
    validator = validator_with_visium
    validator.adata.obs["in_tissue"] = 0
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
    ]


def test_one_random_in_tissue_zero(validator_with_visium):
    validator = validator_with_visium
    validator.adata.obs.replace({"cell_type_ontology_term_id": {"unknown": "CL:0002601"}}, inplace=True)
    validator.adata.obs["in_tissue"] = 1
    random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
    validator.adata.obs.loc[validator.adata.obs.index[random_index], "in_tissue"] = 0
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
    ]


def test_slide_seq_in_tissue(validator_with_slide_seq_adatas):
    validator = validator_with_slide_seq_adatas
    validator.adata.obs["in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


def test_non_spatial_in_tissue(validator_with_non_spatial_adata):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


def test_change_in_tissue_to_zero(validator_with_visium):
    validator = validator_with_visium
    validator.adata.obs.replace({"in_tissue": {1: 0}}, inplace=True)
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
    ]


@pytest.mark.parametrize(
    "change_dict", (
        {"cell_type_ontology_term_id": {"unknown": "CL:0002601"}},
        {"in_tissue": {1: 0}}
    )
)
def test_is_single_false(validator_with_visium, change_dict):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["is_single"] = False
    validator.adata.obs.replace(change_dict, inplace=True)
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: uns['spatial'][library_id] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
        "ERROR: obs['array_col'] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
        "ERROR: obs['array_row'] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
        "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


# should pass, rare but sometimes cell type is unknown
def test_in_tissue_one_unknown(validator_with_visium):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs["cell_type_ontology_term_id"] == "unknown", "in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
