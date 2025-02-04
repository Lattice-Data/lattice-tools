"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1196
https://github.com/chanzuckerberg/single-cell-curation/pull/1230
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_visium,
    validator_with_non_visium_adatas
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"

def test_array_col_passes(validator_with_visium):
    validator = validator_with_visium
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []

@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_array_col_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_col"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_col'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_array_row_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_row'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_in_tissue_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "in_tissue"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['in_tissue'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_array_col_all_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[:], "array_col"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_col'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_array_row_all_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[:], "array_row"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_row'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_in_tissue_all_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[:], "in_tissue"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['in_tissue'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_multiple_contain_nan(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
    validator.adata.obs.loc[validator.adata.obs.index[0], "in_tissue"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[:2] == [
        f"ERROR: obs['array_row'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
        f"ERROR: obs['in_tissue'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


@pytest.mark.parametrize(
    "obs_col", ("array_col", "array_row", "in_tissue")
)
@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_all_contain_nan(validator_with_visium, obs_value, obs_col):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], obs_col] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[:3] == [
        f"ERROR: obs['{obs_col}'] cannot have missing or NaN values when obs['assay_ontology_term_id'] is "
        "a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


def test_array_col_missing(validator_with_visium):
    validator = validator_with_visium
    del validator.adata.obs['array_col']
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_col'] is required for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]

def test_array_row_missing(validator_with_visium):
    validator = validator_with_visium
    del validator.adata.obs['array_row']
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['array_row'] is required for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]

def test_in_tissue_missing(validator_with_visium):
    validator = validator_with_visium
    del validator.adata.obs['in_tissue']
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['in_tissue'] is required for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]

@pytest.mark.parametrize(
    "obs_value", (None, np.nan, float('nan'), float('NaN'), float('NAN'))
)

def test_obsm_visium_is_single_false_with_nans(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.uns["spatial"]["is_single"] = False
    del validator.adata.uns["spatial"][LIBRARY_ID]
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_col"] = obs_value
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
    validator.adata.obs.loc[validator.adata.obs.index[0], "in_tissue"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors[:3] == [
        f"ERROR: obs['array_col'] is only allowed for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.", 
        f"ERROR: obs['array_row'] is only allowed for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.", 
        f"ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] is a descendant of "
        "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]
