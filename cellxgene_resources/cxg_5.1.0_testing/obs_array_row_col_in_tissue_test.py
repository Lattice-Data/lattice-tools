"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/829
https://github.com/chanzuckerberg/single-cell-curation/pull/861
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_all_adatas,
    validator_with_slide_seq_adatas,
    validator_with_spatial_adatas,
    validator_with_visium,
    validator_with_non_spatial_adata,
)

LIBRARY_ID = "spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna"

def test_array_col_passes(validator_with_visium):
    validator = validator_with_visium
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize(
    "obs_value", (None, -1, -1.0, 128, 128.0, "string", True, False)
)
def test_array_col_values(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_col"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False


@pytest.mark.parametrize(
    "obs_value", (None, -1, -1.0, 78, 78.0, "string", True, False)
)
def test_array_row_values(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False


@pytest.mark.parametrize(
    "obs_value", (None, -1, -1.0, 2, 2.0, "string", True, False)
)
def test_in_tissue_values(validator_with_visium, obs_value):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], "in_tissue"] = obs_value
    validator.validate_adata()
    assert validator.is_valid is False


@pytest.mark.parametrize(
    "obs_col", ("array_row", "array_col", "in_tissue")
)
def test_array_col_names(validator_with_visium, obs_col):
    validator = validator_with_visium
    del validator.adata.obs[obs_col]
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['{obs_col}'] is required for obs['assay_ontology_term_id'] 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


# in-range floats converted to int and pass
@pytest.mark.parametrize(
    "obs_col", ("array_row", "array_col", "in_tissue")
)
def test_obs_float_values_pass(validator_with_visium, obs_col):
    validator = validator_with_visium
    validator.adata.obs.loc[validator.adata.obs.index[0], obs_col] = 1.0
    validator.validate_adata()
    assert validator.is_valid is True
    assert validator.errors == []


@pytest.mark.parametrize(
    "obs_col", ("array_row", "array_col", "in_tissue")
)
def test_array_col_values_is_single_false(validator_with_visium, obs_col):
    validator = validator_with_visium
    del validator.adata.obs[obs_col]
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.is_valid is False
    assert len(validator.errors) == 3


@pytest.mark.parametrize(
    "obs_col,obs_type", 
        (
            ("array_row", "float32"),
            ("array_row", "float64"),
            ("array_row", "str"),
            ("array_row", "object"),
            ("array_row", "bool"),
            ("array_col", "float32"),
            ("array_col", "float64"),
            ("array_col", "str"),
            ("array_col", "object"),
            ("array_col", "bool"),
            ("in_tissue", "float32"),
            ("in_tissue", "float64"),
            ("in_tissue", "str"),
            ("in_tissue", "object"),
            ("in_tissue", "bool"),
    )
)
def test_non_int_col_row_in_tissue(validator_with_visium, obs_col, obs_type):
    validator = validator_with_visium
    validator.adata.obs[obs_col] = validator.adata.obs[obs_col].astype(obs_type)   
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        f"ERROR: obs['{obs_col}'] must be of int type, it is {validator.adata.obs[obs_col].dtype}."
    ]
