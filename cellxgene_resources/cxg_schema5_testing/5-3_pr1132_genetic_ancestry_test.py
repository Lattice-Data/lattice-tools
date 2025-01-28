"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1114
https://github.com/chanzuckerberg/single-cell-curation/pull/1132/
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    validator_with_all_visiums,
    validator_with_spatial_adatas,
    validator_with_non_spatial_adata,
    validator_with_all_adatas,
    validator_with_non_visium_adatas,
    validator_with_human_adatas,
    validator_human_adata,
    validator_mouse_adata
)

ANCESTRY_COLUMNS = [
    "genetic_ancestry_African",
    "genetic_ancestry_East_Asian",
    "genetic_ancestry_European",
    "genetic_ancestry_Indigenous_American",
    "genetic_ancestry_Oceanian",
    "genetic_ancestry_South_Asian",
]

def add_valid_nan_ancestry(adata, fill_value=float("nan")):
    for col in ANCESTRY_COLUMNS:
        adata.obs[col] = fill_value
        adata.obs[col] = adata.obs[col].astype("float")


def test_valid_with_nan_ancestry(validator_human_adata):
    validator = validator_human_adata
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_valid_mouse_nan_ancestry(validator_mouse_adata):
    validator = validator_mouse_adata
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


def test_one_row_different(validator_human_adata):
    validator = validator_human_adata
    values = (0.25, 0.25, 0.25, 0.25, 0, 0)
    for value, ancestry_col in zip(values, ANCESTRY_COLUMNS):
        validator.adata.obs.loc[validator.adata.obs.index[100], ancestry_col] = value
    donor_id = validator.adata.obs.loc[validator.adata.obs.index[100]]["donor_id"]
    validator.validate_adata()
    assert not validator.is_valid
    assert len(validator.errors) == 1
    assert (
        f"ERROR: obs rows with donor ids ['{donor_id}'] have invalid genetic_ancestry_* values. "
        "All observations with the same donor_id must contain the same genetic_ancestry_* values. "
        "If organism_ontolology_term_id is NOT 'NCBITaxon:9606' for Homo sapiens, then all geneticancestry "
        "values MUST be float('nan'). If organism_ontolology_term_id is 'NCBITaxon:9606' for Homo sapiens, "
        "then the value MUST be a float('nan') if unavailable; otherwise, the sum of all genetic_ancestry_* "
        "fields must be equal to 1.0"
    ) in validator.errors


def test_ancestry_is_one(validator_human_adata):
    validator = validator_human_adata
    values = (0.25, 0.25, 0.25, 0.25, 0, 0)
    first_donor = validator.adata.obs["donor_id"].unique()[0]
    for value, ancestry_col in zip(values, ANCESTRY_COLUMNS):
        validator.adata.obs.loc[validator.adata.obs["donor_id"] == first_donor, ancestry_col] = value
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize(
    "new_float",
    (
        pytest.param(0.26),
        pytest.param(0.250000001),
        pytest.param(0.25000001),
        pytest.param(0.2500001),
        pytest.param(0.250001),
        pytest.param(0.25001),
        pytest.param(0.00000001),
        pytest.param(0.0000001),
        pytest.param(0.000001),
        pytest.param(0.00001),
        pytest.param(0.0001),
    )
)
def test_donor_values_differ(validator_human_adata, new_float):
    validator = validator_human_adata
    values = (0.25, 0.25, 0.25, 0.25, 0, 0)
    first_donor = validator.adata.obs["donor_id"].unique()[0]
    first_donor_indices = validator.adata.obs[validator.adata.obs["donor_id"] == first_donor].index
    
    for value, ancestry_col in zip(values, ANCESTRY_COLUMNS):
        validator.adata.obs.loc[validator.adata.obs["donor_id"] == first_donor, ancestry_col] = value

    validator.adata.obs.loc[validator.adata.obs.index == first_donor_indices[1], ANCESTRY_COLUMNS[0]] = new_float
    
    validator.validate_adata()
    assert not validator.is_valid
    assert len(validator.errors) == 1
    assert (
        f"ERROR: obs rows with donor ids ['{first_donor}'] have invalid genetic_ancestry_* values. "
        "All observations with the same donor_id must contain the same genetic_ancestry_* values. "
        "If organism_ontolology_term_id is NOT 'NCBITaxon:9606' for Homo sapiens, then all geneticancestry "
        "values MUST be float('nan'). If organism_ontolology_term_id is 'NCBITaxon:9606' for Homo sapiens, "
        "then the value MUST be a float('nan') if unavailable; otherwise, the sum of all genetic_ancestry_* "
        "fields must be equal to 1.0"
    ) in validator.errors


@pytest.mark.parametrize(
    "nan_value,expected",
    (
        pytest.param(float("nan"), True, id="float nan"),
        pytest.param(None, True, id="None"),
        pytest.param(np.nan, True, id="np.nan"),
        pytest.param(pd.NA, False, id="pd.NA"),
    )
)
def test_different_nan_values(validator_human_adata, nan_value, expected):
    validator = validator_human_adata
    add_valid_nan_ancestry(validator.adata, nan_value)
    
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == []


@pytest.mark.parametrize(
    "value",
    (
        pytest.param(1.0000000001),
        pytest.param(1.000000001),
        pytest.param(1.00000001),
        pytest.param(1.0000001),
        pytest.param(1.000001),
        pytest.param(1.00001),
    )
)
def test_value_over_one(validator_human_adata, value):
    validator = validator_human_adata

    validator.adata.obs.loc[validator.adata.obs.index[0], ANCESTRY_COLUMNS[0]] = value
    first_donor = validator.adata.obs.loc[validator.adata.obs.index[0]]["donor_id"]
        
    validator.validate_adata()
    assert not validator.is_valid
    assert len(validator.errors) == 2
    errors = [
        (
            f"ERROR: obs rows with donor ids ['{first_donor}'] have invalid genetic_ancestry_* values. "
            "All observations with the same donor_id must contain the same genetic_ancestry_* values. "
            "If organism_ontolology_term_id is NOT 'NCBITaxon:9606' for Homo sapiens, then all geneticancestry "
            "values MUST be float('nan'). If organism_ontolology_term_id is 'NCBITaxon:9606' for Homo sapiens, "
            "then the value MUST be a float('nan') if unavailable; otherwise, the sum of all genetic_ancestry_* "
            "fields must be equal to 1.0"
        ),
        (
            f"ERROR: Column 'genetic_ancestry_African' in obs contains invalid values: [{value}]. "
            "Valid values are floats between 0 and 1 or float('nan')."
        )
    ]
    for error in errors:
        assert error in validator.errors
