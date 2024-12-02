"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1114
https://github.com/chanzuckerberg/single-cell-curation/pull/1132/
"""

import numpy as np
import pytest
from fixtures.valid_adatas import (
    validator_with_all_visiums,
    validator_with_spatial_adatas,
    validator_with_non_spatial_adata,
    validator_with_all_adatas,
    validator_with_non_visium_adatas,
    validator_with_human_adatas
)

ANCESTORY_COLUMNS = [
    "genetic_ancestry_African",
    "genetic_ancestry_East_Asian",
    "genetic_ancestry_European",
    "genetic_ancestry_Indigenous_American",
    "genetic_ancestry_Oceanian",
    "genetic_ancestry_South_Asian",
]

def add_valid_ancestory(adata):
    fill_value = float("nan")
    for col in ANCESTORY_COLUMNS:
        adata.obs[col] = fill_value
        adata.obs[col] = adata.obs[col].astype("float")


def test_valid_with_nan_ancestory(validator_with_human_adatas):
    validator = validator_with_human_adatas
    add_valid_ancestory(validator.adata)
    validator.validate_adata()
    assert validator.is_valid


def test_ancestory_is_one(validator_with_human_adatas):
    validator = validator_with_human_adatas
    add_valid_ancestory(validator.adata)
    values = (0.25, 0.25, 0.25, 0.25, 0, 0)
    first_donor = validator.adata.obs["donor_id"].unique()[0]
    for value, ancestory_col in zip(values, ANCESTORY_COLUMNS):
        validator.adata.obs.loc[validator.adata.obs["donor_id"] == first_donor, ancestory_col] = value
    validator.validate_adata()
    assert validator.is_valid
