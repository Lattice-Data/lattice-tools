"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/728
https://github.com/chanzuckerberg/single-cell-curation/pull/753
"""

import numpy as np
import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
)

H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
]


@pytest.mark.parametrize("test_h5ads", H5ADS)
@pytest.mark.parametrize(
    "uns_key,empty_value,expected",
    (
        pytest.param('title', 'Here is my title', True, id='Normal string title'),
        pytest.param('batch_condition', ['donor_id'], True, id='Normal batch_condition list'),
        pytest.param('log1p', np.array([1, 2, 3, 4]), True, id='Normal author np.array'),
        pytest.param('log1p', True, True, id='Value set to True'),
        pytest.param('log1p', False, True, id='Value set to False'),
        pytest.param('log1p', None, True, id='Value set to None'),
    )
)
def test_uns_zero_length_passes(validator_with_adatas, uns_key, empty_value, expected):
    validator = validator_with_adatas
    validator.adata.uns[uns_key] = empty_value
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", H5ADS)
@pytest.mark.parametrize(
    "uns_key,empty_value,expected",
    (
        pytest.param('title', '', False, id='Empty string for title'),
        pytest.param('batch_condition', [], False, id='Empty list for batch_condition'),
        pytest.param('log1p', np.array([]), False, id='Empty numpy array, author key'),
        pytest.param('log1p', [], False, id='Empty list, author key'),
        pytest.param('log1p', "", False, id='Empty string, author key'),
        pytest.param('log1p', {}, False, id='Empty dictionary, author key'),
        pytest.param('log1p', (), False, id='Empty set, author key'),
    )
)
def test_uns_zero_length_fails(validator_with_adatas, uns_key, empty_value, expected):
    validator = validator_with_adatas
    validator.adata.uns[uns_key] = empty_value
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [f"ERROR: uns['{uns_key}'] cannot be an empty value."]


@pytest.mark.parametrize("test_h5ads", H5ADS)
@pytest.mark.parametrize(
    "uns_key,empty_value,expected",
    (
        pytest.param('default_embedding', [], False, id='Empty list for batch_condition'),
        pytest.param('QC_colors', np.array([]), False, id='Empty numpy array for {col}_colors'),
        pytest.param('X_approximate_distribution', '', False, id='Empty string for X_approximate_distribution'),
    )
)
def test_uns_zero_length_fails_many_errors(validator_with_adatas, uns_key, empty_value, expected):
    validator = validator_with_adatas
    validator.adata.uns[uns_key] = empty_value
    validator.validate_adata()
    assert validator.is_valid is expected
    assert f"ERROR: uns['{uns_key}'] cannot be an empty value." in validator.errors 


@pytest.mark.parametrize("test_h5ads", H5ADS)
@pytest.mark.parametrize(
    "bool_value,expected",
    (
        pytest.param(True, False, id='title set to True'),
        pytest.param(False, False, id='title set to False'),
        pytest.param(None, False, id='title set to None'),
    )
)
def test_uns_title_bool_none_fails(validator_with_adatas, bool_value, expected):
    validator = validator_with_adatas
    validator.adata.uns['title'] = bool_value
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{bool_value}' in 'uns['title']' is not valid, it must be a string."
    ]
