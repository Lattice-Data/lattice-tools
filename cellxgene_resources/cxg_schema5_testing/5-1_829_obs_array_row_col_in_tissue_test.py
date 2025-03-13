"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/829
https://github.com/chanzuckerberg/single-cell-curation/pull/861
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    ALL_H5ADS,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
NON_VISIUM_H5ADS = [file for file in ALL_H5ADS if not "visium" in file]


@pytest.mark.parametrize("test_h5ads", NON_VISIUM_H5ADS)
class TestObsVisiumColsNonVisium:
    @pytest.mark.parametrize(
        "obs_col", ("array_col", "array_row", "in_tissue")
    )
    def test_non_visium_obs_cols(self, validator_with_adatas, obs_col):
        validator = validator_with_adatas
        validator.adata.obs[obs_col] = np.random.randint(0, 1, validator.adata.obs.shape[0])
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: obs['{obs_col}'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
        ]

    def test_non_visium_all_obs_cols(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["in_tissue"] = 1
        validator.adata.obs["array_col"] = 1
        validator.adata.obs["array_row"] = 1
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: obs['array_col'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
            "ERROR: obs['array_row'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
            "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
        ]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestObsVisiumCols:
    def test_array_col_passes(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    # fails in lots of ways beyond visium-specific obs cols
    def test_change_assay_to_non_spatial(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0011025")
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0011025"
        validator.validate_adata()
        assert validator.is_valid is False

    @pytest.mark.parametrize(
        "obs_value", (None, -1, -1.0, 2, 2.0, "string", True, False)
    )
    def test_in_tissue_values(self, validator_with_adatas, obs_value):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], "in_tissue"] = obs_value
        validator.validate_adata()
        assert validator.is_valid is False

    @pytest.mark.parametrize(
        "obs_col", ("array_row", "array_col", "in_tissue")
    )
    def test_array_col_names(self, validator_with_adatas, obs_col):
        validator = validator_with_adatas
        del validator.adata.obs[obs_col]
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: obs['{obs_col}'] is required for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
        ]

    # in-range floats converted to int and pass
    @pytest.mark.parametrize(
        "obs_col", ("array_row", "array_col", "in_tissue")
    )
    def test_obs_float_values_pass(self, validator_with_adatas, obs_col):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], obs_col] = 1.0
        validator.validate_adata()
        assert validator.is_valid is True
        assert validator.errors == []

    @pytest.mark.parametrize(
        "obs_col", ("array_row", "array_col", "in_tissue")
    )
    def test_array_col_values_is_single_false(self, validator_with_adatas, obs_col):
        validator = validator_with_adatas
        del validator.adata.obs[obs_col]
        validator.adata.uns["spatial"]["is_single"] = False
        validator.validate_adata()
        assert validator.is_valid is False
        assert len(validator.errors) == 4

    @pytest.mark.parametrize("obs_col", ["array_row", "array_col", "in_tissue"])
    @pytest.mark.parametrize("obs_type", ["float32", "float64", "str", "object", "bool"])
    def test_non_int_col_row_in_tissue(self, validator_with_adatas, obs_col, obs_type):
        validator = validator_with_adatas
        validator.adata.obs[obs_col] = validator.adata.obs[obs_col].astype(obs_type)   
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            f"ERROR: obs['{obs_col}'] must be of int type, it is {validator.adata.obs[obs_col].dtype}."
        ]


class TestArrayColValues:
    @pytest.mark.parametrize("test_h5ads", ["visium_human_some_spots.h5ad", "visium_human_all_spots.h5ad"])
    @pytest.mark.parametrize(
        "obs_value", (None, -1, -1.0, 128, 128.0, "string", True, False)
    )
    def test_array_col_values_v1(self, validator_with_adatas, obs_value):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], "array_col"] = obs_value
        validator.validate_adata()
        assert validator.is_valid is False

    @pytest.mark.parametrize("test_h5ads", ["visium_human_some_spots.h5ad", "visium_human_all_spots.h5ad"])
    @pytest.mark.parametrize(
        "obs_value", (None, -1, -1.0, 78, 78.0, "string", True, False)
    )
    def test_array_row_values_v1(self, validator_with_adatas, obs_value):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
        validator.validate_adata()
        assert validator.is_valid is False

    @pytest.mark.parametrize("test_h5ads", ["visium_v2_11mm_human.h5ad"])
    @pytest.mark.parametrize(
        "obs_value", (None, -1, -1.0, 224, 224.0, "string", True, False)
    )
    def test_array_col_values_v2(self, validator_with_adatas, obs_value):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], "array_col"] = obs_value
        validator.validate_adata()
        assert validator.is_valid is False

    @pytest.mark.parametrize("test_h5ads", ["visium_v2_11mm_human.h5ad"])
    @pytest.mark.parametrize(
        "obs_value", (None, -1, -1.0, 128, 128.0, "string", True, False)
    )
    def test_array_row_values_v2(self, validator_with_adatas, obs_value):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs.index[0], "array_row"] = obs_value
        validator.validate_adata()
        assert validator.is_valid is False
