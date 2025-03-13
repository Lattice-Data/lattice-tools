"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/830
https://github.com/chanzuckerberg/single-cell-curation/pull/864
"""

import pytest
import numpy as np
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    ALL_H5ADS,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
def test_slide_seq_in_tissue(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_non_spatial_in_tissue(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["in_tissue"] = 1
    validator.validate_adata()
    assert validator.is_valid is False
    assert validator.errors == [
        "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
        "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
    ]


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
class TestInTissueCLVisium:
    def test_change_unknown_to_cell_type(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.replace({"cell_type_ontology_term_id": {"unknown": "CL:0002601"}}, inplace=True)
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] is a descendant of "
            "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
        ]

    def test_change_all_in_tissue_to_zero(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["in_tissue"] = 0
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] is a descendant of "
            "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
        ]

    def test_one_random_in_tissue_zero(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.replace({"cell_type_ontology_term_id": {"unknown": "CL:0002601"}}, inplace=True)
        validator.adata.obs["in_tissue"] = 1
        random_index = np.random.randint(0, (validator.adata.obs.shape[0] - 1))
        validator.adata.obs.loc[validator.adata.obs.index[random_index], "in_tissue"] = 0
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] is a descendant of "
            "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
        ]

    def test_change_in_tissue_to_zero(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.replace({"in_tissue": {1: 0}}, inplace=True)
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors == [
            "ERROR: obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] is a descendant of "
            "'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."
        ]

    @pytest.mark.parametrize(
        "change_dict", (
            {"cell_type_ontology_term_id": {"unknown": "CL:0002601"}},
            {"in_tissue": {1: 0}}
        )
    )
    def test_is_single_false(self, validator_with_adatas, change_dict):
        validator = validator_with_adatas
        validator.adata.uns["spatial"]["is_single"] = False
        validator.adata.obs.replace(change_dict, inplace=True)
        validator.validate_adata()
        assert validator.is_valid is False
        assert validator.errors[:4] == [
            "ERROR: uns['spatial'][library_id] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
            "ERROR: obs['array_col'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
            "ERROR: obs['array_row'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True.",
            "ERROR: obs['in_tissue'] is only allowed for obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' "
            "(Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True."
        ]

    # should pass, rare but sometimes cell type is unknown
    def test_in_tissue_one_unknown(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs.loc[validator.adata.obs["cell_type_ontology_term_id"] == "unknown", "in_tissue"] = 1
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []
