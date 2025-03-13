"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/873
https://github.com/chanzuckerberg/single-cell-curation/pull/885
"""

import numpy as np
import pandas as pd
import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
    NON_SPATIAL_H5ADS,
    SPATIAL_H5ADS,
)


VISIUM_H5ADS = [file for file in SPATIAL_H5ADS if "visium" in file]
SLIDE_SEQ_H5ADS = [file for file in SPATIAL_H5ADS if "slide_seq" in file]
ASSAY_ERROR_MESSAGE = (
    "ERROR: When obs['assay_ontology_term_id'] is either a descendant of 'EFO:0010961' "
    "(Visium Spatial Gene Expression) or 'EFO:0030062' (Slide-seqV2), all observations "
    "must contain the same value."
)


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_mixed_assays_spatial(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert ASSAY_ERROR_MESSAGE in validator.errors


@pytest.mark.parametrize("test_h5ads", SLIDE_SEQ_H5ADS)
def test_mixed_assays_slide_seq_is_single_false(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert ASSAY_ERROR_MESSAGE in validator.errors


@pytest.mark.parametrize("test_h5ads", VISIUM_H5ADS)
def test_mixed_assays_visium_is_single_false(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0009922")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0009922"
    validator.adata.uns["spatial"]["is_single"] = False
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() == 2
    assert validator.is_valid is False
    assert ASSAY_ERROR_MESSAGE in validator.errors


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
def test_mixed_assays_non_spatial(validator_with_adatas):
    validator = validator_with_adatas
    validator.adata.obs["assay_ontology_term_id"] = validator.adata.obs["assay_ontology_term_id"].cat.add_categories("EFO:0022857")
    validator.adata.obs.loc[validator.adata.obs.index[1], "assay_ontology_term_id"] = "EFO:0022857"
    validator.validate_adata()
    assert validator.adata.obs["assay_ontology_term_id"].nunique() >= 2
    assert validator.is_valid is False
    assert ASSAY_ERROR_MESSAGE in validator.errors
