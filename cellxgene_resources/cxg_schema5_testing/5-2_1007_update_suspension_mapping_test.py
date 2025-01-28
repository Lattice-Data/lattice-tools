"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/953
https://github.com/chanzuckerberg/single-cell-curation/pull/1007/
"""

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


@pytest.mark.parametrize(
    "assay, suspension",
    (
        pytest.param("EFO:0030060", "cell", id="mCT-seq with cell"),
        pytest.param("EFO:0030060", "nucleus", id="mCT-seq with nucleus"),
        pytest.param("EFO:0008992", "na", id="MERFISH with na"),
        pytest.param("EFO:0002761", "nucleus", id="methylation profiling by high throughput sequencing with nucleus"),
        pytest.param("EFO:0022490", "cell", id="ScaleBio single cell RNA sequencing with cell"),
        pytest.param("EFO:0022490", "nucleus", id="ScaleBio single cell RNA sequencing with nucleus"),
        pytest.param("EFO:0030028", "cell", id="sci-RNA-seq3 with cell"),
        pytest.param("EFO:0030028", "nucleus", id="sci-RNA-seq3 with nucleus"),
    )
)
def test_added_mappings(validator_with_non_spatial_adata, assay, suspension):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["assay_ontology_term_id"] = assay
    validator.adata.obs["suspension_type"] = suspension
    for col in ["assay_ontology_term_id", "suspension_type"]:
        validator.adata.obs[col] = validator.adata.obs[col].astype("category")
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize(
    "assay, suspension",
    (
        pytest.param("EFO:0030060", "na", id="mCT-seq with na"),
        pytest.param("EFO:0008992", "cell", id="MERFISH with cell"),
        pytest.param("EFO:0002761", "cell", id="methylation profiling by high throughput sequencing with cell"),
        pytest.param("EFO:0022490", "na", id="ScaleBio single cell RNA sequencing with na"),
        pytest.param("EFO:0030028", "na", id="sci-RNA-seq3 with na"),
    )
)
def test_mismatched_mappings(validator_with_non_spatial_adata, assay, suspension):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["assay_ontology_term_id"] = assay
    validator.adata.obs["suspension_type"] = suspension
    for col in ["assay_ontology_term_id", "suspension_type"]:
        validator.adata.obs[col] = validator.adata.obs[col].astype("category")
    validator.validate_adata()
    assert not validator.is_valid
    assert len(validator.errors) == 1


@pytest.mark.parametrize(
    "assay, suspension",
    (
        pytest.param("EFO:0009294", "cell", id="CITE-seq with cell"),
        pytest.param("EFO:0009294", "nucleus", id="CITE-seq with nucleus"),
        pytest.param("EFO:0009294", "na", id="CITE-seq with na"),
        pytest.param("EFO:0009918", "na", id="smFISH  with na"),
        pytest.param("EFO:0009918", "cell", id="smFISH  with cell"),
        pytest.param("EFO:0009918", "nucleus", id="smFISH  with nucleus"),
        pytest.param("EFO:0008939", "nucleus", id="snmC-seq with nucleus"),
        pytest.param("EFO:0008939", "cell", id="snmC-seq with cell"),
        pytest.param("EFO:0008939", "na", id="snmC-seq with na"),
        pytest.param("EFO:0030027", "nucleus", id="snmC-seq2 with nucleus"),
        pytest.param("EFO:0700000", "na", id="spatial proteomics and its descendants with na"),
        pytest.param("EFO:0700000", "cell", id="spatial proteomics and its descendants with cell"),
        pytest.param("EFO:0700000", "nucleus", id="spatial proteomics and its descendants with nucleus"),
    )
)
def test_deleted_mappings(validator_with_non_spatial_adata, assay, suspension):
    validator = validator_with_non_spatial_adata
    validator.adata.obs["assay_ontology_term_id"] = assay
    validator.adata.obs["suspension_type"] = suspension
    for col in ["assay_ontology_term_id", "suspension_type"]:
        validator.adata.obs[col] = validator.adata.obs[col].astype("category")
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
