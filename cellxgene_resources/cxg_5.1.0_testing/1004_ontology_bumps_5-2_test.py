"""
QA testing for this issue:
https://github.com/chanzuckerberg/cellxgene-ontology-guide/pull/217
https://github.com/chanzuckerberg/single-cell-curation/pull/1004/
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
    "col, term",
    (
        pytest.param("cell_type_ontology_term_id", "CL:4033072", id="cycling gamma-delta T cell"),
        pytest.param("disease_ontology_term_id", "MONDO:1012443", id="malignant melanoma, lion"),
        pytest.param("disease_ontology_term_id", "MONDO:1040020", id="methicillin-susceptible staphylococcus aureus infectious disease"),
        pytest.param("disease_ontology_term_id", "MONDO:1012474", id="retinal degeneration, Smoky Joe, chicken "),
        pytest.param("assay_ontology_term_id", "EFO:0022839", id="STORM-seq"),
        pytest.param("tissue_ontology_term_id", "UBERON:0000201", id="endothelial blood brain barrier"),
        pytest.param("self_reported_ethnicity_ontology_term_id", "HANCESTRO:0597", id="Singaporean Malay"),
    )
)
def test_added_mappings(validator_with_non_spatial_adata, col, term):
    validator = validator_with_non_spatial_adata
    validator.adata.obs[col] = term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
