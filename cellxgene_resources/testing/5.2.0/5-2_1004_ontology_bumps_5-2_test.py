"""
QA testing for this issue:
https://github.com/chanzuckerberg/cellxgene-ontology-guide/pull/217
https://github.com/chanzuckerberg/single-cell-curation/pull/1004/
"""

import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
)

H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
]

class TestAddedMappings:
    @pytest.mark.parametrize("test_h5ads", H5ADS)
    @pytest.mark.parametrize(
        "col, term",
        (
            pytest.param("cell_type_ontology_term_id", "CL:4033072", id="cycling gamma-delta T cell"),
            pytest.param("disease_ontology_term_id", "MONDO:1012443", id="malignant melanoma, lion"),
            pytest.param("disease_ontology_term_id", "MONDO:1040020", id="methicillin-susceptible staphylococcus aureus infectious disease"),
            pytest.param("disease_ontology_term_id", "MONDO:1012474", id="retinal degeneration, Smoky Joe, chicken "),
            pytest.param("assay_ontology_term_id", "EFO:0022839", id="STORM-seq"),
            pytest.param("tissue_ontology_term_id", "UBERON:0000201", id="endothelial blood brain barrier"),
        )
    )
    def test_added_mappings(self, validator_with_adatas, col, term):
        validator = validator_with_adatas
        validator.adata.obs[col] = term
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    @pytest.mark.parametrize("test_h5ads", ["valid_human.h5ad"])
    def test_added_ethnicity(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = "HANCESTRO:0597"
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []
