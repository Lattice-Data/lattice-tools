"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/901
https://github.com/chanzuckerberg/single-cell-curation/pull/902
"""

import numpy as np
import pandas as pd
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


# couldn't find an applicable EFO NTR that we might use, didn't seem to be any new assay terms
@pytest.mark.parametrize(
    "col, term",
    (
        pytest.param("cell_type_ontology_term_id", "CL:4042002", id="TAC3-positive medium spiny neuron"),
        pytest.param("disease_ontology_term_id", "MONDO:0958349", id="dorsal spinal cord lipoma"),
        pytest.param("tissue_ontology_term_id", "UBERON:8600048", id="apex of urinary bladder"),
        # clonal hematopoiesis terms applicable to one published collection
        # https://github.com/chanzuckerberg/single-cell-curation/issues/825
        pytest.param("disease_ontology_term_id", "MONDO:0100542", id="clonal hematopoiesis"),
        pytest.param("disease_ontology_term_id", "MONDO:0100544", id="age-related clonal hematopoiesis"),
        pytest.param("disease_ontology_term_id", "MONDO:0100543", id="clonal hematopoiesis of indeterminate potential"),
    )
)
def test_new_ontology_terms(validator_with_non_spatial_adata, col, term):
    validator = validator_with_non_spatial_adata
    validator.adata.obs[col] = term
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
