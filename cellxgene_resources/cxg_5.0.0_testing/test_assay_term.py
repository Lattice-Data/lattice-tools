"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/717
"""

import pytest
from fixtures.valid_adatas import validator_with_adata


@pytest.mark.parametrize(
    "at_id,expected",
    (
        pytest.param('EFO:0001456', True, id='Valid EFO:0001456, descendant of assay by molecule'),
        pytest.param('EFO:0008995', True, id='Valid EFO:0008995, descendant of single cell library construction'),
        pytest.param('EFO:0009900', True, id='Valid EFO:0009900, current 4.0.0 pinned ontology'),
        pytest.param('EFO:0030003', True, id='Valid EFO:0030003, 10x 3prime transcription profiling parent term'),
        pytest.param('EFO:0022490', True, id='New valid EFO:0022490, new term in pinned ontology bump'),
        pytest.param('EFO:0022492', True, id='New valid EFO:0022492, new term in pinned ontology bump'),
    )
)
def test_assay_ontology_term_id_valid(validator_with_adata, at_id, expected):
    validator = validator_with_adata
    validator.adata.obs['assay_ontology_term_id'] = at_id
    validator.validate_adata()
    assert validator.is_valid is expected


@pytest.mark.parametrize(
    "at_id,expected",
    (
        pytest.param('EFO:0000001', False, id='Not allowed EFO:0000001, top EFO term'),
        pytest.param('EFO:0002772', False, id='Not allowed EFO:0002772, assay by molecule'),
        pytest.param('EFO:0010183', False, id='Not allowed EFO:0010183, single cell library construction'),
    )
)
def test_assay_ontology_term_id_not_allowed(validator_with_adata, at_id, expected):
    validator = validator_with_adata
    validator.adata.obs['assay_ontology_term_id'] = at_id
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{at_id}' in 'assay_ontology_term_id' is not an allowed term id. "
        "Only descendant terms of either 'EFO:0002772' or 'EFO:0010183' are allowed for assay_ontology_term_id"
    ]


@pytest.mark.parametrize(
    "at_id,expected",
    (
        pytest.param('', False, id='Empty string'),
        pytest.param('PATO:0000461', False, id='PATO:0000461, wrong ontology'),
    )
)
def test_assay_ontology_term_id_not_valid(validator_with_adata, at_id, expected):
    validator = validator_with_adata
    validator.adata.obs['assay_ontology_term_id'] = at_id
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{at_id}' in 'assay_ontology_term_id' is not a valid ontology term id of 'EFO'. "
        "Only descendant terms of either 'EFO:0002772' or 'EFO:0010183' are allowed for assay_ontology_term_id"
    ]


@pytest.mark.parametrize(
    "at_id,expected",
    (
        pytest.param('EFO:0700005', True, id='Valid EFO:0700005, warning for suspension schema'),
    )
)
def test_assay_ontology_term_id_suspension_warning(validator_with_adata, at_id, expected):
    validator = validator_with_adata
    validator.adata.obs['assay_ontology_term_id'] = at_id
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == []
    assert validator.warnings == [
        "WARNING: Data contains assay(s) that are not represented in the 'suspension_type' schema definition table. "
        "Ensure you have selected the most appropriate value for the assay(s) between 'cell', 'nucleus', and 'na'. "
        "Please contact cellxgene@chanzuckerberg.com during submission so that the assay(s) can be added to the schema "
        "definition document."
    ]
