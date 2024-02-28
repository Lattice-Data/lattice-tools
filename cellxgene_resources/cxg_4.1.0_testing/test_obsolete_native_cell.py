"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/735
"""

import pytest
from fixtures.valid_adatas import validator_with_adata


@pytest.mark.parametrize(
    "cell_type_term,tissue_type,expected",
    (
        pytest.param('unknown', 'tissue', True, id='cell_type_term is unknown'),
        pytest.param('CL:0000034', 'tissue', True, id='CL:0000034, child of CL:0000548, animal cell, tissue_type == tissue'),
        pytest.param('CL:0000463', 'tissue', True, id='CL:0000463, child of CL:0000255, eukaryotic cell, tissue_type == tissue'),
        pytest.param('unknown', 'cell culture', True, id='cell_type_term is unknown, tissue_type == cell culture'),
        pytest.param('CL:0000034', 'cell culture', True, id='CL:0000034, child of CL:0000548, animal cell, tissue_type == cell culture'),
        pytest.param('CL:0000463', 'cell culture', True, id='CL:0000463, child of CL:0000255, eukaryotic cell, tissue_type == cell culture'),
        pytest.param('unknown', 'organoid', True, id='cell_type_term is unknown, tissue_type == organoid'),
        pytest.param('CL:0000034', 'organoid', True, id='CL:0000034, child of CL:0000548, animal cell, tissue_type == organoid'),
        pytest.param('CL:0000463', 'organoid', True, id='CL:0000463, child of CL:0000255, eukaryotic cell, tissue_type == organoid'),
    )
)
def test_cell_ontology_term_is_valid(validator_with_adata, cell_type_term, tissue_type, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')

    # set to correct col here, will set to wrong col in other tests
    if tissue_type == 'cell culture':
        validator.adata.obs['tissue_ontology_term_id'] = cell_type_term
    else:
        validator.adata.obs['cell_type_ontology_term_id'] = cell_type_term

    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == []


@pytest.mark.parametrize(
    "cell_type_term,tissue_type,expected",
    (
        pytest.param('CL:0000003', 'tissue', False, id='CL:0000003, native cell, tissue_type == tissue'),
        pytest.param('CL:0000548', 'tissue', False, id='CL:0000548, animal cell, tissue_type == tissue'),
        pytest.param('CL:0000003', 'organoid', False, id='CL:0000003, native cell, tissue_type == organoid'),
        pytest.param('CL:0000548', 'organoid', False, id='CL:0000548, animal cell, tissue_type == organoid'),
    )
)
def test_cell_ontology_term_fails_deprecated(validator_with_adata, cell_type_term, tissue_type, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['cell_type_ontology_term_id'] = cell_type_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{cell_type_term}' in 'cell_type_ontology_term_id' is a deprecated term id of 'CL'."
    ]


@pytest.mark.parametrize(
    "cell_type_term,tissue_type,expected",
    (
        pytest.param('CL:0000255', 'tissue', False, id='CL:0000255, eukaryotic cell, tissue_type == tissue'),
        pytest.param('CL:0000255', 'organoid', False, id='CL:0000255, eukaryotic cell, tissue_type == organoid'),
    )
)
def test_cell_ontology_term_fails_not_allowed(validator_with_adata, cell_type_term, tissue_type, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['cell_type_ontology_term_id'] = cell_type_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{cell_type_term}' in 'cell_type_ontology_term_id' is not allowed."
    ]


@pytest.mark.parametrize(
    "cell_type_term,tissue_type,expected",
    (
        pytest.param('CL:0000255', 'cell culture', False, id='CL:0000255, eukaryotic cell, tissue_type == cell culture'),
    )
)
def test_cell_ontology_term_cell_culture_fails(validator_with_adata, cell_type_term, tissue_type, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['tissue_ontology_term_id'] = cell_type_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{cell_type_term}' in 'tissue_ontology_term_id' is not allowed. When 'tissue_type' is 'cell culture', "
        "'tissue_ontology_term_id' MUST be either a CL term (excluding 'CL:0000255' (eukaryotic cell), 'CL:0000257' "
        "(Eumycetozoan cell), and 'CL:0000548' (animal cell)) or 'unknown'."
    ]


@pytest.mark.parametrize(
    "cell_type_term,tissue_type,expected",
    (
        pytest.param('CL:0000003', 'cell culture', False, id='CL:0000003, native cell, tissue_type == cell culture'),
        pytest.param('CL:0000548', 'cell culture', False, id='CL:0000548, animal cell, tissue_type == cell culture'),
    )
)
def test_cell_ontology_term_cell_culture_fails_deprecated(validator_with_adata, cell_type_term, tissue_type, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['tissue_ontology_term_id'] = cell_type_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{cell_type_term}' in 'tissue_ontology_term_id' is a deprecated term id of 'CL'. "
        "When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST be either a CL term (excluding 'CL:0000255'"
        " (eukaryotic cell), 'CL:0000257' (Eumycetozoan cell), and 'CL:0000548' (animal cell)) or 'unknown'."
    ]


@pytest.mark.parametrize(
    "tissue_type,ontology_term,expected",
    (
        pytest.param('cell culture', 'UBERON:0001062', False, id='cell culture with UBERON term'),
        pytest.param('cell culture', '', False, id='Empty string, cell culture'),
        pytest.param('cell culture', 'EFO:0000001', False, id='cell culture with EFO term'),
        pytest.param('cell culture', 'UBERON:0001062', False, id='UBERON:0001062, animal cell, tissue_type == cell culture'),
    )
)
def test_cell_culture_wrong_mismatched_tissue_ontology_fails(validator_with_adata, tissue_type, ontology_term, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['tissue_ontology_term_id'] = ontology_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{ontology_term}' in 'tissue_ontology_term_id' is not a valid ontology term id of 'CL'. "
        "When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST be either a CL term (excluding 'CL:0000255' (eukaryotic cell), "
        "'CL:0000257' (Eumycetozoan cell), and 'CL:0000548' (animal cell)) or 'unknown'."
    ]


@pytest.mark.parametrize(
    "tissue_type,ontology_term,expected",
    (
        pytest.param('tissue', '', False, id='Empty string, tissue'),
        pytest.param('tissue', 'CL:0000034', False, id='tissue with CL term'),
        pytest.param('tissue', 'EFO:0000001', False, id='tissue with EFO term'),
        pytest.param('organoid', '', False, id='Empty string, organoid'),
        pytest.param('organoid', 'EFO:0000001', False, id='organoid with EFO term'),
        pytest.param('organoid', 'CL:0000034', False, id='organoid with CL term'),
    )
)
def test_tissue_organoid_wrong_mismatched_tissue_ontology_fails(validator_with_adata, tissue_type, ontology_term, expected):
    validator = validator_with_adata
    validator.adata.obs['tissue_type'] = tissue_type
    validator.adata.obs['tissue_type'] = validator.adata.obs['tissue_type'].astype('category')
    validator.adata.obs['tissue_ontology_term_id'] = ontology_term
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{ontology_term}' in 'tissue_ontology_term_id' is not a valid ontology term id of 'UBERON'. "
        "When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' MUST be a child term id of 'UBERON:0001062' (anatomical entity)."
    ]
