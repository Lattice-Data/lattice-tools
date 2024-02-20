"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/719
"""

import os
import sys
import anndata as ad
import pytest

scc_repo_loc = os.path.expanduser('~/GitClones/CZI/')

sys.path.append(os.path.abspath(scc_repo_loc + 'single-cell-curation/cellxgene_schemea_cli/'))

from cellxgene_schema.validate import Validator


@pytest.fixture
def validator_with_adata() -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad('fixtures/valid.h5ad')
    return validator


@pytest.mark.parametrize(
    "dt_id,expected",
    (
        pytest.param('PATO:0000461', True, id='Valid PATO:0000461, healthy'),
        pytest.param('MONDO:0700096', True, id='Valid MONDO:0700096, human disease parent term'),
        pytest.param('MONDO:0005491', True, id='Valid MONDO:0005491, human disease child term'),
        pytest.param('MONDO:0021178', True, id='Valid MONDO:0021178, injury parent term'),
        pytest.param('MONDO:0015796', True, id='Valid MONDO:0015796, injury child term'),
        pytest.param('MONDO:0005583', True, id='Valid MONDO:0005583, non-human disease parent term'),
        pytest.param('MONDO:1011335', True, id='Valid MONDO:1011335, non-human disease child term'),
        pytest.param('MONDO:1011336', True, id='Valid MONDO:1011336, nervous system disorder, non-human animal'),
        pytest.param('MONDO:1010239', True, id='Valid MONDO:1010239, peripheral neuropathy, non-human animal'),
        pytest.param('MONDO:1010003', True, id='Valid MONDO:1010003, narcolepsy non-human animal'),
        pytest.param('MONDO:1010421', True, id='Valid MONDO:1010421, animal disease, sheep narcolepsy'),
        pytest.param('MONDO:0800478', True, id='Valid MONDO:0800478, new term with pinned release, trigeminal trophic syndrome'),
    )
)
def test_disease_ontology_term_id_true(validator_with_adata, dt_id, expected):
    validator = validator_with_adata
    validator.adata.obs['disease_ontology_term_id'] = dt_id
    validator.validate_adata()
    assert validator.is_valid is expected


@pytest.mark.parametrize(
    "dt_id,expected",
    (
        pytest.param('MONDO:0042489', False, id='Invalid MONDO:0042489, disease susceptiblity parent term'),
        pytest.param('PATO:0001894', False, id='Invalid PATO:0001894 term id'),
        pytest.param('MONDO:0021125', False, id='Invalid MONDO:0021125, disease characteristic'),
        pytest.param('MONDO:0000001', False, id='Invalid MONDO:0000001, disease parent term'),
        pytest.param('MONDO:0012153', False, id='Invalid MONDO:0012153, disease susceptiblity child term'),
        pytest.param('MONDO:0021135', False, id='Invalid MONDO:0021135, disease characteristic child'),
    )
)
def test_disease_ontology_term_false(validator_with_adata, dt_id, expected):
    validator = validator_with_adata
    validator.adata.obs['disease_ontology_term_id'] = dt_id
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{dt_id}' in 'disease_ontology_term_id' is not an allowed term id. "
        "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or children terms thereof, "
        "or children terms of 'MONDO:0000001' (disease) are allowed"
    ]


@pytest.mark.parametrize(
    "dt_id,expected",
    (
        pytest.param('', False, id='Empty str'),
        pytest.param('EFO:000001', False, id='EFO:000001, Invalid ontology'),
        pytest.param('NCIT:C158547', False, id='NCIT:C158547, Invalid ontology'),
        pytest.param('MONDO:0100535', False, id='Invalid MONDO:0100535, not approved as of 1-25-24'),
        pytest.param('MONDO:0100536', False, id='Invalid child of MONDO:0100535, not approved as of 1-25-24'),
    )
)
def test_disease_ontology_term_not_valid(validator_with_adata, dt_id, expected):
    validator = validator_with_adata
    validator.adata.obs['disease_ontology_term_id'] = dt_id
    validator.validate_adata()
    assert validator.is_valid is expected
    assert validator.errors == [
        f"ERROR: '{dt_id}' in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'. "
        "Only 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or children terms thereof, "
        "or children terms of 'MONDO:0000001' (disease) are allowed"
    ]
