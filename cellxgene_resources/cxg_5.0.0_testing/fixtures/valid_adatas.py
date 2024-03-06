import anndata as ad
import os
import pytest
import sys


# pytest can now discover and successfully run tests from any directory level of repo
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

SCC_REPO_LOC = os.path.expanduser('~/GitClones/CZI/')
sys.path.append(os.path.abspath(SCC_REPO_LOC + 'single-cell-curation/cellxgene_schemea_cli/'))

from cellxgene_schema.validate import Validator


# fixture exported to other tests, returns and therefor tests with each h5ad
@pytest.fixture(params=['valid_human.h5ad', 'valid_mouse.h5ad'])
def validator_with_adata(request, fixtures_root=FIXTURES_ROOT) -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad(f'{fixtures_root}/{request.param}')
    return validator
