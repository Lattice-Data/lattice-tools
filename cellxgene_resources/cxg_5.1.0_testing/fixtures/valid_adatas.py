import anndata as ad
import os
import pytest
import sys


# pytest can now discover and successfully run tests from any directory level of repo
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

SCC_REPO_LOC = os.path.expanduser("~/GitClones/CZI/")
sys.path.append(
    os.path.abspath(SCC_REPO_LOC + "single-cell-curation/cellxgene_schemea_cli/")
)

from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender


H5ADS = [
    "valid_human.h5ad", 
    "visium_human_all_spots.h5ad",
    "slide_seq_image_human.h5ad",
    "slide_seq_no_image_human.h5ad",
]

# fixture exported to other tests, returns and therefor tests with each h5ad
@pytest.fixture(params=H5ADS)
def validator_with_all_adatas(request) -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    return validator


@pytest.fixture(params=H5ADS[1:])
def validator_with_spatial_adatas(request) -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    return validator


@pytest.fixture
def validator_with_visium() -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad(f"{FIXTURES_ROOT}/visium_human_all_spots.h5ad")
    return validator


@pytest.fixture
def label_writer(validator_with_all_adatas: Validator) -> AnnDataLabelAppender:
    validator = validator_with_all_adatas
    validator.validate_adata()
    return AnnDataLabelAppender(validator)


@pytest.fixture
def validator_with_non_spatial_adata() -> Validator:
    validator = Validator()
    validator.adata = ad.read_h5ad(f"{FIXTURES_ROOT}/valid_human.h5ad")
    return validator
