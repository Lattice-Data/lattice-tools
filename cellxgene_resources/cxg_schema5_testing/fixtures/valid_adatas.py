import anndata as ad
import gc
import os
import pytest
import sys
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad

# pytest can now discover and successfully run tests from any directory level of repo
FIXTURES_ROOT = os.path.join(os.path.dirname(__file__))

H5ADS = [
    "valid_human.h5ad", 
    "valid_mouse.h5ad", 
    "slide_seq_image_human.h5ad",
    "slide_seq_no_image_human.h5ad",
    "visium_v2_11mm_human.h5ad",    # originally CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_curated.h5ad from Jenny
    "visium_human_all_spots.h5ad",
    "visium_human_some_spots.h5ad",
]

MULTISPECIES_H5ADS = [
    "valid_fly.h5ad",
    "valid_worm.h5ad",
    "valid_zebrafish.h5ad",
    "valid_lemur.h5ad",
    "valid_rat.h5ad"
]

# will add better check for file, maybe to automatically download as well
if not os.path.isfile(f"{FIXTURES_ROOT}/visium_v2_11mm_human.h5ad"):
    raise FileNotFoundError('This file lives in S3, please download before running tests')

# fixture exported to other tests, returns and therefor tests with each h5ad
@pytest.fixture(params=H5ADS)
def validator_with_all_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture()
def validator_human_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{H5ADS[0]}")
    yield validator


@pytest.fixture()
def validator_mouse_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{H5ADS[1]}")
    yield validator


@pytest.fixture(params=[file for file in H5ADS if "human" in file])
def validator_with_human_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=H5ADS[2:])
def validator_with_spatial_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=H5ADS[2:4])
def validator_with_slide_seq_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=H5ADS[:3])
def validator_with_non_visium_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture
def validator_with_visium_some() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/visium_human_some_spots.h5ad")
    yield validator


@pytest.fixture
def validator_with_visium() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/visium_human_all_spots.h5ad")
    yield validator


@pytest.fixture(params=H5ADS[4:])
def validator_with_all_visiums(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture
def label_writer(validator_with_all_adatas: Validator) -> AnnDataLabelAppender:
    gc.collect()
    validator = validator_with_all_adatas
    validator.validate_adata()
    yield AnnDataLabelAppender(validator.adata)


@pytest.fixture
def label_writer_multispecies(validator_with_multispecies_adatas: Validator) -> AnnDataLabelAppender:
    gc.collect()
    validator = validator_with_multispecies_adatas
    validator.validate_adata()
    yield AnnDataLabelAppender(validator.adata)


@pytest.fixture(params=H5ADS[:2])
def validator_with_non_spatial_adata(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=MULTISPECIES_H5ADS)
def validator_with_multispecies_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture
def validator_with_fly_adata(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_fly.h5ad")
    yield validator


@pytest.fixture
def validator_with_worm_adata(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_worm.h5ad")
    yield validator


@pytest.fixture
def validator_with_zebrafish_adata(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_zebrafish.h5ad")
    yield validator


@pytest.fixture
def validator_with_lemur_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_lemur.h5ad")
    yield validator

@pytest.fixture
def validator_with_rat_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_rat.h5ad")
    yield validator