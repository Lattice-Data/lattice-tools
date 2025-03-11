import anndata as ad
import gc
import os
import pytest
import sys
from pathlib import Path
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad

# pytest can now discover and successfully run tests from any directory level of repo
FIXTURES_ROOT = Path(__file__).absolute().parent

# few standard lists to select fixtures, can also whitelist for any given test
ALL_H5ADS = [
    f.name for f in FIXTURES_ROOT.iterdir() 
        if f.suffix == ".h5ad"
]
SPATIAL_H5ADS = [
    f for f in ALL_H5ADS 
        if any(included in f for included in [
            "slide_seq",
            "visium"
        ])
]
MULTISPECIES_H5ADS = [
    f for f in ALL_H5ADS 
        if not any(excluded in f for excluded in [
            "human", 
            "mouse"
        ])
]
MODEL_ORGANISM_H5ADS = [
    f for f in MULTISPECIES_H5ADS 
        if any(included in f for included in [
            "fly", 
            "worm", 
            "zebrafish"
        ])
]

print(ALL_H5ADS)
print(SPATIAL_H5ADS)
print(MULTISPECIES_H5ADS)
print(MODEL_ORGANISM_H5ADS)
# will add better check for file, maybe to automatically download as well
if not os.path.isfile(FIXTURES_ROOT / "visium_v2_11mm_human.h5ad"):
    raise FileNotFoundError('This file lives in S3, please download before running tests')


# base fixture to get any given h5ad
@pytest.fixture(params=ALL_H5ADS)
def test_h5ads(request):
    yield request.param


# fixture exported to other tests, returns and therefor tests with each h5ad
@pytest.fixture
def validator_with_adatas(test_h5ads) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(FIXTURES_ROOT / test_h5ads)
    yield validator


@pytest.fixture()
def validator_human_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{ALL_H5ADS[0]}")
    yield validator


@pytest.fixture()
def validator_mouse_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{ALL_H5ADS[1]}")
    yield validator


@pytest.fixture(params=[file for file in ALL_H5ADS if "human" in file])
def validator_with_human_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=ALL_H5ADS[2:])
def validator_with_spatial_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=ALL_H5ADS[2:4])
def validator_with_slide_seq_adatas(request) -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    yield validator


@pytest.fixture(params=ALL_H5ADS[:3])
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


@pytest.fixture(params=ALL_H5ADS[4:])
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


@pytest.fixture(params=ALL_H5ADS[:2])
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


@pytest.fixture
def validator_with_rabbit_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_rabbit.h5ad")
    yield validator

@pytest.fixture
def validator_with_gorilla_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_gorilla.h5ad")
    yield validator

@pytest.fixture
def validator_with_marmoset_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_marmoset.h5ad")
    yield validator

@pytest.fixture
def validator_with_pig_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_domestic_pig.h5ad")
    yield validator

@pytest.fixture
def validator_with_chimp_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_chimp.h5ad")
    yield validator

@pytest.fixture
def validator_with_rhesus_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_rhesus.h5ad")
    yield validator


@pytest.fixture
def validator_with_crab_eating_macaque_adata() -> Validator:
    gc.collect()
    validator = Validator()
    validator.adata = read_h5ad(f"{FIXTURES_ROOT}/valid_crab_eating_macaque.h5ad")
    yield validator
