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
NON_SPATIAL_H5ADS = [
    f for f in ALL_H5ADS 
        if not any(excluded in f for excluded in [
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

def get_library_id(adata):
    return [key for key in adata.uns['spatial'].keys() if 'is_single' not in key][0]


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


@pytest.fixture
def label_writer(validator_with_adatas: Validator) -> AnnDataLabelAppender:
    gc.collect()
    validator = validator_with_adatas
    validator.validate_adata()
    yield AnnDataLabelAppender(validator.adata)
