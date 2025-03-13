"""
Fixture structure framework:

For any given test file, need to at least import the following:
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
)

By default, validator_with_adatas will yield a Validator instance with attached adata
for every h5ad in ALL_H5ADS.

def test_all(validator_with_adatas):
    validator = validator_with_adatas 
    ...
    assert ...

This tests everything.

To limit which fixtures to use, provide a decorator for "test_h5ads" with an iterable:

@pytest.mark.parameterize("test_h5ads", SPATIAL_H5ADS)
def test_spatial(validator_with_adatas):
    validator = validator_with_adatas 
    ...
    assert ...

Now only the h5ads in SPATIAL_H5ADS will be tested

Tests can be grouped into classes to use the same h5ad subset for multiple tests without
having to repeat the decoration for each test method/function:

@pytest.mark.parameterize("test_h5ads", SPATIAL_H5ADS)
class TestSpatial:
    def test_spatial_1(self, validator_with_adatas):
        validator = validator_with_adatas 
        ...
        assert ...

    def test_spatial_2(self, validator_with_adatas):
        validator = validator_with_adatas 
        ...
        assert ...

Remember to include the self argument for each method in the test class

Can create your own custom interables/lists for any given test 
Or just import some of the constant lists from below
"""
import anndata as ad
import gc
import os
import pytest
from pathlib import Path
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad

# returns absolute path of fixtures directory
FIXTURES_ROOT = Path(__file__).absolute().parent

# few standard lists to select fixtures, can also whitelist for any given test
ALL_H5ADS = [
    file.name for file in FIXTURES_ROOT.iterdir() 
        if file.suffix == ".h5ad"
]
SPATIAL_H5ADS = [
    file for file in ALL_H5ADS 
        if any(included in file for included in [
            "slide_seq",
            "visium"
        ])
]
NON_SPATIAL_H5ADS = [
    file for file in ALL_H5ADS 
        if not any(excluded in file for excluded in [
            "slide_seq",
            "visium"
        ])
]
MULTISPECIES_H5ADS = [
    file for file in ALL_H5ADS 
        if not any(excluded in file for excluded in [
            "human", 
            "mouse"
        ])
]
MODEL_ORGANISM_H5ADS = [
    file for file in MULTISPECIES_H5ADS 
        if any(included in file for included in [
            "fly", 
            "worm", 
            "zebrafish"
        ])
]

# helper function for visium datasets to get library_id
def get_library_id(adata):
    return [key for key in adata.uns['spatial'].keys() if 'is_single' not in key][0]


# will add better check for file, maybe to automatically download as well
if not os.path.isfile(FIXTURES_ROOT / "visium_v2_11mm_human.h5ad"):
    raise FileNotFoundError('This file lives in S3, please download before running tests')


# base fixture to allow for all h5ads by default, whitelist with decorator
@pytest.fixture(params=ALL_H5ADS)
def test_h5ads(request):
    yield request.param


# fixture exported to other tests, yields and therefor tests with each h5ad
@pytest.fixture
def validator_with_adatas(test_h5ads) -> Validator:
    gc.collect()
    validator = Validator()
    # read_h5ad from cellxgene-schema; lazily loads matrices/layers
    # see dask_test.py or 5-3_pr1169_csr_matrix.py for how to load matrix and manipulate for testing
    validator.adata = read_h5ad(FIXTURES_ROOT / test_h5ads)
    yield validator


@pytest.fixture
def label_writer(validator_with_adatas: Validator) -> AnnDataLabelAppender:
    gc.collect()
    validator = validator_with_adatas
    validator.validate_adata()
    yield AnnDataLabelAppender(validator.adata)
