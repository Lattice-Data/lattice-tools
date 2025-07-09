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

@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
def test_spatial(validator_with_adatas):
    validator = validator_with_adatas 
    ...
    assert ...

Now only the h5ads in SPATIAL_H5ADS will be tested

Tests can be grouped into classes to use the same h5ad subset for multiple tests without
having to repeat the decoration for each test method/function:

@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
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
import pandas as pd
import pytest
from dataclasses import dataclass
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
    return [key for key in adata.uns["spatial"].keys() if "is_single" not in key][0]


# will add better check for file, maybe to automatically download as well
if not os.path.isfile(FIXTURES_ROOT / "visium_v2_11mm_human.h5ad"):
    raise FileNotFoundError(
        "This 2.3 GB file lives on Google Drive, see the QA Fixture Info sheet;"
        " please download before running tests"
    )


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


"""
Fixture setup for ATAC/Multiome testing.
Need to provide paths of fragment file and h5ad, so this adds further setup
beyond just referencing the in-memory AnnData object, like for most other 
validator testing

AtacTestData holds in-memory Anndata object and fragment df, and names of 
each file

Each test should take in an AtacTestData object, modify AnnData and/or 
fragment file for the test setup, then call to_temp_files() to save these
two objects to a temp dir

to_temp_files() saves objects in temp dir and returns temp paths; these paths
can then be used by process_fragment() to run the test

See 5.3.0/5-3_pr1095_cxg602_atac_fragment_test.py::test_mock() as an example
of test setup
"""

@dataclass
class AtacTestData:
    h5ad_file_name: str | Path
    adata: ad.AnnData
    fragment_df: pd.DataFrame
    fragment_file_name: str | Path


ATAC_H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
]


def bundle_atac_test_data(h5ad_file_name: str) -> AtacTestData:
    fragment_file_name = h5ad_file_name.replace(".h5ad", "_fragments.tsv.gz")

    adata = read_h5ad(FIXTURES_ROOT / h5ad_file_name)
    # set to default valid for atac since fixtures used for non-atac tests
    adata.obs["assay_ontology_term_id"] = "EFO:0030059"
    adata.obs["is_primary_data"] = True

    # will see how this try block works for mixing/matching h5ads with valid fragment
    # fixtures; need to set appropriate fragment within the test
    # otherwise, will get various AttributeErrors as other functions attempt to access
    # None attributes
    try:
        fragments = pd.read_csv(
            FIXTURES_ROOT / fragment_file_name,
            sep="\t",
            header=None,
        )
    except FileNotFoundError:
        fragments = None

    return AtacTestData(
        adata=adata,
        h5ad_file_name=h5ad_file_name,
        fragment_df=fragments,
        fragment_file_name=fragment_file_name
    )


@pytest.fixture(params=ATAC_H5ADS)
def yield_atac_h5ads(request):
    yield request.param


@pytest.fixture
def yield_atac_fixture_data(yield_atac_h5ads) -> AtacTestData:
    gc.collect()
    atac_fixture_data = bundle_atac_test_data(yield_atac_h5ads)
    yield atac_fixture_data


def _to_anndata_file(atac_data: AtacTestData, tmp_path: Path) -> str | Path:
    tmp_file_name = tmp_path / atac_data.h5ad_file_name 
    atac_data.adata.write(tmp_file_name, compression="gzip")
    return tmp_file_name

    
def _to_fragment_file(atac_data: AtacTestData, tmp_path: Path) -> str | Path:
    tmp_file_name = tmp_path / atac_data.fragment_file_name
    atac_data.fragment_df.to_csv(
        tmp_file_name, 
        sep="\t", 
        header=False, 
        index=False, 
        compression="gzip"
    )
    return tmp_file_name


def to_temp_files(test_data: AtacTestData, tmp_path: Path | str) -> dict:
    return {
        "anndata_file": _to_anndata_file(test_data, tmp_path), 
        "fragment_file": _to_fragment_file(test_data, tmp_path)
    }
