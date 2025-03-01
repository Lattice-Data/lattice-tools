"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell/issues/721
https://github.com/chanzuckerberg/single-cell-curation/pull/1095/

For fragment validation, need location of h5ad and fragment file, so need
to provide temp file location of files before calling process-fragment

tmpfile is magic item for pytest that will use tmp dir for test then breakdown/
delete after

current fixture setup will by default use all fixtures in h5ads list
can specify specific fixtures with this:
@pytest.mark.parametrize("atac_h5ads", [{list specific h5ads}])
should apply this to other fixtures

TODO:
    maybe enum for files?
    flush out other tests
"""

import anndata as ad
import gc
import pandas as pd
import pytest
from cellxgene_schema.atac_seq import process_fragment
from dataclasses import dataclass
from pathlib import Path
from fixtures.valid_adatas import (
    FIXTURES_ROOT,
    read_h5ad
)

@dataclass
class AtacTestData:
    h5ad_file_name: str | Path
    adata: ad.AnnData
    fragment_df: pd.DataFrame
    fragment_file_name: str | Path


h5ads = [
    "valid_human.h5ad",
    "valid_mouse.h5ad"
]


def bundle_atac_test_data(h5ad_file_name) -> AtacTestData:
    fragment_file_name = h5ad_file_name.replace(".h5ad", "_fragments.tsv.gz")

    adata = read_h5ad(f"{FIXTURES_ROOT}/{h5ad_file_name}")

    fragments = pd.read_csv(
        f"{FIXTURES_ROOT}/{fragment_file_name}",
        sep="\t",
        header=None,
    )

    return AtacTestData(
        adata=adata,
        h5ad_file_name=h5ad_file_name,
        fragment_df=fragments,
        fragment_file_name=fragment_file_name
    )


@pytest.fixture(params=h5ads)
def atac_h5ads(request):
    return request.param


@pytest.fixture
def yeild_atac_fixture_data(atac_h5ads) -> AtacTestData:
    gc.collect()
    atac_fixture_data = bundle_atac_test_data(atac_h5ads)
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


# only mouse is valid at the moment, human fragment with duplicates
@pytest.mark.parametrize("atac_h5ads", ["valid_mouse.h5ad"])
def test_mock(yeild_atac_fixture_data, tmpdir):
    """
    get fixture
    modify
    call to_temp_files()
    run process_fragment
    assert
    """
    test_data = yeild_atac_fixture_data

    temp_files = to_temp_files(test_data, tmpdir)
    results = process_fragment(**temp_files)

    assert results == []


# human fails at the moment
@pytest.mark.parametrize("atac_h5ads", ["valid_human.h5ad"])
def test_mock_fails(yeild_atac_fixture_data, tmpdir):
    """
    get fixture
    modify
    call to_temp_files()
    run process_fragment
    assert
    """
    test_data = yeild_atac_fixture_data

    temp_files = to_temp_files(test_data, tmpdir)
    results = process_fragment(**temp_files)

    assert len(results) > 0


# edit human files to make valid
@pytest.mark.parametrize("atac_h5ads", ["valid_human.h5ad"])
def test_mock_with_edits(yeild_atac_fixture_data, tmpdir):
    """
    get fixture
    modify
    call to_temp_files()
    run process_fragment
    assert
    """
    # arrange
    test_data = yeild_atac_fixture_data

    test_data.adata.obs['is_primary_data'] = True
    test_data.fragment_df = test_data.fragment_df.drop_duplicates(keep="first")

    temp_files = to_temp_files(test_data, tmpdir)

    # act
    results = process_fragment(**temp_files)

    # assert
    assert results == []


# only human dataset with errors
@pytest.mark.parametrize("atac_h5ads", ["valid_human.h5ad"])
def test_against_errors(yeild_atac_fixture_data, tmpdir):
    """
    get fixture
    modify
    call to_temp_files()
    run process_fragment
    assert
    """
    test_data = yeild_atac_fixture_data

    temp_files = to_temp_files(test_data, tmpdir)
    results = process_fragment(**temp_files)

    assert "Fragment file has duplicate rows." in results
    assert "Anndata.obs.is_primary_data must all be True." in results
