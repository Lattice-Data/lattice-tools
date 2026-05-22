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
import numpy as np
import os
import pandas as pd
import pytest
import random
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Generator
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad

# returns absolute path of fixtures directory
FIXTURES_ROOT = Path(__file__).absolute().parent
GUIDE_CSV_PATH = FIXTURES_ROOT / "guide_csvs"
DELIMITER = " || "
random.seed(1116)


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
GUIDE_H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
    "valid_zebrafish.h5ad",
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
def validator_with_adatas(test_h5ads) -> Generator[Validator, None, None]:
    gc.collect()
    validator = Validator()
    # read_h5ad from cellxgene-schema; lazily loads matrices/layers
    # see dask_test.py or 5-3_pr1169_csr_matrix.py for how to load matrix and manipulate for testing
    validator.adata = read_h5ad(FIXTURES_ROOT / test_h5ads)
    yield validator


@pytest.fixture
def label_writer(validator_with_adatas: Validator) -> Generator[AnnDataLabelAppender, None, None]:
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
def yield_atac_fixture_data(yield_atac_h5ads) -> Generator[AtacTestData, None, None]:
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

"""
Setup for fixtures with guide metadata for genetic pertubations:
csvs for uns metadata stored in guide_csvs/ directory:
    human_guides.csv based off of Marson Order 1 AN00024294
    mouse_guides.csv based off of Trapnell Scale GENE12-R117 NVUS202101701-66
    zebrafish_guides.csv based off of Weissman Order 3 AN00025549

"""

class GeneticPerturbationStrategy(Enum):
    na = "no perturbations"
    control = "control"
    activation_screen = "CRISPR activation screen"
    interference_screen = "CRISPR interference screen"
    knockout_mutant = "CRISPR knockout mutant"
    knockout_screen = "CRISPR knockout screen"


def make_uns_dict(df: pd.DataFrame) -> dict[str, Any]:
    """
    Take guide csv and create valid uns["genetic_perturbations"] dict
    """
    perturb_dict = {}
    for row in df.itertuples():
        perturb_dict[row.guide_id] = {
            "role": row.guide_role,
            "protospacer_sequence": row.guide_protospacer,
            "protospacer_adjacent_motif": row.guide_PAM,
        }
        if row.guide_target_gene_id not in {np.nan, "Ignore"}:
            perturb_dict[row.guide_id]["intended_features"] = {row.guide_target_gene_id: ""}

    return perturb_dict


# can have 1-10 guides per cell, weighting to skew towards 0
num_choices = [num for num in range(1, 11)]
weights = [100, 88, 78, 36, 25, 4, 1, 1, 1, 1]


def get_num_choices(
    choice_range: list[int] = num_choices, 
    weights: list[int] = weights,
) -> int:
    """
    Pick number of guides per cell within the choice range based on weights
    """
    return random.choices(choice_range, weights=weights)[0]


def create_genetic_perturbation_id(adata: ad.AnnData) -> str:
    """
    Assign genetic perturbation id to cells based on uns metadata
    Returns str formatted with CXG delimiter and sorted
    """
    perturbation_ids = adata.uns["genetic_perturbations"]
    if isinstance(perturbation_ids, dict):
        perturbation_ids = list(perturbation_ids.keys())
    k = get_num_choices()
    choices = sorted(random.sample(perturbation_ids, k))
    return DELIMITER.join(choices)



def get_perturbation_strategy(str_input: str, adata: ad.AnnData) -> str:
    perturb_roles: list[GeneticPerturbationStrategy] = [
        GeneticPerturbationStrategy.interference_screen,
        GeneticPerturbationStrategy.activation_screen,
        GeneticPerturbationStrategy.knockout_mutant,
        GeneticPerturbationStrategy.knockout_screen
    ]
    # need to provide list, otherwise will be set of characters
    control_set = set(["control"])

    if str_input is np.nan:
        return GeneticPerturbationStrategy.na.value

    input_list = str_input.split(DELIMITER)
    role_set = {adata.uns["genetic_perturbations"][genetic_id]["role"] for genetic_id in input_list}
    
    if role_set == control_set:
        return GeneticPerturbationStrategy.control.value
    
    # most real data probalby just one type, using all other types for testing
    return random.sample(perturb_roles, 1)[0].value


@pytest.fixture(params=GUIDE_H5ADS)
def yield_guide_file_name(request):
    yield request.param


@pytest.fixture
def yield_guide_validator(yield_guide_file_name) -> Generator[Validator, None, None]:
    gc.collect()
    # get species name from file
    species = yield_guide_file_name.split("_")[1].split(".")[0]
    adata = read_h5ad(FIXTURES_ROOT / yield_guide_file_name)
    validator = Validator()

    guide_df = pd.read_csv(GUIDE_CSV_PATH / f"{species}_guides.csv")
    guide_dict = make_uns_dict(guide_df)
    adata.uns["genetic_perturbations"] = guide_dict

    adata.obs["genetic_perturbation_id"] = "na"
    adata.obs["genetic_perturbation_id"] = adata.obs["genetic_perturbation_id"].apply(lambda x: create_genetic_perturbation_id(adata))
    adata.obs["genetic_perturbation_strategy"] = adata.obs["genetic_perturbation_id"].apply(lambda x: get_perturbation_strategy(x, adata))

    for column in ["genetic_perturbation_strategy", "genetic_perturbation_id"]:
        adata.obs[column] = adata.obs[column].astype("category")

    validator.adata = adata

    yield validator
