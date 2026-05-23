# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
import pytest
import sys
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    SPATIAL_H5ADS,              # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    MULTISPECIES_H5ADS,         # noqa: F401
    GUIDE_H5ADS,
    DELIMITER,
    GeneticPerturbationStrategy,# noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_guide_validator,      # noqa: F401
    yield_guide_file_name,      # noqa: F401
)

GENETIC_OBS_COLS = [
    "genetic_perturbation_id",
    "genetic_perturbation_strategy",
]
NON_GUIDE_H5ADS = [
    file
    for file in ALL_H5ADS
        if not any(excluded in file for excluded in [
            "human",
            "mouse",
            "zebrafish",
        ])
]

@pytest.mark.parametrize("yield_guide_file_name", GUIDE_H5ADS)
class TestGeneticPerturbationPasses:
    @pytest.fixture(autouse=True)
    def setup(self, yield_guide_validator):
        self.validator = yield_guide_validator


    def test_guide_h5ads_passes(self):

        # h5ads with guide metadata are valid, assert to make sure metadata is attached

        assert "genetic_perturbations" in self.validator.adata.uns
        assert len(self.validator.adata.uns["genetic_perturbations"]) > 1
        assert "genetic_perturbation_id" in self.validator.adata.obs.columns
        assert "genetic_perturbation_strategy" in self.validator.adata.obs.columns
        assert "control" in self.validator.adata.obs["genetic_perturbation_strategy"].unique()

        print(self.validator.adata.obs[["genetic_perturbation_id", "genetic_perturbation_strategy"]].head(), file=sys.__stderr__)
        print(self.validator.adata.obs[["genetic_perturbation_id", "genetic_perturbation_strategy"]].tail(), file=sys.__stderr__)

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_guide_h5ads_with_na_passes(self):

        # one "na" is valid, along with "no perturbations"

        for column, value in zip(GENETIC_OBS_COLS, ["na", "no perturbations"]):
            self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(value)

        first_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_id"] = "na"
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.na.value

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("yield_guide_file_name", NON_GUIDE_H5ADS)
class TestGeneticPerturbationNonGuideSpeciesFails:
    @pytest.fixture(autouse=True)
    def setup(self, yield_guide_validator):
        self.validator = yield_guide_validator

    def test_non_guides_fail(self):

        # fixture set up to attach human guide csv if no file found
        # all these should fail with other species guides attached to not human, mouse, or zebrafish

        species = self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: uns['organism_ontology_term_id'] is '{species}' but feature_ids "
            "are from [<SupportedOrganisms.HOMO_SAPIENS: 'NCBITaxon:9606'>]."
        ) in self.validator.errors


@pytest.mark.parametrize("yield_guide_file_name", GUIDE_H5ADS)
class TestGeneticPerturbationFails:
    @pytest.fixture(autouse=True)
    def setup(self, yield_guide_validator):
        self.validator = yield_guide_validator


    @pytest.mark.parametrize("column", GENETIC_OBS_COLS)
    def test_obs_colum_deletion_fails(self, column):
        
        # deleting one of two required obs columns fails

        self.validator.adata.obs.drop(columns=column, inplace=True)

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) >= 2
        assert (
            "ERROR: When adata.uns['genetic_perturbations'] is present, "
            "obs must contain 'genetic_perturbation_id' and 'genetic_perturbation_strategy'."
        ) in self.validator.errors

        assert (
            f"ERROR: Dataframe 'obs' is missing column '{column}'."
        ) in self.validator.errors


    def test_both_obs_columns_deleted_fails(self):
        
        # deleting both required obs columns fails

        self.validator.adata.obs.drop(columns=GENETIC_OBS_COLS, inplace=True)

        expected_errors = [
            (
                "ERROR: Dataframe 'obs' is missing column 'genetic_perturbation_id'."
            ),
            (
                "ERROR: When adata.uns['genetic_perturbations'] is present, obs must contain "
                "'genetic_perturbation_id' and 'genetic_perturbation_strategy'."
            ),
        ]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        for error in expected_errors:
            assert error in self.validator.errors


    def test_remove_uns_perturb_dict_fails(self):
        
        # deleting uns["genetic_perturbations"] fails

        del self.validator.adata.uns["genetic_perturbations"]

        expected_errors = [
            (
                "ERROR: Column 'genetic_perturbation_id' in dataframe 'obs' must not be "
                "present when uns['genetic_perturbations'] is not present."
            ),
            (
                "ERROR: obs['genetic_perturbation_id'] is present but adata.uns['genetic_perturbations'] "
                "is missing. When genetic perturbation columns exist in obs, "
                "uns['genetic_perturbations'] must be present."
            ),
        ]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        for error in expected_errors:
            assert error in self.validator.errors


    def test_all_nas_fails(self):
        
        # CURRENTLY PASSES
        # all nas for perturb_id should fail

        self.validator.adata.obs["genetic_perturbation_id"] = "na"
        self.validator.adata.obs["genetic_perturbation_strategy"] = GeneticPerturbationStrategy.na.value
        for column in GENETIC_OBS_COLS:
            self.validator.adata.obs[column] = self.validator.adata.obs[column].astype("category")

        assert self.validator.adata.obs["genetic_perturbation_strategy"].unique() == [GeneticPerturbationStrategy.na.value]
        assert self.validator.adata.obs["genetic_perturbation_id"].unique() == ["na"]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("delimiter", [",", "|", "||", "-", "_", " "])
    @pytest.mark.parametrize("num_ids", [num for num in range(2, 6)])
    def test_other_delimiters_fail(self, delimiter, num_ids):
        
        # other delimiters should fail
        # error message might be too non-specific?

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        keys = keys[:num_ids]
        new_id = delimiter.join(keys)
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        first_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' does not match any key in uns['genetic_perturbations']."
        ) in self.validator.errors


    @pytest.mark.parametrize("num_repeats", [num for num in range(2, 6)])
    def test_duplicate_id_fails(self, num_repeats):
        
        # id with delimiter and two duplicate ids should fail

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        first_key = keys[0]
        repeats = [first_key] * num_repeats
        new_id = DELIMITER.join(repeats)
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        first_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' contains duplicates."
        ) in self.validator.errors


    def test_id_out_of_order_fails(self):
        
        # id not in order fails

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        first_keys = keys[:5]
        new_id = DELIMITER.join(sorted(first_keys, reverse=True))
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        first_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[first_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' is not in ascending lexical order."
        ) in self.validator.errors
