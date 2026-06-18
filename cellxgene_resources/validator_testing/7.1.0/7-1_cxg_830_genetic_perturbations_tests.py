# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
import pytest
import sys
from typing import Any
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


def get_first_key(dictionary: dict[str, Any]) -> str:
    return [key for key in dictionary][0]


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


    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_guide_h5ads_with_na_passes(self, index):

        # one "na" is valid, along with "no perturbations"

        for column, value in zip(GENETIC_OBS_COLS, ["na", "no perturbations"]):
            self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(value)

        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = "na"
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.na.value

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
class TestGeneticPerturbationObsFails:
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


    def test_all_nas_fails(self):
        
        # CURRENTLY PASSES
        # all nas for perturb_id should fail

        self.validator.adata.obs["genetic_perturbation_id"] = "na"
        self.validator.adata.obs["genetic_perturbation_strategy"] = GeneticPerturbationStrategy.na.value
        for column in GENETIC_OBS_COLS:
            self.validator.adata.obs[column] = self.validator.adata.obs[column].astype("category")

        assert self.validator.adata.obs["genetic_perturbation_strategy"].unique() == [GeneticPerturbationStrategy.na.value]
        assert self.validator.adata.obs["genetic_perturbation_id"].unique() == ["na"]
        assert (self.validator.adata.obs["genetic_perturbation_id"] == "na").all()

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


    @pytest.mark.parametrize("delimiter", [",", "|", "||", "-", "_", " "])
    @pytest.mark.parametrize("num_ids", [num for num in range(2, 6)])
    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_other_delimiters_fail(self, delimiter, num_ids, index):
        
        # other delimiters should fail
        # error message might be too non-specific?

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        keys = keys[:num_ids]
        new_id = delimiter.join(keys)
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' does not match any key in uns['genetic_perturbations']."
        ) in self.validator.errors


    @pytest.mark.parametrize("num_repeats", [num for num in range(2, 6)])
    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_duplicate_id_fails(self, num_repeats, index):
        
        # id with delimiter and two duplicate ids should fail

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        first_key = keys[0]
        repeats = [first_key] * num_repeats
        new_id = DELIMITER.join(repeats)
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' contains duplicates."
        ) in self.validator.errors


    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_id_out_of_order_fails(self, index):
        
        # id not in order fails

        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        first_keys = keys[:5]
        new_id = DELIMITER.join(sorted(first_keys, reverse=True))
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' is not in ascending lexical order."
        ) in self.validator.errors


    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_id_not_in_uns_fails(self, index):
        
        # id not in uns fails

        new_id = "test"
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{new_id}' in 'genetic_perturbation_id' does not match any key in uns['genetic_perturbations']."
        ) in self.validator.errors


    @pytest.mark.parametrize("index", [0, 10, 19])
    def test_non_id_in_multi_id_fails(self, index):
        
        # id not in uns fails when part of a multi-id

        new_id = "test"
        keys = list(self.validator.adata.uns["genetic_perturbations"].keys())
        several_ids = [keys[0], new_id, keys[1]]
        new_id = DELIMITER.join(sorted(several_ids, reverse=True))
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].cat.add_categories(new_id)
        obs_index = self.validator.adata.obs.index[index]
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_id"] = new_id
        self.validator.adata.obs.loc[obs_index, "genetic_perturbation_strategy"] = GeneticPerturbationStrategy.interference_screen.value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'test' in 'genetic_perturbation_id' does not match any key in uns['genetic_perturbations']."
        ) in self.validator.errors


    def test_just_id_as_na_fails(self):
        
        # OTHER ERRORS RELATING TO UNS AND MISSING COLUMNS, NOT ALL NAs
        # with just obs genetic perturbation id all as na, should fail

        self.validator.adata.obs["genetic_perturbation_id"] = "na"
        self.validator.adata.obs["genetic_perturbation_id"] = self.validator.adata.obs["genetic_perturbation_id"].astype("category")
        self.validator.adata.obs.drop(columns="genetic_perturbation_strategy", inplace=True)

        assert self.validator.adata.obs["genetic_perturbation_id"].unique() == ["na"]
        assert (self.validator.adata.obs["genetic_perturbation_id"] == "na").all()
        assert "genetic_perturbation_strategy" not in self.validator.adata.obs.columns

        expected_errors = [
            (
                "ERROR: Dataframe 'obs' is missing column 'genetic_perturbation_strategy'."
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


    @pytest.mark.parametrize("value", ["test", "control guide", "interference"])
    def test_other_strategy_fails(self, value):
        
        # strategies outside of the enum should fail

        self.validator.adata.obs["genetic_perturbation_strategy"] = value
        self.validator.adata.obs["genetic_perturbation_strategy"] = self.validator.adata.obs["genetic_perturbation_strategy"].astype("category")


        self.validator.validate_adata()
        start_of_error = (
            "ERROR: Column 'genetic_perturbation_strategy' in dataframe 'obs' "
            f"contains invalid values '['{value}']'. Values must be one of ['control', "
            "'CRISPR activation screen', 'CRISPR interference screen', 'CRISPR knockout mutant', "
            "'CRISPR knockout screen'] when 'genetic_perturbation_id' is in"
        )
        assert not self.validator.is_valid
        assert self.validator.errors[0].startswith(start_of_error)


@pytest.mark.parametrize("yield_guide_file_name", GUIDE_H5ADS)
class TestGeneticUnsPerturbationsFail:
    @pytest.fixture(autouse=True)
    def setup(self, yield_guide_validator):
        self.validator = yield_guide_validator


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


    @pytest.mark.parametrize("second_key", ["protospacer_sequence", "protospacer_adjacent_motif", "role"])
    @pytest.mark.parametrize(
        "value", [
            True,
            False,
            None,
            [],
            23,
            23.3,
        ]
    )
    def test_required_key_non_string_fails(self, value, second_key):
        
        # non-string values fail for required keys

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key][second_key] = value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{value}' in 'uns['genetic_perturbations']['{first_key}']['{second_key}']' "
            "is not valid, it must be a string."
        ) in self.validator.errors


    def test_na_as_key_fails(self):
        
        # na should not be allowed in dict

        key_list = [key for key in self.validator.adata.uns["genetic_perturbations"]]
        self.validator.adata.uns["genetic_perturbations"]["na"] = self.validator.adata.uns["genetic_perturbations"][key_list[0]]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: Key 'na' must not be present in 'uns.genetic_perturbations'."
        ) in self.validator.errors


    @pytest.mark.parametrize("char", ["'", "/", "\\" ",", " ", '"'])
    def test_not_allowed_id_characters_fails(self, char):
        
        # CURRENTLY PASSES
        # other characters should not be allowed in ids

        key_list = [key for key in self.validator.adata.uns["genetic_perturbations"]]
        new_id = f"test{char}"
        print(new_id, file=sys.__stderr__)
        self.validator.adata.uns["genetic_perturbations"][new_id] = self.validator.adata.uns["genetic_perturbations"][key_list[0]]

        assert new_id in self.validator.adata.uns["genetic_perturbations"]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: Key 'na' must not be present in 'uns.genetic_perturbations'."
        ) in self.validator.errors


    @pytest.mark.parametrize("value", ["Control", "na", "non-targeting", "Targeting", "interference"])
    def test_non_valid_role_fails(self, value):
        
        # only "control" and "targeting" should pass

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key]["role"] = value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{value}' in 'uns['genetic_perturbations']['{first_key}']['role']' "
            "is not valid. Allowed terms: ['control', 'targeting']."
        ) in self.validator.errors


    @pytest.mark.parametrize("key", ["targeted_feature", "sequence", "PAM", "guide"])
    def test_additional_keys_fail(self, key):
        
        # extra keys should not pass

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key][key] = "some data"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: 'uns['genetic_perturbations']['{first_key}']' contains unexpected key(s) ['{key}']. "
            "Additional key-value pairs MUST NOT be present."
        ) in self.validator.errors


    @pytest.mark.parametrize("key", ["role", "protospacer_sequence", "protospacer_adjacent_motif"])
    def test_missing_keys_fails(self, key):
        
        # missing required keys should not pass

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        del self.validator.adata.uns["genetic_perturbations"][first_key][key]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{key}' in 'uns['genetic_perturbations']['{first_key}']' is not present."
        ) in self.validator.errors


    @pytest.mark.parametrize("key", ["derived_genomic_regions", "derived_features"])
    def test_annotated_keys_fails(self, key):
        
        # annotated keys by CXG Discover should fail

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key][key] = "some data"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Key '{key}' is a reserved key name of 'uns.genetic_perturbations.{first_key}'. "
            "Remove it from h5ad and try again."
        ) in self.validator.errors


    @pytest.mark.parametrize(
        "sequence", [
            "AAAAAAAAAAAAA", 
            "RAAAAAAAAAAAAAA", 
            "AAAAAAAAAAAAAAAAAAAAAAA",
        ]
    )
    def test_protospacer_sequence_fails(self, sequence):
        
        # various sequences fail

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key]["protospacer_sequence"] = sequence

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{sequence}' in 'uns['genetic_perturbations']['{first_key}']['protospacer_sequence']' "
            "does not match the required pattern."
        ) in self.validator.errors


    @pytest.mark.parametrize(
        "sequence", [
            "SDFS", 
            "ZZZ", 
            "3' ZZZ",
            "3'",
            "Z",
            "",
            "5' NGG",
            "NGG",
        ]
    )
    def test_pam_fails(self, sequence):
        
        # various sequences fail

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key]["protospacer_adjacent_motif"] = sequence

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{sequence}' in 'uns['genetic_perturbations']['{first_key}']['protospacer_adjacent_motif']' "
            "does not match the required pattern."
        ) in self.validator.errors


    @pytest.mark.parametrize(
        "value", [
            True,
            False,
            None,
            [],
            23,
            23.3,
        ]
    )
    def test_intended_feature_non_dict_fails(self, value):
        
        # non-dict values fail for intended features

        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"]["intended_features"] = value

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{value}' in 'uns['genetic_perturbations']['intended_features']' "
            "is not valid, it must be a dictionary."
        ) in self.validator.errors


    @pytest.mark.parametrize("new_key", ["gene_id", "my_gene", "Oct2", "GFP"])
    def test_additional_intended_features_keys_fail(self, new_key):
        
        # non-dict values fail for intended features

        # first keys for all fixtures have intended features
        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        self.validator.adata.uns["genetic_perturbations"][first_key]["intended_features"][new_key] = "test_data"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Could not infer organism from feature ID '{new_key}' in "
            f"'uns['genetic_perturbations']['{first_key}']['intended_features']', make sure it is a valid ID."
        ) in self.validator.errors


    def test_ensembl_version_fail(self):
        
        # ensembl version number in intended_features fails
        # same error as wrong ensembl ID, this is probably alright

        # first keys for all fixtures have intended features
        first_key = get_first_key(self.validator.adata.uns["genetic_perturbations"])
        first_feature = get_first_key(self.validator.adata.uns["genetic_perturbations"][first_key]["intended_features"])
        version_id_key = f"{first_feature}.4"
        self.validator.adata.uns["genetic_perturbations"][first_key]["intended_features"][version_id_key] = ""

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Could not infer organism from feature ID '{version_id_key}' in "
            f"'uns['genetic_perturbations']['{first_key}']['intended_features']', make sure it is a valid ID."
        ) in self.validator.errors
