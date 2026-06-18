# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
import pytest
import tempfile
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    SPATIAL_H5ADS,              # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    MULTISPECIES_H5ADS,         # noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_atac_fixture_data,    # noqa: F401
    yield_atac_h5ads,           # noqa: F401
    _to_anndata_file,           # noqa: F401
    to_temp_files               # noqa: F401
)
from cellxgene_schema.write_labels import AnnDataLabelAppender


DELIMITER = " || "
TERM_LIST = [
    "uniprot:P05112", 
    "anti-uniprot:Q99467", 
    "EFO:0002757", 
    "CHEBI:16412", 
    "EFO:0001702",
    "CHEBI:41774",
]
CORRECT_LABEL = "lipopolysaccharide || tamoxifen || temperature || high fat diet || anti-CD180_HUMAN || IL4_HUMAN"

def make_sorted_delimiter_str(input_terms = list[str], delimiter=DELIMITER) -> str:
    return delimiter.join(sorted(input_terms))


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestExperimentalConditionPasses:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize(
        "condition_value",
        [
            "uniprot:P05112", 
            "anti-uniprot:Q99467", 
            "EFO:0002757", 
            "CHEBI:16412", 
            "EFO:0001702",
            "CHEBI:41774",
        ]
    )
    def test_new_ontologies_pass(self, condition_value):

        # new experimental_condition passes on all fixtures
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = condition_value
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_one_na_passes(self):

        # na is allowed
    
        # set to all valid value
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = "CHEBI:16412"
        first_index = self.validator.adata.obs.index[0]
        # add just one na
        self.validator.adata.obs.loc[first_index, "experimental_condition_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_sorted_delimiter_str_passes(self):

        # lexical sorting with correct delimiter works

        multi_term_str = make_sorted_delimiter_str(TERM_LIST)
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = multi_term_str
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize(
        "condition_value",
        [
            "EFO:0002755", 
            "EFO:0803365", 
            "EFO:0002757", 
            "EFO:0009371", 
            "EFO:0001702",
            "EFO:0002758",
            "EFO:0002756", 
        ]
    )
    def test_efo_diet_fasting_temp_passes(self, condition_value):

        # new experimental_condition passes on all fixtures
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = condition_value
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestExperimentalPerturbationFails:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_all_na_fails(self):

        # all na not allowed
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0
        assert self.validator.errors[0].startswith(
            "ERROR: Column 'experimental_condition_ontology_term_id' in "
            "dataframe 'obs' MUST NOT be present when all values are 'na'."
        )


    def test_sorted_wrong_delimiter_str_fails(self):

        # lexical sorting with wrong delimiter fails

        multi_term_str = make_sorted_delimiter_str(TERM_LIST, ",")
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = multi_term_str
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors[0].startswith(
            "ERROR: 'CHEBI:16412,CHEBI:41774,EFO:0001702,EFO:0002757,anti-uniprot:Q99467,uniprot:P05112' "
            "in 'experimental_condition_ontology_term_id' is not a valid ontology term id of 'CHEBI, EFO, UniProt'."
        )


    def test_duplicate_multiterm_fails(self):

        # lexical sorting with duplication fails

        multi_term_str = make_sorted_delimiter_str([TERM_LIST[0], TERM_LIST[0]])
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = multi_term_str
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors[0].startswith(
            "ERROR: 'uniprot:P05112 || uniprot:P05112' in 'experimental_condition_ontology_term_id' contains duplicates."
        )


    def test_unsorted_multiterm_fails(self):

        # unsorted multiterm string fails

        multi_term_str = DELIMITER.join(TERM_LIST)
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = multi_term_str
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors[0].startswith(
            "ERROR: 'uniprot:P05112 || anti-uniprot:Q99467 || EFO:0002757 || CHEBI:16412 || "
            "EFO:0001702 || CHEBI:41774' in 'experimental_condition_ontology_term_id' is not "
            "in ascending lexical order."
        )


    @pytest.mark.parametrize(
        "condition_value",
        [
            "CHEBI:23367",
            "CHEBI:24431",
            "CHEBI:24835",
            "CHEBI:24867",
            "CHEBI:24870",
            "CHEBI:25212",
            "CHEBI:25367",
            "CHEBI:25699",
            "CHEBI:33238",
            "CHEBI:33259",
            "CHEBI:33497",
            "CHEBI:33595",
            "CHEBI:33674",
            "CHEBI:33675",
            "CHEBI:36342",
            "CHEBI:36357",
            "CHEBI:36358",
            "CHEBI:36914",
            "CHEBI:37577",
            "CHEBI:50906",
        ]
    )
    def test_chebi_restricted_fails(self, condition_value):

        # restricted CHEBI terms fail
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = condition_value
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0
        assert self.validator.errors[0].startswith(f"ERROR: '{condition_value}' in 'experimental_condition_ontology_term_id' is not allowed.")


    @pytest.mark.parametrize(
        "condition_value",
        [
            "CHEBI:36368",  # subatomic particle descendant, strange quark
            "CHEBI:50913",  # role descendant, fixative
        ]
    )
    def test_chebi_descendants_fail(self, condition_value):

        # restricted CHEBI terms fail
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = condition_value
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0
        assert self.validator.errors[0].startswith(f"ERROR: '{condition_value}' in 'experimental_condition_ontology_term_id' is not allowed.")


    @pytest.mark.parametrize(
        "condition_value",
        [
            "EFO:0022604",      # 10x 3' v4
            "EFO:0005316",      # sample pooling
        ]
    )
    def test_random_efo_fails(self, condition_value):

        # other efo terms should fail, including those allowed for assay_ontology_term_id
    
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = condition_value
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert self.validator.errors[0].startswith(
            f"ERROR: '{condition_value}' in 'experimental_condition_ontology_term_id' is not an allowed term id"
        )


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestLabelWriting:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_labels_are_correct(self):
        multi_term_str = make_sorted_delimiter_str(TERM_LIST)
        self.validator.adata.obs["experimental_condition_ontology_term_id"] = multi_term_str
        first_index = self.validator.adata.obs.index[0]
        second_index = self.validator.adata.obs.index[1]
        self.validator.adata.obs.loc[first_index, "experimental_condition_ontology_term_id"] = "na"

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []

        with tempfile.TemporaryDirectory() as temp_dir:
            label_class = AnnDataLabelAppender(self.validator.adata)
            labels_path = temp_dir + "labels.h5ad"
            write_successful = label_class.write_labels(labels_path)

            assert label_class.adata.obs.loc[first_index, "experimental_condition"] == "na"
            assert label_class.adata.obs.loc[second_index, "experimental_condition"] == CORRECT_LABEL
            assert write_successful
            assert not label_class.errors
