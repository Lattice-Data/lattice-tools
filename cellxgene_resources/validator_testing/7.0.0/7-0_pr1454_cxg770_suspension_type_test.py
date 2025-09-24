"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1454

Testing conditions:
Should pass
(Y) - new suspension_type rules pass on default fixtures
(Y) - new suspension_type combos work across non-spatial fixtures

Should not pass
(N) - invalid suspension_type combos are appropriate across non-spatial fixtures
(Y) - EFO:0010550 sci-rna-seq removed, all suspensions should pass
"""

import pytest
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from cellxgene_ontology_guide.ontology_parser import OntologyParser 
from fixtures.valid_adatas import (
    ALL_H5ADS,
    NON_SPATIAL_H5ADS,
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
)

ONTOLOGY_PARSER = OntologyParser()
SUSPENSION_VALUES_SET = {
    "na",
    "cell",
    "nucleus",
}
# parallel test running needs consistent ordering, sets are unordered and different each run
SUSPENSION_VALUES_LIST = sorted(list(SUSPENSION_VALUES_SET))


ASSAY_SUSPENSION_DICT = {
    "EFO:0008679": {"cell"},                # CEL-seq
    "EFO:0008877": {"cell"},                # Quartz-seq
    "EFO:0010010": {"cell", "nucleus"},     # CEL-seq2 and descendants
    "EFO:0030074": {"cell", "nucleus"},     # CEL-seq2 and descendants: SORT-seq
    "EFO:0009919": {"cell", "nucleus"},     # SPLiT-seq and descendants
    "EFO:0022603": {"cell", "nucleus"},     # SPLiT-seq and descendants: Parse Biosciences technology
    "EFO:0022600": {"cell", "nucleus"},     # SPLiT-seq and descendants: Parse Evercode Whole Transcriptome v1
    "EFO:0022601": {"cell", "nucleus"},     # SPLiT-seq and descendants: Parse Evercode Whole Transcriptome v2
    "EFO:0022602": {"cell", "nucleus"},     # SPLiT-seq and descendants: Parse Evercode Whole Transcriptome v3
    "EFO:0008953": {"cell"},                # STRT-seq and descendants
    "EFO:0022846": {"cell"},                # STRT-seq and descendants: 5' STRT-seq
    "EFO:0022845": {"cell"},                # STRT-seq and descendants: modified STRT-seq
}


def generate_tuple_list(create_invalid: bool = False) -> list[tuple[str, str]]:
    """
    Use the above dictionary to return list of tuples: (assay, suspension_type)
    create_invalid uses set difference to find invalid terms
    """
    results = []
    for assay, values_set in ASSAY_SUSPENSION_DICT.items():

        # creates non-valid set based on valid values in input set
        values = SUSPENSION_VALUES_SET.difference(values_set) if create_invalid else values_set
        for value in values:
            results.append((assay, value))
    return results


# use this for tests to get assay term and allowed values per assay term
VALID_ASSAY_SUSPENSION_VALUES = generate_tuple_list()
INVALID_ASSAY_SUSPENSION_VALUES = generate_tuple_list(create_invalid=True)


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestSuspensionTypePasses:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_suspension_type_passes(self):

        # new suspension_type rules pass on default fixtures

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


# will use non spatial to avoid visium rules with suspension
@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestSuspensionTypeValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize("assay_suspension", VALID_ASSAY_SUSPENSION_VALUES)
    def test_suspension_type_new_rules_pass(self, assay_suspension):

        # new suspension_type combos work across non-spatial fixtures

        assay, suspension_type = assay_suspension
        self.validator.adata.obs["assay_ontology_term_id"] = assay
        self.validator.adata.obs["suspension_type"] = suspension_type
        for col in ["assay_ontology_term_id", "suspension_type"]:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].astype("category")
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("assay_suspension", INVALID_ASSAY_SUSPENSION_VALUES)
    def test_suspension_type_new_rules_fails(self, assay_suspension):

        # invalid suspension_type combos are appropriate across non-spatial fixtures
        # currently STRT-seq descendants are valid with nucleus and na

        assay, invalid_suspension_type = assay_suspension
        self.validator.adata.obs["assay_ontology_term_id"] = assay
        self.validator.adata.obs["suspension_type"] = invalid_suspension_type
        for col in ["assay_ontology_term_id", "suspension_type"]:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].astype("category")
            self.validator.adata.obs[col] = self.validator.adata.obs[col].cat.remove_unused_categories()
        self.validator.validate_adata()

        # grab valid suspension types, sort set to list casting to match error string
        suspension_types = sorted(list(ASSAY_SUSPENSION_DICT[assay]))
        assert not self.validator.is_valid
        assert (
            f"ERROR: Column 'suspension_type' in dataframe 'obs' contains invalid values "
            f"'['{invalid_suspension_type}']'. Values must be one of {suspension_types} when "
            f"'assay_ontology_term_id' is in ['{assay}']"
        ) in self.validator.errors


    @pytest.mark.parametrize("suspension", SUSPENSION_VALUES_LIST)
    def test_removed_assay_passes(self, suspension):

        # EFO:0010550 sci-rna-seq removed, all suspensions should pass

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0010550"
        self.validator.adata.obs["suspension_type"] = suspension
        for col in ["assay_ontology_term_id", "suspension_type"]:
            self.validator.adata.obs[col] = self.validator.adata.obs[col].astype("category")
            self.validator.adata.obs[col] = self.validator.adata.obs[col].cat.remove_unused_categories()
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []
