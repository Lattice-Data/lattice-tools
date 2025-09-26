"""
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1454

Testing conditions:
Should not pass

Should pass
"""

import numpy as np
import pytest
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from testing_internals.utils import (
    make_valid_cell_line_fixture,
    NA_COLUMNS                  # noqa: F401
)
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
)

NON_HUMAN_H5ADS = [
    file for file in ALL_H5ADS 
        if not any(excluded in file for excluded in [
            "human", 
        ])
]
HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]
NON_CELL_LINE_TISSUES = [
    "organoid",
    "primary cell culture",
    "tissue",
]
ALL_TISSUE_TYPES = [
    "organoid",
    "primary cell culture",
    "tissue",
    "cell line",
]


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestEthnicityFixturesWork:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_tissue_cell_line_ethnicity_na_passes(self):

        # tissue_type == "cell line" with na for ethnicity term

        self.validator.adata = make_valid_cell_line_fixture(self.validator.adata)
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", NON_HUMAN_H5ADS)
class TestEthnicityValidationNonHuman:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_non_human_na_passes(self):

        # non-human fixtures must be na
        # default fixtures should be ok, will more thoroughly assert due to refactoring
        # done on this PR

        assert self.validator.adata.uns["organism_ontology_term_id"] != "NCBITaxon:9606"
        assert self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"].unique()[0] == "na"
        assert len(self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"].unique()) == 1
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", HUMAN_H5ADS)
class TestEthnicityValidationHuman:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_unknown_passes(self):

        # all unknown should pass

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize(
        "ethnicity_term", 
        [
            "HANCESTRO:0861",   # Black British
            "HANCESTRO:0847",   # Asian
            "HANCESTRO:0852",   # Middle Eastern
            "HANCESTRO:0846",   # Native American
            "HANCESTRO:0862",   # Asian American
            "AfPO:0000365",     # Egyptian, this term is on OLS, might make it to CXG with COG update
        ]
    )
    def test_current_cxg_forbidden_passes(self, ethnicity_term):

        # not currently allowed but are geography or ethnicity category terms
        # AfPO term does not pass, unclear if this will make it with latest HANCESTRO COG update

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = ethnicity_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_na_fails_non_cell_line(self, tissue_type):

        # all na fails with non cell line tissues

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id "
            "of 'HANCESTRO, AFPO'. When 'organism_ontology_term_id' is 'NCBITaxon:9606' (Homo sapiens) "
            "and 'tissue_type' is not 'cell line', self_reported_ethnicity_ontology_term_id MUST be "
            "formatted as one or more AfPO or HANCESTRO terms that are descendants of 'HANCESTRO:0601' "
            "for ethnicity category or 'HANCESTRO:0602' for geography-based population category, in "
            "ascending lexical order with the delimiter ` || `, or 'unknown' if unavailable."
        ) in self.validator.errors


    @pytest.mark.parametrize(
        "ethnicity_term", 
        [
            "HANCESTRO:0004",   # ancestory category
            "HANCESTRO:0005",   # European ancestory
            "HANCESTRO:0306",   # admixed ancestory
            "HANCESTRO:0632",   # reference population
            "HANCESTRO:0730",   # Ami in Taiwan
        ]
    )
    def test_non_valid_ethnicity_child_fails(self, ethnicity_term):

        # child terms outside of ethnicity and geography based population categories fail

        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = ethnicity_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{ethnicity_term}' in 'self_reported_ethnicity_ontology_term_id' is not an allowed term id. "
            "When 'organism_ontology_term_id' is 'NCBITaxon:9606' (Homo sapiens) "
            "and 'tissue_type' is not 'cell line', self_reported_ethnicity_ontology_term_id MUST be "
            "formatted as one or more AfPO or HANCESTRO terms that are descendants of 'HANCESTRO:0601' "
            "for ethnicity category or 'HANCESTRO:0602' for geography-based population category, in "
            "ascending lexical order with the delimiter ` || `, or 'unknown' if unavailable."
        ) in self.validator.errors


    @pytest.mark.parametrize("tissue_type", NON_CELL_LINE_TISSUES)
    def test_tissue_not_cellline_ethnicity_one_na_invalid(self, tissue_type):

        # tissue_type != cell line and ethnicity with one random 'na' + normal ids

        self.validator.adata.obs["tissue_type"] = tissue_type
        self.validator.adata.obs["tissue_type"] = self.validator.adata.obs["tissue_type"].astype("category")

        if tissue_type == "primary cell culture":
            self.validator.adata.obs["tissue_ontology_term_id"] = "CL:0000617"

        random_index = np.random.randint(0, (self.validator.adata.obs.shape[0] - 1))
        self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = self.validator.adata.obs["self_reported_ethnicity_ontology_term_id"].cat.add_categories(["na"])
        self.validator.adata.obs.loc[self.validator.adata.obs.index[random_index], "self_reported_ethnicity_ontology_term_id"] = "na"
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: 'na' in 'self_reported_ethnicity_ontology_term_id' is not a valid ontology term id "
            "of 'HANCESTRO, AFPO'. When 'organism_ontology_term_id' is 'NCBITaxon:9606' (Homo sapiens) "
            "and 'tissue_type' is not 'cell line', self_reported_ethnicity_ontology_term_id MUST be "
            "formatted as one or more AfPO or HANCESTRO terms that are descendants of 'HANCESTRO:0601' "
            "for ethnicity category or 'HANCESTRO:0602' for geography-based population category, in "
            "ascending lexical order with the delimiter ` || `, or 'unknown' if unavailable."
        ) in self.validator.errors
