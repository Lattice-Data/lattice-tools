"""
Will likely be several PRs for this issue
See CXG-831 for more info on ontologies with new and deprecated terms

Testing conditions:

Should pass:
( ) cell_type_ontology_term_id is absent when uns['is_pre_analysis'] is True
( ) no obsm keys when uns['is_pre_analysis'] is True
( ) uns['default_embedding'] can be present when uns['is_pre_analysis'] is False
( ) uns['is_pre_analysis'] value is bool (may not be needed if this is added by CELLxGENE Discover)
( ) pre-analysis checks work with correct entry point and not changing Validator flag outside of normal paths

Shouldn't pass: 
( ) cell_type_ontology_term_id is absent when uns['is_pre_analysis'] is False
( ) obsm keys when uns['is_pre_analysis'] is True
( ) uns['default_embedding'] present when uns['is_pre_analysis'] is False
( ) uns['is_pre_analysis'] value is other type beside bool (may not be needed if this is added by CELLxGENE Discover)
"""

from cellxgene_schema.validate import Validator
import pytest
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    SPATIAL_H5ADS,              # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_atac_fixture_data,    # noqa: F401
    yield_atac_h5ads,           # noqa: F401
    _to_anndata_file,           # noqa: F401
    to_temp_files               # noqa: F401
)

def set_to_preanalysis_state(validator: Validator, remove_obsm: bool = True) -> Validator:
    """
    Set flag and remove obs and uns fields to get to pre-analysis state
    """
    validator.pre_analysis_check_flag = True
    validator.adata.obs.drop(columns="cell_type_ontology_term_id", inplace=True)

    # option to remove to test for schema fields
    if remove_obsm:
        validator.adata.obsm = None

    # should always be removed for pre_analysis, purely a data portal flag for Explorer
    if "default_embedding" in validator.adata.uns:
        del validator.adata.uns["default_embedding"]

    return validator


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestIsPreAnalysisNonSpatialPasses:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_is_pre_analysis_passes(self):

        # new is_pre_analysis changes pass on default fixtures
    
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_remove_fields_passes(self):

        # FAILS: visium datasets raise KeyError at line 2310;
        # spatial checks look for cell_type_ontology_term_id but this is not present
        # debating bigger question if pre-analysis spatial data should be a thing

        self.validator = set_to_preanalysis_state(self.validator)

        assert "cell_type_ontology_term_id" not in self.validator.adata.obs.columns
        assert "default_embedding" not in self.validator.adata.uns
        # obsm is AxisArray object, when keys cleared, length should be 0
        assert len(self.validator.adata.obsm) == 0

        self.validator.validate_adata()

        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_keep_obsm_passes(self):

        # FAILS: visium datasets raise KeyError at line 2310;
        # spatial checks look for cell_type_ontology_term_id but this is not present
        # debating bigger question if pre-analysis spatial data should be a thing

        self.validator = set_to_preanalysis_state(self.validator, remove_obsm=False)

        assert "cell_type_ontology_term_id" not in self.validator.adata.obs.columns
        assert "default_embedding" not in self.validator.adata.uns
        assert len(self.validator.adata.obsm) > 0

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", SPATIAL_H5ADS)
class TestIsPreAnalysisSpatialPasses:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_is_pre_analysis_passes(self):

        # new is_pre_analysis changes pass on default fixtures

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_remove_all_fields_passes(self):

        # FAILS: visium datasets raise KeyError at line 2310;
        # spatial checks look for cell_type_ontology_term_id but this is not present
        # debating bigger question if pre-analysis spatial data should be a thing

        self.validator = set_to_preanalysis_state(self.validator)

        assert "cell_type_ontology_term_id" not in self.validator.adata.obs.columns
        assert "default_embedding" not in self.validator.adata.uns
        # obsm is AxisArray object, when keys cleared, length should be 0
        assert len(self.validator.adata.obsm) == 0

        self.validator.validate_adata()

        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_remove_some_fields_passes(self):

        # FAILS: visium datasets raise KeyError at line 2310;
        # spatial checks look for cell_type_ontology_term_id but this is not present
        # debating bigger question if pre-analysis spatial data should be a thing

        self.validator = set_to_preanalysis_state(self.validator, remove_obsm=False)

        assert "cell_type_ontology_term_id" not in self.validator.adata.obs.columns
        assert "default_embedding" not in self.validator.adata.uns
        assert len(self.validator.adata.obsm) > 0

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestIsPreAnalysisFails:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_is_pre_analysis_key_in_uns_fails(self):

        # add uns key to True and should fail

        self.validator.adata.uns["is_pre_analysis"] = True
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: Key 'is_pre_analysis' is a reserved key name of 'uns'. "
            "Remove it from h5ad and try again."
        ) in self.validator.errors


# will use non spatial to avoid visium rules with assay
# @pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
# class TestOntologyUpdateValidation:
#     @pytest.fixture(autouse=True)
#     def setup(self, validator_with_adatas):
#         self.validator = validator_with_adatas
#
#
#     @pytest.mark.parametrize("column_ontology_terms", NEW_TERMS)
#     def test_new_terms_pass(self, column_ontology_terms):
#
#         # new ontology combos work across non-spatial fixtures
#
#         column, ontology_term = column_ontology_terms
#         self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(ontology_term)
#         self.validator.adata.obs[column] = ontology_term
#
#         # new SHARE-seq term is ATAC and needs nucleus suspension
#         if ontology_term == "EFO:0022962":
#             self.validator.adata.obs["suspension_type"] = "nucleus"
#             self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
#
#         # if spatial assay, set suspension to na, EFO term is 'spatial transcriptomics'
#         if ontology_term in ONTOLOGY_PARSER.get_term_descendants("EFO:0008994"):
#             self.validator.adata.obs["suspension_type"] = "na"
#             self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")
#
#         self.validator.validate_adata()
#         assert self.validator.is_valid
#         assert self.validator.errors == []
#
#
#     @pytest.mark.parametrize("column_ontology_terms", DEPRECATED_TERMS)
#     def test_deprecated_terms_fail(self, column_ontology_terms):
#
#         # deprecated terms are not valid
#
#         column, ontology_term = column_ontology_terms
#         self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(ontology_term)
#         self.validator.adata.obs[column] = ontology_term
#
#         ontology_name = ONTOLOGY_PARSER._parse_ontology_name(ontology_term)
#         organism = self.validator.adata.uns["organism_ontology_term_id"]
#
#         # human fixtures report deprecation error instead of not valid value
#         if ontology_name == "HANCESTRO" and organism != "NCBITaxon:9606":
#             prefix_end = f"not a valid value of '{column}'."
#         else:
#             prefix_end = f"a deprecated term id of '{ontology_name}'."
#         error_prefix = f"ERROR: '{ontology_term}' in '{column}' is {prefix_end}"
#
#         self.validator.validate_adata()
#         assert not self.validator.is_valid
#         for error in self.validator.errors:
#             assert error.startswith(error_prefix)
#
#
# @pytest.mark.parametrize("test_h5ads", ["valid_fly.h5ad"])
# class TestFlyOntologyNewTerms:
#     @pytest.fixture(autouse=True)
#     def setup(self, validator_with_adatas):
#         self.validator = validator_with_adatas
#
#
#     @pytest.mark.parametrize(
#         "fbbt_term", 
#         [
#             "FBbt:20011413",    # adult lobula intrinsic neuron Li31
#             "FBbt:20011378",    # adult central medulla intrinsic neuron Cm29
#             "FBbt:00053545",    # adult olfactory receptor neuron Or35a ac3II
#             "FBbt:00053548",    # adult olfactory receptor neuron Or85b ab11
#         ]
#     )
#     def test_new_fly_terms(self, fbbt_term):
#
#         # new fly terms should pass on fly fixture
#
#         self.validator.adata.obs["cell_type_ontology_term_id"] = fbbt_term
#         self.validator.validate_adata()
#         assert self.validator.is_valid
#         assert self.validator.errors == []
#
#
# def test_atac_fragment_with_share_seq(yield_atac_fixture_data, tmpdir):
#
#     # SHARE-seq should be valid with fragment files
#
#     test_data = yield_atac_fixture_data
#     test_data.adata.obs["assay_ontology_term_id"] = "EFO:0022962"
#
#     temp_files = to_temp_files(test_data, tmpdir)
#     results = process_fragment(**temp_files)
#
#     assert results == []
