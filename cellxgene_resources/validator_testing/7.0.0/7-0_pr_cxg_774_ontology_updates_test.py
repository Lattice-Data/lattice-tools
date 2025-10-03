"""
PR for this issue: no PR for single-cell-curation, updates on COG repo
See CXG-774 for more info on ontologies with new and deprecated terms

Testing conditions:

Should pass:
(Y) - new ontology updates pass on default fixtures
(Y) - new ontology combos work across non-spatial fixtures
(Y) - new fly terms should pass on fly fixture
(Y) - SHARE-seq should be valid with fragment files

Should not pass:
(Y) - deprecated terms are not valid
"""

import pytest
# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
from cellxgene_ontology_guide.ontology_parser import OntologyParser 
from cellxgene_schema.atac_seq import process_fragment
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_atac_fixture_data,    # noqa: F401
    yield_atac_h5ads,           # noqa: F401
    _to_anndata_file,           # noqa: F401
    to_temp_files               # noqa: F401
)

ONTOLOGY_PARSER = OntologyParser()
NEW_TERMS_DICT = {
    "assay_ontology_term_id": [
        "EFO:0900000",          # particle-templated instant partition sequencing
        "EFO:0900001",          # Asteria Single-cell RNA-seq Benchtop Kit
        "EFO:0900002",          # HIVE CLX Single-Cell RNAseq Solution
        "EFO:0022962",          # SHARE-seq
    ],
    "cell_type_ontology_term_id": [
        "CL:4072003",           # mature myelinating oligodendrocyte
        "CL:4033148",           # nodose ganglion TRPV1 neuron
        "CL:4033163",           # myenteric ganglion of small intestine nNOS/VIP neuron
    ],
    "disease_ontology_term_id": [
        "MONDO:7770001",        # argyrophilic grain disease
        "MONDO:1060151",        # schizoaffective bipolar disorder
        "MONDO:1060152",        # schizoaffective depressive disorder
    ],
    "tissue_ontology_term_id": [
        "UBERON:8920023",       # transverse pancreatic artery
        "UBERON:8920048",       # right gastroepiploic vein
        "UBERON:8920030",       # anterior inferior pancreaticoduodenal vein
    ],
}
DEPRECATED_TERMS_DICT = {
    "assay_ontology_term_id": [
        "EFO:1000398",          # Non-Functional Pancreatic Neuroendocrine Tumor
    ],
    "cell_type_ontology_term_id": [
        "CL:0000452",           # thyroid hormone secreting cell
        "CL:0000054",           # bone matrix secreting cell
        "CL:0000627",           # transporting cell
    ],
    "disease_ontology_term_id": [
        "MONDO:0035474",        # PIEZO1-related generalized lymphatic dysplasia with non-immune hydrops fetalis
        "MONDO:0007828",        # indifference to pain, congenital, autosomal dominant
    ],
    "self_reported_ethnicity_ontology_term_id": [
        "HANCESTRO:0518",       # Somali
        "HANCESTRO:0332",       # Sahrawi
    ],
}


def generate_tuple_list(term_dict: dict[str, list[str]]) -> list[tuple[str, str]]:
    """
    Use input dictionary to return list of tuples: (ontology column, ontology term)

    term_dict: dict of str ontology column name and values of list of ontology terms
    """
    results = []
    for column, ontology_terms in term_dict.items():
        for value in ontology_terms:
            results.append((column, value))
    return results


NEW_TERMS = generate_tuple_list(NEW_TERMS_DICT)
DEPRECATED_TERMS = generate_tuple_list(DEPRECATED_TERMS_DICT)


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestOntologyUpdatePasses:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_ontology_update_passes(self):

        # new ontology updates pass on default fixtures

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


# will use non spatial to avoid visium rules with assay
@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestOntologyUpdateValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize("column_ontology_terms", NEW_TERMS)
    def test_new_terms_pass(self, column_ontology_terms):

        # new ontology combos work across non-spatial fixtures

        column, ontology_term = column_ontology_terms
        self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(ontology_term)
        self.validator.adata.obs[column] = ontology_term

        # new SHARE-seq term is ATAC and needs nucleus suspension
        if ontology_term == "EFO:0022962":
            self.validator.adata.obs["suspension_type"] = "nucleus"
            self.validator.adata.obs["suspension_type"] = self.validator.adata.obs["suspension_type"].astype("category")

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("column_ontology_terms", DEPRECATED_TERMS)
    def test_deprecated_terms_fail(self, column_ontology_terms):

        # deprecated terms are not valid

        column, ontology_term = column_ontology_terms
        self.validator.adata.obs[column] = self.validator.adata.obs[column].cat.add_categories(ontology_term)
        self.validator.adata.obs[column] = ontology_term

        ontology_name = ONTOLOGY_PARSER._parse_ontology_name(ontology_term)
        organism = self.validator.adata.uns["organism_ontology_term_id"]

        # human fixtures report deprecation error instead of not valid value
        if ontology_name == "HANCESTRO" and organism != "NCBITaxon:9606":
            prefix_end = f"not a valid value of '{column}'."
        else:
            prefix_end = f"a deprecated term id of '{ontology_name}'."
        error_prefix = f"ERROR: '{ontology_term}' in '{column}' is {prefix_end}"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        for error in self.validator.errors:
            assert error.startswith(error_prefix)


@pytest.mark.parametrize("test_h5ads", ["valid_fly.h5ad"])
class TestFlyOntologyNewTerms:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize(
        "fbbt_term", 
        [
            "FBbt:20011413",    # adult lobula intrinsic neuron Li31
            "FBbt:20011378",    # adult central medulla intrinsic neuron Cm29
            "FBbt:00053545",    # adult olfactory receptor neuron Or35a ac3II
            "FBbt:00053548",    # adult olfactory receptor neuron Or85b ab11
        ]
    )
    def test_new_fly_terms(self, fbbt_term):

        # new fly terms should pass on fly fixture

        self.validator.adata.obs["cell_type_ontology_term_id"] = fbbt_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


def test_atac_fragment_with_share_seq(yield_atac_fixture_data, tmpdir):

    # SHARE-seq should be valid with fragment files

    test_data = yield_atac_fixture_data
    test_data.adata.obs["assay_ontology_term_id"] = "EFO:0022962"

    temp_files = to_temp_files(test_data, tmpdir)
    results = process_fragment(**temp_files)

    assert results == []
