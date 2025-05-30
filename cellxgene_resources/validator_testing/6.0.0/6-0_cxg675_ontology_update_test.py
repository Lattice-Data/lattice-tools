"""
QA testing for this issue:
PR for this issue:

Should pass:
(N) CL: new terms ["CL:4052040", "CL:4070016"]
(N) MONDO: new terms ["MONDO:0700303", "MONDO:0700293"]
(N) UBERON: new terms ["UBERON:8600126", "UBERON:8920001"]

Should not pass:
(N) CL: deprecated terms ["CL:0000215", "CL:0000217"]
(N) MONDO: deprecated terms ["MONDO:0022577", "MONDO:0800112", "MONDO:0000414"]

"""

import pytest
import anndata as ad
from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.utils import read_h5ad
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)


NEW_TERMS = {
    "CL":["CL:4052040", "CL:4070016"],
    "MONDO":["MONDO:0700303", "MONDO:0700293"],
    "UBERON":["UBERON:8600126", "UBERON:8920001"]
}


DEPRECATED_TERMS = {
    "CL":["CL:0000215", "CL:0000217"],
    "MONDO":["MONDO:0022577", "MONDO:0800112", "MONDO:0000414"]
}

VISIUM_ASSAYS = ["EFO:0022857","EFO:0022858","EFO:0022860","EFO:0022859"]

@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestNewOntologyTerms:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_passes_new(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("new_term", NEW_TERMS["CL"])
    def test_new_CL_terms(self, new_term):

        # new CL terms -> valid   # currently not passing

        self.validator.adata.obs["cell_type_ontology_term_id"] = self.validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories(new_term)

        if  self.validator.adata.obs["assay_ontology_term_id"].unique()[0] in VISIUM_ASSAYS and self.validator.adata.uns["spatial"]["is_single"] == True:
            self.validator.adata.obs.loc[self.validator.adata.obs["in_tissue"] == 1, "cell_type_ontology_term_id"] = new_term

        else:
            self.validator.adata.obs["cell_type_ontology_term_id"] = new_term

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("new_term", NEW_TERMS["MONDO"])
    def test_new_MONDO_terms(self, new_term):

        # new MONDO terms -> valid   # currently not passing

        self.validator.adata.obs["disease_ontology_term_id"] = new_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("new_term", NEW_TERMS["UBERON"])
    def test_new_UBERON_terms(self, new_term):

        # new UBERON terms -> valid   # currently not passing

        self.validator.adata.obs["tissue_ontology_term_id"] = new_term
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestDeprecatedOntologyTerms:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_passes_dep(self):
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("dep_term", DEPRECATED_TERMS["CL"])
    def test_dep_CL_terms(self, dep_term):

        # deprecated CL terms -> invalid   # currently valid

        self.validator.adata.obs["cell_type_ontology_term_id"] = self.validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories(dep_term)

        if  self.validator.adata.obs["assay_ontology_term_id"].unique()[0] in VISIUM_ASSAYS and self.validator.adata.uns["spatial"]["is_single"] == True:
            self.validator.adata.obs.loc[self.validator.adata.obs["in_tissue"] == 1, "cell_type_ontology_term_id"] = dep_term

        else:
            self.validator.adata.obs["cell_type_ontology_term_id"] = dep_term

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{dep_term}' in 'cell_type_ontology_term_id' is not a valid ontology term id of 'CL, ZFA, FBbt, WBbt'.'"
            ) in self.validator.errors


    @pytest.mark.parametrize("dep_term", DEPRECATED_TERMS["MONDO"])
    def test_dep_MONDO_terms(self, dep_term):

        # deprecated MONDO terms -> invalid   # currently valid

        self.validator.adata.obs["disease_ontology_term_id"] = dep_term
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{dep_term}' in 'disease_ontology_term_id' is not a valid ontology term id of 'MONDO, PATO'. Individual terms 'PATO:0000461' (normal), 'MONDO:0021178' (injury) or descendant terms thereof, or descendant terms of 'MONDO:0000001' (disease) are allowed. Multiple terms are supported if in ascending lexical order with the delimiter ` || ` if all terms are valid MONDO terms."
        ) in self.validator.errors