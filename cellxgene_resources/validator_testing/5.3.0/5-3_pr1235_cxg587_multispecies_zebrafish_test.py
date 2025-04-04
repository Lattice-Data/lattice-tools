"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1103
https://github.com/chanzuckerberg/single-cell-curation/pull/1235/
"""

import pytest
from fixtures.valid_adatas import (
    MULTISPECIES_H5ADS,
    test_h5ads,
    validator_with_adatas
)
 
DEV_ERROR_SUFFIX = (
    "When 'organism_ontology_term_id' is 'NCBITaxon:7955' (Danio rerio), "
    "'development_stage_ontology_term_id' MUST be the most accurate descendant of 'ZFS:0100000' and it "
    "MUST NOT be 'ZFS:0000000' for Unknown. The str 'unknown' is acceptable."
)

CELL_CULTURE_ERROR_SUFFIX = (
    "When 'tissue_type' is 'cell culture', 'tissue_ontology_term_id' MUST "
    "follow the validation rules for cell_type_ontology_term_id."
)

# fixtures are valid
@pytest.mark.parametrize("test_h5ads", MULTISPECIES_H5ADS)
def test_all_valid(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []


@pytest.mark.parametrize("test_h5ads", ["valid_zebrafish.h5ad"])
class TestZebraFishOntologyTerms:
    @pytest.mark.parametrize(
        "cell_type_term,error_suffix",
        (
            pytest.param("CL:0000255", "in 'cell_type_ontology_term_id' is not allowed.", id="CL eukaryotic cell forbidden"),
            pytest.param("CL:0000257", "in 'cell_type_ontology_term_id' is not allowed.", id="CL Eumycetozoan cell forbidden"),
            pytest.param("CL:0000548", "in 'cell_type_ontology_term_id' is a deprecated term id of 'CL'.", id="CL animal cell forbidden"),      # error says this is a deprecated term
            pytest.param("ZFA:0009000", "in 'cell_type_ontology_term_id' is not an allowed term id.", id="ZFA cell term forbidden"),   # error slightly different, "not an allowed term id"
        )
    )
    def test_cl_forbidden_zebrafish(self, validator_with_adatas, cell_type_term, error_suffix):
        validator = validator_with_adatas
        validator.adata.obs["cell_type_ontology_term_id"] = cell_type_term 
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: '{cell_type_term}' {error_suffix}"
        ) in validator.errors

    @pytest.mark.parametrize(
        "cell_type_term",
        (
            pytest.param("FBbt:00058005", id="Fly term not allowed for zebrafish"),
            pytest.param("WBbt:0008611", id="Worm term not allowed for zebrafish"),
        )
    )
    def test_non_cl_forbidden_for_zebrafish(self, validator_with_adatas, cell_type_term):
        validator = validator_with_adatas
        validator.adata.obs["cell_type_ontology_term_id"] = cell_type_term 
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: cell_type_ontology_term_id must be a valid CL term. If organism is NCBITaxon:6239,"
            " it can be a valid CL term or a valid WBbt term. If organism is NCBITaxon:7955, it can be a "
            "valid CL term or a valid ZFA term. If organism is NCBITaxon:7227, it can be a valid CL term or "
            "a valid FBbt term."
        ) in validator.errors

    @pytest.mark.parametrize(
        "cell_type_term",
        (
            pytest.param("CL:0000169", id="CL term for type B pancreatic cell"),
            pytest.param("ZFA:0009217", id="ZFA term for digestive enzyme secreting cell"),
        )
    )
    def test_cl_in_tissue_term_for_cell_culture(self, validator_with_adatas, cell_type_term):
        validator = validator_with_adatas
        validator.adata.obs["tissue_type"] = "cell culture"
        validator.adata.obs["tissue_type"] = validator.adata.obs["tissue_type"].astype("category")
        validator.adata.obs["tissue_ontology_term_id"] = cell_type_term 
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    @pytest.mark.parametrize(
        "cell_type_term,error_variable",
        (
            pytest.param("CL:0000255", "is not allowed.", id="CL eukaryotic cell forbidden"),
            pytest.param("CL:0000257", "is not allowed.", id="CL Eumycetozoan cell forbidden"),
            pytest.param("CL:0000548", "is a deprecated term id of 'CL'.", id="CL animal cell forbidden"),      # error says this is a deprecated term
            pytest.param("ZFA:0009000", "is not an allowed term id.", id="ZFA cell term forbidden"),   # error slightly different, "not an allowed term id"
        )
    )
    def test_cl_forbidden_cell_culture(self, validator_with_adatas, cell_type_term, error_variable):
        validator = validator_with_adatas
        validator.adata.obs["tissue_type"] = "cell culture"
        validator.adata.obs["tissue_ontology_term_id"] = cell_type_term 
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: '{cell_type_term}' in 'tissue_ontology_term_id' {error_variable} {CELL_CULTURE_ERROR_SUFFIX}"
        ) in validator.errors

    @pytest.mark.parametrize(
        "cell_type_term",
        (
            pytest.param("FBbt:00058005", id="Fly term not allowed for zebrafish"),
            pytest.param("WBbt:0008611", id="Worm term not allowed for zebrafish"),
        )
    )
    def test_non_cl_forbidden_for_zebrafish_cell_culture(self, validator_with_adatas, cell_type_term):
        validator = validator_with_adatas
        validator.adata.obs["tissue_type"] = "cell culture"
        validator.adata.obs["tissue_ontology_term_id"] = cell_type_term 
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. "
            "If organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. "
            "If organism is NCBITaxon:7955, it can be a valid UBERON term or a valid ZFA term. If organism is "
            "NCBITaxon:7227, it can be a valid UBERON term or a valid FBbt term. When tissue_type is cell culture, "
            "tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id."
        ) in validator.errors

    @pytest.mark.parametrize(
        "dev_term,error_variable",
        (
            pytest.param("ZFS:0100000", "is not an allowed term id.", id="Forbidden zebrafish dev term for zebrafish stage"),
            pytest.param("ZFS:0000000", "is not allowed.", id="Forbidden zebrafish dev term for zebrafish 'Unknown'"),
            # typo in error suffix for ZFA instead of ZFS
            pytest.param("FBdv:00007014", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term for adult age in days"),
            pytest.param("FBdv:00005259", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term for development stage"),
            pytest.param("FBdv:00007012", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term for life stage"),
            pytest.param("HsapDv:0000002", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, using HsapDv"),
            pytest.param("MmusDv:0000002", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, using MmusDv"),
            pytest.param("UBERON:0000105", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, UBERON forbidden parent"),
            pytest.param("UBERON:0000071", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, UBERON excluded death term"),
            pytest.param("UBERON:0034920", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, UBERON infant stage"),
            pytest.param("WBls:0000669", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, using worm term"),
            pytest.param("WBls:0000803", "is not a valid ontology term id of 'ZFA'.", id="Forbidden zebrafish dev term, using worm forbidden term"),
        )
    )
    def test_dev_terms_for_zebrafish(self, validator_with_adatas, dev_term, error_variable):
        validator = validator_with_adatas
        validator.adata.obs["development_stage_ontology_term_id"] = dev_term 
        validator.validate_adata()
        assert not validator.is_valid
        # looks like yaml schema definition file adds weird formatting into error message with line breaks
        assert (
            f"ERROR: '{dev_term}' in 'development_stage_ontology_term_id' {error_variable} "
            f"{DEV_ERROR_SUFFIX}"
       ) in validator.errors

    @pytest.mark.parametrize(
        "tissue_term,error_variable",
        (
            pytest.param("FBbt:10000000", "is not an allowed term id.", id="Zebrafish tissue term, forbidden parent"),
            pytest.param("FBbt:00007002", "is not allowed.", id="Zebrafish tissue term, forbidden cell and its descendants term"),
            pytest.param("UBERON:0001062", "is not an allowed term id.", id="Zebrafish tissue term, forbidden UBERON descendant term"),
            pytest.param("ZFA:0100000", "is not an allowed term id.", id="Zebrafish tissue term, forbidden zebrafish parent term"),
            pytest.param("ZFA:0001093", "is not allowed.", id="Zebrafish tissue term, forbidden zebrafish excluded term"),
            pytest.param("ZFA:0009000", "is not allowed.", id="Zebrafish tissue term, forbidden zebrafish cell term"),
            pytest.param("WBbt:0005766", "is not an allowed term id.", id="Zebrafish tissue term, forbidden worm anatomy parent term"),
            pytest.param("WBbt:0007849", "is not allowed.", id="Zebrafish tissue term, forbidden worm anatomy hermaphrodite"),
            pytest.param("WBbt:0007850", "is not allowed.", id="Zebrafish tissue term, forbidden worm anatomy male"),
            pytest.param("WBbt:0008595", "is not allowed.", id="Zebrafish tissue term, forbidden worm anatomy female"),
            pytest.param("WBbt:0004017", "is not allowed.", id="Zebrafish tissue term, forbidden worm anatomy cell and descendants"),
            pytest.param("WBbt:0006803", "is not allowed.", id="Zebrafish tissue term, forbidden worm anatomy nucleus and descendants"),
        )
    )
    def test_tissue_terms_zebrafish(self, validator_with_adatas, tissue_term, error_variable):
        validator = validator_with_adatas
        validator.adata.obs["tissue_ontology_term_id"] = tissue_term 
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            f"ERROR: '{tissue_term}' in 'tissue_ontology_term_id' {error_variable} "
            "When 'tissue_type' is 'tissue' or 'organoid', 'tissue_ontology_term_id' "
            "must be a valid UBERON, ZFA, FBbt, or WBbt term."
       ) in validator.errors
