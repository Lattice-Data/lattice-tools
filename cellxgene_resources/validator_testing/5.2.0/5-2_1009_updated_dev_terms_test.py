"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/988
https://github.com/chanzuckerberg/single-cell-curation/pull/1009/
"""

import pytest
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas,
)

H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
]

@pytest.mark.parametrize("test_h5ads", H5ADS)
class TestUpdatedDevTerms:
    @pytest.mark.parametrize(
        "dev_term, organism",
        (
            pytest.param("HsapDv:0000000", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000001", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000004", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000082", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000235", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000080", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000174", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000256", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000083", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000084", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000081", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000085", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000236", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000086", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000204", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000089", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000088", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000087", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000090", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000092", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000091", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000094", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000093", "NCBITaxon:9606"),
            pytest.param("MmusDv:0000000", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000001", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000016", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000030", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000037", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000038", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000039", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000040", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000041", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000044", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000045", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000046", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000047", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000048", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000049", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000050", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000051", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000052", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000053", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000054", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000055", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000056", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000057", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000058", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000059", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000061", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000065", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000066", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000067", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000068", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000070", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000071", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000072", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000073", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000074", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000075", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000076", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000096", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000097", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000098", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000099", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000100", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000101", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000102", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000112", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000113", "NCBITaxon:10090"),
            pytest.param("CL:0000234", "NCBITaxon:10090"),
            pytest.param("EFO:0009922", "NCBITaxon:10090"),
            pytest.param("MONDO:0005178", "NCBITaxon:10090"),
            pytest.param("PATO:0000384", "NCBITaxon:10090"),
            pytest.param("HANCESTRO:0027", "NCBITaxon:10090"),
            pytest.param("UBERON:0000411", "NCBITaxon:10090"),
        )
    )
    def test_old_dev_terms_should_fail(self, validator_with_adatas, dev_term, organism):
        validator = validator_with_adatas
        validator.adata.obs["development_stage_ontology_term_id"] = dev_term
        validator.adata.obs["organism_ontology_term_id"] = organism
        for col in ["development_stage_ontology_term_id", "organism_ontology_term_id"]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category")
        validator.validate_adata()
        assert not validator.is_valid


    @pytest.mark.parametrize(
        "dev_term, organism",
        (
            pytest.param("HsapDv:0010000", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000262", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000260", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000264", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000273", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000261", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000265", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000270", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000271", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000268", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000266", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000226", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000258", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000267", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000272", "NCBITaxon:9606"),
            pytest.param("HsapDv:0000227", "NCBITaxon:9606"),
            pytest.param("MmusDv:0000114", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000115", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000116", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000117", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000118", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000119", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000120", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000121", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000122", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000123", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000124", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000125", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000126", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000127", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000128", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000129", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000130", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000131", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000132", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000133", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000134", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000135", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000136", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000137", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000138", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000139", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000140", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000141", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000142", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000143", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000144", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000145", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000146", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000147", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000148", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000149", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000150", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000151", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000152", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000153", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000154", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000155", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000156", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000157", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000158", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000159", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000160", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000161", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000162", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000163", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000164", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000165", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000166", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000167", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000168", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000169", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000170", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000171", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000172", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000173", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000174", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000175", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000176", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000177", "NCBITaxon:10090"),
            pytest.param("MmusDv:0000178", "NCBITaxon:10090"),
        )
    )
    def test_new_dev_terms(self, validator_with_adatas, dev_term, organism):
        validator = validator_with_adatas
        validator.adata.obs["development_stage_ontology_term_id"] = dev_term
        validator.adata.obs["organism_ontology_term_id"] = organism

        if organism == "NCBITaxon:9606":
            validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
        else:
            validator.adata.obs["self_reported_ethnicity_ontology_term_id"] = "na"

        for col in [
            "development_stage_ontology_term_id",
            "organism_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id"
        ]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category")

        validator.validate_adata()
        assert validator.is_valid
