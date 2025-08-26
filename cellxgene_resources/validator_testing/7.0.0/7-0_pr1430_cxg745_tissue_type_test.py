"""
QA testing for PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1430

Should pass:
- tissue_type == "cell line" → Note that since we don't have the cellosaurus changes in yet, the "cell line" requirements for tissue_ontology_term_id are not yet implemented.
- tissue_type == "primary cell culture"
- tissue_type == "organoid" + tissue_ontology_term_id == descendant of "UBERON:0001062"
- tissue_type == "tissue" + tissue_ontology_term_id == descendant of "UBERON:0001062"

Should not pass:
- tissue_type == "cell culture"
- tissue_type == "organoid" + tissue_ontology_term_id == "UBERON:0000922" for embryo
- tissue_type == "organoid" + tissue_ontology_term_id != descendant of "UBERON:0001062"
- tissue_type == "tissue" + tissue_ontology_term_id != descendant of "UBERON:0001062"
"""

import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

MODEL_ORG_CL_TERMS = {
     "NCBITaxon:6239":"WBbt:0006925",
     "NCBITaxon:7955":"ZFA:0009129",
     "NCBITaxon:7227":"FBbt:00053322"
}

@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestTissueTypeValidation:
     @pytest.fixture(autouse=True)
     def setup(self, validator_with_adatas):
          self.validator = validator_with_adatas


     def test_valid(self):
          self.validator.validate_adata()
          assert self.validator.is_valid


     def test_tissue_type_cell_line_pass(self):

          # tissue_type == cell line -> pass

          self.validator.adata.obs['tissue_type'] = 'cell line'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.validate_adata()
          assert self.validator.is_valid


     def test_tissue_type_primarycellculture_pass(self):

          # tissue_type == "primary cell culture" -> pass

          self.validator.adata.obs['tissue_type'] = 'primary cell culture'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          if self.validator.adata.uns["organism_ontology_term_id"] in MODEL_ORG_CL_TERMS.keys():
               self.validator.adata.obs['tissue_ontology_term_id'] = MODEL_ORG_CL_TERMS[self.validator.adata.uns["organism_ontology_term_id"]]
          else:
               self.validator.adata.obs['tissue_ontology_term_id'] = 'CL:0000738'
          self.validator.validate_adata()
          assert self.validator.is_valid






     def test_cell_culture_fail(self):

          assert "ERROR: When tissue_type is tissue or organoid, tissue_ontology_term_id must be a valid UBERON term. "
          "If organism is NCBITaxon:6239, it can be a valid UBERON term or a valid WBbt term. "
          "If organism is NCBITaxon:7955, it can be a valid UBERON term or a valid ZFA term. "
          "If organism is NCBITaxon:7227, it can be a valid UBERON term or a valid FBbt term. "
          "When tissue_type is primary cell culture, tissue_ontology_term_id must follow the validation rules "
          "for cell_type_ontology_term_id."