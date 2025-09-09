"""
QA testing for PRs: https://github.com/chanzuckerberg/single-cell-curation/pull/1430,
https://github.com/chanzuckerberg/single-cell-curation/pull/1443

→ Note that since we don't have the cellosaurus changes in yet, the "cell line" requirements for tissue_ontology_term_id are not yet implemented.

Should pass:
(Y) tissue_type == "cell line" + development_stage_ontology_term_id == "na" + tissue_ontology_term_id == cellosaurus term
(Y) tissue_type == "primary cell culture" + valid tissue_ontology_term_id CL term
(Y) tissue_type == "organoid" + tissue_ontology_term_id == descendant of "UBERON:0001062"
(Y) tissue_type == "tissue" + tissue_ontology_term_id == descendant of "UBERON:0001062"

Should not pass:
(Y) tissue_type == "cell line" + development_stage_ontology_term_id == "na" + tissue_ontology_term_id != cellosaurus term
(Y) tissue_type == "cell culture"
(Y) tissue_type == "organoid" + tissue_ontology_term_id == "UBERON:0000922" for embryo
(Y) tissue_type == "organoid" + tissue_ontology_term_id != descendant of "UBERON:0001062"
(Y) tissue_type == "tissue" + tissue_ontology_term_id != descendant of "UBERON:0001062"
(Y) tissue_type == "primary cell culture" + invalid tissue_ontology_term_id term
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
          assert self.validator.errors == []


     def test_tissue_type_cell_line_valid(self):

          # tissue_type == cell line + development_stage_ontology_term_id == 'na' + tissue_ontology_term_id == cellosaurus term

          self.validator.adata.obs['tissue_type'] = 'cell line'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['development_stage_ontology_term_id'] = 'na'
          self.validator.adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'
          self.validator.validate_adata()
          assert self.validator.is_valid
          assert self.validator.errors == []


     def test_tissue_type_primarycellculture_valid(self):

          # tissue_type == "primary cell culture" + valid tissue_ontology_term_id CL term

          self.validator.adata.obs['tissue_type'] = 'primary cell culture'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          if self.validator.adata.uns['organism_ontology_term_id'] in MODEL_ORG_CL_TERMS.keys():
               self.validator.adata.obs['tissue_ontology_term_id'] = MODEL_ORG_CL_TERMS[self.validator.adata.uns['organism_ontology_term_id']]
          else:
               self.validator.adata.obs['tissue_ontology_term_id'] = 'CL:0000738'
          self.validator.validate_adata()
          assert self.validator.is_valid
          assert self.validator.errors == []


     @pytest.mark.parametrize('types', ['organoid','tissue'])
     def test_tissue_type_organoid_tissue_valid(self,types):

          # tissue_type == "organoid" OR tissue_type == "tissue" + tissue_ontology_term_id == descendant of "UBERON:0001062"

          self.validator.adata.obs['tissue_type'] = types
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['tissue_ontology_term_id'] = 'UBERON:0002048'
          self.validator.validate_adata()
          assert self.validator.is_valid
          assert self.validator.errors == []


     def test_tissue_type_cell_line_invalid(self):

          # tissue_type == cell line + development_stage_ontology_term_id == 'na' + tissue_ontology_term_id != cellosaurus term

          self.validator.adata.obs['tissue_type'] = 'cell line'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['development_stage_ontology_term_id'] = 'na'
          self.validator.adata.obs['tissue_ontology_term_id'] = 'UBERON:0002048'
          self.validator.validate_adata()
          assert not self.validator.is_valid
          assert (
               "ERROR: 'UBERON:0002048' in 'tissue_ontology_term_id' is not a valid ontology term id of 'CVCL'. "
               "When 'tissue_type' is 'cell line', 'tissue_ontology_term_id' must be a valid CVCL term."
               ) in self.validator.errors

     def test_cell_culture_invalid(self):

          # tissue_type == "cell culture"

          self.validator.adata.obs['tissue_type'] = 'cell culture'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          if self.validator.adata.uns['organism_ontology_term_id'] in MODEL_ORG_CL_TERMS.keys():
               self.validator.adata.obs['tissue_ontology_term_id'] = MODEL_ORG_CL_TERMS[self.validator.adata.uns['organism_ontology_term_id']]
          else:
               self.validator.adata.obs['tissue_ontology_term_id'] = 'CL:0000738'
          self.validator.validate_adata()
          assert not self.validator.is_valid
          assert (
               "ERROR: Column 'tissue_type' in dataframe 'obs' contains invalid values '['cell culture']'. "
               "Values must be one of ['primary cell culture', 'organoid', 'tissue', 'cell line']"
               ) in self.validator.errors


     def test_tissue_type_organoid_embryo_invalid(self):

          # tissue_type == "organoid" + tissue_ontology_term_id == "UBERON:0000922" for embryo

          self.validator.adata.obs['tissue_type'] = 'organoid'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['tissue_ontology_term_id'] = 'UBERON:0000922'
          self.validator.validate_adata()
          assert not self.validator.is_valid
          assert (
               "ERROR: 'UBERON:0000922' in 'tissue_ontology_term_id' is not allowed. "
               "When 'tissue_type' is 'organoid', 'tissue_ontology_term_id' must be a valid UBERON, ZFA, FBbt, or WBbt term."
               ) in self.validator.errors


     @pytest.mark.parametrize('types', ['organoid','tissue'])
     def test_tissue_type_organoid_tissue_invalid(self,types):

          # tissue_type == "organoid" OR tissue_type == "tissue" + tissue_ontology_term_id != descendant of "UBERON:0001062"

          self.validator.adata.obs['tissue_type'] = types
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['tissue_ontology_term_id'] = 'UBERON:0000104'
          self.validator.validate_adata()
          assert not self.validator.is_valid
          assert (
               "ERROR: 'UBERON:0000104' in 'tissue_ontology_term_id' is not an allowed term id. "
               f"When 'tissue_type' is '{types}', 'tissue_ontology_term_id' must be a valid UBERON, ZFA, FBbt, or WBbt term."
               ) in self.validator.errors


     def test_tissue_type_primarycellculture_invalid(self):

          # tissue_type == "primary cell culture" + invalid tissue_ontology_term_id CL term

          self.validator.adata.obs['tissue_type'] = 'primary cell culture'
          self.validator.adata.obs['tissue_type'] = self.validator.adata.obs['tissue_type'].astype('category')
          self.validator.adata.obs['tissue_ontology_term_id'] = 'UBERON:8600016'
          self.validator.validate_adata()
          assert not self.validator.is_valid
          for error in self.validator.errors:
               cleaned_error = error.replace("'", "").replace("MUST", "must")
               assert cleaned_error.endswith(
                    "When tissue_type is primary cell culture, tissue_ontology_term_id must follow the validation rules for cell_type_ontology_term_id.")
