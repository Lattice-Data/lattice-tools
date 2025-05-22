"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1357

Testing conditions:
Should not pass
(Q) - absent uns.organism_ontology_term_id -> * this doesn't pass but there is no specific error message
(N) - present obs.organism -> * currently allowed
(N) - present obs.organism_ontology_term_id -> * this passes if uns.organism_ontology_term_id is also present, do we want to control for this?
(N) - present uns.organism -> * currently allowed
(Y) - present uns.organism_ontology_term_id_colors
(Y) - present uns.organism_ontology_colors
(Q) - uns.organism_ontology_term_id not on accepted organisms list -> * this doesn't pass but this error message could also be more specific

Should pass
(Y) - uns.organism_ontology_term_id - from list of accepted organisms
"""

import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

ACCEPTED_ORGANISMS = {
    "NCBITaxon:6239": "Caenorhabditis elegans",
    "NCBITaxon:9483": "Callithrix jacchus",
    "NCBITaxon:7955": "Danio rerio",
    "NCBITaxon:7227": "Drosophila melanogaster",
    "NCBITaxon:9595": "Gorilla gorilla gorilla",
    "NCBITaxon:9606": "Homo sapiens",
    "NCBITaxon:9541": "Macaca fascicularis",
    "NCBITaxon:9544": "Macaca mulatta",
    "NCBITaxon:30608": "Microcebus murinus",
    "NCBITaxon:10090": "Mus musculus",
    "NCBITaxon:9986": "Oryctolagus cuniculus",
    "NCBITaxon:9598": "Pan troglodytes",
    "NCBITaxon:10116": "Rattus norvegicus",
    "NCBITaxon:9823": "Sus scrofa",
    "NCBITaxon:1747270": "Macaca fascicularis aureus",
    "NCBITaxon:1654737": "Macaca mulatta brevicaudus",
    "NCBITaxon:947985": "Mus musculus albula",
    "NCBITaxon:230741": "Oryctolagus cuniculus algirus",
    "NCBITaxon:756884": "Pan troglodytes ellioti",
    "NCBITaxon:947987": "Rattus norvegicus albus",
    "NCBITaxon:1170810": "Sus scrofa affinis",
    "NCBITaxon:9825": "Sus scrofa domesticus"
}

EXEMPT_ORGANISMS = {
    "Severe acute respiratory syndrome coronavirus 2":"NCBITaxon:2697049",
    "synthetic construct":"NCBITaxon:32630"
}


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestOrganismValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_organism_term_in_uns_pass(self):

        # fixtures were updated with organism_ontology_term_id in uns

        assert "organism_ontology_term_id" not in self.validator.adata.obs.columns
        assert "organism_ontology_term_id" in self.validator.adata.uns.keys()
        assert self.validator.adata.uns["organism_ontology_term_id"] in ACCEPTED_ORGANISMS.keys()
        self.validator.validate_adata()
        assert self.validator.is_valid


    def test_organism_term_not_in_uns(self):

        # uns.organism_ontology_term_id not present -> fail

        '''
        Question: Fails but currently only getting dependency errors for missing uns.organism_ontology_term_id.
        -> Error message for catching this could be more straight forward. Like ERROR: 'organism_ontology_term_id' is required in uns.

            1)ERROR: Checking values with dependencies failed for adata.obs['sex_ontology_term_id'], likely due to missing column or uns key.
            2)ERROR: Checking values with dependencies failed for adata.obs['self_reported_ethnicity_ontology_term_id'], likely due to missing column or uns key.
            3)ERROR: 'HANCESTRO:0022' in 'self_reported_ethnicity_ontology_term_id' is not a valid value of 'self_reported_ethnicity_ontology_term_id'. When
            'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'.
            4)ERROR: Checking values with dependencies failed for adata.obs['development_stage_ontology_term_id'], likely due to missing column or uns key.
            5)ERROR: 'HsapDv:0000266' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. When 'organism_ontology_term_id'-specific
            requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105'
            excluding 'UBERON:0000071', or unknown.
        '''

        del self.validator.adata.uns["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: 'organism_ontology_term_id' is required in uns") in self.validator.errors


    def test_organism_in_obs(self):

        # organism in obs -> fail

        '''
        Issue: DOES NOT FAIL! Organism in obs is allowed.

        '''

        organism_term_id = self.validator.adata.uns["organism_ontology_term_id"]
        organism_name = ACCEPTED_ORGANISMS[organism_term_id]
        self.validator.adata.obs["organism"] = organism_name
        assert "organism" in self.validator.adata.obs.columns
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: '{organism_name}' {error}") in self.validator.errors

    def test_organism_term_in_obs(self):

        # organism_ontology_term_id in obs -> fail

        '''
        Issue: DOES NOT FAIL if organism_ontology_term_id is in BOTH uns and obs; if term is removed from uns, it fails with unspecific error message like test_organism_term_not_in_uns

        '''
        organism_term_id = self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.adata.obs["organism_ontology_term_id"] = organism_term_id
        assert "organism_ontology_term_id" in self.validator.adata.obs.columns
        del self.validator.adata.uns["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: '{organism_name}' {error}") in self.validator.errors


    def test_organism_in_uns(self):

        # organism in uns -> fail

        '''
        Issue: DOES NOT FAIL! Organism in uns is allowed.
        '''

        organism_term_id = self.validator.adata.uns["organism_ontology_term_id"]
        assert organism_term_id in ACCEPTED_ORGANISMS.keys()
        organism_name = ACCEPTED_ORGANISMS[organism_term_id]
        self.validator.adata.uns["organism"] = organism_name
        assert "organism" in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f'ERROR: "{organism_name}" {error}') in self.validator.errors

    def test_organism_term_id_colors_in_uns(self):

        # present uns.organism_ontology_term_id_colors -> fail

        self.validator.adata.uns["organism_ontology_term_id_colors"] = ["aqua"]
        assert "organism_ontology_term_id_colors" in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Colors field uns[organism_ontology_term_id_colors] does not have a corresponding categorical field in obs"
        ) in self.validator.errors


    def test_organism_colors_in_uns(self):

        # present uns.organism_colors -> fail

        self.validator.adata.uns["organism_colors"] = ["aqua"]
        assert "organism_colors" in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Colors field uns[organism_colors] does not have a corresponding categorical field in obs. Annotate organism_ontology_term_id_colors instead"
        ) in self.validator.errors

    @pytest.mark.parametrize("organism_term", EXEMPT_ORGANISMS.values())
    def test_exempt_organisms(self, organism_term):

        # uns.organism_ontology_term_id - not on accepted organisms (esp COVID) -> fail

        '''
        Question: This error message should probably be more specific?
        '''

        self.validator.adata.uns["organism_ontology_term_id"] = organism_term
        assert self.validator.adata.uns["organism_ontology_term_id"] in EXEMPT_ORGANISMS.values()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{organism_term}' in 'organism_ontology_term_id' is not a valid ontology term id of 'NCBITaxon'."
        ) in self.validator.errors