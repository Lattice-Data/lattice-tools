"""
QA testing for this issue: https://github.com/chanzuckerberg/single-cell-curation/issues/1074
PR for this issue: https://github.com/chanzuckerberg/single-cell-curation/pull/1357

Testing conditions:
Should not pass
(Y, but no specific error message) absent uns.organism_ontology_term_id
(N) present obs.organism - is currently allowed
(Y, but no specific error message)present obs.organism_ontology_term_id
(N) present uns.organism - is currently allowed
(Y) present uns.organism_ontology_term_id_colors
(Y) present uns.organism_ontology_colors
(Q) uns.organism_ontology_term_id - not on accepted organisms list

Should pass
(Y) uns.organism_ontology_term_id - from list of accepted organisms
"""

import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)

#HUMAN_H5ADS = [file for file in ALL_H5ADS if "human" in file]

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
    "synthetic construct":"NCBITaxon:32630" # spike-ins
}


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
class TestOrganismValidation:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    def test_fixtures_invalid(self):

        # fixtures are starting with organism_ontology_term_id in obs -> fail

        assert "organism_ontology_term_id" in self.validator.adata.obs.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid


    def test_organism_term_not_in_uns(self):

        # uns.organism_ontology_term_id not present -> fail

        '''
        Fails but currently getting only dependencies errors for missing uns.organism_ontology_term_id:

            1)ERROR: Checking values with dependencies failed for adata.obs['sex_ontology_term_id'], likely due to missing column or uns key.
            2)ERROR: Checking values with dependencies failed for adata.obs['self_reported_ethnicity_ontology_term_id'], likely due to missing column or uns key.
            3)ERROR: 'HANCESTRO:0022' in 'self_reported_ethnicity_ontology_term_id' is not a valid value of 'self_reported_ethnicity_ontology_term_id'. When
            'organism_ontology_term_id' is NOT 'NCBITaxon:9606' (Homo sapiens), self_reported_ethnicity_ontology_term_id MUST be 'na'.
            4)ERROR: Checking values with dependencies failed for adata.obs['development_stage_ontology_term_id'], likely due to missing column or uns key.
            5)ERROR: 'HsapDv:0000266' in 'development_stage_ontology_term_id' is not a valid ontology term id of 'UBERON'. When 'organism_ontology_term_id'-specific
            requirements are not defined in the schema definition, 'development_stage_ontology_term_id' MUST be a descendant term id of 'UBERON:0000105'
            excluding 'UBERON:0000071', or unknown.

            -> Error message for catching this could be more straight forward. Like ERROR: 'organism_ontology_term_id' is required in uns."
        '''

        del self.validator.adata.obs["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.uns.keys()
        assert "organism_ontology_term_id" not in self.validator.adata.obs.columns
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: 'organism_ontology_term_id' is required in uns") in self.validator.errors


    def test_organism_in_obs(self):

        # organism in obs -> fail

        '''
        Issue: Organism in obs is allowed.

        '''

        organism_term_id = self.validator.adata.obs["organism_ontology_term_id"].unique()[0]
        del self.validator.adata.obs["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.obs.columns
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term_id
        assert organism_term_id == self.validator.adata.uns["organism_ontology_term_id"]
        organism_name = ACCEPTED_ORGANISMS[organism_term_id]
        self.validator.adata.obs["organism"] = organism_name
        assert "organism" in self.validator.adata.obs.columns
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: '{organism_name}' {error}") in self.validator.errors

    def test_organism_term_in_obs(self):

        # organism_ontology_term_id in obs -> fail

        '''
        Issue: Fails but same issue as test_organism_term_not_in_uns -> error message is not specific
        - this test is also somewhat redundant with first test: test_fixtures_invalid

        '''
        assert "organism_ontology_term_id" in self.validator.adata.obs.columns
        assert "organism_ontology_term_id" not in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (f"ERROR: '{organism_name}' {error}") in self.validator.errors

    def test_organism_term_in_uns(self):

        # uns.organism_ontology_term_id from list of accepted organisms -> pass

        organism_term_id = self.validator.adata.obs["organism_ontology_term_id"].unique()[0]
        assert organism_term_id in ACCEPTED_ORGANISMS.keys()
        del self.validator.adata.obs["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.obs.columns
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term_id
        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    def test_organism_in_uns(self):

        # organism in uns -> fail

        '''
        Issue: organism in uns is allowed.
        '''

        organism_term_id = self.validator.adata.obs["organism_ontology_term_id"].unique()[0]
        del self.validator.adata.obs["organism_ontology_term_id"]
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term_id
        assert organism_term_id == self.validator.adata.uns["organism_ontology_term_id"]
        assert "organism_ontology_term_id" not in self.validator.adata.obs.columns
        assert organism_term_id in ACCEPTED_ORGANISMS.keys()
        organism_name = ACCEPTED_ORGANISMS[organism_term_id]
        self.validator.adata.uns["organism"] = organism_name
        assert "organism" in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        #assert (
        #    f'ERROR: "{organism_name}" {error}'
        #)

    def test_organism_term_id_colors_in_uns(self):

        # present uns.organism_ontology_term_id_colors -> fail

        organism_term_id = self.validator.adata.obs["organism_ontology_term_id"].unique()[0]
        del self.validator.adata.obs["organism_ontology_term_id"]
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term_id
        assert organism_term_id == self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.adata.uns["organism_ontology_term_id_colors"] = ["aqua"]
        assert "organism_ontology_term_id_colors" in self.validator.adata.uns.keys()
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: Colors field uns[organism_ontology_term_id_colors] does not have a corresponding categorical field in obs"
        ) in self.validator.errors


    def test_organism_colors_in_uns(self):

        # present uns.organism_colors -> fail

        organism_term_id = self.validator.adata.obs["organism_ontology_term_id"].unique()[0]
        del self.validator.adata.obs["organism_ontology_term_id"]
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term_id
        assert organism_term_id == self.validator.adata.uns["organism_ontology_term_id"]
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
        Isssue: this error message could also be more specific?
        '''

        del self.validator.adata.obs["organism_ontology_term_id"]
        self.validator.adata.uns["organism_ontology_term_id"] = organism_term
        assert organism_term == self.validator.adata.uns["organism_ontology_term_id"]
        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            f"ERROR: '{organism_term}' in 'organism_ontology_term_id' is not a valid ontology term id of 'NCBITaxon'."
        ) in self.validator.errors