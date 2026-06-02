# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
import json
import pytest
from pathlib import Path
from fixtures.ontology_updates import get_valid_ontology_mappings 
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    SPATIAL_H5ADS,              # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    MULTISPECIES_H5ADS,         # noqa: F401
    GeneticPerturbationStrategy,# noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_guide_validator,      # noqa: F401
    yield_guide_file_name,      # noqa: F401
)


JSON_PATH = Path("./7.1.0/")

with open(JSON_PATH.absolute() / "ontology_diff_summary.json", "r") as f:
    ontology_updates = json.load(f)

OBS_ONTOLOGY_COLUMNS = [column for column in ontology_updates]


@pytest.mark.parametrize("test_h5ads", NON_SPATIAL_H5ADS)
class TestOntologyUpdates:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    # want all columns except assay, use another test for new spatial assays
    @pytest.mark.parametrize("column", [column for column in OBS_ONTOLOGY_COLUMNS if column != "assay"])
    @pytest.mark.parametrize("index", [num for num in range(4)])
    def test_new_ontologies_pass(self, column, index):


        valid_ontologies = get_valid_ontology_mappings(self.validator.adata.uns["organism_ontology_term_id"])
        changed_terms = ontology_updates.get(column)

        new_terms = changed_terms.get("new")

        if new_terms == []:
            pytest.skip(f"No new terms for {column}")

        try:
            term = new_terms[index]
        except IndexError:
            pytest.skip("No more new terms in new terms list")

        ontology = term.replace("_", ":").split(":")[0]
        if ontology not in getattr(valid_ontologies, column):
            pytest.skip(f"Fixture organism does not contain valid ontology: {ontology}")

        full_column = f"{column}_ontology_term_id"

        obs_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs[full_column] = self.validator.adata.obs[full_column].cat.add_categories(term)
        self.validator.adata.obs.loc[obs_index, full_column] = term

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("column", OBS_ONTOLOGY_COLUMNS)
    @pytest.mark.parametrize("index", [num for num in range(4)])
    def test_deprecated_ontologies_fail(self, column, index):


        valid_ontologies = get_valid_ontology_mappings(self.validator.adata.uns["organism_ontology_term_id"])
        changed_terms = ontology_updates.get(column)

        deprecated_terms = changed_terms.get("deprecated")

        if deprecated_terms == []:
            pytest.skip(f"No deprecated terms for {column}")

        try:
            term = deprecated_terms[index]
        except IndexError:
            pytest.skip("No more deprecated terms in deprecated terms list")

        ontology = term.replace("_", ":").split(":")[0]
        if ontology not in getattr(valid_ontologies, column):
            pytest.skip(f"Fixture organism does not support current ontology: {ontology}")

        full_column = f"{column}_ontology_term_id"

        obs_index = self.validator.adata.obs.index[0]
        self.validator.adata.obs[full_column] = self.validator.adata.obs[full_column].cat.add_categories(term)
        self.validator.adata.obs.loc[obs_index, full_column] = term

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert len(self.validator.errors) > 0


class TestAssayOntologyUpdates:
    @pytest.fixture(autouse=True)
    def setup(self, validator_with_adatas):
        self.validator = validator_with_adatas


    @pytest.mark.parametrize("test_h5ads", ["visium_human_all_spots.h5ad", "visium_human_some_spots.h5ad"])
    def test_visium_hd_pass(self):

        # new visium hd term passes with 6.5 mm visium fixtures
        
        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0920058"

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("test_h5ads", ["visium_v2_11mm_human.h5ad"])
    def test_visium_hd_with_11mm_fails(self):

        # new visium hd term fails with 11 mm visium fixture
        
        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0920058"

        visium_65_mm_errors = [
            (
                "ERROR: The largest dimension of uns['spatial'][library_id]['images']['hires'] "
                "must be 2000 pixels, it has a largest dimension of 4000 pixels."
            ),
            (
                "ERROR: obs['array_col'] must be between 0 and 127, the min and max are 0 and 223. "
                "This must be the value of the column tissue_positions_in_tissue from the "
                "tissue_positions_list.csv or tissue_positions.csv."
            ),
            (
                "ERROR: obs['array_row'] must be between 0 and 77, the min and max are 0 and 127. "
                "This must be the value of the column tissue_positions_in_tissue from the "
                "tissue_positions_list.csv or tissue_positions.csv."
            ),
        ]

        self.validator.validate_adata()
        assert not self.validator.is_valid
        for error in visium_65_mm_errors:
            assert error in self.validator.errors


    @pytest.mark.parametrize("test_h5ads", ["slide_seq_image_human.h5ad", "slide_seq_no_image_human.h5ad"])
    @pytest.mark.parametrize("term", ["EFO:0920002", "EFO:0920003"])
    def test_curio_seek_pass(self, term):

        # new curio seeker terms pass with slide-seq fixtures

        self.validator.adata.obs["assay_ontology_term_id"] = term

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []


    @pytest.mark.parametrize("test_h5ads", ["slide_seq_image_human.h5ad", "slide_seq_no_image_human.h5ad"])
    def test_curio_seek_parent_term_fails(self):

        # parent term of slide-seq and curio-seeker should fail

        self.validator.adata.obs["assay_ontology_term_id"] = "EFO:0920001"

        self.validator.validate_adata()
        assert not self.validator.is_valid
        assert (
            "ERROR: uns['spatial'] is only allowed when obs['assay_ontology_term_id'] "
            "is either a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) or "
            "a descendant of 'EFO:0920001' (bead-based spatial transcriptomics)"
        ) in self.validator.errors
