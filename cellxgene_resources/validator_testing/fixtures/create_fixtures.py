import anndata as ad
import numpy as np
import pandas as pd
from cellxgene_ontology_guide.ontology_parser import OntologyParser 
from dataclasses import dataclass
from enum import Enum
from scipy import sparse
from urllib.parse import quote


GOOGLE_SHEET_FIXTURES_ID = "1mpwDN1GuCBQFqAqlUpuzBMJWHGCXaox_A-5NQA2Qy6k"
GOOGLE_SHEET_URL = f"https://docs.google.com/spreadsheets/d/{GOOGLE_SHEET_FIXTURES_ID}/gviz/tq?tqx=out:csv&sheet="
VAR_META_DF = pd.read_csv(f"{GOOGLE_SHEET_URL}{quote('ensembl_ids')}")

ONTOLOGY_PARSER = OntologyParser()
NEW_ONTOLOGY_ORGANISMS = {
    "Caenorhabditis elegans",
    "Danio rerio",
    "Drosophila melanogaster"
}
STATIC_OBS_COLUMNS = {
    "disease_ontology_term_id": "PATO:0000461",
    "genetic_ancestry_African": float("nan"),
    "genetic_ancestry_East_Asian": float("nan"),
    "genetic_ancestry_European": float("nan"),
    "genetic_ancestry_Indigenous_American": float("nan"),
    "genetic_ancestry_Oceanian": float("nan"),
    "genetic_ancestry_South_Asian": float("nan"),
    "is_primary_data": True,
    "self_reported_ethnicity_ontology_term_id": "na",
}

@dataclass
class OrganismMeta:
    label: str
    
    def __post_init__(self):
        self.term_id = ONTOLOGY_PARSER.get_term_id_by_label(self.label, "NCBITaxon")


class Organism(Enum):
    worm = OrganismMeta("Caenorhabditis elegans")
    zebrafish = OrganismMeta("Danio rerio")
    fly = OrganismMeta("Drosophila melanogaster")
    marmoset = OrganismMeta("Callithrix jacchus")
    gorilla = OrganismMeta("Gorilla gorilla gorilla")
    rhesus = OrganismMeta("Macaca mulatta")
    chimp = OrganismMeta("Pan troglodytes")
    domestic_pig = OrganismMeta("Sus scrofa domesticus")
    lemur = OrganismMeta("Microcebus murinus")
    rabbit = OrganismMeta("Oryctolagus cuniculus")
    rat = OrganismMeta("Rattus norvegicus")
    crab_eating_macaque = OrganismMeta("Macaca fascicularis")

    @property
    def term_id(self):
        return self.value.term_id

    @property
    def label(self):
        return self.value.label

    @property
    def obs_key(self):
        if self.label in NEW_ONTOLOGY_ORGANISMS:
            return f"{self.label.lower().replace(' ', '_')}_meta"
        else:
            return "common_ontology_meta"


def mock_matrix(cell_num: int, gene_num: int, sparsity: float) -> sparse.csr_matrix:
    """
    Create sparse csr matrix in float32 format for test fixtures.
    Copy this matrix and change data array to float values to make a
    normalized matrix with values in the same index locations
    """
    shape = (cell_num, gene_num)
    dense_values = np.prod(shape)
    num_non_zeros = int((1 - sparsity) * dense_values)

    values_per_row = round(num_non_zeros / cell_num)
    index_matrix = np.random.randint(
        low=0,
        high=(gene_num - 1),
        size=(cell_num, values_per_row)
    )
    
    zeros_matrix = np.zeros(shape)
    for index_row, zeros_row in zip(index_matrix, zeros_matrix):
        for index in index_row:
            zeros_row[index] = np.random.randint(1, 1000)
        
    return sparse.csr_matrix(
        zeros_matrix,
        dtype=np.float32
    )


def create_fixture(organism: Organism) -> ad.AnnData:
    """
    Create fixture from metadata on a google sheet and mocked matrices. Returns Anndata object
    
    organism: Organism Enum
    """

    obs_meta_df = pd.read_csv(f"{GOOGLE_SHEET_URL}{quote(organism.obs_key)}")

    obs_df = obs_meta_df.copy()
    obs_df.set_index("obs_index", inplace=True)
    obs_df["organism_ontology_term_id"] = organism.term_id
    for column, fill_value in STATIC_OBS_COLUMNS.items():
        obs_df[column] = fill_value

    assert organism.label in VAR_META_DF["organism"].unique(), f"Ensembl IDs not added to google sheet for {organism.name}, {organism.label}"
    var_df = VAR_META_DF[VAR_META_DF["organism"] == organism.label].copy()
    var_df.set_index("ensembl_id", inplace=True)
    var_df.drop(columns="organism", inplace=True)

    obsm = {"X_umap": np.zeros([obs_df.shape[0], 2])}
    uns = {
        "title": f"{organism.name.capitalize()} Fixture",
        "default_embedding": "X_umap",
    }

    raw_matrix = mock_matrix(
        cell_num=obs_df.shape[0],
        gene_num=var_df.shape[0],
        sparsity=0.9,
    )
    
    # fill copy of raw matrix with random float values in same indices as raw matrix
    normalized_matrix = raw_matrix.copy()
    normalized_matrix.data = np.random.uniform(low=1.0, high=100.0, size=normalized_matrix.data.shape)
    # uniform creates float64 array, no dtype argument, so convert back to float32
    normalized_matrix.data = normalized_matrix.data.astype(np.float32)

    adata = ad.AnnData(
        X=normalized_matrix,
        obs=obs_df,
        var=var_df,
        uns=uns,
        obsm=obsm,
        raw=ad.AnnData(raw_matrix, var=var_df),
    )
    adata.var["feature_is_filtered"] = False

    return adata


def get_all_fixtures(dry_run = False):
    """
    Call to iterate through all Organism enums and save h5ad files to current directory

    dry_run will create adata and print out adata meta without saving to file
    """
    for organism in Organism:
        adata = create_fixture(organism)
        filename = f"valid_{organism.name}.h5ad"
        if not dry_run:
            adata.write(filename=filename, compression="gzip")
            print(f"Created {filename} in local directory")
        else:
            print(f"Dry run generation for {organism.name}")
            print(adata)


if __name__ == "__main__":
    get_all_fixtures()
