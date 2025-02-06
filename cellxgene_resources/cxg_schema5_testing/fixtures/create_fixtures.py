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


def mock_matrix(
    cell_num: int, gene_num: int, sparsity: float, float_fill: float | int | None = None
) -> sparse.csr_matrix:
    """
    Create sparse csr matrix in float32 format for test fixtures.
    Defaults to raw count matrix, will fill with specified float value if
    that parameter is provided
    """
    shape = (cell_num, gene_num)
    int_fill = 1

    # int cast rounds toward 0
    step = int(1 / (1 - sparsity))
    if sparsity < 0.5:
        ones = np.ones(shape, dtype=np.float32)
        step = int(1 / sparsity)
        int_fill = 0
        ones[:, ::step] = float_fill if float_fill else int_fill
        matrix = sparse.csr_matrix(ones, dtype=np.float32)
    else:
        matrix = sparse.csr_matrix(shape, dtype=np.float32)
        matrix[:, ::step] = float_fill if float_fill else int_fill

    return matrix


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
    adata = ad.AnnData(
        X=mock_matrix(
            cell_num=obs_df.shape[0],
            gene_num=var_df.shape[0],
            sparsity=0.9,
            float_fill=2.3,
        ),
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
