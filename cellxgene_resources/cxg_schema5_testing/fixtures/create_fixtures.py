import anndata as ad
import numpy as np
import pandas as pd
from enum import Enum
from scipy import sparse
from urllib.parse import quote


class Organism(Enum):
    worm = "Caenorhabditis elegans"
    zebrafish = "Danio rerio"
    fly = "Drosophila melanogaster"


GOOGLE_SHEET_FIXTURES_ID = "1mpwDN1GuCBQFqAqlUpuzBMJWHGCXaox_A-5NQA2Qy6k"
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
    sheet_key = organism.value.lower().replace(" ", "_")
    tab_name = f"{sheet_key}_meta"
    url = f"https://docs.google.com/spreadsheets/d/{GOOGLE_SHEET_FIXTURES_ID}/gviz/tq?tqx=out:csv&sheet={quote(tab_name)}"
    meta_df = pd.read_csv(url)

    obs_df = meta_df.loc[:, "obs_index":"tissue_type"].copy()
    obs_df.set_index("obs_index", inplace=True)
    for column, fill_value in STATIC_OBS_COLUMNS.items():
        obs_df[column] = fill_value

    var_df = pd.DataFrame(meta_df.loc[:, "ensembl_id"])
    var_df.set_index("ensembl_id", inplace=True)

    obsm = {"X_umap": np.zeros([obs_df.shape[0], 2])}
    uns = {
        "title": meta_df["title"].dropna().unique()[0],
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


def get_all_fixtures():
    """
    Call to iterate through all Organism enums and save h5ad files to current directory
    """
    for organism in Organism:
        adata = create_fixture(organism)
        adata.write(filename=f"valid_{organism.name}.h5ad", compression="gzip")
