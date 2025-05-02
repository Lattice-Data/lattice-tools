import anndata as ad
import multiprocessing
import scanpy as sc
from cellxgene_schema.utils import read_h5ad
from cellxgene_schema.matrix_utils import compute_column_sums
from pathlib import Path
from cellxgene_mods import revise_cxg


LOCATION = Path("/Users/brianmott/Downloads/visium_fif")
files = [file.name for file in LOCATION.iterdir() if ".h5ad" in file.name and "_revised.h5ad" not in file.name]

def assign_feature_is_filtered(row):
    # all 0 X and non-zero raw X need to be fif True
    if row['X_column_sums'] == 0 and row['raw_X_column_sums'] != 0:
        value = True
    # all 0s for both matrices can be True or False, assume original value is correct
    elif (
        row['X_column_sums'] == 0 and
        row['raw_X_column_sums'] == 0
    ):
        value = row['feature_is_filtered']
    # other values not allowed
    else:
        value = False
    return value


def worker(file):
    adata = sc.read_h5ad(LOCATION / file)
    adata = revise_cxg(adata)
    organism = adata.obs["organism_ontology_term_id"].unique()[0]
    adata.uns["organism_ontology_term_id"] = organism
    adata.obs.drop(columns="organism_ontology_term_id", inplace=True)

    # correct fif error
    adata.var['X_column_sums'] = compute_column_sums(adata.X)
    adata.var['raw_X_column_sums'] = compute_column_sums(adata.raw.X)
    adata.var['feature_is_filtered'] = adata.var.apply(assign_feature_is_filtered, axis=1)
    adata.var.drop(columns=['X_column_sums', 'raw_X_column_sums'], inplace=True)

    print(f"Saving {file} to {LOCATION}")
    adata.write(filename=LOCATION / file)
    

if __name__ == "__main__":
    with multiprocessing.Pool() as pool:
        pool.map(worker, files)
