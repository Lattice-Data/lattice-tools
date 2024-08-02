import anndata as ad
from flattener import main
import gc
import lattice
import os
import sys


SCC_REPO_LOC = os.path.expanduser("~/GitClones/CZI/")
sys.path.append(
    os.path.abspath(SCC_REPO_LOC + "single-cell-curation/cellxgene_schemea_cli/")
)

from cellxgene_schema.validate import Validator

mode = "demo"
files = [
    "LATDF190KNY",
    "LATDF519CUZ",
    "LATDF808TGB",
]

if __name__ == "__main__":
    connection = lattice.Connection(mode)
    for file in files:
        print(f"Flattening file {file}...")
        main(file, connection)
        gc.collect()

    adatas = [f for f in os.listdir() if ".h5ad" in f]

    for file in adatas:
        validator = Validator()
        print(f"Validating file {file}...")
        validator.adata = ad.read_h5ad(file)
        validator.validate_adata()
        if validator.is_valid:
            print(f"File {file} IS VALID")
        else:
            print(f"File {file} NOT VALID")
