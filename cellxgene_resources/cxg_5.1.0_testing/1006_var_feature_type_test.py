"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/991
https://github.com/chanzuckerberg/single-cell-curation/pull/1006/
"""

import anndata as ad
import numpy as np
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import tempfile
from scipy import sparse
from fixtures.valid_adatas import (
    FIXTURES_ROOT,
    Validator,
    AnnDataLabelAppender,
)

H5ADS = [
    "valid_human.h5ad",
    "valid_mouse.h5ad",
]

var_index = [
    "ERCC-00016",
    "ENSG00000284662",
    "ENSMUSG00000057715",
    "ENSSASG00005000002",
]

starting_var = pd.DataFrame(
    [
        [False],
        [False],
        [False],
        [False],
    ],
    index=var_index,
    columns=["feature_is_filtered"],
)

# Expected var, this is what the obs above should look like after adding the necessary columns with the validator,
# these columns are defined in the schema
var_expected = pd.DataFrame(
    [
        [False, "ERCC-00016 (spike-in control)", "NCBITaxon:32630", "spike-in", 844, "synthetic"],
        [False, "OR4F16", "NCBITaxon:9606", "gene", 939, "protein_coding"],
        [False, "A830018L16Rik", "NCBITaxon:10090", "gene", 1856, "protein_coding"],
        [False, "ORF1ab_ENSSASG00005000002", "NCBITaxon:2697049", "gene", 21290, "protein_coding"]
    ],
    index=var_index,
    columns=[
        "feature_is_filtered",
        "feature_name",
        "feature_reference",
        "feature_biotype",
        "feature_length",
        "feature_type",
    ],
)
for col in ["feature_name", "feature_reference", "feature_biotype", "feature_type", "feature_length"]:
    var_expected[col] = var_expected[col].astype("category")


@pytest.fixture(params=H5ADS)
def all_adatas(request) -> ad.AnnData:
    adata = ad.read_h5ad(f"{FIXTURES_ROOT}/{request.param}")
    return adata


def test_label_write_version_and_reference_is_correct(all_adatas):
    adata = all_adatas
    zero_matrix = np.zeros(shape=(adata.obs.shape[0], starting_var.shape[0]), dtype=np.float32)
    for row in zero_matrix:
        idx = np.random.randint(0, 3)
        row[idx] = 16
    new_adata = ad.AnnData(
        obs=adata.obs,
        var=starting_var,
        X=sparse.csr_matrix(zero_matrix),
        obsm=adata.obsm,
        uns=adata.uns
    )
    validator = Validator()
    validator.adata = new_adata
    validator.validate_adata()
    assert validator.is_valid
    with tempfile.TemporaryDirectory() as temp_dir:
        labels_path = temp_dir + "labels.h5ad"
        labels = AnnDataLabelAppender(validator)
        labels.write_labels(labels_path)
        var = labels.adata.var

        assert labels.was_writing_successful
        assert labels.adata.shape[0] == adata.shape[0]  # cells remain
        assert labels.adata.shape[1] == 4               # var now test var
        assert_frame_equal(var, var_expected)
        assert not labels.errors
