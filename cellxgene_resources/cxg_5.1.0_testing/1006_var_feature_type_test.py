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

# all human gene types, some mouse types
var_index = [
    "ERCC-00016",
    "ENSG00000284662",
    "ENSG00000290825",
    "ENSG00000233750",
    "ENSG00000227232",
    "ENSG00000264293",
    "ENSG00000222623",
    "ENSG00000278267",
    "ENSG00000225931",
    "ENSG00000223972",
    "ENSG00000221083",
    "ENSG00000228463",
    "ENSG00000252956",
    "ENSG00000276674",
    "ENSG00000211451",
    "ENSG00000233999",
    "ENSG00000211693",
    "ENSG00000227182",
    "ENSG00000211687",
    "ENSG00000252830",
    "ENSG00000252404",
    "ENSG00000236597",
    "ENSG00000231202",
    "ENSG00000210049",
    "ENSG00000279493",
    "ENSG00000281133",
    "ENSG00000211593",
    "ENSG00000211592",
    "ENSG00000254017",
    "ENSG00000252350",
    "ENSG00000227191",
    "ENSG00000283367",
    "ENSG00000282431",
    "ENSG00000249446",
    "ENSG00000237111",
    "ENSG00000254893",
    "ENSG00000211459",
    "ENSG00000270123",
    "ENSG00000283571",
    "ENSG00000236824",
    "ENSMUSG00000089699",
    "ENSMUSG00000102851",
    "ENSMUSG00000102693",
    "ENSMUSG00000099473",
    "ENSMUSG00002076544",
    "ENSMUSG00000077244",
    "ENSMUSG00000064842",
    "ENSMUSG00000057715",
    "ENSSASG00005000002",
]

starting_var = pd.DataFrame(
    [
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
        [False],
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
        [False, "DDX11L2_ENSG00000290825", "NCBITaxon:9606", "gene", 1657, "lncRNA"],
        [False, "CICP27", "NCBITaxon:9606", "gene", 3812, "processed_pseudogene"],
        [False, "WASH7P", "NCBITaxon:9606", "gene", 1351, "unprocessed_pseudogene"],
        [False, "RN7SL657P", "NCBITaxon:9606", "gene", 293, "misc_RNA"],
        [False, "RNU6-1100P", "NCBITaxon:9606", "gene", 104, "snRNA"],
        [False, "MIR6859-1", "NCBITaxon:9606", "gene", 68, "miRNA"],
        [False, "ENSG00000225931.3", "NCBITaxon:9606", "gene", 1297, "TEC"],
        [False, "DDX11L1", "NCBITaxon:9606", "gene", 632, "transcribed_unprocessed_pseudogene"],
        [False, "ENSG00000221083.1", "NCBITaxon:9606", "gene", 125, "snoRNA"],
        [False, "ENSG00000228463.10", "NCBITaxon:9606", "gene", 1902, "transcribed_processed_pseudogene"],
        [False, "RNA5SP40", "NCBITaxon:9606", "gene", 110, "rRNA_pseudogene"],
        [False, "IGKV1OR1-1", "NCBITaxon:9606", "gene", 347, "IG_V_pseudogene"],
        [False, "GNRHR2", "NCBITaxon:9606", "gene", 1132, "transcribed_unitary_pseudogene"],
        [False, "IGKV3OR2-268", "NCBITaxon:9606", "gene", 350, "IG_V_gene"],
        [False, "TRGV11", "NCBITaxon:9606", "gene", 463, "TR_V_gene"],
        [False, "VN1R28P", "NCBITaxon:9606", "gene", 896, "unitary_pseudogene"],
        [False, "TRGJ2", "NCBITaxon:9606", "gene", 50, "TR_J_gene"],
        [False, "RNA5SP533", "NCBITaxon:9606", "gene", 110, "rRNA"],
        [False, "ENSG00000252404.1", "NCBITaxon:9606", "gene", 184, "scaRNA"],
        [False, "IGHD7-27", "NCBITaxon:9606", "gene", 11, "IG_D_gene"],
        [False, "TRGVB", "NCBITaxon:9606", "gene", 471, "TR_V_pseudogene"],
        [False, "MT-TF", "NCBITaxon:9606", "gene", 71, "Mt_tRNA"],
        [False, "ENSG00000279493.1", "NCBITaxon:9606", "gene", 513, "artifact"],
        [False, "ENSG00000281133.1", "NCBITaxon:9606", "gene", 96, "pseudogene"],
        [False, "IGKJ5", "NCBITaxon:9606", "gene", 38, "IG_J_gene"],
        [False, "IGKC", "NCBITaxon:9606", "gene", 523, "IG_C_gene"],
        [False, "IGHEP2", "NCBITaxon:9606", "gene", 1256, "IG_C_pseudogene"],
        [False, "RPPH1-3P", "NCBITaxon:9606", "gene", 322, "ribozyme"],
        [False, "TRGC2", "NCBITaxon:9606", "gene", 1013, "TR_C_gene"],
        [False, "ENSG00000283367.1", "NCBITaxon:9606", "gene", 160, "sRNA"],
        [False, "TRBD1", "NCBITaxon:9606", "gene", 12, "TR_D_gene"],
        [False, "TRAJ60", "NCBITaxon:9606", "gene", 57, "TR_J_pseudogene"],
        [False, "IGHJ3P", "NCBITaxon:9606", "gene", 50, "IG_J_pseudogene"],
        [False, "RAP1BL", "NCBITaxon:9606", "gene", 555, "translated_processed_pseudogene"],
        [False, "MT-RNR1", "NCBITaxon:9606", "gene", 954, "Mt_rRNA"],
        [False, "VTRNA2-1", "NCBITaxon:9606", "gene", 128, "vault_RNA"],
        [False, "ENSG00000283571.1", "NCBITaxon:9606", "gene", 306, "IG_pseudogene"],
        [False, "BCYRN1", "NCBITaxon:9606", "gene", 200, "scRNA"],
        [False, "Gm1992", "NCBITaxon:10090", "gene", 250, "lncRNA"],
        [False, "Gm18956", "NCBITaxon:10090", "gene", 480, "processed_pseudogene"],
        [False, "4933401J01Rik", "NCBITaxon:10090", "gene", 1070, "TEC"],
        [False, "Gm18775", "NCBITaxon:10090", "gene", 586, "transcribed_unprocessed_pseudogene"],
        [False, "Gm56304", "NCBITaxon:10090", "gene", 101, "miRNA"],
        [False, "Gm23274", "NCBITaxon:10090", "gene", 130, "snoRNA"],
        [False, "Gm26206", "NCBITaxon:10090", "gene", 110, "snRNA"],
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


def test_var_feature_type_label_write(all_adatas):
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
        assert labels.adata.shape[1] == len(var_index)
        assert_frame_equal(var, var_expected)
        assert not labels.errors
