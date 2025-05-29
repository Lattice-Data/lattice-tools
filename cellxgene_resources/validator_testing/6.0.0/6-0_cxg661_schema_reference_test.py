"""
QA testing for this issue:
PR for this issue:

Check write_labels (add-labels)
- check uns properties

Should pass:
(N) schema_reference == "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md" in uns
(N) schema_version == "6.0.0" in uns

-> currently both are at 5.3.0

"""


import pytest
import tempfile
import shutil
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    read_h5ad,
    validator_with_adatas
)


@pytest.mark.parametrize("test_h5ads", [ALL_H5ADS[0]])
@pytest.mark.parametrize("property", ["organism","schema_version","schema_reference","citation"])
def test_add_labels_uns(validator_with_adatas,property):

    # add_labels check: uns.schema_reference should be modified to 6.0.0 -> pass  # FAILED
    # add_labels check: uns.schema_version should be modified to 6.0.0 -> pass  # FAILED

    TEMP_DIR = tempfile.mkdtemp()
    schema_version = "6.0.0"
    try:
        validator = validator_with_adatas
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []
        assert property not in validator.adata.uns.keys()

        labeler = AnnDataLabelAppender(validator.adata)
        labeler.write_labels(TEMP_DIR + 'valid_adata_with_labels.h5ad')
        labeled_adata = read_h5ad(TEMP_DIR + 'valid_adata_with_labels.h5ad')
        assert property in labeled_adata.uns.keys()
        assert labeled_adata.uns["schema_reference"] == f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{schema_version}/schema.md"
        assert labeled_adata.uns["schema_version"] == schema_version

    finally:
                shutil.rmtree(TEMP_DIR)
