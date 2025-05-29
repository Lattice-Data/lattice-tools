"""
QA testing for this issue:
PR for this issue:

Check add-labels
- check uns properties

Should pass:
- schema_reference == "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md" in uns
- schema_version == "6.0.0" in uns


"""


import pytest
from cellxgene_schema.write_labels import AnnDataLabelAppender
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)


@pytest.mark.parametrize("test_h5ads", ALL_H5ADS)
def test_add_labels_uns(validator_with_adatas):

    # add_labels check: uns.schema_reference should be modified to 6.0.0 -> pass
    # add_labels check: uns.schema_version should be modified to 6.0.0 -> pass

    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid
    assert validator.errors == []
    assert "organism" not in validator.adata.uns.keys()

    labeler = AnnDataLabelAppender(validator.adata)
    labeler._add_labels()

    assert "organism" in labeler.adata.uns.keys()
    assert "schema_version" in labeler.adata.uns.keys()
    assert labeler.adata.uns["schema_version"] == "6.0.0"
    assert "schema_reference" in labeler.adata.uns.keys()
    assert labeler.adata.uns["schema_reference"] == "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md"

