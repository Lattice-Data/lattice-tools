import tempfile
from fixtures.valid_adatas import (
    test_h5ads,
    validator_with_adatas, 
    AnnDataLabelAppender,
    label_writer,
) 

SCHEMA_VERSION = "5.3.0"
SCHEMA_REFERENCE = f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{SCHEMA_VERSION}/schema.md"


def test_validator_schema_version_is_correct(validator_with_adatas):
    validator = validator_with_adatas
    validator.validate_adata()
    assert validator.is_valid is True
    assert validator.schema_version == SCHEMA_VERSION


def test_label_write_version_and_reference_is_correct(label_writer):
    with tempfile.TemporaryDirectory() as temp_dir:
        labels_path = temp_dir + "labels.h5ad"
        write_successful = label_writer.write_labels(labels_path)

        assert label_writer.adata.uns["schema_version"] == SCHEMA_VERSION
        assert label_writer.adata.uns["schema_reference"] == SCHEMA_REFERENCE
        assert write_successful
        assert not label_writer.errors


def test_label_write_general_pass(label_writer):
    with tempfile.TemporaryDirectory() as temp_dir:
        labels_path = temp_dir + "labels.h5ad"
        write_successful = label_writer.write_labels(labels_path)

        assert write_successful
        assert not label_writer.errors
