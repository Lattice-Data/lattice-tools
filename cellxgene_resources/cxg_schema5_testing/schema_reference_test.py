import tempfile
from fixtures.valid_adatas import validator_with_all_adatas, label_writer

SCHEMA_VERSION = "5.1.0"
SCHEMA_REFERENCE = f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{SCHEMA_VERSION}/schema.md"


def test_validator_schema_version_is_correct(validator_with_all_adatas):
    validator = validator_with_all_adatas
    validator.validate_adata()
    assert validator.is_valid is True
    assert validator.schema_version == SCHEMA_VERSION


def test_label_write_version_and_reference_is_correct(label_writer):
    with tempfile.TemporaryDirectory() as temp_dir:
        labels_path = temp_dir + "labels.h5ad"
        label_writer.write_labels(labels_path)

        assert label_writer.adata.uns["schema_version"] == SCHEMA_VERSION
        assert label_writer.adata.uns["schema_reference"] == SCHEMA_REFERENCE
        assert label_writer.was_writing_successful
        assert not label_writer.errors
