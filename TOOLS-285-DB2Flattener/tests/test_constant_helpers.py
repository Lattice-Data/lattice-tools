import pytest
from dataclasses import dataclass
from generate_constants import (
    SchemaIDs,
    camel_to_snake,
    make_plural,
    make_singular,
    snake_to_camel,
)


@dataclass
class SchemaNames:
    api_name: str
    slug: str
    url_prefix: str


NAMES = [
    SchemaNames("AccessKey", "access_key", "access_keys"),
    SchemaNames("CellLine", "cell_line", "cell_lines"),
    SchemaNames("ControlledTerm", "controlled_term", "controlled_terms"),
    SchemaNames("Document", "document", "documents"),
    SchemaNames("DropletBasedLibrary", "droplet_based_library", "droplet_based_libraries"),
    SchemaNames("ExperimentalCondition", "experimental_condition", "experimental_conditions"),
    SchemaNames("GeneticModification", "genetic_modification", "genetic_modifications"),
    SchemaNames("HumanDonor", "human_donor", "human_donors"),
    SchemaNames("Image", "image", "images"),
    SchemaNames("Lab", "lab", "labs"),
    SchemaNames("MatrixFileSet", "matrix_file_set", "matrix_file_sets"),
    SchemaNames("NonHumanDonor", "non_human_donor", "non_human_donors"),
    SchemaNames("Organoid", "organoid", "organoids"),
    SchemaNames("Page", "page", "pages"),
    SchemaNames("PlateBasedLibrary", "plate_based_library", "plate_based_libraries"),
    SchemaNames("PrimaryCellCulture", "primary_cell_culture", "primary_cell_cultures"),
    SchemaNames("ProcessedMatrixFile", "processed_matrix_file", "processed_matrix_files"),
    SchemaNames("RawMatrixFile", "raw_matrix_file", "raw_matrix_files"),
    SchemaNames("SequenceFile", "sequence_file", "sequence_files"),
    SchemaNames("SequenceFileSet", "sequence_file_set", "sequence_file_sets"),
    SchemaNames("Source", "source", "sources"),
    SchemaNames("TabularFile", "tabular_file", "tabular_files"),
    SchemaNames("Tissue", "tissue", "tissues"),
    SchemaNames("Treatment", "treatment", "treatments"),
    SchemaNames("User", "user", "users"),
]


@pytest.mark.parametrize(
    "input,expected",
    [
        (name.slug, name.url_prefix) for name in NAMES
    ]
)
def test_make_plural(input, expected):
    assert make_plural(input) == expected


@pytest.mark.parametrize(
    "input,expected",
    [
        (name.url_prefix, name.slug) for name in NAMES
    ]
)
def test_make_singular(input, expected):
    assert make_singular(input) == expected


@pytest.mark.parametrize(
    "input,expected",
    [
        (name.slug, name.api_name) for name in NAMES
    ]
)
def test_snake_to_camel(input, expected):
    assert snake_to_camel(input) == expected


@pytest.mark.parametrize(
    "input,expected",
    [
        (name.api_name, name.slug) for name in NAMES
    ]
)
def test_camel_to_snake(input, expected):
    assert camel_to_snake(input) == expected


@pytest.mark.parametrize(
    "api_name,slug,url_prefix", 
    [
        (name.api_name, name.slug, name.url_prefix) for name in NAMES
    ]
)
def test_schema_ids(api_name, slug, url_prefix):
    for input in [api_name, slug, url_prefix]:
        schema_id = SchemaIDs(input)
        assert schema_id.input_str == input
        assert schema_id.api_name == api_name
        assert schema_id.slug == slug
        assert schema_id.url_prefix == url_prefix
