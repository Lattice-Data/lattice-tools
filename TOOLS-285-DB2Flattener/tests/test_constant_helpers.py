import pytest
from dataclasses import dataclass
from generate_constants import (
    DEFAULT_DEMO_MODE,
    DEFAULT_PROD_MODE,
    SchemaIDs,
    build_parser,
    camel_to_snake,
    make_plural,
    make_singular,
    resolve_modes,
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


@pytest.mark.parametrize(
    "args,expected",
    [
        # no flags runs both defaults
        ([], [DEFAULT_DEMO_MODE, DEFAULT_PROD_MODE]),
        # a bare flag runs only that instance, on its default mode
        (["--demo"], [DEFAULT_DEMO_MODE]),
        (["-d"], [DEFAULT_DEMO_MODE]),
        (["--prod"], [DEFAULT_PROD_MODE]),
        (["-p"], [DEFAULT_PROD_MODE]),
        (["--demo", "--prod"], [DEFAULT_DEMO_MODE, DEFAULT_PROD_MODE]),
        # a flag with a value overrides the default mode
        (["--demo", "db2_sandbox"], ["db2_sandbox"]),
        (["-d", "db2_sandbox"], ["db2_sandbox"]),
        (["--prod", "db2_myprod"], ["db2_myprod"]),
        (["--demo", "db2_a", "--prod", "db2_b"], ["db2_a", "db2_b"]),
        # mixing a bare flag with a valued one
        (["--demo", "db2_sandbox", "--prod"], ["db2_sandbox", DEFAULT_PROD_MODE]),
    ]
)
def test_resolve_modes_from_args(args, expected):
    parsed = build_parser().parse_args(args)
    assert resolve_modes(parsed.demo, parsed.prod) == expected


def test_parser_defaults_are_none_until_resolved():
    """resolve_modes, not argparse, decides the 'run both' fallback"""
    parsed = build_parser().parse_args([])
    assert parsed.demo is None
    assert parsed.prod is None


def test_parser_rejects_unknown_args():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--nope"])
