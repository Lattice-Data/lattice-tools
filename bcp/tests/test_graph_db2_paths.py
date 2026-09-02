"""Tests for graph_db2.cyto_elements path and colour normalization."""

from __future__ import annotations

import pytest

from graph_db2.cyto_elements import (
    GROUP_SEPARATOR,
    UNMAPPED_COLOR,
    canonical_id,
    color_for,
    css_color,
    edge_element,
    group_id,
    normalize_path,
)
from graph_db2.models import GraphDB2Error, LatticeNode, NodeColor

from tests.graph_db2_helpers import clean_graph_state  # noqa: F401  (autouse)

UUID = "aaaa1111-2222-3333-4444-555555555555"
ALIAS = "test-lab:AVN00_mfs"
CORRECT_TYPE_UUID = f"/matrix_file_sets/{UUID}/"
CORRECT_ALIAS = f"/{ALIAS}/"
CORRECT_UUID = f"/{UUID}/"

# normalize_path() decides nothing about which spelling it was given - only the
# server can - so these are the three shapes it has to pass through untouched
SPELLINGS = ("object_type/uuid", "alias", "uuid")


@pytest.mark.parametrize(
    "raw,expected",
    [
        (f"/matrix_file_sets/{UUID}/", CORRECT_TYPE_UUID),
        (f"matrix_file_sets/{UUID}", CORRECT_TYPE_UUID),
        (f"matrix_file_sets/{UUID}/", CORRECT_TYPE_UUID),
        (f"/matrix_file_sets/{UUID}", CORRECT_TYPE_UUID),
        (f"  /matrix_file_sets/{UUID}/  ", CORRECT_TYPE_UUID),
        (f"\tmatrix_file_sets/{UUID}\n", CORRECT_TYPE_UUID),
        (f"/{ALIAS}/", CORRECT_ALIAS),
        (f"/{ALIAS}", CORRECT_ALIAS),
        (f"{ALIAS}/", CORRECT_ALIAS),
        (f"  /{ALIAS}/  ", CORRECT_ALIAS),
        (f"\t{ALIAS}\n", CORRECT_ALIAS),
        (f"/{UUID}/", CORRECT_UUID),
        (f"/{UUID}", CORRECT_UUID),
        (f"{UUID}/", CORRECT_UUID),
        (f"  /{UUID}/  ", CORRECT_UUID),
        (f"\t{UUID}\n", CORRECT_UUID),
    ],
)
def test_normalize_path_accepts_every_spelling(raw: str, expected: str) -> None:
    assert normalize_path(raw) == expected


def test_normalize_path_is_idempotent() -> None:
    assert (
        normalize_path(normalize_path(f"matrix_file_sets/{UUID}")) == CORRECT_TYPE_UUID
    )
    assert normalize_path(normalize_path(ALIAS)) == CORRECT_ALIAS
    assert normalize_path(normalize_path(UUID)) == CORRECT_UUID


@pytest.mark.parametrize(
    "raw",
    [
        f"{ALIAS}/extra",
        "lab:sample#1",
        "lab:a?b=c",
        "lab:100%pure",
        "lab:my sample",
        "matrix_file_sets",
        "LATSQ000001",
    ],
)
def test_normalize_path_does_not_prejudge_the_spelling(raw: str) -> None:
    """Anything that survives the slash stripping is a question for the server.
    A bare word could be an alias, punctuation is legal inside one, and a
    resolvable path is not knowable without a request - so normalize_path()
    must not reject on shape the way it did when '/type/uuid/' was the only
    accepted form."""
    assert normalize_path(raw) == f"/{raw}/"


@pytest.mark.parametrize("bad", ["", "   ", "\t\n", "/", "//", "///"])
def test_normalize_path_rejects_an_empty_seed(bad: str) -> None:
    """Left alone these become '/', which is a request for the portal root -
    a 200 that resolves to no object at all."""
    with pytest.raises(GraphDB2Error, match="Expected an object path"):
        normalize_path(bad)


@pytest.mark.parametrize(
    "bad",
    [
        f"https://api.data.lattice-data.org/matrix_file_sets/{UUID}/",
        f"http://localhost:8050/matrix_file_sets/{UUID}/",
        f"HTTPS://api.data.lattice-data.org/matrix_file_sets/{UUID}/",
    ],
)
def test_normalize_path_rejects_a_full_url(bad: str) -> None:
    """A URL joined onto the instance server yields a nonsense path; saying so
    up front beats a 404 that names something the user never typed."""
    with pytest.raises(GraphDB2Error, match="not a URL"):
        normalize_path(bad)


@pytest.mark.parametrize("bad", [f"matrix_file_sets//{UUID}", f"a//b//{UUID}"])
def test_normalize_path_rejects_empty_interior_segments(bad: str) -> None:
    with pytest.raises(GraphDB2Error, match="empty path segments"):
        normalize_path(bad)


def test_normalize_path_errors_are_value_errors() -> None:
    """explorer.py catches ValueError around seed handling, so a GraphDB2Error
    that is not one would reach the user as a Dash traceback."""
    assert issubclass(GraphDB2Error, ValueError)


def test_normalize_path_error_names_the_input() -> None:
    with pytest.raises(GraphDB2Error, match=r"got '\\t\\n'"):
        normalize_path("\t\n")


# --------------------------------------------------------------------------
# canonical_id - the guard between a 200 and LatticeNode
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "json_object",
    [
        {"@id": CORRECT_TYPE_UUID},
        {"@id": "/labs/test-lab/"},
        {"@id": "/sequence_files/LATSQ000001/"},
    ],
)
def test_canonical_id_accepts_a_single_object(json_object: dict) -> None:
    assert canonical_id(json_object) == json_object["@id"]


@pytest.mark.parametrize(
    "json_object",
    [
        {},  # a 200 with no '@id' at all
        {"@id": None},
        {"@id": "/"},  # the portal root
        {"@id": "/matrix_file_sets/"},  # a collection
        {"@id": f"/matrix_file_sets/{UUID}/extra/"},
        {"@id": ["/matrix_file_sets/", UUID]},  # not even a string
    ],
)
def test_canonical_id_rejects_anything_lattice_node_cannot_parse(
    json_object: dict,
) -> None:
    """LatticeNode splits an '@id' into exactly two segments, so every one of
    these reaches it as a bare unpacking ValueError or an AttributeError on
    None. Rejecting here is what keeps the message about the seed."""
    assert canonical_id(json_object) is None


def test_canonical_id_output_is_parseable_by_lattice_node() -> None:
    for json_object in ({"@id": CORRECT_TYPE_UUID}, {"@id": "/labs/test-lab/"}):
        assert LatticeNode(canonical_id(json_object)).uuid_path == json_object["@id"]


def test_normalize_path_matches_lattice_node() -> None:
    """The canonical form has to equal LatticeNode.uuid_path, or elements built
    from a raw string get ids that no node answers to."""
    for raw in (f"matrix_file_sets/{UUID}", f"/matrix_file_sets/{UUID}"):
        assert normalize_path(raw) == LatticeNode(raw).uuid_path


def test_css_color_drops_opaque_alpha() -> None:
    """cytoscape.js rejects 8-digit #rrggbbaa and silently drops the mapping;
    vis.js accepted it, which is why every NodeColor value carries one."""
    assert css_color("#cfe2f3ff") == "#cfe2f3"


def test_css_color_converts_partial_alpha_to_rgba() -> None:
    assert css_color("#11223380") == "rgba(17, 34, 51, 0.502)"


@pytest.mark.parametrize("value", ["#cfe2f3", "rgb(1, 2, 3)", "red", ""])
def test_css_color_passes_through_non_8_digit(value: str) -> None:
    assert css_color(value) == value


def test_every_node_color_survives_css_color() -> None:
    for member in NodeColor:
        converted = css_color(member.value)
        assert len(converted) == 7 and converted.startswith("#"), member.name


def test_color_for_concrete_type() -> None:
    assert color_for(LatticeNode("/sequence_files/x/")) == css_color(
        NodeColor.SequenceFile.value
    )


def test_color_for_maps_through_abstract_type() -> None:
    """Tissue has no NodeColor member; ABSTRACT_MAPPING sends it to Biosample."""
    assert color_for(LatticeNode("/tissues/x/")) == css_color(NodeColor.Biosample.value)


def test_color_for_unmapped_type_falls_back() -> None:
    assert color_for(LatticeNode("/labs/some-lab/")) == UNMAPPED_COLOR


def test_edge_element_id_is_order_independent() -> None:
    """A Library and its FileSet each reference the other, so the same
    relationship gets discovered from both ends and must dedupe to one edge."""
    one, other = "/tissues/a/", "/raw_matrix_files/b/"
    assert (
        edge_element(one, other)["data"]["id"] == edge_element(other, one)["data"]["id"]
    )


def test_edge_element_endpoints_are_the_given_paths() -> None:
    element = edge_element("/tissues/b/", "/raw_matrix_files/a/")["data"]
    assert {element["source"], element["target"]} == {
        "/tissues/b/",
        "/raw_matrix_files/a/",
    }


def test_group_id_is_namespaced_under_its_parent() -> None:
    parent = "/matrix_file_sets/x/"
    assert (
        group_id(parent, "RawMatrixFile") == f"{parent}{GROUP_SEPARATOR}RawMatrixFile"
    )


def test_group_ids_differ_per_parent_and_type() -> None:
    assert group_id("/a/1/", "RawMatrixFile") != group_id("/a/2/", "RawMatrixFile")
    assert group_id("/a/1/", "RawMatrixFile") != group_id("/a/1/", "SequenceFile")
