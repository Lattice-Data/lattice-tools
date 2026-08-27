"""Tests for graph_db2.cyto_elements path and colour normalization."""

from __future__ import annotations

import pytest

from graph_db2.cyto_elements import (
    GROUP_SEPARATOR,
    UNMAPPED_COLOR,
    color_for,
    css_color,
    edge_element,
    group_id,
    normalize_path,
)
from graph_db2.models import LatticeNode, NodeColor

from tests.graph_db2_helpers import clean_graph_state  # noqa: F401  (autouse)

UUID = "aaaa1111-2222-3333-4444-555555555555"
CANONICAL = f"/matrix_file_sets/{UUID}/"


@pytest.mark.parametrize(
    "raw",
    [
        f"/matrix_file_sets/{UUID}/",
        f"matrix_file_sets/{UUID}",
        f"matrix_file_sets/{UUID}/",
        f"/matrix_file_sets/{UUID}",
        f"  /matrix_file_sets/{UUID}/  ",
        f"\tmatrix_file_sets/{UUID}\n",
    ],
)
def test_normalize_path_accepts_every_spelling(raw: str) -> None:
    assert normalize_path(raw) == CANONICAL


def test_normalize_path_is_idempotent() -> None:
    assert normalize_path(normalize_path(f"matrix_file_sets/{UUID}")) == CANONICAL


@pytest.mark.parametrize(
    "bad",
    [
        "",
        "   ",
        "/",
        "//",
        "matrix_file_sets",
        f"/{UUID}/",
        f"a/b/{UUID}",
        f"matrix_file_sets//{UUID}",
    ],
)
def test_normalize_path_rejects_malformed(bad: str) -> None:
    with pytest.raises(ValueError, match="Expected 'object_type/uuid'"):
        normalize_path(bad)


def test_normalize_path_error_names_the_input() -> None:
    with pytest.raises(ValueError, match="got 'matrix_file_sets'"):
        normalize_path("matrix_file_sets")


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
