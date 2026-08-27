"""Tests for graph_db2.cyto_elements element bookkeeping and labelling."""

from __future__ import annotations

from graph_db2.cyto_elements import (
    drop_node,
    edge_element,
    fan_summary,
    group_element,
    label_for,
    member_options,
    merge_elements,
    node_element,
    not_yet_drawn,
    properties_of,
)
from graph_db2.models import LatticeNode

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    TISSUE,
    build_objects,
    clean_graph_state,
    fake_gatherer,
    raw_matrix_file,
    sequence_file,
)

A = "/tissues/aaa/"
B = "/raw_matrix_files/bbb/"
C = "/sequence_files/ccc/"


# --------------------------------------------------------------------------
# label_for
# --------------------------------------------------------------------------


def test_label_uses_alias_local_part() -> None:
    """Aliases are 'lab-name:local-id'; the local id is the informative half."""
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "aliases": ["some-lab:my_thing"]}
    assert label_for(LatticeNode(TISSUE)) == "my_thing"


def test_label_falls_back_to_accession() -> None:
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "accession": "LATTS000001"}
    assert label_for(LatticeNode(TISSUE)) == "LATTS000001"


def test_label_survives_empty_alias_list() -> None:
    """LatticeNode.alias does aliases[0] after an `is not None` check, which
    raises IndexError on []. label_for reads the cache directly to avoid it."""
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "aliases": []}
    assert label_for(LatticeNode(TISSUE)) == LatticeNode(TISSUE).uuid[:8]


def test_label_of_uncached_node_is_a_uuid_stub() -> None:
    assert label_for(LatticeNode(TISSUE)) == LatticeNode(TISSUE).uuid[:8]


def test_label_for_does_not_fetch_uncached_nodes() -> None:
    """A surprise API call from a property access would make rendering
    unpredictably slow; label_for must stay offline."""
    label_for(LatticeNode(TISSUE))
    assert TISSUE not in LatticeNode._cache


# --------------------------------------------------------------------------
# node_element / group_element
# --------------------------------------------------------------------------


def test_node_element_shape() -> None:
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "aliases": ["lab:tis"]}
    data = node_element(LatticeNode(TISSUE))["data"]
    assert data["id"] == TISSUE
    assert data["label"] == "tis"
    assert data["node_type"] == "Tissue"
    assert data["expanded"] is False
    assert data["color"].startswith("#")


def test_group_element_shape() -> None:
    members = [raw_matrix_file(index) for index in range(4)]
    data = group_element(MFS, "RawMatrixFile", members, "db2_test")["data"]
    assert data["is_group"] is True
    assert data["parent_path"] == MFS
    assert data["members"] == members
    assert data["label"] == "RawMatrixFile × 4"
    assert data["expanded"] is False


# --------------------------------------------------------------------------
# merge_elements
# --------------------------------------------------------------------------


def test_merge_adds_nodes_and_edges() -> None:
    merged = merge_elements([], [_node(A), _node(B)], [edge_element(A, B)])
    assert len(merged) == 3


def test_merge_replaces_a_node_with_its_richer_version() -> None:
    """A stub gains its real label once the batch report lands."""
    stub = {"data": {"id": A, "label": "aaa", "node_type": "Tissue"}}
    labelled = {"data": {"id": A, "label": "my_tissue", "node_type": "Tissue"}}
    merged = merge_elements([stub], [labelled], [])
    assert merged[0]["data"]["label"] == "my_tissue"


def test_merge_keeps_the_first_edge_seen() -> None:
    first = edge_element(A, B)
    first["data"]["marker"] = "original"
    merged = merge_elements([first], [], [edge_element(B, A)])
    assert len(merged) == 1
    assert merged[0]["data"]["marker"] == "original"


def test_merge_is_idempotent() -> None:
    nodes, edges = [_node(A), _node(B)], [edge_element(A, B)]
    once = merge_elements([], nodes, edges)
    assert len(merge_elements(once, nodes, edges)) == len(once)


def test_merge_dedupes_reciprocal_edges() -> None:
    """A Library and its FileSet each name the other, so the relationship gets
    discovered twice and must collapse to a single undirected edge."""
    merged = merge_elements([], [], [edge_element(A, B), edge_element(B, A)])
    assert len(merged) == 1


# --------------------------------------------------------------------------
# drop_node
# --------------------------------------------------------------------------


def test_drop_node_removes_node_and_its_edges() -> None:
    elements = merge_elements(
        [], [_node(A), _node(B), _node(C)], [edge_element(A, B), edge_element(B, C)]
    )
    remaining = drop_node(elements, B)
    ids = {element["data"]["id"] for element in remaining}
    assert B not in ids
    assert ids == {A, C}


def test_drop_node_leaves_unrelated_edges() -> None:
    elements = merge_elements([], [_node(A), _node(B), _node(C)], [edge_element(A, B)])
    remaining = drop_node(elements, C)
    assert edge_element(A, B)["data"]["id"] in {
        element["data"]["id"] for element in remaining
    }


def test_drop_node_missing_id_is_a_noop() -> None:
    elements = merge_elements([], [_node(A)], [])
    assert drop_node(elements, "/tissues/nope/") == elements


# --------------------------------------------------------------------------
# not_yet_drawn
# --------------------------------------------------------------------------


def test_not_yet_drawn_filters_present_and_keeps_order() -> None:
    elements = merge_elements([], [_node(B)], [])
    assert not_yet_drawn(elements, [C, B, A]) == [C, A]


def test_not_yet_drawn_empty_when_all_present() -> None:
    elements = merge_elements([], [_node(A), _node(B)], [])
    assert not_yet_drawn(elements, [A, B]) == []


def test_not_yet_drawn_on_empty_canvas_keeps_everything() -> None:
    assert not_yet_drawn([], [A, B]) == [A, B]


# --------------------------------------------------------------------------
# fan_summary
# --------------------------------------------------------------------------


def test_fan_summary_counts_drawn_neighbors() -> None:
    assert fan_summary([_node(MFS), _node(A), _node(B)]) == "2 drawn"


def test_fan_summary_counts_group_members_not_placeholders() -> None:
    """'1 neighbor' for a collapsed fan of 64 reads as a failure - the whole
    reason this function exists."""
    members = [raw_matrix_file(index) for index in range(64)]
    nodes = [_node(MFS), group_element(MFS, "RawMatrixFile", members, "db2_test")]
    assert fan_summary(nodes) == "64 RawMatrixFile grouped"


def test_fan_summary_reports_both_halves() -> None:
    members = [sequence_file(index) for index in range(30)]
    nodes = [
        _node(MFS),
        _node(A),
        _node(B),
        group_element(MFS, "SequenceFile", members, "db2_test"),
    ]
    assert fan_summary(nodes) == "2 drawn, 30 SequenceFile grouped"


def test_fan_summary_for_a_leaf() -> None:
    assert fan_summary([_node(MFS)]) == "no neighbors"


# --------------------------------------------------------------------------
# properties_of
# --------------------------------------------------------------------------


def test_properties_of_drops_noise_keys() -> None:
    LatticeNode(TISSUE).object_json = {
        "@id": TISSUE,
        "@context": "/terms/",
        "@type": ["Tissue"],
        "audit": {"WARNING": []},
        "actions": [{"name": "edit"}],
        "schema_version": 2,
        "uuid": "ccc",
        "status": "current",
    }
    assert set(properties_of(TISSUE)) == {"@id", "status"}


def test_properties_of_drops_empty_values() -> None:
    LatticeNode(TISSUE).object_json = {
        "@id": TISSUE,
        "aliases": [],
        "note": "",
        "extra": None,
        "meta": {},
        "status": "current",
    }
    assert set(properties_of(TISSUE)) == {"@id", "status"}


def test_properties_of_keeps_falsy_scalars() -> None:
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "count": 0, "flag": False}
    assert set(properties_of(TISSUE)) == {"@id", "count", "flag"}


def test_properties_of_is_sorted() -> None:
    LatticeNode(TISSUE).object_json = {"zebra": 1, "alpha": 2, "@id": TISSUE}
    assert list(properties_of(TISSUE)) == ["@id", "alpha", "zebra"]


def test_properties_of_uncached_node_is_empty() -> None:
    assert properties_of(TISSUE) == {}


# --------------------------------------------------------------------------
# member_options
# --------------------------------------------------------------------------


def test_member_options_are_sorted_by_label() -> None:
    objects = build_objects()
    for index in (2, 0, 1):
        path = raw_matrix_file(index)
        LatticeNode(path).object_json = objects[path]

    options = member_options([raw_matrix_file(index) for index in (2, 0, 1)])
    assert [option["label"] for option in options] == sorted(
        option["label"] for option in options
    )
    assert {option["value"] for option in options} == {
        raw_matrix_file(index) for index in range(3)
    }


def test_member_options_without_labels_fall_back_to_stubs() -> None:
    """Documented precondition: call fetch_labels() first or the dropdown is
    unsearchable because every option is a uuid prefix."""
    options = member_options([raw_matrix_file(0)])
    assert options[0]["label"] == LatticeNode(raw_matrix_file(0)).uuid[:8]


def _node(path: str) -> dict:
    return node_element(LatticeNode(path))
