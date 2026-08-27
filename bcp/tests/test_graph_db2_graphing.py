"""Tests for graph_db2.graphing, the whole-graph pyvis path.

build_adjacency_dict walks the graph one node at a time, so the cache is
pre-seeded here the way fetch_all_references_from_cache() would leave it. Any
cache miss becomes a real API request, which the last test pins down.
"""

from __future__ import annotations

import pytest

from graph_db2 import graphing
from graph_db2.graphing import (
    build_adjacency_dict,
    create_nodes_for_pyvis,
    find_edges,
    graph_dict,
)
from graph_db2.models import LatticeNode

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    RAW_MATRIX_FILE_COUNT,
    SEQUENCE_FILE_COUNT,
    TISSUE,
    build_objects,
    clean_graph_state,
    fake_gatherer,
    raw_matrix_file,
    sequence_file,
)


@pytest.fixture
def seeded_cache() -> dict[str, dict]:
    """Every reachable object already cached, as after a batch prefetch."""
    objects = build_objects()
    LatticeNode._cache.update(objects)
    return objects


def test_adjacency_dict_reaches_every_node(seeded_cache: dict[str, dict]) -> None:
    build_adjacency_dict(LatticeNode(MFS))
    expected = (
        1 + RAW_MATRIX_FILE_COUNT + RAW_MATRIX_FILE_COUNT * SEQUENCE_FILE_COUNT + 1
    )
    assert len(graph_dict) == expected


def test_adjacency_dict_is_keyed_by_canonical_path(
    seeded_cache: dict[str, dict],
) -> None:
    build_adjacency_dict(LatticeNode(MFS))
    assert all(key.startswith("/") and key.endswith("/") for key in graph_dict)


def test_adjacency_dict_excludes_labs_and_users(seeded_cache: dict[str, dict]) -> None:
    build_adjacency_dict(LatticeNode(MFS))
    assert not any(key.startswith(("/labs/", "/users/")) for key in graph_dict)


def test_adjacency_dict_terminates_on_cycles(seeded_cache: dict[str, dict]) -> None:
    """RawMatrixFile names its MatrixFileSet and vice versa; without the
    already-visited check the recursion would not stop."""
    build_adjacency_dict(LatticeNode(MFS))
    assert MFS in graph_dict


def test_adjacency_dict_from_a_leaf(seeded_cache: dict[str, dict]) -> None:
    build_adjacency_dict(LatticeNode(sequence_file(0)))
    assert set(graph_dict) == {sequence_file(0)}


def test_graph_dict_is_a_module_global_that_accumulates(
    seeded_cache: dict[str, dict],
) -> None:
    """Documented sharp edge: it is not reset between calls, so two unrelated
    walks merge into one graph."""
    build_adjacency_dict(LatticeNode(sequence_file(0)))
    build_adjacency_dict(LatticeNode(sequence_file(1)))
    assert set(graph_dict) == {sequence_file(0), sequence_file(1)}
    assert graphing.graph_dict is graph_dict


def test_find_edges_pairs_every_neighbor(seeded_cache: dict[str, dict]) -> None:
    build_adjacency_dict(LatticeNode(sequence_file(0)))
    assert find_edges(graph_dict) == []


def test_find_edges_covers_both_directions(seeded_cache: dict[str, dict]) -> None:
    """find_edges emits one tuple per (node, neighbor), so a reciprocal
    reference appears from both ends - pyvis collapses them on add_edges."""
    build_adjacency_dict(LatticeNode(MFS))
    edges = find_edges(graph_dict)
    assert (MFS, raw_matrix_file(0)) in edges
    assert (raw_matrix_file(0), MFS) in edges


def test_find_edges_endpoints_are_all_known_nodes(
    seeded_cache: dict[str, dict],
) -> None:
    build_adjacency_dict(LatticeNode(MFS))
    for source, target in find_edges(graph_dict):
        assert source in graph_dict
        assert target in graph_dict


def test_create_nodes_for_pyvis_shape(seeded_cache: dict[str, dict]) -> None:
    build_adjacency_dict(LatticeNode(MFS))
    nodes = {node["n_id"]: node for node in create_nodes_for_pyvis(graph_dict)}

    assert len(nodes) == len(graph_dict)
    seed = nodes[MFS]
    assert seed["label"] == "test-lab:my_matrixfileset"
    assert seed["color"].endswith("ff")  # pyvis accepts #rrggbbaa; cytoscape does not
    assert seed["title"].startswith("Neighbors: <br>")


def test_create_nodes_for_pyvis_labels_fall_back_to_path(
    seeded_cache: dict[str, dict],
) -> None:
    """SequenceFiles in the fixture have no aliases, so the label falls back to
    the uuid path."""
    build_adjacency_dict(LatticeNode(MFS))
    nodes = {node["n_id"]: node for node in create_nodes_for_pyvis(graph_dict)}
    assert nodes[sequence_file(0)]["label"] == sequence_file(0)


def test_neighbor_titles_tolerate_aliasless_neighbors(
    seeded_cache: dict[str, dict],
) -> None:
    """The hover title joins neighbor labels. LatticeNode.alias returns None
    for an aliasless object, and '<br>'.join() raises TypeError on None - so a
    RawMatrixFile whose SequenceFiles have no aliases used to break the whole
    snapshot."""
    build_adjacency_dict(LatticeNode(MFS))
    nodes = {node["n_id"]: node for node in create_nodes_for_pyvis(graph_dict)}

    title = nodes[raw_matrix_file(0)]["title"]
    assert sequence_file(0) in title
    assert "None" not in title


def test_node_label_survives_empty_alias_list() -> None:
    """LatticeNode.alias does aliases[0] after an `is not None` check, which
    raises IndexError on []."""
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "aliases": []}
    assert graphing.node_label(LatticeNode(TISSUE)) == TISSUE


def test_build_adjacency_dict_requires_a_warm_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Skipping fetch_all_references_from_cache makes the walk request one
    object at a time - roughly 40x slower, per graphing.py's docstring."""
    objects = build_objects()
    requests_made: list[str] = []

    def fake_fetch(self: LatticeNode) -> dict:
        requests_made.append(self.uuid_path)
        return objects[self.uuid_path]

    monkeypatch.setattr(LatticeNode, "fetch_object_json", fake_fetch)
    build_adjacency_dict(LatticeNode(MFS))

    assert len(requests_made) == len(graph_dict)
