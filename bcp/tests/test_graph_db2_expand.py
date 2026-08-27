"""Tests for graph_db2.cyto_elements one-hop expansion and fan grouping."""

from __future__ import annotations

import pytest
import requests

from graph_db2.cyto_elements import (
    expand,
    explode,
    fetch_labels,
    group_id,
    is_expanded,
    merge_elements,
    promote_members,
)
from graph_db2.models import LatticeNode

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    RAW_MATRIX_FILE_COUNT,
    SEQUENCE_FILE_COUNT,
    TISSUE,
    fake_gatherer,
    clean_graph_state,
    dangling_edges,
    edges_of,
    node_ids,
    patched_fetch,
    raw_matrix_file,
    sequence_file,
)

BIG_BUDGET = 500
SMALL_BUDGET = 10


def groups_in(nodes: list[dict]) -> list[dict]:
    return [node["data"] for node in nodes if node["data"].get("is_group")]


# --------------------------------------------------------------------------
# draw budget
# --------------------------------------------------------------------------


def test_draws_whole_fan_under_budget() -> None:
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert groups_in(nodes) == []
    assert len(nodes) == RAW_MATRIX_FILE_COUNT + 1  # + the seed
    assert len(edges) == RAW_MATRIX_FILE_COUNT


def test_groups_oversized_fan() -> None:
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    groups = groups_in(nodes)
    assert [group["node_type"] for group in groups] == ["RawMatrixFile"]
    assert len(groups[0]["members"]) == RAW_MATRIX_FILE_COUNT
    assert len(nodes) == 2 and len(edges) == 1


def test_zero_budget_never_groups() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=0)
    assert groups_in(nodes) == []


def test_grouping_skips_label_fetches() -> None:
    """A grouped type costs nothing until opened - its api name comes from the
    path, so no report is needed to draw the placeholder."""
    gatherer = fake_gatherer()
    expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    assert gatherer.calls == []


def test_drawing_fetches_labels_in_one_call_per_type() -> None:
    gatherer = fake_gatherer()
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert [call[0] for call in gatherer.calls] == ["RawMatrixFile"]


def test_group_placeholder_id_is_namespaced() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    assert groups_in(nodes)[0]["id"] == group_id(MFS, "RawMatrixFile")


def test_group_label_reports_member_count() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    assert groups_in(nodes)[0]["label"] == f"RawMatrixFile × {RAW_MATRIX_FILE_COUNT}"


def test_mixed_fan_groups_only_the_oversized_type() -> None:
    """A RawMatrixFile fans out to 2 SequenceFiles, 1 Tissue and 1 MatrixFileSet;
    at a budget of 3 the whole fan is over budget but only types above
    FAN_THRESHOLD collapse, so with FAN_THRESHOLD=25 nothing groups."""
    gatherer = fake_gatherer()
    nodes, edges = expand(raw_matrix_file(0), gatherer, draw_budget=3)
    assert groups_in(nodes) == []
    assert len(edges) == SEQUENCE_FILE_COUNT + 2


# --------------------------------------------------------------------------
# path normalization - the regression that blanked the canvas
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        MFS,
        MFS.strip("/"),
        f"{MFS.strip('/')}/",
        f"/{MFS.strip('/')}",
    ],
)
def test_every_seed_spelling_expands_identically(raw: str) -> None:
    nodes, edges = expand(raw, fake_gatherer(), draw_budget=BIG_BUDGET)
    elements = merge_elements([], nodes, edges)
    assert node_ids(elements) == {MFS} | {
        raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)
    }
    assert dangling_edges(elements) == []


def test_seed_without_slashes_produces_no_dangling_edges() -> None:
    """expand() used to build node ids from LatticeNode.uuid_path but edges from
    the raw argument. Cytoscape drops the entire graph when an edge names a
    missing endpoint, so this failed as a blank canvas rather than a lost edge."""
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS.strip("/"), gatherer, draw_budget=BIG_BUDGET)
    assert dangling_edges(merge_elements([], nodes, edges)) == []


def test_group_parent_path_is_canonical() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS.strip("/"), gatherer, draw_budget=SMALL_BUDGET)
    assert groups_in(nodes)[0]["parent_path"] == MFS


def test_malformed_seed_raises_before_any_request() -> None:
    gatherer = fake_gatherer()
    with pytest.raises(ValueError, match="Expected 'object_type/uuid'"):
        expand("matrix_file_sets", gatherer)
    assert gatherer.calls == []


def test_missing_object_raises_http_error() -> None:
    gatherer = fake_gatherer()
    with pytest.raises(requests.HTTPError):
        expand("/matrix_file_sets/does-not-exist/", gatherer)


# --------------------------------------------------------------------------
# invariants that hold at every hop
# --------------------------------------------------------------------------


def test_no_dangling_edges_after_multi_hop_walk() -> None:
    gatherer = fake_gatherer()
    elements: list[dict] = []
    for path in (MFS, raw_matrix_file(0), sequence_file(0), TISSUE):
        nodes, edges = expand(path, gatherer, draw_budget=BIG_BUDGET)
        elements = merge_elements(elements, nodes, edges)
        assert dangling_edges(elements) == [], f"after expanding {path}"


def test_every_drawn_node_has_a_label() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert all(node["data"]["label"] for node in nodes)


def test_seed_is_marked_expanded_and_neighbors_are_not() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    by_id = {node["data"]["id"]: node["data"] for node in nodes}
    assert by_id[MFS]["expanded"] is True
    assert by_id[raw_matrix_file(0)]["expanded"] is False


def test_is_expanded_tracks_full_fetches() -> None:
    gatherer = fake_gatherer()
    assert not is_expanded(MFS)
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert is_expanded(MFS)


def test_re_expanding_does_not_refetch_the_node(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gatherer = fake_gatherer()
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)

    def _boom(_node: LatticeNode) -> dict:
        raise AssertionError("fetch_full called for an already-expanded node")

    monkeypatch.setattr("graph_db2.cyto_elements.fetch_full", _boom)
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)


def test_unconfigured_type_still_yields_nodes_and_edges() -> None:
    """chunk_and_fetch returns [] for a type with no OBJECT_CONFIG entry. The
    node has to survive on a uuid stub, or its edge would dangle."""
    nodes, edges = expand(
        MFS,
        fake_gatherer(unconfigured={"RawMatrixFile"}),
        draw_budget=BIG_BUDGET,
    )
    elements = merge_elements([], nodes, edges)
    assert len(edges) == RAW_MATRIX_FILE_COUNT
    assert dangling_edges(elements) == []


# --------------------------------------------------------------------------
# opening a group
# --------------------------------------------------------------------------


def test_explode_returns_every_member() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    member_nodes, member_edges = explode(groups_in(nodes)[0], gatherer)
    assert len(member_nodes) == RAW_MATRIX_FILE_COUNT
    assert len(member_edges) == RAW_MATRIX_FILE_COUNT


def test_explode_edges_attach_to_the_parent() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    _, member_edges = explode(groups_in(nodes)[0], gatherer)
    for edge in member_edges:
        assert MFS in (edge["data"]["source"], edge["data"]["target"])


def test_explode_resolves_labels() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    member_nodes, _ = explode(groups_in(nodes)[0], gatherer)
    assert all("matrix_" in node["data"]["label"] for node in member_nodes)


def test_promote_members_draws_only_the_chosen() -> None:
    gatherer = fake_gatherer()
    picked = [raw_matrix_file(3), raw_matrix_file(7)]
    nodes, edges = promote_members(picked, MFS, gatherer)
    assert {node["data"]["id"] for node in nodes} == set(picked)
    assert len(edges) == 2


def test_promote_members_batches_one_call_per_type() -> None:
    gatherer = fake_gatherer()
    promote_members([raw_matrix_file(1), raw_matrix_file(2)], MFS, gatherer)
    assert len(gatherer.calls) == 1


def test_promote_members_empty_selection_is_a_noop() -> None:
    gatherer = fake_gatherer()
    nodes, edges = promote_members([], MFS, gatherer)
    assert (nodes, edges) == ([], [])
    assert gatherer.calls == []


def test_fetch_labels_skips_already_cached() -> None:
    gatherer = fake_gatherer()
    fetch_labels([raw_matrix_file(0)], gatherer)
    calls_after_first = len(gatherer.calls)
    fetch_labels([raw_matrix_file(0)], gatherer)
    assert len(gatherer.calls) == calls_after_first
