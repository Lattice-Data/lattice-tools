"""Tests for graph_db2.cyto_elements one-hop expansion and fan grouping."""

from __future__ import annotations

import pytest
import requests

from graph_db2 import cyto_elements
from graph_db2.cyto_elements import (
    already_drawn,
    drop_nodes,
    edge_element,
    expand,
    fetch_labels,
    is_fully_fetched,
    merge_elements,
    neighbors_drawn,
    placeholder_id,
    promote_members,
    settle,
)
from graph_db2.models import GraphDB2Error, LatticeNode

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    RAW_MATRIX_FILE_COUNT,
    SEQUENCE_FILE_COUNT,
    TISSUE,
    fake_gatherer,
    clean_graph_state,
    dangling_edges,
    db2_env,
    edges_of,
    node_ids,
    patched_fetch,
    patched_requests,
    raw_matrix_file,
    sequence_file,
)

BIG_BUDGET = 500
SMALL_BUDGET = 10


def settled(nodes: list[dict], edges: list[dict]) -> list[dict]:
    """An expansion as it lands on an empty canvas, placeholders and all"""
    return settle(merge_elements([], nodes, edges))


def placeholder(elements: list[dict], api_name: str) -> dict:
    return next(
        element["data"]
        for element in elements
        if element["data"]["id"] == placeholder_id(api_name)
    )


# --------------------------------------------------------------------------
# draw budget
# --------------------------------------------------------------------------


def test_draws_whole_fan_under_budget() -> None:
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert len(nodes) == RAW_MATRIX_FILE_COUNT + 1  # + the seed
    assert len(edges) == RAW_MATRIX_FILE_COUNT


def test_holds_back_an_oversized_fan() -> None:
    """Nothing stands in for it in the expansion itself - it is waiting in its
    type's placeholder once the canvas settles."""
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    assert [node["data"]["id"] for node in nodes] == [MFS]
    assert edges == []
    held = placeholder(settled(nodes, edges), "RawMatrixFile")
    assert len(held["members"]) == RAW_MATRIX_FILE_COUNT


def test_zero_budget_never_holds_back() -> None:
    gatherer = fake_gatherer()
    nodes, _ = expand(MFS, gatherer, draw_budget=0)
    assert len(nodes) == RAW_MATRIX_FILE_COUNT + 1


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


def test_placeholder_label_counts_the_held_back_members() -> None:
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    held = placeholder(settled(nodes, edges), "RawMatrixFile")
    assert held["label"] == f"RawMatrixFile × {RAW_MATRIX_FILE_COUNT}"


def test_a_fan_already_on_the_canvas_gets_its_edges() -> None:
    """Reached again from a second parent, every member is already drawn: all
    that is missing is the edges."""
    everything = [raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)]
    nodes, edges = expand(
        MFS, fake_gatherer(), draw_budget=SMALL_BUDGET, on_canvas=everything
    )
    assert {edge["data"]["target"] for edge in edges} == set(everything)
    assert placeholder(settled(nodes, edges), "RawMatrixFile")["label"] == (
        "RawMatrixFile × 0"
    )


def test_a_partly_drawn_fan_links_the_drawn_members_and_holds_the_rest() -> None:
    drawn = [raw_matrix_file(0), raw_matrix_file(1)]
    nodes, edges = expand(
        MFS, fake_gatherer(), draw_budget=SMALL_BUDGET, on_canvas=drawn
    )
    assert set(drawn) <= {edge["data"]["target"] for edge in edges}
    # the drawn two come back among the nodes, so they are on this canvas too
    held = placeholder(settled(nodes, edges), "RawMatrixFile")
    # still every member, so the picker shows the drawn two ticked
    assert len(held["members"]) == RAW_MATRIX_FILE_COUNT
    assert held["label"] == f"RawMatrixFile × {RAW_MATRIX_FILE_COUNT - len(drawn)}"


def test_drawn_members_do_not_count_toward_the_budget() -> None:
    """What is left to draw fits the budget, so it is drawn rather than
    held back."""
    drawn = [raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT - 5)]
    nodes, edges = expand(
        MFS, fake_gatherer(), draw_budget=SMALL_BUDGET, on_canvas=drawn
    )
    assert len(edges) == RAW_MATRIX_FILE_COUNT


def test_mixed_fan_groups_only_the_oversized_type() -> None:
    """A RawMatrixFile fans out to 2 SequenceFiles, 1 Tissue and 1 MatrixFileSet;
    at a budget of 3 the whole fan is over budget but only types above
    FAN_THRESHOLD are held back, so with FAN_THRESHOLD=25 nothing is."""
    gatherer = fake_gatherer()
    nodes, edges = expand(raw_matrix_file(0), gatherer, draw_budget=3)
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


def test_placeholder_members_are_canonical() -> None:
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS.strip("/"), gatherer, draw_budget=SMALL_BUDGET)
    assert placeholder(settled(nodes, edges), "MatrixFileSet")["members"] == [MFS]


def test_empty_seed_raises_before_any_request() -> None:
    """A bare type name is now a question for the server - it could be an alias -
    but an empty seed cannot name anything, and left alone it becomes a request
    for the server root."""
    gatherer = fake_gatherer()
    with pytest.raises(ValueError, match="Expected an object path"):
        expand("   ", gatherer)
    assert gatherer.calls == []


def test_missing_object_raises_before_fetching() -> None:
    """Resolution now happens before the full fetch, so a seed that names
    nothing fails there. The 404 is chained onto a GraphDB2Error that names the
    seed, because 'HTTPError' on its own does not say which path was wrong."""
    gatherer = fake_gatherer()
    with pytest.raises(GraphDB2Error, match="does-not-exist") as raised:
        expand("/matrix_file_sets/does-not-exist/", gatherer)
    assert isinstance(raised.value.__cause__, requests.HTTPError)


def test_object_that_resolves_but_cannot_be_fetched_raises_http_error() -> None:
    """fetch_full() is still the layer that reports a profile the caller's key
    cannot read, so its HTTPError has to keep escaping expand()."""
    gatherer = fake_gatherer()
    LatticeNode(MFS).object_json = {"@id": MFS}  # resolves from the cache

    def forbidden(node: LatticeNode) -> dict:
        raise requests.HTTPError(f"403 Client Error for url: {node.uuid_path}")

    with pytest.raises(requests.HTTPError, match="403"):
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(cyto_elements, "fetch_full", forbidden)
            expand(MFS, gatherer)


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


def test_is_fully_fetched_tracks_full_fetches() -> None:
    gatherer = fake_gatherer()
    assert not is_fully_fetched(MFS)
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    assert is_fully_fetched(MFS)


# --------------------------------------------------------------------------
# reloading the seed
# --------------------------------------------------------------------------


def test_reload_after_a_walk_leaves_every_node_expandable() -> None:
    """The regression: `expanded` used to read the process-wide fetch cache, so
    after a Load emptied the canvas every node clicked on the first pass still
    claimed to be expanded and the explorer refused to redraw its edges."""
    gatherer = fake_gatherer()

    elements: list[dict] = []
    for path in (MFS, raw_matrix_file(0), sequence_file(0)):
        nodes, edges = expand(path, gatherer, draw_budget=BIG_BUDGET)
        elements = merge_elements(elements, nodes, edges)

    # what pressing Load does: same seed, canvas thrown away
    nodes, edges = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    reloaded = merge_elements([], nodes, edges)

    assert neighbors_drawn(reloaded, MFS)
    assert not neighbors_drawn(reloaded, raw_matrix_file(0))

    # and the second-pass expansion still brings its edges with it
    nodes, edges = expand(raw_matrix_file(0), gatherer, draw_budget=BIG_BUDGET)
    reloaded = merge_elements(reloaded, nodes, edges)
    assert edges, "re-expanded node produced no edges"
    assert dangling_edges(reloaded) == []
    assert neighbors_drawn(reloaded, raw_matrix_file(0))


def test_expanded_flag_survives_being_re_emitted_as_a_neighbor() -> None:
    """Expanding a neighbor re-sends the already-expanded node as a stub; if the
    merge dropped its flag the node would offer a second, edge-free expansion."""
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    elements = merge_elements([], nodes, edges)

    nodes, edges = expand(raw_matrix_file(0), gatherer, draw_budget=BIG_BUDGET)
    elements = merge_elements(elements, nodes, edges)

    assert neighbors_drawn(elements, MFS)


def test_neighbors_drawn_is_false_for_a_node_not_on_the_canvas() -> None:
    assert not neighbors_drawn([], MFS)


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
# drawing from a placeholder
# --------------------------------------------------------------------------


def test_promote_members_draws_only_the_chosen() -> None:
    gatherer = fake_gatherer()
    picked = [raw_matrix_file(3), raw_matrix_file(7)]
    nodes = promote_members(picked, gatherer)
    assert {node["data"]["id"] for node in nodes} == set(picked)


def test_promote_members_resolves_labels() -> None:
    nodes = promote_members([raw_matrix_file(3)], fake_gatherer())
    assert "matrix_" in nodes[0]["data"]["label"]


def test_promote_members_batches_one_call_per_type() -> None:
    gatherer = fake_gatherer()
    promote_members([raw_matrix_file(1), raw_matrix_file(2)], gatherer)
    assert len(gatherer.calls) == 1


def test_promote_members_empty_selection_is_a_noop() -> None:
    gatherer = fake_gatherer()
    assert promote_members([], gatherer) == []
    assert gatherer.calls == []


def _element_ids(elements: list[dict]) -> set[str]:
    return {element["data"]["id"] for element in elements}


def _held_back() -> tuple[list[dict], dict]:
    """MFS expanded under a small budget: every RawMatrixFile held back"""
    gatherer = fake_gatherer()
    nodes, edges = expand(MFS, gatherer, draw_budget=SMALL_BUDGET)
    elements = settled(nodes, edges)
    return elements, placeholder(elements, "RawMatrixFile")


def _pick(elements: list[dict], paths: list[str]) -> list[dict]:
    """What the picker does with newly ticked members"""
    return settle(merge_elements(elements, promote_members(paths, fake_gatherer()), []))


def _untick(elements: list[dict], paths: list[str]) -> list[dict]:
    return settle(drop_nodes(elements, paths))


def test_picker_starts_with_nothing_ticked() -> None:
    elements, held = _held_back()
    assert already_drawn(elements, held["members"]) == []


def test_ticking_members_shows_them_as_drawn() -> None:
    elements, held = _held_back()
    picked = [held["members"][3], held["members"][7]]
    elements = _pick(elements, picked)
    assert already_drawn(elements, held["members"]) == picked


def test_a_picked_member_is_linked_to_what_references_it() -> None:
    """A placeholder has no parent to hang an edge off; the expanded node that
    listed the member is where it belongs."""
    elements, held = _held_back()
    elements = _pick(elements, [held["members"][3]])
    assert edges_of(elements) == [edge_element(MFS, held["members"][3])]


def test_a_placeholder_has_no_edges() -> None:
    elements, held = _held_back()
    elements = _pick(elements, held["members"][:2])
    touching = [
        edge
        for edge in edges_of(elements)
        if held["id"] in (edge["data"]["source"], edge["data"]["target"])
    ]
    assert touching == []


def test_unticking_every_member_restores_the_held_back_graph() -> None:
    """The picker applies the difference between its selection and the canvas,
    so tick-then-untick has to land back exactly where it started."""
    before, held = _held_back()
    elements = _pick(before, [held["members"][3], held["members"][7]])
    assert _element_ids(elements) > _element_ids(before)

    elements = _untick(elements, already_drawn(elements, held["members"]))
    assert _element_ids(elements) == _element_ids(before)
    assert dangling_edges(elements) == []


def test_unticking_one_of_two_leaves_the_other_drawn() -> None:
    elements, held = _held_back()
    kept, dropped = held["members"][3], held["members"][7]
    elements = _untick(_pick(elements, [kept, dropped]), [dropped])
    assert already_drawn(elements, held["members"]) == [kept]
    assert dangling_edges(elements) == []


def test_placeholder_survives_unticking_its_members() -> None:
    """The placeholder is not one of its own members, so pruning the fan back to
    nothing must leave it clickable rather than delete it."""
    elements, held = _held_back()
    elements = _untick(_pick(elements, [held["members"][0]]), [held["members"][0]])
    assert held["id"] in _element_ids(elements)


def test_drawing_every_member_keeps_the_placeholder_at_zero() -> None:
    """It is the type's, not the fan's: it stays so the picker can prune back."""
    elements, held = _held_back()
    elements = _pick(elements, held["members"])
    assert already_drawn(elements, held["members"]) == held["members"]
    assert placeholder(elements, "RawMatrixFile")["label"] == "RawMatrixFile × 0"


def test_fetch_labels_skips_already_cached() -> None:
    gatherer = fake_gatherer()
    fetch_labels([raw_matrix_file(0)], gatherer)
    calls_after_first = len(gatherer.calls)
    fetch_labels([raw_matrix_file(0)], gatherer)
    assert len(gatherer.calls) == calls_after_first
