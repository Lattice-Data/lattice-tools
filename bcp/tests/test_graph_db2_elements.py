"""Tests for graph_db2.cyto_elements element bookkeeping and labelling."""

from __future__ import annotations

from math import hypot

import pytest

from graph_db2.cyto_elements import (
    COLUMN_GAP,
    NODE_PITCH,
    RING_RADIUS,
    ROW_PITCH,
    UNPLACED_COLUMN,
    already_drawn,
    anchor_position,
    column_of,
    column_positions,
    drop_node,
    drop_nodes,
    edge_element,
    fan_summary,
    label_for,
    member_options,
    merge_elements,
    node_element,
    not_yet_drawn,
    place_expansion,
    place_placeholders,
    placeholder_element,
    placeholder_id,
    position_of,
    properties_of,
    ring_offsets,
    seed_positions,
    settle,
)
from graph_db2 import cyto_elements
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
# node_element / placeholder_element
# --------------------------------------------------------------------------


def test_node_element_shape() -> None:
    LatticeNode(TISSUE).object_json = {"@id": TISSUE, "aliases": ["lab:tis"]}
    data = node_element(LatticeNode(TISSUE))["data"]
    assert data["id"] == TISSUE
    assert data["label"] == "tis"
    assert data["node_type"] == "Tissue"
    assert data["expanded"] is False
    assert data["color"].startswith("#")


def test_placeholder_element_shape() -> None:
    members = [raw_matrix_file(index) for index in range(4)]
    data = placeholder_element("RawMatrixFile", members, 3)["data"]
    assert data["id"] == placeholder_id("RawMatrixFile")
    assert data["is_group"] is True
    assert data["members"] == members
    assert data["label"] == "RawMatrixFile\n3 out of 4"
    assert data["expanded"] is False


# --------------------------------------------------------------------------
# settle: placeholders and known edges
# --------------------------------------------------------------------------


def _fetched(path: str, neighbors: list[str]) -> None:
    """Register a full profile for `path`, as expand() would leave it"""
    LatticeNode(path).object_json = {"@id": path, "links": neighbors}
    cyto_elements._fully_fetched.add(path)


def _expanded(path: str) -> dict:
    return node_element(LatticeNode(path), expanded=True)


def _placeholder(elements: list[dict], api_name: str) -> dict:
    return next(
        element
        for element in elements
        if element["data"]["id"] == placeholder_id(api_name)
    )


def _members() -> list[str]:
    return [raw_matrix_file(index) for index in range(4)]


def test_every_type_on_the_canvas_gets_a_placeholder() -> None:
    """Not just the ones with something held back - it is also where a drawn
    node is found and cleared."""
    elements = settle([_node(A), _node(B)])
    assert _placeholder(elements, "Tissue")["data"]["label"] == "Tissue\n1 out of 1"
    assert _placeholder(elements, "RawMatrixFile")["data"]["members"] == [B]


def test_an_expanded_node_s_references_are_members() -> None:
    _fetched(MFS, _members())
    held = _placeholder(settle([_expanded(MFS)]), "RawMatrixFile")["data"]
    assert held["members"] == _members()
    assert held["label"] == "RawMatrixFile\n0 out of 4"


def test_an_unexpanded_node_s_references_are_not() -> None:
    """Its profile may be a batch report's partial one, so its neighbor list
    cannot be trusted - and reading it must not cost a request."""
    _fetched(MFS, _members())
    elements = settle([_node(MFS)])
    assert placeholder_id("RawMatrixFile") not in {e["data"]["id"] for e in elements}


def test_placeholder_counts_up_as_members_are_drawn() -> None:
    _fetched(MFS, _members())
    elements = settle(
        merge_elements([_expanded(MFS)], [_node(p) for p in _members()[:2]], [])
    )
    assert _placeholder(elements, "RawMatrixFile")["data"]["label"] == (
        "RawMatrixFile\n2 out of 4"
    )


def test_placeholder_counts_down_as_members_come_off() -> None:
    _fetched(MFS, _members())
    elements = settle(
        merge_elements([_expanded(MFS)], [_node(p) for p in _members()[:2]], [])
    )
    elements = settle(drop_nodes(elements, [_members()[0]]))
    assert _placeholder(elements, "RawMatrixFile")["data"]["label"] == (
        "RawMatrixFile\n1 out of 4"
    )


def test_two_parents_referencing_the_same_fan_share_one_placeholder() -> None:
    _fetched(MFS, _members())
    _fetched(TISSUE, _members())
    elements = settle([_expanded(MFS), _expanded(TISSUE)])
    placeholders = [e for e in elements if e["data"].get("is_group")]
    assert (
        len([p for p in placeholders if p["data"]["node_type"] == "RawMatrixFile"]) == 1
    )
    assert _placeholder(elements, "RawMatrixFile")["data"]["members"] == _members()


def test_a_drawn_member_is_linked_to_every_expanded_node_that_references_it() -> None:
    _fetched(MFS, _members())
    _fetched(TISSUE, _members())
    elements = settle([_expanded(MFS), _expanded(TISSUE), _node(_members()[0])])
    assert {edge["data"]["id"] for edge in elements if "source" in edge["data"]} == {
        edge_element(MFS, _members()[0])["data"]["id"],
        edge_element(TISSUE, _members()[0])["data"]["id"],
    }


def test_settling_twice_changes_nothing() -> None:
    _fetched(MFS, _members())
    once = settle(merge_elements([_expanded(MFS)], [_node(_members()[0])], []))
    assert settle(once) == once


def test_a_placeholder_keeps_its_position() -> None:
    elements = settle([_node(B)])
    moved = [
        {**element, "position": {"x": 5, "y": 7}}
        if element["data"].get("is_group")
        else element
        for element in elements
    ]
    assert _placeholder(settle(moved), "RawMatrixFile")["position"] == {"x": 5, "y": 7}


def test_a_new_placeholder_goes_above_the_top_of_its_column() -> None:
    """Only matters while the layout is held, when nothing else will place it."""
    nodes = [
        {**_node(B), "position": {"x": 300, "y": -40}},
        {**_node(raw_matrix_file(1)), "position": {"x": 300, "y": 60}},
    ]
    held = placeholder_element("RawMatrixFile", [B], 0)
    (placed,) = place_placeholders([held], nodes)
    assert placed["position"] == {"x": 300, "y": -40 - ROW_PITCH}


def test_placeholders_sharing_a_column_stack_rather_than_overlap() -> None:
    nodes = [{**_node(A), "position": {"x": 0, "y": 0}}]
    first = {
        **placeholder_element("Tissue", [A], 0),
        "position": {"x": 0, "y": -ROW_PITCH},
    }
    second = placeholder_element("CellLine", ["/cell_lines/x/"], 1)
    placed = place_placeholders([first, second], nodes)
    assert placed[1]["position"] == {"x": 0, "y": -2 * ROW_PITCH}


def test_a_placeholder_with_no_column_yet_starts_one_on_the_right() -> None:
    nodes = [{**_node(A), "position": {"x": 0, "y": 0}}]
    held = placeholder_element("RawMatrixFile", [B], 1)
    (placed,) = place_placeholders([held], nodes)
    assert placed["position"]["x"] == COLUMN_GAP


def test_nothing_to_place_against_leaves_placeholders_to_the_layout() -> None:
    held = placeholder_element("RawMatrixFile", [B], 1)
    assert place_placeholders([held], [_node(A)]) == [held]


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


def test_merge_keeps_the_position_of_a_replaced_node() -> None:
    """Positions are canvas state the browser owns and a freshly built node has
    none. Dropping it teleports an already-drawn node to the origin the moment
    the layout is held and cannot put it back."""
    placed = {**_node(A), "position": {"x": 40, "y": -12}}
    merged = merge_elements([placed], [_node(A)], [])
    assert merged[0]["position"] == {"x": 40, "y": -12}


def test_merge_invents_no_position() -> None:
    """A node nothing has placed has to stay unplaced, or the layout is handed
    an arrangement it never produced."""
    assert "position" not in merge_elements([], [_node(A)], [])[0]


def test_merge_prefers_a_freshly_seeded_position() -> None:
    """seed_positions() runs before the merge, so a position on the incoming
    node is the deliberate one."""
    placed = {**_node(A), "position": {"x": 40, "y": -12}}
    seeded = {**_node(A), "position": {"x": 0, "y": 300}}
    assert merge_elements([placed], [seeded], [])[0]["position"] == {"x": 0, "y": 300}


def test_merge_still_keeps_expanded_alongside_a_position() -> None:
    previous = {
        "data": {**_node(A)["data"], "expanded": True},
        "position": {"x": 1, "y": 2},
    }
    merged = merge_elements([previous], [_node(A)], [])
    assert merged[0]["data"]["expanded"] is True
    assert merged[0]["position"] == {"x": 1, "y": 2}


# --------------------------------------------------------------------------
# held-layout placement
# --------------------------------------------------------------------------


def test_ring_offsets_places_every_node() -> None:
    assert len(ring_offsets(73)) == 73


def test_ring_offsets_of_nothing_is_empty() -> None:
    assert ring_offsets(0) == []


def test_ring_offsets_keeps_a_small_fan_on_one_ring() -> None:
    radii = {round(hypot(x, y)) for x, y in ring_offsets(6)}
    assert radii == {RING_RADIUS}


def test_ring_offsets_spills_a_wide_fan_onto_further_rings() -> None:
    """A single ring either overlaps or sits so far out that the node the fan
    belongs to is off screen."""
    radii = sorted({round(hypot(x, y)) for x, y in ring_offsets(60)})
    assert len(radii) > 1
    assert radii == [RING_RADIUS * ring for ring in range(1, len(radii) + 1)]


def test_ring_offsets_leaves_room_between_neighbors() -> None:
    offsets = ring_offsets(12)
    gaps = [
        hypot(one[0] - other[0], one[1] - other[1])
        for one, other in zip(offsets, offsets[1:])
    ]
    assert min(gaps) >= NODE_PITCH


def test_seed_positions_places_new_nodes_around_the_anchor() -> None:
    anchor = {"x": 100.0, "y": -50.0}
    placed = seed_positions([_node(B), _node(C)], [_node(A)], anchor)
    for element in placed:
        offset = hypot(
            element["position"]["x"] - anchor["x"],
            element["position"]["y"] - anchor["y"],
        )
        assert round(offset) == RING_RADIUS


def test_seed_positions_leaves_drawn_nodes_alone() -> None:
    """Moving what is already on screen is exactly what holding the layout is
    meant to prevent."""
    drawn = {**_node(A), "position": {"x": 7, "y": 9}}
    placed = seed_positions([_node(A), _node(B)], [drawn], {"x": 0.0, "y": 0.0})
    assert "position" not in placed[0]
    assert placed[1]["position"] != {"x": 7, "y": 9}


def test_seed_positions_gives_every_new_node_a_distinct_spot() -> None:
    """Cytoscape drops a positionless node at the origin, so a whole expansion
    landing in one pile is the failure this replaces."""
    nodes = [_node(path) for path in (A, B, C)]
    spots = {
        tuple(element["position"].values())
        for element in seed_positions(nodes, [], {"x": 0.0, "y": 0.0})
    }
    assert len(spots) == 3


def test_position_of_reads_a_placed_node() -> None:
    placed = {**_node(A), "position": {"x": 3, "y": 4}}
    assert position_of([placed], A) == {"x": 3, "y": 4}


def test_position_of_an_unplaced_or_absent_node_is_none() -> None:
    assert position_of([_node(A)], A) is None
    assert position_of([], A) is None


def test_anchor_prefers_the_live_tap_position() -> None:
    """`elements` only catches up on add, remove and drag, so after a layout
    run it is an arrangement behind; tapNode is what the user just clicked."""
    stale = {**_node(A), "position": {"x": 0, "y": 0}}
    tap_node = {"data": {"id": A}, "position": {"x": 250, "y": 80}}
    assert anchor_position([stale], tap_node, A) == {"x": 250, "y": 80}


def test_anchor_ignores_a_tap_on_a_different_node() -> None:
    placed = {**_node(A), "position": {"x": 5, "y": 6}}
    tap_node = {"data": {"id": B}, "position": {"x": 250, "y": 80}}
    assert anchor_position([placed], tap_node, A) == {"x": 5, "y": 6}


def test_anchor_falls_back_to_the_origin() -> None:
    assert anchor_position([], None, A) == {"x": 0.0, "y": 0.0}


# --------------------------------------------------------------------------
# columns by type
# --------------------------------------------------------------------------

PIPELINE = [
    "HumanDonor",
    "Tissue",
    "PlateBasedLibrary",
    "SequenceFileSet",
    "SequenceFile",
    "RawMatrixFile",
    "MatrixFileSet",
]


def _typed(node_id: str, node_type: str, label: str = "") -> dict:
    """A node element with only the fields column_positions() reads."""
    return {"data": {"id": node_id, "node_type": node_type, "label": label or node_id}}


def test_columns_run_left_to_right_in_pipeline_order() -> None:
    """Donors, then biosamples, then libraries, then sequence files, then
    their sets, then matrix files, then matrix file sets."""
    assert [column_of(name) for name in PIPELINE] == sorted(
        column_of(name) for name in PIPELINE
    )
    assert len({column_of(name) for name in PIPELINE}) == len(PIPELINE)


@pytest.mark.parametrize(
    "one,other",
    [
        ("HumanDonor", "NonHumanDonor"),
        ("Tissue", "CellLine"),
        ("Tissue", "Organoid"),
        ("PlateBasedLibrary", "DropletBasedLibrary"),
    ],
)
def test_subtypes_share_their_legend_column(one: str, other: str) -> None:
    """A column is a legend bucket, so a Tissue and a CellLine line up the
    same way they share a colour - listing every subtype by hand would go
    stale the first time a schema is added."""
    assert column_of(one) == column_of(other)


def test_an_unmapped_type_gets_its_own_column_on_the_end() -> None:
    """Better a visible extra column than silently stacked on the donors."""
    assert column_of("Lab") == UNPLACED_COLUMN
    assert column_of(None) == UNPLACED_COLUMN


def test_column_positions_places_every_node_and_no_edge() -> None:
    elements = merge_elements(
        [], [_typed(A, "Tissue"), _typed(B, "RawMatrixFile")], [edge_element(A, B)]
    )
    assert set(column_positions(elements)) == {A, B}


def test_column_positions_on_an_empty_canvas() -> None:
    assert column_positions([]) == {}


def test_column_positions_puts_the_pipeline_in_order() -> None:
    elements = [_typed(name, name) for name in PIPELINE]
    positions = column_positions(elements)
    assert [positions[name]["x"] for name in PIPELINE] == sorted(
        positions[name]["x"] for name in PIPELINE
    )


def test_columns_are_packed() -> None:
    """A graph with no libraries in it should have no empty gutter where they
    would have gone - the two columns present sit next to each other."""
    positions = column_positions([_typed(A, "HumanDonor"), _typed(B, "MatrixFileSet")])
    assert abs(positions[B]["x"] - positions[A]["x"]) == COLUMN_GAP


def test_a_column_is_spread_and_centred() -> None:
    files = [
        _typed(f"/raw_matrix_files/{index}/", "RawMatrixFile") for index in range(4)
    ]
    heights = sorted(position["y"] for position in column_positions(files).values())
    assert [round(one - other) for one, other in zip(heights[1:], heights)] == (
        [ROW_PITCH] * 3
    )
    assert round(sum(heights)) == 0


def test_one_node_in_a_column_sits_on_the_centre_line() -> None:
    assert column_positions([_typed(A, "Tissue")])[A]["y"] == 0


def test_a_placeholder_sits_at_the_top_of_its_type_s_column() -> None:
    files = [_typed(raw_matrix_file(index), "RawMatrixFile") for index in range(3)]
    held = placeholder_element("RawMatrixFile", [raw_matrix_file(0)], 0)
    positions = column_positions([held, *files])
    column = [positions[f["data"]["id"]] for f in files]
    assert positions[held["data"]["id"]]["x"] == column[0]["x"]
    assert positions[held["data"]["id"]]["y"] == min(p["y"] for p in column) - ROW_PITCH


def test_placeholders_do_not_move_the_rows() -> None:
    files = [_typed(raw_matrix_file(index), "RawMatrixFile") for index in range(3)]
    held = placeholder_element("RawMatrixFile", [raw_matrix_file(0)], 0)
    with_it = column_positions([held, *files])
    without = column_positions(files)
    assert {key: with_it[key] for key in without} == without


def test_a_type_with_nothing_drawn_is_a_column_of_its_placeholder() -> None:
    held = placeholder_element("RawMatrixFile", [raw_matrix_file(0)], 1)
    positions = column_positions([held, _typed(MFS, "MatrixFileSet")])
    assert positions[held["data"]["id"]]["x"] < positions[MFS]["x"]


def test_placeholders_sharing_a_column_stack_in_type_order() -> None:
    tissue = placeholder_element("Tissue", [A], 0)
    cell_line = placeholder_element("CellLine", ["/cell_lines/x/"], 1)
    positions = column_positions([tissue, cell_line, _typed(A, "Tissue")])
    # CellLine sorts first, so it is the higher of the two
    assert (
        positions[cell_line["data"]["id"]]["y"] < positions[tissue["data"]["id"]]["y"]
    )
    assert positions[tissue["data"]["id"]]["y"] == positions[A]["y"] - ROW_PITCH


def test_rows_do_not_cross_the_edges_to_the_next_column() -> None:
    """What the barycentre passes are for. Sorted by label alone these two
    edges cross; a column of 30 files against a column of donors in an
    unrelated order crosses 30 times and is unreadable at any zoom."""
    donors = [
        _typed("/human_donors/1/", "HumanDonor", "a_donor"),
        _typed("/human_donors/2/", "HumanDonor", "b_donor"),
    ]
    # labelled so that the alphabet puts them the wrong way round
    tissues = [
        _typed("/tissues/1/", "Tissue", "z_tissue"),
        _typed("/tissues/2/", "Tissue", "a_tissue"),
    ]
    elements = merge_elements(
        [],
        donors + tissues,
        [
            edge_element("/human_donors/1/", "/tissues/1/"),
            edge_element("/human_donors/2/", "/tissues/2/"),
        ],
    )
    positions = column_positions(elements)

    # whichever way round the pair ended up, each donor is level with its own
    # tissue rather than the other one's
    donor_order = (
        positions["/human_donors/1/"]["y"] < positions["/human_donors/2/"]["y"]
    )
    tissue_order = positions["/tissues/1/"]["y"] < positions["/tissues/2/"]["y"]
    assert donor_order == tissue_order


def test_an_isolated_node_does_not_migrate_to_the_top() -> None:
    """With no neighbors to line up with, a node keeps the row it had."""
    tissues = [
        _typed("/tissues/1/", "Tissue", "a_tissue"),
        _typed("/tissues/2/", "Tissue", "b_tissue"),
        _typed("/tissues/3/", "Tissue", "c_tissue"),
    ]
    positions = column_positions(tissues)
    heights = [positions[element["data"]["id"]]["y"] for element in tissues]
    assert heights == sorted(heights)


def test_column_positions_is_deterministic() -> None:
    """The same graph twice has to come out identical, or every callback that
    re-emits the layout reshuffles the canvas for no reason."""
    elements = merge_elements(
        [],
        [_typed(A, "Tissue"), _typed(B, "RawMatrixFile"), _typed(C, "SequenceFile")],
        [edge_element(A, B), edge_element(B, C)],
    )
    assert column_positions(elements) == column_positions(list(reversed(elements)))


def test_place_expansion_hangs_new_nodes_off_the_tapped_node() -> None:
    tap_node = {"data": {"id": A}, "position": {"x": 300.0, "y": 0.0}}
    placed = place_expansion([_node(A), _node(B)], [_node(A)], tap_node, A)
    assert "position" not in placed[0]
    offset = hypot(placed[1]["position"]["x"] - 300.0, placed[1]["position"]["y"])
    assert round(offset) == RING_RADIUS


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
# drop_nodes
# --------------------------------------------------------------------------


def test_drop_nodes_removes_every_named_node_and_its_edges() -> None:
    elements = merge_elements(
        [], [_node(A), _node(B), _node(C)], [edge_element(A, B), edge_element(A, C)]
    )
    remaining = drop_nodes(elements, [B, C])
    assert {element["data"]["id"] for element in remaining} == {A}


def test_drop_nodes_with_nothing_to_drop_returns_the_same_list() -> None:
    elements = merge_elements([], [_node(A), _node(B)], [edge_element(A, B)])
    assert drop_nodes(elements, []) == elements


def test_drop_nodes_does_not_cascade_to_orphans() -> None:
    """Unticking one member of a group must not silently delete a subtree the
    user expanded from it - the orphan stays, visibly disconnected."""
    elements = merge_elements(
        [], [_node(A), _node(B), _node(C)], [edge_element(A, B), edge_element(B, C)]
    )
    remaining = drop_nodes(elements, [B])
    assert {element["data"]["id"] for element in remaining} == {A, C}


# --------------------------------------------------------------------------
# not_yet_drawn / already_drawn
# --------------------------------------------------------------------------


def test_not_yet_drawn_filters_present_and_keeps_order() -> None:
    elements = merge_elements([], [_node(B)], [])
    assert not_yet_drawn(elements, [C, B, A]) == [C, A]


def test_not_yet_drawn_empty_when_all_present() -> None:
    elements = merge_elements([], [_node(A), _node(B)], [])
    assert not_yet_drawn(elements, [A, B]) == []


def test_not_yet_drawn_on_empty_canvas_keeps_everything() -> None:
    assert not_yet_drawn([], [A, B]) == [A, B]


def test_already_drawn_keeps_only_present_paths_in_order() -> None:
    elements = merge_elements([], [_node(C), _node(A)], [])
    assert already_drawn(elements, [A, B, C]) == [A, C]


def test_already_drawn_and_not_yet_drawn_partition_the_paths() -> None:
    """The group picker's tick state is one half and its to-fetch list the
    other, so a member must land in exactly one of them."""
    elements = merge_elements([], [_node(B)], [])
    paths = [A, B, C]
    assert sorted(already_drawn(elements, paths) + not_yet_drawn(elements, paths)) == (
        sorted(paths)
    )


def test_already_drawn_on_empty_canvas_is_empty() -> None:
    assert already_drawn([], [A, B]) == []


# --------------------------------------------------------------------------
# fan_summary
# --------------------------------------------------------------------------


def test_fan_summary_counts_drawn_neighbors() -> None:
    assert fan_summary([_node(MFS), _node(A), _node(B)]) == "2 drawn"


def test_fan_summary_counts_what_was_held_back() -> None:
    """'no neighbors' for a held-back fan of 64 reads as a failure - the whole
    reason this function exists."""
    _fetched(MFS, [raw_matrix_file(index) for index in range(64)])
    assert fan_summary([_node(MFS)]) == "64 RawMatrixFile in its placeholder"


def test_fan_summary_reports_both_halves() -> None:
    _fetched(MFS, [A, B, *(sequence_file(index) for index in range(30))])
    nodes = [_node(MFS), _node(A), _node(B)]
    assert fan_summary(nodes) == "2 drawn, 30 SequenceFile in its placeholder"


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
