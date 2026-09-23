"""Tests for graph_db2.explorer presentation helpers.

Only the pure helpers are exercised here - the callbacks are closures inside
build_app() and need a live Dash app to reach.
"""

from __future__ import annotations

from math import hypot

import pytest

from graph_db2.cyto_elements import (
    RING_RADIUS,
    css_color,
    edge_element,
    group_element,
    merge_elements,
    node_element,
)
from graph_db2.explorer import (
    BASE_STYLESHEET,
    COLUMN_LAYOUT,
    DEFAULT_LAYOUT,
    HOLD_LAYOUT,
    KEEP_VIEW,
    LAYOUTS,
    SAMPLE_SEED,
    computed_here,
    detail_panel,
    format_value,
    layout_for,
    legend,
    status_text,
    suggest_layout,
)
from graph_db2.models import LatticeNode, NodeColor

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    TISSUE,
    built_app,
    clean_graph_state,
    click_node,
    TEST_MODE,
    TEST_SERVER,
    db2_env,
    fire_callback,
    node_ids,
    patched_fetch,
    patched_requests,
    pick_layout,
    press_load,
    raw_matrix_file,
)


def _node(path: str) -> dict:
    return node_element(LatticeNode(path))


def _star(spokes: int) -> list[dict]:
    """One hub joined to `spokes` leaves."""
    paths = [raw_matrix_file(index) for index in range(spokes)]
    return merge_elements(
        [],
        [_node(MFS)] + [_node(path) for path in paths],
        [edge_element(MFS, path) for path in paths],
    )


def _chain(length: int) -> list[dict]:
    paths = [raw_matrix_file(index) for index in range(length)]
    return merge_elements(
        [],
        [_node(path) for path in paths],
        [edge_element(paths[i], paths[i + 1]) for i in range(length - 1)],
    )


# --------------------------------------------------------------------------
# suggest_layout
# --------------------------------------------------------------------------


def test_suggest_layout_picks_concentric_for_a_star() -> None:
    """dagre lays a 64-wide fan out as one ~2900px column, which reads as
    nothing having loaded."""
    assert suggest_layout(_star(64)) == "concentric"


def test_suggest_layout_picks_dagre_for_a_chain() -> None:
    assert suggest_layout(_chain(20)) == DEFAULT_LAYOUT


def test_suggest_layout_defaults_for_small_graphs() -> None:
    """Below the size threshold the shape does not matter enough to override."""
    assert suggest_layout(_star(3)) == DEFAULT_LAYOUT


def test_suggest_layout_on_empty_canvas() -> None:
    assert suggest_layout([]) == DEFAULT_LAYOUT


def test_suggest_layout_on_isolated_nodes() -> None:
    nodes = [_node(raw_matrix_file(index)) for index in range(15)]
    assert suggest_layout(nodes) == DEFAULT_LAYOUT


def test_suggest_layout_returns_a_known_layout() -> None:
    for elements in (_star(64), _chain(20), _star(3), []):
        assert suggest_layout(elements) in LAYOUTS


def test_every_layout_fits_on_run() -> None:
    """Without fit the first layout leaves the graph parked off-centre."""
    for name, options in LAYOUTS.items():
        assert options.get("fit") is True, name
        assert options.get("animate") is False, name


# --------------------------------------------------------------------------
# layout_for - the hold-view switch
# --------------------------------------------------------------------------


def arrangement(layout: dict) -> dict:
    """A layout minus the parts that are not the arrangement itself."""
    return {
        key: value
        for key, value in layout.items()
        if key not in ("fit", "nonce", "positions")
    }


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_layout_for_fits_by_default(name: str) -> None:
    """Unticked has to mean exactly the old behaviour."""
    assert arrangement(layout_for(name)) == arrangement(LAYOUTS[name])
    assert layout_for(name, keep_view=False)["fit"] is True


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_layout_for_holding_the_view_turns_fit_off(name: str) -> None:
    """`fit` is the whole mechanism: every element change re-renders the
    component, react-cytoscapejs sees a new layout object and re-runs it, and a
    run with fit on ends in cy.fit() - the zoom-out that loses your place."""
    assert layout_for(name, keep_view=True)["fit"] is False


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_layout_for_changes_nothing_but_fit(name: str) -> None:
    """The layout the user picked still has to be the layout that runs -
    holding the view is not a different arrangement."""
    held = layout_for(name, keep_view=True)
    assert arrangement(held) == arrangement(LAYOUTS[name])


def test_layout_for_does_not_mutate_the_preset() -> None:
    """LAYOUTS is module state shared by every session on the process."""
    layout_for(DEFAULT_LAYOUT, keep_view=True)
    assert LAYOUTS[DEFAULT_LAYOUT]["fit"] is True


# --------------------------------------------------------------------------
# the hold-view control, wired up
# --------------------------------------------------------------------------


def toolbar_control(app, control_id: str):
    for component in app.layout.children[0].children:
        if getattr(component, "id", None) == control_id:
            return component
        # the checklist sits inside a Div that carries its tooltip
        child = getattr(component, "children", None)
        if getattr(child, "id", None) == control_id:
            return child
    raise AssertionError(f"no {control_id} in the toolbar")


def test_hold_view_box_starts_unticked() -> None:
    """Holding the view changes how the canvas behaves, so it has to be opt in:
    an untouched toolbar is the behaviour that shipped before."""
    assert toolbar_control(built_app(MFS), "keep-view").value == []


def test_hold_view_box_offers_exactly_the_keep_value() -> None:
    options = toolbar_control(built_app(MFS), "keep-view").options
    assert [option["value"] for option in options] == [KEEP_VIEW]


def test_hold_view_box_is_explained() -> None:
    """A checkbox labelled 'hold view' says nothing about what it holds."""
    rendered = str(built_app(MFS).layout.children[0])
    assert "zoom" in rendered.lower() and "untick" in rendered.lower()


def test_starting_layout_fits() -> None:
    app = built_app(MFS)
    assert app.layout.children[1].children[0].children.layout["fit"] is True


def test_ticking_the_box_turns_fit_off() -> None:
    assert pick_layout(built_app(MFS), DEFAULT_LAYOUT, [KEEP_VIEW])["fit"] is False


def test_unticking_the_box_turns_fit_back_on() -> None:
    """Unticking is the 'show me the whole graph again' gesture - the callback
    re-emits the layout, which re-runs it and fits."""
    assert pick_layout(built_app(MFS), DEFAULT_LAYOUT, [])["fit"] is True


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_picking_a_layout_while_holding_the_view_keeps_that_layout(name: str) -> None:
    layout = pick_layout(built_app(MFS), name, [KEEP_VIEW])
    assert layout["name"] == LAYOUTS[name]["name"]
    assert layout["fit"] is False


# --------------------------------------------------------------------------
# layout_for - the nonce
# --------------------------------------------------------------------------


def test_layout_for_carries_no_nonce_by_default() -> None:
    assert "nonce" not in layout_for(DEFAULT_LAYOUT)


def test_layout_for_adds_nothing_but_the_nonce() -> None:
    """react-cytoscapejs compares the layout prop key by key, so a re-layout
    the user asked for needs the dict to differ - and differ in nothing that
    changes the arrangement."""
    numbered = layout_for(DEFAULT_LAYOUT, nonce=7)
    assert numbered["nonce"] == 7
    assert {key: value for key, value in numbered.items() if key != "nonce"} == (
        layout_for(DEFAULT_LAYOUT)
    )


# --------------------------------------------------------------------------
# columns by type
# --------------------------------------------------------------------------


def test_the_column_layout_is_a_preset() -> None:
    """Cytoscape has no layout that groups by an attribute, so this one is
    computed here and handed over as positions."""
    assert LAYOUTS[COLUMN_LAYOUT]["name"] == "preset"
    assert computed_here(COLUMN_LAYOUT) is True


@pytest.mark.parametrize("name", [n for n in LAYOUTS if n != COLUMN_LAYOUT])
def test_every_other_layout_is_cytoscape_s_own(name: str) -> None:
    assert computed_here(name) is False


def test_the_column_layout_carries_positions() -> None:
    elements = _star(4)
    positions = layout_for(COLUMN_LAYOUT, elements=elements)["positions"]
    assert set(positions) == node_ids(elements)


def test_the_column_layout_without_a_canvas_is_empty_but_valid() -> None:
    """The dropdown can be changed before anything is loaded."""
    assert layout_for(COLUMN_LAYOUT)["positions"] == {}


@pytest.mark.parametrize("name", [n for n in LAYOUTS if n != COLUMN_LAYOUT])
def test_no_other_layout_carries_positions(name: str) -> None:
    """dagre works the arrangement out for itself; sending positions with it
    would be dead weight on every callback."""
    assert "positions" not in layout_for(name, elements=_star(4))


def test_the_canvas_changing_re_columns_the_graph() -> None:
    """A preset is only ever the positions it was built from, so a node added
    after it was built has no place in it."""
    layout = pick_layout(
        built_app(MFS),
        COLUMN_LAYOUT,
        [],
        [],
        triggered="graph.elements",
        elements=_star(4),
    )
    assert layout["positions"]


def test_dragging_a_node_does_not_re_column_the_graph() -> None:
    """A drag pushes the elements back too, with the same nodes on them; a
    fresh preset would snap the dragged node straight back to its column."""
    elements = _star(4)
    layout = pick_layout(
        built_app(MFS),
        COLUMN_LAYOUT,
        [],
        [],
        triggered="graph.elements",
        elements=elements,
        current=layout_for(COLUMN_LAYOUT, elements=elements),
    )
    assert layout is None


def test_a_node_drawn_after_a_drag_still_re_columns_the_graph() -> None:
    before = _star(3)
    after = _star(4)
    layout = pick_layout(
        built_app(MFS),
        COLUMN_LAYOUT,
        [],
        [],
        triggered="graph.elements",
        elements=after,
        current=layout_for(COLUMN_LAYOUT, elements=before),
    )
    assert set(layout["positions"]) == node_ids(after)


@pytest.mark.parametrize("name", [n for n in LAYOUTS if n != COLUMN_LAYOUT])
def test_the_canvas_changing_leaves_cytoscape_s_own_layouts_alone(name: str) -> None:
    """dagre is re-run by dash-cytoscape on add and remove; re-emitting it
    from here as well would run it twice for every click."""
    assert pick_layout(built_app(MFS), name, [], [], triggered="graph.elements") is None


def test_the_canvas_changing_does_not_re_column_a_held_graph() -> None:
    """Hold Layout outranks it: the whole point is that drawing a node moves
    nothing."""
    layout = pick_layout(
        built_app(MFS),
        COLUMN_LAYOUT,
        [],
        [HOLD_LAYOUT],
        triggered="graph.elements",
        elements=_star(4),
    )
    assert layout is None


def test_load_re_columns_a_held_graph() -> None:
    """Load changes the elements and bumps relayout in one action. Reading
    only the first of the two triggers would drop the run it asked for, and a
    freshly loaded graph would sit in a pile on the origin."""
    layout = pick_layout(
        built_app(MFS),
        COLUMN_LAYOUT,
        [],
        [HOLD_LAYOUT],
        triggered="graph.elements,relayout.data",
        bump=1,
        elements=_star(4),
    )
    assert layout["positions"]


# --------------------------------------------------------------------------
# the hold-layout control, wired up
# --------------------------------------------------------------------------


def auto_refresh(app, keep_layout: list[str]) -> bool:
    """Whether dash-cytoscape re-runs the layout on every add and remove."""
    response = fire_callback(
        app,
        "graph.autoRefreshLayout",
        keep_layout,
        outputs={"id": "graph", "property": "autoRefreshLayout"},
        triggered="keep-layout.value",
    )
    return response["graph"]["autoRefreshLayout"]


def test_hold_layout_box_starts_unticked() -> None:
    assert toolbar_control(built_app(MFS), "keep-layout").value == []


def test_hold_layout_box_offers_exactly_the_hold_value() -> None:
    options = toolbar_control(built_app(MFS), "keep-layout").options
    assert [option["value"] for option in options] == [HOLD_LAYOUT]


def test_hold_layout_box_is_explained() -> None:
    """'Hold Layout' next to 'Hold View' says nothing about which is which."""
    rendered = str(built_app(MFS).layout.children[0]).lower()
    assert "ring" in rendered and "untick" in rendered


def test_the_graph_refreshes_its_layout_until_told_not_to() -> None:
    """Unticked has to mean exactly the old behaviour."""
    app = built_app(MFS)
    assert app.layout.children[1].children[0].children.autoRefreshLayout is True
    assert auto_refresh(app, []) is True


def test_ticking_hold_layout_stops_the_automatic_relayout() -> None:
    """The whole mechanism: dash-cytoscape re-runs the layout on `add remove`,
    and dagre re-flows every node when one arrives. Off, nodes stay put."""
    assert auto_refresh(built_app(MFS), [HOLD_LAYOUT]) is False


def test_ticking_hold_layout_does_not_re_arrange_the_graph() -> None:
    """Ticking it says 'leave what is on screen alone', so it must not be the
    one thing that moves it."""
    layout = pick_layout(
        built_app(MFS),
        DEFAULT_LAYOUT,
        [],
        [HOLD_LAYOUT],
        triggered="keep-layout.value",
    )
    assert layout is None


def test_unticking_hold_layout_re_arranges_the_graph() -> None:
    """Releasing it is the 'tidy this back up' gesture."""
    layout = pick_layout(
        built_app(MFS), DEFAULT_LAYOUT, [], [], triggered="keep-layout.value"
    )
    assert layout["name"] == LAYOUTS[DEFAULT_LAYOUT]["name"]


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_picking_a_layout_while_holding_it_still_runs(name: str) -> None:
    """The dropdown is the escape hatch: a held graph the user wants
    re-arranged once."""
    layout = pick_layout(
        built_app(name and DEFAULT_LAYOUT),
        name,
        [],
        [HOLD_LAYOUT],
        triggered="layout-choice.value",
    )
    assert layout["name"] == LAYOUTS[name]["name"]


def test_a_relayout_bump_runs_the_layout() -> None:
    """Load's way of asking for one run on a graph that has no positions to
    hold yet."""
    app = built_app(MFS)
    layout = pick_layout(
        app, DEFAULT_LAYOUT, [], [HOLD_LAYOUT], triggered="relayout.data", bump=1
    )
    assert layout["name"] == LAYOUTS[DEFAULT_LAYOUT]["name"]


def test_consecutive_layouts_are_never_the_same_dict() -> None:
    """An identical dict is not a new prop, and react-cytoscapejs only runs the
    layout when the prop changes - so a second request for the same layout
    would do nothing at all."""
    app = built_app(MFS)
    first = pick_layout(app, DEFAULT_LAYOUT, [])
    second = pick_layout(app, DEFAULT_LAYOUT, [])
    assert first != second
    assert first["nonce"] != second["nonce"]


def test_the_two_hold_boxes_are_independent() -> None:
    """Holding the layout says nothing about where the viewport goes."""
    layout = pick_layout(
        built_app(MFS),
        DEFAULT_LAYOUT,
        [KEEP_VIEW],
        [HOLD_LAYOUT],
        "layout-choice.value",
    )
    assert layout["fit"] is False
    layout = pick_layout(
        built_app(MFS), DEFAULT_LAYOUT, [], [HOLD_LAYOUT], "layout-choice.value"
    )
    assert layout["fit"] is True


# --------------------------------------------------------------------------
# holding the layout through an expansion
# --------------------------------------------------------------------------


def drawn_canvas(app) -> list[dict]:
    """
    The elements the browser is holding after a Load, as it reports them back.

    dash-cytoscape pushes the canvas up on every add, remove and drag, with a
    position on every node - so this is what the next callback's State sees,
    not the positionless list the server sent down.
    """
    elements = press_load(app, MFS)["graph"]["elements"]
    return [
        element
        if "source" in element["data"]
        else {**element, "position": {"x": 0.0, "y": float(index * 30)}}
        for index, element in enumerate(elements)
    ]


def positions_in(elements: list[dict]) -> dict[str, dict]:
    return {
        element["data"]["id"]: element["position"]
        for element in elements
        if "position" in element
    }


def expand_first_file(app, elements: list[dict], keep_layout: list[str]) -> list[dict]:
    """Tap a drawn RawMatrixFile, with the canvas reporting where it sits."""
    target = raw_matrix_file(0)
    tapped = next(
        element["data"] for element in elements if element["data"]["id"] == target
    )
    response = click_node(
        app,
        tapped,
        elements,
        keep_layout=keep_layout,
        tap_node={"data": {"id": target}, "position": {"x": 500.0, "y": 25.0}},
    )
    return response["graph"]["elements"]


def test_a_held_expansion_leaves_every_drawn_node_where_it_was() -> None:
    """The point of the box. A node re-emitted by the expansion - the one that
    was tapped, and any neighbour already on screen - is rebuilt from the fetch
    cache and arrives without a position; keeping it is what stops the canvas
    jumping."""
    app = built_app("", TEST_MODE)
    before = drawn_canvas(app)
    after = expand_first_file(app, before, [HOLD_LAYOUT])
    kept = positions_in(after)
    assert positions_in(before).items() <= kept.items()


def test_a_held_expansion_rings_its_new_nodes_round_the_tapped_node() -> None:
    """Cytoscape drops a positionless node at the origin, so without this the
    whole expansion lands in one pile in the corner."""
    app = built_app("", TEST_MODE)
    before = drawn_canvas(app)
    after = expand_first_file(app, before, [HOLD_LAYOUT])

    fresh = set(positions_in(after)) - set(positions_in(before))
    assert fresh
    for node_id in fresh:
        position = positions_in(after)[node_id]
        offset = hypot(position["x"] - 500.0, position["y"] - 25.0)
        assert round(offset) == RING_RADIUS


def test_an_unheld_expansion_places_its_new_nodes_too() -> None:
    """Not only while holding: cytoscape drops a positionless node at the
    origin, and the layout that will move it does not run until 100ms after
    the add - long enough to watch the whole expansion pile into the corner.
    A running layout is free to overrule the ring."""
    app = built_app("", TEST_MODE)
    before = drawn_canvas(app)
    after = expand_first_file(app, before, [])
    assert set(positions_in(after)) > set(positions_in(before))


def test_loading_a_seed_while_holding_the_layout_asks_for_one_run() -> None:
    """A fresh graph has no positions to hold, so the held layout still has to
    run once for it - otherwise every node of it is stacked on the origin."""
    response = press_load(
        built_app("", TEST_MODE), MFS, keep_layout=[HOLD_LAYOUT], runs=3
    )
    assert response["relayout"]["data"] == 4


def test_loading_a_seed_without_holding_the_layout_asks_for_nothing() -> None:
    """autoRefreshLayout already runs it; a bump would run it a second time."""
    assert "relayout" not in press_load(built_app("", TEST_MODE), MFS)


def test_layout_has_a_single_writer() -> None:
    """graph.layout is written by one callback on purpose. A second writer -
    say grow_graph fitting on Load - would race it, and which fit won would
    depend on callback ordering."""
    app = built_app(MFS)
    writers = [key for key in app.callback_map if "graph.layout" in key]
    assert writers == ["graph.layout"]


# --------------------------------------------------------------------------
# status_text
# --------------------------------------------------------------------------


def test_status_text_ok_is_grey() -> None:
    assert status_text("all good").style["color"] == "#555"


def test_status_text_failure_is_red_and_bold() -> None:
    """An empty canvas is indistinguishable from a silent 404, so failures have
    to look different."""
    style = status_text("boom", ok=False).style
    assert style["color"] == "#b00020"
    assert style["fontWeight"] == "600"


def test_status_text_carries_the_message() -> None:
    assert status_text("64 drawn").children == "64 drawn"


# --------------------------------------------------------------------------
# format_value
# --------------------------------------------------------------------------


def test_format_value_short_list_is_joined() -> None:
    assert format_value(["a", "b"]) == "a, b"


def test_format_value_long_list_is_counted() -> None:
    assert format_value(list(range(10))) == "10 items"


def test_format_value_truncates_long_strings() -> None:
    formatted = format_value("x" * 400)
    assert formatted.endswith("...")
    assert len(formatted) == 140


def test_format_value_leaves_short_strings_alone() -> None:
    assert format_value("current") == "current"


def test_format_value_stringifies_scalars() -> None:
    assert format_value(0) == "0"
    assert format_value(False) == "False"


# --------------------------------------------------------------------------
# detail_panel
# --------------------------------------------------------------------------


def test_detail_panel_without_selection_prompts() -> None:
    assert "Click a node" in str(detail_panel(None, "db2_test"))


def test_detail_panel_for_a_group_offers_a_picker() -> None:
    members = [raw_matrix_file(index) for index in range(30)]
    data = group_element(MFS, "RawMatrixFile", members, TEST_MODE)["data"]
    rendered = str(detail_panel(data, TEST_MODE))
    assert "member-pick" in rendered
    assert "fan-out" in rendered
    assert "Fan out all 30" in rendered


def _picker(children: list):
    """The member-pick dropdown out of a group detail panel."""
    return next(
        child for child in children if getattr(child, "id", None) == "member-pick"
    )


def test_group_picker_ticks_nothing_when_no_member_is_drawn() -> None:
    members = [raw_matrix_file(index) for index in range(30)]
    data = group_element(MFS, "RawMatrixFile", members, TEST_MODE)["data"]
    assert _picker(detail_panel(data, TEST_MODE)).value == []


def test_group_picker_ticks_exactly_the_drawn_members() -> None:
    """The tick state is the canvas, so the panel takes it as an argument rather
    than keeping its own record of what was clicked."""
    members = [raw_matrix_file(index) for index in range(30)]
    data = group_element(MFS, "RawMatrixFile", members, TEST_MODE)["data"]
    drawn = [members[2], members[11]]
    picker = _picker(detail_panel(data, TEST_MODE, drawn))
    assert picker.value == drawn
    # every member stays selectable, drawn or not
    assert len(picker.options) == 30


def test_group_picker_explains_that_unticking_removes() -> None:
    members = [raw_matrix_file(index) for index in range(30)]
    data = group_element(MFS, "RawMatrixFile", members, TEST_MODE)["data"]
    assert "untick" in str(detail_panel(data, TEST_MODE))


def test_detail_panel_for_a_node_lists_properties() -> None:
    LatticeNode(TISSUE).object_json = {
        "@id": TISSUE,
        "aliases": ["lab:my_tissue"],
        "status": "current",
    }
    rendered = str(detail_panel(node_element(LatticeNode(TISSUE))["data"], TEST_MODE))
    assert "my_tissue" in rendered
    assert "current" in rendered
    assert "open in API" in rendered


def test_detail_panel_for_an_uncached_node_says_so() -> None:
    rendered = str(detail_panel(node_element(LatticeNode(TISSUE))["data"], TEST_MODE))
    assert "No cached properties" in rendered


def test_detail_panel_api_link_points_at_the_mode_server() -> None:
    rendered = str(detail_panel(node_element(LatticeNode(TISSUE))["data"], TEST_MODE))
    assert f"{TEST_SERVER.rstrip('/')}{TISSUE}" in rendered


# --------------------------------------------------------------------------
# legend and stylesheet
# --------------------------------------------------------------------------


def test_legend_lists_every_node_color() -> None:
    rendered = str(legend())
    for member in NodeColor:
        assert member.name in rendered


def test_legend_swatches_use_css_safe_colors() -> None:
    """The 8-digit enum values must not reach the DOM unconverted - the browser
    ignores #rrggbbaa in a `background` shorthand the same way cytoscape does."""
    rendered = str(legend())
    for member in NodeColor:
        assert member.value not in rendered, member.name
        assert css_color(member.value) in rendered, member.name


def test_stylesheet_marks_unexpanded_nodes() -> None:
    selectors = [rule["selector"] for rule in BASE_STYLESHEET]
    assert "node[!expanded]" in selectors


def test_stylesheet_styles_group_placeholders() -> None:
    selectors = [rule["selector"] for rule in BASE_STYLESHEET]
    assert "node[?is_group]" in selectors


def test_sample_seed_is_a_canonical_path() -> None:
    from graph_db2.cyto_elements import normalize_path

    assert normalize_path(SAMPLE_SEED) == SAMPLE_SEED


@pytest.mark.parametrize("name", list(LAYOUTS))
def test_layout_names_have_a_cytoscape_name(name: str) -> None:
    assert LAYOUTS[name]["name"]
