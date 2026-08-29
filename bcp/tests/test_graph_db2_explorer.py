"""Tests for graph_db2.explorer presentation helpers.

Only the pure helpers are exercised here - the callbacks are closures inside
build_app() and need a live Dash app to reach.
"""

from __future__ import annotations

import pytest

from graph_db2.cyto_elements import (
    css_color,
    edge_element,
    group_element,
    merge_elements,
    node_element,
)
from graph_db2.explorer import (
    BASE_STYLESHEET,
    DEFAULT_LAYOUT,
    LAYOUTS,
    SAMPLE_SEED,
    detail_panel,
    format_value,
    legend,
    status_text,
    suggest_layout,
)
from graph_db2.models import LatticeNode, NodeColor

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    MFS,
    TISSUE,
    clean_graph_state,
    TEST_MODE,
    TEST_SERVER,
    db2_env,
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
