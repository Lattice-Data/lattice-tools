from collections import Counter

import dash_cytoscape as cyto
import requests
from dash import Dash, Input, Output, State, ctx, dcc, html, no_update

from .cyto_elements import (
    DRAW_BUDGET,
    FAN_THRESHOLD,
    css_color,
    drop_node,
    expand,
    explode,
    fan_summary,
    fetch_labels,
    make_gatherer,
    member_options,
    merge_elements,
    neighbors_drawn,
    normalize_path,
    not_yet_drawn,
    object_url,
    promote_members,
    properties_of,
)
from .models import LatticeNode, NodeColor

cyto.load_extra_layouts()  # dagre, cose-bilkent, klay, cola

FIT = {"animate": False, "fit": True, "padding": 30}
LAYOUTS = {
    "dagre (left to right)": {
        "name": "dagre",
        "rankDir": "LR",
        "nodeSep": 18,
        "rankSep": 110,
        **FIT,
    },
    "dagre (top to bottom)": {
        "name": "dagre",
        "rankDir": "TB",
        "nodeSep": 24,
        "rankSep": 90,
        **FIT,
    },
    "breadthfirst": {"name": "breadthfirst", "spacingFactor": 1.1, **FIT},
    "cose-bilkent": {"name": "cose-bilkent", "nodeRepulsion": 9000, **FIT},
    "concentric": {"name": "concentric", "minNodeSpacing": 30, **FIT},
}
DEFAULT_LAYOUT = "dagre (left to right)"
# only valid on DEFAULT_MODE's server; see main()
SAMPLE_SEED = "/matrix_file_sets/f1ef71ee-98d8-4145-84a7-24b68bcc769e/"

BASE_STYLESHEET = [
    {
        "selector": "node",
        "style": {
            "background-color": "data(color)",
            "label": "data(label)",
            "font-size": "9px",
            "color": "#333",
            "text-valign": "bottom",
            "text-halign": "center",
            "text-margin-y": "4px",
            "text-wrap": "ellipsis",
            "text-max-width": "110px",
            "width": 26,
            "height": 26,
            "border-width": 1.5,
            "border-color": "#666",
        },
    },
    # dashed ring = neighbors not resolved yet, i.e. still clickable
    {
        "selector": "node[!expanded]",
        "style": {"border-style": "dashed", "border-color": "#aaa"},
    },
    {
        "selector": "node[?is_group]",
        "style": {
            "shape": "round-rectangle",
            "width": 122,
            "height": 24,
            "border-width": 2,
            "border-color": "#888",
            "font-size": "10px",
            "text-valign": "center",
            "text-margin-y": 0,
            "text-max-width": "116px",
        },
    },
    {
        "selector": "edge",
        "style": {"width": 1.2, "line-color": "#c4c4c4", "curve-style": "bezier"},
    },
    {
        "selector": "node:selected",
        "style": {"border-width": 3.5, "border-color": "#1f6feb"},
    },
]

# a column so the type filter can sit in a pinned footer while the legend and
# the (unbounded) property table scroll above it
PANEL_STYLE = {
    "width": "330px",
    "flex": "0 0 330px",
    "padding": "10px 14px",
    "borderLeft": "1px solid #e2e2e2",
    "display": "flex",
    "flexDirection": "column",
    "height": "100%",
    "boxSizing": "border-box",
    "overflow": "hidden",
    "fontFamily": "system-ui, sans-serif",
    "fontSize": "12px",
}
# minHeight 0 or the flex child refuses to shrink below its content and the
# footer gets pushed off the panel instead of the content scrolling
PANEL_SCROLL_STYLE = {"flex": "1 1 auto", "overflowY": "auto", "minHeight": 0}
PANEL_FOOTER_STYLE = {
    "flex": "0 0 auto",
    "borderTop": "1px solid #e2e2e2",
    "paddingTop": "8px",
    "marginTop": "8px",
    # its own scroll, so a graph with many types can't crowd out the details
    "maxHeight": "38%",
    "overflowY": "auto",
}


def status_text(message: str, ok: bool = True) -> html.Span:
    """Status line; failures are red, because an empty canvas looks identical
    to a load that silently 404'd."""
    return html.Span(
        message,
        style={"color": "#555"} if ok else {"color": "#b00020", "fontWeight": "600"},
    )


def suggest_layout(elements: list[dict]) -> str:
    """
    Pick a layout to match the graph's shape.

    dagre reads well for lineage chains and badly for stars - a 64-wide fan
    becomes one ~2900px column, which looks like nothing loaded. Concentric puts
    the same fan in a readable ring. Only used for a freshly loaded seed; once
    the user has picked a layout, expansions leave it alone.
    """
    nodes = [element for element in elements if "source" not in element["data"]]
    if len(nodes) < 12:
        return DEFAULT_LAYOUT

    degree = Counter()
    for element in elements:
        if "source" in element["data"]:
            degree[element["data"]["source"]] += 1
            degree[element["data"]["target"]] += 1

    hub = max(degree.values(), default=0)
    return "concentric" if hub >= 0.8 * (len(nodes) - 1) else DEFAULT_LAYOUT


def legend() -> html.Div:
    swatches = [
        html.Div(
            [
                html.Span(
                    style={
                        "display": "inline-block",
                        "width": "11px",
                        "height": "11px",
                        "borderRadius": "50%",
                        "border": "1px solid #666",
                        "background": css_color(member.value),
                        "marginRight": "6px",
                        "verticalAlign": "middle",
                    }
                ),
                html.Span(member.name),
            ],
            style={"marginBottom": "3px"},
        )
        for member in NodeColor
    ]
    return html.Div(swatches)


def format_value(value) -> str:
    if isinstance(value, list):
        if len(value) > 4:
            return f"{len(value)} items"
        return ", ".join(str(item) for item in value)
    text = str(value)
    return text if len(text) <= 140 else text[:137] + "..."


def detail_panel(node_data: dict | None, mode: str) -> list:
    if not node_data:
        return [html.P("Click a node to expand it and see its properties.")]

    if node_data.get("is_group"):
        members = node_data["members"]
        return [
            html.H4(node_data["label"], style={"margin": "0 0 4px"}),
            html.P(
                f"{len(members)} {node_data['node_type']} references, collapsed "
                "because the whole fan is over the draw budget. Search and tick "
                "as many as you want - clearing one leaves it on the canvas.",
                style={"marginTop": 0},
            ),
            dcc.Dropdown(
                id="member-pick",
                options=member_options(members, mode),
                placeholder=f"search {len(members)} {node_data['node_type']}...",
                multi=True,
                optionHeight=44,
            ),
            html.Button(
                f"Fan out all {len(members)}",
                id="fan-out",
                n_clicks=0,
                style={"marginTop": "8px"},
            ),
        ]

    uuid_path = node_data["id"]
    properties = properties_of(uuid_path)
    rows = [
        html.Tr(
            [
                html.Td(
                    key,
                    style={
                        "verticalAlign": "top",
                        "paddingRight": "8px",
                        "color": "#666",
                        "whiteSpace": "nowrap",
                    },
                ),
                html.Td(format_value(value), style={"wordBreak": "break-word"}),
            ]
        )
        for key, value in properties.items()
    ]

    return [
        html.H4(node_data["label"], style={"margin": "0 0 2px"}),
        html.Div(
            node_data["node_type"], style={"color": "#666", "marginBottom": "6px"}
        ),
        html.A(
            "open in API",
            href=object_url(uuid_path, mode),
            target="_blank",
            style={"fontSize": "11px"},
        ),
        html.Hr(style={"margin": "8px 0"}),
        html.Table(rows, style={"borderCollapse": "collapse", "width": "100%"})
        if rows
        else html.P("No cached properties - click the node to fetch it."),
    ]


def build_app(seed: str, mode: str, fetch_new: bool) -> Dash:
    gatherer = make_gatherer(mode, fetch_new)
    LatticeNode.mode = mode

    elements: list[dict] = []
    if not seed:
        status = status_text(f"On {mode}. Paste an object path and press Load.")
        notice = [html.P(f"Connected to {mode}. Load an object path to start.")]
    else:
        try:
            # so the input box shows the canonical form of whatever was typed
            seed = normalize_path(seed)
            seed_nodes, seed_edges = expand(seed, gatherer, mode=mode)
            elements = merge_elements([], seed_nodes, seed_edges)
            status = status_text(f"{seed} - {fan_summary(seed_nodes)}")
            notice = detail_panel(None, mode)
        except (requests.HTTPError, ValueError, KeyError) as error:
            # an empty canvas is indistinguishable from a silent failure, so say
            # so in the panel too rather than only in the toolbar
            status = status_text(f"Could not load {seed}", ok=False)
            notice = [
                html.H4(
                    "Nothing loaded", style={"margin": "0 0 6px", "color": "#b00020"}
                ),
                html.P(f"{seed} could not be fetched from {mode}:"),
                html.Pre(
                    str(error),
                    style={
                        "whiteSpace": "pre-wrap",
                        "wordBreak": "break-word",
                        "background": "#f6f6f6",
                        "padding": "6px",
                        "fontSize": "11px",
                    },
                ),
                html.P("A path from one deployment will not resolve on another."),
            ]
    initial_layout = suggest_layout(elements)

    # the member picker only exists once a group node is selected
    app = Dash(__name__, title="DB2 graph explorer", suppress_callback_exceptions=True)
    app.layout = html.Div(
        [
            html.Div(
                [
                    dcc.Input(
                        id="seed",
                        value=seed,
                        debounce=True,
                        style={"width": "440px", "fontFamily": "monospace"},
                    ),
                    html.Button("Load", id="load", n_clicks=0),
                    html.Label(
                        [
                            "draw up to ",
                            dcc.Input(
                                id="fan",
                                type="number",
                                value=DRAW_BUDGET,
                                min=1,
                                step=10,
                                style={"width": "66px"},
                            ),
                            " neighbors",
                        ],
                        title=(
                            "Fans this size or smaller are drawn in full. Past it, "
                            f"any type with more than {FAN_THRESHOLD} members "
                            "collapses into a placeholder you can search."
                        ),
                        style={"whiteSpace": "nowrap", "color": "#555"},
                    ),
                    dcc.Dropdown(
                        id="layout-choice",
                        options=list(LAYOUTS),
                        value=initial_layout,
                        clearable=False,
                        style={"width": "220px"},
                    ),
                    html.Span(
                        status, id="status", style={"color": "#555", "fontSize": "12px"}
                    ),
                ],
                style={
                    "display": "flex",
                    "flex": "0 0 auto",
                    "gap": "10px",
                    "alignItems": "center",
                    "padding": "8px 12px",
                    "borderBottom": "1px solid #e2e2e2",
                    "fontFamily": "system-ui, sans-serif",
                    "fontSize": "12px",
                },
            ),
            html.Div(
                [
                    # Cytoscape sizes itself from its parent, so the flex child
                    # needs a resolved height - not an intermediate wrapper that
                    # collapses to 0.
                    html.Div(
                        cyto.Cytoscape(
                            id="graph",
                            elements=elements,
                            layout=LAYOUTS[initial_layout],
                            stylesheet=BASE_STYLESHEET,
                            style={"width": "100%", "height": "100%"},
                            boxSelectionEnabled=True,
                            # first layout runs before the flex child has its
                            # final width, so re-fit once the container settles
                            responsive=True,
                        ),
                        style={"flex": "1", "minWidth": "0", "height": "100%"},
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div("Legend", style={"color": "#666"}),
                                    legend(),
                                    html.Hr(style={"margin": "10px 0"}),
                                    html.Div(notice, id="details"),
                                ],
                                style=PANEL_SCROLL_STYLE,
                            ),
                            html.Div(
                                [
                                    html.Div("Show types", style={"color": "#666"}),
                                    dcc.Checklist(
                                        id="type-filter", options=[], value=[]
                                    ),
                                ],
                                style=PANEL_FOOTER_STYLE,
                            ),
                        ],
                        style=PANEL_STYLE,
                    ),
                ],
                # minHeight 0 lets this row shrink to the space the toolbar
                # leaves, instead of forcing its content height onto the page
                style={"display": "flex", "flex": "1 1 auto", "minHeight": 0},
            ),
        ],
        # Fixed to the viewport rather than sized with calc(100vh - toolbar):
        # that magic number ignored the body's default margin, so the panel
        # (and its pinned footer) hung below the fold. Colours are pinned light
        # so the inline text colours stay legible whatever the browser prefers.
        style={
            "position": "fixed",
            "top": 0,
            "left": 0,
            "right": 0,
            "bottom": 0,
            "display": "flex",
            "flexDirection": "column",
            "background": "#ffffff",
            "color": "#222",
            "colorScheme": "light",
        },
    )

    @app.callback(
        Output("graph", "elements"),
        Output("status", "children"),
        Output("layout-choice", "value"),
        Input("load", "n_clicks"),
        Input("graph", "tapNodeData"),
        State("seed", "value"),
        State("fan", "value"),
        State("graph", "elements"),
        prevent_initial_call=True,
    )
    def grow_graph(_clicks, tapped, seed_value, fan, elements):
        budget = int(fan) if fan else 0
        elements = elements or []
        loading = ctx.triggered_id == "load"

        if loading:
            if not seed_value:
                return (
                    no_update,
                    status_text("Enter an object path to load.", ok=False),
                    no_update,
                )
            try:
                target = normalize_path(seed_value)
            except ValueError as error:
                return no_update, status_text(str(error), ok=False), no_update
            elements = []
        else:
            if not tapped:
                return no_update, no_update, no_update
            target = tapped["id"]
            # a group tap only opens the side panel; fanning out hundreds of
            # nodes needs the explicit button there
            if tapped.get("is_group"):
                return (
                    no_update,
                    status_text(f"{tapped['label']} - choose from the panel."),
                    no_update,
                )
            # canvas state, not the fetch cache: a reload empties the canvas
            # while leaving every profile cached, and those nodes still need
            # their edges drawn again
            if neighbors_drawn(elements, target):
                return (
                    no_update,
                    status_text(f"{tapped['label']} is already expanded."),
                    no_update,
                )

        try:
            nodes, edges = expand(target, gatherer, mode=mode, draw_budget=budget)
        except (requests.HTTPError, ValueError, KeyError) as error:
            return (
                no_update,
                status_text(f"Could not load {target}: {error}", ok=False),
                no_update,
            )
        action = f"{target} - {fan_summary(nodes)}"

        merged = merge_elements(elements, nodes, edges)
        node_count = sum(1 for element in merged if "source" not in element["data"])
        return (
            merged,
            status_text(
                f"{action} - {node_count} nodes, {len(merged) - node_count} edges"
            ),
            # a new seed replaces the graph, so re-pick the layout for its shape;
            # an expansion builds on what the user is already looking at
            suggest_layout(merged) if loading else no_update,
        )

    @app.callback(
        Output("type-filter", "options"),
        Output("type-filter", "value"),
        Input("graph", "elements"),
        State("type-filter", "value"),
        State("type-filter", "options"),
    )
    def sync_type_filter(elements, selected, previous_options):
        present = sorted(
            {
                element["data"]["node_type"]
                for element in elements or []
                if "node_type" in element["data"]
            }
        )
        # types that just appeared default to visible; keep existing choices
        newly_seen = [name for name in present if name not in (previous_options or [])]
        return present, [
            name for name in present if name in (selected or []) or name in newly_seen
        ]

    @app.callback(
        Output("graph", "stylesheet"),
        Input("type-filter", "value"),
        State("type-filter", "options"),
    )
    def filter_by_type(selected, options):
        hidden = [name for name in (options or []) if name not in (selected or [])]
        return BASE_STYLESHEET + [
            {
                "selector": f'node[node_type = "{name}"]',
                "style": {"display": "none"},
            }
            for name in hidden
        ]

    @app.callback(Output("graph", "layout"), Input("layout-choice", "value"))
    def choose_layout(choice):
        return LAYOUTS[choice]

    @app.callback(
        Output("graph", "elements", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("member-pick", "value"),
        Input("fan-out", "n_clicks"),
        State("graph", "tapNodeData"),
        State("graph", "elements"),
        prevent_initial_call=True,
    )
    def draw_group_members(picked, _clicks, tapped, elements):
        if not tapped or not tapped.get("is_group"):
            return no_update, no_update

        elements = elements or []
        if ctx.triggered_id == "fan-out":
            nodes, edges = explode(tapped, gatherer, mode=mode)
            # the placeholder is gone once its members are all on the canvas
            elements = drop_node(elements, tapped["id"])
            action = f"fanned out {len(nodes)} {tapped['node_type']}"
        else:
            # the dropdown reports its whole selection on every change, so only
            # the members not already on the canvas are worth fetching
            fresh = not_yet_drawn(elements, picked or [])
            if not fresh:
                return no_update, no_update
            nodes, edges = promote_members(fresh, tapped["parent_path"], gatherer, mode)
            action = (
                f"added {nodes[0]['data']['label']}"
                if len(fresh) == 1
                else f"added {len(fresh)} {tapped['node_type']}"
            )

        merged = merge_elements(elements, nodes, edges)
        node_count = sum(1 for element in merged if "source" not in element["data"])
        return (
            merged,
            status_text(
                f"{action} - {node_count} nodes, {len(merged) - node_count} edges"
            ),
        )

    # prevent_initial_call keeps the startup notice; without it this fires once
    # with tapNodeData=None and overwrites any load error with the generic hint
    @app.callback(
        Output("details", "children"),
        Input("graph", "tapNodeData"),
        prevent_initial_call=True,
    )
    def show_details(tapped):
        if tapped and tapped.get("is_group"):
            # one batched report so the picker is searchable by alias rather
            # than by uuid prefix
            fetch_labels(tapped["members"], gatherer)
        return detail_panel(tapped, mode)

    return app
