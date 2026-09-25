"""
Cytoscape element construction and one-hop expansion.

graphing.py materializes the whole reachable graph up front so pyvis can write
it into a static file. This module walks outward one node at a time instead:
explorer.py asks for a node's neighbors only when a user clicks it, so the
payload stays small no matter how large the reachable graph is.

Only genuinely wide fans are held back. A fan of DRAW_BUDGET or fewer is drawn
in full, so a 64-file MatrixFileSet looks like its 64 files. Past that, any type
with more than FAN_THRESHOLD members is left off the canvas - a 500-wide fan is
both unreadable and 500 label fetches.

Every type discovered so far gets one placeholder, unconnected, at the top of
its column. It lists every node of that type the canvas knows of - drawn, or
referenced by an expanded node - and its picker draws and clears them. See
sync_placeholders().

Edges are undirected. LatticeNode.get_ids() walks the JSON on both sides of a
reference, so a Library and its FileSet each name the other; picking a direction
from discovery order would draw arrows that don't mean lineage.
"""

from collections import Counter, defaultdict
from collections.abc import Collection
from math import cos, pi, sin
from urllib.parse import quote, urljoin

import requests
from .connection import get_connection
from .constants import DEFAULT_MODE, FETCH_NEW
from .models import (
    LatticeNode,
    GraphDB2Error,
    NodeColor,
    batch_request_chunk_and_fetch,
    group_batch_request,
)
from .schema import create_configs

from db2_flattener.gather.gatherer import DB2Gatherer

UNMAPPED_COLOR = "#e8e8e8"
# neighbors drawn individually before any are held back
DRAW_BUDGET = 25
# once over budget, a type this size or smaller is still drawn rather than held back
FAN_THRESHOLD = 25
# no object path starts with a colon, so a placeholder id cannot collide
PLACEHOLDER_PREFIX = "::type:"
# Where a new node goes before any layout has had a say: RING_RADIUS out from
# the node it was expanded from, NODE_PITCH of room each along the ring. See
# seed_positions().
RING_RADIUS = 170
NODE_PITCH = 70

# "columns by type" geometry: the gap between columns, the room each node gets
# down one, and how many ordering passes to make over them. See
# column_positions().
COLUMN_GAP = 240
ROW_PITCH = 46
ORDERING_SWEEPS = 2
# Left to right, roughly the order an experiment happens in, so lineage reads
# across the screen. Each entry is one column, named by the legend's buckets
# rather than by api_name - Tissue, CellLine and Organoid share the Biosample
# column the same way they share its colour. A type the legend does not map
# lands in a trailing column of its own.
COLUMN_TYPES = (
    ("Donor",),
    ("Biosample",),
    (
        "GeneticModification",
        "Treatment",
        "ExperimentalCondition",
        "TabularFile",
    ),
    ("Library",),
    ("SequenceFileSet",),
    ("SequenceFile",),
    (
        "RawMatrixFile",
        "ProcessedMatrixFile",
    ),
    ("MatrixFileSet",),
)
_COLUMN_OF = {name: index for index, names in enumerate(COLUMN_TYPES) for name in names}
UNPLACED_COLUMN = len(COLUMN_TYPES)

# Paths whose complete profile has been fetched. Batch reports only carry the
# fields in OBJECT_CONFIG, so a node known only from a batch response can have
# an incomplete neighbor set and must be re-fetched before it is expanded.
_fully_fetched: set[str] = set()


def make_gatherer(mode: str = DEFAULT_MODE, fetch_new: bool = FETCH_NEW) -> DB2Gatherer:
    """DB2Gatherer for the server behind `mode`, configured from constants.yaml"""
    connection = get_connection(mode)
    config = create_configs(mode, fetch_new)
    return DB2Gatherer(connection, config)


def fetch_full(node: LatticeNode) -> dict:
    """
    Authenticated GET of a node's complete profile.

    LatticeNode.fetch_object_json() sends no credentials, so it only resolves
    objects released to the public; expansion needs whatever the caller's key
    can see, or the neighbor set comes back short.
    """
    connection = node.connection
    response = requests.get(
        urljoin(connection.server, node.uuid_path),
        auth=connection.auth,
        headers=connection.headers,
    )
    response.raise_for_status()
    return response.json()


def normalize_path(path: str) -> str:
    """
    Canonicalize a seed string to a URL path with leading and trailing slashes.

    Three spellings address an object on a DB2 instance, and this accepts all of
    them without deciding which it was handed:
        /{object_type}/{uuid}/
        /{alias}/
        /{uuid}/
    Only the server can say which one it is, so the two checks here are the ones
    that need no request: something has to survive the slash stripping, and a
    full URL is not a path on the instance.
    """
    cleaned = path.strip().strip("/")
    if not cleaned:
        raise GraphDB2Error(f"Expected an object path, alias or uuid, got {path!r}")
    if cleaned.lower().startswith(("http://", "https://")):
        raise GraphDB2Error(f"Expected a path on the instance, not a URL: {path!r}")
    # '//' means an empty segment, which no spelling has - and left alone it
    # turns a seed into a request for a different path than the user typed
    if "//" in cleaned:
        raise GraphDB2Error(f"Expected no empty path segments, got {path!r}")
    return f"/{cleaned}/"


def canonical_id(json_object: dict) -> str | None:
    """
    The '@id' of a JSON response, if it names a single object.

    A collection ('/matrix_file_sets/') and the portal root both answer 200
    without naming one, and LatticeNode cannot parse either - it splits an '@id'
    into exactly two segments.
    """
    type_and_uuid = json_object.get("@id")
    if not isinstance(type_and_uuid, str):
        return None
    segments = [segment for segment in type_and_uuid.split("/") if segment]
    return type_and_uuid if len(segments) == 2 else None


def normalize_to_type_and_uuid(path: str, mode: str = DEFAULT_MODE) -> str:
    """
    Resolve any accepted seed spelling to the canonical '/{object_type}/{uuid}/'.

    An alias and a bare uuid each address an object but neither is usable as a
    node id, and only the server can say which object they name. So the path is
    requested once: the profile it answers with lands in the LatticeNode cache,
    and its '@id' is the canonical form everything downstream keys off.

    A path already in the cache is already an '@id' - every key in that dict came
    from one - so it resolves to itself without a request. That is what keeps
    tapping an already-fetched node from costing an extra round trip.
    """
    if path in LatticeNode._cache:
        return path

    # outside the try: a missing credential exits, and naming the server in the
    # failure message needs the connection to have been built
    connection = get_connection(mode)
    # an alias is free-form text, and an unencoded '#' or '?' in one silently
    # truncates the request path rather than reaching the object
    url = urljoin(connection.server, quote(path, safe="/:@"))
    try:
        response = requests.get(
            url,
            auth=connection.auth,
            headers=connection.headers,
            timeout=60,
        )
        response.raise_for_status()
        json_object = response.json()
    except requests.RequestException as exc:
        # covers the 404 for an unknown alias, the 403 for an object this key
        # cannot see, a dead server, and a 200 that is not JSON
        raise GraphDB2Error(
            f"Could not resolve seed {path!r} on {connection.server}: {exc}"
        ) from exc

    type_and_uuid = canonical_id(json_object)
    if type_and_uuid is None:
        raise GraphDB2Error(
            f"Seed {path!r} resolved to {json_object.get('@id')!r} on "
            f"{connection.server}, which is not a single object - expected "
            "'/{object_type}/{uuid}/'"
        )

    LatticeNode._cache[type_and_uuid] = json_object
    return type_and_uuid


def resolve_seed(seed: str, mode: str = DEFAULT_MODE) -> str:
    """
    Whatever the user typed, as a node id.

    The two steps are never useful apart, and calling them separately is how a
    caller ends up resolving a db2_demo seed against db2_prod: the mode has to
    be threaded through the second one, and there is nothing in the first to
    suggest it.
    """
    return normalize_to_type_and_uuid(normalize_path(seed), mode)


def object_url(uuid_path: str, mode: str = DEFAULT_MODE) -> str:
    return urljoin(get_connection(mode).server, uuid_path)


def is_fully_fetched(uuid_path: str) -> bool:
    """
    Whether this node's complete profile is in the cache.

    A process-wide fact, not a canvas one: pressing Load wipes the elements but
    leaves the cache warm, so this must not be used to decide whether a node's
    edges are on screen. See neighbors_drawn().
    """
    return uuid_path in _fully_fetched


def neighbors_drawn(elements: list[dict], node_id: str) -> bool:
    """
    Whether this node's neighbors are already on the given canvas.

    The `expanded` flag lives on the element, so it resets with the canvas.
    Reading the fetch cache here instead made a re-Load unrecoverable: every
    node clicked before the reload was still "expanded", so the second pass
    refused to redraw any of their edges.
    """
    for element in elements:
        if element["data"]["id"] == node_id:
            return bool(element["data"].get("expanded"))
    return False


def css_color(value: str) -> str:
    """
    NodeColor values are 8-digit #rrggbbaa. vis.js accepts that; cytoscape.js
    rejects it and silently drops the mapping, so normalize to something CSS
    understands.
    """
    if len(value) != 9 or not value.startswith("#"):
        return value

    red, green, blue, alpha = (
        int(value[index : index + 2], 16) for index in (1, 3, 5, 7)
    )
    if alpha == 255:
        return value[:7]
    return f"rgba({red}, {green}, {blue}, {round(alpha / 255, 3)})"


def color_for(node: LatticeNode) -> str:
    """NodeColor with a fallback: _missing_ returns None for unmapped types"""
    try:
        return css_color(NodeColor(node.schema_ids.api_name).value)
    except ValueError:
        return UNMAPPED_COLOR


def label_for(node: LatticeNode) -> str:
    """
    Best available human label, without triggering a request.

    Reads the cache directly rather than going through LatticeNode.alias, which
    would fetch on a miss and raise IndexError on an empty aliases list.
    """
    json_object = LatticeNode._cache.get(node.uuid_path)
    if json_object is None:
        return node.uuid[:8]

    aliases = json_object.get("aliases") or []
    if aliases:
        # aliases are 'lab-name:local-id'; the local id is the informative half
        return aliases[0].split(":", 1)[-1]
    return json_object.get("accession") or node.uuid[:8]


def node_element(node: LatticeNode, expanded: bool = False) -> dict:
    """
    One cytoscape node. `expanded` means "its neighbors are on this canvas", so
    only the centre of an expansion sets it - a node re-emitted as someone
    else's neighbor arrives False and merge_elements() keeps the older True.
    """
    return {
        "data": {
            "id": node.uuid_path,
            "label": label_for(node),
            "node_type": node.schema_ids.api_name,
            "color": color_for(node),
            "expanded": expanded,
        }
    }


def edge_element(one_path: str, other_path: str) -> dict:
    """Undirected edge with an order-independent id, so A-B and B-A dedupe"""
    source, target = sorted((one_path, other_path))
    return {"data": {"id": f"{source}--{target}", "source": source, "target": target}}


def placeholder_id(api_name: str) -> str:
    return f"{PLACEHOLDER_PREFIX}{api_name}"


def group_label(api_name: str, drawn: int, total: int) -> str:
    # two lines, so a long type name's ellipsis cannot eat the count
    return f"{api_name}\n{drawn} out of {total}"


def placeholder_element(api_name: str, members: list[str], drawn: int) -> dict:
    """
    The one placeholder for a type. `members` is every node of it the canvas
    knows of; the label says how many of them are drawn.
    """
    return {
        "data": {
            "id": placeholder_id(api_name),
            "label": group_label(api_name, drawn, len(members)),
            "node_type": api_name,
            "color": color_for(LatticeNode(members[0])),
            "expanded": False,
            "is_group": True,
            "members": members,
        }
    }


def fetch_labels(paths: Collection[str], gatherer: DB2Gatherer) -> None:
    """
    Batch report per type for anything not already cached - the same
    group -> chunk_and_fetch -> cache_update path that
    fetch_all_references_from_cache() uses, one hop wide instead of the whole
    reachable graph.
    """
    grouped = group_batch_request(paths)
    if grouped:
        LatticeNode.batch_cache_update(batch_request_chunk_and_fetch(grouped, gatherer))


def expand(
    uuid_path: str,
    gatherer: DB2Gatherer,
    mode: str = DEFAULT_MODE,
    draw_budget: int = DRAW_BUDGET,
    on_canvas: Collection[str] = (),
) -> tuple[list[dict], list[dict]]:
    """
    Resolve one node's neighbors into cytoscape elements.

    A fan of draw_budget or fewer is drawn in full - a 64-file MatrixFileSet is
    exactly what someone opening it wants to see. Only past that are the big
    types held back, since a 512-wide fan lays out ~23,000px tall and is
    unreadable at any zoom. Nothing stands in for them here: once this node is
    expanded they are in its type's placeholder, which sync_placeholders()
    builds from the whole canvas.

    Neighbors already `on_canvas` cost nothing to draw, so they never count
    toward either limit and always get their edge.

    Costs one authenticated GET for the node itself plus one batched report per
    type actually drawn. Held-back types cost nothing until the user draws one,
    since a neighbor's type is readable from its path alone.
    """
    # everything downstream keys off node ids, so work from the canonical form
    # rather than whatever the caller typed. Free for a node already in the
    # cache, which is every node the user can tap.
    uuid_path = resolve_seed(uuid_path, mode)
    node = LatticeNode(uuid_path)
    if node.uuid_path not in _fully_fetched:
        node.object_json = fetch_full(node)  # setter writes through to _cache
        _fully_fetched.add(node.uuid_path)

    neighbors = sorted(node.neighbors)
    by_type = defaultdict(list)
    for neighbor in neighbors:
        by_type[LatticeNode(neighbor).schema_ids.api_name].append(neighbor)

    on_canvas = set(on_canvas)
    drawn: list[str] = []
    if draw_budget and sum(1 for n in neighbors if n not in on_canvas) > draw_budget:
        for paths in by_type.values():
            waiting = [path for path in paths if path not in on_canvas]
            if len(waiting) > FAN_THRESHOLD:
                drawn.extend(path for path in paths if path in on_canvas)
            else:
                drawn.extend(paths)
    else:
        drawn = neighbors

    fetch_labels(drawn, gatherer)

    nodes = [node_element(node, expanded=True)]
    edges = []
    for path in drawn:
        nodes.append(node_element(LatticeNode(path)))
        edges.append(edge_element(uuid_path, path))

    return nodes, edges


def promote_members(paths: Collection[str], gatherer: DB2Gatherer) -> list[dict]:
    """
    Nodes for members picked from a placeholder. Their edges come from
    link_known_edges(), since a placeholder has no parent to hang them off.

    Takes the whole selection at once so a multi-pick costs one batched report
    per type rather than one per node.
    """
    fetch_labels(paths, gatherer)
    return [node_element(LatticeNode(path)) for path in paths]


def member_options(members: list[str], mode: str = DEFAULT_MODE) -> list[dict]:
    """
    Dropdown options for a group's members, sorted by label.

    Call fetch_labels() on the members first, or every option falls back to a
    uuid prefix and the dropdown is unsearchable.
    """
    options = [
        {"label": label_for(LatticeNode(path)), "value": path} for path in members
    ]
    return sorted(options, key=lambda option: option["label"])


def fan_summary(nodes: list[dict]) -> str:
    """
    What an expansion actually produced, counting what it held back too -
    "1 drawn" for a fan of 65 reads as a failure. `nodes` is expand()'s, whose
    first entry is the node that was expanded.
    """
    centre = nodes[0]["data"]["id"]
    shown = {node["data"]["id"] for node in nodes}
    held_back = Counter(
        LatticeNode(path).schema_ids.api_name
        for path in known_neighbors(centre)
        if path not in shown
    )

    parts = [f"{len(nodes) - 1} drawn"] if len(nodes) > 1 else []
    parts += [
        f"{count} {api_name} in its placeholder"
        for api_name, count in sorted(held_back.items())
    ]
    return ", ".join(parts) or "no neighbors"


def known_neighbors(uuid_path: str) -> set[str]:
    """
    A node's neighbors if its full profile is cached, else nothing. Never
    fetches: a batch report's partial profile would give a short list, and a
    miss would cost a request per node on every canvas change.
    """
    if uuid_path not in _fully_fetched or uuid_path not in LatticeNode._cache:
        return set()
    return LatticeNode(uuid_path).neighbors


def not_yet_drawn(elements: list[dict], paths: Collection[str]) -> list[str]:
    """The subset of `paths` not yet on the canvas, in selection order"""
    present = {element["data"]["id"] for element in elements}
    return [path for path in paths if path not in present]


def already_drawn(elements: list[dict], paths: Collection[str]) -> list[str]:
    """
    The subset of `paths` currently on the canvas, in the given order.

    Backs the group picker's tick state: the dropdown shows the canvas rather
    than a separate record of what was clicked, so the two cannot drift.
    """
    present = {element["data"]["id"] for element in elements}
    return [path for path in paths if path in present]


def merge_elements(
    existing: list[dict], new_nodes: list[dict], new_edges: list[dict]
) -> list[dict]:
    """
    Fold an expansion into the elements already on screen.

    Nodes are replaced (a stub gains its real label once fetched), edges are
    kept on first sight. Two things survive replacement. `expanded`, because
    expanding B re-emits its neighbor A as an unexpanded stub and dropping A's
    flag would offer a second, edge-free expansion of a node already opened.
    And `position`, because that is canvas state the browser owns and a node
    rebuilt here has none - dropping it teleports an already-drawn node to the
    origin the moment the layout is held and cannot put it back.
    """
    by_id = {element["data"]["id"]: element for element in existing}

    for element in new_nodes:
        node_id = element["data"]["id"]
        previous = by_id.get(node_id)
        if previous and previous["data"].get("expanded"):
            element = {**element, "data": {**element["data"], "expanded": True}}
        if previous and "position" in previous and "position" not in element:
            element = {**element, "position": previous["position"]}
        by_id[node_id] = element
    for element in new_edges:
        by_id.setdefault(element["data"]["id"], element)

    return list(by_id.values())


def settle(elements: list[dict]) -> list[dict]:
    """
    The canvas after a draw or removal: every known edge between drawn nodes
    in place, and the placeholders rebuilt to match. Each callback that writes
    the elements finishes with this.
    """
    return sync_placeholders(link_known_edges(elements))


def link_known_edges(elements: list[dict]) -> list[dict]:
    """
    An edge from every expanded node to each of its neighbors on the canvas.

    An expansion only links the node that was clicked, and a member picked from
    a placeholder has no parent at all. Without this, a library drawn from one
    CellLine would stay unlinked from the three others that reference it.
    """
    present = {
        element["data"]["id"]
        for element in elements
        if "source" not in element["data"] and not element["data"].get("is_group")
    }
    edges = [
        edge_element(element["data"]["id"], neighbor)
        for element in elements
        if element["data"].get("expanded")
        for neighbor in sorted(known_neighbors(element["data"]["id"]) & present)
    ]
    return merge_elements(elements, [], edges)


def sync_placeholders(elements: list[dict]) -> list[dict]:
    """
    One placeholder per type the canvas knows of, rebuilt from scratch.

    A type's members are its nodes on the canvas plus every one an expanded
    node references, so the picker can clear a drawn node as well as draw a
    held-back one, and a type with nothing held back still has one. Built from
    the whole canvas rather than kept per click: when two CellLines reference
    the same 48 libraries, one list covers both, and it cannot fall out of step
    with what is on screen. A placeholder keeps the position it had.
    """
    kept = [element for element in elements if not element["data"].get("is_group")]
    previous = {
        element["data"]["id"]: element
        for element in elements
        if element["data"].get("is_group")
    }
    present = {
        element["data"]["id"] for element in kept if "source" not in element["data"]
    }
    known = set(present)
    for element in kept:
        if element["data"].get("expanded"):
            known |= known_neighbors(element["data"]["id"])

    by_type: dict[str, list[str]] = defaultdict(list)
    for path in sorted(known):
        by_type[LatticeNode(path).schema_ids.api_name].append(path)

    placeholders = []
    for api_name, members in sorted(by_type.items()):
        drawn = sum(1 for path in members if path in present)
        element = placeholder_element(api_name, members, drawn)
        old = previous.get(element["data"]["id"])
        if old and "position" in old:
            element["position"] = old["position"]
        placeholders.append(element)

    return kept + place_placeholders(placeholders, kept)


def place_placeholders(placeholders: list[dict], elements: list[dict]) -> list[dict]:
    """
    A position for each new placeholder, above the top of its column as it
    stands.

    Only matters while the layout is held: otherwise the column layout re-runs
    and puts them there itself, and every other layout hides them. A type with
    nothing drawn has no column yet, so it starts one to the right of the rest.
    """
    # placeholders that kept a position count too, so a new one stacks above
    # them rather than on top of them
    positioned = [
        (element["data"].get("node_type"), element["position"])
        for element in [*elements, *placeholders]
        if "position" in element and "source" not in element["data"]
    ]
    if not positioned:
        return placeholders

    right = max(position["x"] for _, position in positioned) + COLUMN_GAP
    top = min(position["y"] for _, position in positioned)
    placed = []
    for element in placeholders:
        if "position" in element:
            placed.append(element)
            continue
        column = column_of(element["data"]["node_type"])
        peers = [
            position
            for node_type, position in positioned
            if column_of(node_type) == column
        ]
        if peers:
            highest = min(peers, key=lambda position: position["y"])
            position = {"x": highest["x"], "y": highest["y"] - ROW_PITCH}
        else:
            position = {"x": right, "y": top - ROW_PITCH}
            right += COLUMN_GAP
        positioned.append((element["data"]["node_type"], position))
        placed.append({**element, "position": position})
    return placed


def ring_offsets(count: int) -> list[tuple[float, float]]:
    """
    `count` offsets from a centre, on as many rings as it takes to leave
    NODE_PITCH between neighbors sharing one.

    A single ring does not scale: 60 nodes on it either overlap or sit so far
    out that the node they belong to is off screen. Even-numbered rings are
    turned half a step so a wide fan does not come out as spokes.
    """
    offsets: list[tuple[float, float]] = []
    ring = 1
    while len(offsets) < count:
        radius = RING_RADIUS * ring
        room = max(1, int(2 * pi * radius / NODE_PITCH))
        placing = min(room, count - len(offsets))
        step = 2 * pi / placing
        turn = step / 2 if ring % 2 == 0 else 0.0
        offsets += [
            (radius * cos(index * step + turn), radius * sin(index * step + turn))
            for index in range(placing)
        ]
        ring += 1
    return offsets


def position_of(elements: list[dict], node_id: str) -> dict | None:
    """The canvas position recorded for a node, if it has one"""
    for element in elements:
        if element["data"]["id"] == node_id:
            return element.get("position")
    return None


def anchor_position(elements: list[dict], tap_node: dict | None, node_id: str) -> dict:
    """
    Where to hang an expansion's new nodes.

    tapNode carries the live position of the node the user just clicked;
    `elements` only catches up when dash-cytoscape pushes the canvas back, which
    it does on add, remove and drag but not after a layout run - so it can be
    one arrangement stale. Both beat the origin, which is where cytoscape drops
    a node that arrives without a position.
    """
    if tap_node and (tap_node.get("data") or {}).get("id") == node_id:
        live = tap_node.get("position")
        if live:
            return live
    return position_of(elements, node_id) or {"x": 0.0, "y": 0.0}


def seed_positions(
    new_nodes: list[dict], elements: list[dict], anchor: dict
) -> list[dict]:
    """
    `new_nodes` with a position on each one not already on the canvas.

    Only for a held layout, where nothing is going to place these. Cytoscape
    drops a positionless node at the origin, so without this a whole expansion
    lands in one pile in the corner rather than around what it came from.

    Nodes already drawn come back untouched: merge_elements() carries their
    existing position over, and moving them is the thing holding the layout is
    meant to prevent.
    """
    drawn = {element["data"]["id"] for element in elements}
    fresh = [element for element in new_nodes if element["data"]["id"] not in drawn]
    offsets = dict(
        zip(
            (element["data"]["id"] for element in fresh),
            ring_offsets(len(fresh)),
        )
    )

    placed = []
    for element in new_nodes:
        offset = offsets.get(element["data"]["id"])
        if offset is None:
            placed.append(element)
            continue
        position = {"x": anchor["x"] + offset[0], "y": anchor["y"] + offset[1]}
        placed.append({**element, "position": position})
    return placed


def type_group(api_name: str) -> str | None:
    """
    The legend bucket a type falls in - Tissue and CellLine are both Biosample
    - or None for one the legend does not map.

    Goes through NodeColor rather than ABSTRACT_MAPPING directly, so a node's
    column and its colour can never disagree about what it is.
    """
    if not api_name:
        return None
    try:
        return NodeColor(api_name).name
    except ValueError:
        return None


def column_of(api_name: str | None) -> int:
    """Which column a node type belongs in, left to right"""
    return _COLUMN_OF.get(type_group(api_name or ""), UNPLACED_COLUMN)


def spread(node_ids: list[str]) -> dict[str, float]:
    """One column's nodes as heights, ROW_PITCH apart and centred on zero"""
    middle = (len(node_ids) - 1) / 2
    return {
        node_id: (index - middle) * ROW_PITCH for index, node_id in enumerate(node_ids)
    }


def barycentre(
    node_id: str,
    linked: dict[str, set[str]],
    heights: dict[str, float],
    column: set[str],
) -> float:
    """
    The average height of a node's neighbors in the *other* columns, or its own
    current height where it has none - an isolated node should stay put rather
    than migrate to the top.
    """
    neighbors = [
        heights[other]
        for other in linked.get(node_id, ())
        if other in heights and other not in column
    ]
    return sum(neighbors) / len(neighbors) if neighbors else heights[node_id]


def column_positions(elements: list[dict]) -> dict[str, dict]:
    """
    A position for every node, in vertical columns by type: the positions map
    behind the "columns by type" layout.

    Cytoscape has no layout that groups by an arbitrary attribute - dagre ranks
    by distance from a root, which puts a Tissue and a SequenceFile in the same
    column whenever the path lengths happen to match - so this is a `preset`
    computed here instead. x is the type's column, y the node's row in it.

    The rows are not alphabetical: after an initial sort by label, a couple of
    barycentre passes pull each node level with its neighbors in the columns
    either side. Without them a column of 30 files faces a column of donors in
    an unrelated order, and every edge in the graph crosses every other.

    Columns are packed, so a graph with no Libraries in it has no empty gutter
    where they would have gone.

    Placeholders take no part in the ordering - they have no edges to be pulled
    level with - and stack above the top of their column instead, in type
    order. A type with nothing drawn yet is a column of placeholder alone.
    """
    nodes = [
        element
        for element in elements
        if "source" not in element["data"] and not element["data"].get("is_group")
    ]
    headers: dict[int, list[str]] = defaultdict(list)
    for element in sorted(
        (element for element in elements if element["data"].get("is_group")),
        key=lambda element: element["data"]["node_type"],
    ):
        headers[column_of(element["data"]["node_type"])].append(element["data"]["id"])
    if not nodes and not headers:
        return {}

    columns: dict[int, list[str]] = defaultdict(list)
    for element in sorted(
        nodes,
        key=lambda element: (element["data"].get("label") or "", element["data"]["id"]),
    ):
        column = column_of(element["data"].get("node_type"))
        columns[column].append(element["data"]["id"])

    linked: dict[str, set[str]] = defaultdict(set)
    for element in elements:
        data = element["data"]
        if "source" in data:
            linked[data["source"]].add(data["target"])
            linked[data["target"]].add(data["source"])

    heights: dict[str, float] = {}
    for node_ids in columns.values():
        heights.update(spread(node_ids))

    for sweep in range(ORDERING_SWEEPS):
        # alternating direction, so the pass sees the columns it just moved
        for column in sorted(columns, reverse=bool(sweep % 2)):
            node_ids = columns[column]
            here = set(node_ids)
            node_ids.sort(
                key=lambda node_id: (
                    barycentre(node_id, linked, heights, here),
                    # the height it already had, so the sort is stable and the
                    # columns do not reshuffle between identical graphs
                    heights[node_id],
                )
            )
            heights.update(spread(node_ids))

    for column, header_ids in headers.items():
        top = min(
            (heights[node_id] for node_id in columns.get(column, ())), default=0.0
        )
        for rank, header_id in enumerate(header_ids):
            heights[header_id] = top - (len(header_ids) - rank) * ROW_PITCH

    everything = {**headers}
    for column, node_ids in columns.items():
        everything[column] = everything.get(column, []) + node_ids
    packed = {column: index for index, column in enumerate(sorted(everything))}
    return {
        node_id: {"x": float(packed[column] * COLUMN_GAP), "y": heights[node_id]}
        for column, node_ids in everything.items()
        for node_id in node_ids
    }


def place_expansion(
    new_nodes: list[dict], elements: list[dict], tap_node: dict | None, anchor_id: str
) -> list[dict]:
    """seed_positions() around the node an expansion was clicked on"""
    return seed_positions(
        new_nodes, elements, anchor_position(elements, tap_node, anchor_id)
    )


def drop_nodes(elements: list[dict], node_ids: Collection[str]) -> list[dict]:
    """
    Remove nodes and any edge touching one of them.

    Only the named nodes go. Anything that was reachable only through them is
    left in place as an isolated node rather than cascaded away - unticking one
    member of a group should not silently delete a subtree the user expanded
    from it.
    """
    doomed = set(node_ids)
    if not doomed:
        return elements
    return [
        element
        for element in elements
        if element["data"]["id"] not in doomed
        and doomed.isdisjoint(
            (element["data"].get("source"), element["data"].get("target"))
        )
    ]


def drop_node(elements: list[dict], node_id: str) -> list[dict]:
    """Remove a node and any edge touching it"""
    return drop_nodes(elements, [node_id])


def properties_of(uuid_path: str) -> dict:
    """
    Cached JSON for a node, minus the keys that are noise in a detail panel.
    Returns {} for a node known only as an unfetched stub.
    """
    skipped = {"@context", "@type", "audit", "actions", "schema_version", "uuid"}
    json_object = LatticeNode._cache.get(uuid_path) or {}
    return {
        key: value
        for key, value in sorted(json_object.items())
        if key not in skipped and value not in (None, [], {}, "")
    }
