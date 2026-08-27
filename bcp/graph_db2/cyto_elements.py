"""
Cytoscape element construction and one-hop expansion.

graphing.py materializes the whole reachable graph up front so pyvis can write
it into a static file. This module walks outward one node at a time instead:
explorer.py asks for a node's neighbors only when a user clicks it, so the
payload stays small no matter how large the reachable graph is.

Only genuinely wide fans are grouped. A fan of DRAW_BUDGET or fewer is drawn in
full, so a 64-file MatrixFileSet looks like its 64 files. Past that, any type
with more than FAN_THRESHOLD members collapses into one placeholder the user can
search or fan out - a 500-wide fan is both unreadable and 500 label fetches.

Edges are undirected. LatticeNode.get_ids() walks the JSON on both sides of a
reference, so a Library and its FileSet each name the other; picking a direction
from discovery order would draw arrows that don't mean lineage.
"""

from collections import defaultdict
from collections.abc import Collection
from urllib.parse import urljoin

import requests
from .connection import get_connection
from .constants import DEFAULT_MODE, FETCH_NEW
from .models import (
    LatticeNode,
    NodeColor,
    batch_request_chunk_and_fetch,
    group_batch_request,
)
from .schema import create_configs

from db2_flattener.gather.gatherer import DB2Gatherer

UNMAPPED_COLOR = "#e8e8e8"
# neighbors drawn individually before any grouping kicks in
DRAW_BUDGET = 25
# once over budget, a type this size or smaller is still drawn rather than collapsed
FAN_THRESHOLD = 25
GROUP_SEPARATOR = "::group:"

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
    Canonicalize an object path to the '/object_type/uuid/' form.

    '/tissues/<uuid>/', 'tissues/<uuid>' and 'tissues/<uuid>/' all name the same
    object. LatticeNode normalizes internally, so mixing a raw input string with
    a LatticeNode.uuid_path silently produces edges whose endpoints match no
    node id - and cytoscape drops the whole graph rather than one bad edge.
    """
    cleaned = path.strip().strip("/")
    parts = cleaned.split("/")
    if len(parts) != 2 or not all(parts):
        raise ValueError(f"Expected 'object_type/uuid', got {path!r}")
    return f"/{cleaned}/"


def object_url(uuid_path: str, mode: str = DEFAULT_MODE) -> str:
    return urljoin(get_connection(mode).server, uuid_path)


def is_expanded(uuid_path: str) -> bool:
    return uuid_path in _fully_fetched


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


def node_element(node: LatticeNode) -> dict:
    return {
        "data": {
            "id": node.uuid_path,
            "label": label_for(node),
            "node_type": node.schema_ids.api_name,
            "color": color_for(node),
            "expanded": is_expanded(node.uuid_path),
        }
    }


def edge_element(one_path: str, other_path: str) -> dict:
    """Undirected edge with an order-independent id, so A-B and B-A dedupe"""
    source, target = sorted((one_path, other_path))
    return {"data": {"id": f"{source}--{target}", "source": source, "target": target}}


def group_id(parent_path: str, api_name: str) -> str:
    return f"{parent_path}{GROUP_SEPARATOR}{api_name}"


def group_element(
    parent_path: str, api_name: str, members: list[str], mode: str
) -> dict:
    """Placeholder standing in for a fan too wide to draw"""
    return {
        "data": {
            "id": group_id(parent_path, api_name),
            "label": f"{api_name} × {len(members)}",
            "node_type": api_name,
            "color": color_for(LatticeNode(members[0])),
            "expanded": False,
            "is_group": True,
            "parent_path": parent_path,
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
) -> tuple[list[dict], list[dict]]:
    """
    Resolve one node's neighbors into cytoscape elements.

    A fan of draw_budget or fewer is drawn in full - a 64-file MatrixFileSet is
    exactly what someone opening it wants to see. Only past that do the big
    types collapse into placeholders, since a 512-wide fan lays out ~23,000px
    tall and is unreadable at any zoom.

    Costs one authenticated GET for the node itself plus one batched report per
    type actually drawn. Grouped types cost nothing until the user fans them
    out, since a neighbor's type is readable from its path alone.
    """
    # everything downstream keys off node ids, so work from the canonical form
    # rather than whatever the caller typed
    uuid_path = normalize_path(uuid_path)
    node = LatticeNode(uuid_path)
    if node.uuid_path not in _fully_fetched:
        node.object_json = fetch_full(node)  # setter writes through to _cache
        _fully_fetched.add(node.uuid_path)

    neighbors = sorted(node.neighbors)
    by_type = defaultdict(list)
    for neighbor in neighbors:
        by_type[LatticeNode(neighbor).schema_ids.api_name].append(neighbor)

    drawn: list[str] = []
    grouped: dict[str, list[str]] = {}
    if draw_budget and len(neighbors) > draw_budget:
        for api_name, paths in by_type.items():
            if len(paths) > FAN_THRESHOLD:
                grouped[api_name] = paths
            else:
                drawn.extend(paths)
    else:
        drawn = neighbors

    fetch_labels(drawn, gatherer)

    nodes = [node_element(node)]
    edges = []
    for path in drawn:
        nodes.append(node_element(LatticeNode(path)))
        edges.append(edge_element(uuid_path, path))
    for api_name, paths in sorted(grouped.items()):
        nodes.append(group_element(uuid_path, api_name, paths, mode))
        edges.append(edge_element(uuid_path, group_id(uuid_path, api_name)))

    return nodes, edges


def promote_members(
    paths: Collection[str],
    parent_path: str,
    gatherer: DB2Gatherer,
    mode: str = DEFAULT_MODE,
) -> tuple[list[dict], list[dict]]:
    """
    Lift chosen members of a group onto the canvas.

    Takes the whole selection at once so a multi-pick costs one batched report
    per type rather than one per node.
    """
    fetch_labels(paths, gatherer)
    nodes = [node_element(LatticeNode(path)) for path in paths]
    edges = [edge_element(parent_path, path) for path in paths]
    return nodes, edges


def explode(
    group_data: dict, gatherer: DB2Gatherer, mode: str = DEFAULT_MODE
) -> tuple[list[dict], list[dict]]:
    """Turn a group placeholder into every one of its member nodes"""
    return promote_members(
        group_data["members"], group_data["parent_path"], gatherer, mode
    )


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
    What an expansion actually produced, counting group members rather than
    placeholders - "1 neighbor" for a collapsed fan of 64 reads as a failure.
    """
    drawn = sum(1 for node in nodes if not node["data"].get("is_group")) - 1
    groups = [node["data"] for node in nodes if node["data"].get("is_group")]

    parts = [f"{drawn} drawn"] if drawn else []
    parts += [f"{len(data['members'])} {data['node_type']} grouped" for data in groups]
    return ", ".join(parts) or "no neighbors"


def not_yet_drawn(elements: list[dict], paths: Collection[str]) -> list[str]:
    """The subset of `paths` not yet on the canvas, in selection order"""
    present = {element["data"]["id"] for element in elements}
    return [path for path in paths if path not in present]


def merge_elements(
    existing: list[dict], new_nodes: list[dict], new_edges: list[dict]
) -> list[dict]:
    """
    Fold an expansion into the elements already on screen.

    Nodes are replaced (a stub gains its real label once fetched), edges are
    kept on first sight.
    """
    by_id = {element["data"]["id"]: element for element in existing}

    for element in new_nodes:
        by_id[element["data"]["id"]] = element
    for element in new_edges:
        by_id.setdefault(element["data"]["id"], element)

    return list(by_id.values())


def drop_node(elements: list[dict], node_id: str) -> list[dict]:
    """Remove a node and any edge touching it"""
    return [
        element
        for element in elements
        if element["data"]["id"] != node_id
        and node_id
        not in (element["data"].get("source"), element["data"].get("target"))
    ]


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
