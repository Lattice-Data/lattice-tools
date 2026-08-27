from .models import GraphDict, LatticeNode, NodeColor

graph_dict: GraphDict = {}


def build_adjacency_dict(start_node: LatticeNode):
    """
    Builds adjacency dict for use to get edges.
    Should call connection.fetch_all_references_from_cache() prior to load
    all JSON objects into LatticeNode cache, otherwise will walk graph
    one request at a time and will take ~40X as long

    Currently using global dict outside of scope, can change to pass
    as arg for the recursive function calls
    """
    graph_dict[start_node.uuid_path] = start_node

    if not start_node.neighbors:
        return

    for neighbor in start_node.neighbors:
        node = LatticeNode(neighbor)
        if node.uuid_path in graph_dict:
            continue
        graph_dict[node.uuid_path] = node
        build_adjacency_dict(node)


def find_edges(graph_dict: GraphDict) -> list[tuple[str, str]]:
    """
    Iterate through all neighbors in nodes to create tuple pairings
    of edges
    """
    result = []
    for lattice_id, node in graph_dict.items():
        for neighbor in node.neighbors:
            result.append((lattice_id, neighbor))

    return result


def create_nodes_for_pyvis(graph_dict: GraphDict) -> list[dict]:
    """
    Create input for pyvis visual graphing
    """
    result = []
    for uuid_path, node in graph_dict.items():
        result.append(
            {
                "n_id": uuid_path,
                "label": node.alias if node.alias is not None else uuid_path,
                "color": NodeColor(node.schema_ids.api_name).value,
                "title": "Neighbors: <br>"
                + "<br>".join(
                    [graph_dict[neighbor].alias for neighbor in node.neighbors]
                ),
            }
        )
    return result
