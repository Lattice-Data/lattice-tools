from collections import defaultdict
from functools import cache
from typing import Collection

from constants import DEFAULT_MODE
from db2_flattener.gather.lattice import Connection
from db2_flattener.gather.gatherer import DB2Gatherer
from models import BatchResponse, GroupedBatchRequest, LatticeNode


@cache
def get_connection(mode: str = DEFAULT_MODE) -> Connection:
    """One Connection per mode, built on first use"""
    return Connection(mode)


def configure_connection(mode: str = DEFAULT_MODE, **attributes) -> Connection:
    """Set attributes on connection if needed"""
    connection = get_connection(mode)
    for attribute, value in attributes.items():
        setattr(connection, attribute, value)
    return connection


def group_batch_request(uuid_paths: Collection) -> GroupedBatchRequest:
    """
    DB2Gatherer needs API name and UUID to make batch request
    Returns dict with key of API name and set of UUIDs to request
    """
    grouped_paths = defaultdict(set)
    for path in uuid_paths:
        # skip if already in cache
        if path in LatticeNode._cache:
            continue
        node = LatticeNode(path)
        grouped_paths[node.schema_ids.api_name].add(node.uuid)

    return grouped_paths


def batch_request_chunk_and_fetch(
    grouped_request: GroupedBatchRequest, gatherer: DB2Gatherer
) -> BatchResponse:
    """
    After grouping many objects by type, send off a batch request per type via DB2Gatherer.chunk_and_fetch()
    """
    responses = []
    for api_name, uuid_set in grouped_request.items():
        responses.extend(
            gatherer.chunk_and_fetch(
                obj_type=api_name,
                object_ids=uuid_set,
            )
        )
    return responses


def fetch_all_references_from_cache(node: LatticeNode, gatherer: DB2Gatherer):
    """
    After one API request stores a JSON response in the LatticeNode cache, can
    then run this function to batch request all references not excluded.
    Should call this prior to constructing adjacency dict
    """
    if not LatticeNode._cache:
        print("No items in cache, exiting")
        return

    running = True

    while running:
        cache_ids = set()
        for json_object in LatticeNode._cache.values():
            cache_ids.update(LatticeNode.get_ids(json_object, node.excluded))
        grouped_request = group_batch_request(cache_ids)
        if not grouped_request:
            running = False
            break
        responses = batch_request_chunk_and_fetch(grouped_request, gatherer)
        LatticeNode.batch_cache_update(responses)
