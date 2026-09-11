"""Tests for graph_db2.models: LatticeNode, NodeColor, batch-request helpers."""

from __future__ import annotations

import pytest

from graph_db2.constants import EXCLUDED_SCHEMAS
from graph_db2.models import (
    LatticeNode,
    NodeColor,
    batch_request_chunk_and_fetch,
    fetch_all_references_from_cache,
    group_batch_request,
)

from tests.graph_db2_helpers import (  # noqa: F401  (clean_graph_state is autouse)
    LAB,
    MFS,
    RAW_MATRIX_FILE_COUNT,
    SEQUENCE_FILE_COUNT,
    TISSUE,
    USER,
    FakeGatherer,
    build_objects,
    fake_gatherer,
    clean_graph_state,
    raw_matrix_file,
    sequence_file,
)


def test_parses_path_into_parts() -> None:
    node = LatticeNode(MFS)
    assert node.uuid_path == MFS
    assert node.uuid == MFS.split("/")[2]
    assert node.node_type == "matrix_file_sets"
    assert node.schema_ids.api_name == "MatrixFileSet"


@pytest.mark.parametrize(
    "raw",
    [
        "/matrix_file_sets/abc/",
        "matrix_file_sets/abc",
        "matrix_file_sets/abc/",
        "/matrix_file_sets/abc",
    ],
)
def test_uuid_path_is_canonical_regardless_of_input(raw: str) -> None:
    assert LatticeNode(raw).uuid_path == "/matrix_file_sets/abc/"


def test_mode_is_a_classvar_not_a_field() -> None:
    """build_app sets LatticeNode.mode once for the process; passing mode= to
    the constructor is a TypeError, and the README depends on that."""
    with pytest.raises(TypeError):
        LatticeNode("/tissues/x/", mode="db2_demo")  # type: ignore[call-arg]


def test_object_json_reads_through_cache_without_fetching() -> None:
    node = LatticeNode(TISSUE)
    node.object_json = {"@id": TISSUE, "aliases": ["lab:thing"]}
    assert LatticeNode(TISSUE).object_json["aliases"] == ["lab:thing"]


def test_cache_is_shared_across_instances() -> None:
    LatticeNode(TISSUE).object_json = {"@id": TISSUE}
    assert TISSUE in LatticeNode._cache


def test_get_ids_finds_nested_references() -> None:
    objects = build_objects()
    found = set(LatticeNode.get_ids(objects[raw_matrix_file(0)], EXCLUDED_SCHEMAS))
    assert TISSUE in found
    assert MFS in found
    assert sequence_file(0) in found


def test_get_ids_skips_excluded_schemas() -> None:
    """labs and users are in EXCLUDED_SCHEMAS, which is why a MatrixFileSet
    naming a lab and a submitter can still show only its files."""
    objects = build_objects()
    found = set(LatticeNode.get_ids(objects[MFS], EXCLUDED_SCHEMAS))
    assert LAB not in found
    assert USER not in found


def test_get_ids_ignores_terms_sentinel() -> None:
    payload = {"@context": "/terms/", "sample": "/tissues/x/"}
    assert set(LatticeNode.get_ids(payload, EXCLUDED_SCHEMAS)) == {"/tissues/x/"}


def test_get_ids_ignores_non_path_strings() -> None:
    payload = {"a": "not/a/path", "b": "/no-trailing-slash", "c": "plain"}
    assert set(LatticeNode.get_ids(payload, EXCLUDED_SCHEMAS)) == set()


def test_neighbors_excludes_self() -> None:
    objects = build_objects()
    node = LatticeNode(raw_matrix_file(0))
    node.object_json = objects[raw_matrix_file(0)]
    assert node.uuid_path not in node.neighbors


def test_neighbors_of_matrix_file_set() -> None:
    objects = build_objects()
    node = LatticeNode(MFS)
    node.object_json = objects[MFS]
    assert node.neighbors == {
        raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)
    }


def test_batch_cache_update_keys_on_at_id() -> None:
    LatticeNode.batch_cache_update([{"@id": TISSUE, "status": "current"}])
    assert LatticeNode._cache[TISSUE]["status"] == "current"


def test_batch_cache_update_skips_entries_without_at_id() -> None:
    LatticeNode.batch_cache_update([{"status": "current"}, {"@id": TISSUE}])
    assert list(LatticeNode._cache) == [TISSUE]


def test_group_batch_request_groups_by_api_name() -> None:
    grouped = group_batch_request([raw_matrix_file(0), sequence_file(0), TISSUE])
    assert set(grouped) == {"RawMatrixFile", "SequenceFile", "Tissue"}
    assert grouped["RawMatrixFile"] == {raw_matrix_file(0).split("/")[2]}


def test_group_batch_request_skips_cached_paths() -> None:
    """Already-cached objects are the reason a repeat expansion is nearly free."""
    LatticeNode(TISSUE).object_json = {"@id": TISSUE}
    assert group_batch_request([TISSUE]) == {}


def test_group_batch_request_empty_input() -> None:
    assert group_batch_request([]) == {}


def test_batch_request_chunk_and_fetch_one_call_per_type() -> None:
    gatherer = fake_gatherer()
    responses = batch_request_chunk_and_fetch(
        group_batch_request([raw_matrix_file(0), sequence_file(0)]), gatherer
    )
    assert {call[0] for call in gatherer.calls} == {"RawMatrixFile", "SequenceFile"}
    assert len(gatherer.calls) == 2
    assert {item["@id"] for item in responses} == {raw_matrix_file(0), sequence_file(0)}


def test_fetch_all_references_from_cache_walks_to_completion() -> None:
    """The pyvis path batch-loads everything reachable before walking, which is
    what makes build_adjacency_dict cheap."""
    objects = build_objects()
    gatherer = FakeGatherer(objects)
    LatticeNode(MFS).object_json = objects[MFS]
    fetch_all_references_from_cache(LatticeNode(MFS), gatherer)

    expected = (
        1  # the seed
        + RAW_MATRIX_FILE_COUNT
        + RAW_MATRIX_FILE_COUNT * SEQUENCE_FILE_COUNT
        + 1  # tissue
    )
    assert len(LatticeNode._cache) == expected


def test_fetch_all_references_from_cache_noop_on_empty_cache() -> None:
    gatherer = fake_gatherer()
    fetch_all_references_from_cache(LatticeNode(MFS), gatherer)
    assert gatherer.calls == []


def test_node_color_concrete_member() -> None:
    assert NodeColor("SequenceFile") is NodeColor.SequenceFile


@pytest.mark.parametrize(
    ("api_name", "expected"),
    [
        ("Tissue", NodeColor.Biosample),
        ("CellLine", NodeColor.Biosample),
        ("Organoid", NodeColor.Biosample),
        ("HumanDonor", NodeColor.Donor),
        ("NonHumanDonor", NodeColor.Donor),
        ("DropletBasedLibrary", NodeColor.Library),
        ("PlateBasedLibrary", NodeColor.Library),
    ],
)
def test_node_color_maps_concrete_to_abstract(api_name: str, expected) -> None:
    assert NodeColor(api_name) is expected


@pytest.mark.parametrize("api_name", ["Lab", "User", "Publication", "AnalysisStep"])
def test_node_color_raises_for_unmapped(api_name: str) -> None:
    with pytest.raises(ValueError):
        NodeColor(api_name)
