"""Shared fakes and fixtures for graph_db2 tests.

Nothing here touches the network. `expand()` normally makes one authenticated
GET per node plus a batched report per neighbor type; tests replace the GET with
`fake_fetch_full` and the DB2Gatherer with `FakeGatherer`.
"""

from __future__ import annotations

from typing import Any

import pytest
import requests

from graph_db2 import cyto_elements, graphing
from graph_db2.connection import get_connection
from graph_db2.constants import DEFAULT_MODE
from graph_db2.models import LatticeNode
from graph_db2.schema import create_configs

TEST_MODE = "db2_test"
TEST_SERVER = "https://db2-test.example.org/"

MFS = "/matrix_file_sets/aaaa1111-2222-3333-4444-555555555555/"
TISSUE = "/tissues/cccc1111-2222-3333-4444-555555555555/"
LAB = "/labs/test-lab/"
USER = "/users/dddd1111-2222-3333-4444-555555555555/"

RAW_MATRIX_FILE_COUNT = 30
SEQUENCE_FILE_COUNT = 2


def raw_matrix_file(index: int) -> str:
    return f"/raw_matrix_files/bbbb{index:04d}-2222-3333-4444-555555555555/"


def sequence_file(index: int) -> str:
    return f"/sequence_files/eeee{index:04d}-2222-3333-4444-555555555555/"


def build_objects() -> dict[str, dict]:
    """
    A synthetic DB2 graph:

        MatrixFileSet -> 30 RawMatrixFile -> 1 Tissue + 2 SequenceFile

    plus a lab and a submitter, which live in EXCLUDED_SCHEMAS and must never
    appear as neighbors.
    """
    objects: dict[str, dict] = {
        MFS: {
            "@id": MFS,
            "@type": ["MatrixFileSet", "FileSet", "Item"],
            "uuid": MFS.split("/")[2],
            "aliases": ["test-lab:my_matrixfileset"],
            "status": "current",
            "schema_version": 2,
            "audit": {},
            "raw_matrix_files": [
                raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)
            ],
            "lab": {"@id": LAB, "title": "Test Lab"},
            "submitted_by": {"@id": USER, "title": "Test User"},
        },
        TISSUE: {
            "@id": TISSUE,
            "@type": ["Tissue", "Biosample", "Item"],
            "uuid": TISSUE.split("/")[2],
            "aliases": ["test-lab:my_tissue"],
            "preservation_method": "fixed-frozen",
        },
        LAB: {"@id": LAB, "title": "Test Lab"},
        USER: {"@id": USER, "title": "Test User"},
    }

    for index in range(RAW_MATRIX_FILE_COUNT):
        path = raw_matrix_file(index)
        objects[path] = {
            "@id": path,
            "@type": ["RawMatrixFile", "File", "Item"],
            "uuid": path.split("/")[2],
            "aliases": [f"test-lab:matrix_{index}.h5"],
            "file_format": "h5",
            "matrix_file_sets": [MFS],
            "derived_from": [
                sequence_file(index * SEQUENCE_FILE_COUNT + offset)
                for offset in range(SEQUENCE_FILE_COUNT)
            ],
            "sample": TISSUE,
            "lab": {"@id": LAB, "title": "Test Lab"},
        }

    for index in range(RAW_MATRIX_FILE_COUNT * SEQUENCE_FILE_COUNT):
        path = sequence_file(index)
        objects[path] = {
            "@id": path,
            "@type": ["SequenceFile", "File", "Item"],
            "uuid": path.split("/")[2],
            # deliberately no aliases: label_for must fall back to accession
            "accession": f"LATSQ{index:06d}",
            "file_format": "fastq",
        }

    return objects


class FakeGatherer:
    """
    Stands in for DB2Gatherer.

    Only `chunk_and_fetch` is used, via models.batch_request_chunk_and_fetch.
    `unconfigured` mirrors a type with no OBJECT_CONFIG entry, which the real
    gatherer answers with a warning and an empty list.
    """

    def __init__(
        self, objects: dict[str, dict], unconfigured: set[str] | None = None
    ) -> None:
        self._objects = objects
        self.unconfigured = set(unconfigured or ())
        self.calls: list[tuple[str, frozenset[str]]] = []

    @property
    def fetched_uuids(self) -> set[str]:
        return {uuid for _, uuids in self.calls for uuid in uuids}

    def chunk_and_fetch(
        self, obj_type: str, object_ids: Any, **_kwargs: Any
    ) -> list[dict]:
        uuids = {str(uuid) for uuid in object_ids}
        self.calls.append((obj_type, frozenset(uuids)))
        if obj_type in self.unconfigured:
            return []

        return [
            json_object
            for path, json_object in self._objects.items()
            if path.split("/")[2] in uuids
            and LatticeNode(path).schema_ids.api_name == obj_type
        ]


def fake_fetch_full(objects: dict[str, dict]):
    """Replacement for cyto_elements.fetch_full: cache lookup, no HTTP."""

    def _fetch(node: LatticeNode) -> dict:
        if node.uuid_path not in objects:
            raise requests.HTTPError(
                f"404 Client Error: Not Found for url: {TEST_SERVER}{node.uuid_path}"
            )
        return objects[node.uuid_path]

    return _fetch


def node_ids(elements: list[dict]) -> set[str]:
    return {
        element["data"]["id"] for element in elements if "source" not in element["data"]
    }


def edges_of(elements: list[dict]) -> list[dict]:
    return [element for element in elements if "source" in element["data"]]


def dangling_edges(elements: list[dict]) -> list[str]:
    """
    Edge ids whose endpoints are not both present as nodes.

    Cytoscape refuses to render a graph containing one of these - it drops
    everything rather than the single bad edge - so this is the invariant that
    keeps the canvas from going blank.
    """
    ids = node_ids(elements)
    return [
        edge["data"]["id"]
        for edge in edges_of(elements)
        if not {edge["data"]["source"], edge["data"]["target"]} <= ids
    ]


def fake_gatherer(unconfigured: set[str] | None = None) -> FakeGatherer:
    """
    A gatherer over a fresh copy of the synthetic graph.

    A plain factory rather than a fixture: a fixture named `gatherer` would be
    shadowed by every test parameter of the same name, which ruff flags as F811.
    """
    return FakeGatherer(build_objects(), unconfigured)


@pytest.fixture(autouse=True)
def patched_fetch(monkeypatch: pytest.MonkeyPatch) -> None:
    """Swap the authenticated GET in expand() for a cache lookup."""
    monkeypatch.setattr(cyto_elements, "fetch_full", fake_fetch_full(build_objects()))


@pytest.fixture(autouse=True)
def db2_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Credentials for TEST_MODE.

    db2_flattener's Connection calls sys.exit when any of the three is missing,
    so anything reaching get_connection() needs these set.
    """
    monkeypatch.setenv("DB2_TEST_KEY", "test-key")
    monkeypatch.setenv("DB2_TEST_SECRET", "test-secret")
    monkeypatch.setenv("DB2_TEST_SERVER", TEST_SERVER)


@pytest.fixture(autouse=True)
def clean_graph_state():
    """
    graph_db2 keeps process-wide state: LatticeNode._cache and .mode are class
    attributes, cyto_elements._fully_fetched and graphing.graph_dict are module
    globals, and get_connection/create_configs are @cache'd. Leaking any of it
    between tests makes results order-dependent.
    """

    def _reset() -> None:
        LatticeNode._cache.clear()
        LatticeNode.mode = DEFAULT_MODE
        cyto_elements._fully_fetched.clear()
        graphing.graph_dict.clear()
        get_connection.cache_clear()
        create_configs.cache_clear()

    _reset()
    yield
    _reset()
