"""Shared fakes and fixtures for graph_db2 tests.

Nothing here touches the network. `expand()` normally makes one authenticated
GET per node plus a batched report per neighbor type; tests replace the GET with
`fake_fetch_full` and the DB2Gatherer with `FakeGatherer`.
"""

from __future__ import annotations

import json
from typing import Any
from urllib.parse import unquote, urlparse

import pytest
import requests

from graph_db2 import cyto_elements, graphing
from graph_db2.connection import get_connection
from graph_db2.constants import DEFAULT_MODE
from graph_db2.models import LatticeNode
from graph_db2.schema import create_configs

TEST_MODE = "db2_test"
TEST_SERVER = "https://db2-test.example.org/"
# expand() and normalize_to_type_and_uuid() default to DEFAULT_MODE, so tests
# that do not pass a mode still build a Connection; point it at a host that
# cannot answer, alongside the patched requests.get, so nothing can escape
DEFAULT_SERVER = "https://db2-default.example.org/"

MFS = "/matrix_file_sets/aaaa1111-2222-3333-4444-555555555555/"
TISSUE = "/tissues/cccc1111-2222-3333-4444-555555555555/"
LAB = "/labs/test-lab/"
USER = "/users/dddd1111-2222-3333-4444-555555555555/"

# the three spellings of MFS that --seed accepts
MFS_UUID = MFS.split("/")[2]
MFS_ALIAS = "test-lab:my_matrixfileset"
# an alias whose punctuation has to survive being put in a URL: unencoded, '#'
# truncates the request path and '?' turns the rest into a query string
MFS_AWKWARD_ALIAS = "test-lab:mfs#1?v=2"

RAW_MATRIX_FILE_COUNT = 30
SEQUENCE_FILE_COUNT = 2

# collection endpoints answer 200 with an '@id' that names no single object
COLLECTIONS = ("matrix_file_sets", "raw_matrix_files", "sequence_files", "tissues")


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
            "uuid": MFS_UUID,
            "aliases": [MFS_ALIAS, MFS_AWKWARD_ALIAS],
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


class FakeResponse:
    """The three bits of requests.Response that graph_db2 touches."""

    def __init__(self, status_code: int, payload: Any, url: str = "") -> None:
        self.status_code = status_code
        self.url = url
        self._payload = payload

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(
                f"{self.status_code} Client Error: for url: {self.url}", response=self
            )

    def json(self) -> Any:
        if isinstance(self._payload, Exception):
            raise self._payload
        return self._payload


def fake_requests_get(objects: dict[str, dict]):
    """
    Stand-in for requests.get against a DB2 instance.

    Answers the three seed spellings the API resolves - '/{object_type}/{uuid}/',
    '/{alias}/' and '/{uuid}/' - plus the two 200s that name no single object: a
    collection and the portal root.

    Records the decoded path of every request on `.paths` and the URL exactly as
    built on `.urls`, so a test can assert what was asked for, how often, on
    which server, and how it was encoded.
    """
    by_uuid = {
        json_object["uuid"]: json_object
        for json_object in objects.values()
        if "uuid" in json_object
    }
    by_alias = {
        alias: json_object
        for json_object in objects.values()
        for alias in json_object.get("aliases") or []
    }
    paths: list[str] = []
    urls: list[str] = []

    def _get(url: str, auth: Any = None, **_kwargs: Any) -> FakeResponse:
        # what reaches the server, not what the caller meant: percent-encoding
        # bugs show up here or nowhere
        urls.append(url)
        path = unquote(urlparse(url).path)
        paths.append(path)
        if auth is None:
            # an unauthenticated GET only sees released objects, which is how a
            # neighbor set comes back short rather than erroring
            raise AssertionError(f"unauthenticated request for {path}")

        segments = [segment for segment in path.split("/") if segment]
        if not segments:
            return FakeResponse(200, {"@type": ["Portal"]}, url)
        if len(segments) == 1 and segments[0] in COLLECTIONS:
            return FakeResponse(200, {"@id": f"/{segments[0]}/", "@graph": []}, url)
        if len(segments) == 1:
            found = by_uuid.get(segments[0]) or by_alias.get(segments[0])
        elif len(segments) == 2:
            found = objects.get(f"/{segments[0]}/{segments[1]}/")
        else:
            found = None

        if found is None:
            return FakeResponse(404, {"code": 404, "title": "Not Found"}, url)
        return FakeResponse(200, found, url)

    _get.paths = paths
    _get.urls = urls
    return _get


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
def patched_requests(monkeypatch: pytest.MonkeyPatch, db2_env: None):
    """
    Swap requests.get for the fake DB2 endpoint.

    normalize_to_type_and_uuid() resolves a seed by requesting it, so seed
    handling cannot be tested without either a server or this. Patching the
    module attribute rather than cyto_elements' reference also means a call from
    anywhere in graph_db2 is caught instead of quietly leaving the machine.

    Takes db2_env because resolving needs a Connection, and a module importing
    this fixture should not also have to remember the credentials it implies.
    """
    monkeypatch.setattr(requests, "get", fake_requests_get(build_objects()))


def requested_paths() -> list[str]:
    """
    The decoded paths the patched requests.get has been asked for, in order.

    A plain accessor rather than a fixture returning the fake: a fixture named
    `patched_requests` would be shadowed by a test parameter of the same name,
    which ruff flags as F811 - the same reason fake_gatherer() is a factory.
    """
    return requests.get.paths


def requested_urls() -> list[str]:
    """The URLs the patched requests.get was built with, unencoded as sent."""
    return requests.get.urls


@pytest.fixture(autouse=True)
def db2_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Credentials for TEST_MODE and for DEFAULT_MODE.

    db2_flattener's Connection calls sys.exit when any of the three is missing,
    so anything reaching get_connection() needs these set. DEFAULT_MODE is in
    here because expand() and normalize_to_type_and_uuid() default to it, and a
    test that does not name a mode should still not read the developer's real
    prod credentials out of the environment.
    """
    monkeypatch.setenv("DB2_TEST_KEY", "test-key")
    monkeypatch.setenv("DB2_TEST_SECRET", "test-secret")
    monkeypatch.setenv("DB2_TEST_SERVER", TEST_SERVER)
    for suffix, value in (
        ("KEY", "default-key"),
        ("SECRET", "default-secret"),
        ("SERVER", DEFAULT_SERVER),
    ):
        monkeypatch.setenv(f"{DEFAULT_MODE.upper()}_{suffix}", value)


# the callback that owns Load, keyed by its outputs
GROW_GRAPH = "..graph.elements...status.children...layout-choice.value.."
GROW_GRAPH_OUTPUTS = [
    {"id": "graph", "property": "elements"},
    {"id": "status", "property": "children"},
    {"id": "layout-choice", "property": "value"},
]


def built_app(seed: str, mode: str = TEST_MODE):
    """
    build_app() with the gatherer faked out.

    create_configs() looks the mode's server up in db2_flattener's
    constants.yaml, which has no entry for a test host, so make_gatherer has to
    go before build_app can run at all. Everything else - seed resolution
    included - is the real thing.
    """
    from graph_db2 import explorer

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(
            explorer,
            "make_gatherer",
            lambda mode, fetch_new: FakeGatherer(build_objects()),
        )
        return explorer.build_app(seed, mode, False)


def press_load(app, seed_value: str, fan: int = 500, elements: list | None = None):
    """
    Drive the Load button, as the browser would.

    Dash callbacks are closures inside build_app(), reachable only through
    callback_map. The wrapper builds its own dash.ctx from the callback_context
    kwarg, so `triggered_inputs` is how a test says which Input fired.
    """
    from dash._utils import AttributeDict

    callback = app.callback_map[GROW_GRAPH]["callback"]
    response = callback(
        1,
        None,
        seed_value,
        fan,
        elements if elements is not None else [],
        outputs_list=GROW_GRAPH_OUTPUTS,
        callback_context=AttributeDict(
            {
                "triggered_inputs": [{"prop_id": "load.n_clicks", "value": 1}],
                # the wrapper writes into this on the way out
                "updated_props": {},
            }
        ),
    )
    return json.loads(response)["response"]


def status_of(response: dict) -> str:
    """The status line's text out of a serialized callback response."""
    return str(response["status"]["children"]["props"]["children"])


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
