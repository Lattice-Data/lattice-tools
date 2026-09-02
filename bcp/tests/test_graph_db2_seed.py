"""
Tests for seed resolution - turning any accepted --seed spelling into a node id.

A DB2 instance resolves three spellings of the same object: '/type/uuid/', a
bare alias and a bare uuid. Only the first is usable as a cytoscape node id, and
only the server can say which object the other two name, so
normalize_to_type_and_uuid() requests the path and reads '@id' off the answer.

Nothing here touches the network: the autouse `patched_requests` fixture swaps
requests.get for a fake DB2 endpoint over the same synthetic graph the rest of
the suite uses, and requested_paths()/requested_urls() report what it was asked
for.
"""

from __future__ import annotations

import pytest
import requests

from graph_db2.cyto_elements import (
    expand,
    group_id,
    merge_elements,
    normalize_to_type_and_uuid,
    resolve_seed,
)
from graph_db2.models import GraphDB2Error, LatticeNode

from tests.graph_db2_helpers import (  # noqa: F401  (fixtures + autouse reset)
    DEFAULT_SERVER,
    MFS,
    MFS_ALIAS,
    MFS_AWKWARD_ALIAS,
    MFS_UUID,
    RAW_MATRIX_FILE_COUNT,
    TEST_MODE,
    TEST_SERVER,
    FakeResponse,
    built_app,
    clean_graph_state,
    dangling_edges,
    db2_env,
    fake_gatherer,
    node_ids,
    patched_fetch,
    patched_requests,
    press_load,
    raw_matrix_file,
    requested_paths,
    requested_urls,
    status_of,
)

BIG_BUDGET = 500

# what every caller runs; aliased for brevity here
resolve = resolve_seed


# --------------------------------------------------------------------------
# the three accepted spellings all land on the same node id
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        MFS,
        MFS.strip("/"),
        MFS_ALIAS,
        f"/{MFS_ALIAS}/",
        f"  {MFS_ALIAS}  ",
        MFS_UUID,
        f"/{MFS_UUID}/",
        f"\t{MFS_UUID}\n",
    ],
)
def test_every_spelling_resolves_to_the_canonical_id(raw: str) -> None:
    assert resolve(raw) == MFS


@pytest.mark.parametrize("raw", [MFS, MFS_ALIAS, MFS_UUID])
def test_resolved_id_is_a_usable_node_id(raw: str) -> None:
    """The whole point of resolving: an alias or uuid cannot be a node id
    because LatticeNode cannot parse one, and an edge naming an id no node
    answers to makes cytoscape drop the entire graph."""
    assert LatticeNode(resolve(raw)).uuid_path == MFS


@pytest.mark.parametrize("raw", [MFS_ALIAS, MFS_UUID])
def test_alias_and_uuid_seeds_expand_like_a_path_seed(raw: str) -> None:
    expected = {MFS} | {
        raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)
    }
    nodes, edges = expand(raw, fake_gatherer(), draw_budget=BIG_BUDGET)
    elements = merge_elements([], nodes, edges)
    assert node_ids(elements) == expected
    assert dangling_edges(elements) == []


def test_alias_seed_is_marked_expanded_under_its_canonical_id() -> None:
    """A seed typed as an alias has to come back flagged as expanded, or the
    canvas offers a second, edge-free expansion of the node just drawn."""
    nodes, _ = expand(MFS_ALIAS, fake_gatherer(), draw_budget=BIG_BUDGET)
    by_id = {node["data"]["id"]: node["data"] for node in nodes}
    assert by_id[MFS]["expanded"] is True


# --------------------------------------------------------------------------
# what actually goes over the wire
# --------------------------------------------------------------------------


def test_resolution_requests_the_path_as_typed() -> None:
    assert resolve(MFS_ALIAS) == MFS
    assert requested_paths() == [f"/{MFS_ALIAS}/"]


def test_resolution_is_authenticated(monkeypatch: pytest.MonkeyPatch) -> None:
    """The fake endpoint raises on a missing auth tuple; a private object only
    resolves for a key that can see it, so an anonymous GET would 404 on
    exactly the objects a curator is most likely to seed from."""
    seen: dict = {}

    def _get(url, auth=None, headers=None, timeout=None):
        seen.update(url=url, auth=auth, headers=headers, timeout=timeout)
        return FakeResponse(200, {"@id": MFS}, url)

    monkeypatch.setattr(requests, "get", _get)
    assert resolve(MFS_ALIAS) == MFS
    assert seen["auth"] is not None
    assert seen["headers"]["accept"] == "application/json"
    # LatticeNode.fetch_object_json() passes one; without it a hung instance
    # hangs the Dash callback with it
    assert seen["timeout"]


def test_an_awkward_alias_is_percent_encoded() -> None:
    """Unencoded, the '#' truncates the request path to '/test-lab:mfs' and the
    '?' turns the rest into a query string - so the server is asked for a
    different object than the user typed, and answers 404."""
    assert resolve(MFS_AWKWARD_ALIAS) == MFS
    assert requested_paths() == [f"/{MFS_AWKWARD_ALIAS}/"]


def test_a_colon_is_sent_literally() -> None:
    """Every alias is 'lab-name:local-id', so a colon is the one character that
    is in all of them. It is legal unencoded in a path segment, and that is the
    form the API's own examples use."""
    resolve(MFS_ALIAS)
    assert requested_urls()[0].endswith(f"/{MFS_ALIAS}/")


def test_resolution_uses_the_server_for_the_given_mode() -> None:
    """A seed only means something on the server it came from. Resolving a
    db2_test seed against DEFAULT_MODE either 404s or - worse - finds a
    different object that happens to share the uuid."""
    resolve(MFS_ALIAS, TEST_MODE)
    assert requested_urls() == [f"{TEST_SERVER.rstrip('/')}/{MFS_ALIAS}/"]


def test_default_mode_is_only_used_when_asked_for() -> None:
    resolve(MFS_ALIAS)
    assert requested_urls()[0].startswith(DEFAULT_SERVER)


def test_expand_resolves_against_its_own_mode() -> None:
    """expand() resolves the seed itself, so its mode has to reach the
    resolver too - not just the batched label reports."""
    expand(MFS_ALIAS, fake_gatherer(), mode=TEST_MODE, draw_budget=BIG_BUDGET)
    assert requested_urls()[0].startswith(TEST_SERVER)


def test_resolution_failure_names_the_server() -> None:
    with pytest.raises(GraphDB2Error, match=DEFAULT_SERVER):
        resolve("no-such-lab:no-such-alias")


# --------------------------------------------------------------------------
# the cache short circuit
# --------------------------------------------------------------------------


def test_resolution_caches_the_profile_it_fetched() -> None:
    resolve(MFS_ALIAS)
    assert LatticeNode._cache[MFS]["aliases"] == [MFS_ALIAS, MFS_AWKWARD_ALIAS]


def test_a_canonical_id_in_the_cache_resolves_without_a_request() -> None:
    """Every key in the cache came from an '@id', so it is already canonical.
    Tapping a node calls expand() with an id that is always in the cache, and
    re-requesting it would put a round trip behind every click."""
    resolve(MFS_ALIAS)
    requested_paths().clear()

    assert normalize_to_type_and_uuid(MFS) == MFS
    assert requested_paths() == []


def test_a_second_alias_lookup_still_costs_a_request() -> None:
    """An alias is not a cache key - only '@id's are - so the short circuit
    cannot fire for one. Documented rather than fixed: caching alias -> '@id'
    would be a second mapping to invalidate when an alias is edited."""
    resolve(MFS_ALIAS)
    resolve(MFS_ALIAS)
    assert requested_paths() == [f"/{MFS_ALIAS}/"] * 2


def test_expanding_a_tapped_node_makes_no_resolution_request() -> None:
    """The regression this guards: the cache branch returned
    _cache['@id'] - a KeyError - so every tap after the first died in
    resolution and surfaced as 'Could not load ...'."""
    gatherer = fake_gatherer()
    expand(MFS, gatherer, draw_budget=BIG_BUDGET)
    requested_paths().clear()

    expand(raw_matrix_file(0), gatherer, draw_budget=BIG_BUDGET)
    assert requested_paths() == []


def test_a_failed_resolution_leaves_nothing_in_the_cache() -> None:
    with pytest.raises(GraphDB2Error):
        resolve("no-such-lab:no-such-alias")
    assert LatticeNode._cache == {}


# --------------------------------------------------------------------------
# 200s that name no object
# --------------------------------------------------------------------------


def test_a_collection_seed_is_rejected() -> None:
    """'matrix_file_sets' could have been an alias, so normalize_path() lets it
    through - but the collection endpoint answers 200 with '@id':
    '/matrix_file_sets/', which LatticeNode cannot split into two segments."""
    with pytest.raises(GraphDB2Error, match="not a single object"):
        resolve("matrix_file_sets")


def test_a_response_without_an_id_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    """Unguarded this returned None, cached the body under the None key, and
    then failed in LatticeNode with "'NoneType' object has no attribute
    'strip'" - a traceback that never mentions the seed."""
    monkeypatch.setattr(
        requests, "get", lambda url, **_: FakeResponse(200, {"@type": ["Portal"]}, url)
    )
    with pytest.raises(GraphDB2Error, match="not a single object"):
        resolve(MFS_ALIAS)
    assert LatticeNode._cache == {}


def test_a_non_json_200_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    """An HTML error page or a proxy login form answers 200; response.json()
    then raises, and it has to come back as the same GraphDB2Error a 404 does
    rather than as a bare JSONDecodeError."""
    broken = FakeResponse(200, requests.exceptions.JSONDecodeError("x", "x", 0))
    monkeypatch.setattr(requests, "get", lambda url, **_: broken)
    with pytest.raises(GraphDB2Error, match="Could not resolve seed"):
        resolve(MFS_ALIAS)


def test_an_unreachable_server_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    def _get(url, **_):
        raise requests.ConnectionError("nodename nor servname provided")

    monkeypatch.setattr(requests, "get", _get)
    with pytest.raises(GraphDB2Error, match="Could not resolve seed"):
        resolve(MFS_ALIAS)


@pytest.mark.parametrize(
    "raw", ["no-such-lab:no-such-alias", "ffffffff-2222-3333-4444-555555555555"]
)
def test_an_unknown_alias_or_uuid_is_rejected(raw: str) -> None:
    with pytest.raises(GraphDB2Error, match="Could not resolve seed"):
        resolve(raw)


def test_every_resolution_failure_is_catchable_as_one_type() -> None:
    """explorer.py catches ValueError around seed handling. Each of these used
    to escape as a different type - HTTPError, KeyError, AttributeError,
    JSONDecodeError - and only some of them were caught."""
    for raw in ("matrix_file_sets", "no-such-lab:no-such-alias", "   "):
        with pytest.raises(ValueError):
            resolve(raw)


# --------------------------------------------------------------------------
# the app: --seed at startup and the Load button
# --------------------------------------------------------------------------


def seed_box(app) -> str:
    """The value in the toolbar's seed Input."""
    for component in app.layout.children[0].children:
        if getattr(component, "id", None) == "seed":
            return component.value
    raise AssertionError("no seed input in the toolbar")


@pytest.mark.parametrize("raw", [MFS, MFS_ALIAS, MFS_UUID, MFS_AWKWARD_ALIAS])
def test_startup_seed_box_shows_the_canonical_id(raw: str) -> None:
    """An alias in the box would be resolved again on the next Load, and a
    curator reading it has no way to tell which object it turned out to be."""
    assert seed_box(built_app(raw)) == MFS


@pytest.mark.parametrize("raw", [MFS_ALIAS, MFS_UUID])
def test_startup_seed_draws_the_graph(raw: str) -> None:
    """Startup uses the default draw budget, so the 30-file fan arrives as one
    placeholder - but under the canonical id, which is what makes the group
    tappable at all."""
    app = built_app(raw)
    elements = app.layout.children[1].children[0].children.elements
    assert node_ids(elements) == {MFS, group_id(MFS, "RawMatrixFile")}
    assert dangling_edges(elements) == []


def test_startup_resolves_against_the_mode_server() -> None:
    """build_app() takes a mode; passing the seed on without it resolved every
    startup against DEFAULT_MODE, so --mode db2_demo --seed <demo alias> either
    404'd into an empty canvas or found a prod object with the same uuid."""
    built_app(MFS_ALIAS, TEST_MODE)
    assert requested_urls()[0].startswith(TEST_SERVER)


def test_an_unresolvable_startup_seed_says_so_in_the_panel() -> None:
    """An empty canvas looks the same as a silent failure, so the reason has to
    land somewhere the user is looking."""
    app = built_app("no-such-lab:no-such-alias", TEST_MODE)
    rendered = str(app.layout)
    assert "Nothing loaded" in rendered
    assert "no-such-lab:no-such-alias" in rendered


def test_an_unresolvable_startup_seed_does_not_raise() -> None:
    """Every resolution failure has to be one of the types build_app() catches;
    anything else takes the whole app down before it serves a page."""
    assert built_app("matrix_file_sets", TEST_MODE) is not None
    assert built_app("   ", TEST_MODE) is not None


@pytest.mark.parametrize("raw", [MFS, MFS_ALIAS, MFS_UUID, MFS_AWKWARD_ALIAS])
def test_load_accepts_every_spelling(raw: str) -> None:
    response = press_load(built_app("", TEST_MODE), raw)
    assert node_ids(response["graph"]["elements"]) == {MFS} | {
        raw_matrix_file(index) for index in range(RAW_MATRIX_FILE_COUNT)
    }


@pytest.mark.parametrize("raw", [MFS_ALIAS, MFS_UUID])
def test_load_status_names_the_object_the_seed_resolved_to(raw: str) -> None:
    """Typing an alias and being told only 'test-lab:my_matrixfileset - 30
    drawn' leaves the canonical path - the thing worth copying elsewhere -
    nowhere on screen."""
    status = status_of(press_load(built_app("", TEST_MODE), raw))
    assert status.startswith(MFS)


def test_load_resolves_against_the_mode_server() -> None:
    app = built_app("", TEST_MODE)
    requested_urls().clear()
    press_load(app, MFS_ALIAS)
    assert requested_urls()[0].startswith(TEST_SERVER)


def test_load_reports_an_unresolvable_seed_in_the_status_line() -> None:
    response = press_load(built_app("", TEST_MODE), "no-such-lab:no-such-alias")
    assert "no-such-lab:no-such-alias" in status_of(response)
    assert "graph" not in response  # the canvas is left alone


def test_load_reports_a_collection_seed_in_the_status_line() -> None:
    response = press_load(built_app("", TEST_MODE), "matrix_file_sets")
    assert "not a single object" in status_of(response)


def test_load_reports_a_full_url_in_the_status_line() -> None:
    response = press_load(built_app("", TEST_MODE), f"{TEST_SERVER}{MFS.lstrip('/')}")
    assert "not a URL" in status_of(response)
