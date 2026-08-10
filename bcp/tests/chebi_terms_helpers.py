"""Shared ChEBI HTTP mocks for chebi_terms tests."""

from __future__ import annotations

import json
import re
from copy import deepcopy
from pathlib import Path
from unittest.mock import MagicMock

FIXTURES = Path(__file__).parent / "fixtures" / "chebi_terms"
CHEBI_LIVE = FIXTURES / "chebi_live"

# Recorded fixtures, keyed by the *requested* numeric ID.
ETHANOL = "16236"
EPICATECHIN_SECONDARY = "18484"  # resolves to CHEBI:90
# Markup-heavy name *and* no CAS xref: the no_cas_recorded coverage rests on the
# second property, so do not swap in another markup-heavy ID that has a CAS.
MARKUP_NAME = "44567"


def load_payload(numeric_id: str) -> dict:
    """Load a recorded ChEBI payload by requested numeric ID."""
    path = CHEBI_LIVE / f"{numeric_id}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def payload_copy(numeric_id: str) -> dict:
    """A mutable deep copy, for tests that need to perturb a recorded payload."""
    return deepcopy(load_payload(numeric_id))


def mock_response(status_code: int, payload: dict | None = None) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = payload or {}
    return resp


def route_chebi_get(url: str, *args, **kwargs) -> MagicMock:
    """Serve recorded fixtures for /compound/{id}/, 404 for anything unrecorded."""
    match = re.search(r"/compound/(\d+)/", url)
    if match:
        path = CHEBI_LIVE / f"{match.group(1)}.json"
        if path.exists():
            return mock_response(200, json.loads(path.read_text(encoding="utf-8")))
    return mock_response(404)
