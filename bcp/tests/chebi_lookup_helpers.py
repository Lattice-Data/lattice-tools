"""Shared PubChem HTTP mocks for chebi_lookup tests."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

FIXTURES = Path(__file__).parent / "fixtures" / "chebi_lookup"


def load_json(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def mock_response(status_code: int, payload: dict | None = None) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = payload or {}
    return resp


def route_pubchem_get(url: str, *args, **kwargs):
    if url.endswith("/cids/JSON"):
        return mock_response(200, load_json("cid_list.json"))
    if "/property/" in url:
        return mock_response(200, load_json("properties.json"))
    if "/xrefs/RegistryID/JSON" in url:
        return mock_response(200, load_json("registry_ids_chebi.json"))
    if "/synonyms/JSON" in url:
        return mock_response(200, load_json("synonyms.json"))
    return mock_response(404)
