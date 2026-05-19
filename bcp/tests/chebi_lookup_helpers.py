"""Shared PubChem HTTP mocks for chebi_lookup tests."""

from __future__ import annotations

import json
import re
from pathlib import Path
from unittest.mock import MagicMock

FIXTURES = Path(__file__).parent / "fixtures" / "chebi_lookup"
PUBCHEM_LIVE = FIXTURES / "pubchem_live"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_json_name(name: str) -> dict:
    return load_json(FIXTURES / name)


def mock_response(status_code: int, payload: dict | None = None) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = payload or {}
    return resp


def _live_cid_for_cas(cas: str) -> int | None:
    cids_path = PUBCHEM_LIVE / cas / "cids.json"
    if not cids_path.exists():
        return None
    try:
        cids = load_json(cids_path).get("IdentifierList", {}).get("CID", [])
        return int(cids[0]) if cids else None
    except (TypeError, ValueError, KeyError, IndexError):
        return None


def _route_pubchem_live(url: str) -> MagicMock | None:
    """Return a mock response from recorded pubchem_live fixtures, or None."""
    name_match = re.search(r"/compound/name/([^/]+)/cids", url)
    if name_match:
        cas = name_match.group(1)
        cids_path = PUBCHEM_LIVE / cas / "cids.json"
        if cids_path.exists():
            return mock_response(200, load_json(cids_path))

    cid_match = re.search(r"/compound/cid/(\d+)/", url)
    if cid_match:
        cid = int(cid_match.group(1))
        for cas_dir in sorted(PUBCHEM_LIVE.iterdir()) if PUBCHEM_LIVE.exists() else []:
            if not cas_dir.is_dir():
                continue
            if _live_cid_for_cas(cas_dir.name) != cid:
                continue
            if "/property/" in url:
                path = cas_dir / "properties.json"
            elif "/xrefs/RegistryID/" in url:
                path = cas_dir / "registry_ids.json"
            elif "/synonyms/" in url:
                path = cas_dir / "synonyms.json"
            else:
                return None
            if path.exists():
                return mock_response(200, load_json(path))
    return None


def route_pubchem_get(url: str, *args, **kwargs):
    live = _route_pubchem_live(url)
    if live is not None:
        return live

    if url.endswith("/cids/JSON"):
        return mock_response(200, load_json_name("cid_list.json"))
    if "/property/" in url:
        return mock_response(200, load_json_name("properties.json"))
    if "/xrefs/RegistryID/JSON" in url:
        return mock_response(200, load_json_name("registry_ids_chebi.json"))
    if "/synonyms/JSON" in url:
        return mock_response(200, load_json_name("synonyms.json"))
    return mock_response(404)
