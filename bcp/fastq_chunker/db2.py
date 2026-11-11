"""Read-only DB2 client for ``resolve``: fetch objects by uuid through search.

Credentials follow the repo's ``--mode`` convention (db2_flattener, graph_db2):
``<MODE>_KEY``, ``<MODE>_SECRET`` and ``<MODE>_SERVER``. Nothing here writes.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable

from .plan import PlanError

# uuids per search request; keeps the query string well under URL limits
BATCH = 50
TIMEOUT_S = 120


class Client:
    def __init__(self, mode: str, session=None):
        prefix = mode.upper()
        values = {
            k: os.environ.get(f"{prefix}_{k}") for k in ("KEY", "SECRET", "SERVER")
        }
        missing = [f"{prefix}_{k}" for k, v in values.items() if not v]
        if missing:
            raise PlanError(f"DB2 credentials not set: {', '.join(missing)}")
        self.server = values["SERVER"].rstrip("/")
        if session is None:
            import requests

            session = requests.Session()
        session.auth = (values["KEY"], values["SECRET"])
        session.headers.update({"Accept": "application/json"})
        self.session = session

    def search(self, params: dict) -> list[dict]:
        r = self.session.get(f"{self.server}/search/", params=params, timeout=TIMEOUT_S)
        if r.status_code == 404:
            # snovault answers an empty search with 404 and an empty @graph
            return []
        r.raise_for_status()
        return list(r.json().get("@graph", []))


def fetch_by_uuid(
    search: Callable[[dict], list[dict]], item_type: str, uuids: Iterable[str]
) -> dict[str, dict]:
    """uuid -> object (frame=object) for every uuid the search returns.

    A uuid missing from the result was not found or is not visible to the
    search (deleted objects are hidden by default); the caller reports it.
    """
    wanted = list(dict.fromkeys(uuids))
    found: dict[str, dict] = {}
    for i in range(0, len(wanted), BATCH):
        params = {
            "type": item_type,
            "uuid": wanted[i : i + BATCH],
            "frame": "object",
            "format": "json",
            "limit": "all",
        }
        for obj in search(params):
            found[obj["uuid"]] = obj
    return found
