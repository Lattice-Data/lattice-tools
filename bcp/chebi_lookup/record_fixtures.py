"""Record live PubChem JSON responses as committed test fixtures."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import urllib.parse
from pathlib import Path

import requests

from cas_registry import CAS_INVALID_FORMAT, CAS_MISSING, classify_cas

from .client import BASE, PROPERTIES, REQUEST_DELAY, get_with_retry

log = logging.getLogger(__name__)

# bcp/tests/fixtures/chebi_lookup/pubchem_live/{cas}/
FIXTURES_ROOT = (
    Path(__file__).resolve().parent.parent
    / "tests"
    / "fixtures"
    / "chebi_lookup"
    / "pubchem_live"
)


def _save_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    log.info("Wrote %s", path)


def record_fixtures_for_cas(cas: str, out_root: Path | None = None) -> int | None:
    """
    Fetch PubChem responses for one CAS and write fixture JSON files.

    Returns the resolved CID, or None if CAS did not resolve.

    The recorder validates the CAS the same way production does, then queries and
    names the fixture directory after the normalised value. An operator pasting a
    corrupted number from the sheet would otherwise write a fixture for a URL
    production would never issue, and the live test would pass against a request
    path the shipped code cannot reach.
    """
    queried, cas_class, repairs = classify_cas(cas)
    if cas_class in (CAS_MISSING, CAS_INVALID_FORMAT):
        log.error("Not a CAS Registry Number, not recording: %r (%s)", cas, cas_class)
        return None
    if repairs:
        log.warning(
            "CAS %r repaired to %r before recording (%s)", cas, queried, repairs
        )

    out_dir = (out_root or FIXTURES_ROOT) / queried

    resp = get_with_retry(
        f"{BASE}/compound/name/{urllib.parse.quote(queried, safe='')}/cids/JSON"
    )
    time.sleep(REQUEST_DELAY)
    if resp is None:
        log.error("No CID response for CAS %s", cas)
        return None
    cids_payload = resp.json()
    _save_json(out_dir / "cids.json", cids_payload)

    try:
        cids = cids_payload.get("IdentifierList", {}).get("CID", [])
        cid = cids[0] if cids else None
    except (ValueError, KeyError, AttributeError, TypeError, IndexError):
        cid = None
    if cid is None:
        log.error("Could not parse CID for CAS %s", cas)
        return None

    resp = get_with_retry(f"{BASE}/compound/cid/{cid}/property/{PROPERTIES}/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        log.error("No properties response for CID %s", cid)
        return None
    _save_json(out_dir / "properties.json", resp.json())

    reg_resp = get_with_retry(f"{BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON")
    time.sleep(REQUEST_DELAY)
    if reg_resp is None:
        log.error("No registry ID response for CID %s", cid)
        return None
    # Held in its own name: `resp` is rebound to the synonyms response below, and
    # the xref parse used to read RegistryID back out of that instead, so chebi_id
    # was always "" and the recorder logged "ChEBI --" for every compound.
    registry_payload = reg_resp.json()
    _save_json(out_dir / "registry_ids.json", registry_payload)

    resp = get_with_retry(f"{BASE}/compound/cid/{cid}/synonyms/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        log.error("No synonyms response for CID %s", cid)
        return None
    _save_json(out_dir / "synonyms.json", resp.json())

    chebi_id = ""
    try:
        reg_ids = (
            registry_payload.get("InformationList", {})
            .get("Information", [{}])[0]
            .get("RegistryID", [])
        )
        for rid in reg_ids:
            if rid.upper().startswith("CHEBI:"):
                chebi_id = rid
                break
    except (ValueError, KeyError, AttributeError, TypeError, IndexError) as exc:
        # Logged for the reason client.py gives at the same parse: an unusable xref
        # body would otherwise read as "this compound has no ChEBI cross-reference",
        # and the recorder would write a fixture asserting that.
        log.warning("Unusable registry-ID payload for CID %s: %s", cid, exc)

    log.info(
        "Recorded fixtures for CAS %s (CID %s, ChEBI %s) → %s",
        queried,
        cid,
        chebi_id or "—",
        out_dir,
    )
    return cid


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Record live PubChem JSON fixtures for chebi_lookup tests."
    )
    parser.add_argument(
        "--cas",
        default="64-17-5",
        help="CAS Registry Number to record (default: 64-17-5 ethanol).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=f"Output root (default: {FIXTURES_ROOT}).",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        cid = record_fixtures_for_cas(args.cas, args.out_dir)
    except requests.exceptions.RequestException as exc:
        log.error("Network error: %s", exc)
        sys.exit(1)

    if cid is None:
        sys.exit(1)


if __name__ == "__main__":
    main()
