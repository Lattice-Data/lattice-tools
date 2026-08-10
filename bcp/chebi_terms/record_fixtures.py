"""Record live ChEBI JSON responses as committed test fixtures."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import requests

from .client import fetch_compound, normalize_chebi_id

log = logging.getLogger(__name__)

# bcp/tests/fixtures/chebi_terms/chebi_live/{numeric_id}.json
FIXTURES_ROOT = (
    Path(__file__).resolve().parent.parent
    / "tests"
    / "fixtures"
    / "chebi_terms"
    / "chebi_live"
)

# ethanol, a secondary ID that resolves to CHEBI:90, and a markup-heavy name.
DEFAULT_IDS = ["CHEBI:16236", "CHEBI:18484", "CHEBI:44567"]


def _save_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log.info("Wrote %s", path)


def record_fixture_for_id(chebi_id: str, out_root: Path | None = None) -> str | None:
    """
    Fetch the ChEBI payload for one ID and write it as a fixture file.

    Fixtures are keyed by the *requested* numeric ID so that secondary-ID
    behaviour stays reproducible offline.
    """
    parsed = normalize_chebi_id(chebi_id)
    if parsed is None:
        log.error("Not a ChEBI identifier: %s", chebi_id)
        return None
    numeric, accession = parsed

    payload = fetch_compound(numeric)
    if payload is None:
        log.error("No ChEBI record for %s", accession)
        return None

    out_dir = out_root or FIXTURES_ROOT
    _save_json(out_dir / f"{numeric}.json", payload)
    log.info(
        "Recorded %s (returned %s, %s)",
        accession,
        payload.get("chebi_accession", "—"),
        payload.get("name", "—"),
    )
    return numeric


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Record live ChEBI JSON fixtures for chebi_terms tests."
    )
    parser.add_argument(
        "--chebi",
        action="append",
        default=None,
        help=(f"ChEBI ID to record; repeatable. Default: {' '.join(DEFAULT_IDS)}"),
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

    chebi_ids = args.chebi or DEFAULT_IDS
    recorded = 0
    for chebi_id in chebi_ids:
        try:
            if record_fixture_for_id(chebi_id, args.out_dir) is not None:
                recorded += 1
        except requests.exceptions.RequestException as exc:
            log.error("Network error for %s: %s", chebi_id, exc)

    if recorded != len(chebi_ids):
        sys.exit(1)


if __name__ == "__main__":
    main()
