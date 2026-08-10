"""Tests for chebi_terms.record_fixtures (mocked HTTP)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from chebi_terms.client import ChebiUnavailableError
from chebi_terms.record_fixtures import (
    DEFAULT_IDS,
    record_fixture_for_id,
)
from tests.chebi_terms_helpers import ETHANOL, load_payload


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_record_fixture_writes_payload_keyed_by_requested_id(
    mock_fetch: MagicMock, tmp_path: Path
) -> None:
    mock_fetch.return_value = load_payload(ETHANOL)

    assert record_fixture_for_id("CHEBI:16236", tmp_path) == ETHANOL

    written = json.loads((tmp_path / f"{ETHANOL}.json").read_text(encoding="utf-8"))
    assert written["name"] == "ethanol"


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_record_fixture_keys_secondary_id_by_what_was_asked_for(
    mock_fetch: MagicMock, tmp_path: Path
) -> None:
    """Keyed by the requested ID so secondary-ID behaviour stays reproducible."""
    mock_fetch.return_value = {"chebi_accession": "CHEBI:90", "name": "x"}

    assert record_fixture_for_id("CHEBI:18484", tmp_path) == "18484"

    assert (tmp_path / "18484.json").exists()
    assert not (tmp_path / "90.json").exists()


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_record_fixture_rejects_junk_without_fetching(
    mock_fetch: MagicMock, tmp_path: Path
) -> None:
    assert record_fixture_for_id("not-an-id", tmp_path) is None
    mock_fetch.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_record_fixture_returns_none_when_no_such_compound(
    mock_fetch: MagicMock, tmp_path: Path
) -> None:
    mock_fetch.return_value = None
    assert record_fixture_for_id("CHEBI:99999999", tmp_path) is None
    assert list(tmp_path.iterdir()) == []


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_record_fixture_propagates_unavailable(
    mock_fetch: MagicMock, tmp_path: Path
) -> None:
    """main() catches this to log a clean error; it must not be swallowed here."""
    mock_fetch.side_effect = ChebiUnavailableError("down")
    with pytest.raises(ChebiUnavailableError):
        record_fixture_for_id("CHEBI:16236", tmp_path)


@patch("chebi_terms.record_fixtures.fetch_compound")
def test_main_reports_unreachable_chebi_without_traceback(
    mock_fetch: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    import sys

    import chebi_terms.record_fixtures as rf

    mock_fetch.side_effect = ChebiUnavailableError("down")
    sys.argv = ["record_fixtures", "--chebi", "CHEBI:16236", "--out-dir", str(tmp_path)]

    with pytest.raises(SystemExit) as exc_info:
        rf.main()

    assert exc_info.value.code == 1
    assert "Could not reach ChEBI" in caplog.text


def test_default_ids_are_the_recorded_fixtures() -> None:
    # Guards against the defaults drifting away from what tests load.
    assert "CHEBI:16236" in DEFAULT_IDS
    assert len(DEFAULT_IDS) == 3
