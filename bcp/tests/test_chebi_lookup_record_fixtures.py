"""Tests for chebi_lookup.record_fixtures (mocked HTTP)."""

from __future__ import annotations

import logging

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from chebi_lookup.record_fixtures import record_fixtures_for_cas


@patch("chebi_lookup.record_fixtures.get_with_retry")
def test_recorder_refuses_a_cell_that_is_not_a_registry_number(
    mock_get: MagicMock, tmp_path: Path
) -> None:
    """The same cells cas_to_cid_status refuses must not produce a fixture."""
    assert record_fixtures_for_cas("not a cas", tmp_path) is None
    mock_get.assert_not_called()


@patch("chebi_lookup.record_fixtures.time.sleep")
@patch("chebi_lookup.record_fixtures.get_with_retry")
def test_recorder_queries_and_names_the_directory_after_the_normalised_value(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A repaired CAS must not write a fixture for a URL production never issues.

    `0362-07-02` is 2-methoxyestradiol after a padded check digit and a leading
    zero. Production looks up `362-07-2`; a fixture keyed on the raw cell would
    make the live test pass against a request path the shipped code cannot reach.
    """

    def _json(payload: dict) -> MagicMock:
        resp = MagicMock()
        resp.json.return_value = payload
        return resp

    # The xrefs slot carries a ChEBI ID and the synonyms slot carries none, so the
    # two payloads are distinguishable. With identical empty shapes, reading
    # RegistryID back out of the *synonyms* response -- the bug 24c3399 fixed, which
    # reported "ChEBI --" for every compound -- produced the same result and this
    # test passed either way.
    mock_get.side_effect = [
        _json({"IdentifierList": {"CID": [123]}}),
        _json({"PropertyTable": {"Properties": [{}]}}),
        _json({"InformationList": {"Information": [{"RegistryID": ["CHEBI:16236"]}]}}),
        _json(
            {"InformationList": {"Information": [{"Synonym": ["2-methoxyestradiol"]}]}}
        ),
    ]

    with caplog.at_level(logging.INFO, logger="chebi_lookup.record_fixtures"):
        assert record_fixtures_for_cas("0362-07-02", tmp_path) == 123
    url = mock_get.call_args_list[0][0][0]
    assert "/name/362-07-2/" in url
    assert (tmp_path / "362-07-2" / "cids.json").exists()
    assert not (tmp_path / "0362-07-02").exists()
    # The xref is read from the xrefs response, not from whatever `resp` last held.
    assert "ChEBI CHEBI:16236" in caplog.text
