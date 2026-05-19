"""Unit tests for chebi_lookup client and io (mocked HTTP)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from chebi_lookup.client import (
    OUTPUT_FIELDS_APPENDED,
    cas_to_cid,
    lookup_cas,
)
from chebi_lookup.io import (
    CasMappingError,
    build_single_cas_row,
    emit_single_cas,
    map_cas_file,
)
from tests.chebi_lookup_helpers import (
    FIXTURES,
    load_json_name,
    mock_response,
    route_pubchem_get,
)


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_to_cid_success(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.side_effect = route_pubchem_get
    assert cas_to_cid("64-17-5") == 702


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_to_cid_not_found(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.return_value = mock_response(404)
    assert cas_to_cid("00-00-0") is None


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_lookup_cas_fills_chebi(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.side_effect = route_pubchem_get
    result = lookup_cas("64-17-5")
    assert result["pubchem_cid"] == 702
    assert result["chebi_id"] == "CHEBI:16236"
    assert result["preferred_name"] == "Ethanol"
    assert "ethanol" in result["synonyms"].lower()


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_lookup_cas_no_chebi_xref(mock_get: MagicMock, _sleep: MagicMock) -> None:
    def route_no_chebi(url: str, *args, **kwargs):
        if "/xrefs/RegistryID/JSON" in url:
            return mock_response(200, load_json_name("registry_ids_no_chebi.json"))
        return route_pubchem_get(url, *args, **kwargs)

    mock_get.side_effect = route_no_chebi
    result = lookup_cas("64-17-5")
    assert result["pubchem_cid"] == 702
    assert result["chebi_id"] == ""


def test_lookup_cas_empty() -> None:
    result = lookup_cas("  ")
    assert result == {field: "" for field in OUTPUT_FIELDS_APPENDED}


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_build_single_cas_row(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.side_effect = route_pubchem_get
    row = build_single_cas_row("64-17-5", lookup_cas("64-17-5"))
    assert row["CAS"] == "64-17-5"
    assert row["pubchem_cid"] == 702
    assert row["chebi_id"] == "CHEBI:16236"


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_emit_single_cas_json_file(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    import json

    mock_get.side_effect = route_pubchem_get
    out = tmp_path / "result.json"
    emit_single_cas("64-17-5", out, fmt="json")
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["chebi_id"] == "CHEBI:16236"


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_map_cas_file(mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path) -> None:
    mock_get.side_effect = route_pubchem_get
    input_path = FIXTURES / "sample_input.csv"
    output_path = tmp_path / "out.csv"
    map_cas_file(input_path, "CAS", output_path)

    text = output_path.read_text(encoding="utf-8")
    assert "chebi_id" in text
    assert "CHEBI:16236" in text
    assert "cmp_003" in text


def test_map_cas_file_missing_column(tmp_path: Path) -> None:
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("id,notes\n1,foo\n", encoding="utf-8")
    with pytest.raises(CasMappingError, match="Column 'CAS' not found"):
        map_cas_file(bad_csv, "CAS", tmp_path / "out.csv")


def test_map_cas_file_missing_input(tmp_path: Path) -> None:
    with pytest.raises(CasMappingError, match="Input file not found"):
        map_cas_file(tmp_path / "missing.csv", "CAS", tmp_path / "out.csv")
