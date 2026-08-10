"""CLI tests for chebi_terms.main."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tests.chebi_terms_helpers import route_chebi_get


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    import chebi_terms.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--help"]
        chebi_terms.cli.main()
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "ChEBI" in out
    assert "--expect-name" in out
    assert "--cas-column" in out


def test_cli_requires_input_or_chebi() -> None:
    import chebi_terms.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms"]
        chebi_terms.cli.main()
    assert exc_info.value.code == 2


def test_cli_rejects_both_modes(tmp_path: Path) -> None:
    import chebi_terms.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = [
            "chebi_terms",
            "--input",
            str(tmp_path / "in.csv"),
            "--chebi",
            "CHEBI:16236",
        ]
        chebi_terms.cli.main()
    assert exc_info.value.code == 2


def test_cli_rejects_zero_max_synonyms() -> None:
    import chebi_terms.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--chebi", "CHEBI:16236", "--max-synonyms", "0"]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_id_json(
    mock_get: MagicMock,
    _sleep: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    sys.argv = ["chebi_terms", "--chebi", "CHEBI:16236"]
    chebi_terms.cli.main()

    data = json.loads(capsys.readouterr().out)
    assert data["chebi_id"] == "CHEBI:16236"
    assert data["chebi_accession"] == "CHEBI:16236"
    assert data["chebi_name"] == "ethanol"
    assert data["id_status"] == "ok"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_id_with_expectations(
    mock_get: MagicMock,
    _sleep: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    sys.argv = [
        "chebi_terms",
        "--chebi",
        "18484",
        "--expect-name",
        "(-)-epicatechin",
        "--expect-cas",
        "490-46-0",
    ]
    chebi_terms.cli.main()

    data = json.loads(capsys.readouterr().out)
    assert data["id_status"] == "secondary"
    assert data["chebi_accession"] == "CHEBI:90"
    assert data["name_verdict"] == "match"
    assert data["cas_verdict"] == "confirmed"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_id_csv_format(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get
    out = tmp_path / "one.csv"

    sys.argv = [
        "chebi_terms",
        "--chebi",
        "CHEBI:16236",
        "--format",
        "csv",
        "--output",
        str(out),
    ]
    chebi_terms.cli.main()

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["chebi_name"] == "ethanol"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_with_checks(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    src = tmp_path / "in.csv"
    src.write_text(
        "chebi_id,name,CAS\nCHEBI:16236,ethanol,64-17-5\nCHEBI:16236,caffeine,58-08-2\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.csv"

    sys.argv = [
        "chebi_terms",
        "--input",
        str(src),
        "--chebi-column",
        "chebi_id",
        "--name-column",
        "name",
        "--cas-column",
        "CAS",
        "--output",
        str(out),
    ]
    chebi_terms.cli.main()

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert [r["name_verdict"] for r in rows] == ["match", "mismatch"]
    assert [r["cas_verdict"] for r in rows] == ["confirmed", "not_in_chebi"]


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_default_output_path(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    src = tmp_path / "ids.csv"
    src.write_text("chebi_id\nCHEBI:16236\n", encoding="utf-8")

    sys.argv = ["chebi_terms", "--input", str(src)]
    chebi_terms.cli.main()

    expected = tmp_path / "ids_chebi_checked.csv"
    assert expected.exists()
    assert "ethanol" in expected.read_text(encoding="utf-8")


def test_cli_batch_missing_column_exits_1(tmp_path: Path) -> None:
    import chebi_terms.cli

    src = tmp_path / "in.csv"
    src.write_text("other\nx\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--input", str(src)]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
