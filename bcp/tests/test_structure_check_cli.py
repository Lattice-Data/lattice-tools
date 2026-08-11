"""CLI tests for structure_check.main."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from structure_check.client import REVIEW_INVESTIGATE, REVIEW_OK, empty_result


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    import structure_check.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["structure_check", "--help"]
        structure_check.cli.main()
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "InChIKey" in out
    assert "--chebi-column" in out


def test_cli_requires_input_or_cas() -> None:
    import structure_check.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["structure_check"]
        structure_check.cli.main()
    assert exc_info.value.code == 2


def test_cli_single_needs_something_to_compare_against() -> None:
    import structure_check.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["structure_check", "--cas", "64-17-5"]
        structure_check.cli.main()
    assert exc_info.value.code == 1


@patch("structure_check.io.check_row")
def test_cli_single_json(
    mock_check: MagicMock, capsys: pytest.CaptureFixture[str]
) -> None:
    import structure_check.cli

    mock_check.return_value = {**empty_result(), "review": REVIEW_INVESTIGATE}
    sys.argv = [
        "structure_check",
        "--cas",
        "22573-88-2",
        "--name",
        "Alexidine_dihydrochloride",
        "--chebi",
        "CHEBI:27391",
    ]
    structure_check.cli.main()

    payload = json.loads(capsys.readouterr().out)
    assert payload["review"] == REVIEW_INVESTIGATE
    assert payload["chebi_id"] == "CHEBI:27391"


@patch("structure_check.io.check_row")
def test_cli_batch_writes_review_column(mock_check: MagicMock, tmp_path: Path) -> None:
    import structure_check.cli

    mock_check.return_value = {**empty_result(), "review": REVIEW_OK}
    src = tmp_path / "in.csv"
    src.write_text("Name,CAS\nEthanol,64-17-5\n", encoding="utf-8")
    out = tmp_path / "out.csv"

    sys.argv = [
        "structure_check",
        "--input",
        str(src),
        "--cas-column",
        "CAS",
        "--name-column",
        "Name",
        "--output",
        str(out),
    ]
    structure_check.cli.main()

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["review"] == REVIEW_OK


@patch("structure_check.io.check_row")
def test_cli_batch_default_output_path(mock_check: MagicMock, tmp_path: Path) -> None:
    import structure_check.cli

    mock_check.return_value = empty_result()
    src = tmp_path / "sheet.csv"
    src.write_text("Name,CAS\nEthanol,64-17-5\n", encoding="utf-8")

    sys.argv = ["structure_check", "--input", str(src), "--name-column", "Name"]
    structure_check.cli.main()

    assert (tmp_path / "sheet_structure_checked.csv").exists()


@patch("structure_check.io.check_row")
def test_cli_batch_findings_still_exit_0(mock_check: MagicMock, tmp_path: Path) -> None:
    """Findings are the product; only a broken run is a failure."""
    import structure_check.cli

    mock_check.return_value = {**empty_result(), "review": REVIEW_INVESTIGATE}
    src = tmp_path / "in.csv"
    src.write_text("Name,CAS\nAlexidine,22573-88-2\n", encoding="utf-8")

    sys.argv = [
        "structure_check",
        "--input",
        str(src),
        "--name-column",
        "Name",
        "--output",
        str(tmp_path / "out.csv"),
    ]
    structure_check.cli.main()  # must not raise SystemExit


def test_cli_batch_missing_column_exits_1(tmp_path: Path) -> None:
    import structure_check.cli

    src = tmp_path / "in.csv"
    src.write_text("other\nx\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["structure_check", "--input", str(src), "--name-column", "other"]
        structure_check.cli.main()
    assert exc_info.value.code == 1


def test_cli_batch_without_anything_to_check_exits_1(tmp_path: Path) -> None:
    import structure_check.cli

    src = tmp_path / "in.csv"
    src.write_text("CAS\n64-17-5\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["structure_check", "--input", str(src), "--cas-column", "CAS"]
        structure_check.cli.main()
    assert exc_info.value.code == 1
