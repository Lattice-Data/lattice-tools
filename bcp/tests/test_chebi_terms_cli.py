"""CLI tests for chebi_terms.main."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

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


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_empty_chebi_arg_reports_missing_not_crash(
    mock_get: MagicMock,
    _sleep: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """`--chebi ""` is falsy but must not fall through to batch mode."""
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    sys.argv = ["chebi_terms", "--chebi", ""]
    chebi_terms.cli.main()

    data = json.loads(capsys.readouterr().out)
    assert data["id_status"] == "missing"
    mock_get.assert_not_called()


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_exits_1_when_chebi_unreachable(
    mock_get: MagicMock,
    _sleep: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = requests.exceptions.ConnectionError("down")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--chebi", "CHEBI:16236"]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
    # The JSON still lands on stdout, saying lookup_failed rather than not_found.
    data = json.loads(capsys.readouterr().out)
    assert data["id_status"] == "lookup_failed"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_exits_0_for_an_unwelcome_verdict(
    mock_get: MagicMock,
    _sleep: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A mismatch is a successful run: the tool did its job and said no."""
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    sys.argv = ["chebi_terms", "--chebi", "CHEBI:16236", "--expect-name", "caffeine"]
    chebi_terms.cli.main()

    assert json.loads(capsys.readouterr().out)["name_verdict"] == "mismatch"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_exits_1_when_any_row_lookup_failed(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_terms.cli
    from chebi_terms.client import MAX_RETRIES

    calls = {"n": 0}

    def flaky(url: str, *args, **kwargs):
        calls["n"] += 1
        if calls["n"] <= MAX_RETRIES:
            raise requests.exceptions.ConnectionError("blip")
        return route_chebi_get(url)

    mock_get.side_effect = flaky

    src = tmp_path / "in.csv"
    src.write_text("chebi_id\nCHEBI:16236\nCHEBI:16236\n", encoding="utf-8")
    out = tmp_path / "out.csv"

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--input", str(src), "--output", str(out)]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1

    # The degraded run still leaves a usable CSV behind.
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert [r["id_status"] for r in rows] == ["lookup_failed", "ok"]


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_exits_1_when_every_lookup_missed(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    """The other half of the exit-1 contract: a wholesale 404."""
    import chebi_terms.cli
    from chebi_terms.io import SUSPICIOUS_TOTAL_MISS

    mock_get.return_value = MagicMock(status_code=404)

    src = tmp_path / "in.csv"
    ids = "\n".join(f"CHEBI:{16236 + i}" for i in range(SUSPICIOUS_TOTAL_MISS))
    src.write_text(f"chebi_id\n{ids}\n", encoding="utf-8")
    out = tmp_path / "out.csv"

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--input", str(src), "--output", str(out)]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1

    # Per-row verdicts stay honest; the run-level verdict is what changed.
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert all(row["id_status"] == "not_found" for row in rows)


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_exits_0_when_some_ids_resolve(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    """A genuine not_found alongside a hit is data, not an outage."""
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    src = tmp_path / "in.csv"
    src.write_text("chebi_id\nCHEBI:16236\nCHEBI:99999999\n", encoding="utf-8")
    out = tmp_path / "out.csv"

    sys.argv = ["chebi_terms", "--input", str(src), "--output", str(out)]
    chebi_terms.cli.main()

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert [r["id_status"] for r in rows] == ["ok", "not_found"]


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_exits_1_when_no_row_held_a_usable_id(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    """--chebi-column pointed at the wrong column must not look like success."""
    import chebi_terms.cli
    from chebi_terms.io import SUSPICIOUS_TOTAL_MISS

    mock_get.side_effect = route_chebi_get

    src = tmp_path / "in.csv"
    names = "\n".join(f"compound-{i}" for i in range(SUSPICIOUS_TOTAL_MISS))
    src.write_text(f"compound_name\n{names}\n", encoding="utf-8")
    out = tmp_path / "out.csv"

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = [
            "chebi_terms",
            "--input",
            str(src),
            "--chebi-column",
            "compound_name",
            "--output",
            str(out),
        ]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
    mock_get.assert_not_called()


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_single_exits_1_on_unwritable_output(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = [
            "chebi_terms",
            "--chebi",
            "CHEBI:16236",
            "--output",
            str(tmp_path / "nope" / "out.json"),
        ]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
    mock_get.assert_not_called()


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_cli_batch_exits_1_on_unwritable_output(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    """The batch counterpart to test_cli_single_exits_1_on_unwritable_output."""
    import chebi_terms.cli

    mock_get.side_effect = route_chebi_get

    src = tmp_path / "in.csv"
    src.write_text("chebi_id\nCHEBI:16236\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = [
            "chebi_terms",
            "--input",
            str(src),
            "--output",
            str(tmp_path / "nope" / "out.csv"),
        ]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
    mock_get.assert_not_called()


def test_cli_batch_missing_column_exits_1(tmp_path: Path) -> None:
    import chebi_terms.cli

    src = tmp_path / "in.csv"
    src.write_text("other\nx\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_terms", "--input", str(src)]
        chebi_terms.cli.main()
    assert exc_info.value.code == 1
