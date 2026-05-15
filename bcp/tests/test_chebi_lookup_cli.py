"""CLI tests for chebi_lookup.main."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tests.chebi_lookup_helpers import FIXTURES, route_pubchem_get


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    import chebi_lookup.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_lookup", "--help"]
        chebi_lookup.cli.main()
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "CAS" in out
    assert "PubChem" in out
    assert "ChEBI" in out


def test_cli_missing_input() -> None:
    import chebi_lookup.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_lookup"]
        chebi_lookup.cli.main()
    assert exc_info.value.code == 2


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cli_runs_with_mocked_http(
    mock_get: MagicMock,
    _sleep: MagicMock,
    tmp_path: Path,
) -> None:
    import chebi_lookup.cli

    mock_get.side_effect = route_pubchem_get

    input_copy = tmp_path / "in.csv"
    input_copy.write_text(
        (FIXTURES / "sample_input.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_path = tmp_path / "out.csv"

    sys.argv = [
        "chebi_lookup",
        "--input",
        str(input_copy),
        "--cas-column",
        "CAS",
        "--output",
        str(output_path),
    ]
    chebi_lookup.cli.main()

    content = output_path.read_text(encoding="utf-8")
    header = content.splitlines()[0]
    assert "pubchem_cid" in header
    assert "chebi_id" in header
    assert "CHEBI:16236" in content
