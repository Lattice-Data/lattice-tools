"""CLI tests for chebi_lookup.main."""

from __future__ import annotations

import sys

import pytest


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    import chebi_lookup.cli

    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["chebi_lookup", "--help"]
        chebi_lookup.cli.main()
    assert exc_info.value.code == 0
    assert "ChEBI" in capsys.readouterr().out
