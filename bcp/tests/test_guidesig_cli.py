"""CLI tests for guidesig."""

from __future__ import annotations

from pathlib import Path

import pytest

from guidesig.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


def test_cli_multi_file_prints_signature_and_path(capsys) -> None:
    paths = [
        FIXTURES / "lib_two_layout_a.tsv",
        FIXTURES / "lib_two_layout_b.tsv",
    ]
    main([str(p) for p in paths])
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) == 2
    for line, path in zip(out, paths, strict=True):
        sig, reported_path = line.split("\t", 1)
        assert sig.startswith("gsig1:set:n=60:")
        assert reported_path == str(path)


def test_cli_compare_lib_two_match(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--compare",
                str(FIXTURES / "lib_two_layout_a.tsv"),
                str(FIXTURES / "lib_two_layout_b.tsv"),
            ]
        )
    assert exc.value.code == 0
    assert "match" in capsys.readouterr().out


def test_cli_compare_lib_one_mismatch(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--compare",
                str(FIXTURES / "lib_one.tsv"),
                str(FIXTURES / "lib_one_prepended_g.tsv"),
            ]
        )
    assert exc.value.code != 0
    out = capsys.readouterr().out
    assert "mismatch" in out
    assert "|A|=" in out
    assert "|A - B|=" in out
    assert "A - B examples:" in out
    assert "B - A examples:" in out


def test_cli_error_exits_nonzero(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURES / "lib_one.tsv"), "--column", "missing_column"])
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert "error:" in err
    assert "absent from header" in err
