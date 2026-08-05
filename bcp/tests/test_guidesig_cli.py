"""CLI tests for guidesig."""

from __future__ import annotations

import pytest

from guidesig.cli import main
from tests.guidesig_helpers import FIXTURES


def test_cli_multi_file_prints_signature_and_path(capsys) -> None:
    paths = [
        FIXTURES / "lib_two_layout_a.tsv",
        FIXTURES / "lib_two_layout_b.tsv",
    ]
    main(["--format", "tsv", *[str(p) for p in paths]])
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) == 2
    for line, path in zip(out, paths, strict=True):
        sig, reported_path = line.split("\t", 1)
        assert sig.startswith("gsig1:set:n=60:")
        assert reported_path == str(path)


def test_cli_csv_matches_tsv_signature(capsys) -> None:
    """The same library as CSV and as TSV must report the same signature."""
    main(["--format", "csv", str(FIXTURES / "lib_two_layout_a.csv")])
    csv_sig = capsys.readouterr().out.split("\t", 1)[0]
    main(["--format", "tsv", str(FIXTURES / "lib_two_layout_a.tsv")])
    tsv_sig = capsys.readouterr().out.split("\t", 1)[0]
    assert csv_sig == tsv_sig == "gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78"


def test_cli_compare_lib_two_match(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--format",
                "tsv",
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
                "--format",
                "tsv",
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
        main(
            [
                "--format",
                "tsv",
                str(FIXTURES / "lib_one.tsv"),
                "--column",
                "missing_column",
            ]
        )
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert "error:" in err
    assert "absent from header" in err


def test_cli_format_is_required(capsys) -> None:
    """Omitting --format must fail rather than assuming a format."""
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURES / "lib_one.tsv")])
    assert exc.value.code == 2
    assert "--format" in capsys.readouterr().err


def test_cli_rejects_unknown_format(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "xlsx", str(FIXTURES / "lib_one.tsv")])
    assert exc.value.code == 2
    assert "--format" in capsys.readouterr().err


def test_cli_wrong_format_exits_nonzero(capsys) -> None:
    """Reading the CSV fixture as TSV fails loudly instead of misparsing."""
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", str(FIXTURES / "lib_two_layout_a.csv")])
    assert exc.value.code != 0
    assert "absent from header" in capsys.readouterr().err


def test_cli_compare_across_formats_names_the_format(capsys) -> None:
    """--compare reads both files under one --format, so a mixed pair must say so."""
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--format",
                "csv",
                "--compare",
                str(FIXTURES / "lib_two_layout_a.csv"),
                str(FIXTURES / "lib_two_layout_b.tsv"),
            ]
        )
    assert exc.value.code != 0
    captured = capsys.readouterr()
    assert captured.out == ""
    err = captured.err
    assert "lib_two_layout_b.tsv" in err
    assert "looks like tsv rather than csv" in err


def test_cli_multi_file_failure_emits_no_partial_stdout(tmp_path, capsys) -> None:
    """A later failing file must not leave earlier signatures on stdout."""
    bad = tmp_path / "bad.tsv"
    bad.write_bytes(b"property\tguide_id\r\n\tg1\r\n")  # no protospacer column
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", str(FIXTURES / "lib_one.tsv"), str(bad)])
    assert exc.value.code != 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "error:" in captured.err
