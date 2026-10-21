"""Tests for the guidesig content checks.

Covers the guide_role / guide_target_gene_id consistency check and the
set-wide trailing G / leading C check.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from guidesig import GuideSigError, header_warnings, protospacer_set, signature
from guidesig.cli import main
from tests.guidesig_helpers import FIXTURES, write_delimited

HEADER = [
    "property",
    "guide_id",
    "guide_protospacer",
    "guide_role",
    "guide_target_gene_id",
]
TARGETING = "Targeting a Gene"
GENE = "ENSGFAKE0000001"


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    return write_delimited(path, rows, file_format="tsv")


# --- guide_role / guide_target_gene_id ---------------------------------------


@pytest.mark.parametrize("missing", ["guide_role", "guide_target_gene_id"])
def test_missing_role_or_target_column_raises(tmp_path: Path, missing: str) -> None:
    header = [name for name in HEADER if name != missing]
    row = ["", "g1", "AAAA", "x"]
    path = _write_tsv(tmp_path / "missing.tsv", [header, row])
    with pytest.raises(
        GuideSigError, match=f"required column '{missing}' absent from header"
    ):
        signature(path, file_format="tsv")


def test_control_with_na_target_accepted(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "ok.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", TARGETING, GENE],
            ["", "g2", "CCCA", "control", "NA"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCA"}


def test_control_with_gene_target_raises_with_row_number(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "control_gene.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", TARGETING, GENE],
            ["", "g2", "CCCA", "control", GENE],
        ],
    )
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    message = str(exc.value)
    assert "row 3" in message
    assert f"guide_role is 'control' but guide_target_gene_id is '{GENE}'" in message
    assert "(expected one of 'NA', 'na', 'N/A', 'n/a')" in message


@pytest.mark.parametrize("target", ["NA", "na", "N/A", "n/a", " n/a "])
def test_control_accepts_na_spellings(tmp_path: Path, target: str) -> None:
    path = _write_tsv(
        tmp_path / "control_na_spelling.tsv",
        [HEADER, ["", "g1", "AAAA", "control", target]],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}


@pytest.mark.parametrize("target", ["", "Na", "N/a", "none", "-"])
def test_control_rejects_other_values(tmp_path: Path, target: str) -> None:
    """Only the four NA spellings are accepted; others and empty cells fail."""
    path = _write_tsv(
        tmp_path / "control_not_na.tsv",
        [HEADER, ["", "g1", "AAAA", "control", target]],
    )
    with pytest.raises(GuideSigError, match="row 2"):
        signature(path, file_format="tsv")


@pytest.mark.parametrize("role", [TARGETING, "Targeting a Region"])
def test_non_control_with_na_target_accepted(tmp_path: Path, role: str) -> None:
    """Only control rows are held to NA; other roles may carry NA too."""
    path = _write_tsv(
        tmp_path / "non_control_na.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", "control", "NA"],
            ["", "g2", "CCCA", role, "NA"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCA"}


def test_role_case_and_padding_ignored(tmp_path: Path) -> None:
    """The role matches in any case, and both cells are stripped first."""
    path = _write_tsv(
        tmp_path / "padded_role.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", " Control ", " NA "],
            ["", "g2", "CCCA", "CONTROL", "NA"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCA"}


def test_role_mismatch_in_comment_row_ignored(tmp_path: Path) -> None:
    """Comment rows are skipped before the role check, like every row check."""
    path = _write_tsv(
        tmp_path / "comment_mismatch.tsv",
        [
            HEADER,
            ["#example", "g0", "TTTT", "control", GENE],
            ["", "g1", "AAAA", TARGETING, GENE],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}


def test_role_mismatch_checked_when_protospacer_empty(tmp_path: Path) -> None:
    """The role check runs before rows with no protospacer are skipped."""
    path = _write_tsv(
        tmp_path / "empty_protospacer_mismatch.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", TARGETING, GENE],
            ["", "g2", "", "control", GENE],
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="tsv")


# --- trailing G / leading C ---------------------------------------------------


def _seqs(n: int, *, prefix: str = "", suffix: str = "") -> list[str]:
    """``n`` distinct ACGT protospacers with a fixed start and end."""
    cores = ["".join(bits) for bits in itertools.product("AT", repeat=5)]
    return [f"{prefix}{core}{suffix}" for core in cores[:n]]


def _end_sheet(tmp_path: Path, seqs: list[str], name: str = "ends.tsv") -> Path:
    rows = [["", f"g{i}", seq, TARGETING, GENE] for i, seq in enumerate(seqs)]
    return _write_tsv(tmp_path / name, [HEADER, *rows])


def test_all_trailing_g_raises(tmp_path: Path) -> None:
    path = _end_sheet(tmp_path, _seqs(10, suffix="G"))
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    assert str(exc.value) == (
        f"{path}: Probable trailing G issue, all 10 protospacers end with 'G'"
    )


def test_all_leading_c_raises(tmp_path: Path) -> None:
    path = _end_sheet(tmp_path, _seqs(10, prefix="C"))
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    assert str(exc.value) == (
        f"{path}: Probable leading C issue, all 10 protospacers start with 'C'"
    )


def test_both_conditions_report_trailing_g_first(tmp_path: Path) -> None:
    path = _end_sheet(tmp_path, _seqs(10, prefix="C", suffix="G"))
    with pytest.raises(GuideSigError, match="Probable trailing G issue"):
        signature(path, file_format="tsv")


def test_lowercase_g_and_c_count(tmp_path: Path) -> None:
    """Values are uppercased before the check, so lowercase ends still count."""
    path = _end_sheet(tmp_path, [seq.lower() for seq in _seqs(10, suffix="G")])
    with pytest.raises(GuideSigError, match="Probable trailing G issue"):
        signature(path, file_format="tsv")


def test_count_is_unique_protospacers(tmp_path: Path) -> None:
    """Duplicate rows collapse to the set, so the message counts unique values."""
    seqs = _seqs(10, suffix="G")
    path = _end_sheet(tmp_path, [*seqs, seqs[0], seqs[1]])
    with pytest.raises(GuideSigError, match="all 10 protospacers end with 'G'"):
        signature(path, file_format="tsv")


def test_duplicates_do_not_reach_the_minimum(tmp_path: Path) -> None:
    """Nine unique sequences in twelve rows is below the minimum: not checked."""
    seqs = _seqs(9, suffix="G")
    path = _end_sheet(tmp_path, [*seqs, *seqs[:3]])
    assert protospacer_set(path, file_format="tsv") == set(seqs)


@pytest.mark.parametrize("n", [1, 9])
def test_end_check_skipped_below_minimum(tmp_path: Path, n: int) -> None:
    """Fewer than 10 unique protospacers can all share an end by chance."""
    trailing = _end_sheet(tmp_path, _seqs(n, suffix="G"), "trailing.tsv")
    leading = _end_sheet(tmp_path, _seqs(n, prefix="C"), "leading.tsv")
    assert len(protospacer_set(trailing, file_format="tsv")) == n
    assert len(protospacer_set(leading, file_format="tsv")) == n


def test_mixed_ends_accepted(tmp_path: Path) -> None:
    """Some guides ending in G and some starting with C is not a failure."""
    seqs = [*_seqs(5, suffix="G"), *_seqs(5, prefix="C")]
    path = _end_sheet(tmp_path, seqs)
    assert protospacer_set(path, file_format="tsv") == set(seqs)


def test_cli_reports_end_check_and_prints_no_signature(tmp_path: Path, capsys) -> None:
    path = _end_sheet(tmp_path, _seqs(10, suffix="G"))
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", str(path)])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "error:" in captured.err
    assert "Probable trailing G issue" in captured.err


# --- recommended columns (warning only) ----------------------------------------


def test_header_warnings_name_missing_gene_name_column(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "no_gene_name.tsv", [HEADER, ["", "g1", "AAAA", TARGETING, GENE]]
    )
    assert header_warnings(path, file_format="tsv") == [
        f"{path}: column 'guide_target_gene_name' absent from header"
    ]


def test_header_warnings_empty_when_gene_name_present() -> None:
    assert header_warnings(FIXTURES / "lib_one.tsv", file_format="tsv") == []


def test_header_warnings_ignore_unreadable_file(tmp_path: Path) -> None:
    """The signature read reports unreadable files; the warning pass stays quiet."""
    assert header_warnings(tmp_path / "missing.tsv", file_format="tsv") == []


def test_cli_missing_gene_name_warns_but_still_signs(tmp_path: Path, capsys) -> None:
    path = _write_tsv(
        tmp_path / "no_gene_name.tsv", [HEADER, ["", "g1", "AAAA", TARGETING, GENE]]
    )
    main(["--format", "tsv", str(path)])
    captured = capsys.readouterr()
    assert captured.out.startswith("gsig1:set:n=1:")
    assert (
        f"warning: {path}: column 'guide_target_gene_name' absent from header"
        in captured.err
    )


def test_cli_no_warning_when_gene_name_present(capsys) -> None:
    main(["--format", "tsv", str(FIXTURES / "lib_one.tsv")])
    assert "warning:" not in capsys.readouterr().err


# --- guide_target_gene_name on control rows ------------------------------------

NAME_HEADER = [*HEADER, "guide_target_gene_name"]


def test_control_with_gene_name_raises_with_row_number(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "control_name.tsv",
        [
            NAME_HEADER,
            ["", "g1", "AAAA", TARGETING, GENE, "GENE1"],
            ["", "g2", "CCCA", "control", "NA", "Tfam"],
        ],
    )
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    assert str(exc.value) == (
        f"{path}: row 3: guide_role is 'control' but guide_target_gene_name is "
        "'Tfam' (expected one of 'NA', 'na', 'N/A', 'n/a')"
    )


@pytest.mark.parametrize("name", ["NA", " NA ", "na", "N/A", "n/a"])
def test_control_with_na_gene_name_accepted(tmp_path: Path, name: str) -> None:
    path = _write_tsv(
        tmp_path / "control_name_na.tsv",
        [
            NAME_HEADER,
            ["", "g1", "AAAA", TARGETING, GENE, "GENE1"],
            ["", "g2", "CCCA", "Control", "NA", name],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCA"}


def test_gene_name_rule_skipped_without_column(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "no_name_column.tsv",
        [HEADER, ["", "g1", "AAAA", "control", "NA"]],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}


def test_gene_name_mismatch_in_comment_row_ignored(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "comment_name.tsv",
        [
            NAME_HEADER,
            ["#example", "g0", "TTTT", "control", "NA", "Tfam"],
            ["", "g1", "AAAA", TARGETING, GENE, "GENE1"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}
