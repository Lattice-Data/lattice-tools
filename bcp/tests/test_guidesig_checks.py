"""Tests for the guidesig content checks.

Covers the guide_role / guide_target_gene_id consistency check and the
set-wide trailing G / leading C check.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from guidesig import GuideSigError, protospacer_set, signature
from guidesig.cli import main
from tests.guidesig_helpers import write_delimited

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
    assert "(expected 'NA')" in message


@pytest.mark.parametrize("target", ["", "na", "N/A"])
def test_control_requires_exact_na(tmp_path: Path, target: str) -> None:
    """Only the literal NA is accepted; look-alikes and empty cells fail."""
    path = _write_tsv(
        tmp_path / "control_not_na.tsv",
        [HEADER, ["", "g1", "AAAA", "control", target]],
    )
    with pytest.raises(GuideSigError, match="row 2"):
        signature(path, file_format="tsv")


def test_non_control_with_na_target_raises_with_row_number(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "targeting_na.tsv",
        [
            HEADER,
            ["", "g1", "AAAA", "control", "NA"],
            ["", "g2", "CCCA", TARGETING, "NA"],
        ],
    )
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    message = str(exc.value)
    assert "row 3" in message
    assert f"guide_target_gene_id is 'NA' but guide_role is '{TARGETING}'" in message
    assert "(expected 'control')" in message


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


def test_all_trailing_g_raises(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "trailing_g.tsv",
        [
            HEADER,
            ["", "g1", "AAAG", TARGETING, GENE],
            ["", "g2", "TTTG", TARGETING, GENE],
        ],
    )
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    assert str(exc.value) == (
        f"{path}: Probable trailing G issue, all 2 protospacers end with 'G'"
    )


def test_all_leading_c_raises(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "leading_c.tsv",
        [
            HEADER,
            ["", "g1", "CAAA", TARGETING, GENE],
            ["", "g2", "CTTT", TARGETING, GENE],
        ],
    )
    with pytest.raises(GuideSigError) as exc:
        signature(path, file_format="tsv")
    assert str(exc.value) == (
        f"{path}: Probable leading C issue, all 2 protospacers start with 'C'"
    )


def test_both_conditions_report_trailing_g_first(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "both.tsv",
        [
            HEADER,
            ["", "g1", "CAAG", TARGETING, GENE],
            ["", "g2", "CTTG", TARGETING, GENE],
        ],
    )
    with pytest.raises(GuideSigError, match="Probable trailing G issue"):
        signature(path, file_format="tsv")


def test_lowercase_g_and_c_count(tmp_path: Path) -> None:
    """Values are uppercased before the check, so lowercase ends still count."""
    path = _write_tsv(
        tmp_path / "lowercase.tsv",
        [
            HEADER,
            ["", "g1", "aaag", TARGETING, GENE],
            ["", "g2", "tttg", TARGETING, GENE],
        ],
    )
    with pytest.raises(GuideSigError, match="Probable trailing G issue"):
        signature(path, file_format="tsv")


def test_count_is_unique_protospacers(tmp_path: Path) -> None:
    """Duplicate rows collapse to the set, so the message counts unique values."""
    path = _write_tsv(
        tmp_path / "duplicates.tsv",
        [
            HEADER,
            ["", "g1", "AAAG", TARGETING, GENE],
            ["", "g1b", "AAAG", TARGETING, GENE],
            ["", "g2", "TTTG", TARGETING, GENE],
        ],
    )
    with pytest.raises(GuideSigError, match="all 2 protospacers end with 'G'"):
        signature(path, file_format="tsv")


def test_single_guide_ending_in_g_raises(tmp_path: Path) -> None:
    """There is no minimum set size: one guide ending in G is all of them."""
    path = _write_tsv(
        tmp_path / "single.tsv",
        [HEADER, ["", "g1", "AAAG", TARGETING, GENE]],
    )
    with pytest.raises(GuideSigError, match="all 1 protospacers end with 'G'"):
        signature(path, file_format="tsv")


def test_mixed_ends_accepted(tmp_path: Path) -> None:
    """Some guides ending in G and some starting with C is not a failure."""
    path = _write_tsv(
        tmp_path / "mixed.tsv",
        [
            HEADER,
            ["", "g1", "AAAG", TARGETING, GENE],
            ["", "g2", "CTTT", TARGETING, GENE],
            ["", "g3", "CAAG", TARGETING, GENE],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAG", "CTTT", "CAAG"}


def test_cli_reports_end_check_and_prints_no_signature(tmp_path: Path, capsys) -> None:
    path = _write_tsv(
        tmp_path / "trailing_g.tsv",
        [
            HEADER,
            ["", "g1", "AAAG", TARGETING, GENE],
            ["", "g2", "TTTG", TARGETING, GENE],
        ],
    )
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", str(path)])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "error:" in captured.err
    assert "Probable trailing G issue" in captured.err
