"""CSV-input tests for guidesig.

The parsing rules are shared with TSV, but the reported row numbers and the
quoting behaviour are what users act on, so they are re-asserted here rather
than assumed from the TSV suite.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from guidesig import GuideSigError, protospacer_set, signature
from tests.guidesig_helpers import FIXTURES, write_delimited

LAYOUT_A_DIGEST = "gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78"


def _write_csv(path: Path, rows: list[list[str]], *, bom: bool = True) -> Path:
    return write_delimited(path, rows, file_format="csv", bom=bom)


def test_csv_fixture_matches_tsv_golden_digest() -> None:
    """The CSV conversion of lib_two_layout_a must hit the pinned TSV digest."""
    assert signature(FIXTURES / "lib_two_layout_a.csv", file_format="csv") == (
        LAYOUT_A_DIGEST
    )


def test_same_sequences_in_csv_and_tsv_match(tmp_path: Path) -> None:
    """The signature depends on the protospacer set, not on the file format."""
    rows = [
        ["property", "guide_id", "guide_protospacer"],
        ["", "g1", "AAAA"],
        ["", "g2", "CCCC"],
    ]
    as_csv = _write_csv(tmp_path / "lib.csv", rows)
    as_tsv = write_delimited(tmp_path / "lib.tsv", rows, file_format="tsv")
    assert signature(as_csv, file_format="csv") == signature(as_tsv, file_format="tsv")


def test_quoted_comma_does_not_shift_protospacer_column(tmp_path: Path) -> None:
    """Commas inside quoted cells on both sides of the target column."""
    path = _write_csv(
        tmp_path / "quoted.csv",
        [
            ["property", "guide_id", "guide_protospacer", "gene"],
            ["", "g1, replicate A", "AAAA", "GENE1"],
            ["", "g2", "CCCC", "GENE2, alias"],
        ],
    )
    assert protospacer_set(path, file_format="csv") == {"AAAA", "CCCC"}


def test_quoted_protospacer_padding_and_cr_stripped(tmp_path: Path) -> None:
    clean = _write_csv(
        tmp_path / "clean.csv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    padded = tmp_path / "padded.csv"
    # Quote the field so the embedded CR is data, not a line terminator.
    padded.write_bytes(
        b'\xef\xbb\xbfproperty,guide_protospacer\r\n,  AAAA  \r\n,"CCCC\r"\r\n'
    )
    assert signature(clean, file_format="csv") == signature(padded, file_format="csv")


def test_csv_without_bom_parses(tmp_path: Path) -> None:
    path = _write_csv(
        tmp_path / "nobom.csv",
        [["property", "guide_protospacer"], ["", "AAAA"]],
        bom=False,
    )
    assert not path.read_bytes().startswith(b"\xef\xbb\xbf")
    assert signature(path, file_format="csv").startswith("gsig1:set:n=1:")


def test_csv_comment_rows_filtered(tmp_path: Path) -> None:
    path = _write_csv(
        tmp_path / "comments.csv",
        [
            ["property", "guide_protospacer"],
            ["#comment", "TTTTTTTTTTTTTTTTTTT"],
            ["# example", "GGGGGGGGGGGGGGGGGGG"],
            ["", "AAAA"],
        ],
    )
    assert protospacer_set(path, file_format="csv") == {"AAAA"}


def test_csv_blank_lines_skipped(tmp_path: Path) -> None:
    path = tmp_path / "blank_lines.csv"
    path.write_bytes(
        b"\xef\xbb\xbfproperty,guide_protospacer\r\n,AAAA\r\n\r\n,CCCC\r\n\r\n"
    )
    assert protospacer_set(path, file_format="csv") == {"AAAA", "CCCC"}


def test_csv_ragged_row_raises_with_row_number(tmp_path: Path) -> None:
    path = _write_csv(
        tmp_path / "ragged.csv",
        [
            ["property", "guide_id", "guide_protospacer"],
            ["", "g1", "AAAA"],
            ["", "g2"],
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="csv")


@pytest.mark.parametrize("bad", ["N", "U", "0", "-"])
def test_csv_non_acgt_raises_with_row_number(tmp_path: Path, bad: str) -> None:
    path = _write_csv(
        tmp_path / "bad.csv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", f"AAA{bad}"],
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="csv")


def test_csv_multiline_quoted_cell_keeps_row_numbers_physical(tmp_path: Path) -> None:
    """A spreadsheet cell with an embedded newline must not shift row numbers."""
    path = _write_csv(
        tmp_path / "multiline.csv",
        [
            ["property", "guide_target_region", "guide_protospacer"],
            ["", "chrFAKE:1000-1100\nsecond locus", "AAAA"],
            ["", "chrFAKE:1200-1300", "AAAN"],  # non-ACGT, on physical line 4
        ],
    )
    assert path.read_bytes().count(b"\n") == 4
    with pytest.raises(GuideSigError, match="row 4"):
        signature(path, file_format="csv")


def test_csv_missing_column_raises(tmp_path: Path) -> None:
    path = _write_csv(
        tmp_path / "missing.csv",
        [["property", "guide_id"], ["", "g1"]],
    )
    with pytest.raises(GuideSigError, match="absent from header"):
        signature(path, file_format="csv")


def test_csv_read_as_tsv_reports_missing_column() -> None:
    """Naming the wrong format fails loudly instead of misparsing."""
    with pytest.raises(GuideSigError, match="absent from header"):
        signature(FIXTURES / "lib_two_layout_a.csv", file_format="tsv")


def test_tsv_read_as_csv_reports_missing_column() -> None:
    with pytest.raises(GuideSigError, match="absent from header"):
        signature(FIXTURES / "lib_two_layout_a.tsv", file_format="csv")


def test_csv_error_message_names_csv(tmp_path: Path) -> None:
    """An empty file must not be reported as a TSV problem."""
    path = tmp_path / "empty.csv"
    path.write_bytes(b"")
    with pytest.raises(GuideSigError, match="not CSV-parseable"):
        signature(path, file_format="csv")


def test_csv_alternate_column_honored(tmp_path: Path) -> None:
    """--column selects a different sequence column in CSV as it does in TSV."""
    path = _write_csv(
        tmp_path / "two_columns.csv",
        [
            ["property", "guide_rc_sequence", "guide_protospacer"],
            ["", "TTTT", "AAAA"],
            ["", "GGGG", "CCCC"],
        ],
    )
    default = signature(path, file_format="csv")
    rc = signature(path, file_format="csv", column="guide_rc_sequence")
    assert protospacer_set(path, file_format="csv", column="guide_rc_sequence") == {
        "TTTT",
        "GGGG",
    }
    assert rc != default
