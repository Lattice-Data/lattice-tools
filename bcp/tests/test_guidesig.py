"""Unit and regression tests for guidesig."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from guidesig import GuideSigError, protospacer_set, signature
from tests.guidesig_helpers import FIXTURES, write_delimited

GOLDEN = {
    "lib_one.tsv": "gsig1:set:n=60:62d08f4954bf50b7249c6277e2808b16",
    "lib_one_prepended_g.tsv": "gsig1:set:n=60:8a3af4c57186c0766609fd1841136e7c",
    "lib_two_layout_a.tsv": "gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78",
    "lib_two_layout_b.tsv": "gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78",
}


def _write_tsv(path: Path, rows: list[list[str]], *, bom: bool = True) -> Path:
    """Write a CRLF TSV, optionally with a UTF-8 BOM."""
    return write_delimited(path, rows, file_format="tsv", bom=bom)


def test_row_order_permuted_identical(tmp_path: Path) -> None:
    header = ["property", "guide_id", "guide_protospacer"]
    rows_a = [
        header,
        ["", "g1", "AAAA"],
        ["", "g2", "CCCC"],
        ["", "g3", "GGGG"],
    ]
    rows_b = [
        header,
        ["", "g3", "GGGG"],
        ["", "g1", "AAAA"],
        ["", "g2", "CCCC"],
    ]
    a = _write_tsv(tmp_path / "a.tsv", rows_a)
    b = _write_tsv(tmp_path / "b.tsv", rows_b)
    assert signature(a, file_format="tsv") == signature(b, file_format="tsv")


def test_non_protospacer_columns_irrelevant(tmp_path: Path) -> None:
    a = _write_tsv(
        tmp_path / "a.tsv",
        [
            ["property", "guide_id", "guide_protospacer", "gene"],
            ["", "g1", "AAAA", "GENE1"],
            ["", "g2", "CCCC", "GENE2"],
        ],
    )
    b = _write_tsv(
        tmp_path / "b.tsv",
        [
            ["property", "guide_id", "guide_protospacer", "gene", "pam"],
            ["", "other1", "AAAA", "OTHER", "NGG"],
            ["", "other2", "CCCC", "DIFFERENT", "NGG"],
        ],
    )
    assert signature(a, file_format="tsv") == signature(b, file_format="tsv")


def test_column_reorder_and_add_remove_identical(tmp_path: Path) -> None:
    a = _write_tsv(
        tmp_path / "a.tsv",
        [
            ["property", "guide_protospacer", "guide_id"],
            ["", "AAAA", "g1"],
            ["", "CCCC", "g2"],
        ],
    )
    b = _write_tsv(
        tmp_path / "b.tsv",
        [
            ["guide_id", "extra", "guide_protospacer", "property"],
            ["g1", "x", "AAAA", ""],
            ["g2", "y", "CCCC", ""],
        ],
    )
    assert signature(a, file_format="tsv") == signature(b, file_format="tsv")


def test_duplicates_collapse_to_set(tmp_path: Path) -> None:
    once = _write_tsv(
        tmp_path / "once.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", "CCCC"],
        ],
    )
    many = _write_tsv(
        tmp_path / "many.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", "AAAA"],
            ["", "CCCC"],
            ["", "CCCC"],
            ["", "CCCC"],
        ],
    )
    assert signature(once, file_format="tsv") == signature(many, file_format="tsv")
    assert signature(once, file_format="tsv").startswith("gsig1:set:n=2:")


def test_lowercase_identical_to_uppercase(tmp_path: Path) -> None:
    upper = _write_tsv(
        tmp_path / "upper.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    lower = _write_tsv(
        tmp_path / "lower.tsv",
        [["property", "guide_protospacer"], ["", "aaaa"], ["", "cccc"]],
    )
    assert signature(upper, file_format="tsv") == signature(lower, file_format="tsv")


def test_whitespace_and_stray_cr_stripped(tmp_path: Path) -> None:
    clean = _write_tsv(
        tmp_path / "clean.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    # Quote the field so the embedded CR is data, not a line terminator.
    padded = tmp_path / "padded.tsv"
    padded.write_bytes(
        b'\xef\xbb\xbfproperty\tguide_protospacer\r\n\t  AAAA  \r\n\t"CCCC\r"\r\n'
    )
    assert signature(clean, file_format="tsv") == signature(padded, file_format="tsv")


def test_signature_stable_across_calls_and_reread(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "stable.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    first = signature(path, file_format="tsv")
    second = signature(path, file_format="tsv")
    third = signature(Path(str(path)), file_format="tsv")
    assert first == second == third


def test_one_sequence_changed_different(tmp_path: Path) -> None:
    a = _write_tsv(
        tmp_path / "a.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    b = _write_tsv(
        tmp_path / "b.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "TTTT"]],
    )
    assert signature(a, file_format="tsv") != signature(b, file_format="tsv")


def test_one_sequence_added_increments_n(tmp_path: Path) -> None:
    a = _write_tsv(
        tmp_path / "a.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    b = _write_tsv(
        tmp_path / "b.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", "CCCC"],
            ["", "GGGG"],
        ],
    )
    sig_a = signature(a, file_format="tsv")
    sig_b = signature(b, file_format="tsv")
    assert sig_a != sig_b
    assert ":n=2:" in sig_a
    assert ":n=3:" in sig_b


def test_prepended_g_different_signature(tmp_path: Path) -> None:
    a = _write_tsv(
        tmp_path / "a.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    b = _write_tsv(
        tmp_path / "b.tsv",
        [["property", "guide_protospacer"], ["", "GAAAA"], ["", "GCCCC"]],
    )
    assert signature(a, file_format="tsv") != signature(b, file_format="tsv")
    assert ":n=2:" in signature(a, file_format="tsv")
    assert ":n=2:" in signature(b, file_format="tsv")


def test_sequence_and_reverse_complement_both_retained(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "rc.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAC"],
            ["", "GTTT"],  # reverse complement of AAAC
        ],
    )
    sequences = protospacer_set(path, file_format="tsv")
    assert sequences == {"AAAC", "GTTT"}
    assert signature(path, file_format="tsv").startswith("gsig1:set:n=2:")


def test_comment_and_example_rows_filtered(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "comments.tsv",
        [
            ["property", "guide_protospacer"],
            ["#comment", "TTTTTTTTTTTTTTTTTTT"],
            ["#example", "GGGGGGGGGGGGGGGGGGG"],
            ["", "AAAA"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}


def test_only_example_row_raises_empty_set(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "only_example.tsv",
        [
            ["property", "guide_protospacer"],
            ["#example", "GGGGGGGGGGGGGGGGGGG"],
        ],
    )
    with pytest.raises(GuideSigError, match="zero usable sequences"):
        signature(path, file_format="tsv")


def test_hash_space_after_hash_still_filtered(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "hash_space.tsv",
        [
            ["property", "guide_protospacer"],
            ["# example", "TTTTTTTTTTTTTTTTTTT"],
            ["", "AAAA"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA"}


def test_empty_first_column_retained(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "empty_first.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", "CCCC"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCC"}


def test_blank_lines_skipped(tmp_path: Path) -> None:
    """Trailing/intervening blank lines must not raise or affect the set."""
    path = tmp_path / "blank_lines.tsv"
    path.write_bytes(
        b"\xef\xbb\xbfproperty\tguide_protospacer\r\n\tAAAA\r\n\r\n\tCCCC\r\n\r\n"
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCC"}
    assert signature(path, file_format="tsv").startswith("gsig1:set:n=2:")


def test_bom_and_header_lookup(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "bom.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"]],
        bom=True,
    )
    assert path.read_bytes().startswith(b"\xef\xbb\xbf")
    assert signature(path, file_format="tsv").startswith("gsig1:set:n=1:")


def test_crlf_does_not_leave_cr_in_set(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "crlf.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    assert path.read_bytes().count(b"\r\n") >= 2
    assert all("\r" not in seq for seq in protospacer_set(path, file_format="tsv"))


def test_empty_protospacer_skipped(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "empty_cell.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", ""],
            ["", "CCCC"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCC"}


def test_ragged_row_raises_with_row_number(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "ragged.tsv",
        [
            ["property", "guide_id", "guide_protospacer"],
            ["", "g1", "AAAA"],
            ["", "g2"],  # missing protospacer cell
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="tsv")


def test_ragged_row_with_protospacer_still_raises(tmp_path: Path) -> None:
    """Dropped trailing tabs are a ragged row even when the sequence is present."""
    path = _write_tsv(
        tmp_path / "ragged_trailing_columns.tsv",
        [
            ["property", "guide_protospacer", "guide_id", "gene", "pam"],
            ["", "AAAA", "g1", "GENE1", "NGG"],
            ["", "CCCC"],  # protospacer present, trailing columns were dropped
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="tsv")


def test_multiline_quoted_cell_parses(tmp_path: Path) -> None:
    """A quoted cell holding an embedded newline is one record, not two rows."""
    path = _write_tsv(
        tmp_path / "multiline.tsv",
        [
            ["property", "guide_id", "guide_protospacer"],
            ["", "pasted from\na spreadsheet", "AAAA"],
            ["", "g2", "CCCC"],
        ],
    )
    assert protospacer_set(path, file_format="tsv") == {"AAAA", "CCCC"}


def test_multiline_quoted_cell_keeps_row_numbers_physical(tmp_path: Path) -> None:
    """An earlier embedded newline must not shift the reported row number.

    The ragged row sits on physical line 4 because the preceding record spans
    lines 2 and 3, which is the line a user opening the file would look at.
    """
    path = _write_tsv(
        tmp_path / "multiline_ragged.tsv",
        [
            ["property", "guide_id", "guide_protospacer"],
            ["", "pasted from\na spreadsheet", "AAAA"],
            ["", "g2"],  # missing protospacer cell, on physical line 4
        ],
    )
    assert path.read_bytes().count(b"\n") == 4
    with pytest.raises(GuideSigError, match="row 4"):
        signature(path, file_format="tsv")


def test_multiline_offending_row_reports_its_first_line(tmp_path: Path) -> None:
    """A bad record that itself spans lines is reported at the line it starts on.

    Both records span two lines, so the offending one starts on physical line 4.
    """
    path = _write_tsv(
        tmp_path / "multiline_bad.tsv",
        [
            ["property", "guide_id", "guide_protospacer"],
            ["", "pasted from\na spreadsheet", "AAAA"],
            ["", "also pasted\nfrom a spreadsheet", "AAAN"],
        ],
    )
    assert path.read_bytes().count(b"\n") == 5
    with pytest.raises(GuideSigError, match="row 4"):
        signature(path, file_format="tsv")


@pytest.mark.parametrize("bad", ["N", "U", "0", "-"])
def test_non_acgt_raises_with_row_number(tmp_path: Path, bad: str) -> None:
    path = _write_tsv(
        tmp_path / "bad.tsv",
        [
            ["property", "guide_protospacer"],
            ["", "AAAA"],
            ["", f"AAA{bad}"],
        ],
    )
    with pytest.raises(GuideSigError, match="row 3"):
        signature(path, file_format="tsv")


def test_missing_column_raises(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "missing.tsv",
        [["property", "guide_id"], ["", "g1"]],
    )
    with pytest.raises(GuideSigError, match="absent from header"):
        signature(path, file_format="tsv")


def test_non_utf8_file_raises_guidesig_error(tmp_path: Path) -> None:
    """A Latin-1/Windows-1252 export must raise GuideSigError, not a raw decode error."""
    path = tmp_path / "latin1.tsv"
    # 0xE9 is 'é' in Latin-1 and is an invalid UTF-8 start byte here.
    path.write_bytes(b"property\tguide_protospacer\r\n\tGENE\xe9\tAAAA\r\n")
    with pytest.raises(GuideSigError, match="not valid UTF-8"):
        signature(path, file_format="tsv")


def test_serialization_pinning_join_and_no_trailing_newline(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "pin.tsv",
        [["property", "guide_protospacer"], ["", "AAAA"], ["", "CCCC"]],
    )
    expected_digest = hashlib.sha256(b"AAAA\nCCCC").hexdigest()[:32]
    assert signature(path, file_format="tsv") == f"gsig1:set:n=2:{expected_digest}"


@pytest.mark.parametrize("filename,expected", list(GOLDEN.items()))
def test_golden_digests(filename: str, expected: str) -> None:
    assert signature(FIXTURES / filename, file_format="tsv") == expected


def test_lib_two_layouts_must_match() -> None:
    """Layout, order, case, and padding must not affect the signature."""
    a = signature(FIXTURES / "lib_two_layout_a.tsv", file_format="tsv")
    b = signature(FIXTURES / "lib_two_layout_b.tsv", file_format="tsv")
    assert a == b == GOLDEN["lib_two_layout_a.tsv"]


def test_lib_one_and_prepended_g_must_not_match() -> None:
    """Equal n with a 5' G difference must not match (no length normalization)."""
    one = signature(FIXTURES / "lib_one.tsv", file_format="tsv")
    prepended = signature(FIXTURES / "lib_one_prepended_g.tsv", file_format="tsv")
    assert one != prepended
    assert ":n=60:" in one
    assert ":n=60:" in prepended


def test_alternate_column_differs_from_default() -> None:
    """--column must be honored; RC column is not interchangeable."""
    default = signature(FIXTURES / "lib_two_layout_b.tsv", file_format="tsv")
    rc = signature(
        FIXTURES / "lib_two_layout_b.tsv",
        file_format="tsv",
        column="guide_rc_sequence",
    )
    assert rc == "gsig1:set:n=60:429cecb0d8b27564cddc32dfde6972ff"
    assert rc != default


def test_file_format_is_required() -> None:
    """The format must always be stated; there is no default."""
    with pytest.raises(TypeError):
        signature(FIXTURES / "lib_one.tsv")  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        protospacer_set(FIXTURES / "lib_one.tsv")  # type: ignore[call-arg]


def test_unknown_file_format_raises() -> None:
    with pytest.raises(GuideSigError, match="unknown file format"):
        signature(FIXTURES / "lib_one.tsv", file_format="xlsx")
