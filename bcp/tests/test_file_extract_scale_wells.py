"""Tests for file_extract.scale_wells."""

from __future__ import annotations

import pytest

from file_extract.scale_wells import (
    WellOwnershipError,
    correlate_samples,
    expand_barcodes,
    normalize_rt_index,
    parse_rt_index_cell,
    parse_well,
)


def test_normalize_rt_index_strips_and_flips() -> None:
    assert normalize_rt_index("SCALEQUANT-A11") == "11A"
    assert normalize_rt_index("SCALEQUANT-A1") == "1A"
    assert normalize_rt_index("A1") == "1A"
    assert normalize_rt_index("11A") == "11A"


def test_normalize_rt_index_rejects_bad_values() -> None:
    from file_extract.scale_wells import ScaleExtractError

    with pytest.raises(ScaleExtractError):
        normalize_rt_index("SCALEQUANT-Z1")
    with pytest.raises(ScaleExtractError):
        normalize_rt_index("SCALEQUANT-A13")


def test_parse_rt_index_cell_splits_commas() -> None:
    assert parse_rt_index_cell("SCALEQUANT-A1,SCALEQUANT-B1,SCALEQUANT-H11") == (
        "1A",
        "1B",
        "11H",
    )
    assert parse_rt_index_cell("SCALEQUANT-A1, SCALEQUANT-B1") == ("1A", "1B")
    assert parse_rt_index_cell("") == ()


def test_parse_rt_index_cell_rejects_invalid_token() -> None:
    from file_extract.scale_wells import ScaleExtractError

    with pytest.raises(ScaleExtractError, match="SCALEQUANT-Z1"):
        parse_rt_index_cell("SCALEQUANT-A1,SCALEQUANT-Z1")


def test_expand_barcodes_1a_to_2c() -> None:
    wells = expand_barcodes("1A-2C")
    assert wells == (
        "1A",
        "1B",
        "1C",
        "1D",
        "1E",
        "1F",
        "1G",
        "1H",
        "2A",
        "2B",
        "2C",
    )


def test_expand_barcodes_same_column_and_semicolon() -> None:
    assert expand_barcodes("2A-2G") == ("2A", "2B", "2C", "2D", "2E", "2F", "2G")
    assert expand_barcodes("3A-3G;8H") == (
        "3A",
        "3B",
        "3C",
        "3D",
        "3E",
        "3F",
        "3G",
        "8H",
    )
    assert expand_barcodes("1A") == ("1A",)
    assert expand_barcodes("12H") == ("12H",)


def test_parse_well_rejects_row_col_form() -> None:
    from file_extract.scale_wells import ScaleExtractError

    with pytest.raises(ScaleExtractError):
        parse_well("A1")


def test_correlate_samples_paired_vs_control() -> None:
    result = correlate_samples(
        [
            ("SAMP-01", "1A-2C"),
            ("SAMP-02", "3A-3G"),
            ("CTRL-01", "12H"),
        ],
        sheet_wells={
            "1A",
            "1B",
            "2A",
            "2B",
            "2C",
            "3A",
            "3B",
            "3C",
            "3D",
            "3E",
            "3F",
            "3G",
        },
    )
    assert result.paired == ("SAMP-01", "SAMP-02")
    assert [item.sample for item in result.controls] == ["CTRL-01"]
    assert result.controls[0].barcodes == "12H"


def test_correlate_samples_duplicate_well_is_error() -> None:
    with pytest.raises(WellOwnershipError, match="1A"):
        correlate_samples(
            [("SAMP-01", "1A"), ("SAMP-02", "1A-1B")],
            sheet_wells={"1A"},
        )
