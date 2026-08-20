"""Tests for file_extract.scale_h5ad."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from file_extract.metadata_sheet import parse_sample_template, sheet_wells_from_csv
from file_extract.scale_flags import validate_id_list
from file_extract.scale_h5ad import (
    default_scale_h5ad_output_name,
    extract_scale_h5ad,
    file_belongs_to_sample,
    is_scale_h5ad_key,
    parse_samples_csv,
    sample_from_filename,
)
from file_extract.scale_wells import ScaleExtractError
from tests.file_extract_helpers import FIXTURES, MockS3Client

BUCKET = "czi-cro"
RUNDATE = "proj/ORD01/processed/Run_2001-01-13/"
OTHER = "proj/ORD01/processed/Run_1999-01-01/samples/OTHER.QSR-1_anndata.h5ad"
SAMPLES_CSV = FIXTURES / "scale_samples.csv"
SHEET_CSV = FIXTURES / "scale_sample_template.csv"

SAMP01_MERGED = f"{RUNDATE}samples/SAMP-01_anndata.h5ad"
SAMP01_QSR = f"{RUNDATE}samples/SAMP-01.QSR-1_anndata.h5ad"
SAMP02_QSR = f"{RUNDATE}samples/SAMP-02.QSR-1_anndata.h5ad"
CTRL_QSR = f"{RUNDATE}samples/CTRL-01.QSR-1_anndata.h5ad"
ALLCELLS = f"{RUNDATE}samples/SAMP-01.merged.allCells.csv"
CR_H5 = f"{RUNDATE}samples/sample_filtered_feature_bc_matrix.h5"


def _sheet_text() -> str:
    return SHEET_CSV.read_text(encoding="utf-8")


def _samples_text() -> str:
    return SAMPLES_CSV.read_text(encoding="utf-8")


def test_sample_from_filename() -> None:
    assert sample_from_filename("SAMP-01.QSR-1_anndata.h5ad") == "SAMP-01"
    assert sample_from_filename("SAMP-01_anndata.h5ad") == "SAMP-01_anndata"


def test_file_belongs_to_sample() -> None:
    assert file_belongs_to_sample("SAMP-01_anndata.h5ad", "SAMP-01")
    assert file_belongs_to_sample("SAMP-01.QSR-1_anndata.h5ad", "SAMP-01")
    assert not file_belongs_to_sample("SAMP-02.QSR-1_anndata.h5ad", "SAMP-01")


def test_is_scale_h5ad_key() -> None:
    assert is_scale_h5ad_key(SAMP01_QSR, RUNDATE)
    assert not is_scale_h5ad_key(ALLCELLS, RUNDATE)
    assert not is_scale_h5ad_key(CR_H5, RUNDATE)
    assert not is_scale_h5ad_key(OTHER, RUNDATE)
    assert not is_scale_h5ad_key(f"{RUNDATE}samples.csv", RUNDATE)


def test_default_scale_h5ad_output_name() -> None:
    assert (
        default_scale_h5ad_output_name(RUNDATE) == "Run_2001-01-13_scale_h5ad_info.tsv"
    )


def test_parse_samples_csv() -> None:
    rows = parse_samples_csv(_samples_text())
    assert rows[0] == ("SAMP-01", "1A-2C")
    assert rows[-1] == ("CTRL-01", "12H")


def test_parse_sample_template_normalizes_wells() -> None:
    rows = parse_sample_template(_sheet_text())
    wells = {row["well"] for row in rows}
    assert "1A" in wells
    assert "2G" in wells
    assert "11A" not in wells
    assert sheet_wells_from_csv(_sheet_text()) == wells
    assert not any("Plate only" in (row.get("RT_index") or "") for row in rows)


def test_parse_sample_template_skips_hash_comment_rows() -> None:
    csv_text = (
        "sample_name,RT_index\n"
        "# comment,Plate only: if multiple indices, list using a comma as separator\n"
        "tissue-1A,SCALEQUANT-A1\n"
    )
    rows = parse_sample_template(csv_text)
    assert [row["well"] for row in rows] == ["1A"]
    assert sheet_wells_from_csv(csv_text) == {"1A"}


def test_parse_sample_template_splits_comma_separated_rt_index() -> None:
    csv_text = (
        "sample_name,RT_index\n"
        'tissue-pool,"SCALEQUANT-A1,SCALEQUANT-B1,SCALEQUANT-C12"\n'
    )
    rows = parse_sample_template(csv_text)
    assert [row["well"] for row in rows] == ["1A", "1B", "12C"]
    assert {row["sample_name"] for row in rows} == {"tissue-pool"}


def test_parse_sample_template_invalid_data_row_still_errors() -> None:
    with pytest.raises(ScaleExtractError, match="SCALEQUANT-Z1"):
        parse_sample_template("sample_name,RT_index\ntissue-bad,SCALEQUANT-Z1\n")


def test_validate_id_list_rejects_empty() -> None:
    with pytest.raises(ScaleExtractError):
        validate_id_list([""], flag="--wafers")
    assert validate_id_list([" 426971 ", "441969"], flag="--wafers") == [
        "426971",
        "441969",
    ]


def test_extract_scale_h5ad_writes_non_control_rows(tmp_path: Path) -> None:
    keys = [
        f"{RUNDATE}samples.csv",
        SAMP01_MERGED,
        SAMP01_QSR,
        SAMP02_QSR,
        CTRL_QSR,
        ALLCELLS,
        CR_H5,
        OTHER,
    ]
    client = MockS3Client(
        keys=keys,
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={
            SAMP01_MERGED: "crc-merged",
            SAMP01_QSR: "crc-qsr1",
            SAMP02_QSR: "crc-qsr2",
            CTRL_QSR: "crc-ctrl",
        },
    )
    out = tmp_path / "out.tsv"
    summary = extract_scale_h5ad(
        client,
        BUCKET,
        RUNDATE,
        str(out),
        metadata_gid="sheet-uuid",
        cro_orders=["ORD01"],
        wafers=["426971"],
        show_progress=False,
        sheet_csv=_sheet_text(),
    )

    assert summary.total == 3
    assert summary.crc_ok == 3
    assert any("CTRL-01" in warning for warning in summary.warnings)
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert [row["filename"] for row in rows] == [
        "SAMP-01_anndata.h5ad",
        "SAMP-01.QSR-1_anndata.h5ad",
        "SAMP-02.QSR-1_anndata.h5ad",
    ]
    assert rows[1]["sample"] == "SAMP-01"
    assert rows[1]["s3_uri"] == f"s3://{BUCKET}/{SAMP01_QSR}"
    assert rows[1]["crc64nvme_base64"] == "crc-qsr1"
    assert all(row["filename"] != "CTRL-01.QSR-1_anndata.h5ad" for row in rows)


def test_extract_scale_h5ad_zero_matches_does_not_write(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", CTRL_QSR],
        object_bodies={f"{RUNDATE}samples.csv": "sample,barcodes\nCTRL-01,12H\n"},
        crc_by_key={CTRL_QSR: "crc-ctrl"},
    )
    out = tmp_path / "empty.tsv"
    summary = extract_scale_h5ad(
        client,
        BUCKET,
        RUNDATE,
        str(out),
        metadata_gid="sheet-uuid",
        cro_orders=["ORD01"],
        wafers=["426971"],
        show_progress=False,
        sheet_csv=_sheet_text(),
    )
    assert summary.total == 0
    assert not out.exists()


def test_extract_scale_h5ad_uses_fetch_sheet(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1"},
    )
    fetched: list[str] = []

    def fetch_sheet(sheet_id: str) -> str:
        fetched.append(sheet_id)
        return _sheet_text()

    out = tmp_path / "fetched.tsv"
    extract_scale_h5ad(
        client,
        BUCKET,
        RUNDATE,
        str(out),
        metadata_gid="sheet-uuid",
        cro_orders=["ORD01"],
        wafers=["426971"],
        show_progress=False,
        fetch_sheet=fetch_sheet,
    )
    assert fetched == ["sheet-uuid"]
    assert out.exists()


def test_extract_scale_h5ad_missing_samples_csv() -> None:
    client = MockS3Client(keys=[SAMP01_QSR], crc_by_key={SAMP01_QSR: "crc"})
    with pytest.raises(ScaleExtractError, match="samples.csv"):
        extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            "unused.tsv",
            metadata_gid="sheet-uuid",
            cro_orders=["ORD01"],
            wafers=["426971"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )


def test_metadata_sheet_parse_requires_rt_index() -> None:
    with pytest.raises(ScaleExtractError, match="RT_index"):
        parse_sample_template("sample_name\nfoo\n")
