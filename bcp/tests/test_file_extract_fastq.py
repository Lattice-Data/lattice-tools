"""Tests for file_extract.fastq."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from file_extract.fastq import (
    default_fastq_output_name,
    extract_fastq,
    is_target_file,
    parse_read_lane,
    r1_r2_mismatch_warning,
)
from tests.file_extract_helpers import FIXTURES, MockS3Client

PREFIX = "test-order/NVUS0000000000-01/"
BUCKET = "example-bucket"
R1_KEY = f"{PREFIX}GROUP1/raw/sample_L001_R1_001.fastq.gz"
R2_KEY = f"{PREFIX}GROUP1/raw/sample_L001_R2_001.fastq.gz"


@pytest.mark.parametrize(
    ("key", "require_raw", "expected"),
    [
        (R1_KEY, True, True),
        (R1_KEY.replace("/raw/", "/processed/"), True, False),
        (f"{PREFIX}GROUP1/raw/sample_sample.fastq.gz", True, False),
        (f"{PREFIX}GROUP1/raw/unmatched.fastq.gz", True, False),
        (R1_KEY.replace("/raw/", "/other/"), False, True),
    ],
)
def test_is_target_file(key: str, require_raw: bool, expected: bool) -> None:
    assert is_target_file(key, require_raw=require_raw) is expected


@pytest.mark.parametrize(
    ("fname", "read", "lane"),
    [
        ("sample_L001_R1_001.fastq.gz", "R1", "001"),
        ("sample_L002_R2_003.fastq.gz", "R2", "002"),
        ("sample_I1_001.fastq.gz", "I1", ""),
        ("no_lane_R1.fastq.gz", "", ""),
    ],
)
def test_parse_read_lane(fname: str, read: str, lane: str) -> None:
    assert parse_read_lane(fname) == (read, lane)


def test_default_fastq_output_name() -> None:
    assert default_fastq_output_name("order/prefix/") == "prefix_fastq_info.tsv"


def test_r1_r2_mismatch_warning() -> None:
    assert r1_r2_mismatch_warning({"R1": 2, "R2": 2}) is None
    msg = r1_r2_mismatch_warning({"R1": 3, "R2": 2})
    assert msg is not None
    assert "R1 (3)" in msg


def test_extract_fastq_integration(tmp_path: Path) -> None:
    meta_body = (FIXTURES / "fastq_metadata.json").read_text(encoding="utf-8")
    keys = [R1_KEY, R2_KEY]
    sizes = {k: 1000 for k in keys}
    crc = {k: f"crc-{k[-6:]}" for k in keys}
    bodies = {k + "-metadata.json": meta_body for k in keys}
    client = MockS3Client(
        keys=keys,
        sizes=sizes,
        crc_by_key=crc,
        object_bodies=bodies,
    )

    out = tmp_path / "out.tsv"
    summary = extract_fastq(
        client,
        BUCKET,
        PREFIX,
        str(out),
        workers=2,
        show_progress=False,
        inline=True,
    )

    assert summary.total == 2
    assert summary.crc_ok == 2
    assert summary.enrichment_ok == 2
    assert summary.read_tally["R1"] == 1
    assert summary.read_tally["R2"] == 1

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0] == [
        "filename",
        "s3_uri",
        "read",
        "lane",
        "size_bytes",
        "crc64nvme_base64",
        "read_count",
        "crc_error",
        "metadata_error",
    ]
    assert len(rows) == 3
    data_rows = rows[1:]
    reads = {row[2] for row in data_rows}
    assert reads == {"R1", "R2"}
    assert all(row[6] == "1234567" for row in data_rows)


def test_extract_fastq_missing_metadata(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[R1_KEY],
        sizes={R1_KEY: 100},
        crc_by_key={R1_KEY: "abc"},
        object_bodies={},
    )
    out = tmp_path / "out.tsv"
    summary = extract_fastq(
        client,
        BUCKET,
        PREFIX,
        str(out),
        show_progress=False,
        inline=True,
    )
    assert summary.crc_ok == 1
    assert summary.enrichment_ok == 0
    assert len(summary.failures) == 1
