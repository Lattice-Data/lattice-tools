"""Tests for file_extract.h5."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from file_extract.h5 import (
    default_h5_output_name,
    extract_h5,
    extract_library,
    extract_sample_name,
    h5_columns,
    parse_metrics_cells_from_text,
)
from tests.file_extract_helpers import FIXTURES, MockS3Client

PREFIX = "proj/lib/processed/outs/per_sample_outs/"
BUCKET = "example-bucket"
SAMPLE = "sampleA"
H5_KEY = (
    f"proj/lib/processed/outs/per_sample_outs/{SAMPLE}/"
    "sample_filtered_feature_bc_matrix.h5"
)
METRICS_KEY = f"proj/lib/processed/outs/per_sample_outs/{SAMPLE}/metrics_summary.csv"


def test_extract_sample_name() -> None:
    assert extract_sample_name(H5_KEY, PREFIX) == SAMPLE
    assert extract_sample_name("other/key", PREFIX) == ""


def test_extract_library() -> None:
    assert extract_library(H5_KEY) == "lib"
    assert extract_library("bucket/only/key.h5") == ""


def test_default_h5_output_name() -> None:
    assert default_h5_output_name("a/b/per_sample_outs/") == "b_h5_info.tsv"


def test_h5_columns_variants() -> None:
    base = h5_columns(do_introspect=False, do_genome=False, do_metrics=False)
    assert base == [
        "library",
        "sample",
        "s3_uri",
        "size_bytes",
        "crc64nvme_base64",
        "crc_error",
    ]
    full = h5_columns(do_introspect=True, do_genome=True, do_metrics=True)
    assert "gene_counts_by_genome" in full
    assert "metrics_cells_match" in full


def test_parse_metrics_long_format() -> None:
    text = (FIXTURES / "metrics_summary_long.csv").read_text(encoding="utf-8")
    assert parse_metrics_cells_from_text(text) == 1234


def test_parse_metrics_wide_format() -> None:
    text = (FIXTURES / "metrics_summary_wide.csv").read_text(encoding="utf-8")
    assert parse_metrics_cells_from_text(text) == 1234


def test_parse_metrics_missing_cells() -> None:
    with pytest.raises(RuntimeError, match="cell-count"):
        parse_metrics_cells_from_text("col1,col2\n1,2\n")


def test_extract_h5_checksum_only(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[H5_KEY],
        sizes={H5_KEY: 5000},
        crc_by_key={H5_KEY: "crc-h5"},
    )
    out = tmp_path / "h5.tsv"
    summary = extract_h5(
        client,
        BUCKET,
        PREFIX,
        str(out),
        do_introspect=False,
        show_progress=False,
        workers=1,
    )
    assert summary.total == 1
    assert summary.crc_ok == 1

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[1][4] == "crc-h5"
    assert rows[1][0] == "lib"
    assert rows[1][1] == SAMPLE


def test_extract_h5_with_metrics(tmp_path: Path) -> None:
    metrics_text = (FIXTURES / "metrics_summary_long.csv").read_text(encoding="utf-8")
    client = MockS3Client(
        keys=[H5_KEY],
        sizes={H5_KEY: 5000},
        crc_by_key={H5_KEY: "crc-h5"},
        object_bodies={METRICS_KEY: metrics_text},
    )
    out = tmp_path / "h5_metrics.tsv"
    summary = extract_h5(
        client,
        BUCKET,
        PREFIX,
        str(out),
        do_introspect=False,
        do_metrics=True,
        show_progress=False,
        workers=1,
    )
    assert summary.total == 1
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    metrics_idx = header.index("metrics_cells")
    assert rows[1][metrics_idx] == "1234"
