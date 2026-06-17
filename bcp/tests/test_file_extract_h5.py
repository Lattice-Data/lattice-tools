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
    map_feature_counts,
    parse_metrics_cells_from_text,
)
from tests.file_extract_helpers import FIXTURES, MockS3Client

LATTICE_FEATURE_TYPES = frozenset(
    {
        "antibody capture",
        "gene",
        "guide capture",
        "peak",
        "transcription factor",
    }
)

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
    introspect = h5_columns(do_introspect=True, do_genome=False, do_metrics=False)
    assert introspect == [
        "library",
        "sample",
        "s3_uri",
        "size_bytes",
        "crc64nvme_base64",
        "observation_count",
        "feature_counts",
        "feature_count_total",
        "unmapped_feature_types",
        "crc_error",
        "h5_error",
    ]
    full = h5_columns(do_introspect=True, do_genome=True, do_metrics=True)
    assert "gene_counts_by_genome" in full
    assert "metrics_cells_match" in full


def _assert_lattice_feature_counts(
    lattice_fc: list[dict[str, int | str]],
) -> None:
    for item in lattice_fc:
        assert item["feature_type"] in LATTICE_FEATURE_TYPES


@pytest.mark.parametrize(
    ("raw_counts", "expected_fc", "expected_unmapped"),
    [
        (
            {"Gene Expression": 18129},
            [{"feature_type": "gene", "feature_count": 18129}],
            {},
        ),
        (
            {"Gene Expression": 18129, "CRISPR Guide Capture": 1842},
            [
                {"feature_type": "gene", "feature_count": 18129},
                {"feature_type": "guide capture", "feature_count": 1842},
            ],
            {},
        ),
        (
            {"Gene Expression": 100, "Multiplexing Capture": 4},
            [{"feature_type": "gene", "feature_count": 100}],
            {"Multiplexing Capture": 4},
        ),
        (
            {},
            [],
            {},
        ),
    ],
)
def test_map_feature_counts(
    raw_counts: dict[str, int],
    expected_fc: list[dict[str, int | str]],
    expected_unmapped: dict[str, int],
) -> None:
    lattice_fc, unmapped = map_feature_counts(raw_counts)
    assert lattice_fc == expected_fc
    assert unmapped == expected_unmapped
    _assert_lattice_feature_counts(lattice_fc)


def test_map_feature_counts_strips_whitespace_for_lookup() -> None:
    lattice_fc, unmapped = map_feature_counts({"  Gene Expression  ": 42})
    assert lattice_fc == [{"feature_type": "gene", "feature_count": 42}]
    assert unmapped == {}


def test_map_feature_counts_sums_shared_lattice_type(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from file_extract import h5

    monkeypatch.setattr(
        h5,
        "CR_TO_LATTICE_FEATURE_TYPE",
        {"Type A": "gene", "Type B": "gene"},
    )
    lattice_fc, unmapped = h5.map_feature_counts({"Type A": 10, "Type B": 5})
    assert lattice_fc == [{"feature_type": "gene", "feature_count": 15}]
    assert unmapped == {}


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
