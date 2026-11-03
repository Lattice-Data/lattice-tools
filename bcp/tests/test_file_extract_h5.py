"""Tests for file_extract.h5."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from file_extract.h5 import (
    default_h5_output_name,
    empty_raw_fastq_warning,
    extract_h5,
    extract_library,
    extract_sample_name,
    h5_columns,
    h5_worker_ceiling,
    map_feature_counts,
    missing_processed_warning,
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
H5_BODY = b"\x89HDF\r\n\x1a\n" + bytes(range(256)) * 4


def _md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


def test_extract_sample_name() -> None:
    full_key = (
        "proj/AN00000001/cohort1_batch1_1/processed/cellranger/2024-01-01/"
        "outs/per_sample_outs/sampleA/count/sample_filtered_feature_bc_matrix.h5"
    )
    assert extract_sample_name(H5_KEY) == SAMPLE
    assert extract_sample_name(full_key) == SAMPLE
    assert extract_sample_name("other/key") == ""
    assert (
        extract_sample_name(
            "proj/lib/processed/outs/per_sample_outs/"
            "sample_filtered_feature_bc_matrix.h5"
        )
        == ""
    )


def test_extract_library() -> None:
    assert extract_library(H5_KEY) == "lib"
    assert extract_library("bucket/only/key.h5") == ""


def test_default_h5_output_name() -> None:
    assert (
        default_h5_output_name("a/b/per_sample_outs/") == "per_sample_outs_h5_info.tsv"
    )
    assert default_h5_output_name("proj/AN00000001") == "AN00000001_h5_info.tsv"
    assert default_h5_output_name("") == "output_h5_info.tsv"


def test_h5_worker_ceiling() -> None:
    assert h5_worker_ceiling(do_introspect=True, do_md5=True) == 8
    assert h5_worker_ceiling(do_introspect=True, do_md5=False) == 8
    # Each md5 holds a streaming GET open for a whole h5, so the CRC-only
    # fan-out of 64 applies only when nothing reads file content.
    assert h5_worker_ceiling(do_introspect=False, do_md5=True) == 8
    assert h5_worker_ceiling(do_introspect=False, do_md5=False) == 64
    assert h5_worker_ceiling(do_introspect=True, do_md5=True, workers=3) == 3
    with pytest.raises(ValueError, match="positive integer"):
        h5_worker_ceiling(do_introspect=True, do_md5=True, workers=0)
    with pytest.raises(ValueError, match="positive integer"):
        h5_worker_ceiling(do_introspect=False, do_md5=False, workers=-1)


def test_extract_h5_forwards_pool_size(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[H5_KEY],
        sizes={H5_KEY: 5000},
        crc_by_key={H5_KEY: "crc-h5"},
    )
    out = tmp_path / "h5.tsv"
    with patch(
        "file_extract.h5.introspect_h5", return_value=(10, {}, None)
    ) as mock_intro:
        summary = extract_h5(
            client,
            BUCKET,
            PREFIX,
            str(out),
            lab="example-lab",
            do_introspect=True,
            show_progress=False,
        )
    assert summary.total == 1
    mock_intro.assert_called_once_with(BUCKET, H5_KEY, pool_size=8)


def test_h5_columns_variants() -> None:
    base = h5_columns(
        do_introspect=False, do_genome=False, do_metrics=False, do_md5=False
    )
    assert base == [
        "library",
        "sample",
        "s3_uri",
        "size_bytes",
        "crc64nvme_base64",
        "derived_from",
        "crc_error",
    ]
    hashed = h5_columns(
        do_introspect=False, do_genome=False, do_metrics=False, do_md5=True
    )
    assert hashed == [
        "library",
        "sample",
        "s3_uri",
        "size_bytes",
        "crc64nvme_base64",
        "md5sum",
        "derived_from",
        "crc_error",
        "md5_error",
    ]
    introspect = h5_columns(
        do_introspect=True, do_genome=False, do_metrics=False, do_md5=True
    )
    assert introspect == [
        "library",
        "sample",
        "s3_uri",
        "size_bytes",
        "crc64nvme_base64",
        "md5sum",
        "derived_from",
        "observation_count",
        "feature_counts",
        "feature_count_total",
        "unmapped_feature_types",
        "crc_error",
        "md5_error",
        "h5_error",
    ]
    full = h5_columns(do_introspect=True, do_genome=True, do_metrics=True, do_md5=True)
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
        object_bodies={H5_KEY: H5_BODY},
    )
    out = tmp_path / "h5.tsv"
    summary = extract_h5(
        client,
        BUCKET,
        PREFIX,
        str(out),
        lab="example-lab",
        do_introspect=False,
        show_progress=False,
        workers=1,
    )
    assert summary.total == 1
    assert summary.crc_ok == 1
    assert summary.md5_ok == 1
    assert summary.failures == []
    assert summary.warnings == [empty_raw_fastq_warning("proj/lib/raw/")]

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["crc64nvme_base64"] == "crc-h5"
    assert rows[0]["md5sum"] == _md5(H5_BODY)
    assert rows[0]["md5_error"] == ""
    assert client.get_object_calls == [
        {"Bucket": BUCKET, "Key": H5_KEY, "ChecksumMode": "ENABLED"}
    ]
    assert rows[0]["library"] == "lib"
    assert rows[0]["sample"] == SAMPLE
    assert json.loads(rows[0]["derived_from"]) == []


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
        lab="example-lab",
        do_introspect=False,
        do_metrics=True,
        do_md5=False,
        show_progress=False,
        workers=1,
    )
    assert summary.total == 1
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["metrics_cells"] == "1234"
    assert json.loads(rows[0]["derived_from"]) == []


def _h5_key(parent: str, sample: str) -> str:
    return (
        f"{parent}/processed/cellranger/2024-01-01/outs/"
        f"per_sample_outs/{sample}/sample_filtered_feature_bc_matrix.h5"
    )


def test_extract_h5_derived_from_is_shared_per_library(tmp_path: Path) -> None:
    lib_a = "proj/order/libA"
    lib_b = "proj/order/libB"
    lib_a_other = "proj/other/libA"
    h5_a1 = _h5_key(lib_a, "sampleA")
    h5_a2 = _h5_key(lib_a, "sampleB")
    h5_b = _h5_key(lib_b, "sampleC")
    h5_other = _h5_key(lib_a_other, "sampleD")
    bare = "proj/order/no_processed/sample_filtered_feature_bc_matrix.h5"
    h5_root = _h5_key("", "sampleRoot").removeprefix("/")
    keys = [
        h5_a1,
        h5_a2,
        h5_b,
        h5_other,
        bare,
        h5_root,
        f"{lib_a}/raw/nested/a_R2.fastq.gz",
        f"{lib_a}/raw/a_R1.fastq.gz",
        f"{lib_a}/raw/unmatched_R1.fastq.gz",
        f"{lib_a}/raw/foo_sample_R1.fastq.gz",
        f"{lib_b}/raw/b_R1.fastq.gz",
        f"{lib_a_other}/raw/other_R1.fastq.gz",
        "raw/a.fastq.gz",
    ]
    client = MockS3Client(
        keys=keys,
        crc_by_key={key: "crc" for key in keys if key.endswith(".h5")},
    )
    out = tmp_path / "grouped.tsv"
    summary = extract_h5(
        client,
        BUCKET,
        "",
        str(out),
        lab="/labs/example-lab/",
        do_introspect=False,
        do_md5=False,
        show_progress=False,
        workers=1,
    )
    assert summary.total == 6
    with out.open(encoding="utf-8") as fh:
        rows = {row["s3_uri"]: row for row in csv.DictReader(fh, delimiter="\t")}

    shared = ["example-lab:a_R1.fastq.gz", "example-lab:a_R2.fastq.gz"]
    assert json.loads(rows[f"s3://{BUCKET}/{h5_a1}"]["derived_from"]) == shared
    assert json.loads(rows[f"s3://{BUCKET}/{h5_a2}"]["derived_from"]) == shared
    assert rows[f"s3://{BUCKET}/{h5_a1}"]["sample"] == "sampleA"
    assert rows[f"s3://{BUCKET}/{h5_a1}"]["library"] == "libA"
    assert json.loads(rows[f"s3://{BUCKET}/{h5_b}"]["derived_from"]) == [
        "example-lab:b_R1.fastq.gz"
    ]
    assert json.loads(rows[f"s3://{BUCKET}/{h5_other}"]["derived_from"]) == [
        "example-lab:other_R1.fastq.gz"
    ]
    assert json.loads(rows[f"s3://{BUCKET}/{bare}"]["derived_from"]) == []
    assert rows[f"s3://{BUCKET}/{bare}"]["sample"] == ""
    assert missing_processed_warning(bare) in summary.warnings
    assert json.loads(rows[f"s3://{BUCKET}/{h5_root}"]["derived_from"]) == []
    assert missing_processed_warning(h5_root) in summary.warnings
    assert empty_raw_fastq_warning("raw/") not in summary.warnings


# --- md5sum ------------------------------------------------------------------

H5_KEY_B = H5_KEY.replace(SAMPLE, "sampleB")
H5_BODY_B = b"second matrix" * 10


def _h5_md5_run(
    tmp_path: Path, client: MockS3Client, **kwargs: Any
) -> tuple[Any, dict[str, dict[str, str]]]:
    out = tmp_path / "md5.tsv"
    summary = extract_h5(
        client,
        BUCKET,
        PREFIX,
        str(out),
        lab="example-lab",
        do_introspect=False,
        show_progress=False,
        **kwargs,
    )
    with out.open(encoding="utf-8") as fh:
        rows = {row["sample"]: row for row in csv.DictReader(fh, delimiter="\t")}
    return summary, rows


def _two_h5_client(bodies: dict[str, bytes], cls: type = MockS3Client) -> Any:
    return cls(
        keys=[H5_KEY, H5_KEY_B],
        crc_by_key={H5_KEY: "crc-a", H5_KEY_B: "crc-b"},
        object_bodies=bodies,
    )


def test_extract_h5_md5sum_is_lowercase_hex_per_object(tmp_path: Path) -> None:
    summary, rows = _h5_md5_run(
        tmp_path, _two_h5_client({H5_KEY: H5_BODY, H5_KEY_B: H5_BODY_B})
    )
    assert summary.md5_ok == 2 and not summary.failures
    assert rows[SAMPLE]["md5sum"] == _md5(H5_BODY)
    assert rows["sampleB"]["md5sum"] == _md5(H5_BODY_B)
    for row in rows.values():
        assert len(row["md5sum"]) == 32
        assert row["md5sum"] == row["md5sum"].lower()


def test_extract_h5_failed_md5_leaves_cell_empty_and_is_listed(
    tmp_path: Path,
) -> None:
    summary, rows = _h5_md5_run(tmp_path, _two_h5_client({H5_KEY_B: H5_BODY_B}))

    assert rows[SAMPLE]["md5sum"] == ""
    assert "NoSuchKey" in rows[SAMPLE]["md5_error"]
    assert rows[SAMPLE]["crc64nvme_base64"] == "crc-a"
    assert rows["sampleB"]["md5sum"] == _md5(H5_BODY_B)
    assert summary.crc_ok == 2
    assert summary.md5_ok == 1
    assert [(key, crc) for key, crc, _detail in summary.failures] == [(H5_KEY, "")]
    assert summary.failures[0][2].startswith("md5: ")


def test_extract_h5_failure_detail_keeps_the_h5_error(tmp_path: Path) -> None:
    client = MockS3Client(keys=[H5_KEY], crc_by_key={H5_KEY: "crc-a"})
    with patch(
        "file_extract.h5.introspect_h5", side_effect=RuntimeError("no matrix group")
    ):
        summary = extract_h5(
            client,
            BUCKET,
            PREFIX,
            str(tmp_path / "both.tsv"),
            lab="example-lab",
            do_introspect=True,
            retries=1,
            show_progress=False,
        )
    (detail,) = [detail for _key, _crc, detail in summary.failures]
    assert detail.startswith("no matrix group; md5: ")


def test_extract_h5_no_md5_drops_the_columns_and_reads_nothing(
    tmp_path: Path,
) -> None:
    client = _two_h5_client({H5_KEY: H5_BODY, H5_KEY_B: H5_BODY_B})
    summary, rows = _h5_md5_run(tmp_path, client, do_md5=False)

    assert "md5sum" not in rows[SAMPLE]
    assert "md5_error" not in rows[SAMPLE]
    assert summary.md5_ok == 0 and not summary.failures
    assert client.get_object_calls == []


def test_extract_h5_retries_a_transient_md5_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from botocore.exceptions import ClientError

    class ThrottledOnceClient(MockS3Client):
        """Throttles the first GET of every key, then serves it."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.throttled: set[str] = set()

        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            if Key not in self.throttled:
                self.throttled.add(Key)
                raise ClientError(
                    {"Error": {"Code": "SlowDown", "Message": "slow"}}, "GetObject"
                )
            return super().get_object(Bucket=Bucket, Key=Key, **kwargs)

    monkeypatch.setattr("file_extract.retry.time.sleep", lambda _seconds: None)
    bodies = {H5_KEY: H5_BODY, H5_KEY_B: H5_BODY_B}

    summary, rows = _h5_md5_run(
        tmp_path, _two_h5_client(bodies, ThrottledOnceClient), retries=2
    )
    assert summary.md5_ok == 2 and not summary.failures
    assert rows[SAMPLE]["md5sum"] == _md5(H5_BODY)

    summary, rows = _h5_md5_run(
        tmp_path, _two_h5_client(bodies, ThrottledOnceClient), retries=1
    )
    assert summary.md5_ok == 0 and len(summary.failures) == 2
    assert all("SlowDown" in detail for _k, _c, detail in summary.failures)
    assert rows[SAMPLE]["md5sum"] == ""


def test_extract_h5_md5_error_without_a_message_is_still_listed(
    tmp_path: Path,
) -> None:
    class SilentResetClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            raise ConnectionResetError()

    summary, rows = _h5_md5_run(
        tmp_path, _two_h5_client({}, SilentResetClient), retries=1
    )
    assert rows[SAMPLE]["md5sum"] == ""
    assert rows[SAMPLE]["md5_error"] == "failed without an error message"
    assert summary.md5_ok == 0 and len(summary.failures) == 2


def test_extract_h5_interrupt_cancels_queued_files(tmp_path: Path) -> None:
    import time

    class InterruptingClient(MockS3Client):
        """The first GET is interrupted; later ones are slow enough to observe."""

        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            first = not self.get_object_calls
            result = super().get_object(Bucket=Bucket, Key=Key, **kwargs)
            if first:
                raise KeyboardInterrupt
            time.sleep(0.2)
            return result

    keys = [f"{PREFIX}s{i}/sample_filtered_feature_bc_matrix.h5" for i in range(6)]
    client = InterruptingClient(
        keys=keys,
        crc_by_key={key: "crc" for key in keys},
        object_bodies={key: H5_BODY for key in keys},
    )
    with pytest.raises(KeyboardInterrupt):
        extract_h5(
            client,
            BUCKET,
            PREFIX,
            str(tmp_path / "interrupted.tsv"),
            lab="example-lab",
            do_introspect=False,
            workers=1,
            show_progress=False,
        )
    # The interrupted stream, plus at most the one the worker had already
    # picked up. Without the cancel, all six run to completion first.
    assert len(client.get_object_calls) <= 2
