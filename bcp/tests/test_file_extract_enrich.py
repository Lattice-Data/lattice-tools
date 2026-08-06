"""Tests for file_extract.enrich, the enrichment core shared by fastq and cram."""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from file_extract.enrich import (
    fetch_one,
    fetch_read_count,
    fetch_results,
    process_one,
)
from tests.file_extract_helpers import MockS3Client

BUCKET = "example-bucket"
KEY = "test-order/NVUS0000000000-01/S1/raw/sample_L001_R1_001.fastq.gz"


def _client(
    *, read_count: int | None = 42, crc: str | None = "crc-abc"
) -> MockS3Client:
    bodies = {}
    if read_count is not None:
        bodies[f"{KEY}-metadata.json"] = json.dumps({"read_count": read_count})
    return MockS3Client(
        keys=[KEY],
        sizes={KEY: 100},
        crc_by_key={KEY: crc} if crc is not None else {},
        object_bodies=bodies,
    )


def test_fetch_read_count_reads_companion_json() -> None:
    assert fetch_read_count(_client(read_count=7), BUCKET, KEY) == 7


def test_fetch_read_count_rejects_json_without_read_count() -> None:
    """A companion file that parses but omits the field is a delivery problem."""
    client = MockS3Client(object_bodies={f"{KEY}-metadata.json": json.dumps({"a": 1})})
    with pytest.raises(RuntimeError, match="read_count"):
        fetch_read_count(client, BUCKET, KEY)


def test_fetch_one_returns_crc_and_read_count() -> None:
    result = fetch_one(_client(), BUCKET, KEY, retries=1)
    assert result == {
        "crc": "crc-abc",
        "crc_error": "",
        "read_count": 42,
        "metadata_error": "",
    }


def test_fetch_one_records_missing_metadata_without_raising() -> None:
    """One unreadable companion file must not abandon the rest of an order."""
    result = fetch_one(_client(read_count=None), BUCKET, KEY, retries=1)
    assert result["crc"] == "crc-abc"
    assert result["read_count"] is None
    assert result["metadata_error"]


def test_fetch_one_records_missing_checksum_without_raising() -> None:
    result = fetch_one(_client(crc=None), BUCKET, KEY, retries=1)
    assert result["crc"] is None
    assert result["crc_error"]
    assert result["read_count"] == 42


def test_fetch_results_inline_keys_by_s3_key() -> None:
    results = fetch_results(
        _client(), BUCKET, [KEY], retries=1, show_progress=False, inline=True
    )
    assert list(results) == [KEY]
    assert results[KEY]["read_count"] == 42


def test_process_one_builds_its_own_client() -> None:
    """A boto3 client cannot be pickled into a worker, so each task makes one."""
    client = _client()
    with patch("boto3.client", return_value=client) as mock_client:
        result = process_one(BUCKET, KEY, retries=1)
    mock_client.assert_called_once_with("s3")
    assert result["read_count"] == 42
    assert result["crc"] == "crc-abc"
