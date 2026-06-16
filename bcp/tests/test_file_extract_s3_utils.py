"""Tests for file_extract.s3_utils."""

from __future__ import annotations

import pytest

from file_extract.s3_utils import (
    fetch_crc64nvme,
    list_objects_with_size,
    parse_s3_uri,
    s3_uri_for,
)
from tests.file_extract_helpers import MockS3Client


def test_parse_s3_uri_with_prefix() -> None:
    loc = parse_s3_uri("s3://example-bucket/path/to/order")
    assert loc.bucket == "example-bucket"
    assert loc.prefix == "path/to/order/"


def test_parse_s3_uri_bucket_only() -> None:
    loc = parse_s3_uri("s3://example-bucket")
    assert loc.bucket == "example-bucket"
    assert loc.prefix == ""


def test_parse_s3_uri_invalid() -> None:
    with pytest.raises(ValueError, match="Invalid S3 URI"):
        parse_s3_uri("https://example.com/foo")


def test_s3_uri_for() -> None:
    assert s3_uri_for("b", "k/e.y") == "s3://b/k/e.y"


def test_list_objects_with_size() -> None:
    keys = [
        "order/raw/a.fastq.gz",
        "order/raw/b.fastq.gz",
        "order/other.txt",
    ]
    client = MockS3Client(keys=keys, sizes={k: 100 for k in keys})
    objs = list_objects_with_size(
        client,
        "b",
        "order/",
        predicate=lambda k: k.endswith(".fastq.gz"),
    )
    assert len(objs) == 2
    assert objs[0].size_bytes == 100


def test_fetch_crc64nvme_success() -> None:
    client = MockS3Client(crc_by_key={"k": "abc123"})
    assert fetch_crc64nvme(client, "b", "k") == "abc123"


def test_fetch_crc64nvme_missing() -> None:
    client = MockS3Client()
    with pytest.raises(RuntimeError, match="ChecksumCRC64NVME"):
        fetch_crc64nvme(client, "b", "missing")
