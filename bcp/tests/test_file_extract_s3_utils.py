"""Tests for file_extract.s3_utils."""

from __future__ import annotations

import hashlib
import io
from typing import Any

import pytest
from botocore.exceptions import IncompleteReadError

from file_extract.s3_utils import (
    fetch_crc64nvme,
    list_objects_with_size,
    parse_s3_uri,
    s3_uri_for,
    stream_md5,
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


MD5_BUCKET = "example-bucket"
MD5_KEY_A = "proj/libA/processed/outs/per_sample_outs/s1/matrix.h5"
MD5_KEY_B = "proj/libB/processed/outs/per_sample_outs/s2/matrix.h5"
MD5_BODY_A = b"a" * 1000 + b"tail"
MD5_BODY_B = bytes(range(256)) * 3


def _md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


# --- stream_md5 --------------------------------------------------------------


def test_stream_md5_matches_hashlib_and_counts_bytes() -> None:
    client = MockS3Client(object_bodies={MD5_KEY_A: MD5_BODY_A})
    assert stream_md5(client, MD5_BUCKET, MD5_KEY_A) == (
        _md5(MD5_BODY_A),
        len(MD5_BODY_A),
    )
    assert client.get_object_calls == [
        {"Bucket": MD5_BUCKET, "Key": MD5_KEY_A, "ChecksumMode": "ENABLED"}
    ]


def test_stream_md5_chunked_read_is_the_same_digest() -> None:
    client = MockS3Client(object_bodies={MD5_KEY_B: MD5_BODY_B})
    assert stream_md5(client, MD5_BUCKET, MD5_KEY_B, chunk_bytes=7) == (
        _md5(MD5_BODY_B),
        len(MD5_BODY_B),
    )


def test_stream_md5_empty_object() -> None:
    client = MockS3Client(object_bodies={MD5_KEY_A: b""})
    assert stream_md5(client, MD5_BUCKET, MD5_KEY_A) == (_md5(b""), 0)


def test_stream_md5_short_read_is_an_error() -> None:
    class TruncatingClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {
                "Body": io.BytesIO(MD5_BODY_A[:-3]),
                "ContentLength": len(MD5_BODY_A),
            }

    client = TruncatingClient(object_bodies={MD5_KEY_A: MD5_BODY_A})
    # botocore's own type, so retry_with_backoff treats it as transient too.
    with pytest.raises(
        IncompleteReadError, match=r"1001 read, but total bytes expected is 1004"
    ):
        stream_md5(client, MD5_BUCKET, MD5_KEY_A)


def test_stream_md5_over_long_read_is_an_error() -> None:
    class PaddingClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {
                "Body": io.BytesIO(MD5_BODY_A + b"xyz"),
                "ContentLength": len(MD5_BODY_A),
            }

    with pytest.raises(
        IncompleteReadError, match=r"1007 read, but total bytes expected is 1004"
    ):
        stream_md5(PaddingClient(), MD5_BUCKET, MD5_KEY_A)


def test_stream_md5_without_content_length_trusts_the_stream() -> None:
    class NoLengthClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {"Body": io.BytesIO(MD5_BODY_A)}

    assert stream_md5(NoLengthClient(), MD5_BUCKET, MD5_KEY_A) == (
        _md5(MD5_BODY_A),
        len(MD5_BODY_A),
    )


def test_stream_md5_with_a_real_boto3_client(monkeypatch: pytest.MonkeyPatch) -> None:
    """The mock above cannot vouch for ChecksumMode or StreamingBody; moto can."""
    pytest.importorskip("moto")
    import boto3
    from moto import mock_aws

    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-east-1")
    with mock_aws():
        client = boto3.client("s3", region_name="us-east-1")
        client.create_bucket(Bucket=MD5_BUCKET)
        # Stored with a checksum botocore computes natively, so the response
        # carries x-amz-checksum-crc32 and the validating body is exercised.
        client.put_object(
            Bucket=MD5_BUCKET, Key=MD5_KEY_B, Body=MD5_BODY_B, ChecksumAlgorithm="CRC32"
        )
        assert stream_md5(client, MD5_BUCKET, MD5_KEY_B, chunk_bytes=100) == (
            _md5(MD5_BODY_B),
            len(MD5_BODY_B),
        )
