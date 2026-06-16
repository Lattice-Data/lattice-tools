from __future__ import annotations

from typing import Any, Callable
from urllib.parse import urlparse

from .models import ListedObject, S3Location
from .retry import retry_with_backoff


def parse_s3_uri(uri: str) -> S3Location:
    """Parse s3://bucket/prefix into bucket and normalized prefix (trailing /)."""
    parsed = urlparse(uri)
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError(f"Invalid S3 URI: {uri}")
    bucket = parsed.netloc
    prefix = parsed.path.strip("/")
    if prefix:
        prefix += "/"
    return S3Location(bucket=bucket, prefix=prefix)


def list_objects_with_size(
    s3_client: Any,
    bucket: str,
    prefix: str,
    *,
    predicate: Callable[[str], bool] | None = None,
) -> list[ListedObject]:
    """List objects under prefix, optionally filtered by predicate(key) -> bool."""
    paginator = s3_client.get_paginator("list_objects_v2")
    results: list[ListedObject] = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if predicate is None or predicate(key):
                results.append(ListedObject(key=key, size_bytes=obj["Size"]))
    return results


def fetch_crc64nvme(s3_client: Any, bucket: str, key: str) -> str:
    """Fetch ChecksumCRC64NVME from object attributes."""
    resp = s3_client.get_object_attributes(
        Bucket=bucket,
        Key=key,
        ObjectAttributes=["Checksum"],
    )
    crc = resp.get("Checksum", {}).get("ChecksumCRC64NVME")
    if crc is None:
        raise RuntimeError("No ChecksumCRC64NVME in object attributes")
    return crc


def fetch_crc64nvme_with_retry(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> tuple[str | None, str | None]:
    """Fetch CRC with retry; returns (crc, error_str)."""
    return retry_with_backoff(fetch_crc64nvme, s3_client, bucket, key, retries=retries)


def get_object_bytes(s3_client: Any, bucket: str, key: str) -> bytes:
    """Download object body as bytes."""
    resp = s3_client.get_object(Bucket=bucket, Key=key)
    return resp["Body"].read()


def get_object_text(s3_client: Any, bucket: str, key: str) -> str:
    """Download object body as UTF-8 text."""
    return get_object_bytes(s3_client, bucket, key).decode("utf-8")


def s3_uri_for(bucket: str, key: str) -> str:
    return f"s3://{bucket}/{key}"
