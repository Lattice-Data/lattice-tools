from __future__ import annotations

import hashlib
from typing import Any, Callable
from urllib.parse import urlparse

from botocore.exceptions import IncompleteReadError

from .constants import GEO_MD5_CHUNK_BYTES
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


def parse_s3_object_uri(uri: str) -> tuple[str, str]:
    """Parse ``s3://bucket/key`` into ``(bucket, key)``.

    Unlike :func:`parse_s3_uri` this names one object, so a bare bucket or a
    trailing slash is an error. Split by hand rather than with ``urlparse``: an
    S3 key may legally contain ``?`` or ``#``, which a URL parser would strip.
    """
    text = uri.strip()
    scheme = "s3://"
    if not text.startswith(scheme):
        raise ValueError(f"Invalid S3 URI: {uri}")
    bucket, slash, key = text[len(scheme) :].partition("/")
    if not bucket:
        raise ValueError(f"Invalid S3 URI: {uri}")
    if not slash or not key:
        raise ValueError(f"S3 URI does not name an object: {uri}")
    if key.endswith("/"):
        raise ValueError(f"S3 URI names a prefix, not an object: {uri}")
    if key.startswith("/"):
        raise ValueError(f"S3 key starts with '/': {uri}")
    if "\n" in key or "\r" in key:
        raise ValueError(f"S3 URI holds a line break: {uri!r}")
    return bucket, key


def stream_md5(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    chunk_bytes: int = GEO_MD5_CHUNK_BYTES,
) -> tuple[str, int]:
    """MD5 hex digest of a whole object, and the byte count it covers.

    The body is read in chunks so a multi-GB matrix never sits in memory.
    botocore's ``StreamingBody`` checks the byte count against
    ``ContentLength`` on the final read, and ``ChecksumMode=ENABLED`` asks S3
    for its stored checksum so botocore verifies the stream against that too,
    where it has the algorithm (CRC64NVME needs ``awscrt``; see
    ``requirements.txt``). The count is checked here as well, so a body that
    skips that step still cannot hand GEO a confident md5 of a truncated
    file. A short read and a checksum mismatch both surface as botocore
    exceptions that ``retry.is_transient`` retries.
    """
    resp = s3_client.get_object(Bucket=bucket, Key=key, ChecksumMode="ENABLED")
    expected = resp.get("ContentLength")
    body = resp["Body"]
    digest = hashlib.md5(usedforsecurity=False)
    total = 0
    while True:
        chunk = body.read(chunk_bytes)
        if not chunk:
            break
        digest.update(chunk)
        total += len(chunk)
    if expected is not None and total != expected:
        # The same type botocore raises, so the retry decision is the same.
        raise IncompleteReadError(actual_bytes=total, expected_bytes=expected)
    return digest.hexdigest(), total
