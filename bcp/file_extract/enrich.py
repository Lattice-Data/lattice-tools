"""Shared per-file enrichment for deliverables under an S3 order prefix.

FASTQ and CRAM deliverables are enriched identically: fetch the object's
CRC64NVME, then read ``read_count`` from the companion ``<key>-metadata.json``
the CRO writes beside it. Both retry transient S3 errors the same way and both
fan the work out across a process pool. That orchestration lives here so a
change to retry behaviour, error shaping or the pool pattern lands once instead
of being copied per file format.

File matching, row shaping and per-format tallies deliberately stay in fastq.py
and cram.py: those genuinely differ, and pulling them in here would trade
duplication for a callback-heavy indirection that reads worse than either.
"""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from tqdm import tqdm

from .retry import retry_with_backoff
from .s3_utils import fetch_crc64nvme, get_object_bytes

DEFAULT_WORKERS = 64


def fetch_read_count(s3_client: Any, bucket: str, key: str) -> int:
    """Read read_count from the companion -metadata.json beside an object."""
    metadata_key = key + "-metadata.json"
    data = json.loads(get_object_bytes(s3_client, bucket, metadata_key))
    if "read_count" not in data:
        raise RuntimeError("'read_count' not present in metadata JSON")
    return data["read_count"]


def fetch_one(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Fetch CRC and read_count for a single deliverable key.

    Never raises: a failed fetch leaves its value None and records the reason,
    so one unreadable object cannot abandon a whole order.
    """
    result: dict[str, object] = {
        "crc": None,
        "crc_error": "",
        "read_count": None,
        "metadata_error": "",
    }

    crc, crc_err = retry_with_backoff(
        fetch_crc64nvme, s3_client, bucket, key, retries=retries
    )
    result["crc"] = crc
    result["crc_error"] = crc_err or ""

    rc, rc_err = retry_with_backoff(
        fetch_read_count, s3_client, bucket, key, retries=retries
    )
    if rc_err:
        result["metadata_error"] = rc_err
    else:
        result["read_count"] = rc

    return result


def process_one(
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Enrich one key in a worker process.

    Module-level, and builds its own client, because a boto3 client cannot be
    pickled across a ProcessPoolExecutor boundary.
    """
    import boto3

    s3_client = boto3.client("s3")
    return fetch_one(s3_client, bucket, key, retries=retries)


def fetch_results(
    s3_client: Any,
    bucket: str,
    keys: list[str],
    *,
    retries: int = 5,
    workers: int | None = None,
    show_progress: bool = True,
    inline: bool = False,
) -> dict[str, dict[str, object]]:
    """Enrich every key, returned keyed by S3 key rather than completion order.

    ``inline`` runs in-process, which is what lets tests drive a mock client the
    pool could not pickle.
    """
    if inline:
        return {key: fetch_one(s3_client, bucket, key, retries=retries) for key in keys}

    results: dict[str, dict[str, object]] = {}
    max_workers = min(workers or DEFAULT_WORKERS, len(keys))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_one, bucket, key, retries=retries): key
            for key in keys
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(keys), desc="Fetching")
        for fut in iterator:
            results[futures[fut]] = fut.result()
    return results
