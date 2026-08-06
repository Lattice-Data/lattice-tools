from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any

from tqdm import tqdm

from .constants import CRAM_COLUMNS
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_bytes,
    list_objects_with_size,
    s3_uri_for,
)
from .tsv_writer import TsvWriter


def is_unmatched_cram(key: str) -> bool:
    """Return True if key is an unmatched CRAM deliverable."""
    if not key.endswith(".cram"):
        return False
    fname = key.rsplit("/", 1)[-1]
    return "_unmatched.cram" in fname or "-unmatched.cram" in fname


def is_target_file(key: str, *, require_raw: bool = True) -> bool:
    """Return True if key is a deliverable CRAM under the order prefix."""
    if require_raw and "/raw/" not in key:
        return False
    if not key.endswith(".cram"):
        return False
    if is_unmatched_cram(key):
        return False
    return True


def _fetch_read_count(s3_client: Any, bucket: str, key: str) -> int:
    metadata_key = key + "-metadata.json"
    data = json.loads(get_object_bytes(s3_client, bucket, metadata_key))
    if "read_count" not in data:
        raise RuntimeError("'read_count' not present in metadata JSON")
    return data["read_count"]


def fetch_one_cram(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Fetch CRC and read_count for a single CRAM key."""
    result: dict[str, object] = {
        "crc": None,
        "crc_error": "",
        "read_count": "",
        "metadata_error": "",
    }

    crc, crc_err = retry_with_backoff(
        fetch_crc64nvme, s3_client, bucket, key, retries=retries
    )
    result["crc"] = crc
    result["crc_error"] = crc_err or ""

    rc, rc_err = retry_with_backoff(
        _fetch_read_count, s3_client, bucket, key, retries=retries
    )
    if rc_err:
        result["metadata_error"] = rc_err
    else:
        result["read_count"] = rc

    return result


def process_one_cram(
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Process a single CRAM key (module-level for ProcessPoolExecutor)."""
    import boto3

    s3_client = boto3.client("s3")
    return fetch_one_cram(s3_client, bucket, key, retries=retries)


def default_cram_output_name(prefix: str) -> str:
    order_name = prefix.rstrip("/").rsplit("/", 1)[-1] if prefix else "output"
    return f"{order_name}_cram_info.tsv"


def cram_columns() -> list[str]:
    return list(CRAM_COLUMNS)


@dataclass(frozen=True)
class CramListingWarnings:
    """Counts from a prefix scan used for CLI guardrails."""

    unmatched_cram_count: int = 0
    ucram_count: int = 0


def scan_cram_listing_warnings(
    s3_client: Any,
    bucket: str,
    prefix: str,
) -> CramListingWarnings:
    """Scan prefix for unmatched CRAM and .ucram keys."""
    unmatched = 0
    ucram = 0
    paginator = s3_client.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith(".ucram"):
                ucram += 1
            elif is_unmatched_cram(key):
                unmatched += 1
    return CramListingWarnings(
        unmatched_cram_count=unmatched,
        ucram_count=ucram,
    )


def only_unmatched_cram_warning(unmatched_count: int) -> str | None:
    """Return warning when only unmatched CRAMs were found under the prefix."""
    if unmatched_count <= 0:
        return None
    return (
        f"Found {unmatched_count} unmatched .cram file(s) but no deliverable "
        "CRAMs under /raw/ (unmatched files are excluded)"
    )


def ucram_found_warning(ucram_count: int) -> str | None:
    """Return warning when .ucram files are present under the prefix."""
    if ucram_count <= 0:
        return None
    return f"Found {ucram_count} .ucram file(s) under the prefix (not processed)"


def extract_cram(
    s3_client: Any,
    bucket: str,
    prefix: str,
    output_path: str,
    *,
    require_raw: bool = True,
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
    inline: bool = False,
) -> RunSummary:
    """List CRAM files, enrich with CRC and read_count, write TSV."""
    targets = list_objects_with_size(
        s3_client,
        bucket,
        prefix,
        predicate=lambda k: is_target_file(k, require_raw=require_raw),
    )
    summary = RunSummary(total=len(targets))
    if not targets:
        return summary

    writer = TsvWriter(output_path, cram_columns())
    size_by_key = {obj.key: obj.size_bytes for obj in targets}
    max_workers = min(workers or 64, len(targets))

    def _handle_result(key: str, r: dict[str, object]) -> None:
        nonlocal summary
        fname = key.rsplit("/", 1)[-1]

        crc_err = str(r["crc_error"])
        meta_err = str(r["metadata_error"])
        if not crc_err:
            summary.crc_ok += 1
        if not crc_err and not meta_err:
            summary.enrichment_ok += 1
        if crc_err or meta_err:
            summary.failures.append((key, crc_err, meta_err))

        writer.append_row(
            [
                fname,
                s3_uri_for(bucket, key),
                size_by_key[key],
                r["crc"] if r["crc"] is not None else "",
                r["read_count"],
                crc_err,
                meta_err,
            ]
        )

    if inline:
        for obj in targets:
            r = fetch_one_cram(s3_client, bucket, obj.key, retries=retries)
            _handle_result(obj.key, r)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_one_cram, bucket, obj.key, retries=retries
                ): obj.key
                for obj in targets
            }
            iterator = as_completed(futures)
            if show_progress:
                iterator = tqdm(iterator, total=len(targets), desc="Fetching")

            for fut in iterator:
                key = futures[fut]
                _handle_result(key, fut.result())

    return summary
