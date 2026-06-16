from __future__ import annotations

import json
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from tqdm import tqdm

from .constants import FASTQ_COLUMNS
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_bytes,
    list_objects_with_size,
    s3_uri_for,
)
from .tsv_writer import TsvWriter

_READ_RE = re.compile(r"_(R[12]|I[12])_")
_LANE_RE = re.compile(r"_L(\d{3})_")


def is_target_file(key: str, *, require_raw: bool = True) -> bool:
    """Return True if key is a deliverable FASTQ.gz under the order prefix."""
    if require_raw and "/raw/" not in key:
        return False
    if not key.endswith(".fastq.gz"):
        return False
    fname = key.rsplit("/", 1)[-1].lower()
    if "_sample" in fname or "unmatched" in fname:
        return False
    return True


def parse_read_lane(fname: str) -> tuple[str, str]:
    """Parse read (R1/R2/I1/I2) and lane from a FASTQ basename."""
    read_match = _READ_RE.search(fname)
    lane_match = _LANE_RE.search(fname)
    return (
        (read_match.group(1) if read_match else ""),
        (lane_match.group(1) if lane_match else ""),
    )


def _fetch_read_count(s3_client: Any, bucket: str, key: str) -> int:
    metadata_key = key + "-metadata.json"
    data = json.loads(get_object_bytes(s3_client, bucket, metadata_key))
    if "read_count" not in data:
        raise RuntimeError("'read_count' not present in metadata JSON")
    return data["read_count"]


def fetch_one_fastq(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Fetch CRC and read_count for a single FASTQ key."""
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


def process_one_fastq(
    bucket: str,
    key: str,
    *,
    retries: int = 5,
) -> dict[str, object]:
    """Process a single FASTQ key (module-level for ProcessPoolExecutor)."""
    import boto3

    s3_client = boto3.client("s3")
    return fetch_one_fastq(s3_client, bucket, key, retries=retries)


def default_fastq_output_name(prefix: str) -> str:
    order_name = prefix.rstrip("/").rsplit("/", 1)[-1] if prefix else "output"
    return f"{order_name}_fastq_info.tsv"


def fastq_columns() -> list[str]:
    return list(FASTQ_COLUMNS)


def extract_fastq(
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
    """List FASTQ.gz files, enrich with CRC and read_count, write TSV."""
    targets = list_objects_with_size(
        s3_client,
        bucket,
        prefix,
        predicate=lambda k: is_target_file(k, require_raw=require_raw),
    )
    summary = RunSummary(total=len(targets))
    if not targets:
        return summary

    writer = TsvWriter(output_path, fastq_columns())
    size_by_key = {obj.key: obj.size_bytes for obj in targets}
    max_workers = min(workers or 64, len(targets))
    read_tally: Counter[str] = Counter()

    def _handle_result(key: str, r: dict[str, object]) -> None:
        nonlocal summary
        fname = key.rsplit("/", 1)[-1]
        read, lane = parse_read_lane(fname)
        read_tally[read or "(none)"] += 1

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
                read,
                lane,
                size_by_key[key],
                r["crc"] if r["crc"] is not None else "",
                r["read_count"],
                crc_err,
                meta_err,
            ]
        )

    if inline:
        for obj in targets:
            r = fetch_one_fastq(s3_client, bucket, obj.key, retries=retries)
            _handle_result(obj.key, r)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_one_fastq, bucket, obj.key, retries=retries
                ): obj.key
                for obj in targets
            }
            iterator = as_completed(futures)
            if show_progress:
                iterator = tqdm(iterator, total=len(targets), desc="Fetching")

            for fut in iterator:
                key = futures[fut]
                _handle_result(key, fut.result())

    summary.read_tally = dict(read_tally)
    return summary


def r1_r2_mismatch_warning(read_tally: dict[str, int]) -> str | None:
    """Return warning message when R1 and R2 counts differ."""
    r1 = read_tally.get("R1", 0)
    r2 = read_tally.get("R2", 0)
    if r1 != r2:
        return f"R1 ({r1}) and R2 ({r2}) counts differ -- possible unpaired files"
    return None
