from __future__ import annotations

import json
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from tqdm import tqdm

from .constants import FASTQ_COLUMNS, SHEET_HELPER_COLUMNS
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_bytes,
    list_objects_with_size,
    s3_uri_for,
)
from .sheets import (
    SequenceFileRecord,
    SheetOptions,
    build_fastq_record,
    enrich_record,
    group_records,
    parse_fastq_slot,
    sample_dir_for,
    validate_aliases,
    write_sheets,
)
from .tsv_writer import TsvWriter

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
    """Parse read (R1/R2/R3/I1/I2) and lane from a FASTQ basename.

    The read half defers to the set-stem parser so a file's read designator is
    the same value here, in the read tally, and in its SequenceFileSet slot.
    """
    _, slot, _ = parse_fastq_slot(fname)
    lane_match = _LANE_RE.search(fname)
    return (slot, (lane_match.group(1) if lane_match else ""))


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
        "read_count": None,
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
    return list(FASTQ_COLUMNS) + list(SHEET_HELPER_COLUMNS)


def _fetch_results(
    s3_client: Any,
    bucket: str,
    keys: list[str],
    *,
    retries: int,
    workers: int | None,
    show_progress: bool,
    inline: bool,
) -> dict[str, dict[str, object]]:
    """Enrich every key, returned keyed by S3 key rather than completion order."""
    if inline:
        return {
            key: fetch_one_fastq(s3_client, bucket, key, retries=retries)
            for key in keys
        }

    results: dict[str, dict[str, object]] = {}
    max_workers = min(workers or 64, len(keys))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_one_fastq, bucket, key, retries=retries): key
            for key in keys
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(keys), desc="Fetching")
        for fut in iterator:
            results[futures[fut]] = fut.result()
    return results


def _diagnostic_row(
    record: SequenceFileRecord,
    result: dict[str, object],
    *,
    namespace: str | None,
) -> list[object]:
    _, lane = parse_read_lane(record.filename)
    set_alias = ""
    if record.set_stem:
        set_alias = record.set_alias(namespace) if namespace else record.set_stem
    return [
        record.filename,
        record.s3_uri,
        record.slot,
        lane,
        record.file_size,
        result["crc"] if result["crc"] is not None else "",
        result["read_count"] if result["read_count"] is not None else "",
        result["crc_error"],
        result["metadata_error"],
        sample_dir_for(record.key),
        set_alias,
    ]


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
    sheets: SheetOptions | None = None,
) -> RunSummary:
    """List FASTQ.gz files, enrich with CRC and read_count, write TSVs.

    Rows are buffered and emitted in S3-key order rather than in completion
    order: the submission sheets need R1 next to its R2, and a diff between two
    runs of the same order should be empty.
    """
    targets = list_objects_with_size(
        s3_client,
        bucket,
        prefix,
        predicate=lambda k: is_target_file(k, require_raw=require_raw),
    )
    summary = RunSummary(total=len(targets))
    if not targets:
        return summary

    plan = sorted(
        (
            build_fastq_record(
                key=obj.key,
                s3_uri=s3_uri_for(bucket, obj.key),
                file_size=obj.size_bytes,
            )
            for obj in targets
        ),
        key=lambda record: record.sort_key,
    )
    namespace = sheets.lab.namespace if sheets is not None else None
    if namespace is not None:
        # Before spending a request per file: a collision here is unsubmittable.
        validate_aliases(plan, namespace=namespace)

    results = _fetch_results(
        s3_client,
        bucket,
        [record.key for record in plan],
        retries=retries,
        workers=workers,
        show_progress=show_progress,
        inline=inline,
    )

    writer = TsvWriter(output_path, fastq_columns())
    read_tally: Counter[str] = Counter()
    records: list[SequenceFileRecord] = []
    for record in plan:
        result = results[record.key]
        crc_err = str(result["crc_error"])
        meta_err = str(result["metadata_error"])
        read_tally[record.slot or "(none)"] += 1
        if not crc_err:
            summary.crc_ok += 1
        if not crc_err and not meta_err:
            summary.enrichment_ok += 1
        if crc_err or meta_err:
            summary.failures.append((record.key, crc_err, meta_err))

        writer.append_row(_diagnostic_row(record, result, namespace=namespace))
        records.append(
            enrich_record(
                record,
                crc=result["crc"] if isinstance(result["crc"], str) else None,
                read_count=(
                    result["read_count"]
                    if isinstance(result["read_count"], int)
                    else None
                ),
            )
        )

    summary.read_tally = dict(read_tally)

    if sheets is not None:
        groups, warnings = group_records(records)
        write_sheets(records, groups, options=sheets)
        summary.set_count = len(groups)
        summary.warnings.extend(warnings)

    return summary


def r1_r2_mismatch_warning(read_tally: dict[str, int]) -> str | None:
    """Return warning message when R1 and R2 counts differ."""
    r1 = read_tally.get("R1", 0)
    r2 = read_tally.get("R2", 0)
    if r1 != r2:
        return f"R1 ({r1}) and R2 ({r2}) counts differ -- possible unpaired files"
    return None
