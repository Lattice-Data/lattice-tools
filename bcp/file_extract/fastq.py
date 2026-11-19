from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Any, Sequence

from .constants import FASTQ_COLUMNS, SHEET_HELPER_COLUMNS, SOURCE_ORDER_COLUMN
from .models import ListedObject, RunSummary
from .sequence_sets import OutputTarget, plan_and_fetch
from .sheets import (
    SequenceFileRecord,
    SheetOptions,
    build_fastq_record,
    enrich_record,
    group_records,
    sample_dir_for,
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


def parse_lane(fname: str) -> str:
    """Parse the lane from a FASTQ basename, empty when absent.

    Lane only: a file's read designator comes from parse_fastq_slot via
    SequenceFileRecord.slot, so it is not re-derived here.
    """
    lane_match = _LANE_RE.search(fname)
    return lane_match.group(1) if lane_match else ""


def default_fastq_output_name(prefix: str) -> str:
    order_name = prefix.rstrip("/").rsplit("/", 1)[-1] if prefix else "output"
    return f"{order_name}_fastq_info.tsv"


@dataclass(frozen=True)
class FastqListing:
    """One traversal of an order prefix: deliverables, and FASTQs left out.

    ``outside_raw`` holds the keys of FASTQs that would be deliverables but are
    not under a /raw/ folder. They are reported so they are not dropped
    silently; it stays empty under --no-require-raw, which includes them.
    """

    targets: tuple[ListedObject, ...] = ()
    outside_raw: tuple[str, ...] = ()


def scan_fastq_listing(
    s3_client: Any,
    bucket: str,
    prefix: str,
    *,
    require_raw: bool = True,
) -> FastqListing:
    """List deliverable FASTQs and the ones skipped for being outside /raw/."""
    targets: list[ListedObject] = []
    outside_raw: list[str] = []
    paginator = s3_client.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if is_target_file(key, require_raw=require_raw):
                targets.append(ListedObject(key=key, size_bytes=obj["Size"]))
            elif require_raw and is_target_file(key, require_raw=False):
                outside_raw.append(key)
    return FastqListing(targets=tuple(targets), outside_raw=tuple(outside_raw))


def fastq_columns(*, source_order: bool = False) -> list[str]:
    columns = list(FASTQ_COLUMNS) + list(SHEET_HELPER_COLUMNS)
    if source_order:
        columns.append(SOURCE_ORDER_COLUMN)
    return columns


def _diagnostic_row(
    record: SequenceFileRecord,
    result: dict[str, object],
    *,
    namespace: str | None,
    source_order: bool = False,
) -> list[object]:
    lane = parse_lane(record.filename)
    set_alias = ""
    if record.set_stem:
        set_alias = record.set_alias(namespace) if namespace else record.set_stem
    row: list[object] = [
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
    if source_order:
        row.append(record.cro_order)
    return row


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
    listing: FastqListing | None = None,
) -> RunSummary:
    """List FASTQ.gz files, enrich with CRC and read_count, write TSVs.

    Rows are buffered and emitted in S3-key order rather than in completion
    order: the submission sheets need R1 next to its R2, and a diff between two
    runs of the same order should be empty.

    Pass ``listing`` to reuse a scan the caller already made -- the CLI does, so
    the prefix is walked once for both the outside-raw/ warning and the work.
    """
    if listing is None:
        listing = scan_fastq_listing(s3_client, bucket, prefix, require_raw=require_raw)
    if not listing.targets:
        return RunSummary(total=0)
    target = OutputTarget(
        output_path=output_path,
        sheets=sheets,
        files=[(obj, "") for obj in listing.targets],
    )
    return extract_fastq_sets(
        s3_client,
        bucket,
        [target],
        workers=workers,
        retries=retries,
        show_progress=show_progress,
        inline=inline,
    )[0]


def extract_fastq_sets(
    s3_client: Any,
    bucket: str,
    targets: Sequence[OutputTarget],
    *,
    source_order_column: bool = False,
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
    inline: bool = False,
) -> list[RunSummary]:
    """Enrich and write every output set of one run, one summary per set.

    All sets are validated together and fetched in one pass; see
    sequence_sets.plan_and_fetch.
    """
    plans, results = plan_and_fetch(
        s3_client,
        bucket,
        targets,
        build_fastq_record,
        retries=retries,
        workers=workers,
        show_progress=show_progress,
        inline=inline,
    )
    return [
        _write_fastq_set(
            plan,
            results,
            target,
            source_order_column=source_order_column,
        )
        for plan, target in zip(plans, targets)
    ]


def _write_fastq_set(
    plan: Sequence[SequenceFileRecord],
    results: dict[str, dict[str, object]],
    target: OutputTarget,
    *,
    source_order_column: bool,
) -> RunSummary:
    sheets = target.sheets
    summary = RunSummary(total=len(plan))
    namespace = sheets.lab.namespace if sheets is not None else None
    writer = TsvWriter(
        target.output_path, fastq_columns(source_order=source_order_column)
    )
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

        writer.append_row(
            _diagnostic_row(
                record,
                result,
                namespace=namespace,
                source_order=source_order_column,
            )
        )
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
