from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from .constants import CRAM_COLUMNS, SHEET_HELPER_COLUMNS, SOURCE_ORDER_COLUMN
from .models import ListedObject, RunSummary
from .sequence_sets import OutputTarget, plan_and_fetch
from .sheets import (
    SequenceFileRecord,
    SheetOptions,
    build_cram_record,
    enrich_record,
    group_records,
    sample_dir_for,
    write_sheets,
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


def default_cram_output_name(prefix: str) -> str:
    order_name = prefix.rstrip("/").rsplit("/", 1)[-1] if prefix else "output"
    return f"{order_name}_cram_info.tsv"


def cram_columns(*, source_order: bool = False) -> list[str]:
    columns = list(CRAM_COLUMNS) + list(SHEET_HELPER_COLUMNS)
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
    set_alias = record.set_alias(namespace) if namespace else record.set_stem
    row: list[object] = [
        record.filename,
        record.s3_uri,
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


@dataclass(frozen=True)
class CramListingWarnings:
    """Counts from a prefix scan used for CLI guardrails."""

    unmatched_cram_count: int = 0
    ucram_count: int = 0


@dataclass(frozen=True)
class CramListing:
    """One traversal of an order prefix, classified into targets and warnings.

    ``outside_raw`` holds the keys of CRAMs that would be deliverables but are
    not under a /raw/ folder. They are reported so they are not dropped
    silently; it stays empty under --no-require-raw, which includes them.
    """

    targets: tuple[ListedObject, ...] = ()
    warnings: CramListingWarnings = field(default_factory=CramListingWarnings)
    outside_raw: tuple[str, ...] = ()


def scan_cram_listing(
    s3_client: Any,
    bucket: str,
    prefix: str,
    *,
    require_raw: bool = True,
) -> CramListing:
    """Classify every key under the prefix in a single traversal.

    Deliverables and the guardrail counts come from the same listing so they
    cannot disagree, and the prefix is walked once rather than once per concern.

    The counts are deliberately not scoped to /raw/: a .ucram or an unmatched
    CRAM anywhere under the order is worth reporting, because it means the
    delivery holds more than this run processed.
    """
    targets: list[ListedObject] = []
    outside_raw: list[str] = []
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
            elif is_target_file(key, require_raw=require_raw):
                targets.append(ListedObject(key=key, size_bytes=obj["Size"]))
            elif require_raw and is_target_file(key, require_raw=False):
                outside_raw.append(key)
    return CramListing(
        targets=tuple(targets),
        warnings=CramListingWarnings(
            unmatched_cram_count=unmatched,
            ucram_count=ucram,
        ),
        outside_raw=tuple(outside_raw),
    )


def only_unmatched_cram_warning(
    unmatched_count: int,
    *,
    require_raw: bool = True,
) -> str | None:
    """Return warning when only unmatched CRAMs were found under the prefix.

    The /raw/ clause is dropped under --no-require-raw: claiming nothing was
    found there would name a location the run never restricted itself to.
    """
    if unmatched_count <= 0:
        return None
    location = " under /raw/" if require_raw else ""
    return (
        f"Found {unmatched_count} unmatched .cram file(s) but no deliverable "
        f"CRAMs{location} (unmatched files are excluded)"
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
    sheets: SheetOptions | None = None,
    listing: CramListing | None = None,
) -> RunSummary:
    """Enrich CRAM deliverables with CRC and read_count, write TSVs.

    Rows are buffered and emitted in S3-key order so repeated runs of the same
    order produce identical files.

    Pass ``listing`` to reuse a scan the caller already made -- the CLI does, so
    the prefix is walked once for both the guardrail warnings and the work.
    """
    if listing is None:
        listing = scan_cram_listing(s3_client, bucket, prefix, require_raw=require_raw)
    if not listing.targets:
        return RunSummary(total=0)
    target = OutputTarget(
        output_path=output_path,
        sheets=sheets,
        files=[(obj, "") for obj in listing.targets],
    )
    return extract_cram_sets(
        s3_client,
        bucket,
        [target],
        workers=workers,
        retries=retries,
        show_progress=show_progress,
        inline=inline,
    )[0]


def extract_cram_sets(
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
        build_cram_record,
        retries=retries,
        workers=workers,
        show_progress=show_progress,
        inline=inline,
    )
    return [
        _write_cram_set(
            plan,
            results,
            target,
            source_order_column=source_order_column,
        )
        for plan, target in zip(plans, targets)
    ]


def _write_cram_set(
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
        target.output_path, cram_columns(source_order=source_order_column)
    )
    records: list[SequenceFileRecord] = []
    for record in plan:
        result = results[record.key]
        crc_err = str(result["crc_error"])
        meta_err = str(result["metadata_error"])
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

    if sheets is not None:
        groups, warnings = group_records(records)
        write_sheets(records, groups, options=sheets)
        summary.set_count = len(groups)
        summary.warnings.extend(warnings)

    return summary
