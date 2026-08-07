from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .constants import CRAM_COLUMNS, SHEET_HELPER_COLUMNS
from .enrich import fetch_results
from .models import ListedObject, RunSummary
from .s3_utils import s3_uri_for
from .sheets import (
    SequenceFileRecord,
    SheetOptions,
    build_cram_record,
    enrich_record,
    group_records,
    sample_dir_for,
    validate_plan,
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


def cram_columns() -> list[str]:
    return list(CRAM_COLUMNS) + list(SHEET_HELPER_COLUMNS)


def _diagnostic_row(
    record: SequenceFileRecord,
    result: dict[str, object],
    *,
    namespace: str | None,
) -> list[object]:
    set_alias = record.set_alias(namespace) if namespace else record.set_stem
    return [
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


@dataclass(frozen=True)
class CramListingWarnings:
    """Counts from a prefix scan used for CLI guardrails."""

    unmatched_cram_count: int = 0
    ucram_count: int = 0


@dataclass(frozen=True)
class CramListing:
    """One traversal of an order prefix, classified into targets and warnings."""

    targets: tuple[ListedObject, ...] = ()
    warnings: CramListingWarnings = field(default_factory=CramListingWarnings)


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
    return CramListing(
        targets=tuple(targets),
        warnings=CramListingWarnings(
            unmatched_cram_count=unmatched,
            ucram_count=ucram,
        ),
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
    targets = listing.targets
    summary = RunSummary(total=len(targets))
    if not targets:
        return summary

    plan = sorted(
        (
            build_cram_record(
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
        # Before spending a request per file: these listings are unsubmittable.
        validate_plan(plan, namespace=namespace)

    results = fetch_results(
        s3_client,
        bucket,
        [record.key for record in plan],
        retries=retries,
        workers=workers,
        show_progress=show_progress,
        inline=inline,
    )

    writer = TsvWriter(output_path, cram_columns())
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

    if sheets is not None:
        groups, warnings = group_records(records)
        write_sheets(records, groups, options=sheets)
        summary.set_count = len(groups)
        summary.warnings.extend(warnings)

    return summary
