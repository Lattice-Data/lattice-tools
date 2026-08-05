"""
Reconciliation of a SeaHub trimmed upload against its untrimmed vendor source.

Consumes the two ``(wafer, UG)`` indexes built by :mod:`qa_seahub_source` and
classifies every well.  See that module for why ``(wafer, UG)`` is the identity
key and how each side is indexed.

Read counts
-----------
The per-well failure CSV declares a different ``total read count`` per format
block for the same well, and the failure rows are largely shared between blocks,
so there is no defensible single "failed reads" figure for a well.  The check is
therefore framed as "did trimming consume the whole delivered file": the vendor
CRAM's ``read_count`` must equal one of the input totals the trimmer declared.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from qa_seahub_source import (
    IdentityKey,
    SourceCoverage,
    SourceEntry,
    TrimmedEntry,
    UntrimmedSources,
    finding_row as _row,
)

__all__ = [
    "TrimmingReport",
    "reconcile_trimming",
]


@dataclass
class TrimmingReport:
    """Outcome of reconciling a trimmed upload against its vendor source."""

    rows: list[dict[str, Any]] = field(default_factory=list)
    source_total: int = 0
    trimmed_total: int = 0
    matched: int = 0
    source_coverage: list[SourceCoverage] = field(default_factory=list)
    # Wafers present in the trimmed upload that no listed source delivered --
    # the signature of a forgotten vendor order.
    unsourced_wafers: list[str] = field(default_factory=list)

    @property
    def counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for row in self.rows:
            counts[row["category"]] = counts.get(row["category"], 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))

    def verdict(self) -> str:
        counts = self.counts
        verdict = (
            f"{self.source_total} vendor CRAMs, {self.trimmed_total} trimmed wells, "
            f"{self.matched} matched, "
            f"{counts.get('not_trimmed', 0)} not trimmed, "
            f"{counts.get('orphan_trimmed', 0)} orphan trimmed, "
            f"{counts.get('read_count_mismatch', 0)} read-count mismatches"
        )
        # Appended only when non-zero, so the five clauses above stay stable.
        for category, label in (
            ("size_not_reduced", "not smaller than source"),
            ("duplicate_source_well", "duplicate source wells"),
            ("duplicate_trimmed_well", "duplicate trimmed wells"),
        ):
            if counts.get(category):
                verdict += f", {counts[category]} {label}"
        if self.unsourced_wafers:
            verdict += (
                f"; {len(self.unsourced_wafers)} wafer(s) have no listed source "
                f"({', '.join(self.unsourced_wafers)}) -- an untrimmed order is "
                "probably missing from untrimmed_s3_paths"
            )
        return verdict


def _reconcile_read_counts(
    source: SourceEntry,
    trimmed: TrimmedEntry,
    fail_counts: dict[str, dict[str, dict[str, int]]],
) -> dict[str, Any] | None:
    """Check the vendor read count against the trimmer's declared input totals."""
    if source.read_count is None:
        return _row(
            "metadata_unavailable",
            wafer=source.wafer,
            ug=source.ug,
            sublibrary=trimmed.sublibrary,
            source_key=source.cram_key,
            trimmed_stem=trimmed.stem,
            detail="vendor .cram-metadata.json missing or has no read_count",
        )
    per_format = fail_counts.get(trimmed.stem)
    if not per_format:
        return _row(
            "metadata_unavailable",
            wafer=source.wafer,
            ug=source.ug,
            sublibrary=trimmed.sublibrary,
            source_key=source.cram_key,
            trimmed_stem=trimmed.stem,
            detail="no per-well trimmer failure counts parsed for this stem",
        )
    totals = {int(stats["total"]) for stats in per_format.values()}
    if source.read_count in totals:
        return None
    return _row(
        "read_count_mismatch",
        wafer=source.wafer,
        ug=source.ug,
        sublibrary=trimmed.sublibrary,
        source_key=source.cram_key,
        trimmed_stem=trimmed.stem,
        detail=(
            f"vendor read_count {source.read_count} matches none of the trimmer "
            f"input totals {sorted(totals)}, so trimming did not consume the "
            "whole delivered file"
        ),
    )


def _reconcile_size(
    source: SourceEntry, trimmed: TrimmedEntry
) -> dict[str, Any] | None:
    """Trimming removes reads, so the output must be smaller than its input.

    Silent when either size is unknown (0): sizes are only collected in S3 mode.
    A trimmed CRAM at or above its source size is the signature of untrimmed data
    landing in the trimmed bucket -- the inverse of the size ratio that confirmed
    REF3's uploads really were post-trim (14.4 GB against a 35-45 GB source).
    """
    if not source.size_bytes or not trimmed.size_bytes:
        return None
    if trimmed.size_bytes < source.size_bytes:
        return None
    return _row(
        "size_not_reduced",
        wafer=source.wafer,
        ug=source.ug,
        sublibrary=trimmed.sublibrary,
        source_key=source.cram_key,
        trimmed_stem=trimmed.stem,
        detail=(
            f"trimmed CRAM is {trimmed.size_bytes} bytes against a vendor source "
            f"of {source.size_bytes}; trimming should reduce the file, so this "
            "may be untrimmed data in the trimmed bucket"
        ),
    )


def _reconcile_identity(
    source: SourceEntry, trimmed: TrimmedEntry
) -> list[dict[str, Any]]:
    """Compare the fields that are not part of the identity key."""
    rows: list[dict[str, Any]] = []
    if source.barcode != trimmed.barcode:
        rows.append(
            _row(
                "identity_mismatch",
                wafer=source.wafer,
                ug=source.ug,
                sublibrary=trimmed.sublibrary,
                source_key=source.cram_key,
                trimmed_stem=trimmed.stem,
                detail=(
                    f"barcode differs: vendor {source.barcode} vs trimmed "
                    f"{trimmed.barcode}"
                ),
            )
        )
    if (
        source.assay is not None
        and trimmed.assay is not None
        and source.assay != trimmed.assay
    ):
        rows.append(
            _row(
                "identity_mismatch",
                wafer=source.wafer,
                ug=source.ug,
                sublibrary=trimmed.sublibrary,
                source_key=source.cram_key,
                trimmed_stem=trimmed.stem,
                detail=(
                    f"sublibrary type differs: vendor {source.assay} vs "
                    f"trimmed {trimmed.assay}"
                ),
            )
        )
    return rows


def reconcile_trimming(
    untrimmed: UntrimmedSources | dict[IdentityKey, SourceEntry],
    trimmed_index: dict[IdentityKey, TrimmedEntry],
    fail_counts: dict[str, dict[str, dict[str, int]]] | None = None,
    trimmed_findings: list[dict[str, Any]] | None = None,
) -> TrimmingReport:
    """Compare one or more vendor deliveries against a trimmed upload, well by well.

    Accepts either an :class:`UntrimmedSources` container or a bare index. Only
    the container carries per-prefix coverage and the indexers' own findings, so a
    bare index cannot report a forgotten order as anything but orphans.
    """
    fail_counts = fail_counts or {}
    if isinstance(untrimmed, UntrimmedSources):
        untrimmed_index = untrimmed.index
        seeded = list(untrimmed.findings)
        coverage = list(untrimmed.coverage)
    else:
        untrimmed_index = untrimmed
        seeded = []
        coverage = []

    report = TrimmingReport(
        source_total=len(untrimmed_index),
        trimmed_total=len(trimmed_index),
        source_coverage=coverage,
    )
    report.rows.extend(seeded)
    report.rows.extend(trimmed_findings or [])

    for identity in sorted(untrimmed_index):
        source = untrimmed_index[identity]
        trimmed = trimmed_index.get(identity)
        if trimmed is None or not trimmed.has_cram:
            report.rows.append(
                _row(
                    "not_trimmed",
                    wafer=source.wafer,
                    ug=source.ug,
                    source_key=source.cram_key,
                    detail=(
                        "no trimmed CRAM for this well in the upload"
                        if trimmed is None
                        else f"well present as {trimmed.stem} but its CRAM is absent"
                    ),
                )
            )
            continue

        report.matched += 1
        report.rows.extend(_reconcile_identity(source, trimmed))
        size_row = _reconcile_size(source, trimmed)
        if size_row is not None:
            report.rows.append(size_row)
        read_row = _reconcile_read_counts(source, trimmed, fail_counts)
        if read_row is not None:
            report.rows.append(read_row)
        for entry in report.source_coverage:
            if entry.source_uri == source.source_uri:
                entry.matched += 1

    sourced_wafers = {e.wafer for e in untrimmed_index.values()}
    report.unsourced_wafers = sorted(
        {t.wafer for t in trimmed_index.values() if t.wafer not in sourced_wafers}
    )

    for identity in sorted(trimmed_index):
        if identity in untrimmed_index:
            continue
        trimmed = trimmed_index[identity]
        report.rows.append(
            _row(
                "orphan_trimmed",
                wafer=trimmed.wafer,
                ug=trimmed.ug,
                sublibrary=trimmed.sublibrary,
                trimmed_stem=trimmed.stem,
                detail="trimmed well has no matching vendor CRAM in the source prefix",
            )
        )
    return report
