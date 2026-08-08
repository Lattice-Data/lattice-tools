"""
QA data gathering: S3 traversal and file collection for the QA pipeline.

Extracts the data-gathering logic (notebook Cell 2) into testable classes so
the notebook shrinks to a single function call.  Pure validation lives in
qa_checks; parsing in qa_mods.
"""

from __future__ import annotations

import json
import tempfile

import s3fs
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import Any
import re

from qa_mods import (
    QARunContext,
    finalize_merged_wafer_stats,
    apply_seahub_trim_fail_blocks,
    parse_seahub_trim_fail_csv,
    grab_trimmer_stats,
    grab_trimmer_failure_codes_wafer_metrics,
    ingest_merged_trimmer_from_s3,
    merge_partial_wafer_stats,
    parse_seahub_raw_path,
    seahub_file_stem,
    seahub_stem_and_family,
    seahub_trimmer_failure_storage_key,
    seahub_trimmer_group_storage_key,
    trimmer_failure_storage_key,
    trimmer_group_storage_key,
    is_order_level_processed_folder,
    is_valid_cellranger_run_dir_name,
    load_files_from_manifest,
    parse_pct_pf_q30_from_lines,
    parse_raw_filename,
    parse_scale_samples_csv,
    parse_scale_workflow_info,
    valid_assays,
)
from qa_constants import (
    SEAHUB_FAIL_SUFFIXES,
    SEAHUB_PLATE_SIZES,
    SEAHUB_SUBLIBRARY_TYPES,
    SEAHUB_UNTYPED_LABEL,
)
from qa_seahub_sop import sop_violation_summary, validate_seahub_stems

METADATA_DOWNLOAD_MAX_WORKERS = 16
TRIM_FAIL_DOWNLOAD_MAX_WORKERS = 16
METADATA_DOWNLOAD_PROGRESS_INTERVAL = 250
PCT_Q30_READ_MAX_WORKERS = 16
PCT_Q30_READ_PROGRESS_INTERVAL = 250
_RAW_ASSAYS_SKIP_PCT_Q30_ERRORS = frozenset(
    {"scale", "sci_jumbo", "sci_plex", "seahub_sci"}
)


def _normalize_group_id_for_compare(group_id: str) -> str:
    """Normalize group IDs for equality checks.

    Group identifiers may be represented with either hyphens or underscores
    across directory names and raw filenames (for example ``CUIMC-001`` vs
    ``CUIMC_001``). For validation, treat these as equivalent while preserving
    original values in logs and gathered output.
    """
    return group_id.replace("-", "_").upper()


@dataclass
class QAGatheredData:
    """All data collected during the gathering phase of QA."""

    all_raw_files: list[str] = field(default_factory=list)
    all_proc_files: dict[str, list[str]] = field(default_factory=dict)
    fastq_log: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    read_metadata: dict[str, dict] = field(default_factory=dict)
    trimmer_failure_stats: dict[str, dict[str, list[float]]] = field(
        default_factory=dict
    )
    group_failure_stats: dict[str, dict[str, list[float]]] = field(default_factory=dict)
    exp_to_run_map: dict[str, str] = field(default_factory=dict)
    merged_wafer_stats: dict[str, dict] = field(default_factory=dict)
    scale_workflow_infos: dict[str, dict] = field(default_factory=dict)
    scale_samples_info: dict[str, dict] = field(default_factory=dict)
    scale_proc_files: dict[str, list[str]] = field(default_factory=dict)
    pct_q30_values: dict[str, float] = field(default_factory=dict)
    has_raw: bool = False
    has_processed: bool = False
    gathering_errors: list[str] = field(default_factory=list)
    gathering_warnings: list[str] = field(default_factory=list)
    # SeaHub-only: every wafer folder seen under raw/{sublibrary}/{wafer}/, so
    # the wafer table can list wafers that produced no parsable failure CSV.
    discovered_wafers: set[str] = field(default_factory=set)
    # SeaHub-only: SOP path/filename violations, one entry per stem+rule.
    sop_violations: list[dict[str, str]] = field(default_factory=list)
    # SeaHub-only: absolute per-format trimmer failed/total counts keyed by stem.
    # Keyed by (raw_dir, stem): a bare stem is not unique across folders.
    seahub_fail_counts: dict[tuple[str, str], dict[str, dict[str, int]]] = field(
        default_factory=dict
    )
    # SeaHub-only, S3 mode only: object sizes by key, for the trimmed-vs-vendor
    # size check.  Empty in manifest mode, where the manifest supplies only URIs;
    # a missing size reads as "unknown" and the check stays silent.
    raw_file_sizes: dict[str, int] = field(default_factory=dict)


class QADataGatherer:
    """Collects QA data from S3 or a manifest file.

    Usage::

        gatherer = QADataGatherer(ctx, s3_client)
        data = gatherer.gather()
    """

    def __init__(self, ctx: QARunContext, s3_client: Any) -> None:
        self.ctx = ctx
        self.s3 = s3_client
        self.paginator = s3_client.get_paginator("list_objects")
        self.bucket = ctx.bucket
        self.raw_assay = ctx.raw_assay
        self._data = QAGatheredData()
        self._metadata_lock = Lock()
        self._pct_q30_lock = Lock()
        # Collected during the walk and fetched as one batch afterwards; see
        # _download_seahub_trim_fail_batch.
        self._seahub_trim_fail_pending: list[str] = []
        s3_fs_kwargs: dict[str, Any] = {}
        region = getattr(getattr(s3_client, "meta", None), "region_name", None)
        if region:
            s3_fs_kwargs["client_kwargs"] = {"region_name": region}
        self._s3_fs = s3fs.S3FileSystem(**s3_fs_kwargs)

    def gather(self) -> QAGatheredData:
        """Run the full gathering pipeline and return collected data."""
        if self.ctx.allow_truncated_stats_name:
            self._data.gathering_warnings.append(
                "RELAXED NAMING: allow_truncated_stats_name=True; "
                "*_stats.csv accepted as alias of *_trimmer-stats.csv"
            )
        if self.ctx.data_source == "manifest":
            self._gather_from_manifest()
        elif self.ctx.data_source == "s3":
            self._gather_from_s3()
        else:
            raise ValueError(f"Invalid data_source: {self.ctx.data_source}")
        self._gather_pct_q30()
        finalize_merged_wafer_stats(self._data.merged_wafer_stats)
        return self._data

    # ------------------------------------------------------------------
    # Manifest mode
    # ------------------------------------------------------------------

    def _gather_from_manifest(self) -> None:
        all_raw, all_proc, _bucket = load_files_from_manifest(
            manifest_path=self.ctx.manifest_path,
            delimiter=self.ctx.manifest_delimiter,
            s3_column=self.ctx.manifest_s3_column,
            has_header=self.ctx.manifest_has_header,
        )
        if _bucket != self.ctx.bucket:
            raise ValueError(
                f"Manifest bucket {_bucket!r} does not match "
                f"context ({self.ctx.bucket!r})."
            )
        self._data.all_raw_files = all_raw
        self._data.all_proc_files = all_proc
        self._data.has_raw = len(all_raw) > 0
        self._data.has_processed = len(all_proc) > 0
        if self.raw_assay == "seahub_sci" and self._data.has_raw:
            self._enrich_seahub_raw_files(self.ctx.order)

    # ------------------------------------------------------------------
    # S3 mode — top-level
    # ------------------------------------------------------------------

    def _gather_from_s3(self) -> None:
        if self.raw_assay == "seahub_sci":
            self._gather_from_s3_seahub()
            return
        o = self.ctx.listing_prefix
        r_order = self.s3.list_objects(Bucket=self.bucket, Prefix=o, Delimiter="/")

        if "CommonPrefixes" not in r_order:
            return

        groups = [e["Prefix"] for e in r_order["CommonPrefixes"]]
        for g in groups:
            if is_order_level_processed_folder(o, g):
                self._data.gathering_warnings.append(
                    f"WARNING: order-level `processed/` folder found at `{g}` — skipping."
                )
                continue

            group_name = g.replace(o, "").rstrip("/")
            self._data.fastq_log[group_name] = {}

            r_group = self.s3.list_objects(Bucket=self.bucket, Prefix=g, Delimiter="/")
            subdirs = [e["Prefix"] for e in r_group.get("CommonPrefixes", [])]

            if len(subdirs) > 2:
                self._data.gathering_warnings.append(f"EXTRA subdirs {g}")

            if f"{g}raw/" not in subdirs:
                self._data.gathering_warnings.append(f"raw/ MISSING {g}")
            else:
                self._gather_group_raw(g, group_name)

            if f"{g}processed/" not in subdirs:
                self._data.gathering_warnings.append(f"processed/ MISSING {g}")
            else:
                self._gather_group_processed(g, group_name)

        if self.raw_assay in ("10x", "10x_cram", "10x_viral_ORF"):
            self._gather_order_level_merged_trimmers()

    def _gather_from_s3_seahub(self) -> None:
        """SeaHub layout: {proj}/{ExperimentID}/raw/{sublibrary}/{wafer}/.

        Everything under ``raw/`` is collected recursively, at whatever depth it
        sits. Walking the delimiter tree and keeping only the wafer-level
        ``Contents`` dropped exactly the objects the SOP wants to hear about: an
        object directly in ``raw/`` or in a sublibrary folder never entered
        ``all_raw_files``, so ``bad_path_depth`` could only ever fire for keys too
        *deep*, never too shallow, and S3 mode disagreed with manifest mode --
        whose rule is simply ``"/raw/" in key`` -- about the bucket's contents.
        The delimiter walk is kept, but only for ``discovered_wafers``, which is a
        question about folders rather than objects.
        """
        o = self.ctx.listing_prefix
        experiment_id = self.ctx.order
        # No fastq_log seed here. _process_raw_file keys that by sublibrary and
        # creates the entry it needs, so an ExperimentID key only ever stayed
        # empty -- and manifest mode, which reaches the same wells, never made
        # one.

        r_top = self._list_folders(o, "the experiment folder")
        top_subdirs = [e["Prefix"] for e in r_top.get("CommonPrefixes", [])]
        if f"{o}processed/" in top_subdirs:
            self._data.has_processed = True
            self._data.gathering_warnings.append(
                "processed/ present for seahub_sci; processed validation is skipped"
            )

        raw_prefix = f"{o}raw/"
        raw_files: list[str] = []
        for page in self.paginator.paginate(Bucket=self.bucket, Prefix=raw_prefix):
            for content in page.get("Contents", []):
                key = content["Key"]
                raw_files.append(key)
                self._data.raw_file_sizes[key] = int(content.get("Size", 0) or 0)

        self._discover_seahub_wafers(raw_prefix)

        if not raw_files:
            # Emptiness is now a statement about objects, not about folders: a
            # raw/ holding only loose objects used to be reported as missing.
            self._data.gathering_warnings.append(f"raw/ MISSING or empty at {o}")
            return

        self._data.has_raw = True
        self._data.all_raw_files = raw_files
        self._enrich_seahub_raw_files(experiment_id)

    def _list_folders(self, prefix: str, what: str) -> dict:
        """One delimiter listing, saying so if S3 cut it short.

        ``list_objects`` returns at most 1000 entries and reports the fact only
        in ``IsTruncated``, which nothing was reading: past the cap the folders
        beyond it simply did not exist as far as QA was concerned, and the run
        still looked clean. Unreachable with anything observed so far -- REF3 has
        7 sublibraries and GENE7 has 9 -- but silence is the wrong failure mode
        for a limit nobody would think to check.

        Only folder enumeration goes through here. Objects are paginated, so a
        truncated listing now costs ``discovered_wafers`` entries or the
        ``processed/`` notice, never a raw file.
        """
        result = self.s3.list_objects(Bucket=self.bucket, Prefix=prefix, Delimiter="/")
        if result.get("IsTruncated"):
            self._data.gathering_warnings.append(
                f"LISTING TRUNCATED: {what} under {prefix} returned only the "
                "first 1000 entries; wafer discovery and the processed/ check "
                "are incomplete for this run"
            )
        return result

    def _discover_seahub_wafers(self, raw_prefix: str) -> None:
        """Record the wafer folders under ``raw/{sublibrary}/``.

        A wafer is a directory, so this is the one part of the walk that has to
        ask about ``CommonPrefixes`` rather than keys: a wafer folder holding no
        objects still counts as discovered.
        """
        r_raw = self._list_folders(raw_prefix, "the sublibrary folders")
        for sublib_entry in r_raw.get("CommonPrefixes", []):
            r_wafers = self._list_folders(sublib_entry["Prefix"], "the wafer folders")
            for wafer_entry in r_wafers.get("CommonPrefixes", []):
                wafer = wafer_entry["Prefix"].rstrip("/").split("/")[-1]
                if wafer:
                    self._data.discovered_wafers.add(wafer)

    def _enrich_seahub_raw_files(self, experiment_id: str) -> None:
        """Parse SeaHub trim artifacts, metadata, and plate-size warnings."""
        metadata_files: list[str] = []
        plate_counts: dict[tuple[str, str], set[str]] = {}
        cram_stems: dict[str, str] = {}
        metadata_stems: set[str] = set()
        wrong_experiment: list[tuple[str, str]] = []

        if not experiment_id:
            self._data.gathering_warnings.append(
                "EXPERIMENT UNKNOWN: no ExperimentID resolved for this run, so the "
                "cross-experiment check did not run; objects from another "
                "ExperimentID would not be reported"
            )

        for rf in self._data.all_raw_files:
            name = rf.split("/")[-1]
            path_info = parse_seahub_raw_path(rf)
            if path_info is not None:
                if experiment_id and path_info["experiment_id"] != experiment_id:
                    wrong_experiment.append((path_info["experiment_id"], rf))
                parsed_stem = seahub_stem_and_family(name)
                if parsed_stem is not None:
                    stem, suffix, _family = parsed_stem
                    key = (path_info["sublibrary"], path_info["wafer"])
                    plate_counts.setdefault(key, set()).add(stem)
                    stem_key = "/".join(rf.split("/")[:-1]) + "/" + stem
                    if suffix.endswith(".cram-metadata.json"):
                        metadata_stems.add(stem_key)
                    elif suffix.endswith(".cram"):
                        # Keep the delivered spelling: both families reach here,
                        # so the warning must not name a ".trim.*" file that a
                        # bare-family well never had.
                        cram_stems[stem_key] = suffix

            if self._should_download_metadata_json(rf):
                metadata_files.append(rf)
                continue
            self._process_raw_file(rf, experiment_id, is_10x=False)

        self._download_metadata_json_batch(metadata_files)
        self._download_seahub_trim_fail_batch(self._seahub_trim_fail_pending)
        self._append_seahub_wrong_experiment_errors(experiment_id, wrong_experiment)
        self._append_seahub_plate_warnings(plate_counts)
        self._append_seahub_missing_metadata_warnings(cram_stems, metadata_stems)
        self._append_seahub_sop_violations()

    def _append_seahub_sop_violations(self) -> None:
        """Validate the listing against the SOP, one entry per stem and rule.

        Kept out of ``gathering_warnings`` deliberately: a wholly misnamed
        upload produces one entry per well, which would bury every other
        warning.  Only a summary line goes to the warnings list; the detail is
        reported as its own table.
        """
        violations = validate_seahub_stems(self.bucket, self._data.all_raw_files)
        if not violations:
            return
        self._data.sop_violations = [v.as_dict() for v in violations]
        summary = sop_violation_summary(violations)
        breakdown = ", ".join(f"{name}={count}" for name, count in summary.items())
        self._data.gathering_warnings.append(
            f"SOP VIOLATIONS: {len(violations)} across "
            f"{len(summary)} rule(s): {breakdown} "
            f"— per-rule examples in the SOP violations cell below; full detail in "
            f"{self.ctx.output_label}_raw_sop_violations.csv"
        )

    def _append_seahub_wrong_experiment_errors(
        self, experiment_id: str, wrong: list[tuple[str, str]]
    ) -> None:
        """Report objects from another ExperimentID once, not once per object.

        A prefix holding the wrong experiment holds *all* of it, so the per-object
        form wrote thousands of copies of one fact into ``{order}_errors.txt`` and
        buried everything else — the same reason
        :meth:`_append_seahub_sop_violations` summarises rather than enumerates.
        The count and the offending ExperimentIDs are what identify the mixup;
        two keys are enough to point at it.
        """
        if not wrong:
            return
        found = sorted({eid for eid, _ in wrong})
        examples = [key for _, key in wrong[:2]]
        detail = "; ".join(examples)
        if len(wrong) > len(examples):
            detail += f"; and {len(wrong) - len(examples)} more"
        self._data.gathering_errors.append(
            f"WRONG EXPERIMENT: {len(wrong)} object(s) belong to "
            f"{', '.join(found)}, not {experiment_id}: {detail}"
        )

    def _append_seahub_plate_warnings(
        self, plate_counts: dict[tuple[str, str], set[str]]
    ) -> None:
        for (sublibrary, wafer), stems in sorted(plate_counts.items()):
            count = len(stems)
            if count not in SEAHUB_PLATE_SIZES:
                self._data.gathering_warnings.append(
                    f"PLATE SIZE: {sublibrary}/{wafer} has {count} wells "
                    f"(expected one of {sorted(SEAHUB_PLATE_SIZES)})"
                )

    def _append_seahub_missing_metadata_warnings(
        self, cram_stems: dict[str, str], metadata_stems: set[str]
    ) -> None:
        """Warn once for the CRAMs lacking their metadata sidecar.

        The metadata sidecar is classified as optional (its absence is not a
        hard inventory failure), but a missing sidecar means read-count QA
        cannot run for that well, so surface it as an informational warning.

        One warning, not one per well. CZI generates these sidecars for the
        whole upload, so when they are absent they are absent for every well at
        once: the per-well form printed 336 identical lines on REF3 and 864 on
        GENE7, which is the "one fact reported N times" the SOP ``scope``
        machinery exists to prevent.

        ``cram_stems`` carries the suffix each well was delivered with, so the
        examples name ``.cram`` / ``.cram-metadata.json`` for a bare-family well
        rather than ``.trim.*`` files it never had. Both spellings can appear in
        one upload, so the count is broken down by suffix.
        """
        missing = sorted(set(cram_stems) - metadata_stems)
        if not missing:
            return
        by_suffix = Counter(cram_stems[stem] for stem in missing)
        breakdown = ", ".join(
            f"{count}x {suffix}-metadata.json"
            for suffix, count in sorted(by_suffix.items())
        )
        examples = "; ".join(f"{stem}{cram_stems[stem]}" for stem in missing[:2])
        if len(missing) > 2:
            examples += f"; and {len(missing) - 2} more"
        self._data.gathering_warnings.append(
            f"METADATA MISSING: {len(missing)} CRAM(s) have no matching "
            f"metadata sidecar ({breakdown}); read-count QA cannot run for "
            f"them: {examples}"
        )

    # ------------------------------------------------------------------
    # S3 mode — raw files
    # ------------------------------------------------------------------

    def _gather_group_raw(self, group_prefix: str, group_name: str) -> None:
        r_raw: dict[str, list] = {"Contents": [], "CommonPrefixes": []}
        for page in self.paginator.paginate(
            Bucket=self.bucket, Prefix=f"{group_prefix}raw/", Delimiter="/"
        ):
            r_raw["Contents"].extend(page.get("Contents", []))
            r_raw["CommonPrefixes"].extend(page.get("CommonPrefixes", []))

        raw_files: list[str] = []
        is_10x = False

        if r_raw["CommonPrefixes"]:
            non_10x_runs = [e["Prefix"] for e in r_raw["CommonPrefixes"]]
            for run in non_10x_runs:
                for page in self.paginator.paginate(
                    Bucket=self.bucket, Prefix=run, Delimiter="/"
                ):
                    raw_files.extend([c["Key"] for c in page.get("Contents", [])])
            is_10x = False

        # Flat files take precedence (preserves original notebook behaviour)
        if r_raw["Contents"]:
            raw_files = [c["Key"] for c in r_raw["Contents"]]
            is_10x = True

        if not raw_files:
            return

        self._data.has_raw = True
        self._data.all_raw_files.extend(raw_files)

        metadata_files: list[str] = []
        for rf in raw_files:
            if self._should_download_metadata_json(rf):
                metadata_files.append(rf)
                continue
            self._process_raw_file(rf, group_name, is_10x)
        self._download_metadata_json_batch(metadata_files)

    def _process_raw_file(self, rf: str, group_name: str, is_10x: bool) -> None:
        parsed = parse_raw_filename(rf, self.raw_assay)
        if parsed is not None:
            _run, group, assay, _ug, _barcode = parsed
            # SeaHub allows a narrower vocabulary than valid_assays, and a
            # missing type is already reported as invalid_sublibrary_type by the
            # SOP validator, so it must not also raise a WRONG ASSAY error.
            if self.raw_assay == "seahub_sci":
                allowed_assays: tuple[str, ...] | list[str] = SEAHUB_SUBLIBRARY_TYPES
                assay_is_wrong = assay is not None and assay not in allowed_assays
            else:
                assay_is_wrong = assay not in valid_assays
            if assay_is_wrong:
                self._data.gathering_errors.append(f"WRONG ASSAY: {assay} {rf}")
            if (
                _normalize_group_id_for_compare(group)
                != _normalize_group_id_for_compare(group_name)
                and is_10x
            ):
                self._data.gathering_errors.append(
                    f"WRONG GROUP: {group} {group_name} {rf}"
                )
            if self._should_count_for_fastq_log(rf, assay):
                fl = self._data.fastq_log
                if group not in fl:
                    fl[group] = {}
                assay_key = assay if assay is not None else SEAHUB_UNTYPED_LABEL
                fl[group].setdefault(assay_key, []).append(rf.split("/")[-1])

        if self._should_download_metadata_json(rf):
            self._download_metadata_json(rf)
        elif rf.endswith(
            ("trimmer-failure_codes.csv", "trimmer-failure-codes.csv")
        ) and not rf.endswith("merged_trimmer-failure_codes.csv"):
            self._download_trimmer_failure_codes(rf)
        elif self.raw_assay == "seahub_sci" and self._is_seahub_trim_fail(rf):
            # Deferred, not fetched: one per well, and serially that is one
            # round-trip per well against S3.
            self._seahub_trim_fail_pending.append(rf)
        elif _is_merged_trimmer_file(rf):
            ingest_merged_trimmer_from_s3(
                self.bucket, rf, self._data.merged_wafer_stats, self.s3
            )

    def _is_seahub_trim_fail(self, rf: str) -> bool:
        """Is this a per-well trim failure CSV whose numbers QA should ingest?

        The suffix alone is not enough. ``_fail.csv`` is generic, so the bare
        family is only recognised when the stem before it parses -- the same
        gate ``seahub_stem_and_family`` applies, and the reason a 480-well
        upload did not report as 960. Without it here, a file the SOP table
        calls ``unexpected_suffix`` still fed its numbers into the wafer and
        sublibrary trimmer-fail histograms.

        ``.trim_fail.csv`` stays ungated, deliberately: ``.trim.*`` is
        distinctive enough that a malformed stem is a stem defect rather than a
        different kind of file. ``seahub_file_stem`` already draws that line.
        """
        return rf.endswith(SEAHUB_FAIL_SUFFIXES) and seahub_file_stem(rf) is not None

    def _should_download_metadata_json(self, rf: str) -> bool:
        """Return True for raw metadata sidecars that QA parses."""
        if rf.endswith("fastq.gz-metadata.json") and not rf.endswith(
            "_sample.fastq.gz-metadata.json"
        ):
            return True
        if (
            rf.endswith(".cram-metadata.json")
            and self.raw_assay in ("scale", "10x_cram", "sci_plex", "seahub_sci")
            and "-unmatched.cram-metadata.json" not in rf
            and "_unmatched.cram-metadata.json" not in rf
        ):
            return True
        return False

    def _download_metadata_json_batch(self, metadata_files: list[str]) -> None:
        """Download raw metadata sidecars, using concurrency for large Scale runs."""
        if not metadata_files:
            return

        total = len(metadata_files)
        if total == 1:
            self._download_metadata_json(metadata_files[0])
            return

        if total >= METADATA_DOWNLOAD_PROGRESS_INTERVAL:
            print(f"Downloading {total} raw metadata JSON file(s)...")

        workers = min(METADATA_DOWNLOAD_MAX_WORKERS, total)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_path = {
                executor.submit(self._download_metadata_json, rf): rf
                for rf in metadata_files
            }
            for completed, future in enumerate(as_completed(future_to_path), start=1):
                future.result()
                if total >= METADATA_DOWNLOAD_PROGRESS_INTERVAL and (
                    completed % METADATA_DOWNLOAD_PROGRESS_INTERVAL == 0
                    or completed == total
                ):
                    print(f"  downloaded {completed}/{total} metadata JSON file(s)")

    def _should_count_for_fastq_log(self, rf: str, assay: str | None) -> bool:
        if rf.endswith(".fastq.gz") and not rf.endswith("_sample.fastq.gz"):
            return True
        if rf.endswith(".cram") and assay == "viral_ORF":
            return True
        if (
            rf.endswith(".cram")
            and self.raw_assay == "scale"
            and not rf.endswith("_unmatched.cram")
        ):
            return True
        if rf.endswith(".cram") and self.raw_assay == "10x_cram":
            return False
        # Both SeaHub families count: bare ".cram" here is a trim output whose
        # ".trim" infix was dropped, not vendor data, so excluding it would hide
        # most of a misnamed upload from the per-sublibrary counts.
        if rf.endswith(".cram") and self.raw_assay == "seahub_sci":
            return True
        return False

    # ------------------------------------------------------------------
    # S3 mode — processed files (dispatch)
    # ------------------------------------------------------------------

    def _gather_group_processed(self, group_prefix: str, group_name: str) -> None:
        if self.raw_assay == "scale":
            self._gather_group_processed_scale(group_prefix, group_name)
        else:
            self._gather_group_processed_cellranger(group_prefix, group_name)

    # ------------------------------------------------------------------
    # S3 mode — Scale processed
    # ------------------------------------------------------------------

    def _gather_group_processed_scale(self, group_prefix: str, group_name: str) -> None:
        r_proc = self.s3.list_objects(
            Bucket=self.bucket,
            Prefix=f"{group_prefix}processed/",
            Delimiter="/",
        )
        proc_subdirs = [e["Prefix"] for e in r_proc.get("CommonPrefixes", [])]

        for psd in proc_subdirs:
            r_psd = self.s3.list_objects(Bucket=self.bucket, Prefix=psd, Delimiter="/")
            psd_files = [c["Key"] for c in r_psd.get("Contents", [])]
            wf_keys = [k for k in psd_files if k.endswith("workflow_info.json")]
            if not wf_keys:
                continue

            self._data.has_processed = True

            with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as tf:
                local_wf = tf.name
            try:
                self.s3.download_file(self.bucket, wf_keys[0], local_wf)
                self._data.scale_workflow_infos[group_name] = parse_scale_workflow_info(
                    local_wf
                )
            finally:
                Path(local_wf).unlink(missing_ok=True)

            csv_keys = [k for k in psd_files if k.endswith("samples.csv")]
            if csv_keys:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tf:
                    local_csv = tf.name
                try:
                    self.s3.download_file(self.bucket, csv_keys[0], local_csv)
                    self._data.scale_samples_info[group_name] = parse_scale_samples_csv(
                        local_csv
                    )
                finally:
                    Path(local_csv).unlink(missing_ok=True)

            samples_files: list[str] = []
            for page in self.paginator.paginate(
                Bucket=self.bucket, Prefix=f"{psd}samples/"
            ):
                if "Contents" in page:
                    for obj in page["Contents"]:
                        samples_files.append(obj["Key"])
            self._data.scale_proc_files[group_name] = samples_files

    # ------------------------------------------------------------------
    # S3 mode — CellRanger processed
    # ------------------------------------------------------------------

    def _gather_group_processed_cellranger(
        self, group_prefix: str, group_name: str
    ) -> None:
        proc_prefix = f"{group_prefix}processed/"
        r_proc = self.s3.list_objects(
            Bucket=self.bucket, Prefix=proc_prefix, Delimiter="/"
        )

        if "Contents" in r_proc:
            self._data.gathering_warnings.append(f"UNEXPECTED FILES {proc_prefix}")

        cr_prefix = f"{proc_prefix}cellranger/"
        if cr_prefix not in [e["Prefix"] for e in r_proc.get("CommonPrefixes", [])]:
            self._data.gathering_warnings.append(f"cellranger/ MISSING {group_prefix}")
            return

        self._data.has_processed = True

        r_cr = self.s3.list_objects(Bucket=self.bucket, Prefix=cr_prefix, Delimiter="/")
        if "Contents" in r_cr:
            self._data.gathering_warnings.append(f"UNEXPECTED FILES {cr_prefix}")

        run_dates = [e["Prefix"] for e in r_cr.get("CommonPrefixes", [])]
        for rd in run_dates:
            date = rd.split("/")[-2]
            if not is_valid_cellranger_run_dir_name(date):
                self._data.gathering_errors.append(
                    f"INCORRECT DATE FORMAT: {date} {rd}"
                )

            r_date = self.s3.list_objects(Bucket=self.bucket, Prefix=rd, Delimiter="/")
            if "Contents" in r_date:
                for c in r_date["Contents"]:
                    self._data.gathering_warnings.append(f"UNEXPECTED FILES {c['Key']}")

            outsdirs = [e["Prefix"] for e in r_date.get("CommonPrefixes", [])]
            if len(outsdirs) > 1:
                self._data.gathering_warnings.append(f"EXTRA subdirs {rd}")
            if f"{rd}outs/" not in outsdirs:
                self._data.gathering_warnings.append(f"NO outs/ {rd}")
            else:
                files: list[str] = []
                for page in self.paginator.paginate(
                    Bucket=self.bucket, Prefix=f"{rd}outs/"
                ):
                    if "Contents" in page:
                        for obj in page["Contents"]:
                            files.append(obj["Key"])
                self._data.all_proc_files[group_name] = files

    # ------------------------------------------------------------------
    # S3 mode — order-level merged trimmers (10x)
    # ------------------------------------------------------------------

    def _gather_order_level_merged_trimmers(self) -> None:
        o = self.ctx.listing_prefix
        for page in self.paginator.paginate(Bucket=self.bucket, Prefix=o):
            for obj in page.get("Contents", []):
                key = obj["Key"]
                if not key.startswith(o) or key[len(o) :].count("/") != 0:
                    continue
                name = key.split("/")[-1]
                if (
                    "merged_trimmer-failure_codes.csv" in name
                    or "merged_trimmer-failure-codes.csv" in name
                    or "merged_trimmer-stats.csv" in name
                ):
                    ingest_merged_trimmer_from_s3(
                        self.bucket,
                        key,
                        self._data.merged_wafer_stats,
                        self.s3,
                    )

    # ------------------------------------------------------------------
    # Sublibrary Q30 stats (PCT_PF_Q30_bases)
    # ------------------------------------------------------------------

    def _is_sublibrary_stats_csv(self, key: str) -> bool:
        """Return True for aggregate sublibrary stats CSVs (PCT_PF_Q30_bases source)."""
        name = key.split("/")[-1]
        if not name.endswith(".csv"):
            return False
        if _is_merged_trimmer_file(key):
            return False
        # SeaHub per-well artifacts (".csv" / "_fail.csv", with or without the
        # ".trim" infix) are single-well trimmer output, never sublibrary
        # aggregates, and carry no PCT_PF_Q30_bases.  Compliant names happened to
        # be excluded by the hash_oligo check below, but that only holds for
        # sublibraries whose type token contains "hash_oligo".
        if self.raw_assay == "seahub_sci" and seahub_stem_and_family(name) is not None:
            return False
        if re.search("hash_oligo", name):
            return False
        excluded = (
            "_trimmer-stats.csv",
            "_trimmer-failure_codes.csv",
            "-trimmer-failure-codes.csv",
            "_unmatched.csv",
            "_S1_L001_R1_001.csv",
            "_S1_L001_R2_001.csv",
            "_stats.csv",
            "-stats.csv",
            ".failure_codes.csv",
            ".trimmer_stats.csv",
        )
        if any(name.endswith(suf) for suf in excluded):
            return False
        return parse_raw_filename(key, self.raw_assay) is not None

    def _read_pct_q30_from_s3(self, key: str) -> None:
        """Fetch one sublibrary stats CSV and record its PCT_PF_Q30_bases value."""
        filename = key.split("/")[-1]
        try:
            with self._s3_fs.open(f"{self.bucket}/{key}", "r") as fobj:
                value = parse_pct_pf_q30_from_lines(fobj)
        except Exception as exc:
            if self.raw_assay not in _RAW_ASSAYS_SKIP_PCT_Q30_ERRORS:
                with self._pct_q30_lock:
                    self._data.gathering_errors.append(
                        f"PCT_Q30 read failed for {filename}: {exc}"
                    )
            return

        with self._pct_q30_lock:
            if value is None:
                self._data.gathering_warnings.append(
                    f"PCT_PF_Q30_bases missing in {filename}"
                )
            else:
                self._data.pct_q30_values[filename] = value

    def _gather_pct_q30(self) -> None:
        """Stream sublibrary stats CSVs from S3 and populate pct_q30_values."""
        keys = [
            key
            for key in self._data.all_raw_files
            if self._is_sublibrary_stats_csv(key)
        ]
        if not keys:
            return

        total = len(keys)
        if total == 1:
            self._read_pct_q30_from_s3(keys[0])
            return

        if total >= PCT_Q30_READ_PROGRESS_INTERVAL:
            print(f"Reading {total} sublibrary stats CSV(s) for PCT_PF_Q30_bases...")

        workers = min(PCT_Q30_READ_MAX_WORKERS, total)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(self._read_pct_q30_from_s3, key) for key in keys]
            for completed, future in enumerate(as_completed(futures), start=1):
                future.result()
                if total >= PCT_Q30_READ_PROGRESS_INTERVAL and (
                    completed % PCT_Q30_READ_PROGRESS_INTERVAL == 0
                    or completed == total
                ):
                    print(f"  read {completed}/{total} sublibrary stats CSV(s) for Q30")

    # ------------------------------------------------------------------
    # S3 helpers — file downloads
    # ------------------------------------------------------------------

    def _download_metadata_json(self, rf: str) -> None:
        with tempfile.NamedTemporaryFile(
            mode="w+b", delete=False, suffix=".json"
        ) as tf:
            local = tf.name
        try:
            self.s3.download_file(self.bucket, rf, local)
            with open(local) as fh:
                metadata = json.load(fh)
            reported_filename = str(metadata.get("filename", ""))
            actual_filename = self._canonical_metadata_filename(
                source_key=rf,
                reported_filename=reported_filename,
            )
            metadata["__actual_filename"] = actual_filename
            metadata["__reported_filename"] = reported_filename
            metadata["__source_key"] = rf
            with self._metadata_lock:
                self._data.read_metadata[actual_filename] = metadata
        finally:
            Path(local).unlink(missing_ok=True)

    def _canonical_metadata_filename(
        self, *, source_key: str, reported_filename: str
    ) -> str:
        """Choose a stable metadata key derived from actual S3 object location.

        If the payload reports an s3 URI in ``filename``, keep that shape by
        using ``s3://{bucket}/{source_key}``. Otherwise, use basename form to
        preserve behavior for local/test-style metadata payloads.
        """
        data_key = (
            source_key[: -len("-metadata.json")]
            if source_key.endswith("-metadata.json")
            else source_key
        )
        if reported_filename.startswith("s3://"):
            return f"s3://{self.bucket}/{data_key}"
        return data_key.split("/")[-1]

    def _download_trimmer_failure_codes(self, rf: str) -> None:
        storage_key, run_id = trimmer_failure_storage_key(rf)
        group_key = trimmer_group_storage_key(rf)
        if run_id is not None:
            self._data.exp_to_run_map[run_id] = run_id
            self._data.exp_to_run_map[group_key] = run_id
        with tempfile.NamedTemporaryFile(mode="w+b", delete=False, suffix=".csv") as tf:
            local = tf.name
        try:
            self.s3.download_file(self.bucket, rf, local)
            grab_trimmer_stats(self._data.trimmer_failure_stats, storage_key, local)
            grab_trimmer_stats(self._data.group_failure_stats, group_key, local)
            if run_id is not None:
                rsq_metrics = grab_trimmer_failure_codes_wafer_metrics(local)
                if rsq_metrics:
                    merge_partial_wafer_stats(
                        self._data.merged_wafer_stats, run_id, rsq_metrics
                    )
        finally:
            Path(local).unlink(missing_ok=True)

    def _download_seahub_trim_fail_batch(self, keys: list[str]) -> None:
        """Fetch the per-well SeaHub trim failure CSVs, then fold them in.

        Handles both ``*.trim_fail.csv`` (SOP) and ``*_fail.csv`` (the same
        artifact without the ``.trim`` infix); without the second form, uploads
        that dropped the infix rendered as blank wafer rows.

        SeaHub uploads start from RSQ-passing reads and carry no merged
        trimmer files, so no wafer-level RSQ/TT metrics are derived here; the
        per-well, per-modality trimmer-fail fractions flow into
        ``trimmer_failure_stats`` (keyed by wafer) and ``group_failure_stats``
        (keyed by ExperimentID/sublibrary) for the histograms and the wafer
        sample-level summary via ``build_wafer_failure_stats``.

        There is one of these per well, so fetching them inside the walk meant
        336 sequential round-trips on REF3 and 864 on GENE7 -- while the
        ``.cram-metadata.json`` sidecars sitting beside them went through a
        16-worker pool. They use the same pool now.

        Downloading and parsing happen in the pool; *applying* happens here, in
        listing order. That is what makes the concurrency free of locks: the
        three structures these feed are appended to, so a parallel apply would
        both race and reorder the distributions run to run. The workers touch
        nothing shared -- each returns its blocks and its own warnings, and this
        loop merges them in the order the listing gave, so the output is
        identical to the serial version rather than merely equivalent.
        """
        if not keys:
            return

        total = len(keys)
        if total == 1:
            results = [self._fetch_seahub_trim_fail(keys[0])]
        else:
            if total >= METADATA_DOWNLOAD_PROGRESS_INTERVAL:
                print(f"Downloading {total} per-well trim failure CSV(s)...")
            workers = min(TRIM_FAIL_DOWNLOAD_MAX_WORKERS, total)
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [
                    executor.submit(self._fetch_seahub_trim_fail, rf) for rf in keys
                ]
                results = []
                # In submission order, not as_completed: the apply below has to
                # be ordered anyway, and this makes which failure surfaces first
                # deterministic too.
                for done, future in enumerate(futures, start=1):
                    results.append(future.result())
                    if total >= METADATA_DOWNLOAD_PROGRESS_INTERVAL and (
                        done % METADATA_DOWNLOAD_PROGRESS_INTERVAL == 0 or done == total
                    ):
                        print(f"  read {done}/{total} trim failure CSV(s)")

        for rf, blocks, warnings in results:
            self._apply_seahub_trim_fail(rf, blocks, warnings)
        keys.clear()

    def _fetch_seahub_trim_fail(self, rf: str) -> tuple[str, list[dict], list[str]]:
        """Download and parse one trim failure CSV. Runs in a worker thread.

        Deliberately touches no gatherer state: it returns what it parsed, and
        its own warnings, for the caller to merge in order.
        """
        warnings: list[str] = []
        with tempfile.NamedTemporaryFile(mode="w+b", delete=False, suffix=".csv") as tf:
            local = tf.name
        try:
            self.s3.download_file(self.bucket, rf, local)
            # Parsed once, applied twice. Calling the combined helper for each
            # distribution read_csv'd the same file twice, on an upload that has
            # one of these per well.
            return rf, parse_seahub_trim_fail_csv(local, warnings=warnings), warnings
        finally:
            Path(local).unlink(missing_ok=True)

    def _apply_seahub_trim_fail(
        self, rf: str, blocks: list[dict], warnings: list[str]
    ) -> None:
        """Fold one parsed CSV into both distributions. Single-threaded."""
        storage_key, run_id = seahub_trimmer_failure_storage_key(rf)
        group_key = seahub_trimmer_group_storage_key(rf)
        if run_id is not None:
            self._data.exp_to_run_map[run_id] = run_id
            self._data.exp_to_run_map[group_key] = run_id
        self._data.gathering_warnings.extend(warnings)
        apply_seahub_trim_fail_blocks(
            blocks,
            self._data.trimmer_failure_stats,
            storage_key,
            fail_counts=self._data.seahub_fail_counts,
            # Qualified by folder. On the bare stem, two wells in different
            # sublibrary folders that happen to share one would have merged
            # their per-format counts, and the reconciliation reads this to
            # decide whether the trimmer saw the whole delivered file.
            # index_trimmed_upload computes TrimmedEntry.raw_dir the same way,
            # which is what the lookup uses.
            stem_key=("/".join(rf.split("/")[:-1]), seahub_file_stem(rf)),
        )
        apply_seahub_trim_fail_blocks(blocks, self._data.group_failure_stats, group_key)


def _is_merged_trimmer_file(rf: str) -> bool:
    name = rf.split("/")[-1]
    return (
        "merged_trimmer-failure_codes.csv" in name
        or "merged_trimmer-failure-codes.csv" in name
        or "_merged_trimmer-failure_codes.csv" in name
        or "_merged_trimmer-failure-codes.csv" in name
        or "merged_trimmer-stats.csv" in name
        or "_merged_trimmer-stats.csv" in name
    )


def gather_qa_data(ctx: QARunContext, s3_client: Any) -> QAGatheredData:
    """
    Gather all QA data from S3 or manifest.

    Main entry point: pass a resolved ``QARunContext`` and an S3 client,
    receive a ``QAGatheredData`` with all file listings, metadata, and
    trimmer statistics ready for downstream validation cells.
    """
    return QADataGatherer(ctx, s3_client).gather()
