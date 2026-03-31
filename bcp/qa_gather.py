"""
QA data gathering: S3 traversal and file collection for the QA pipeline.

Extracts the data-gathering logic (notebook Cell 2) into testable classes so
the notebook shrinks to a single function call.  Pure validation lives in
qa_checks; parsing in qa_mods.
"""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from qa_mods import (
    QARunContext,
    extract_run_id_from_trimmer_filename,
    grab_trimmer_stats,
    ingest_merged_trimmer_from_s3,
    is_order_level_processed_folder,
    is_valid_cellranger_run_dir_name,
    load_files_from_manifest,
    parse_raw_filename,
    parse_scale_samples_csv,
    parse_scale_workflow_info,
    valid_assays,
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
    exp_to_run_map: dict[str, str] = field(default_factory=dict)
    merged_wafer_stats: dict[str, dict] = field(default_factory=dict)
    scale_workflow_infos: dict[str, dict] = field(default_factory=dict)
    scale_samples_info: dict[str, dict] = field(default_factory=dict)
    scale_proc_files: dict[str, list[str]] = field(default_factory=dict)
    has_raw: bool = False
    has_processed: bool = False
    gathering_errors: list[str] = field(default_factory=list)
    gathering_warnings: list[str] = field(default_factory=list)


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

    def gather(self) -> QAGatheredData:
        """Run the full gathering pipeline and return collected data."""
        if self.ctx.data_source == "manifest":
            self._gather_from_manifest()
        elif self.ctx.data_source == "s3":
            self._gather_from_s3()
        else:
            raise ValueError(f"Invalid data_source: {self.ctx.data_source}")
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

    # ------------------------------------------------------------------
    # S3 mode — top-level
    # ------------------------------------------------------------------

    def _gather_from_s3(self) -> None:
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

        if self.raw_assay in ("10x", "10x_viral_ORF"):
            self._gather_order_level_merged_trimmers()

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
                r_run = self.s3.list_objects(
                    Bucket=self.bucket, Prefix=run, Delimiter="/"
                )
                if "Contents" in r_run:
                    raw_files.extend([c["Key"] for c in r_run["Contents"]])
            is_10x = False

        # Flat files take precedence (preserves original notebook behaviour)
        if r_raw["Contents"]:
            raw_files = [c["Key"] for c in r_raw["Contents"]]
            is_10x = True

        if not raw_files:
            return

        self._data.has_raw = True
        self._data.all_raw_files.extend(raw_files)

        for rf in raw_files:
            self._process_raw_file(rf, group_name, is_10x)

    def _process_raw_file(self, rf: str, group_name: str, is_10x: bool) -> None:
        parsed = parse_raw_filename(rf, self.raw_assay)
        if parsed is not None:
            _run, group, assay, _ug, _barcode = parsed
            if assay not in valid_assays:
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
                fl[group].setdefault(assay, []).append(rf.split("/")[-1])

        if rf.endswith("fastq.gz-metadata.json") and not rf.endswith(
            "_sample.fastq.gz-metadata.json"
        ):
            self._download_metadata_json(rf)
        elif (
            rf.endswith(".cram-metadata.json")
            and self.raw_assay == "scale"
            and "-unmatched.cram-metadata.json" not in rf
        ):
            self._download_metadata_json(rf)
        elif rf.endswith(
            ("trimmer-failure_codes.csv", "trimmer-failure-codes.csv")
        ) and not rf.endswith("merged_trimmer-failure_codes.csv"):
            self._download_trimmer_stats(rf)
        elif _is_merged_trimmer_file(rf):
            ingest_merged_trimmer_from_s3(
                self.bucket, rf, self._data.merged_wafer_stats, self.s3
            )

    def _should_count_for_fastq_log(self, rf: str, assay: str) -> bool:
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
                    or "merged_trimmer-stats.csv" in name
                ):
                    ingest_merged_trimmer_from_s3(
                        self.bucket,
                        key,
                        self._data.merged_wafer_stats,
                        self.s3,
                    )

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

    def _download_trimmer_stats(self, rf: str) -> None:
        exp = "/".join(rf.split("/")[1:3])
        run_id = extract_run_id_from_trimmer_filename(rf)
        if run_id is not None:
            self._data.exp_to_run_map[exp] = run_id
        with tempfile.NamedTemporaryFile(mode="w+b", delete=False, suffix=".csv") as tf:
            local = tf.name
        try:
            self.s3.download_file(self.bucket, rf, local)
            grab_trimmer_stats(self._data.trimmer_failure_stats, exp, local)
        finally:
            Path(local).unlink(missing_ok=True)


def _is_merged_trimmer_file(rf: str) -> bool:
    name = rf.split("/")[-1]
    return (
        "merged_trimmer-failure_codes.csv" in name
        or "_merged_trimmer-failure_codes.csv" in name
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
