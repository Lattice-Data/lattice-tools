from __future__ import annotations

import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

from tqdm import tqdm

from .constants import (
    DEFAULT_H5_TARGET_FILENAME,
    H5_BASE_COLUMNS,
    H5_GENOME_COLUMN,
    H5_INTROSPECT_COLUMNS,
    H5_METRICS_COLUMNS,
)
from .h5_introspect import introspect_h5
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_text,
    list_objects_with_size,
    s3_uri_for,
)
from .tsv_writer import TsvWriter


def extract_sample_name(key: str, per_sample_outs_prefix: str) -> str:
    """First path segment after the prefix (prefix-relative)."""
    base = per_sample_outs_prefix.rstrip("/") + "/"
    if not key.startswith(base):
        return ""
    remainder = key[len(base) :]
    return remainder.split("/", 1)[0] if "/" in remainder else ""


def extract_library(key: str) -> str:
    """Best-effort library/GEM-well folder: segment before 'processed' or 'outs'."""
    parts = key.split("/")
    for anchor in ("processed", "outs"):
        if anchor in parts:
            idx = parts.index(anchor)
            if idx > 0:
                return parts[idx - 1]
    return ""


def parse_metrics_cells(bucket: str, key: str, s3_client: Any) -> int:
    """
    Parse sibling per_sample_outs/<sample>/metrics_summary.csv and return cell count.
    """
    if "per_sample_outs/" not in key:
        raise RuntimeError("key not under per_sample_outs; metrics cross-check N/A")
    head, _, tail = key.partition("per_sample_outs/")
    sample = tail.split("/", 1)[0]
    metrics_key = f"{head}per_sample_outs/{sample}/metrics_summary.csv"

    text = get_object_text(s3_client, bucket, metrics_key)
    fields = csv.DictReader(io.StringIO(text)).fieldnames or []

    if "Metric Name" in fields and "Metric Value" in fields:
        for row in csv.DictReader(io.StringIO(text)):
            if (row.get("Metric Name") or "").strip() == "Cells":
                return int(
                    float((row.get("Metric Value") or "").replace(",", "").strip())
                )
        raise RuntimeError("'Cells' metric not found in metrics_summary.csv")

    for row in csv.DictReader(io.StringIO(text)):
        for col in ("Estimated Number of Cells", "Cells"):
            val = row.get(col)
            if val:
                return int(float(val.replace(",", "").strip()))
        break
    raise RuntimeError("could not locate a cell-count column in metrics_summary.csv")


def parse_metrics_cells_from_text(text: str) -> int:
    """Parse metrics_summary.csv content (for unit tests)."""
    fields = csv.DictReader(io.StringIO(text)).fieldnames or []

    if "Metric Name" in fields and "Metric Value" in fields:
        for row in csv.DictReader(io.StringIO(text)):
            if (row.get("Metric Name") or "").strip() == "Cells":
                return int(
                    float((row.get("Metric Value") or "").replace(",", "").strip())
                )
        raise RuntimeError("'Cells' metric not found in metrics_summary.csv")

    for row in csv.DictReader(io.StringIO(text)):
        for col in ("Estimated Number of Cells", "Cells"):
            val = row.get(col)
            if val:
                return int(float(val.replace(",", "").strip()))
        break
    raise RuntimeError("could not locate a cell-count column in metrics_summary.csv")


def _matches_target_filename(key: str, target_filename: str) -> bool:
    return key.rsplit("/", 1)[-1] == target_filename


def process_one_h5(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    do_introspect: bool,
    do_metrics: bool,
    do_genome: bool,
    retries: int,
) -> dict[str, object]:
    """Enrich a single h5 key with CRC, optional introspection and metrics."""
    result: dict[str, object] = {
        "crc": None,
        "crc_error": "",
        "observation_count": "",
        "feature_counts": "",
        "feature_count_total": "",
        "gene_counts_by_genome": "",
        "h5_error": "",
        "metrics_cells": "",
        "metrics_cells_match": "",
        "metrics_error": "",
    }

    crc, crc_err = retry_with_backoff(
        fetch_crc64nvme, s3_client, bucket, key, retries=retries
    )
    result["crc"] = crc
    result["crc_error"] = crc_err or ""

    if do_introspect:
        intro, h5_err = retry_with_backoff(introspect_h5, bucket, key, retries=retries)
        if h5_err:
            result["h5_error"] = h5_err
        else:
            obs, type_counts, genome_counts = intro  # type: ignore[misc]
            result["observation_count"] = obs
            fc = [
                {"feature_type": k, "feature_count": v}
                for k, v in sorted(type_counts.items())
            ]
            result["feature_counts"] = json.dumps(fc)
            result["feature_count_total"] = sum(type_counts.values())
            if do_genome and genome_counts is not None:
                result["gene_counts_by_genome"] = json.dumps(genome_counts)

    if do_metrics:
        cells, m_err = retry_with_backoff(
            parse_metrics_cells, bucket, key, s3_client, retries=retries
        )
        if m_err:
            result["metrics_error"] = m_err
        elif cells is not None:
            result["metrics_cells"] = cells
            if result["observation_count"] != "":
                result["metrics_cells_match"] = str(
                    cells == result["observation_count"]
                )

    return result


def h5_columns(
    *,
    do_introspect: bool,
    do_genome: bool,
    do_metrics: bool,
) -> list[str]:
    cols = list(H5_BASE_COLUMNS)
    if do_introspect:
        cols.extend(H5_INTROSPECT_COLUMNS)
        if do_genome:
            cols.append(H5_GENOME_COLUMN)
    if do_metrics:
        cols.extend(H5_METRICS_COLUMNS)
    cols.append("crc_error")
    if do_introspect:
        cols.append("h5_error")
    if do_metrics:
        cols.append("metrics_error")
    return cols


def default_h5_output_name(prefix: str) -> str:
    segments = prefix.rstrip("/").split("/")
    run_or_dir = (
        segments[-2] if len(segments) >= 2 else (segments[-1] if segments else "output")
    )
    return f"{run_or_dir}_h5_info.tsv"


def extract_h5(
    s3_client: Any,
    bucket: str,
    prefix: str,
    output_path: str,
    *,
    target_filename: str = DEFAULT_H5_TARGET_FILENAME,
    do_introspect: bool = True,
    do_genome: bool = False,
    do_metrics: bool = False,
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
) -> RunSummary:
    """List h5 matrices, enrich metadata, write TSV."""
    targets = list_objects_with_size(
        s3_client,
        bucket,
        prefix,
        predicate=lambda k: _matches_target_filename(k, target_filename),
    )
    summary = RunSummary(total=len(targets))
    if not targets:
        return summary

    columns = h5_columns(
        do_introspect=do_introspect,
        do_genome=do_genome,
        do_metrics=do_metrics,
    )
    writer = TsvWriter(output_path, columns)
    size_by_key = {obj.key: obj.size_bytes for obj in targets}
    default_workers = 16 if do_introspect else 64
    max_workers = min(workers or default_workers, len(targets))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_one_h5,
                s3_client,
                bucket,
                obj.key,
                do_introspect=do_introspect,
                do_metrics=do_metrics,
                do_genome=do_genome,
                retries=retries,
            ): obj.key
            for obj in targets
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(targets), desc="Processing")

        for fut in iterator:
            key = futures[fut]
            r = fut.result()
            row: list[object] = [
                extract_library(key),
                extract_sample_name(key, prefix),
                s3_uri_for(bucket, key),
                size_by_key[key],
                r["crc"] if r["crc"] is not None else "",
            ]
            if do_introspect:
                row.extend(
                    [
                        r["observation_count"],
                        r["feature_counts"],
                        r["feature_count_total"],
                    ]
                )
                if do_genome:
                    row.append(r["gene_counts_by_genome"])
            if do_metrics:
                row.extend([r["metrics_cells"], r["metrics_cells_match"]])
            row.append(r["crc_error"])
            if do_introspect:
                row.append(r["h5_error"])
            if do_metrics:
                row.append(r["metrics_error"])

            writer.append_row(row)

            crc_err = str(r["crc_error"])
            h5_err = str(r["h5_error"]) if do_introspect else ""
            if not crc_err:
                summary.crc_ok += 1
            if do_introspect and not h5_err and not crc_err:
                summary.enrichment_ok += 1
            if crc_err or (do_introspect and h5_err):
                summary.failures.append((key, crc_err, h5_err))

    return summary
