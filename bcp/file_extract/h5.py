from __future__ import annotations

import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Sequence

from tqdm import tqdm

from .constants import (
    CR_TO_LATTICE_FEATURE_TYPE,
    DEFAULT_H5_TARGET_FILENAME,
    H5_BASE_COLUMNS,
    H5_GENOME_COLUMN,
    H5_INTROSPECT_COLUMNS,
    H5_METRICS_COLUMNS,
)
from .fastq import is_target_file
from .h5_introspect import introspect_h5
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_text,
    list_objects_with_size,
    s3_uri_for,
)
from .sheets import LabIdentity
from .tsv_writer import TsvWriter


def extract_sample_name(key: str) -> str:
    """Directory immediately under ``per_sample_outs/``, or empty."""
    _head, marker, tail = key.partition("per_sample_outs/")
    if not marker:
        return ""
    sample, slash, _rest = tail.partition("/")
    if not slash or not sample:
        return ""
    return sample


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


def map_feature_counts(
    raw_counts: dict[str, int],
) -> tuple[list[dict[str, int | str]], dict[str, int]]:
    """Map Cell Ranger feature_type counts to the Lattice feature_counts schema."""
    lattice_sums: dict[str, int] = {}
    unmapped: dict[str, int] = {}

    for raw_type, count in raw_counts.items():
        lattice_type = CR_TO_LATTICE_FEATURE_TYPE.get(raw_type.strip())
        if lattice_type is not None:
            lattice_sums[lattice_type] = lattice_sums.get(lattice_type, 0) + count
        else:
            unmapped[raw_type] = count

    lattice_fc = [
        {"feature_type": feature_type, "feature_count": lattice_sums[feature_type]}
        for feature_type in sorted(lattice_sums)
    ]
    return lattice_fc, unmapped


def h5_worker_ceiling(*, do_introspect: bool, workers: int | None = None) -> int:
    """Thread count, also the S3 connection-pool size those threads share.

    Introspection opens each h5 with many range reads, so 8 stays inside one
    pool without flooding S3. Checksum-only is one small request per file.
    """
    if workers is not None:
        return workers
    return 8 if do_introspect else 64


def process_one_h5(
    s3_client: Any,
    bucket: str,
    key: str,
    *,
    do_introspect: bool,
    do_metrics: bool,
    do_genome: bool,
    retries: int,
    pool_size: int,
) -> dict[str, object]:
    """Enrich a single h5 key with CRC, optional introspection and metrics."""
    result: dict[str, object] = {
        "crc": None,
        "crc_error": "",
        "observation_count": "",
        "feature_counts": "",
        "feature_count_total": "",
        "unmapped_feature_types": "",
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
        intro, h5_err = retry_with_backoff(
            introspect_h5, bucket, key, pool_size=pool_size, retries=retries
        )
        if h5_err:
            result["h5_error"] = h5_err
        else:
            obs, type_counts, genome_counts = intro  # type: ignore[misc]
            result["observation_count"] = obs
            lattice_fc, unmapped = map_feature_counts(type_counts)
            result["feature_counts"] = json.dumps(lattice_fc)
            result["feature_count_total"] = sum(type_counts.values())
            result["unmapped_feature_types"] = json.dumps(unmapped) if unmapped else ""
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


def raw_prefix_for_h5(key: str) -> str | None:
    """Sibling ``raw/`` of the first ``processed/`` segment, or None.

    ``proj/order/LIB/processed/.../matrix.h5`` maps to ``proj/order/LIB/raw/``.
    The full parent path is the grouping key, so two folders that share a
    basename stay separate. A key with no directory before ``processed``
    returns None.
    """
    parts = key.split("/")
    try:
        idx = parts.index("processed")
    except ValueError:
        return None
    parent = "/".join(parts[:idx])
    if not parent:
        return None
    return f"{parent}/raw/"


def fastq_aliases_for_raw_prefix(
    s3_client: Any,
    bucket: str,
    raw_prefix: str,
    namespace: str,
) -> list[str]:
    """Sorted unique ``{lab}:{filename}`` aliases of selected FASTQs."""
    objects = list_objects_with_size(
        s3_client,
        bucket,
        raw_prefix,
        predicate=lambda key: is_target_file(key, require_raw=True),
    )
    return sorted({f"{namespace}:{obj.key.rsplit('/', 1)[-1]}" for obj in objects})


def missing_processed_warning(key: str) -> str:
    return f"derived_from is empty for {key!r}; key has no processed/ segment"


def empty_raw_fastq_warning(raw_prefix: str) -> str:
    return (
        f"derived_from is empty for h5 files paired with {raw_prefix!r}; "
        "no selected FASTQs in that directory"
    )


def derived_from_cells(
    s3_client: Any,
    bucket: str,
    keys: Sequence[str],
    namespace: str,
    summary: RunSummary,
) -> dict[str, str]:
    """JSON alias lists for each h5 key, listing each sibling raw/ once."""
    aliases_by_raw: dict[str, list[str]] = {}
    warned_raw: set[str] = set()
    cells: dict[str, str] = {}
    for key in keys:
        raw_prefix = raw_prefix_for_h5(key)
        if raw_prefix is None:
            summary.warnings.append(missing_processed_warning(key))
            cells[key] = json.dumps([])
            continue
        if raw_prefix not in aliases_by_raw:
            aliases_by_raw[raw_prefix] = fastq_aliases_for_raw_prefix(
                s3_client, bucket, raw_prefix, namespace
            )
        aliases = aliases_by_raw[raw_prefix]
        if not aliases and raw_prefix not in warned_raw:
            warned_raw.add(raw_prefix)
            summary.warnings.append(empty_raw_fastq_warning(raw_prefix))
        cells[key] = json.dumps(aliases)
    return cells


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
    order_name = prefix.rstrip("/").rsplit("/", 1)[-1] if prefix else "output"
    return f"{order_name}_h5_info.tsv"


def extract_h5(
    s3_client: Any,
    bucket: str,
    prefix: str,
    output_path: str,
    *,
    lab: str,
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

    namespace = LabIdentity.parse(lab).name
    derived = derived_from_cells(
        s3_client,
        bucket,
        [obj.key for obj in targets],
        namespace,
        summary,
    )
    columns = h5_columns(
        do_introspect=do_introspect,
        do_genome=do_genome,
        do_metrics=do_metrics,
    )
    writer = TsvWriter(output_path, columns)
    size_by_key = {obj.key: obj.size_bytes for obj in targets}
    pool_size = h5_worker_ceiling(do_introspect=do_introspect, workers=workers)
    max_workers = min(pool_size, len(targets))

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
                pool_size=pool_size,
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
                extract_sample_name(key),
                s3_uri_for(bucket, key),
                size_by_key[key],
                r["crc"] if r["crc"] is not None else "",
                derived[key],
            ]
            if do_introspect:
                row.extend(
                    [
                        r["observation_count"],
                        r["feature_counts"],
                        r["feature_count_total"],
                        r["unmapped_feature_types"],
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
