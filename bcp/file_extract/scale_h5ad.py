"""Extract Scale processed AnnData (.h5ad) metadata from a rundate prefix."""

from __future__ import annotations

import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Callable, Sequence

from tqdm import tqdm

from .constants import SCALE_H5AD_COLUMNS
from .metadata_sheet import fetch_sample_template, parse_sample_template
from .models import ListedObject, RunSummary
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_text,
    list_objects_with_size,
    s3_uri_for,
)
from .scale_wells import (
    ScaleExtractError,
    correlate_samples,
)
from .sheets import LabIdentity
from .tsv_writer import TsvWriter

SAMPLES_CSV_NAME = "samples.csv"
SAMPLES_DIR = "samples"
H5AD_SUFFIX = ".h5ad"
MERGED_ANNDATA_SUFFIX = "_anndata.h5ad"


def format_samples_column(sample_names: Sequence[str], lab: str) -> str:
    """JSON list of sheet sample_name values prefixed with ``{lab}:``."""
    prefix = LabIdentity.parse(lab).name
    return json.dumps([f"{prefix}:{name}" for name in sample_names])


def default_scale_h5ad_output_name(prefix: str) -> str:
    rundate = prefix.rstrip("/").rsplit("/", 1)[-1] or "output"
    return f"{rundate}_scale_h5ad_info.tsv"


def samples_csv_key(prefix: str) -> str:
    return prefix.rstrip("/") + f"/{SAMPLES_CSV_NAME}"


def samples_dir_prefix(prefix: str) -> str:
    return prefix.rstrip("/") + f"/{SAMPLES_DIR}/"


def is_scale_h5ad_key(key: str, rundate_prefix: str) -> bool:
    """True for .h5ad files under {rundate}/samples/ whose name contains QSR."""
    base = samples_dir_prefix(rundate_prefix)
    if not key.startswith(base):
        return False
    name = key[len(base) :]
    return (
        bool(name) and "/" not in name and name.endswith(H5AD_SUFFIX) and "QSR" in name
    )


def sample_from_filename(filename: str) -> str:
    return filename.split(".", 1)[0]


def file_belongs_to_sample(filename: str, sample: str) -> bool:
    return filename == f"{sample}{MERGED_ANNDATA_SUFFIX}" or filename.startswith(
        f"{sample}."
    )


def is_control_file(filename: str, control_samples: set[str] | frozenset[str]) -> bool:
    return any(file_belongs_to_sample(filename, sample) for sample in control_samples)


def parse_samples_csv(text: str) -> list[tuple[str, str]]:
    reader = csv.DictReader(io.StringIO(text))
    fields = reader.fieldnames or []
    if "sample" not in fields or "barcodes" not in fields:
        raise ScaleExtractError("samples.csv must have sample and barcodes columns")
    rows: list[tuple[str, str]] = []
    for row in reader:
        sample = (row.get("sample") or "").strip()
        barcodes = (row.get("barcodes") or "").strip()
        if not sample:
            continue
        rows.append((sample, barcodes))
    return rows


def control_warning(sample: str, barcodes: str) -> str:
    return (
        f"control sample {sample!r} (barcodes {barcodes!r}) has no RT_index "
        "pairing; excluding from output"
    )


def _fetch_crc(
    s3_client: Any, bucket: str, key: str, *, retries: int
) -> tuple[str | None, str]:
    crc, err = retry_with_backoff(
        fetch_crc64nvme, s3_client, bucket, key, retries=retries
    )
    return crc, err or ""


def extract_scale_h5ad(
    s3_client: Any,
    bucket: str,
    prefix: str,
    output_path: str,
    *,
    metadata_gid: str,
    lab: str,
    cro_orders: Sequence[str],
    wafers: Sequence[str],
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
    sheet_csv: str | None = None,
    fetch_sheet: Callable[[str], str] | None = None,
) -> RunSummary:
    """List non-control samples/*.h5ad, fetch CRC64NVME, write the TSV.

    ``cro_orders`` and ``wafers`` are accepted for the Scale CLI contract and
    later scale_cram reuse; they are not written to this TSV.
    """
    del cro_orders, wafers
    lab_name = LabIdentity.parse(lab).name

    csv_key = samples_csv_key(prefix)
    try:
        samples_text = get_object_text(s3_client, bucket, csv_key)
    except Exception as exc:
        raise ScaleExtractError(
            f"could not read {s3_uri_for(bucket, csv_key)}: {exc}"
        ) from exc
    sample_rows = parse_samples_csv(samples_text)

    if sheet_csv is None:
        loader = fetch_sheet or fetch_sample_template
        sheet_csv = loader(metadata_gid)
    template_rows = parse_sample_template(sheet_csv)
    sheet_wells = {row["well"] for row in template_rows}
    sheet_names = [(row["well"], row["sample_name"]) for row in template_rows]
    correlation = correlate_samples(sample_rows, sheet_wells, sheet_names)
    control_names = set(correlation.control_set)

    summary = RunSummary()
    for control in correlation.controls:
        summary.warnings.append(control_warning(control.sample, control.barcodes))

    listed = list_objects_with_size(
        s3_client,
        bucket,
        samples_dir_prefix(prefix),
        predicate=lambda key: is_scale_h5ad_key(key, prefix),
    )
    targets: list[ListedObject] = []
    for obj in listed:
        filename = obj.key.rsplit("/", 1)[-1]
        if is_control_file(filename, control_names):
            continue
        targets.append(obj)

    summary.total = len(targets)
    if not targets:
        return summary

    writer = TsvWriter(output_path, SCALE_H5AD_COLUMNS)
    max_workers = min(workers or 64, len(targets))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _fetch_crc, s3_client, bucket, obj.key, retries=retries
            ): obj
            for obj in targets
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(targets), desc="Processing")

        # Keep TSV rows in listing order, not completion order.
        results: dict[str, tuple[str | None, str]] = {}
        for fut in iterator:
            obj = futures[fut]
            results[obj.key] = fut.result()

    for obj in targets:
        filename = obj.key.rsplit("/", 1)[-1]
        crc, crc_err = results[obj.key]
        sample = sample_from_filename(filename)
        writer.append_row(
            [
                filename,
                s3_uri_for(bucket, obj.key),
                crc if crc is not None else "",
                sample,
                format_samples_column(
                    correlation.sample_names.get(sample, ()), lab_name
                ),
            ]
        )
        if not crc_err:
            summary.crc_ok += 1
        else:
            summary.failures.append((obj.key, crc_err, ""))

    return summary
