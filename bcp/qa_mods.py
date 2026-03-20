"""
QA parsing helpers. Constants live in qa_constants; re-exported here for backward compatibility.
"""

from __future__ import annotations

import json
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

import pandas as pd
from bs4 import BeautifulSoup

from qa_constants import (
    ALLOWED_RAW_ASSAYS,
    cellranger_expected,
    chemistries,
    raw_expected,
    raw_optional,
    valid_assays,
)

__all__ = [
    "ALLOWED_RAW_ASSAYS",
    "QARunContext",
    "cellranger_expected",
    "chemistries",
    "ingest_merged_trimmer_from_s3",
    "normalize_raw_assay",
    "raw_expected",
    "raw_optional",
    "resolve_qa_run_context",
    "valid_assays",
    "parse_met_summ",
    "parse_web_summ",
    "grab_trimmer_stats",
    "parse_raw_filename",
    "load_files_from_manifest",
    "extract_run_id_from_trimmer_filename",
    "extract_run_id_from_merged_trimmer_path",
    "grab_merged_trimmer_stats",
    "grab_merged_trimmer_q30",
]


def normalize_raw_assay(value: str | None) -> str:
    """
    Strip whitespace, fix common casing for 10x assays, and validate against ALLOWED_RAW_ASSAYS.

    Raises:
        ValueError: if missing or not a supported raw assay type.
    """
    if value is None or not str(value).strip():
        raise ValueError(
            "ERROR: raw_assay is not specified. "
            "Set it to one of: '10x', '10x_viral_ORF', 'sci_jumbo', 'sci_plex', 'scale'."
        )
    s = str(value).strip()
    lower = s.lower()
    if lower == "10x":
        s = "10x"
    elif lower == "10x_viral_orf":
        s = "10x_viral_ORF"
    if s not in ALLOWED_RAW_ASSAYS:
        raise ValueError(
            f"HUGE ERROR: raw_assay='{s}' is not recognized. "
            "Update the parameter to one of: "
            "'10x', '10x_viral_ORF', 'sci_jumbo', 'sci_plex', 'scale'."
        )
    return s


def _split_s3_uri(uri: str) -> tuple[str | None, str]:
    """Return (bucket, key) for s3:// URIs; else (None, uri as key)."""
    u = str(uri).strip()
    if not u.startswith("s3://"):
        return None, u
    rest = u[5:]
    if "/" not in rest:
        return rest, ""
    bucket, key = rest.split("/", 1)
    return bucket or None, key


def _manifest_buckets_from_column(
    manifest_path: str,
    delimiter: str,
    s3_column: int,
    has_header: bool,
) -> frozenset[str]:
    df = pd.read_csv(manifest_path, sep=delimiter, header=0 if has_header else None)
    buckets: set[str] = set()
    for cell in df.iloc[:, s3_column]:
        if cell is None or (isinstance(cell, float) and pd.isna(cell)):
            continue
        b, _ = _split_s3_uri(str(cell).strip())
        if b:
            buckets.add(b)
    return frozenset(buckets)


@dataclass(frozen=True)
class QARunContext:
    """Resolved QA run identity: S3 bucket, listing prefix, output file stem, normalized assay."""

    data_source: Literal["s3", "manifest"]
    raw_assay: str
    bucket: str
    provider: str
    proj: str
    order: str
    output_label: str
    listing_prefix: str
    manifest_path: str = ""
    manifest_delimiter: str = "\t"
    manifest_s3_column: int = 0
    manifest_has_header: bool = False


def resolve_qa_run_context(
    *,
    data_source: str,
    raw_assay: str,
    s3_path: str = "",
    provider: str = "",
    proj: str = "",
    order: str = "",
    run_label: str = "",
    manifest_path: str = "",
    manifest_delimiter: str = "\t",
    manifest_s3_column: int = 0,
    manifest_has_header: bool = False,
) -> QARunContext:
    """
    Validate notebook parameters and return a single context for gathering and output paths.

    * **s3** mode: require either ``s3://czi-*/proj/order`` or ``provider`` + ``proj`` + ``order``.
      Output files use ``run_label`` if set, else ``order``.
    * **manifest** mode: require ``manifest_path`` and non-empty ``run_label`` for output names.
      ``bucket`` is inferred from ``s3://czi-*`` URIs in the manifest column (single bucket).
    """
    ds = str(data_source).strip().lower()
    if ds not in ("s3", "manifest"):
        raise ValueError(
            f"Invalid data_source: {data_source!r}. Must be 'manifest' or 's3'."
        )

    assay = normalize_raw_assay(raw_assay)

    if ds == "manifest":
        mp = str(manifest_path).strip()
        if not mp:
            raise ValueError("manifest mode requires a non-empty manifest_path.")
        rl = str(run_label).strip()
        if not rl:
            raise ValueError(
                "manifest mode requires run_label (used for output files, e.g. *_errors.txt)."
            )
        buckets = _manifest_buckets_from_column(
            mp, manifest_delimiter, manifest_s3_column, manifest_has_header
        )
        if not buckets:
            raise ValueError(
                "No s3:// URIs with a bucket found in the manifest column; cannot infer bucket."
            )
        if len(buckets) > 1:
            raise ValueError(
                f"Manifest lists multiple S3 buckets {sorted(buckets)}; expected a single bucket."
            )
        bucket = next(iter(buckets))
        if not bucket.startswith("czi-"):
            raise ValueError(f"Invalid bucket in manifest (expected czi-*): {bucket!r}")
        provider_name = bucket[len("czi-") :]
        return QARunContext(
            data_source="manifest",
            raw_assay=assay,
            bucket=bucket,
            provider=provider_name,
            proj="",
            order="",
            output_label=rl,
            listing_prefix="",
            manifest_path=mp,
            manifest_delimiter=manifest_delimiter,
            manifest_s3_column=manifest_s3_column,
            manifest_has_header=manifest_has_header,
        )

    # --- s3 ---
    sp = str(s3_path).strip()
    prov = str(provider).strip()
    p = str(proj).strip()
    ord_seg = str(order).strip()

    if sp:
        u = urlparse(sp)
        if u.scheme != "s3" or not u.netloc:
            raise ValueError(f"Invalid s3_path (expected s3://...): {sp}")
        bucket = u.netloc
        if not bucket.startswith("czi-"):
            raise ValueError(f"Invalid bucket in s3_path (expected czi-*): {bucket}")
        provider_name = bucket[len("czi-") :]
        parts = [x for x in u.path.split("/") if x]
        if len(parts) < 2:
            raise ValueError(f"S3 path must include proj + order segments: {sp}")
        p = parts[-2]
        ord_seg = parts[-1]
    else:
        if not prov or not p or not ord_seg:
            raise ValueError(
                "s3 mode: set s3_path, or set all of provider, proj, and order."
            )
        bucket = f"czi-{prov}"
        provider_name = prov

    out = str(run_label).strip() or ord_seg
    listing = f"{p}/{ord_seg}/"
    return QARunContext(
        data_source="s3",
        raw_assay=assay,
        bucket=bucket,
        provider=provider_name,
        proj=p,
        order=ord_seg,
        output_label=out,
        listing_prefix=listing,
    )


def parse_met_summ(f):
    df = pd.read_csv(f)
    if len(df) == 1:
        report = {"GEX_reads": int(df["Number of Reads"].iloc[0].replace(",", ""))}

        return report

    lib_reads = df[
        (df["Metric Name"].isin(["Number of reads", "Number of short reads skipped"]))
        & (df["Grouped By"] == "Fastq ID")
    ]
    # Keep numeric values in a separate series to avoid assigning int into string column
    metric_value_int = lib_reads["Metric Value"].str.replace(",", "").astype(int)

    gex_mask = lib_reads["Library Type"] == "Gene Expression"
    report = {"GEX_reads": int(metric_value_int[gex_mask].sum())}

    if "CRISPR Guide Capture" in lib_reads["Library Type"].unique():
        cri_mask = lib_reads["Library Type"] == "CRISPR Guide Capture"
        report["CRI_reads"] = int(metric_value_int[cri_mask].sum())

    return report


def parse_web_summ(f):
    report = {"extra": [], "alerts": []}
    data = None
    with open(f) as html_doc:
        soup = BeautifulSoup(html_doc, "html.parser")
    for x in soup.find_all("script"):
        if not x.string:
            continue
        match = re.search("const data = ", x.string)
        if match:
            end = match.end()
            data = json.loads(x.string[end:])
            break

    if data is None:
        raise ValueError(
            f"No embedded Cell Ranger JSON ('const data = ...') found in {f!r}"
        )

    if "summary" in data:
        report["sub"] = data["summary"]["sample"]["subcommand"]
        summary_rows = data["summary"]["summary_tab"]["pipeline_info_table"]["rows"]
        gex_tab = {row[0].lower(): row[1] for row in summary_rows}
        for info in [
            row for row in summary_rows if row[0] in ("Transcriptome", "Probe Set Name")
        ]:
            report[info[0]] = info[1]
    else:
        gex_content = data["library"]["data"]["gex_tab"]["content"]
        summary_rows = gex_content["parameters_table"]["rows"]
        gex_tab = {row[0].lower(): row[1] for row in summary_rows}
        for info in [
            row for row in summary_rows if row[0] in ("Transcriptome", "Probe Set Name")
        ]:
            report[info[0]] = info[1]
        if "sample" in data:
            report["sub"] = data["sample"]["subcommand"]

    chem = gex_tab["chemistry"]
    report["chem"] = chemistries.get(chem, chem)

    if "library" in data:
        if (
            "crispr_tab" in data["library"]["data"]
            and data["library"]["data"]["crispr_tab"]
        ):
            report["extra"].append("CRISPR")
        if "antibody_tab" in data["library"]["data"].keys():
            if data["library"]["data"]["antibody_tab"]:
                report["extra"].append("Antibody")
    if "experimental_design" in data:
        for line in data["experimental_design"]["csv"].split("\n"):
            if line.startswith("["):
                cat = line
                if cat == "[samples]":
                    report["multiplex"] = True
            if "," in line:
                path = line.strip().split(",")
                if path[0] == "skip-cell-annotation" and path[1] == "false":
                    report["extra"].append("CellAnnotate")
                elif path[0] == "min-crispr-umi":
                    report["min-crispr-umi"] = path[1]
                elif path[0] == "create-bam":
                    report["create-bam"] = path[1]
                elif path[0] == "reference" and cat == "[gene-expression]":
                    report["ref"] = path[1]

    # Transcriptome from pipeline table is the canonical QA reference (overrides
    # experimental_design [gene-expression] reference when both exist).
    report["ref"] = gex_tab["transcriptome"]
    if chem != "Flex Gene Expression":
        report["incl_int"] = gex_tab["include introns"].lower()

    if "pipeline version" in gex_tab:
        report["software"] = gex_tab["pipeline version"]
    elif "pipeline_version" in data:
        report["software"] = "cellranger-" + data["pipeline_version"]
    if "library" in data:
        if "software" not in report:
            header_info = data["library"]["data"]["header_info"]
            report["software"] = header_info["Pipeline Version"]
        report["gex_alerts"] = data["library"]["data"]["gex_tab"]["alerts"]
        if (
            "crispr_tab" in data["library"]["data"]
            and data["library"]["data"]["crispr_tab"]
        ):
            report["crispr_alerts"] = data["library"]["data"]["crispr_tab"]["alerts"]

    if "summary" in data:
        report["gex_alerts"] = data["summary"]["alarms"]["alarms"]

    return report


def grab_trimmer_stats(
    trimmer_failure_stats: dict,
    exp: str,
    csv_path: str | Path,
) -> None:
    """
    Parse a trimmer-failure_codes CSV and update trimmer_failure_stats in place.
    Caller is responsible for downloading the file from S3 (if needed) and cleanup.
    exp is typically "/".join(s3_key.split("/")[1:3]) for the experiment identifier.
    """
    if exp not in trimmer_failure_stats:
        trimmer_failure_stats[exp] = {"rsq": [], "trimmer_fail": []}
    stats_df = pd.read_csv(csv_path)
    stats_df.columns = stats_df.columns.str.replace(" ", "_")
    trimmer_fail = 0
    total_reads = 0
    for row in stats_df.itertuples():
        total_reads = row.total_read_count
        if row.reason == "rsq file":
            trimmer_failure_stats[exp]["rsq"].append(
                100 * row.failed_read_count / total_reads
            )
        else:
            trimmer_fail += row.failed_read_count
    if total_reads > 0:
        trimmer_fail_pct = trimmer_fail / total_reads
        trimmer_failure_stats[exp]["trimmer_fail"].append(100 * trimmer_fail_pct)


def extract_run_id_from_trimmer_filename(filename: str) -> str | None:
    """
    Extract the RunID / wafer identifier from a trimmer statistics filename.

    Filenames follow the same leading convention as raw FASTQs/CRAMs:
    `{RunID}-{GroupID}_{Assay}-..._trimmer-failure_codes.csv`
    or `{RunID}-..._trimmer-stats.csv`.

    RunID is usually 6 digits but may grow to up to 8 digits. We therefore
    accept any numeric prefix of length 6–8 characters as a valid RunID.

    Args:
        filename: Basename or full path to the trimmer CSV.

    Returns:
        The RunID string if it can be parsed, otherwise None.
    """
    name = filename.split("/")[-1]
    parts = name.split("-")
    if not parts:
        return None
    run_id = parts[0]
    if run_id.isdigit() and 6 <= len(run_id) <= 8:
        return run_id
    return None


def extract_run_id_from_merged_trimmer_path(s3_key: str) -> str | None:
    """
    Extract RunID (wafer id) from an S3 key pointing to a merged trimmer file.

    - Scale/sci: path like .../raw/438761/merged_trimmer-failure_codes.csv
      → RunID = path parent (438761). Parent segment must be 6–8 digit numeric.
    - 10x: path under order folder (no /raw/), e.g. .../order/438761_merged_...
      → RunID = 6–8 digit prefix from filename.

    Args:
        s3_key: S3 object key (path without bucket).

    Returns:
        RunID string if determinable, else None.
    """
    parts = s3_key.split("/")
    name = parts[-1] if parts else ""

    # Must look like a merged trimmer file
    if (
        "merged_trimmer-failure_codes.csv" not in name
        and "merged_trimmer-stats.csv" not in name
    ):
        return None

    # 10x: merged files under {order}/ — no "raw" in path; filename has optional run_id prefix
    if "/raw/" not in s3_key:
        # Filename: 438761_merged_trimmer-failure_codes.csv or merged_trimmer-failure_codes.csv
        if name.startswith("merged_"):
            return None  # Single-wafer, no prefix
        prefix = name.split("_merged_")[0]
        if prefix.isdigit() and 6 <= len(prefix) <= 8:
            return prefix
        return None

    # Scale/sci: .../raw/{run_id}/merged_trimmer-*.csv — run_id = parent dir
    if len(parts) < 2:
        return None
    run_id = parts[-2]
    if run_id.isdigit() and 6 <= len(run_id) <= 8:
        return run_id
    return None


def grab_merged_trimmer_stats(csv_path: str | Path) -> dict | None:
    """
    Parse a merged trimmer-failure_codes CSV and return TT and RSQ pass/fail stats.

    TT (read_group == "TT"):
      - tt_total_reads: total_read_count from first TT row.
      - tt_failed_reads: sum of failed_read_count across ALL TT rows (all reasons).
      - tt_pass_reads, tt_fail_pct, tt_pass_pct derived from above.

    RSQ (reason == "rsq file", all read groups):
      - rsq_total_reads: sum of total_read_count across all rows with reason "rsq file"
        (each read group has its own total).
      - rsq_failed_reads: sum of failed_read_count across all rows with reason "rsq file".
      - rsq_pass_reads, rsq_fail_pct, rsq_pass_pct derived from above.

    Returns:
        Dict with tt and rsq metrics, or None if no TT row.
    """
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.replace(" ", "_")
    if "read_group" not in df.columns or "reason" not in df.columns:
        return None
    tt = df[df["read_group"] == "TT"]
    if tt.empty:
        return None

    # TT: all rows with read_group == "TT", all reasons
    tt_total = int(tt["total_read_count"].iloc[0])
    if tt_total <= 0:
        return None
    tt_failed = int(tt["failed_read_count"].sum())
    tt_pass = tt_total - tt_failed

    # RSQ: all rows with reason == "rsq file" across all read groups
    rsq_rows = df[df["reason"] == "rsq file"]
    rsq_total = int(rsq_rows["total_read_count"].sum()) if not rsq_rows.empty else 0
    rsq_failed = int(rsq_rows["failed_read_count"].sum()) if not rsq_rows.empty else 0
    rsq_pass = rsq_total - rsq_failed if rsq_total > 0 else 0

    return {
        "tt_total_reads": tt_total,
        "tt_failed_reads": tt_failed,
        "tt_pass_reads": tt_pass,
        "tt_fail_pct": 100.0 * tt_failed / tt_total,
        "tt_pass_pct": 100.0 * tt_pass / tt_total,
        "rsq_total_reads": rsq_total,
        "rsq_failed_reads": rsq_failed,
        "rsq_pass_reads": rsq_pass,
        "rsq_fail_pct": 100.0 * rsq_failed / rsq_total if rsq_total > 0 else None,
        "rsq_pass_pct": 100.0 * rsq_pass / rsq_total if rsq_total > 0 else None,
    }


def grab_merged_trimmer_q30(csv_path: str | Path) -> float | None:
    """
    Parse a merged trimmer-stats CSV and return TT-level Q30 percentage.

    Uses the "low quality bases" segment for read_group "TT":
    Q30 % = 100 * num_matched_bases / (num_matched_bases + num_failures).

    Args:
        csv_path: Path to the merged trimmer-stats CSV.

    Returns:
        Q30 percentage (0–100) or None if TT row or segment missing.
    """
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.replace(" ", "_")
    for col in ("read_group", "segment_label", "num_matched_bases", "num_failures"):
        if col not in df.columns:
            return None
    tt = df[df["read_group"] == "TT"]
    if tt.empty:
        return None
    low_qual = tt[tt["segment_label"] == "low quality bases"]
    if low_qual.empty:
        return None
    row = low_qual.iloc[0]
    matched = row["num_matched_bases"]
    failures = row["num_failures"]
    total = matched + failures
    if total <= 0:
        return None
    return 100.0 * matched / total


def ingest_merged_trimmer_from_s3(
    bucket: str,
    key: str,
    merged_wafer_stats: dict[str, Any],
    s3_client: Any,
) -> None:
    """
    Download a merged trimmer CSV under ``key``, parse TT/RSQ or Q30, and update ``merged_wafer_stats``.
    No-op if the key is not a merged trimmer file or RunID cannot be determined.
    """
    run_id = extract_run_id_from_merged_trimmer_path(key)
    if not run_id:
        return
    name = key.split("/")[-1]
    if (
        "merged_trimmer-failure_codes.csv" not in name
        and "merged_trimmer-stats.csv" not in name
    ):
        return

    with tempfile.NamedTemporaryFile(mode="w+b", delete=False, suffix=".csv") as tf:
        local = tf.name
    try:
        s3_client.download_file(bucket, key, local)
        if "merged_trimmer-failure_codes.csv" in name:
            merged = grab_merged_trimmer_stats(local)
            if merged:
                merged_wafer_stats[run_id] = merged_wafer_stats.get(run_id, {})
                merged_wafer_stats[run_id].update(merged)
        else:
            q30 = grab_merged_trimmer_q30(local)
            if q30 is not None:
                merged_wafer_stats.setdefault(run_id, {})
                merged_wafer_stats[run_id]["tt_q30_pct"] = q30
    finally:
        Path(local).unlink(missing_ok=True)


def parse_raw_filename(f, raw_assay):
    """
    For scale data, use regex for determining assay, and "group" is replaced by "experiment", and there is no "barcode". This is because
    for scale data, all cram files within a single experiment will be used as input for a single run of SeqSuite.
    Otherwise, parse filename, split by '-'
    """
    filename = f.split("/")[-1]
    path = filename.split("-")
    if len(path) < 3:
        return None
    run = path[0]

    if raw_assay == "scale":
        group = f.split("/")[2]
        if re.search("GEX_hash_oligo", filename):
            assay = "GEX_hash_oligo"
        elif re.search("hash_oligo", filename):
            assay = "hash_oligo"
        elif re.search("GEX", filename):
            assay = "GEX"
        ug = filename.split("_")[-2]
        barcode = None

    else:
        group_assay = path[1]
        if group_assay.endswith("GEX_hash_oligo"):
            assay = "GEX_hash_oligo"
        else:
            match = False
            for v_a in valid_assays:
                if group_assay.endswith(v_a):
                    assay = v_a
                    match = True
            if not match:
                assay = group_assay.split("_")[-1]
        group = group_assay.replace(f"_{assay}", "")
        ug = path[2]
        barcode = path[3].split("_")[0].split(".")[0]

    return run, group, assay, ug, barcode


def load_files_from_manifest(
    manifest_path: str, delimiter: str, s3_column: int, has_header: bool = False
):
    """
    Load S3 file paths from a CSV/TSV manifest file.

    Args:
        manifest_path: Path to the manifest file
        delimiter: Field delimiter (',' for CSV, '\t' for TSV)
        s3_column: Column index (0-based) containing S3 URIs
        has_header: Whether the file has a header row to skip

    Returns:
        Tuple of:
        - all_raw_files: List of raw file S3 keys (paths containing '/raw/')
        - all_proc_files: Dict of {group: [processed file keys]} (paths containing '/processed/')
        - manifest_bucket: Bucket from ``s3://`` URIs, or None if none found
    """
    df = pd.read_csv(manifest_path, sep=delimiter, header=0 if has_header else None)
    s3_uris = df.iloc[:, s3_column].tolist()

    all_raw_files = []
    all_proc_files = {}
    buckets_seen: set[str] = set()

    for uri in s3_uris:
        if uri is None or (isinstance(uri, float) and pd.isna(uri)):
            continue
        uri = str(uri).strip()
        # Strip 's3://bucket-name/' prefix to get just the key
        if uri.startswith("s3://"):
            b, key = _split_s3_uri(uri)
            if b:
                buckets_seen.add(b)
            if not key:
                continue
        else:
            key = uri  # Assume it's already a key

        # Separate into raw vs processed
        if "/raw/" in key:
            all_raw_files.append(key)
        elif "/processed/" in key:
            # Extract group name from path: .../GROUP/processed/...
            path_before_processed = key.split("/processed/")[0]
            group = path_before_processed.split("/")[-1]
            if group not in all_proc_files:
                all_proc_files[group] = []
            all_proc_files[group].append(key)

    if len(buckets_seen) > 1:
        raise ValueError(
            f"Manifest mixes S3 buckets {sorted(buckets_seen)}; expected a single bucket."
        )
    manifest_bucket = next(iter(buckets_seen)) if buckets_seen else None

    return all_raw_files, all_proc_files, manifest_bucket
