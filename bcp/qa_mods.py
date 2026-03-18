"""
QA parsing helpers. Constants live in qa_constants; re-exported here for backward compatibility.
"""

import json
import re
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup

from qa_constants import (
    cellranger_expected,
    chemistries,
    raw_expected,
    raw_optional,
    valid_assays,
)

__all__ = [
    "cellranger_expected",
    "chemistries",
    "raw_expected",
    "raw_optional",
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
    with open(f) as html_doc:
        soup = BeautifulSoup(html_doc, "html.parser")
    for x in soup.find_all("script"):
        match = re.search("const data = ", x.string)
        if match:
            end = match.end()
            data = json.loads(x.string[end:])

    if "summary" in data:
        report["sub"] = data["summary"]["sample"]["subcommand"]
        summary_rows = data["summary"]["summary_tab"]["pipeline_info_table"]["rows"]
        gex_tab = {row[0].lower(): row[1] for row in summary_rows}
        for info in [
            row for row in summary_rows if row[0] in ("Transcriptome", "Probe Set Name")
        ]:
            report[info[0]] = info[1]
    else:
        summary_rows = data["library"]["data"]["gex_tab"]["content"][
            "parameters_table"
        ]["rows"]
        gex_tab = {row[0].lower(): row[1] for row in summary_rows}
        for info in [
            row for row in summary_rows if row[0] in ("Transcriptome", "Probe Set Name")
        ]:
            report[info[0]] = info[1]
        if "sample" in data:
            report["sub"] = data["sample"]["subcommand"]

    chem = gex_tab["chemistry"]
    report["chem"] = chemistries.get(chem, chem)

    report["extra"] = []
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

    # location of some additional info to QA
    report["ref"] = gex_tab["transcriptome"]
    if chem != "Flex Gene Expression":
        report["incl_int"] = gex_tab["include introns"].lower()

    if "pipeline version" in gex_tab:
        report["software"] = gex_tab["pipeline version"]
    elif "pipeline_version" in data:
        report["software"] = "cellranger-" + data["pipeline_version"]
    if "library" in data:
        if "software" not in report:
            report["software"] = data["library"]["data"]["header_info"][
                "Pipeline Version"
            ]
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


# TT failure reasons used for "TT failing" (failed TT reads) in merged file.
TT_FAIL_REASONS = ("rsq file", "sequence was too short", "sequence was too long")


def grab_merged_trimmer_stats(csv_path: str | Path) -> dict | None:
    """
    Parse a merged trimmer-failure_codes CSV and return wafer-level RSQ and TT pass/fail
    percentages and read counts.

    - Wafer total: total_read_count from the first TT row (not summed across rows).
    - RSQ failed: sum of failed_read_count over ALL rows where reason == "rsq file".
    - TT failed: sum of failed_read_count over TT rows only where reason is one of
      "rsq file", "sequence was too short", "sequence was too long".

    Returns:
        Dict with wafer_total_reads, rsq_failed_reads, rsq_fail_pct, rsq_pass_pct,
        rsq_pass_reads, tt_failed_reads, tt_fail_pct, tt_pass_pct, tt_pass_reads,
        or None if no TT row.
    """
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.replace(" ", "_")
    if "read_group" not in df.columns or "reason" not in df.columns:
        return None
    tt = df[df["read_group"] == "TT"]
    if tt.empty:
        return None
    # Wafer total from first TT row (total is per read group, do not sum across rows)
    wafer_total = int(tt["total_read_count"].iloc[0])
    if wafer_total <= 0:
        return None

    # RSQ: sum failed over ALL rows (all read groups) with reason "rsq file"
    rsq_failed = int(df.loc[df["reason"] == "rsq file", "failed_read_count"].sum())
    rsq_fail_pct = 100.0 * rsq_failed / wafer_total
    rsq_pass_pct = 100.0 * (wafer_total - rsq_failed) / wafer_total
    rsq_pass_reads = wafer_total - rsq_failed

    # TT failed: TT rows only, reasons rsq file + sequence too short/long
    tt_mask = (df["read_group"] == "TT") & (df["reason"].isin(TT_FAIL_REASONS))
    tt_failed = int(df.loc[tt_mask, "failed_read_count"].sum())
    tt_fail_pct = 100.0 * tt_failed / wafer_total
    tt_pass_pct = 100.0 * (wafer_total - tt_failed) / wafer_total
    tt_pass_reads = wafer_total - tt_failed

    return {
        "wafer_total_reads": wafer_total,
        "rsq_failed_reads": rsq_failed,
        "rsq_fail_pct": rsq_fail_pct,
        "rsq_pass_pct": rsq_pass_pct,
        "rsq_pass_reads": rsq_pass_reads,
        "tt_failed_reads": tt_failed,
        "tt_fail_pct": tt_fail_pct,
        "tt_pass_pct": tt_pass_pct,
        "tt_pass_reads": tt_pass_reads,
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
    """
    df = pd.read_csv(manifest_path, sep=delimiter, header=0 if has_header else None)
    s3_uris = df.iloc[:, s3_column].tolist()

    all_raw_files = []
    all_proc_files = {}

    for uri in s3_uris:
        # Strip 's3://bucket-name/' prefix to get just the key
        if uri.startswith("s3://"):
            # s3://bucket-name/path/to/file -> path/to/file
            parts = uri[5:].split("/", 1)  # Remove 's3://', split on first '/'
            if len(parts) > 1:
                key = parts[1]  # The key (path after bucket)
            else:
                continue  # Invalid URI, skip
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

    return all_raw_files, all_proc_files
