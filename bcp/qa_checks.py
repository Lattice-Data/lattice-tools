"""
QA validation logic extracted from the QA notebook (cells 4-8).
Pure functions: no S3, no file I/O. Callers pass in gathered data and get back
structured results (errors, missing/extra file lists, etc.).
"""

from __future__ import annotations

import logging
from typing import Any

from qa_mods import (
    cellranger_expected,
    parse_raw_filename,
    raw_expected,
    raw_optional,
)


VALID_PROBES = [
    "Chromium Mouse Transcriptome Probe Set v1.1.1",
    "Chromium Human Transcriptome Probe Set v1.1.0",
]


def validate_fastq_counts(
    fastq_log: dict[str, dict[str, list[str]]], raw_assay: str
) -> list[str]:
    """
    Check for mismatching number of files across modalities per sample.

    - scale: GEX vs hash_oligo must have equal counts when both present.
    - sci_plex: GEX vs hash_oligo when both present; GEX_hash_oligo only → no check.
    - sci_jumbo: No validation; logs that sci_jumbo is not validated.
    - 10x: GEX–CRI and GEX–ATAC pairs must have equal counts when present;
      GEX-only or ATAC-only → no check. Logs when CRI+ATAC present (future QA).
    - 10x_viral_ORF: No validation (legacy); logs that it is not validated.

    Returns list of error messages.
    """
    logger = logging.getLogger(__name__)
    errors: list[str] = []

    if raw_assay == "sci_jumbo":
        logger.warning("Not validating fastq counts for sci_jumbo at the moment.")
        return errors

    if raw_assay == "10x_viral_ORF":
        logger.warning(
            "Not validating fastq counts for 10x_viral_ORF (legacy/outlier)."
        )
        return errors

    for sample, v in fastq_log.items():
        if raw_assay == "scale":
            gex_list = v.get("GEX", [])
            hash_list = v.get("hash_oligo", [])
            if gex_list and hash_list and len(gex_list) != len(hash_list):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(hash_list)} hash_oligo"
                )
        elif raw_assay == "sci_plex":
            gex_list = v.get("GEX", [])
            hash_list = v.get("hash_oligo", [])
            if gex_list and hash_list and len(gex_list) != len(hash_list):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(hash_list)} hash_oligo"
                )
            # GEX_hash_oligo only or single modality → no comparison
        elif raw_assay == "10x":
            cri_list = v.get("CRI", [])
            atac_list = v.get("ATAC", [])
            gex_list = v.get("GEX", [])
            if cri_list and atac_list:
                logger.warning(
                    "Unexpected: CRI+ATAC present for sample %s; "
                    "QA expansion needed in the future.",
                    sample,
                )
            if gex_list and cri_list and len(gex_list) != len(cri_list):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(cri_list)} CRI"
                )
            if gex_list and atac_list and len(gex_list) != len(atac_list):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(atac_list)} ATAC"
                )
            # GEX-only or ATAC-only → no comparison

    return errors


def validate_read_metadata(
    read_metadata: dict[str, Any], raw_assay: str
) -> tuple[dict[str, dict[str, int]], list[str]]:
    """
    Build group read counts from metadata and check R1/R2 consistency.
    Returns (group_read_counts, errors).
    group_read_counts: {group: {assay: total_reads}}
    """
    errors: list[str] = []
    group_read_counts: dict[str, dict[str, int]] = {}

    for f, meta in read_metadata.items():
        if "_R2_" in f:
            continue
        if meta.get("errors"):
            errors.append(
                f"METADATA.JSON ERROR: {f} has error in metadata.json:{meta['errors']}"
            )
            continue
        reads = meta.get("read_count")
        if reads is None:
            continue
        parsed = parse_raw_filename(f, raw_assay)
        if parsed is not None:
            _run, group, assay, _ug, _barcode = parsed
            if group not in group_read_counts:
                group_read_counts[group] = {assay: reads}
            elif assay not in group_read_counts[group]:
                group_read_counts[group][assay] = reads
            else:
                group_read_counts[group][assay] += reads
        if "_R1_" in f:
            r2file = f.replace("_R1_", "_R2_")
            if r2file in read_metadata:
                r2meta = read_metadata[r2file]
                if r2meta.get("errors"):
                    errors.append(
                        f"METADATA.JSON ERROR: {r2file} has error in "
                        f"metadata.json:{r2meta['errors']}"
                    )
                    continue
                r2reads = r2meta.get("read_count")
                if reads != r2reads:
                    errors.append(f"READ COUNT ERROR:{f}-{reads},{r2file}-{r2reads}")

    return group_read_counts, errors


def check_expected_raw_files(
    all_raw_files: list[str], raw_assay: str
) -> tuple[int, list[dict[str, Any]], list[str]]:
    """
    Check that expected raw file endings are present for each "beginning".
    Returns (all_good_count, raw_lost list of dicts, raw_found list of paths).
    """
    beginnings: dict[str, dict[str, Any]] = {}
    for fullpath in all_raw_files:
        parsed = parse_raw_filename(fullpath, raw_assay)
        if parsed is None:
            continue
        run, group, assay, ug, barcode = parsed
        b = f"{run}-{group}_{assay}-{ug}-{barcode}"
        if b not in beginnings:
            raw_dir = "/".join(fullpath.split("/")[:-1])
            endings = list(raw_expected.get(raw_assay, []))
            if raw_assay == "10x_viral_ORF" and assay == "GEX":
                endings = list(raw_expected.get("10x", []))
            beginnings[b] = {"raw_dir": raw_dir, "endings": endings}

    all_good = 0
    raw_lost: list[dict[str, Any]] = []
    raw_found: list[str] = []
    raw_found_set = set(all_raw_files)

    for b, v in beginnings.items():
        temp_missing: dict[str, Any] = {"path": b}
        for e in v["endings"]:
            f = f"{v['raw_dir']}/{b}{e}"
            if f not in raw_found_set:
                if e.endswith("-metadata.json") and (
                    f.replace("-metadata.json", "") not in raw_found_set
                ):
                    continue
                temp_missing[e] = f
            else:
                raw_found.append(f)
        if len(temp_missing) > 1:
            raw_lost.append(temp_missing)
        else:
            all_good += 1

    return all_good, raw_lost, raw_found


def check_extra_raw_files(
    all_raw_files: list[str],
    raw_found: list[str],
    raw_assay: str,
) -> list[str]:
    """
    Identify raw files that are not in the expected set (optionals allowed).
    Returns list of extra file paths.
    """
    raw_found_set = set(raw_found)
    optional_endings = raw_optional.get(raw_assay, [])
    if raw_assay == "10x_viral_ORF":
        optional_endings = raw_optional.get("10x", [])

    extra: list[str] = []
    for f in all_raw_files:
        if f in raw_found_set:
            continue
        if f.endswith("-metadata.json") and (
            f.replace("-metadata.json", "") in set(all_raw_files)
        ):
            continue
        if (raw_assay in raw_optional or raw_assay == "10x_viral_ORF") and (
            optional_endings
        ):
            parsed = parse_raw_filename(f, raw_assay)
            if parsed is not None:
                run, group, assay, ug, barcode = parsed
                b = f"{run}-{group}_{assay}-{ug}-{barcode}"
                raw_dir = "/".join(f.split("/")[:-1])
                suffix = f.replace(f"{raw_dir}/{b}", "")
                endings = (
                    raw_optional["10x"]
                    if (raw_assay == "10x_viral_ORF" and assay == "GEX")
                    else optional_endings
                )
                if suffix in endings:
                    continue
        extra.append(f)
    return extra


def validate_processed_group(
    group_name: str,
    proc_files: list[str],
    report: dict[str, Any],
    group_read_counts: dict[str, dict[str, int]],
) -> dict[str, Any]:
    """
    Validate one group's processed (cellranger) outputs against report and expected files.
    Caller provides report from parse_web_summ + parse_met_summ.
    Returns dict with keys: errors, alerts, proc_missing (list of dicts), process_extra (list).
    """
    errors: list[str] = []
    alerts: list[dict[str, Any]] = []
    proc_missing: list[dict[str, Any]] = []
    process_extra: list[str] = []

    chem = report.get("chem")
    extra = report.get("extra", [])
    software = report.get("software", "")

    if software == "cellranger-9.0.1":
        sub = report.get("sub", "")
    elif software == "cellranger-10.0.0":
        sub = "multi"
    else:
        errors.append(
            f"CR ERROR: {group_name} version is {software} but should be 9.0.1 or 10.0.0"
        )
        sub = ""

    if sub:
        if "min-crispr-umi" in report:
            cri_umi = report["min-crispr-umi"]
            if cri_umi != "3":
                errors.append(
                    f"CR ERROR: {group_name} min-crispr-umi is {cri_umi} but should be 3"
                )

        if "incl_int" in report:
            intron = report["incl_int"]
            if intron != "true":
                errors.append(
                    f"CR ERROR: {group_name} include-introns is {intron} but should be true"
                )

        if "create-bam" in report:
            bam = report["create-bam"]
            if chem != "flex" and bam != "true":
                errors.append(
                    f"CR ERROR: {group_name} create-bam is {bam} but should be true"
                )

    for a in report.get("gex_alerts", []):
        alert = {"group": group_name, "modality": "GEX", **a}
        alerts.append(alert)
    for a in report.get("crispr_alerts", []):
        alert = {"group": group_name, "modality": "CRI", **a}
        alerts.append(alert)

    if sub == "multi" and software in cellranger_expected:
        if chem == "flex":
            expected = list(cellranger_expected[software]["flex"]["outs"])
            per_samp_expected = list(
                cellranger_expected[software]["flex"]["per_sample"]
            )
        else:
            expected = list(cellranger_expected[software]["nonflex"]["outs"])
            per_samp_expected = list(
                cellranger_expected[software]["nonflex"]["per_sample"]
            )
    elif sub == "count" and "count" in cellranger_expected:
        expected = list(cellranger_expected["count"]["outs"])
        per_samp_expected = []
    else:
        expected = []
        per_samp_expected = []

    if "CRISPR" in extra or "Antibody" in extra:
        if software != "cellranger-10.0.0":
            expected.append("multi/count/feature_reference.csv")
            per_samp_expected.append("count/feature_reference.csv")
        if "CRISPR" in extra and software != "cellranger-10.0.0":
            per_samp_expected.append("count/crispr_analysis.tar.gz")
        if "Antibody" in extra:
            per_samp_expected.append("count/antibody_analysis.tar.gz")
    if "CellAnnotate" in extra:
        per_samp_expected.append("count/cell_types.tar.gz")

    if sub == "multi" and report.get("multiplex") and software != "cellranger-10.0.0":
        expected.append("multi/multiplexing_analysis.tar.gz")

    actual = [
        f.split("/outs/", 1)[1]
        for f in proc_files
        if "/outs/" in f and f.split("/")[-1] != "curated.h5ad"
    ]
    if not actual:
        return {
            "errors": errors,
            "alerts": alerts,
            "proc_missing": proc_missing,
            "process_extra": process_extra,
        }

    missing = [f for f in expected if f not in actual]
    if missing:
        temp_missing: dict[str, Any] = {"group": group_name}
        for m in missing:
            temp_missing[m] = "Y"
        proc_missing.append(temp_missing)

    if per_samp_expected:
        # actual entries are paths under outs, e.g. per_sample_outs/SAMPLE/count/...
        samples = set(
            p.split("/")[1]
            for p in actual
            if p.startswith("per_sample_outs/") and len(p.split("/")) > 1
        )
        for s in samples:
            expected_samp = [f"per_sample_outs/{s}/{e}" for e in per_samp_expected]
            missing_samp = [e for e in expected_samp if e not in actual]
            if missing_samp:
                temp_missing = {"group": f"{group_name}/{s}"}
                for m in missing_samp:
                    key = m.replace(f"per_sample_outs/{s}/", "")
                    temp_missing[key] = "Y"
                proc_missing.append(temp_missing)
            expected.extend(expected_samp)

    extra_files = [
        f for f in actual if f not in expected and not f.endswith("manifest.json")
    ]
    if extra_files:
        process_extra.extend(extra_files)

    for k, v in report.items():
        if not k.endswith("_reads"):
            continue
        assay = k.replace("_reads", "")
        if group_name in group_read_counts:
            if assay in group_read_counts[group_name]:
                v2 = group_read_counts[group_name][assay]
                if v != v2:
                    errors.append(
                        f"READ COUNT ERROR: {group_name} {assay} {v} from proc,"
                        f"{v2} from raw,{v - v2} diff"
                    )
        else:
            pass  # WARNING only in notebook

    if "Probe Set Name" in report:
        if report["Probe Set Name"] not in VALID_PROBES:
            errors.append(
                f"WARNING: Invalid probe set for {group_name}: "
                f"{report['Probe Set Name']}"
            )

    return {
        "errors": errors,
        "alerts": alerts,
        "proc_missing": proc_missing,
        "process_extra": process_extra,
    }
