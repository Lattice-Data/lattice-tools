"""
QA validation logic extracted from the QA notebook (cells 4-8).
Pure functions: no S3, no file I/O. Callers pass in gathered data and get back
structured results (errors, missing/extra file lists, etc.).
"""

from __future__ import annotations

import logging
from typing import Any

from qa_constants import (
    SCALE_SAMPLES_FORBIDDEN_COLUMNS,
    SCALE_WORKFLOW_REQUIRED_PARAMS,
)
from qa_mods import (
    cellranger_expected,
    extract_read_indicator,
    make_read_partner,
    parse_raw_filename,
    raw_expected,
    raw_optional,
)


VALID_PROBES = [
    "Chromium Mouse Transcriptome Probe Set v1.1.1",
    "Chromium Human Transcriptome Probe Set v1.1.0",
]


def _fastq_count_mode(raw_assay: str) -> str:
    """Internal: assay policy for fastq count validation and summaries."""
    if raw_assay in ("sci_jumbo", "10x_viral_ORF"):
        return "skip"
    if raw_assay in ("scale", "sci_plex"):
        return "gex_hash"
    if raw_assay == "10x":
        return "10x"
    return "unknown"


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

    mode = _fastq_count_mode(raw_assay)
    if mode == "skip":
        if raw_assay == "sci_jumbo":
            logger.warning("Not validating fastq counts for sci_jumbo at the moment.")
        else:
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
            # Only compare when both modalities have at least one file (GEX-only or CRI-only → no check)
            if (
                len(gex_list) > 0
                and len(cri_list) > 0
                and len(gex_list) != len(cri_list)
            ):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(cri_list)} CRI"
                )
            if (
                len(gex_list) > 0
                and len(atac_list) > 0
                and len(gex_list) != len(atac_list)
            ):
                errors.append(
                    f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                    f"{len(atac_list)} ATAC"
                )

    return errors


def summarize_fastq_count_validation(
    fastq_log: dict[str, dict[str, list[str]]],
    raw_assay: str,
    errors: list[str],
) -> str:
    """
    Create a short, single-line report for fastq count validation.

    This is intentionally summary-only (no per-file lists). It reports both:
    - whether any mismatches were found (via `errors`)
    - whether comparisons were actually applicable (via `fastq_log`)
    """

    mismatches = len(errors)
    mode = _fastq_count_mode(raw_assay)

    if mode == "skip":
        return (
            f"Fastq count validation ({raw_assay}): not validated by design; "
            f"mismatches: {mismatches}. No comparisons performed."
        )

    if not fastq_log:
        return (
            f"Fastq count validation ({raw_assay}): checked 0; mismatches: "
            f"{mismatches}. No fastq_log data (nothing to compare)."
        )

    if raw_assay in ("scale", "sci_plex"):
        checked = sum(
            1
            for _sample, v in fastq_log.items()
            if v.get("GEX") and v.get("hash_oligo")
        )
        gex_only = sum(
            1
            for _sample, v in fastq_log.items()
            if v.get("GEX") and not v.get("hash_oligo")
        )
        hash_only = sum(
            1
            for _sample, v in fastq_log.items()
            if v.get("hash_oligo") and not v.get("GEX")
        )

        if checked > 0 and mismatches == 0:
            result = "All matched."
        elif checked == 0 and mismatches == 0:
            if gex_only > 0 and hash_only == 0:
                result = "Only GEX present (missing hash_oligo)."
            elif hash_only > 0 and gex_only == 0:
                result = "Only hash_oligo present (missing GEX)."
            else:
                result = "No comparable pairs (missing modality pairs)."
        elif mismatches > 0:
            result = (
                "Mismatches found (see preceding MISMATCH FQ COUNTS lines in stdout "
                "and the errors log)."
            )
        else:
            result = "No comparisons performed."

        return (
            f"Fastq count validation ({raw_assay}): checked {checked} "
            f"(GEX<->hash_oligo); mismatches: {mismatches}. {result}"
        )

    if raw_assay == "10x":
        checked_gex_cri = sum(
            1 for _sample, v in fastq_log.items() if v.get("GEX") and v.get("CRI")
        )
        checked_gex_atac = sum(
            1 for _sample, v in fastq_log.items() if v.get("GEX") and v.get("ATAC")
        )
        checked_total = checked_gex_cri + checked_gex_atac

        if checked_total > 0 and mismatches == 0:
            result = "All matched."
        elif checked_total == 0 and mismatches == 0:
            gex_present = sum(1 for _sample, v in fastq_log.items() if v.get("GEX"))
            cri_or_atac_present = sum(
                1 for _sample, v in fastq_log.items() if v.get("CRI") or v.get("ATAC")
            )
            if gex_present > 0 and cri_or_atac_present == 0:
                result = "Only GEX present (missing CRI/ATAC)."
            elif gex_present == 0 and cri_or_atac_present > 0:
                result = "Missing GEX (CRI/ATAC present)."
            elif gex_present > 0 and cri_or_atac_present > 0:
                result = "No comparable pairs (GEX groups differ from CRI/ATAC groups)."
            else:
                result = "No comparable pairs (missing modality pairs)."
        elif mismatches > 0:
            result = (
                "Mismatches found (see preceding MISMATCH FQ COUNTS lines in stdout "
                "and the errors log)."
            )
        else:
            result = "No comparisons performed."

        return (
            "Fastq count validation (10x): checked "
            f"{checked_total} (GEX<->CRI: {checked_gex_cri}, GEX<->ATAC: {checked_gex_atac}); "
            f"mismatches: {mismatches}. {result}"
        )

    return (
        f"Fastq count validation ({raw_assay}): checked 0; mismatches: "
        f"{mismatches}. Unsupported assay."
    )


def validate_read_metadata(
    read_metadata: dict[str, Any],
    raw_assay: str,
    *,
    print_success: bool = False,
    success_print_limit: int = 5,
    pairing_paths_print_limit: int = 100,
) -> tuple[dict[str, dict[str, int]], list[str], dict[str, list[str]]]:
    """
    Build group read counts from metadata and check R1/R2 consistency.

    Returns ``(group_read_counts, errors, pairing)`` where ``pairing`` has:

    - ``r1_without_r2_metadata``: R1 paths with ``read_count`` but no R2 key
    - ``r2_without_r1_metadata``: R2 paths with no corresponding R1 key

    Those issues are also appended to ``errors`` with a ``READ METADATA PAIRING:`` prefix
    so notebook callers can log them like other messages.

    R1/R2 classification uses ``extract_read_indicator()`` which matches the
    Illumina read indicator at the *end* of the filename (e.g. ``_R1_001.fastq.gz``),
    so ``_R1_`` / ``_R2_`` embedded in a group ID (e.g. ``q_pcf_R2``) is ignored.
    """
    errors: list[str] = []
    group_read_counts: dict[str, dict[str, int]] = {}

    # Optional instrumentation for stdout-only reporting.
    matched_examples: list[str] = []
    compared_pairs = 0
    matched_pairs = 0
    mismatched_pairs = 0
    skipped_no_r2 = 0
    r2_metadata_error = 0
    r1_without_r2_metadata: list[str] = []

    for f, meta in read_metadata.items():
        if extract_read_indicator(f) == "R2":
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
        if extract_read_indicator(f) == "R1":
            r2file = make_read_partner(f, "R1", "R2")
            if r2file in read_metadata:
                r2meta = read_metadata[r2file]
                if r2meta.get("errors"):
                    r2_metadata_error += 1
                    errors.append(
                        f"METADATA.JSON ERROR: {r2file} has error in "
                        f"metadata.json:{r2meta['errors']}"
                    )
                    continue
                r2reads = r2meta.get("read_count")
                compared_pairs += 1
                if reads != r2reads:
                    errors.append(f"READ COUNT ERROR:{f}-{reads},{r2file}-{r2reads}")
                    mismatched_pairs += 1
                else:
                    matched_pairs += 1
                    if print_success and len(matched_examples) < success_print_limit:
                        matched_examples.append(
                            f"MATCH: {f} ({reads}) == {r2file} ({r2reads})"
                        )
            else:
                skipped_no_r2 += 1
                r1_without_r2_metadata.append(f)

    r2_without_r1_metadata: list[str] = []
    for f in read_metadata:
        if extract_read_indicator(f) != "R2":
            continue
        r1file = make_read_partner(f, "R2", "R1")
        if r1file not in read_metadata:
            r2_without_r1_metadata.append(f)

    for path in r1_without_r2_metadata:
        errors.append(f"READ METADATA PAIRING: R1 present, no R2 metadata key: {path}")
    for path in r2_without_r1_metadata:
        errors.append(f"READ METADATA PAIRING: R2 present, no R1 metadata key: {path}")

    pairing: dict[str, list[str]] = {
        "r1_without_r2_metadata": r1_without_r2_metadata,
        "r2_without_r1_metadata": r2_without_r1_metadata,
    }

    if print_success:
        print(
            f"validate_read_metadata({raw_assay}): "
            f"r1_r2_pairs_compared={compared_pairs} "
            f"(matched={matched_pairs}, mismatched={mismatched_pairs}); "
            f"r1_missing_r2_metadata={skipped_no_r2}, "
            f"r2_missing_r1_metadata={len(r2_without_r1_metadata)}, "
            f"r2_metadata_errors={r2_metadata_error}"
        )
        for line in matched_examples:
            print(line)
        _print_pairing_path_lists(
            r1_without_r2_metadata,
            r2_without_r1_metadata,
            pairing_paths_print_limit,
        )

    return group_read_counts, errors, pairing


def _print_pairing_path_lists(
    r1_missing_r2: list[str],
    r2_missing_r1: list[str],
    limit: int,
) -> None:
    """Print paths for R1/R2 orphan metadata keys (truncated per ``limit`` per side)."""
    if r1_missing_r2:
        print(f"R1 paths with no R2 metadata entry ({len(r1_missing_r2)}):")
        for p in r1_missing_r2[:limit]:
            print(f"  {p}")
        if len(r1_missing_r2) > limit:
            print(f"  ... and {len(r1_missing_r2) - limit} more")
    if r2_missing_r1:
        print(f"R2 paths with no R1 metadata entry ({len(r2_missing_r1)}):")
        for p in r2_missing_r1[:limit]:
            print(f"  {p}")
        if len(r2_missing_r1) > limit:
            print(f"  ... and {len(r2_missing_r1) - limit} more")


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
    all_raw_set = set(all_raw_files)
    optional_endings = raw_optional.get(raw_assay, [])
    if raw_assay == "10x_viral_ORF":
        optional_endings = raw_optional.get("10x", [])

    extra: list[str] = []
    for f in all_raw_files:
        if f in raw_found_set:
            continue
        if f.endswith("-metadata.json") and (
            f.replace("-metadata.json", "") in all_raw_set
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

    actual_set = set(actual)

    # Cell Ranger v7+ may output either a BAI or a CSI index depending on
    # reference genome size. We require exactly one of the two files.
    bai_filename = "possorted_genome_bam.bam.bai"
    csi_filename = "possorted_genome_bam.bam.csi"
    bai_present = bai_filename in actual_set
    csi_present = csi_filename in actual_set
    bai_csi_case: str | None = None

    if bai_filename in expected:
        if bai_present and not csi_present:
            bai_csi_case = "bai_only"
        elif csi_present and not bai_present:
            bai_csi_case = "csi_only"
        elif not bai_present and not csi_present:
            bai_csi_case = "neither"
        else:
            bai_csi_case = "both"

    missing = [f for f in expected if f not in actual_set]

    # If Cell Ranger produced CSI instead of the expected BAI, don't report
    # BAI as missing (but still don't treat CSI as an "extra").
    if bai_csi_case == "csi_only":
        missing = [f for f in missing if f != bai_filename]

    if bai_csi_case == "neither":
        errors.append(
            f"CR ERROR: {group_name} missing both {bai_filename} and "
            f"{csi_filename}; expected exactly one"
        )

    if bai_csi_case == "both":
        errors.append(
            f"CR ERROR: {group_name} has both {bai_filename} and {csi_filename}; "
            "expected exactly one"
        )

    if missing:
        temp_missing: dict[str, Any] = {"group": group_name}
        for m in missing:
            temp_missing[m] = "Y"

        # Surface the BAI/CSI protocol failure in the notebook table too.
        if bai_csi_case == "neither":
            temp_missing[csi_filename] = "Y"
        elif bai_csi_case == "both":
            temp_missing[bai_filename] = "Y"
            temp_missing[csi_filename] = "Y"

        proc_missing.append(temp_missing)
    elif bai_csi_case == "both":
        # No other missing outputs, but we still want the notebook to show
        # BAI/CSI conflict.
        proc_missing.append({"group": group_name, bai_filename: "Y", csi_filename: "Y"})

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

    expected_set_for_extra = set(expected)
    # Suppress CSI from being flagged as an extra when BAI was expected.
    if bai_filename in expected and csi_present:
        expected_set_for_extra.add(csi_filename)

    extra_files = [
        f
        for f in actual
        if f not in expected_set_for_extra and not f.endswith("manifest.json")
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


def build_wafer_failure_stats(
    trimmer_failure_stats: dict[str, dict[str, list[float]]],
    exp_to_run_map: dict[str, str],
) -> dict[str, dict[str, list[float]]]:
    """
    Aggregate experiment-level trimmer failure statistics into wafer-level stats.

    Existing logic accumulates trimmer_failure_stats keyed by an experiment
    identifier (typically "/".join(s3_key.split("/")[1:3]) in the notebook).
    For wafer-level QA, we instead want to summarize by RunID (wafer id).

    This helper is pure: callers provide the mapping from experiment id to
    RunID, which is constructed at data-gathering time (e.g. when iterating
    over S3 keys).

    Args:
        trimmer_failure_stats:
            Dict mapping experiment id -> {"rsq": [...], "trimmer_fail": [...]}.
        exp_to_run_map:
            Dict mapping the same experiment ids to RunID strings.

    Returns:
        Dict mapping RunID (wafer) -> {"rsq": [...], "trimmer_fail": [...]}.
    """
    wafer_failure_stats: dict[str, dict[str, list[float]]] = {}

    for exp, stats in trimmer_failure_stats.items():
        run_id = exp_to_run_map.get(exp)
        if not run_id:
            continue
        if run_id not in wafer_failure_stats:
            wafer_failure_stats[run_id] = {"rsq": [], "trimmer_fail": []}

        rsq_vals = stats.get("rsq", [])
        trim_vals = stats.get("trimmer_fail", [])

        wafer_failure_stats[run_id]["rsq"].extend(rsq_vals)
        wafer_failure_stats[run_id]["trimmer_fail"].extend(trim_vals)

    return wafer_failure_stats


def validate_scale_workflow_info(params: dict) -> tuple[list[str], list[str]]:
    """
    Validate Scale workflow parameters extracted by ``parse_scale_workflow_info``.

    Returns ``(errors, info_messages)``.  Each required parameter in
    ``SCALE_WORKFLOW_REQUIRED_PARAMS`` is checked against the parsed dict;
    genome value is reported as an informational message.
    """
    errors: list[str] = []
    info_messages: list[str] = []

    for key, expected in SCALE_WORKFLOW_REQUIRED_PARAMS.items():
        actual = params.get(key)
        if actual != expected:
            errors.append(
                f"SCALE WORKFLOW ERROR: {key} is '{actual}' but should be '{expected}'"
            )

    info_messages.append(f"SCALE GENOME: {params.get('genome', 'N/A')}")

    return errors, info_messages


def validate_scale_samples_csv(columns: list[str]) -> list[str]:
    """
    Validate that a Scale ``samples.csv`` does not contain forbidden columns.

    Returns list of errors (one per forbidden column found).
    """
    errors: list[str] = []
    for col in SCALE_SAMPLES_FORBIDDEN_COLUMNS:
        if col in columns:
            errors.append(
                f"SCALE SAMPLES ERROR: forbidden column '{col}' is present in samples.csv"
            )
    return errors


def validate_scale_cb_tag(read_metadata: dict[str, Any]) -> list[str]:
    """
    Validate that all Scale CRAM files have ``cb_tag=True`` in their metadata.

    Only inspects entries whose filename ends with ``.cram`` (not ``.cram-metadata.json``).
    Unmatched CRAMs (``*-unmatched.cram``) are skipped — they inherently lack
    cell barcodes so ``cb_tag=False`` is expected.
    Returns list of errors for files where ``cb_tag`` is not ``True``.
    """
    errors: list[str] = []
    for filename, metadata in read_metadata.items():
        if not filename.endswith(".cram"):
            continue
        if "-unmatched.cram" in filename:
            continue
        cb_tag = metadata.get("cb_tag")
        if cb_tag is not True:
            errors.append(f"SCALE CB_TAG ERROR: {filename} has cb_tag={cb_tag}")
    return errors


def validate_scale_processed_files(
    proc_files: list[str], samples_info: dict
) -> dict[str, Any]:
    """
    Validate that expected Scale processed per-sample files are present.

    Checks under the ``samples/`` directory within ``proc_files`` for:
    - ``{sample}_anndata.h5ad`` (merged anndata per sample)
    - ``{sample}.merged.allCells.csv`` (merged cell metrics per sample)
    - ``{sample}.{sublib}_anndata.h5ad`` (per-sublibrary anndata)

    Returns ``{"errors": [...], "missing_files": [...]}``.
    """
    errors: list[str] = []
    missing_files: list[dict[str, str]] = []

    samples_dir_files: set[str] = set()
    for f in proc_files:
        if "/samples/" in f:
            samples_dir_files.add(f.split("/samples/")[-1])

    samples = samples_info.get("samples", [])
    sublibraries = samples_info.get("sublibraries", {})

    for sample in samples:
        merged_anndata = f"{sample}_anndata.h5ad"
        if merged_anndata not in samples_dir_files:
            errors.append(f"SCALE PROCESSED ERROR: missing {merged_anndata}")
            missing_files.append({"sample": sample, "file": merged_anndata})

        allcells = f"{sample}.merged.allCells.csv"
        if allcells not in samples_dir_files:
            errors.append(f"SCALE PROCESSED ERROR: missing {allcells}")
            missing_files.append({"sample": sample, "file": allcells})

        for sublib in sublibraries.get(sample, []):
            sublib_anndata = f"{sample}.{sublib}_anndata.h5ad"
            if sublib_anndata not in samples_dir_files:
                errors.append(f"SCALE PROCESSED ERROR: missing {sublib_anndata}")
                missing_files.append({"sample": sample, "file": sublib_anndata})

    return {"errors": errors, "missing_files": missing_files}
