"""
QA validation logic extracted from the QA notebook (cells 4-8).
Pure functions: no S3, no file I/O. Callers pass in gathered data and get back
structured results (errors, missing/extra file lists, etc.).
"""

from __future__ import annotations

import logging
import os
import re
from typing import Any

from qa_constants import (
    CELLRANGER_FLAT_LAYOUT,
    GEX_PAIRED_MODALITIES,
    ILLUMINA_DEFAULT_PLANNED_READS,
    ILLUMINA_FASTQ_RE,
    ILLUMINA_FLOWCELL_FILE_SUFFIXES,
    ILLUMINA_FLOWCELL_ID_RE,
    ILLUMINA_MODALITIES,
    ILLUMINA_RUN_PARAMETERS_SUFFIX,
    CELLRANGER_SUPPORTED_VERSIONS,
    CELLRANGER_VERSION_ALIASES,
    SCALE_AGGREGATE_FILE_RE,
    SCALE_AGGREGATE_FILE_RELAXED_RE,
    SCALE_RT_FILE_RE,
    SCALE_SAMPLES_FORBIDDEN_COLUMNS,
    SCALE_WAFER_MISC_RE,
    SCALE_WORKFLOW_REQUIRED_PARAMS,
    SEAHUB_TRIM_TO_BARE_SUFFIX,
    SEAHUB_BARE_OPTIONAL,
)
from qa_mods import (
    cellranger_expected,
    cellranger_ignore,
    extract_read_indicator,
    is_trimmer_stats_basename,
    make_read_partner,
    parse_raw_filename,
    parse_seahub_raw_path,
    raw_expected,
    raw_optional,
    resolve_wafer_run_id,
    seahub_stem_and_family,
)


VALID_PROBES = [
    "Chromium Mouse Transcriptome Probe Set v1.1.1",
    "Chromium Human Transcriptome Probe Set v1.1.0",
    "Chromium Human Transcriptome Probe Set v2.0.0",
    "Chromium Mouse Transcriptome Probe Set v2.0.0",
]
MIN_METADATA_READ_COUNT = 1_000_000
MIN_PCT_Q30_THRESHOLD = 60.0

# Read indicators that identify a sample index rather than a library read.
# Excluded from group read totals; see validate_read_metadata.
INDEX_READ_INDICATORS = frozenset({"I1", "I2"})


def validate_pct_q30(
    pct_q30_values: dict[str, float],
    threshold: float = MIN_PCT_Q30_THRESHOLD,
) -> list[str]:
    """
    Return error strings for sublibrary stats files with PCT_PF_Q30_bases below threshold.

    Caller passes pct_q30_values from QAGatheredData (populated during gather).
    """
    errors: list[str] = []
    for filename, q30 in sorted(pct_q30_values.items()):
        if q30 < threshold:
            errors.append(
                f"LOW Q30: {filename} PCT_PF_Q30_bases={q30:.2f} (< {threshold})"
            )
    return errors


def _fastq_count_mode(raw_assay: str) -> str:
    """Internal: assay policy for fastq count validation and summaries."""
    if raw_assay in ("sci_jumbo", "10x_viral_ORF", "10x_cram", "seahub_sci"):
        return "skip"
    if raw_assay in ("scale", "sci_plex"):
        return "gex_hash"
    if raw_assay == "10x":
        return "10x"
    return "unknown"


def _illumina_library_keys(names: list[str]) -> set[tuple[str, str, str, str]]:
    """Distinct libraries behind a list of Illumina fastq basenames.

    An Illumina library contributes one fastq per planned read per lane, so the
    raw file count is ``reads x lanes x libraries``. Comparing those counts
    across modalities reports a mismatch whenever GEX and CRI ran on a
    different number of lanes -- which is normal, since lanes are accepted as
    they arrive. Collapsing to ``(flowcell, group, assay, S-number)`` restores
    what the comparison has always meant under Ultima, where every library
    contributes exactly two files: the same number of libraries per modality.
    """
    keys: set[tuple[str, str, str, str]] = set()
    for name in names:
        m = ILLUMINA_FASTQ_RE.match(name)
        if m is not None:
            keys.add(
                (
                    m.group("flowcell"),
                    m.group("group"),
                    m.group("assay"),
                    m.group("sample"),
                )
            )
    return keys


def validate_fastq_counts(
    fastq_log: dict[str, dict[str, list[str]]],
    raw_assay: str,
    *,
    sequencing_type: str = "Ultima",
) -> list[str]:
    """
    Check for mismatching number of files across modalities per sample.

    - scale: GEX vs hash_oligo must have equal counts when both present.
    - sci_plex: GEX vs hash_oligo when both present; GEX_hash_oligo only → no check.
    - sci_jumbo: No validation; logs that sci_jumbo is not validated.
    - 10x: every modality in ``GEX_PAIRED_MODALITIES[sequencing_type]`` --
      CRI, ATAC, ADT, HTO, ADT_HTO, plus TCR and BCR on Illumina -- must have
      as many libraries as GEX in the same group, when both are present.
      GEX-only, or a paired modality without GEX, is not checked. Ultima counts
      files (each library is exactly one R1 + one R2); Illumina counts distinct
      libraries, since a library has one file per read per lane. Logs when
      CRI+ATAC are both present (future QA).
    - 10x_viral_ORF: No validation (legacy); logs that it is not validated.
    - seahub_sci: No validation. A SeaHub sublibrary carries one combined
      GEX_hash_oligo type token rather than separate GEX/hash_oligo files, so
      there is nothing to compare; logs that it is not validated.

    Returns list of error messages.
    """
    logger = logging.getLogger(__name__)
    errors: list[str] = []

    mode = _fastq_count_mode(raw_assay)
    if mode == "skip":
        if raw_assay == "sci_jumbo":
            logger.warning("Not validating fastq counts for sci_jumbo at the moment.")
        elif raw_assay == "10x_cram":
            logger.warning(
                "Not validating fastq counts for 10x_cram (CRAM-only raw layout)."
            )
        elif raw_assay == "seahub_sci":
            # A SeaHub sublibrary carries one combined GEX_hash_oligo type token,
            # never separate GEX and hash_oligo files, so this check's premise --
            # two modalities to compare -- never holds. Grouping it with scale and
            # sci_plex's gex_hash mode made every clean upload report "checked 0,
            # no comparable pairs", which reads as a defect on a compliant upload.
            logger.warning(
                "Not validating fastq counts for seahub_sci (sublibraries carry a "
                "single combined type token, not separate GEX/hash_oligo files)."
            )
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

            def _units(modality: str) -> list:
                files = v.get(modality, [])
                if sequencing_type == "Illumina":
                    # Compare libraries, not files -- see _illumina_library_keys.
                    return sorted(_illumina_library_keys(files))
                return files

            gex_list = _units("GEX")
            if v.get("CRI") and v.get("ATAC"):
                logger.warning(
                    "Unexpected: CRI+ATAC present for sample %s; "
                    "QA expansion needed in the future.",
                    sample,
                )
            # Only compare when both have at least one file: GEX-only, or a
            # paired modality without GEX, has nothing to compare against.
            for modality in GEX_PAIRED_MODALITIES[sequencing_type]:
                other = _units(modality)
                if gex_list and other and len(gex_list) != len(other):
                    errors.append(
                        f"MISMATCH FQ COUNTS: {sample}: {len(gex_list)} GEX, "
                        f"{len(other)} {modality}"
                    )

    return errors


# Column order for the libraries-by-modality table; anything else follows,
# alphabetically.
_MODALITY_DISPLAY_ORDER = (
    "GEX",
    "CRI",
    "ATAC",
    "viral_ORF",
    "ADT",
    "HTO",
    "ADT_HTO",
    "TCR",
    "BCR",
    "GEX_hash_oligo",
    "hash_oligo",
)


def libraries_by_modality(
    fastq_log: dict[str, dict[str, list[str]]],
) -> list[dict[str, Any]]:
    """Which modalities each GroupID has: one row per GroupID, one column each.

    A library is unique to its GroupID and modality, so every file of that
    modality in the group -- every wafer, flowcell, lane, S-number, UG or
    barcode -- counts as one library. Each count is therefore 0 or 1, the same
    for Ultima and Illumina.

    Every GroupID in ``fastq_log`` gets a row, including one with no files, and
    every modality seen anywhere gets a column (0 where a group has none).
    """
    counts: dict[str, dict[str, int]] = {
        group: {modality: 1 if files else 0 for modality, files in by_modality.items()}
        for group, by_modality in fastq_log.items()
    }

    seen = {modality for row in counts.values() for modality in row}
    columns = [m for m in _MODALITY_DISPLAY_ORDER if m in seen] + sorted(
        seen - set(_MODALITY_DISPLAY_ORDER)
    )
    return [
        {"group_id": group, **{m: counts[group].get(m, 0) for m in columns}}
        for group in sorted(counts)
    ]


def summarize_fastq_count_validation(
    fastq_log: dict[str, dict[str, list[str]]],
    raw_assay: str,
    errors: list[str],
    *,
    sequencing_type: str = "Ultima",
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
        if raw_assay == "10x_cram":
            return (
                "Fastq count validation (10x_cram): not applicable for CRAM-only raw "
                f"inputs; mismatches: {mismatches}. No comparisons performed."
            )
        if raw_assay == "seahub_sci":
            return (
                "Fastq count validation (seahub_sci): not applicable -- SeaHub "
                "sublibraries carry a single combined GEX_hash_oligo type token, "
                f"not separate GEX/hash_oligo files; mismatches: {mismatches}. "
                "No comparisons performed."
            )
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
        paired = GEX_PAIRED_MODALITIES[sequencing_type]
        checked_by_modality = {
            modality: sum(
                1
                for _sample, v in fastq_log.items()
                if v.get("GEX") and v.get(modality)
            )
            for modality in paired
        }
        checked_total = sum(checked_by_modality.values())

        if checked_total > 0 and mismatches == 0:
            result = "All matched."
        elif checked_total == 0 and mismatches == 0:
            gex_present = sum(1 for _sample, v in fastq_log.items() if v.get("GEX"))
            paired_present = sum(
                1
                for _sample, v in fastq_log.items()
                if any(v.get(modality) for modality in paired)
            )
            if gex_present > 0 and paired_present == 0:
                result = "Only GEX present (no paired modality to compare)."
            elif gex_present == 0 and paired_present > 0:
                result = (
                    "Missing GEX (paired modalities present, nothing to compare to)."
                )
            elif gex_present > 0 and paired_present > 0:
                result = (
                    "No comparable pairs (GEX groups differ from the groups with "
                    "paired modalities)."
                )
            else:
                result = "No comparable pairs (missing modality pairs)."
        elif mismatches > 0:
            result = (
                "Mismatches found (see preceding MISMATCH FQ COUNTS lines in stdout "
                "and the errors log)."
            )
        else:
            result = "No comparisons performed."

        # Only the pairs that actually occurred, so an 8-modality Illumina
        # list does not print seven zeros.
        breakdown = ", ".join(
            f"GEX<->{modality}: {n}" for modality, n in checked_by_modality.items() if n
        )
        return (
            f"Fastq count validation (10x): checked {checked_total}"
            + (f" ({breakdown})" if breakdown else "")
            + f"; mismatches: {mismatches}. {result}"
        )

    return (
        f"Fastq count validation ({raw_assay}): checked 0; mismatches: "
        f"{mismatches}. Unsupported assay."
    )


def validate_read_metadata(
    read_metadata: dict[str, Any],
    raw_assay: str,
    *,
    sequencing_type: str = "Ultima",
    print_success: bool = False,
    success_print_limit: int = 5,
    pairing_paths_print_limit: int = 100,
) -> tuple[dict[str, dict[str, int]], list[str], dict[str, list[str]]]:
    """
    Build group read counts from metadata and check R1/R2 consistency.

    Returns ``(group_read_counts, errors, pairing)`` where ``pairing`` has:

    - ``r1_without_r2_metadata``: R1 paths with ``read_count`` but no R2 key
    - ``r2_without_r1_metadata``: R2 paths with no corresponding R1 key

    Pairing issues are appended to ``errors`` with a ``READ METADATA PAIRING:``
    prefix. Metadata filename/content mismatches are appended as
    ``METADATA FILENAME WARNING: ...`` entries so callers can log them as
    warnings without failing count comparisons.

    R1/R2 classification uses ``extract_read_indicator()`` which matches the
    Illumina read indicator at the *end* of the filename (e.g. ``_R1_001.fastq.gz``),
    so ``_R1_`` / ``_R2_`` embedded in a group ID (e.g. ``q_pcf_R2``) is ignored.

    Index reads (I1, I2) are compared too. Every read in one library-lane set
    should carry the same count, so each index read is compared against that
    set's R1, or its R2 when there is no R1; an I2 with neither is compared
    against its I1. A mismatch is a ``READ COUNT ERROR`` in the same format as
    an R1/R2 one. A missing index file is not reported here -- which reads
    should exist is the expected-files check's job.
    """
    errors: list[str] = []
    group_read_counts: dict[str, dict[str, int]] = {}

    skip_r1_r2_pairing = raw_assay in ("10x_cram", "seahub_sci")

    # Optional instrumentation for stdout-only reporting.
    matched_examples: list[str] = []
    compared_pairs = 0
    matched_pairs = 0
    mismatched_pairs = 0
    skipped_no_r2 = 0
    r2_metadata_error = 0
    r1_without_r2_metadata: list[str] = []
    filename_mismatch_warnings = 0
    # Per-check tallies for the printed summary.
    filename_checked = 0
    errored_files = 0
    no_read_count = 0
    min_checked = 0
    min_failed = 0

    for f, meta in read_metadata.items():
        reported_filename = str(
            meta.get("__reported_filename", meta.get("filename", ""))
        )
        if reported_filename:
            filename_checked += 1
        if reported_filename and reported_filename != f:
            filename_mismatch_warnings += 1
            errors.append(
                "METADATA FILENAME WARNING: "
                f"metadata filename does not match source object; "
                f"reported={reported_filename} actual={f}"
            )
        if meta.get("errors"):
            errored_files += 1
            errors.append(
                f"METADATA.JSON ERROR: {f} has error in metadata.json:{meta['errors']}"
            )
            continue
        reads = meta.get("read_count")
        if reads is None:
            no_read_count += 1
            continue
        min_checked += 1
        if reads < MIN_METADATA_READ_COUNT:
            min_failed += 1
            errors.append(
                f"READ COUNT ERROR: {f} read_count {reads} is below minimum "
                f"{MIN_METADATA_READ_COUNT}"
            )
        indicator = extract_read_indicator(f)
        # An index read is a sequencing artifact, not a library read. Summing
        # I1 and I2 into the group total inflates raw reads about threefold on
        # a dual-indexed Illumina delivery and makes every group fail the
        # raw-vs-CellRanger comparison in validate_processed_group. Not gated
        # on pairing mode: no layout has a reason to count one, and no Ultima
        # assay ships index sidecars, so this is a no-op for every delivery QA
        # has seen. R3 is deliberately left alone -- ATAC's third read is a
        # library read and is counted today.
        if indicator in INDEX_READ_INDICATORS:
            continue
        if not skip_r1_r2_pairing and indicator == "R2":
            continue
        parsed = parse_raw_filename(f, raw_assay, sequencing_type=sequencing_type)
        if parsed is not None:
            _run, group, assay, _ug, _barcode = parsed
            if group not in group_read_counts:
                group_read_counts[group] = {assay: reads}
            elif assay not in group_read_counts[group]:
                group_read_counts[group][assay] = reads
            else:
                group_read_counts[group][assay] += reads
        if not skip_r1_r2_pairing and indicator == "R1":
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
    r2_checked = 0
    if not skip_r1_r2_pairing:
        for f in read_metadata:
            if extract_read_indicator(f) != "R2":
                continue
            r2_checked += 1
            r1file = make_read_partner(f, "R2", "R1")
            if r1file not in read_metadata:
                r2_without_r1_metadata.append(f)

    def _usable_count(meta: Any) -> Any:
        if not isinstance(meta, dict) or meta.get("errors"):
            return None
        return meta.get("read_count")

    index_compared = 0
    index_mismatched = 0
    if not skip_r1_r2_pairing:
        for f, meta in read_metadata.items():
            indicator = extract_read_indicator(f)
            if indicator not in INDEX_READ_INDICATORS:
                continue
            reads = _usable_count(meta)
            if reads is None:
                continue
            references = ("R1", "R2") + (("I1",) if indicator == "I2" else ())
            for ref_indicator in references:
                ref_file = make_read_partner(f, indicator, ref_indicator)
                ref_reads = _usable_count(read_metadata.get(ref_file))
                if ref_reads is None:
                    continue
                index_compared += 1
                if reads != ref_reads:
                    index_mismatched += 1
                    errors.append(
                        f"READ COUNT ERROR:{ref_file}-{ref_reads},{f}-{reads}"
                    )
                break

    for path in r1_without_r2_metadata:
        errors.append(f"READ METADATA PAIRING: R1 present, no R2 metadata key: {path}")
    for path in r2_without_r1_metadata:
        errors.append(f"READ METADATA PAIRING: R2 present, no R1 metadata key: {path}")

    pairing: dict[str, list[str]] = {
        "r1_without_r2_metadata": r1_without_r2_metadata,
        "r2_without_r1_metadata": r2_without_r1_metadata,
    }

    if print_success:
        rows = [
            ("metadata.json has no errors", len(read_metadata), errored_files),
            ("read_count present", len(read_metadata) - errored_files, no_read_count),
            (f"read_count >= {MIN_METADATA_READ_COUNT:,}", min_checked, min_failed),
            (
                "filename matches its S3 object",
                filename_checked,
                filename_mismatch_warnings,
            ),
        ]
        if not skip_r1_r2_pairing:
            rows += [
                (
                    "R1 has R2 metadata",
                    compared_pairs + r2_metadata_error + skipped_no_r2,
                    skipped_no_r2,
                ),
                ("R2 has R1 metadata", r2_checked, len(r2_without_r1_metadata)),
                ("R1 count == R2 count", compared_pairs, mismatched_pairs),
                ("I1/I2 count == R1 (or R2)", index_compared, index_mismatched),
            ]
        _print_read_metadata_summary(raw_assay, len(read_metadata), rows)
        if skip_r1_r2_pairing:
            print(
                "  R1/R2 and index read comparisons: not applicable (CRAM-only layout)"
            )
        for line in matched_examples:
            print(line)
        if not skip_r1_r2_pairing:
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


def _print_read_metadata_summary(
    raw_assay: str, total_files: int, rows: list[tuple[str, int, int]]
) -> None:
    """Print one line per read-metadata check: files checked, passed, failed.

    Which files failed is in the error lines (``READ COUNT ERROR``,
    ``METADATA.JSON ERROR``, ...), not here.
    """
    print(f"validate_read_metadata({raw_assay}): {total_files} metadata files")
    width = max(len(name) for name, _, _ in rows)
    print(f"  {'Check':<{width}}  {'Checked':>7}  {'Pass':>7}  {'Fail':>7}")
    for name, checked, failed in rows:
        print(f"  {name:<{width}}  {checked:>7,}  {checked - failed:>7,}  {failed:>7,}")


def _check_expected_illumina_raw_files(
    all_raw_files: list[str],
    planned_reads: dict[str, list[str]] | None,
) -> tuple[int, list[dict[str, Any]], list[str]]:
    """
    Completeness for an Illumina delivery. Same contract as its Ultima sibling:
    ``(all_good_count, raw_lost, raw_found)``.

    An Illumina ``raw/`` directory is checked on two axes, because it holds two
    kinds of thing:

    * **Per flowcell** -- the eight run files
      (``{FC}_CopyComplete.txt`` ... ``{FC}_Logs.zip``). One set per flowcell,
      not one per library, so these cannot be expressed as endings on a
      library "beginning" the way the Ultima tables are.
    * **Per library and lane** -- one fastq for each planned read, plus that
      fastq's ``-metadata.json`` sidecar. A library here is a distinct
      ``(flowcell, group, assay, S-number)``; the lanes it arrived on are
      whatever is present, since lanes are reported rather than enforced.

    Which reads are required comes from ``planned_reads`` (parsed from each
    flowcell's RunParameters.xml), NOT from a fixed R1/R2 assumption: a
    single-index run plans three reads and a dual-indexed one plans four. A
    flowcell missing from ``planned_reads`` falls back to
    ``ILLUMINA_DEFAULT_PLANNED_READS``; the gather step has already warned.
    """
    planned_reads = planned_reads or {}
    raw_found_set = set(all_raw_files)

    # raw_dir -> flowcell -> {(group, assay, sample, lane, chunk)}
    flowcells: dict[str, set[str]] = {}
    libraries: dict[tuple[str, str], set[tuple[str, str, str, str, str]]] = {}

    for fullpath in all_raw_files:
        raw_dir = "/".join(fullpath.split("/")[:-1])
        basename = fullpath.split("/")[-1]

        m = ILLUMINA_FASTQ_RE.match(basename)
        if m is not None:
            fc = m.group("flowcell")
            flowcells.setdefault(raw_dir, set()).add(fc)
            libraries.setdefault((raw_dir, fc), set()).add(
                (
                    m.group("group"),
                    m.group("assay"),
                    m.group("sample"),
                    m.group("lane"),
                    m.group("chunk"),
                )
            )
            continue

        for suffix in ILLUMINA_FLOWCELL_FILE_SUFFIXES:
            if basename.endswith(suffix):
                candidate = basename[: -len(suffix)]
                # A flowcell ID is alphanumeric, which is what lets the fastq
                # regex split it off at the first underscore. Without this
                # guard a stray "FOO-BAR_Logs.zip" would register flowcell
                # "FOO-BAR" and invent seven missing artifacts for it; it
                # belongs in the extra-files list instead.
                if ILLUMINA_FLOWCELL_ID_RE.match(candidate):
                    flowcells.setdefault(raw_dir, set()).add(candidate)
                break

    all_good = 0
    raw_lost: list[dict[str, Any]] = []
    raw_found: list[str] = []

    def _close(row: dict[str, Any]) -> None:
        nonlocal all_good
        if len(row) > 1:
            raw_lost.append(row)
        else:
            all_good += 1

    for raw_dir in sorted(flowcells):
        for fc in sorted(flowcells[raw_dir]):
            row: dict[str, Any] = {"path": f"{raw_dir}/{fc}"}
            for suffix in ILLUMINA_FLOWCELL_FILE_SUFFIXES:
                f = f"{raw_dir}/{fc}{suffix}"
                if f in raw_found_set:
                    raw_found.append(f)
                else:
                    row[suffix] = f
            _close(row)

            reads = list(planned_reads.get(fc) or ILLUMINA_DEFAULT_PLANNED_READS)
            # One row per library+lane: a lane that lost a read is the unit an
            # operator has to go back to the vendor about.
            seen = libraries.get((raw_dir, fc), set())
            units = sorted({(g, a, s, lane, chunk) for g, a, s, lane, chunk in seen})
            for group, assay, sample, lane, chunk in units:
                stem = f"{fc}_{group}_{assay}_S{sample}_L{lane}"
                lib_row: dict[str, Any] = {"path": f"{raw_dir}/{stem}"}
                for read in reads:
                    fastq = f"{raw_dir}/{stem}_{read}_{chunk}.fastq.gz"
                    sidecar = f"{fastq}-metadata.json"
                    if fastq in raw_found_set:
                        raw_found.append(fastq)
                        # Only chase the sidecar when its fastq arrived;
                        # otherwise one absent read is reported twice.
                        if sidecar in raw_found_set:
                            raw_found.append(sidecar)
                        else:
                            lib_row[f"{read}-metadata.json"] = sidecar
                    else:
                        lib_row[read] = fastq
                _close(lib_row)

    return all_good, raw_lost, raw_found


def _illumina_r1_read_counts(
    read_metadata: dict[str, Any] | None,
) -> dict[tuple[str, str], int]:
    """Sum metadata.json ``read_count`` over every R1 fastq, by (raw_dir, flowcell).

    R1 only: each library-lane has exactly one R1, and its count is the number
    of reads in that library-lane. Adding R2 or the index reads would count the
    same reads again.

    Each entry is located by ``__source_key``, the sidecar's full S3 key set by
    the gather step -- the dict key itself is a basename or an s3:// URI
    depending on the payload, so it cannot be relied on for the directory.
    Entries with errors or no ``read_count`` are skipped, as in
    validate_read_metadata.
    """
    totals: dict[tuple[str, str], int] = {}
    for key, meta in (read_metadata or {}).items():
        if not isinstance(meta, dict) or meta.get("errors"):
            continue
        reads = meta.get("read_count")
        if not isinstance(reads, int) or isinstance(reads, bool):
            continue
        source = str(meta.get("__source_key") or key)
        if source.startswith("s3://"):
            source = source.split("/", 3)[-1]
        if source.endswith("-metadata.json"):
            source = source[: -len("-metadata.json")]
        m = ILLUMINA_FASTQ_RE.match(source.split("/")[-1])
        if m is None or m.group("read") != "R1":
            continue
        raw_dir = "/".join(source.split("/")[:-1])
        k = (raw_dir, m.group("flowcell"))
        totals[k] = totals.get(k, 0) + reads
    return totals


def format_read_count(value: Any) -> str | None:
    """Human-readable read count for the flowcell table.

    Below 100,000 the exact number is kept. From 100,000 up it is rounded to
    three significant figures and shown in millions (M) or billions (B):
    12,345,678 -> "12.3M", 250,000 -> "0.250M", 1,234,567,890 -> "1.23B".
    Rounding happens before the unit is chosen, so 999,999,999 becomes
    "1.00B" rather than "1000M". None (or NaN) stays None.
    """
    if value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    if x != x:  # NaN
        return None
    if abs(x) < 100_000:
        return str(int(x))
    rounded = float(f"{x:.3g}")
    unit, suffix = (1e9, "B") if abs(rounded) >= 1e9 else (1e6, "M")
    v = rounded / unit
    mag = abs(v)
    decimals = 0 if mag >= 100 else 1 if mag >= 10 else 2 if mag >= 1 else 3
    return f"{v:.{decimals}f}{suffix}"


# Display order for read types in the flowcell table: the library reads, then
# the index reads.
ILLUMINA_READ_TYPE_ORDER = ("R1", "R2", "R3", "I1", "I2")


def _illumina_observed_read_lengths(
    read_metadata: dict[str, Any] | None,
) -> dict[tuple[str, str], dict[str, set[int]]]:
    """Longest read in each fastq, by (raw_dir, flowcell) and read type.

    Read from each sidecar's ``read_lengths`` histogram (``{"151": n, ...}``).
    The longest length is used because trimming only shortens reads: an
    untrimmed file's longest read is the cycle count, and a trimmed one's
    spread of shorter reads must not read as a mismatch. Kept as a set of
    per-file values so one short file is not hidden by its siblings. Entries
    with errors, or without a usable histogram, are skipped. Located by
    ``__source_key``, as in ``_illumina_r1_read_counts``.
    """
    observed: dict[tuple[str, str], dict[str, set[int]]] = {}
    for key, meta in (read_metadata or {}).items():
        if not isinstance(meta, dict) or meta.get("errors"):
            continue
        histogram = meta.get("read_lengths")
        if not isinstance(histogram, dict):
            continue
        lengths = []
        for length, count in histogram.items():
            try:
                if int(count) > 0:
                    lengths.append(int(length))
            except (TypeError, ValueError):
                continue
        if not lengths:
            continue
        source = str(meta.get("__source_key") or key)
        if source.startswith("s3://"):
            source = source.split("/", 3)[-1]
        if source.endswith("-metadata.json"):
            source = source[: -len("-metadata.json")]
        m = ILLUMINA_FASTQ_RE.match(source.split("/")[-1])
        if m is None:
            continue
        k = ("/".join(source.split("/")[:-1]), m.group("flowcell"))
        observed.setdefault(k, {}).setdefault(m.group("read"), set()).add(max(lengths))
    return observed


def summarize_illumina_flowcells(
    all_raw_files: list[str],
    planned_reads: dict[str, list[str]] | None = None,
    read_metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """
    One row per GroupID: what the RunParameters.xml planned, and what actually
    arrived.

    A GroupID sequenced on several flowcells gets one row: ``flowcell`` lists
    them, ``library_count`` and ``read_count`` are summed, and every other
    column is the unique values across them. ``missing_read_type`` is still
    checked per flowcell, against that flowcell's own planned reads, before
    being combined.

    This is the Illumina counterpart of the Ultima wafer table -- the
    per-delivery summary an operator reads before looking at any error file.
    Lanes are REPORTED, never enforced: a library may legitimately arrive on
    any number of lanes, so the set is information rather than a check.

    ``modality`` is the distinct modalities (GEX, CRI, ADT_HTO, ...) named in
    the row's fastq filenames, in ILLUMINA_MODALITIES order.

    ``observed_read_types`` lists every read type (R1, R2, I1, I2) that has at
    least one fastq for the row. ``missing_read_type`` is the one check in the
    table: each library-lane set that arrived -- every distinct
    ``{FC}_{GroupID}_{Assay}_S{n}_L00#`` stem -- should have a fastq for each
    planned read, and any planned read type absent from any set is listed.
    "none" when every set is complete. A flowcell with no fastqs at all lists
    every planned read, since none arrived. Which files are missing is in
    ``{order}_raw_missing.csv``.

    ``observed_read_length`` is the longest read in each fastq, from its
    sidecar's ``read_lengths`` histogram -- one value per read type when every
    file agrees, all distinct values ("I1=8/10") when not.
    """
    planned_reads = planned_reads or {}
    r1_reads = _illumina_r1_read_counts(read_metadata)
    observed_lengths = _illumina_observed_read_lengths(read_metadata)
    lanes: dict[tuple[str, str], set[str]] = {}
    libraries: dict[tuple[str, str], set[tuple[str, str, str]]] = {}
    # (raw_dir, flowcell) -> {(stem, chunk)}: the library-lane sets that arrived
    lane_sets: dict[tuple[str, str], set[tuple[str, str]]] = {}
    # (raw_dir, flowcell) -> read types with at least one fastq
    found_types: dict[tuple[str, str], set[str]] = {}
    # (raw_dir, flowcell) -> modalities named in the fastq filenames
    modalities: dict[tuple[str, str], set[str]] = {}
    present: set[str] = set()
    seen: set[tuple[str, str]] = set()

    for fullpath in all_raw_files:
        raw_dir = "/".join(fullpath.split("/")[:-1])
        basename = fullpath.split("/")[-1]

        m = ILLUMINA_FASTQ_RE.match(basename)
        if m is not None:
            key = (raw_dir, m.group("flowcell"))
            seen.add(key)
            lanes.setdefault(key, set()).add(f"L{m.group('lane')}")
            libraries.setdefault(key, set()).add(
                (m.group("group"), m.group("assay"), m.group("sample"))
            )
            stem = (
                f"{m.group('flowcell')}_{m.group('group')}_{m.group('assay')}"
                f"_S{m.group('sample')}_L{m.group('lane')}"
            )
            lane_sets.setdefault(key, set()).add((stem, m.group("chunk")))
            found_types.setdefault(key, set()).add(m.group("read"))
            modalities.setdefault(key, set()).add(m.group("assay"))
            present.add(f"{raw_dir}/{basename}")
            continue

        for suffix in ILLUMINA_FLOWCELL_FILE_SUFFIXES:
            if basename.endswith(suffix):
                candidate = basename[: -len(suffix)]
                if ILLUMINA_FLOWCELL_ID_RE.match(candidate):
                    seen.add((raw_dir, candidate))
                break

    def _ordered(types: set[str]) -> str:
        return ", ".join(t for t in ILLUMINA_READ_TYPE_ORDER if t in types)

    def _lengths_text(by_type: dict[str, set[int]]) -> str:
        return ", ".join(
            f"{t}={'/'.join(str(n) for n in sorted(by_type[t]))}"
            for t in ILLUMINA_READ_TYPE_ORDER
            if t in by_type
        )

    def _missing_read_types(key: tuple[str, str], reads: list[str]) -> set[str]:
        sets = lane_sets.get(key, set())
        if not sets:
            return set(reads)
        raw_dir = key[0]
        return {
            read
            for stem, chunk in sets
            for read in reads
            if f"{raw_dir}/{stem}_{read}_{chunk}.fastq.gz" not in present
        }

    # One row per GroupID folder; a group sequenced on several flowcells has
    # each flowcell's values combined -- listed, uniqued or summed.
    by_dir: dict[str, list[str]] = {}
    for raw_dir, flowcell in seen:
        by_dir.setdefault(raw_dir, []).append(flowcell)

    rows: list[dict[str, Any]] = []
    for raw_dir in sorted(by_dir):
        flowcells = sorted(by_dir[raw_dir])
        keys = [(raw_dir, fc) for fc in flowcells]
        # {proj}/{Order#}/{GroupID}/raw -> GroupID.
        parts = raw_dir.split("/")
        group_id = parts[-2] if len(parts) >= 2 and parts[-1] == "raw" else ""

        def _union(table: dict, keys=keys) -> set:
            return set().union(*(table.get(k, set()) for k in keys))

        # Each flowcell is checked against its own planned reads (all four
        # read types when its RunParameters.xml is unusable, as in the
        # expected-files check).
        missing: set[str] = set()
        for key in keys:
            reads = planned_reads.get(key[1]) or ILLUMINA_DEFAULT_PLANNED_READS
            missing |= _missing_read_types(key, list(reads))

        expected = set().union(*(planned_reads.get(fc) or () for fc in flowcells))
        unknown = [fc for fc in flowcells if not planned_reads.get(fc)]
        if not unknown:
            expected_text = _ordered(expected)
        else:
            note = (
                f"no readable {ILLUMINA_RUN_PARAMETERS_SUFFIX.lstrip('_')} for "
                f"{', '.join(unknown)}; assumed "
                + ", ".join(ILLUMINA_DEFAULT_PLANNED_READS)
            )
            expected_text = (
                f"{_ordered(expected)}; UNKNOWN ({note})"
                if expected
                else f"UNKNOWN ({note})"
            )

        observed: dict[str, set[int]] = {}
        for key in keys:
            for read_type, values in observed_lengths.get(key, {}).items():
                observed.setdefault(read_type, set()).update(values)

        group_modalities = _union(modalities)
        rows.append(
            {
                "group_id": group_id,
                "flowcell": ", ".join(flowcells),
                "modality": ", ".join(
                    mod for mod in ILLUMINA_MODALITIES if mod in group_modalities
                ),
                "expected_read_types": expected_text,
                "observed_read_length": _lengths_text(observed),
                "observed_read_types": _ordered(_union(found_types)),
                "missing_read_type": _ordered(missing) or "none",
                "lanes": ", ".join(sorted(_union(lanes))),
                "library_count": sum(len(libraries.get(k, set())) for k in keys),
                "read_count": sum(r1_reads.get(k, 0) for k in keys),
            }
        )
    return rows


def total_reads_by_flowcell(
    all_raw_files: list[str],
    read_metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Total R1 read count per flowcell, across every GroupID on it.

    Sums the R1 ``read_count`` of every fastq ``metadata.json`` by flowcell, so
    a flowcell shared by several GroupIDs gets one total. A flowcell named by a
    fastq or a run file but with no readable R1 count is listed with 0. Sorted
    by flowcell ID.
    """
    totals: dict[str, int] = {}
    for fullpath in all_raw_files:
        basename = fullpath.split("/")[-1]
        m = ILLUMINA_FASTQ_RE.match(basename)
        if m is not None:
            totals.setdefault(m.group("flowcell"), 0)
            continue
        for suffix in ILLUMINA_FLOWCELL_FILE_SUFFIXES:
            if basename.endswith(suffix):
                candidate = basename[: -len(suffix)]
                if ILLUMINA_FLOWCELL_ID_RE.match(candidate):
                    totals.setdefault(candidate, 0)
                break
    for (_raw_dir, flowcell), n in _illumina_r1_read_counts(read_metadata).items():
        totals[flowcell] = totals.get(flowcell, 0) + n
    return [{"flowcell": fc, "read_count": n} for fc, n in sorted(totals.items())]


# Missing raw files with these endings are sequence data; everything else a
# delivery is expected to carry (metadata.json sidecars, trimmer stats, run
# logs and flowcell files, ...) is reported as a log file.
RAW_SEQUENCE_SUFFIXES = (".fastq.gz", ".cram")


def split_missing_raw_files(
    raw_lost: list[dict[str, Any]],
) -> tuple[list[str], list[str]]:
    """Split the missing paths in ``raw_lost`` into (sequence files, log files).

    A sequence file ends in ``.fastq.gz`` or ``.cram``; a ``-metadata.json``
    sidecar of one is a log file. Each row of ``raw_lost`` is
    ``{"path": <beginning>, <ending>: <full missing path>, ...}`` from
    ``check_expected_raw_files``.
    """
    sequence: list[str] = []
    logs: list[str] = []
    for row in raw_lost:
        for key, path in row.items():
            if key == "path":
                continue
            (sequence if str(path).endswith(RAW_SEQUENCE_SUFFIXES) else logs).append(
                path
            )
    return sequence, logs


def check_expected_raw_files(
    all_raw_files: list[str],
    raw_assay: str,
    *,
    allow_truncated_stats_name: bool = False,
    sequencing_type: str = "Ultima",
    planned_reads: dict[str, list[str]] | None = None,
) -> tuple[int, list[dict[str, Any]], list[str]]:
    """
    Check that expected raw file endings are present for each "beginning".
    Returns (all_good_count, raw_lost list of dicts, raw_found list of paths).

    When ``allow_truncated_stats_name`` is True, a missing
    ``<beginning>_trimmer-stats.csv`` is treated as fulfilled if the sibling
    ``<beginning>_stats.csv`` is present in ``all_raw_files``; that sibling
    is then added to ``raw_found`` so it is not double-flagged as extra.

    ``sequencing_type="Illumina"`` delegates to
    ``_check_expected_illumina_raw_files``: an Illumina delivery is not shaped
    like a set of endings hung off a library beginning, so it gets its own
    traversal rather than an entry in ``raw_expected``.
    """
    if sequencing_type == "Illumina":
        return _check_expected_illumina_raw_files(all_raw_files, planned_reads)

    # Keyed by (raw_dir, beginning) for seahub_sci, and by the beginning alone
    # for everything else.
    #
    # A SeaHub beginning is unique inside a wafer folder and not across them, and
    # the value carries the folder every expected path is built from -- so with
    # one key the second folder's copy of a well kept the *first* folder's
    # raw_dir and its artifacts matched nothing. Measured on a well uploaded
    # under two sublibrary folders: two complete wells reported as one complete
    # well plus five unexpected files, because the second copy's required five
    # never reached raw_found and check_extra_raw_files takes everything not
    # there. That upload is already reported as duplicate_trimmed_well, and the
    # roll-up calls the well UNKNOWN; the inventory was the one output
    # describing it wrongly rather than merely differently.
    #
    # Only for seahub_sci. Per-directory keying was tried for every assay and is
    # not a no-op: `scale` derives its group from the key rather than the folder,
    # and _gather_group_raw unions several run subfolders under one group, so
    # sibling folders sharing a beginning would each become a well. Those layouts
    # are outside this mode and have no SOP notion of a well directory to check
    # against.
    beginnings: dict[Any, dict[str, Any]] = {}
    for fullpath in all_raw_files:
        raw_dir = "/".join(fullpath.split("/")[:-1])
        if raw_assay == "seahub_sci":
            parsed_stem = seahub_stem_and_family(fullpath.split("/")[-1])
            if parsed_stem is None:
                continue
            # An object outside a wafer folder belongs to no well. Registering
            # one would invent a well at a path that must not exist -- reporting
            # its other four artifacts as missing from there, and taking the
            # stray itself out of the extra list, which is the only place it
            # should appear. The recursive raw/ walk collects these deliberately
            # so bad_path_depth can see them, and roll_up_wells counts them as
            # unaccounted for the same reason.
            if parse_seahub_raw_path(fullpath) is None:
                continue
            b, _suffix, _family = parsed_stem
            # Completeness asks only whether each of the five artifact *kinds*
            # arrived; the spelling is a separate axis, reported once per stem by
            # the SOP validator as missing_trim_infix. Judging by family instead
            # let one optional sidecar with the other family's name decide the
            # whole requirement set, so a complete well could be reported as
            # missing five files that were never meant to exist -- and this path
            # and roll_up_wells disagreed about which five.
            beginnings.setdefault(
                (raw_dir, b),
                {
                    "raw_dir": raw_dir,
                    "beginning": b,
                    "endings": list(raw_expected[raw_assay]),
                },
            )
            continue

        parsed = parse_raw_filename(fullpath, raw_assay)
        if parsed is None:
            continue
        run, group, assay, ug, barcode = parsed
        b = f"{run}-{group}_{assay}-{ug}-{barcode}"
        if b not in beginnings:
            endings = list(raw_expected.get(raw_assay, []))
            if raw_assay == "10x_viral_ORF" and assay == "GEX":
                endings = list(raw_expected.get("10x", []))
            beginnings[b] = {"raw_dir": raw_dir, "beginning": b, "endings": endings}

    all_good = 0
    raw_lost: list[dict[str, Any]] = []
    raw_found: list[str] = []
    raw_found_set = set(all_raw_files)

    for v in beginnings.values():
        raw_dir, b = v["raw_dir"], v["beginning"]
        # "path" stays the beginning alone, so the CSV keeps its shape -- but for
        # seahub_sci it is no longer unique across rows, since two folders can
        # each contribute one. The per-ending columns hold full paths, which is
        # what tells those rows apart.
        temp_missing: dict[str, Any] = {"path": b}
        for e in v["endings"]:
            f = f"{raw_dir}/{b}{e}"
            # A SeaHub artifact counts under either spelling: the SOP name or the
            # same file with the ".trim" infix dropped. Only the SOP name is
            # reported as missing, since that is what should exist.
            if raw_assay == "seahub_sci" and f not in raw_found_set:
                bare = SEAHUB_TRIM_TO_BARE_SUFFIX.get(e)
                if bare and f"{raw_dir}/{b}{bare}" in raw_found_set:
                    raw_found.append(f"{raw_dir}/{b}{bare}")
                    continue
            if f not in raw_found_set:
                if e.endswith("-metadata.json") and (
                    f.replace("-metadata.json", "") not in raw_found_set
                ):
                    continue
                if allow_truncated_stats_name and e == "_trimmer-stats.csv":
                    alias = f.replace("_trimmer-stats.csv", "_stats.csv")
                    if alias in raw_found_set:
                        raw_found.append(alias)
                        continue
                temp_missing[e] = f
            else:
                raw_found.append(f)
        if len(temp_missing) > 1:
            raw_lost.append(temp_missing)
        else:
            all_good += 1

    return all_good, raw_lost, raw_found


def _is_known_scale_raw_file(
    filepath: str, *, allow_truncated_stats_name: bool = False
) -> bool:
    """Return True if *filepath* matches a recognised Scale raw file pattern.

    Scale filenames don't follow the ``{beginning}{suffix}`` layout used by
    10x / sci assays, so they are validated here with compiled regexes
    (see ``qa_constants.SCALE_*_RE``).

    Recognised categories:
    * Per-RT files  – one per reaction tube / well (e.g. ``_5B.cram`` or ``-5B.cram``)
    * Aggregate files – trimmer stats, failure codes, unmatched reads
    * Wafer-level misc – ``SequencingInfo.json``, ``LibraryInfo.xml``,
      ``merged_trimmer-*.csv``

    ``.cram-metadata.json`` sidecars are handled separately by the existing
    metadata-sidecar logic in :func:`check_extra_raw_files`.

    When ``allow_truncated_stats_name`` is True, the relaxed aggregate regex
    additionally accepts ``..._stats.csv`` / ``...-stats.csv`` as aliases of
    ``..._trimmer-stats.csv``.
    """
    filename = filepath.split("/")[-1]
    if SCALE_RT_FILE_RE.search(filename):
        return True
    aggregate_re = (
        SCALE_AGGREGATE_FILE_RELAXED_RE
        if allow_truncated_stats_name
        else SCALE_AGGREGATE_FILE_RE
    )
    # Aggregate matcher must not absorb wafer-level ``merged_*`` filenames; the
    # wafer regex is the canonical authority for those (and intentionally
    # rejects truncated aliases like ``merged_stats.csv``).
    if not filename.startswith(("merged_", "merged-")) and aggregate_re.search(
        filename
    ):
        return True
    if SCALE_WAFER_MISC_RE.match(filename):
        return True
    return False


def count_files_by_type(paths: list[str]) -> list[tuple[str, int]]:
    """Count files by type, most common first (ties alphabetical).

    A file's type is everything from the first "." of its name, so
    ``x_R1_001.fastq.gz`` is ``.fastq.gz`` and its sidecar is
    ``.fastq.gz-metadata.json``. A name with no "." is ``(no extension)``.
    """
    counts: dict[str, int] = {}
    for path in paths:
        name = path.split("/")[-1]
        kind = name[name.index(".") :] if "." in name else "(no extension)"
        counts[kind] = counts.get(kind, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def check_extra_raw_files(
    all_raw_files: list[str],
    raw_found: list[str],
    raw_assay: str,
    *,
    allow_truncated_stats_name: bool = False,
    sequencing_type: str = "Ultima",
) -> list[str]:
    """
    Identify raw files that are not in the expected set (optionals allowed).
    Returns list of extra file paths.

    When ``allow_truncated_stats_name`` is True, files ending in
    ``_stats.csv`` / ``-stats.csv`` are accepted: for Scale via the relaxed
    aggregate regex, and for non-Scale assays only when a sibling
    ``<beginning>_trimmer-stats.csv`` is part of the expected set (i.e. the
    same beginning appears elsewhere in ``all_raw_files``).
    """
    raw_found_set = set(raw_found)
    all_raw_set = set(all_raw_files)

    if sequencing_type == "Illumina":
        # There are no optional Illumina endings to reconstruct: the Illumina
        # completeness pass already put every flowcell file, fastq and sidecar
        # it matched into raw_found, so whatever is left really is unexpected.
        # The sidecar rule still applies, for the same reason it does below --
        # a sidecar beside a file we know about is not itself a stray.
        return [
            f
            for f in all_raw_files
            if f not in raw_found_set
            and not (
                f.endswith("-metadata.json")
                and f.replace("-metadata.json", "") in all_raw_set
            )
        ]

    optional_endings = raw_optional.get(raw_assay, [])
    if raw_assay == "10x_viral_ORF":
        optional_endings = raw_optional.get("10x", [])

    # Beginnings registered via parseable raw files (used for non-Scale alias).
    registered_beginnings: set[str] = set()
    if allow_truncated_stats_name and raw_assay != "scale":
        for path in all_raw_files:
            parsed = parse_raw_filename(path, raw_assay)
            if parsed is None:
                continue
            run, group, assay, ug, barcode = parsed
            registered_beginnings.add(f"{run}-{group}_{assay}-{ug}-{barcode}")

    extra: list[str] = []
    for f in all_raw_files:
        if f in raw_found_set:
            continue
        if f.endswith("-metadata.json") and (
            f.replace("-metadata.json", "") in all_raw_set
        ):
            continue
        if raw_assay == "scale" and _is_known_scale_raw_file(
            f, allow_truncated_stats_name=allow_truncated_stats_name
        ):
            continue
        if raw_assay == "seahub_sci":
            # Optional endings are family-specific, and the generic
            # parse_raw_filename path below cannot rebuild a SeaHub stem (its
            # group is the folder sublibrary, not the filename token), so the
            # SeaHub decision is made entirely here.
            parsed_stem = seahub_stem_and_family(f.split("/")[-1])
            if parsed_stem is not None:
                _stem, suffix, family = parsed_stem
                family_optional = (
                    raw_optional.get(raw_assay, [])
                    if family == "trim"
                    else SEAHUB_BARE_OPTIONAL
                )
                if suffix in family_optional:
                    continue
            extra.append(f)
            continue
        if (raw_assay in raw_optional or raw_assay == "10x_viral_ORF") and (
            optional_endings
        ):
            parsed = parse_raw_filename(f, raw_assay)
            if parsed is not None:
                run, group, assay, ug, barcode = parsed
                b = f"{run}-{group}_{assay}-{ug}-{barcode}"
                raw_dir = "/".join(f.split("/")[:-1])
                # The reconstructed beginning is not always the literal filename
                # prefix: parse_raw_filename's fallback can duplicate a token
                # (`439047-G1-Z0169-ACGT` parses to `439047-G1_G1-Z0169-ACGT`).
                # Slicing by length would then take a misaligned tail that can
                # coincidentally match an optional ending — silently accepting a
                # file that belongs in the extra list.
                prefix = f"{raw_dir}/{b}"
                if not f.startswith(prefix):
                    extra.append(f)
                    continue
                suffix = f[len(prefix) :]
                endings = (
                    raw_optional["10x"]
                    if (raw_assay == "10x_viral_ORF" and assay == "GEX")
                    else optional_endings
                )
                if suffix in endings:
                    continue
        if (
            allow_truncated_stats_name
            and raw_assay != "scale"
            and is_trimmer_stats_basename(f, allow_truncated=True)
            and not is_trimmer_stats_basename(f, allow_truncated=False)
        ):
            basename = f.split("/")[-1]
            sep = "_" if basename.endswith("_stats.csv") else "-"
            beginning = basename[: -len(f"{sep}stats.csv")]
            if beginning and beginning in registered_beginnings:
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

    # There is a specific cellranger version (10.0.0-5-g8638ac84de) that was run in the cloud that had no impact other than online documentation
    software = CELLRANGER_VERSION_ALIASES.get(software, software)
    if software == "cellranger-9.0.1":
        sub = report.get("sub", "")
    elif software in CELLRANGER_FLAT_LAYOUT:
        sub = "multi"
    else:
        expected_versions = ", ".join(
            v.replace("cellranger-", "") for v in CELLRANGER_SUPPORTED_VERSIONS
        )
        errors.append(
            f"CR ERROR: {group_name} version is {software} but should be "
            f"{expected_versions}"
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

    flat_layout = software in CELLRANGER_FLAT_LAYOUT

    if "CRISPR" in extra or "Antibody" in extra:
        if not flat_layout:
            expected.append("multi/count/feature_reference.csv")
            per_samp_expected.append("count/feature_reference.csv")
        else:
            expected.append("feature_reference.csv")
        if "CRISPR" in extra and not flat_layout:
            per_samp_expected.append("count/crispr_analysis.tar.gz")
        if "Antibody" in extra:
            per_samp_expected.append("count/antibody_analysis.tar.gz")

    if sub == "multi" and report.get("multiplex") and not flat_layout:
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

    # If tar is in actual_set and not found in expected and per_samp_expected:
    # - add tar files to expected and
    # - drop trio of files from expected
    trio = [
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
    ]
    prefixes = ["filtered_feature_bc_matrix", "raw_feature_bc_matrix"]
    for prefix in prefixes:
        if prefix + ".tar.gz" in actual_set and prefix + ".tar.gz" not in expected:
            expected = [
                f for f in expected if f not in [os.path.join(prefix, i) for i in trio]
            ]
            expected.append(prefix + ".tar.gz")
            sample_prefix = "sample_" + prefix
            per_samp_expected = [
                f
                for f in per_samp_expected
                if f not in [os.path.join(sample_prefix, i) for i in trio]
            ]
            per_samp_expected.append(sample_prefix + ".tar.gz")

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

    # Create list of extra files, but ignore any file that contains any of the strings in cellranger_ignore
    extra_files = [
        f
        for f in actual
        if f not in expected_set_for_extra
        and not f.endswith("manifest.json")
        and not any(sub in f for sub in cellranger_ignore)
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
        elif chem == "Flex Gene Expression" and not re.search(
            "v1.1", report["Probe Set Name"]
        ):
            errors.append(
                f"WARNING: Invalid probe set for this assay version {group_name}: "
                f"{report['Probe Set Name']}"
            )
        elif chem == "GEM-X Flex v2" and not re.search(
            "v2.0.0", report["Probe Set Name"]
        ):
            errors.append(
                f"WARNING: Invalid probe set for this assay version {group_name}: "
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
    Aggregate trimmer failure statistics into wafer-level stats keyed by RunID.

    ``trimmer_failure_stats`` may be keyed either by legacy experiment id
    (``"/".join(s3_key.split("/")[1:3])``) or directly by RunID when the
    gatherer uses ``trimmer_failure_storage_key`` (multiple wafers per
    sublibrary).

    Args:
        trimmer_failure_stats:
            Dict mapping stats key -> {"rsq": [...], "trimmer_fail": [...]}.
        exp_to_run_map:
            Dict mapping legacy experiment ids (and optionally RunIDs) to RunID
            strings. Used via ``resolve_wafer_run_id`` for backward
            compatibility.

    Returns:
        Dict mapping RunID (wafer) -> {"rsq": [...], "trimmer_fail": [...]}.
    """
    wafer_failure_stats: dict[str, dict[str, list[float]]] = {}

    for stats_key, stats in trimmer_failure_stats.items():
        run_id = resolve_wafer_run_id(stats_key, exp_to_run_map)
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
        if "-unmatched.cram" in filename or "_unmatched.cram" in filename:
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
