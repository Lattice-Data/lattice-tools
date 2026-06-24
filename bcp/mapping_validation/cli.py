from __future__ import annotations

import argparse
import sys
from typing import List, Set

from .completeness import (
    validate_10x_raw_fastq_completeness,
    validate_scale_raw_completeness,
    validate_sci_raw_completeness,
)
from .parsing import parse_mapping_file
from .reporting import (
    DEFAULT_MAX_EXAMPLES,
    cap,
    configure_reporting,
    count_issue_type as _count_issue_type,
    print_issue_examples as _print_issue_examples,
    report_sif_comparison as _report_sif_comparison,
)
from .sif_io import (
    _normalize_sif_groupid,
    load_sif_group_assays,
    load_sif_library_assays,
)
from .uniqueness import validate_uniqueness
from .validators import (
    compare_groupid_assays,
    find_unmatched_sif_paths_10x,
    validate_10x_raw_file_modalities,
    validate_library_assay_consistency,
    validate_local_paths_scale_raw,
    validate_local_paths_sci_raw,
    validate_s3_10x_cram_raw,
    validate_s3_10x_processed,
    validate_10x_multiome_processed_outs,
    validate_s3_10x_raw,
    validate_10x_raw_fastq_read_mates,
    validate_s3_local_consistency_10x_processed,
    validate_s3_local_consistency_scale,
    validate_s3_local_consistency_sci,
    validate_s3_seahub_raw,
    validate_sif_completeness_10x_processed,
    validate_sif_completeness_seahub,
)


def main() -> None:
    """
    CLI for validating mapping CSV/TSV files against SOP rules.

    Required arguments:
        --mapping PATH       Mapping file with two columns: S3, local
        --provider {novogene,psomagen}
        --data {raw,processed}
        --assay {10x,10x_cram,sci,scale}

    For Scale raw validation, a SIF file is strongly recommended:
        --sif PATH

    For 10x processed Multiome (Cell Ranger ARC) mappings, also pass:
        --tenx-profile multiome

    To control how many example issues are printed per check:
        --max-examples N     (default 5; 0 means show all)
        --verbose / -v       (show every error and warning)
    """

    parser = argparse.ArgumentParser(
        description="Validate S3/local mapping files against SOP rules."
    )
    parser.add_argument("--mapping", required=True, help="Path to mapping CSV/TSV file")
    parser.add_argument(
        "--sif", help="Path to SIF CSV/XLSX for completeness checks (Scale/sci)"
    )
    parser.add_argument(
        "--provider",
        required=True,
        choices=["novogene", "psomagen"],
        help="Sequencing provider (affects SOP rules)",
    )
    parser.add_argument(
        "--data",
        required=True,
        choices=["raw", "processed"],
        help="Which kind of data the mapping represents",
    )
    parser.add_argument(
        "--assay",
        required=True,
        choices=["10x", "10x_cram", "sci", "scale"],
        help="High-level assay family for SOP selection",
    )
    parser.add_argument(
        "--tenx-profile",
        choices=["default", "multiome"],
        default="default",
        help=(
            "For --assay 10x --data processed: 'multiome' enforces Cell Ranger ARC-style "
            "core outs (ATAC + feature matrices + summary) per GroupID."
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=(
            "Print every error and warning instead of a capped sample. "
            "Equivalent to --max-examples 0."
        ),
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=DEFAULT_MAX_EXAMPLES,
        metavar="N",
        help=(
            "Maximum number of example errors/warnings to show per check "
            f"(default: {DEFAULT_MAX_EXAMPLES}; 0 means show all). "
            "Ignored when --verbose is set."
        ),
    )

    args = parser.parse_args()

    configure_reporting(verbose=args.verbose, max_examples=args.max_examples)

    mappings = parse_mapping_file(args.mapping, provider=args.provider)
    if not mappings:
        print(
            "No mappings parsed from file (check delimiter and content).",
            file=sys.stderr,
        )
        raise SystemExit(1)

    exit_code = 0
    fail_reasons: List[str] = []

    # 1. Uniqueness
    uniq = validate_uniqueness(mappings)
    dup_local_count = len(uniq["duplicate_locals"])
    dup_s3_count = len(uniq["duplicate_s3"])
    print(
        f"Uniqueness: {uniq['total']} mappings, "
        f"{dup_local_count} duplicate locals, {dup_s3_count} duplicate S3 paths"
    )
    if dup_local_count or dup_s3_count:
        exit_code = 1
        fail_reasons.append(
            f"duplicate mappings ({dup_local_count} local, {dup_s3_count} S3)"
        )

    # 2. Mode-specific validation
    provider = args.provider
    if args.data == "raw" and args.assay == "10x":
        exit_code |= _validate_10x_raw(provider, mappings, args.sif, fail_reasons)

    elif args.data == "raw" and args.assay == "10x_cram":
        exit_code |= _validate_10x_cram_raw(provider, mappings, args.sif, fail_reasons)

    elif (
        args.data == "raw" and args.assay in ("scale", "sci") and provider == "novogene"
    ):
        exit_code |= _validate_seahub_raw(args.assay, mappings, args.sif, fail_reasons)

    elif args.data == "processed" and args.assay == "10x":
        exit_code |= _validate_10x_processed(
            provider,
            mappings,
            args.sif,
            fail_reasons,
            tenx_profile=args.tenx_profile,
        )

    else:
        print(
            f"Mode (provider={provider}, data={args.data}, assay={args.assay}) "
            "is not implemented yet.",
            file=sys.stderr,
        )
        exit_code = 1
        fail_reasons.append("unsupported validation mode")

    # Final verdict
    print()
    if exit_code == 0:
        print("VERDICT: PASS — mapping validates successfully against SOP rules.")
    else:
        print("VERDICT: FAIL — mapping has issues that need attention:")
        for reason in fail_reasons:
            print(f"  - {reason}")

    raise SystemExit(exit_code)


def _validate_10x_raw(
    provider: str,
    mappings: list,
    sif_path: str | None,
    fail_reasons: List[str],
) -> int:
    """Run all 10x raw validation checks, return non-zero on failure."""
    exit_code = 0

    modality = validate_10x_raw_file_modalities(mappings)
    if modality["fastq_count"] == 0:
        print(
            "10x raw modality: 0 FASTQ rows found; 10x mode requires SOP FASTQ files."
        )
        exit_code = 1
        fail_reasons.append("10x raw requires FASTQ rows but none were found")
    if modality["strict_cram_count"]:
        print(
            "10x raw modality: unexpected CRAM rows found "
            f"({modality['strict_cram_count']}; unmatched CRAM artifacts are the only allowed exception)."
        )
        shown_cram = cap(modality["strict_cram_rows"], 5)
        for row in shown_cram:
            print(f"  line {row.line_num}: {row.s3_path}")
        remaining_cram = modality["strict_cram_count"] - len(shown_cram)
        if remaining_cram > 0:
            print(f"  ... and {remaining_cram} more (use --verbose to see all)")
        exit_code = 1
        fail_reasons.append(
            f"{modality['strict_cram_count']} unexpected CRAM rows in 10x FASTQ mode"
        )

    res = validate_s3_10x_raw(provider, mappings)
    meta_count = res.get("metadata_files", 0)
    s3_ga: dict[str, Set[str]] = res.get("group_assays", {})
    s3_run_ids: set[str] = res.get("run_ids", set())
    print(
        f"10x raw SOP: matched {res['matched']} S3 paths, "
        f"{len(res['errors'])} errors, {len(res['warnings'])} warnings"
        + (f", {meta_count} run-metadata files" if meta_count else "")
    )
    _print_issue_examples("10x raw SOP", res["errors"], "errors")
    _print_issue_examples("10x raw SOP", res["warnings"], "warnings")
    parse_miss_count = _count_issue_type(res["warnings"], "parse_miss")
    if res["matched"] == 0:
        print(
            "  WARNING: none of the S3 paths matched the 10x raw pattern. "
            "Check that the S3 paths follow the SOP naming convention."
        )
        fail_reasons.append("no S3 paths matched the 10x raw pattern (0 matched)")
        exit_code = 1
    if parse_miss_count:
        exit_code = 1
        fail_reasons.append(f"{parse_miss_count} S3 paths do not follow 10x raw SOP")
    if res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(res['errors'])} 10x S3 SOP errors")

    if s3_run_ids:
        print(f"Unique RunIDs (wafer identifiers): {len(s3_run_ids)}")
        for rid in sorted(s3_run_ids):
            print(f"  {rid}")

    if s3_ga:
        print(f"S3 GroupIDs found: {len(s3_ga)}")
        for gid in sorted(s3_ga):
            assays_str = ", ".join(sorted(s3_ga[gid]))
            print(f"  {gid}: {assays_str}")

    mate_res = validate_10x_raw_fastq_read_mates(provider, mappings)
    n_skip_illum = mate_res.get("skipped_non_illumina", 0)
    print(
        f"10x raw FASTQ read mates: {mate_res['groups_checked']} Illumina groups, "
        f"{len(mate_res['errors'])} errors, {len(mate_res['warnings'])} warnings"
        + (f", {n_skip_illum} non-Illumina FASTQ tails skipped" if n_skip_illum else "")
    )
    _print_issue_examples("10x FASTQ read mates", mate_res["errors"], "errors")
    _print_issue_examples("10x FASTQ read mates", mate_res["warnings"], "warnings")
    if (
        mate_res["fastq_rows_seen"] > 0
        and mate_res["groups_checked"] == 0
        and not mate_res["errors"]
    ):
        exit_code = 1
        fail_reasons.append(
            "10x FASTQ rows were present but none matched the SOP Illumina tail pattern"
        )
    if mate_res["skipped_non_illumina"]:
        exit_code = 1
        fail_reasons.append(
            f"{mate_res['skipped_non_illumina']} FASTQ rows do not match SOP Illumina naming"
        )
    if mate_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(mate_res['errors'])} FASTQ read-mate errors")

    completeness = validate_10x_raw_fastq_completeness(provider, mappings)
    print(
        f"10x raw FASTQ completeness: matched {completeness['matched']} suffix rows, "
        f"{completeness['sample_prefixes_checked']} prefixes, "
        f"{completeness['tails_checked']} Illumina tails, "
        f"{len(completeness['errors'])} errors, "
        f"{len(completeness['warnings'])} warnings"
    )
    _print_issue_examples(
        "10x raw FASTQ completeness", completeness["errors"], "errors"
    )
    _print_issue_examples(
        "10x raw FASTQ completeness", completeness["warnings"], "warnings"
    )
    if completeness["errors"]:
        exit_code = 1
        fail_reasons.append(
            f"{len(completeness['errors'])} 10x FASTQ completeness errors"
        )

    if sif_path:
        exit_code |= _validate_10x_sif(
            provider, mappings, sif_path, s3_ga, fail_reasons
        )
    else:
        print("10x mode: no --sif provided, skipping SIF completeness checks.")

    return exit_code


def _validate_10x_cram_raw(
    provider: str,
    mappings: list,
    sif_path: str | None,
    fail_reasons: List[str],
) -> int:
    """Run all 10x_cram raw validation checks, return non-zero on failure."""
    exit_code = 0

    modality = validate_10x_raw_file_modalities(mappings)
    if modality["fastq_count"]:
        print(
            "10x_cram modality: FASTQ rows found "
            f"({modality['fastq_count']}); FASTQ files are forbidden in 10x_cram mode."
        )
        shown_fastq = cap(modality["fastq_rows"], 5)
        for row in shown_fastq:
            print(f"  line {row.line_num}: {row.s3_path}")
        remaining_fastq = modality["fastq_count"] - len(shown_fastq)
        if remaining_fastq > 0:
            print(f"  ... and {remaining_fastq} more (use --verbose to see all)")
        exit_code = 1
        fail_reasons.append(
            f"{modality['fastq_count']} FASTQ rows found in 10x_cram mode"
        )

    res = validate_s3_10x_cram_raw(provider, mappings)
    s3_ga: dict[str, Set[str]] = res.get("group_assays", {})
    s3_run_ids: set[str] = res.get("run_ids", set())
    print(
        f"10x_cram raw SOP: matched {res['matched']} sample paths, "
        f"{len(res['errors'])} errors, {len(res['warnings'])} warnings"
        + (
            f", {res.get('metadata_files', 0)} run-metadata files"
            if res.get("metadata_files", 0)
            else ""
        )
    )
    _print_issue_examples("10x_cram raw SOP", res["errors"], "errors")
    _print_issue_examples("10x_cram raw SOP", res["warnings"], "warnings")
    if res["matched"] == 0:
        fail_reasons.append("no S3 paths matched the 10x_cram raw sample pattern")
        exit_code = 1
    if res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(res['errors'])} 10x_cram SOP errors")

    if s3_run_ids:
        print(f"Unique RunIDs (wafer identifiers): {len(s3_run_ids)}")
        for rid in sorted(s3_run_ids):
            print(f"  {rid}")
    if s3_ga:
        print(f"S3 GroupIDs found: {len(s3_ga)}")
        for gid in sorted(s3_ga):
            assays_str = ", ".join(sorted(s3_ga[gid]))
            print(f"  {gid}: {assays_str}")

    if sif_path:
        exit_code |= _validate_10x_sif(
            provider, mappings, sif_path, s3_ga, fail_reasons
        )
    else:
        print("10x_cram mode: no --sif provided, skipping SIF completeness checks.")

    return exit_code


def _validate_10x_sif(
    provider: str,
    mappings: list,
    sif_path: str,
    s3_ga: dict[str, Set[str]],
    fail_reasons: List[str],
) -> int:
    """Run SIF-based checks for 10x, return non-zero on failure."""
    exit_code = 0

    sif_ga = load_sif_group_assays(sif_path, provider=provider)
    sif_norm: dict[str, set[str]] = {
        _normalize_sif_groupid(k): v for k, v in sif_ga.items()
    }

    # Build actual group_assays with lower-cased assay names for comparison
    actual_ga: dict[str, set[str]] = {
        gid: {a.lower() for a in assays} for gid, assays in s3_ga.items()
    }

    cmp = compare_groupid_assays(sif_norm, actual_ga)
    exit_code |= _report_sif_comparison(cmp, fail_reasons)

    # Library-level cross-checks (assay + GroupID consistency)
    lib_assays = load_sif_library_assays(sif_path, provider=provider)
    if lib_assays:
        exit_code |= _report_library_consistency(
            mappings, lib_assays, provider, fail_reasons
        )

    # Per-path SIF coverage
    sif_ids = set(sif_norm.keys())
    exit_code |= _report_sif_path_coverage(mappings, sif_ids, provider, fail_reasons)

    return exit_code


def _report_library_consistency(
    mappings: list,
    lib_assays: dict[str, str],
    provider: str,
    fail_reasons: List[str],
) -> int:
    """Run and report library-assay consistency checks."""
    exit_code = 0
    lib_res = validate_library_assay_consistency(mappings, lib_assays, provider)
    n_assay = len(lib_res["assay_mismatches"])
    n_gid = len(lib_res["groupid_mismatches"])
    print(
        f"Library consistency: checked {lib_res['checked']} paths, "
        f"{n_assay} assay mismatches, {n_gid} GroupID mismatches"
        + (
            f", {lib_res['skipped']} skipped (no library name found)"
            if lib_res["skipped"]
            else ""
        )
    )
    if n_gid:
        print("  GroupID mismatches (library in local path not in S3 GroupID):")
        shown_gid = cap(lib_res["groupid_mismatches"], 10)
        for item in shown_gid:
            print(
                f"    line {item['line']}: local has '{item['library']}' "
                f"but S3 GroupID is '{item['s3_groupid']}'"
            )
        if n_gid > len(shown_gid):
            print(
                f"    ... and {n_gid - len(shown_gid)} more (use --verbose to see all)"
            )
        exit_code = 1
        fail_reasons.append(
            f"{n_gid} GroupID mismatches (local library name not found in S3 GroupID)"
        )
    if n_assay:
        print("  Assay mismatches (library in local path vs assay in S3):")
        shown_assay = cap(lib_res["assay_mismatches"], 10)
        for item in shown_assay:
            print(
                f"    line {item['line']}: local has '{item['library']}' "
                f"(SIF expects {item['expected_assay']}) "
                f"but S3 assay is {item['s3_assay']}"
            )
        if n_assay > len(shown_assay):
            print(
                f"    ... and {n_assay - len(shown_assay)} more (use --verbose to see all)"
            )
        exit_code = 1
        fail_reasons.append(
            f"{n_assay} library-assay mismatches "
            "(local library name doesn't match S3 assay per SIF)"
        )
    return exit_code


def _report_sif_path_coverage(
    mappings: list,
    sif_ids: set[str],
    provider: str,
    fail_reasons: List[str],
) -> int:
    """Run and report per-path SIF coverage for 10x."""
    exit_code = 0
    sif_cov = find_unmatched_sif_paths_10x(mappings, sif_ids, provider)
    n_unmatched_groups = len(sif_cov["unmatched_by_group"])
    n_unmatched_paths = sum(
        len(rows) for rows in sif_cov["unmatched_by_group"].values()
    )
    n_unparsed = len(sif_cov["unparsed"])
    print(
        f"SIF path coverage: {sif_cov['matched_sif']} paths matched SIF, "
        f"{n_unmatched_paths} paths in {n_unmatched_groups} extra GroupIDs, "
        f"{n_unparsed} unparsed"
        + (f", {sif_cov['metadata']} metadata skipped" if sif_cov["metadata"] else "")
    )
    if sif_cov["unmatched_by_group"]:
        print("  S3 paths with GroupIDs not in SIF:")
        for gid in sorted(sif_cov["unmatched_by_group"]):
            rows = sif_cov["unmatched_by_group"][gid]
            print(f"    {gid}: {len(rows)} files")
            shown_rows = cap(rows, 3)
            for r in shown_rows:
                print(f"      {r.s3_path}")
            if len(rows) > len(shown_rows):
                print(
                    f"      ... and {len(rows) - len(shown_rows)} more "
                    "(use --verbose to see all)"
                )
        exit_code = 1
        fail_reasons.append(
            f"{n_unmatched_paths} S3 paths in {n_unmatched_groups} "
            "GroupIDs not found in SIF"
        )
    if sif_cov["unparsed"]:
        print("  S3 paths that could not be parsed (not metadata):")
        shown_unparsed = cap(sif_cov["unparsed"], 5)
        for r in shown_unparsed:
            print(f"    line {r.line_num}: {r.s3_path}")
        if n_unparsed > len(shown_unparsed):
            print(
                f"    ... and {n_unparsed - len(shown_unparsed)} more "
                "(use --verbose to see all)"
            )
        exit_code = 1
        fail_reasons.append(f"{n_unparsed} S3 paths could not be parsed")
    return exit_code


_SEAHUB_LOCAL_VALIDATORS: dict = {
    "scale": validate_local_paths_scale_raw,
    "sci": validate_local_paths_sci_raw,
}

_SEAHUB_CONSISTENCY_VALIDATORS: dict = {
    "scale": validate_s3_local_consistency_scale,
    "sci": validate_s3_local_consistency_sci,
}


_SEAHUB_COMPLETENESS_VALIDATORS: dict = {
    "scale": validate_scale_raw_completeness,
    "sci": validate_sci_raw_completeness,
}


def _validate_seahub_raw(
    assay_family: str,
    mappings: list,
    sif_path: str | None,
    fail_reasons: List[str],
) -> int:
    """Run all seahub (Scale or sci) raw validation checks."""
    family = assay_family.lower()
    exit_code = 0

    # S3 SOP checks
    s3_res = validate_s3_seahub_raw(family, mappings)
    meta_count = s3_res.get("metadata_files", 0)
    s3_ga: dict[str, Set[str]] = s3_res.get("group_assays", {})
    s3_run_ids: set[str] = s3_res.get("run_ids", set())
    print(
        f"{family} raw (S3): matched {s3_res['matched']} paths, "
        f"{len(s3_res['errors'])} errors, {len(s3_res['warnings'])} warnings"
        + (f", {meta_count} run-metadata files" if meta_count else "")
    )
    _print_issue_examples(f"{family} raw (S3)", s3_res["errors"], "errors")
    _print_issue_examples(f"{family} raw (S3)", s3_res["warnings"], "warnings")
    if s3_res["matched"] == 0:
        print(
            f"  WARNING: none of the S3 paths matched the {family} raw pattern. "
            "Check that the S3 paths follow the SOP naming convention."
        )
        fail_reasons.append(f"no S3 paths matched the {family} raw pattern (0 matched)")
        exit_code = 1
    if s3_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(s3_res['errors'])} {family} S3 SOP errors")

    if s3_run_ids:
        print(f"Unique RunIDs (wafer identifiers): {len(s3_run_ids)}")
        for rid in sorted(s3_run_ids):
            print(f"  {rid}")

    if s3_ga:
        print(f"S3 GroupIDs found: {len(s3_ga)}")
        for gid in sorted(s3_ga):
            assays_str = ", ".join(sorted(s3_ga[gid]))
            print(f"  {gid}: {assays_str}")

    # Local path sanity (family-specific)
    local_validator = _SEAHUB_LOCAL_VALIDATORS[family]
    local_res = local_validator(mappings)
    print(
        f"{family} raw (local): matched {local_res['matched']} paths, "
        f"{len(local_res['errors'])} errors, {len(local_res['warnings'])} warnings"
    )
    _print_issue_examples(f"{family} raw (local)", local_res["errors"], "errors")
    _print_issue_examples(f"{family} raw (local)", local_res["warnings"], "warnings")
    if local_res["matched"] == 0:
        print(
            f"  WARNING: none of the local paths matched any recognised {family} pattern. "
            "Local-path consistency checks were NOT applied."
        )
        fail_reasons.append(f"no local paths matched {family} pattern (0 matched)")
        exit_code = 1
    if local_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(local_res['errors'])} local-path errors")

    # SIF completeness
    if sif_path:
        sif_res = validate_sif_completeness_seahub(family, mappings, sif_path)
        exit_code |= _report_sif_comparison(sif_res, fail_reasons)
    else:
        print(f"{family} mode: no --sif provided, skipping SIF completeness checks.")

    # S3/local consistency (family-specific)
    consistency_validator = _SEAHUB_CONSISTENCY_VALIDATORS[family]
    cons_res = consistency_validator(mappings)
    print(
        f"{family} S3/local consistency: matched {cons_res['matched']} pairs, "
        f"{len(cons_res['errors'])} errors, {len(cons_res['warnings'])} warnings"
    )
    _print_issue_examples(
        f"{family} S3/local consistency", cons_res["errors"], "errors"
    )
    _print_issue_examples(
        f"{family} S3/local consistency", cons_res["warnings"], "warnings"
    )
    if cons_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(cons_res['errors'])} S3/local consistency errors")

    # SOP-driven per-prefix completeness (Scale: per-well + per-UG; sci: per-prefix)
    completeness_validator = _SEAHUB_COMPLETENESS_VALIDATORS[family]
    comp_res = completeness_validator(mappings)
    if family == "scale":
        print(
            f"{family} raw completeness: matched {comp_res['matched']} suffix rows, "
            f"{comp_res.get('well_prefixes_checked', 0)} per-well bundles, "
            f"{comp_res.get('ug_aggregates_checked', 0)} per-UG aggregates, "
            f"{len(comp_res['errors'])} errors, "
            f"{len(comp_res['warnings'])} warnings"
        )
    else:
        print(
            f"{family} raw completeness: matched {comp_res['matched']} suffix rows, "
            f"{comp_res.get('sample_prefixes_checked', 0)} prefixes, "
            f"{len(comp_res['errors'])} errors, "
            f"{len(comp_res['warnings'])} warnings"
        )
    _print_issue_examples(f"{family} raw completeness", comp_res["errors"], "errors")
    _print_issue_examples(
        f"{family} raw completeness", comp_res["warnings"], "warnings"
    )
    if comp_res["errors"]:
        exit_code = 1
        fail_reasons.append(
            f"{len(comp_res['errors'])} {family} raw completeness errors"
        )

    return exit_code


# ---------------------------------------------------------------------------
# 10x processed orchestrator
# ---------------------------------------------------------------------------


def _validate_10x_processed(
    provider: str,
    mappings: list,
    sif_path: str | None,
    fail_reasons: List[str],
    *,
    tenx_profile: str = "default",
) -> int:
    """Run all 10x processed validation checks, return non-zero on failure."""
    exit_code = 0

    # S3 SOP checks
    s3_res = validate_s3_10x_processed(provider, mappings)
    group_ids: set[str] = s3_res.get("group_ids", set())
    pipelines: set[str] = s3_res.get("pipelines", set())
    run_dates: set[str] = s3_res.get("run_dates", set())
    print(
        f"10x processed (S3): matched {s3_res['matched']} paths, "
        f"{len(s3_res['errors'])} errors, {len(s3_res['warnings'])} warnings"
    )
    _print_issue_examples("10x processed (S3)", s3_res["errors"], "errors")
    _print_issue_examples("10x processed (S3)", s3_res["warnings"], "warnings")
    if s3_res["matched"] == 0:
        print(
            "  WARNING: none of the S3 paths matched the 10x processed pattern. "
            "Check that the S3 paths follow the SOP naming convention."
        )
        fail_reasons.append("no S3 paths matched the 10x processed pattern (0 matched)")
        exit_code = 1
    if s3_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(s3_res['errors'])} 10x processed S3 SOP errors")

    if pipelines:
        print(f"Pipelines found: {', '.join(sorted(pipelines))}")
    if run_dates:
        print(f"Run dates found: {', '.join(sorted(run_dates))}")
    if group_ids:
        print(f"Unique GroupIDs found: {len(group_ids)}")
        for gid in sorted(group_ids):
            print(f"  {gid}")

    # SIF completeness
    if sif_path:
        sif_res = validate_sif_completeness_10x_processed(provider, mappings, sif_path)
        sif_count = sif_res["sif_count"]
        s3_count = sif_res["s3_count"]
        missing = sif_res["missing"]
        extra = sif_res["extra"]
        print(
            f"SIF completeness: expected={sif_count} identifiers, "
            f"found={s3_count} in S3, "
            f"missing={len(missing)}, extra={len(extra)}"
        )
        if missing:
            print(f"  Missing from S3 (in SIF but not mapping): {', '.join(missing)}")
            exit_code = 1
            fail_reasons.append(
                f"{len(missing)} SIF identifiers missing from S3: {', '.join(missing)}"
            )
        if extra:
            print(f"  Extra in S3 (in mapping but not SIF): {', '.join(extra)}")
    else:
        print(
            "10x processed mode: no --sif provided, skipping SIF completeness checks."
        )

    # S3/local consistency
    cons_res = validate_s3_local_consistency_10x_processed(provider, mappings)
    normalized_groupid_count = _count_issue_type(
        cons_res["warnings"], "group_id_normalized_match"
    )
    print(
        f"10x processed S3/local consistency: matched {cons_res['matched']} pairs, "
        f"{len(cons_res['errors'])} errors, {len(cons_res['warnings'])} warnings"
    )
    if normalized_groupid_count:
        print(
            "  NOTE: "
            f"{normalized_groupid_count} GroupID warning(s) matched after '-'/'_' "
            "normalization."
        )
    _print_issue_examples(
        "10x processed S3/local consistency", cons_res["errors"], "errors"
    )
    _print_issue_examples(
        "10x processed S3/local consistency", cons_res["warnings"], "warnings"
    )
    if cons_res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(cons_res['errors'])} S3/local consistency errors")

    if tenx_profile == "multiome":
        mo_res = validate_10x_multiome_processed_outs(provider, mappings)
        n_miss = len(mo_res["missing_by_group"])
        skip = mo_res["skipped_unparsed"]
        print(
            f"10x multiome processed outs: checked {mo_res['groups_checked']} GroupIDs, "
            f"{n_miss} with missing required files"
            + (f", {skip} S3 rows skipped (not processed paths)" if skip else "")
        )
        for err in mo_res["errors"]:
            print(f"  {err['group_id']}: {err['detail']}")
        if mo_res["errors"]:
            exit_code = 1
            fail_reasons.append(
                f"{n_miss} GroupIDs missing required multiome processed outs"
            )

    return exit_code


__all__ = ["main"]
