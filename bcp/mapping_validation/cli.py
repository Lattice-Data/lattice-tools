from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from typing import List, Set

from .constants import CANONICAL_ASSAY
from .parsing import parse_mapping_file
from .sif_io import (
    _normalize_sif_groupid,
    load_sif_group_assays,
    load_sif_library_assays,
)
from .uniqueness import validate_uniqueness
from .validators import (
    compare_groupid_assays,
    find_unmatched_sif_paths_10x,
    validate_library_assay_consistency,
    validate_local_paths_scale_raw,
    validate_local_paths_sci_raw,
    validate_s3_10x_raw,
    validate_s3_local_consistency_scale,
    validate_s3_local_consistency_sci,
    validate_s3_seahub_raw,
    validate_sif_completeness_seahub,
)


def _print_issue_examples(label: str, issues: List[dict], kind: str, max_examples: int = 5) -> None:
    """
    Print a short, grouped summary of errors or warnings with a few examples.

    `kind` is either 'errors' or 'warnings' for labeling.
    """
    if not issues:
        return

    by_type: dict[str, int] = defaultdict(int)
    for item in issues:
        by_type[item.get("type", "unknown")] += 1

    print(f"{label} {kind}:")
    for t, count in sorted(by_type.items(), key=lambda kv: kv[0]):
        print(f"  - {t}: {count}")

    print(f"  Showing up to {max_examples} example {kind}:")
    for item in issues[: max_examples]:
        line = item.get("line", "?")
        s3_path = item.get("s3_path")
        local_path = item.get("local_path")
        detail = item.get("detail")
        print(f"    line {line}:")
        if s3_path:
            print(f"      S3:    {s3_path}")
        if local_path:
            print(f"      local: {local_path}")
        if detail:
            print(f"      detail: {detail}")


def _report_sif_comparison(
    cmp: dict,
    fail_reasons: List[str],
) -> int:
    """Print SIF GroupID/assay comparison results and return exit code delta.

    Shared reporting logic for both 10x and Scale SIF completeness checks.
    Returns 1 if any issues warrant failure, 0 otherwise.
    """
    exit_code = 0
    missing = cmp["missing_groupids"]
    extra = cmp["extra_groupids"]
    missing_assays = cmp.get("missing_assays", {})
    extra_assays = cmp.get("extra_assays", {})
    expected_ga = cmp.get("expected_group_assays", {})

    print(
        f"SIF completeness: expected={len(cmp['expected_groupids'])} GroupIDs, "
        f"found={len(cmp['actual_groupids'])}, "
        f"missing={len(missing)}, extra={len(extra)}"
    )

    if missing:
        print(
            "  Missing GroupIDs from S3 (in SIF but not S3): "
            + ", ".join(sorted(missing))
        )
        exit_code = 1
        fail_reasons.append(
            f"{len(missing)} SIF GroupIDs missing from S3: "
            + ", ".join(sorted(missing))
        )
    if extra:
        print(
            "  Extra GroupIDs in S3 (not in SIF): " + ", ".join(sorted(extra))
        )
        exit_code = 1
        fail_reasons.append(
            f"{len(extra)} extra GroupIDs in S3 not in SIF: "
            + ", ".join(sorted(extra))
        )

    if missing_assays:
        print("  GroupIDs with missing assay types in S3:")
        for gid in sorted(missing_assays):
            print(f"    {gid}: missing {sorted(missing_assays[gid])}")
        exit_code = 1
        fail_reasons.append(
            f"{len(missing_assays)} GroupIDs have missing assay types in S3"
        )
    if extra_assays:
        print("  GroupIDs with unexpected assay types in S3 (not in SIF):")
        for gid in sorted(extra_assays):
            print(f"    {gid}: unexpected {sorted(extra_assays[gid])}")

    if expected_ga:
        unique_combos: set[frozenset[str]] = set()
        for assays in expected_ga.values():
            unique_combos.add(frozenset(assays))
        combos_desc = [
            " + ".join(CANONICAL_ASSAY.get(a, a) for a in sorted(c))
            for c in sorted(unique_combos, key=sorted)
        ]
        print(f"  SIF assay combinations: {', '.join(combos_desc)}")

    return exit_code


def main() -> None:
    """
    CLI for validating mapping CSV/TSV files against SOP rules.

    Required arguments:
        --mapping PATH       Mapping file with two columns: S3, local
        --provider {novogene,psomagen}
        --data {raw,processed}
        --assay {10x,sci,scale}

    For Scale raw validation, a SIF file is strongly recommended:
        --sif PATH
    """

    parser = argparse.ArgumentParser(description="Validate S3/local mapping files against SOP rules.")
    parser.add_argument("--mapping", required=True, help="Path to mapping CSV/TSV file")
    parser.add_argument("--sif", help="Path to SIF CSV/XLSX for completeness checks (Scale/sci)")
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
        choices=["10x", "sci", "scale"],
        help="High-level assay family for SOP selection",
    )

    args = parser.parse_args()

    mappings = parse_mapping_file(args.mapping)
    if not mappings:
        print("No mappings parsed from file (check delimiter and content).", file=sys.stderr)
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
        fail_reasons.append(f"duplicate mappings ({dup_local_count} local, {dup_s3_count} S3)")

    # 2. Mode-specific validation
    provider = args.provider
    if args.data == "raw" and args.assay == "10x":
        exit_code |= _validate_10x_raw(provider, mappings, args.sif, fail_reasons)

    elif args.data == "raw" and args.assay in ("scale", "sci") and provider == "novogene":
        exit_code |= _validate_seahub_raw(args.assay, mappings, args.sif, fail_reasons)

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

    res = validate_s3_10x_raw(provider, mappings)
    meta_count = res.get("metadata_files", 0)
    s3_ga: dict[str, Set[str]] = res.get("group_assays", {})
    print(
        f"10x raw SOP: matched {res['matched']} S3 paths, "
        f"{len(res['errors'])} errors, {len(res['warnings'])} warnings"
        + (f", {meta_count} run-metadata files" if meta_count else "")
    )
    _print_issue_examples("10x raw SOP", res["errors"], "errors")
    _print_issue_examples("10x raw SOP", res["warnings"], "warnings")
    if res["matched"] == 0:
        print(
            "  WARNING: none of the S3 paths matched the 10x raw pattern. "
            "Check that the S3 paths follow the SOP naming convention."
        )
        fail_reasons.append("no S3 paths matched the 10x raw pattern (0 matched)")
        exit_code = 1
    if res["errors"]:
        exit_code = 1
        fail_reasons.append(f"{len(res['errors'])} 10x S3 SOP errors")

    if s3_ga:
        print(f"S3 GroupIDs found: {len(s3_ga)}")
        for gid in sorted(s3_ga):
            assays_str = ", ".join(sorted(s3_ga[gid]))
            print(f"  {gid}: {assays_str}")

    if sif_path:
        exit_code |= _validate_10x_sif(provider, mappings, sif_path, s3_ga, fail_reasons)
    else:
        print("10x mode: no --sif provided, skipping SIF completeness checks.")

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

    sif_ga = load_sif_group_assays(sif_path)
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
    lib_assays = load_sif_library_assays(sif_path)
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
        for item in lib_res["groupid_mismatches"][:10]:
            print(
                f"    line {item['line']}: local has '{item['library']}' "
                f"but S3 GroupID is '{item['s3_groupid']}'"
            )
        if n_gid > 10:
            print(f"    ... and {n_gid - 10} more")
        exit_code = 1
        fail_reasons.append(
            f"{n_gid} GroupID mismatches "
            "(local library name not found in S3 GroupID)"
        )
    if n_assay:
        print("  Assay mismatches (library in local path vs assay in S3):")
        for item in lib_res["assay_mismatches"][:10]:
            print(
                f"    line {item['line']}: local has '{item['library']}' "
                f"(SIF expects {item['expected_assay']}) "
                f"but S3 assay is {item['s3_assay']}"
            )
        if n_assay > 10:
            print(f"    ... and {n_assay - 10} more")
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
    n_unmatched_paths = sum(len(rows) for rows in sif_cov["unmatched_by_group"].values())
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
            for r in rows[:3]:
                print(f"      {r.s3_path}")
            if len(rows) > 3:
                print(f"      ... and {len(rows) - 3} more")
        exit_code = 1
        fail_reasons.append(
            f"{n_unmatched_paths} S3 paths in {n_unmatched_groups} "
            "GroupIDs not found in SIF"
        )
    if sif_cov["unparsed"]:
        print("  S3 paths that could not be parsed (not metadata):")
        for r in sif_cov["unparsed"][:5]:
            print(f"    line {r.line_num}: {r.s3_path}")
        if n_unparsed > 5:
            print(f"    ... and {n_unparsed - 5} more")
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

    return exit_code


__all__ = ["main"]
