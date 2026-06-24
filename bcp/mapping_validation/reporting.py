"""
Console reporting helpers for the mapping_validation CLI.

These functions render validator result dicts (``errors``/``warnings`` lists
and SIF comparison summaries) to stdout.  Output verbosity is controlled by a
small module-level configuration set once from the CLI via
:func:`configure_reporting`:

- default: each grouped example list is capped at ``max_examples`` (default 5)
  and the per-group detail lists keep their original short caps;
- verbose: every error/warning and every per-group detail is printed.

Keeping this logic in one place lets the CLI orchestrators stay focused on
*which* validators to run, while presentation rules (caps, ``--verbose``)
live here and apply uniformly.
"""

from __future__ import annotations

from collections import defaultdict
from typing import List, Sequence, TypeVar

from .constants import CANONICAL_ASSAY

_T = TypeVar("_T")

DEFAULT_MAX_EXAMPLES = 5

# Module-level reporting configuration, set once by the CLI.
_VERBOSE = False
_MAX_EXAMPLES = DEFAULT_MAX_EXAMPLES


def configure_reporting(
    *, verbose: bool = False, max_examples: int = DEFAULT_MAX_EXAMPLES
) -> None:
    """Set the global reporting verbosity used by the print helpers.

    ``verbose`` overrides ``max_examples`` and prints every item.  A
    ``max_examples`` of 0 or less is also treated as unlimited.
    """
    global _VERBOSE, _MAX_EXAMPLES
    _VERBOSE = verbose
    _MAX_EXAMPLES = max_examples


def is_verbose() -> bool:
    """Return whether unlimited (verbose) output is currently enabled."""
    return _VERBOSE or _MAX_EXAMPLES <= 0


def cap(items: Sequence[_T], limit: int) -> Sequence[_T]:
    """Truncate ``items`` to ``limit`` unless verbose mode is active.

    A ``limit`` of 0 or less means "no limit".  Used by the CLI orchestrators
    for their per-group detail lists so ``--verbose`` lifts every short cap.
    """
    if _VERBOSE or limit <= 0:
        return items
    return items[:limit]


def print_issue_examples(label: str, issues: List[dict], kind: str) -> None:
    """Print a grouped summary of errors or warnings with examples.

    ``kind`` is either ``'errors'`` or ``'warnings'`` for labeling.  By
    default only the first ``max_examples`` items are shown with a trailing
    "... and N more" hint; ``--verbose`` prints every item.
    """
    if not issues:
        return

    by_type: dict[str, int] = defaultdict(int)
    for item in issues:
        by_type[item.get("type", "unknown")] += 1

    print(f"{label} {kind}:")
    for t, count in sorted(by_type.items(), key=lambda kv: kv[0]):
        print(f"  - {t}: {count}")

    shown = cap(issues, _MAX_EXAMPLES)
    if len(shown) == len(issues):
        print(f"  Showing all {len(issues)} {kind}:")
    else:
        print(f"  Showing up to {_MAX_EXAMPLES} example {kind}:")

    for item in shown:
        _print_issue(item)

    remaining = len(issues) - len(shown)
    if remaining > 0:
        print(f"    ... and {remaining} more {kind} (use --verbose to see all)")


def _print_issue(item: dict) -> None:
    """Render a single error/warning dict as indented detail lines."""
    line = item.get("line", "?")
    s3_path = item.get("s3_path")
    local_path = item.get("local_path")
    detail = item.get("detail")
    present = item.get("present")
    source_lines = item.get("lines")
    missing_files = item.get("missing_files")

    print(f"    line {line}:")
    if s3_path:
        print(f"      S3:    {s3_path}")
    if local_path:
        print(f"      local: {local_path}")
    if detail:
        print(f"      detail: {detail}")
    if present:
        print(f"      present: {', '.join(present)}")
    if missing_files:
        print("      missing files (SOP-expected, not found in mapping):")
        for name in missing_files:
            print(f"        {name}")
    if source_lines:
        print(f"      source lines: {', '.join(str(n) for n in source_lines)}")


def count_issue_type(issues: List[dict], issue_type: str) -> int:
    """Count issues of a specific type."""
    return sum(1 for item in issues if item.get("type") == issue_type)


def report_sif_comparison(cmp: dict, fail_reasons: List[str]) -> int:
    """Print SIF GroupID/assay comparison results and return exit-code delta.

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
        print("  Extra GroupIDs in S3 (not in SIF): " + ", ".join(sorted(extra)))
        exit_code = 1
        fail_reasons.append(
            f"{len(extra)} extra GroupIDs in S3 not in SIF: " + ", ".join(sorted(extra))
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


__all__ = [
    "DEFAULT_MAX_EXAMPLES",
    "configure_reporting",
    "is_verbose",
    "cap",
    "print_issue_examples",
    "count_issue_type",
    "report_sif_comparison",
]
