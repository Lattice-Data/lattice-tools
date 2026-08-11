"""CSV batch I/O and single-row output for structure cross-checking."""

from __future__ import annotations

import csv
import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Literal, NamedTuple

from .client import (
    OUTPUT_FIELDS_APPENDED,
    UNRESOLVED_VERDICTS as _UNRESOLVED,
    REVIEW_CHECK,
    REVIEW_INVESTIGATE,
    REVIEW_OK,
    REVIEW_UNVERIFIED,
    UNRESOLVED_VERDICTS,
    check_row,
)

log = logging.getLogger(__name__)

# Reported in this order, worst first, so the tail of the log is the to-do list.
REVIEW_LEVELS = (REVIEW_INVESTIGATE, REVIEW_CHECK, REVIEW_UNVERIFIED, REVIEW_OK)

# Below this many rows, "nothing resolved" is as likely to be a two-row spot check
# as a broken run, so the degraded verdict stays quiet.
MIN_ROWS_FOR_DEGRADED = 5

SINGLE_INPUT_FIELDS = ["name", "cas", "chebi_id"]
SINGLE_OUTPUT_FIELDS = [*SINGLE_INPUT_FIELDS, *OUTPUT_FIELDS_APPENDED]


class StructureCheckError(Exception):
    """Raised when input file or column validation fails."""


class RunSummary(NamedTuple):
    """Outcome of a batch run."""

    review_counts: Counter[str]
    verdict_counts: Counter[str]
    rows: int
    # Per-side, because a check can fail run-wide while the other side sails through.
    id_unresolved_rows: int = 0
    name_unresolved_rows: int = 0
    id_requested: bool = False
    name_requested: bool = False

    @property
    def needs_attention(self) -> int:
        return self.review_counts[REVIEW_INVESTIGATE] + self.review_counts[REVIEW_CHECK]

    @property
    def dead_checks(self) -> tuple[str, ...]:
        """
        Requested checks that came back unresolved on every single row.

        A per-row verdict is honest about this, but the run-level verdict was not:
        if ChEBI is unreachable while the name column resolves fine, every row gets
        chebi_unresolved *and* a name match, review_level says `ok`, and the sheet
        reads as clean — with the ID question, the one a curator actually asked,
        never answered on a single row.
        """
        dead = []
        if self.id_requested and self.id_unresolved_rows == self.rows:
            dead.append("ChEBI ID vs CAS")
        if self.name_requested and self.name_unresolved_rows == self.rows:
            dead.append("name vs CAS")
        return tuple(dead)

    @property
    def degraded(self) -> bool:
        """
        True when the run cannot be read as a verdict on the sheet.

        Either nothing was compared at all, or a check that was asked for never
        succeeded once. Both produce a complete-looking CSV whose emptiness means
        nothing about the sheet — the single most misleading thing this tool could
        produce.
        """
        if self.rows < MIN_ROWS_FOR_DEGRADED:
            return False
        if self.review_counts[REVIEW_UNVERIFIED] == self.rows:
            return True
        return bool(self.dead_checks)


def check_output_path(output_path: Path | None, input_path: Path | None = None) -> None:
    """
    Fail on an unusable output before spending any requests.

    Refuses to write over the input: the whole file is read before the output is
    opened, so in-place enrichment looks like it works right up until a crash
    leaves a truncated spreadsheet.
    """
    if output_path is None:
        return
    if output_path.is_dir():
        raise StructureCheckError(f"Output path is a directory: {output_path}")
    if not output_path.parent.is_dir():
        raise StructureCheckError(
            f"Output directory does not exist: {output_path.parent}"
        )
    if input_path is not None and output_path.resolve() == input_path.resolve():
        raise StructureCheckError(
            f"Output path is the input file, which would be overwritten: {output_path}"
        )


def emit_single_row(
    output_path: Path | None,
    *,
    name: str | None = None,
    cas: str | None = None,
    chebi_id: str | None = None,
    fmt: Literal["json", "csv"] = "json",
) -> dict[str, Any]:
    """Cross-check one row and write JSON or CSV to a file or stdout."""
    check_output_path(output_path)
    result = check_row(cas=cas, chebi_id=chebi_id, name=name)

    row: dict[str, Any] = {
        "name": (name or "").strip(),
        "cas": (cas or "").strip(),
        "chebi_id": (chebi_id or "").strip(),
    }
    row.update({field: result.get(field, "") for field in OUTPUT_FIELDS_APPENDED})

    log.info("review: %s", result["review"])
    log.info("  ChEBI ID vs CAS : %s", result["id_cas_verdict"])
    log.info("  name vs CAS     : %s", result["name_cas_verdict"])
    if result["name_query"]:
        log.info("  name resolved as: %r", result["name_query"])

    if fmt == "json":
        payload = json.dumps(row, indent=2, ensure_ascii=False)
        if output_path is None:
            print(payload)
        else:
            output_path.write_text(payload + "\n", encoding="utf-8")
            log.info("Output: %s", output_path)
        return result

    out_fh = (
        sys.stdout
        if output_path is None
        else open(output_path, "w", newline="", encoding="utf-8")
    )
    try:
        writer = csv.DictWriter(out_fh, fieldnames=SINGLE_OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerow({k: "" if v is None else str(v) for k, v in row.items()})
    finally:
        if output_path is not None:
            out_fh.close()
            log.info("Output: %s", output_path)
    return result


def _resolve_columns(
    fieldnames: list[str],
    cas_column: str,
    chebi_column: str | None,
    name_column: str | None,
) -> list[str]:
    """Validate the requested columns and return the input columns to carry through."""
    if not chebi_column and not name_column:
        raise StructureCheckError(
            "Nothing to check against the CAS column: pass --chebi-column, "
            "--name-column, or both."
        )
    for label, column in (
        ("--chebi-column", chebi_column),
        ("--name-column", name_column),
    ):
        if column and column == cas_column:
            # Resolving the same string through two PubChem endpoints agrees with
            # itself, so every row would report match/ok — a false all-clear.
            raise StructureCheckError(
                f"Column '{column}' ({label}) is the same as --cas-column, so every "
                f"row would be compared against itself."
            )

    requested = [("--cas-column", cas_column)]
    requested += [("--chebi-column", chebi_column), ("--name-column", name_column)]
    for label, column in requested:
        if not column:
            continue
        if column not in fieldnames:
            raise StructureCheckError(
                f"Column '{column}' ({label}) not found. Available: {fieldnames}"
            )
        if column in OUTPUT_FIELDS_APPENDED:
            raise StructureCheckError(
                f"Column '{column}' ({label}) collides with a column this tool "
                f"writes. Rename it in the input. Reserved: {OUTPUT_FIELDS_APPENDED}"
            )

    clobbered = [c for c in fieldnames if c in OUTPUT_FIELDS_APPENDED]
    if clobbered:
        log.warning("Input columns overwritten: %s", ", ".join(clobbered))
    return [c for c in fieldnames if c not in OUTPUT_FIELDS_APPENDED]


def check_file(
    input_path: Path,
    output_path: Path,
    *,
    cas_column: str,
    chebi_column: str | None = None,
    name_column: str | None = None,
) -> RunSummary:
    """Cross-check every row of a CSV and write the enriched result."""
    if not input_path.is_file():
        raise StructureCheckError(f"Input file not found (or not a file): {input_path}")
    check_output_path(output_path, input_path)

    with open(input_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        carried = _resolve_columns(fieldnames, cas_column, chebi_column, name_column)
        rows = list(reader)

    output_fields = carried + OUTPUT_FIELDS_APPENDED
    total = len(rows)
    # One lookup per distinct triple: the same compound often repeats across rows.
    cache: dict[tuple[str, str, str], dict[str, Any]] = {}
    review_counts: Counter[str] = Counter()
    verdict_counts: Counter[str] = Counter()
    id_unresolved_rows = 0
    name_unresolved_rows = 0

    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=output_fields)
        writer.writeheader()

        for i, row in enumerate(rows, 1):
            cas = (row.get(cas_column) or "").strip()
            chebi_id = (row.get(chebi_column) or "").strip() if chebi_column else ""
            name = (row.get(name_column) or "").strip() if name_column else ""

            key = (cas, chebi_id, name)
            if key not in cache:
                cache[key] = check_row(cas=cas, chebi_id=chebi_id, name=name)
            result = cache[key]

            out_row = {c: row.get(c, "") for c in carried}
            out_row.update({f: result.get(f, "") for f in OUTPUT_FIELDS_APPENDED})
            writer.writerow(out_row)

            review = result["review"]
            review_counts[review] += 1
            for field in ("id_cas_verdict", "name_cas_verdict"):
                verdict_counts[result[field]] += 1
            if result["id_cas_verdict"] in _UNRESOLVED:
                id_unresolved_rows += 1
            if result["name_cas_verdict"] in _UNRESOLVED:
                name_unresolved_rows += 1

            label = name or chebi_id or cas or "(blank row)"
            if review in (REVIEW_INVESTIGATE, REVIEW_CHECK):
                log.warning(
                    "[%s/%s] %s — %s: id_cas=%s name_cas=%s",
                    i,
                    total,
                    label,
                    review.upper() if review == REVIEW_INVESTIGATE else review,
                    result["id_cas_verdict"],
                    result["name_cas_verdict"],
                )
            else:
                log.info("[%s/%s] %s — %s", i, total, label, review)

    summary = RunSummary(
        review_counts=review_counts,
        verdict_counts=verdict_counts,
        rows=total,
        id_unresolved_rows=id_unresolved_rows,
        name_unresolved_rows=name_unresolved_rows,
        id_requested=bool(chebi_column),
        name_requested=bool(name_column),
    )

    log.info("=" * 58)
    log.info("Done. %s rows (%s distinct lookups).", total, len(cache))
    for level in REVIEW_LEVELS:
        if review_counts[level]:
            log.info("  %-12s: %s", level, review_counts[level])
    unresolved = {
        v: verdict_counts[v] for v in UNRESOLVED_VERDICTS if verdict_counts[v]
    }
    if unresolved:
        log.info(
            "  Not compared : %s",
            ", ".join(f"{k}={v}" for k, v in unresolved.items()),
        )
    if summary.dead_checks:
        log.error(
            "%s never succeeded on any of the %s rows, so this output says nothing "
            "about that question. Check the column and that the upstream API is "
            "reachable.",
            " and ".join(summary.dead_checks),
            total,
        )
    elif summary.degraded:
        log.error(
            "Nothing was compared on any of the %s rows. Check that --cas-column "
            "is the CAS column and that PubChem is reachable — this output says "
            "nothing about the sheet.",
            total,
        )
    log.info("Output: %s", output_path)
    return summary
