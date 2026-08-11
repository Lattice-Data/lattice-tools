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
    IDENTIFIER_MISSING,
    IDENTIFIER_RESOLVED,
    IDENTIFIER_UNREACHABLE,
    OUTPUT_FIELDS_APPENDED,
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

# Below this many rows, "nothing was compared" is as likely to be a two-row spot
# check as a broken run, so the degraded verdict stays quiet.
MIN_ROWS_FOR_DEGRADED = 5

# If a column produced nothing usable across at least this many *distinct* values,
# the likelier explanation is that it is not the column it was said to be than that
# every value in it is independently bad. Counted per distinct value, not per row:
# 500 rows of one junk string is a single mistake, not 500 pieces of evidence.
# Mirrors chebi_terms.SUSPICIOUS_TOTAL_MISS, which draws the same line.
SUSPICIOUS_DISTINCT_MISSES = 5

# A run against a wholly unreachable upstream is decided by its first few rows but
# would otherwise grind through every remaining one at ~14s of backoff per request.
# Mirrors chebi_terms.MAX_CONSECUTIVE_FAILURES.
MAX_CONSECUTIVE_OUTAGE_ROWS = 5


class SideTally(NamedTuple):
    """
    What happened to one identifier column over a run.

    Values and rows are counted separately because the two questions this feeds
    need different units. "Is this column what it claims to be?" is evidence about
    *values*, and repeating one bad value does not make it stronger. "Did enough of
    this run actually happen?" is about *rows*, because a row is what goes
    unanswered.
    """

    resolved_values: int = 0
    missing_values: int = 0
    resolved_rows: int = 0
    outage_rows: int = 0


SINGLE_INPUT_FIELDS = ["name", "cas", "chebi_id"]
SINGLE_OUTPUT_FIELDS = [*SINGLE_INPUT_FIELDS, *OUTPUT_FIELDS_APPENDED]


class StructureCheckError(Exception):
    """Raised when input file or column validation fails."""


class RunSummary(NamedTuple):
    """
    Outcome of a batch run.

    `degraded` answers "can this output be read as a verdict on the sheet?", and it
    is deliberately built from three *separate* signals rather than one predicate.
    A single boolean here was wrong four times running, each fix correct for the
    input it targeted and blind to the next, because it was being asked two
    unrelated questions at once: "was upstream up?" and "is this the right column?".
    Those have different evidence, different units, and different remedies — a
    re-run fixes one and cannot fix the other — so they are answered apart.
    """

    review_counts: Counter[str]
    verdict_counts: Counter[str]
    rows: int
    cas: SideTally = SideTally()
    chebi: SideTally = SideTally()
    name: SideTally = SideTally()

    @property
    def needs_attention(self) -> int:
        return self.review_counts[REVIEW_INVESTIGATE] + self.review_counts[REVIEW_CHECK]

    @property
    def _sides(self) -> tuple[tuple[str, str, SideTally], ...]:
        """(human label, the flag that selects it, its tally)."""
        return (
            ("CAS", "--cas-column", self.cas),
            ("ChEBI ID", "--chebi-column", self.chebi),
            ("name", "--name-column", self.name),
        )

    @property
    def outages(self) -> tuple[str, ...]:
        """
        Sides where more rows went unanswered by an outage than were answered.

        Not "any outage at all", which is what the sibling chebi_terms can afford:
        it makes one request per ID, while this tool makes 3 to 44 per row against
        an API that rate-limits, and chebi_lookup cannot tell a throttle that
        outlasted its retries from a genuine outage. Firing on a single occurrence
        would exit 1 on healthy large runs, and an exit code that cries wolf stops
        being read — which would cost more than the case it was added to catch.

        A majority needs no tuned constant: once most of what was asked went
        unanswered, the output is not a verdict on the sheet whatever the rest says.
        Rows below that line are still marked individually, in the `unasked` column.
        """
        return tuple(
            label
            for label, _flag, tally in self._sides
            if tally.outage_rows > tally.resolved_rows and tally.outage_rows
        )

    @property
    def suspect_columns(self) -> tuple[str, ...]:
        """
        Columns where nothing resolved across enough distinct values to look wrong.

        The cheap mistake, and so the likely one: point --chebi-column at a notes
        column and every row fails for a reason that has nothing to do with ChEBI.

        Disqualified per side by that side's *own* upstream being out, never
        globally. A ChEBI column full of free text is rejected by a regex without a
        single request, so a PubChem outage is logically incapable of explaining it;
        letting one suppress the other would hide a real misconfiguration behind an
        unrelated blip. An all-blank column is not suspect either — a partly-filled
        sheet is normal input — because blanks never enter these counts.
        """
        return tuple(
            label
            for label, _flag, tally in self._sides
            if not tally.resolved_values
            and not tally.outage_rows
            and tally.missing_values >= SUSPICIOUS_DISTINCT_MISSES
        )

    @property
    def nothing_compared(self) -> bool:
        """
        Every row came out unverified: this run produced no comparison at all.

        Kept as its own term because it is the only one that reasons about
        *comparisons performed* rather than about identifiers resolved, and the two
        genuinely diverge. A ChEBI column of class terms resolves perfectly — ChEBI
        holds every record — and still compares nothing, because a class term
        carries no structure. Both other signals stay silent on that sheet.
        """
        return self.rows >= MIN_ROWS_FOR_DEGRADED and (
            self.review_counts[REVIEW_UNVERIFIED] == self.rows
        )

    @property
    def degraded(self) -> bool:
        """True when the run cannot be read as a verdict on the sheet."""
        return bool(self.outages or self.suspect_columns or self.nothing_compared)


def _has_outage(result: dict[str, Any]) -> bool:
    """True when any identifier in this row could not be reached."""
    return any(
        result.get(field) == IDENTIFIER_UNREACHABLE
        for field in ("cas_status", "chebi_status", "name_status")
    )


def _is_total_outage_row(result: dict[str, Any]) -> bool:
    """
    An outage row that yielded no comparison at all.

    Not "no identifier resolved": with ChEBI down and only a ChEBI column asked
    for, every CAS still resolves perfectly and every row still compares nothing.
    What makes a run worth abandoning is that it is spending requests and
    producing no verdicts, which is exactly `unverified` plus an unreachable
    upstream.
    """
    return _has_outage(result) and result.get("review") == REVIEW_UNVERIFIED


class _SideAccumulator:
    """
    Per-column tallies, kept in the units each question needs.

    Reads the per-identifier statuses check_row reports rather than re-deriving
    them from the sheet cell. That distinction is the whole point: the cell says
    what was *asked for*, only the resolver knows what *happened*, and every
    previous version of the run verdict was wrong because it divided one by the
    other. A status cannot disagree with itself.
    """

    SIDES = ("cas", "chebi", "name")

    def __init__(self) -> None:
        # value -> best status seen for it, so one transient failure does not
        # condemn a value that resolved on another row.
        self._values: dict[str, dict[str, str]] = {s: {} for s in self.SIDES}
        self._resolved_rows: Counter[str] = Counter()
        self._outage_rows: Counter[str] = Counter()

    # RESOLVED beats MISSING beats UNREACHABLE: an answer outranks a silence.
    _RANK = {IDENTIFIER_RESOLVED: 3, IDENTIFIER_MISSING: 2, IDENTIFIER_UNREACHABLE: 1}

    def add(
        self, result: dict[str, Any], *, cas: str, chebi_id: str, name: str
    ) -> None:
        for side, value in (("cas", cas), ("chebi", chebi_id), ("name", name)):
            status = result.get(f"{side}_status")
            if not value or status not in self._RANK:
                # Blank cells and sides that were never checked are not evidence
                # about anything: a partly-filled sheet is normal input.
                continue
            if status == IDENTIFIER_RESOLVED:
                self._resolved_rows[side] += 1
            elif status == IDENTIFIER_UNREACHABLE:
                self._outage_rows[side] += 1
            best = self._values[side].get(value)
            if best is None or self._RANK[status] > self._RANK[best]:
                self._values[side][value] = status

    def _tally(self, side: str) -> SideTally:
        statuses = self._values[side].values()
        return SideTally(
            resolved_values=sum(1 for s in statuses if s == IDENTIFIER_RESOLVED),
            missing_values=sum(1 for s in statuses if s == IDENTIFIER_MISSING),
            resolved_rows=self._resolved_rows[side],
            outage_rows=self._outage_rows[side],
        )

    def summary(
        self, *, review_counts: Counter[str], verdict_counts: Counter[str], rows: int
    ) -> RunSummary:
        return RunSummary(
            review_counts=review_counts,
            verdict_counts=verdict_counts,
            rows=rows,
            cas=self._tally("cas"),
            chebi=self._tally("chebi"),
            name=self._tally("name"),
        )


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

    requested = [
        ("--cas-column", cas_column),
        ("--chebi-column", chebi_column),
        ("--name-column", name_column),
    ]
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
    # Tracked apart from the cache, which no longer holds outages: "0 distinct
    # lookups" on a run that made hundreds of requests would be a lie.
    attempted_keys: set[tuple[str, str, str]] = set()
    review_counts: Counter[str] = Counter()
    verdict_counts: Counter[str] = Counter()
    tallies = _SideAccumulator()
    consecutive_outage_rows = 0

    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=output_fields)
        writer.writeheader()

        for i, row in enumerate(rows, 1):
            cas = (row.get(cas_column) or "").strip()
            chebi_id = (row.get(chebi_column) or "").strip() if chebi_column else ""
            name = (row.get(name_column) or "").strip() if name_column else ""

            key = (cas, chebi_id, name)
            if key not in cache:
                attempted_keys.add(key)
                result = check_row(cas=cas, chebi_id=chebi_id, name=name)
                # Outages are not cached. A transient failure on one triple must
                # not be replayed as a settled answer for every row that repeats
                # it — the same rule parent_inchikey already applies to its cache,
                # and the reason chebi_terms never caches a failure either.
                if not _has_outage(result):
                    cache[key] = result
            else:
                result = cache[key]

            out_row = {c: row.get(c, "") for c in carried}
            out_row.update({f: result.get(f, "") for f in OUTPUT_FIELDS_APPENDED})
            writer.writerow(out_row)

            review = result["review"]
            review_counts[review] += 1
            for field in ("id_cas_verdict", "name_cas_verdict"):
                verdict_counts[result[field]] += 1
            tallies.add(result, cas=cas, chebi_id=chebi_id, name=name)

            if _is_total_outage_row(result):
                consecutive_outage_rows += 1
                if consecutive_outage_rows >= MAX_CONSECUTIVE_OUTAGE_ROWS:
                    raise StructureCheckError(
                        f"{consecutive_outage_rows} rows in a row resolved nothing "
                        f"because an upstream could not be reached, at row {i} of "
                        f"{total}. Stopping rather than spending the rest of the "
                        f"sheet on backoff. The partial output is at {output_path}."
                    )
            else:
                consecutive_outage_rows = 0

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

    summary = tallies.summary(
        review_counts=review_counts, verdict_counts=verdict_counts, rows=total
    )
    _report(summary, total=total, distinct=len(attempted_keys), output_path=output_path)
    return summary


def _report(
    summary: RunSummary, *, total: int, distinct: int, output_path: Path
) -> None:
    """Log the tally, then say plainly why the run cannot be trusted, if it cannot."""
    log.info("=" * 58)
    log.info("Done. %s rows (%s distinct lookups).", total, distinct)
    for level in REVIEW_LEVELS:
        if summary.review_counts[level]:
            log.info("  %-12s: %s", level, summary.review_counts[level])
    unresolved = {
        v: summary.verdict_counts[v]
        for v in UNRESOLVED_VERDICTS
        if summary.verdict_counts[v]
    }
    if unresolved:
        log.info(
            "  Not compared : %s",
            ", ".join(f"{k}={v}" for k, v in unresolved.items()),
        )

    # An outage below the majority line does not degrade the run, but it still left
    # rows unanswered, and those rows are only distinguishable from genuine answers
    # if someone says so. Reported whether or not the run is degraded.
    for label, _flag, tally in summary._sides:
        if tally.outage_rows and label not in summary.outages:
            log.warning(
                "%s went unanswered on %s of %s rows because the upstream could not "
                "be reached. Those rows are marked in the 'unasked' column; a re-run "
                "would answer them.",
                label,
                tally.outage_rows,
                total,
            )

    for label, _flag, tally in summary._sides:
        if label in summary.outages:
            log.error(
                "%s could not be resolved on %s of %s rows — more than were answered "
                "(%s), so this output is not a verdict on that question. The upstream "
                "API was unreachable, not merely unhelpful; re-run once it is back.",
                label,
                tally.outage_rows,
                total,
                tally.resolved_rows,
            )
    for label, flag, tally in summary._sides:
        if label in summary.suspect_columns:
            log.error(
                "%s (%s) produced nothing usable across %s distinct values, with no "
                "sign of an outage. That is the signature of a column holding "
                "something other than %s values — check %s.",
                label,
                flag,
                tally.missing_values,
                label,
                flag,
            )
    if summary.nothing_compared:
        log.error(
            "Nothing was compared on any of the %s rows, so this output says nothing "
            "about the sheet. Every row resolved at most one side, which is what a "
            "sheet of structureless ChEBI entries or a missing CAS column looks like.",
            total,
        )
    log.info("Output: %s", output_path)
