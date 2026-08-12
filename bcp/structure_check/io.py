"""CSV batch I/O and single-row output for structure cross-checking."""

from __future__ import annotations

import csv
import json
import logging
import sys
from collections import Counter
from collections.abc import Collection
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
# check as a broken run, so the degraded verdict stays quiet. Lifted when an
# outage is in play — see RunSummary.nothing_compared.
MIN_ROWS_FOR_DEGRADED = 5

# Below this many unanswered rows, an outage is not yet evidence about the run.
# Separate from the row floor above because it counts a different thing: not how
# big the sheet was, but how much of it an unreachable upstream actually ate.
# Without it the majority test in RunSummary.outages goes degenerate — against
# zero resolved rows, one outage is a "majority" — and a sparse or wrong column
# would be condemned by a single throttled request.
MIN_OUTAGE_ROWS_FOR_DEGRADED = 5

# If a column produced nothing usable across at least this many *distinct* values,
# the likelier explanation is that it is not the column it was said to be than that
# every value in it is independently bad. Counted per distinct value, not per row:
# 500 rows of one junk string is a single mistake, not 500 pieces of evidence.
# Mirrors chebi_terms.SUSPICIOUS_TOTAL_MISS, which draws the same line.
SUSPICIOUS_DISTINCT_MISSES = 5

# A side that has failed this many attempts in a row is treated as down, and stops
# being queried: the verdict is already settled and every further request costs
# ~14s of backoff to confirm it. Mirrors chebi_terms.MAX_CONSECUTIVE_FAILURES.
MAX_CONSECUTIVE_OUTAGE_ROWS = 5

# Coupled to MIN_OUTAGE_ROWS_FOR_DEGRADED above: a side that trips has by
# construction already met that floor, so past a trip only the majority test is
# still doing work. Lowering the trip threshold below the floor would let a side
# be dropped before its outage counts as evidence.
#
# ...but "down" is an inference, so it is re-tested this often. Long enough that a
# genuine outage costs few probes, short enough that a sparse column which tripped
# on scattered blips recovers well inside a normal sheet.
OUTAGE_PROBE_INTERVAL = 20


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
    missing_rows: int = 0
    outage_rows: int = 0

    @property
    def answered_rows(self) -> int:
        """
        Rows where the upstream said something, whether or not it was useful.

        "No such compound" is an answer. Counting only the rows that *resolved*
        would let five throttled rows outweigh a hundred honest 404s, which reads
        as an outage on precisely the sheet — a column of internal codes PubChem
        does not carry — where the answer is that the column is wrong.
        """
        return self.resolved_rows + self.missing_rows


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
    def sides(self) -> tuple[tuple[str, str, SideTally], ...]:
        """(human label, the flag that selects it, its tally), in report order."""
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

        The majority is against every row the upstream *answered*, not just the
        ones that resolved: a definite "no such compound" is an answer, and
        leaving it out let five throttled rows outweigh a hundred honest 404s.
        That misfires exactly where it costs most — against a wrong column
        nothing resolves by construction, so the outage branch would fire and
        suppress the diagnosis that would have named the real problem.

        A majority needs no tuned constant: once most of what was asked went
        unanswered, the output is not a verdict on the sheet whatever the rest says.
        Rows below that line are still marked individually, in the `unasked` column.

        The majority alone is not enough, though, because it goes degenerate when
        nothing resolved: against zero answers a single outage is a majority. That
        is the normal shape of a *sparse* column (two ChEBI IDs on a 117-row sheet,
        both unlucky) and of a *wrong* one (nothing resolves by definition), so
        without a floor one transient failure would condemn the run in exactly the
        two cases this rewrite exists to stop condemning. The floor is on the
        outage count itself, not on the sheet: what has to be substantial is the
        evidence of an outage, not the size of the input.
        """
        return tuple(
            label
            for label, _flag, tally in self.sides
            if tally.outage_rows > tally.answered_rows
            and tally.outage_rows >= MIN_OUTAGE_ROWS_FOR_DEGRADED
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

        Disqualified by that side appearing in `outages`, not merely by its outage
        count being non-zero. A wrong column resolves nothing by definition, so any
        single throttled request on it would otherwise suppress this diagnosis and
        report an outage instead — telling the operator to re-run a multi-minute
        sheet that will fail again for the reason nobody named. One unanswered row
        beside a hundred answered ones is not a competing explanation for the
        hundred; a genuine outage is, and that is what `outages` already decides.
        """
        out = self.outages
        return tuple(
            label
            for label, _flag, tally in self.sides
            if not tally.resolved_values
            and label not in out
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

        The row floor lifts when an outage is in play, which is what keeps a small
        batch from being quieter than the same input in single mode: four rows
        against a dead ChEBI sit under every other threshold, while `--cas X
        --chebi Y` on one of those same rows exits 1. An outage is positive
        evidence that we failed to ask, so it needs no volume to be believed.
        Absence of results *without* one is only evidence in bulk — two rows of
        compounds PubChem happens not to know is a plausible spot check, not a
        broken run — so that case keeps the floor.

        This does not reopen the sparse-column hole the outage floor closed. Two
        unlucky ChEBI IDs on a 117-row sheet leave the other 115 rows comparing
        fine, so `nothing_compared` is false on the count alone.
        """
        if not self.rows or self.review_counts[REVIEW_UNVERIFIED] != self.rows:
            return False
        if any(tally.outage_rows for _label, _flag, tally in self.sides):
            return True
        return self.rows >= MIN_ROWS_FOR_DEGRADED

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


class _OutageBreaker:
    """
    Stops querying an upstream that has stopped answering, per side.

    Tracked per side rather than per row because a whole-row test cannot see the
    case that costs the most: with ChEBI dead run-wide and a live name column
    whose names match, every row is `ok`, so the streak resets on every row and
    the sheet walks all of itself paying full retry backoff to reach a verdict
    `outages` settled by row 5.

    Once a side has failed MAX_CONSECUTIVE_OUTAGE_ROWS attempts in a row it is
    tripped: its rows are marked unreachable without a request, which is what
    they are, and the run finishes the columns that still work rather than
    aborting. The alternative — killing the whole run — throws away the name
    check on 400 remaining rows because ChEBI is down.

    A tripped side is a *guess* that the upstream is down, so it is re-checked
    every OUTAGE_PROBE_INTERVAL attempts. Without that, five scattered blips on a
    sparsely mapped column would poison every remaining row with an `unreachable`
    no request was ever made for — the exact conflation of an answer with a
    silence this module exists to prevent. A successful probe clears the trip.

    A blank cell neither increments nor resets a streak: it says nothing about
    whether the upstream is up, and a sparsely mapped column is normal input.
    """

    SIDES = ("cas", "chebi", "name")

    def __init__(self) -> None:
        self._streaks: Counter[str] = Counter()
        self._tripped: set[str] = set()
        self._since_probe: Counter[str] = Counter()

    def skip(self) -> frozenset[str]:
        """Sides to leave unqueried on the next row."""
        return frozenset(
            side
            for side in self._tripped
            if self._since_probe[side] < OUTAGE_PROBE_INTERVAL
        )

    def note(
        self,
        result: dict[str, Any],
        skipped: frozenset[str],
        *,
        present: Collection[str] = (),
    ) -> tuple[list[str], list[str]]:
        """
        Fold one row in. Returns (sides newly tripped, sides newly recovered).

        `present` names the sides whose cell was non-blank, so the probe interval
        counts attempts rather than rows: on a sparsely mapped column, rows with
        nothing in them are not tries and must not bring the next probe forward.
        """
        tripped, recovered = [], []
        for side in self.SIDES:
            if side not in present:
                continue
            if side in skipped:
                self._since_probe[side] += 1
                continue
            status = result.get(f"{side}_status")
            if status == IDENTIFIER_UNREACHABLE:
                self._streaks[side] += 1
                if side in self._tripped:
                    self._since_probe[side] = 0  # the probe failed; wait again
                elif self._streaks[side] >= MAX_CONSECUTIVE_OUTAGE_ROWS:
                    self._tripped.add(side)
                    self._since_probe[side] = 0
                    tripped.append(side)
            elif status in (IDENTIFIER_RESOLVED, IDENTIFIER_MISSING):
                self._streaks[side] = 0
                if side in self._tripped:
                    self._tripped.discard(side)
                    recovered.append(side)
        return tripped, recovered


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
        self._missing_rows: Counter[str] = Counter()
        self._outage_rows: Counter[str] = Counter()

    # RESOLVED beats MISSING beats UNREACHABLE: an answer outranks a silence.
    _RANK = {IDENTIFIER_RESOLVED: 3, IDENTIFIER_MISSING: 2, IDENTIFIER_UNREACHABLE: 1}

    def add(
        self,
        result: dict[str, Any],
        *,
        cas: str,
        chebi_id: str,
        name: str,
        skipped: Collection[str] = (),
    ) -> None:
        for side, value in (("cas", cas), ("chebi", chebi_id), ("name", name)):
            status = result.get(f"{side}_status")
            if not value or status not in self._RANK:
                # Blank cells and sides that were never checked are not evidence
                # about anything: a partly-filled sheet is normal input.
                continue
            if side in skipped:
                # The breaker inferred this row's `unreachable` from the rows that
                # tripped it; it is not independent evidence, and folding it in
                # would let the inference confirm itself. Five real failures
                # would become twenty-five outage rows and flip a run that
                # recovered into exit 1 — manufacturing the very verdict the
                # majority test exists to require evidence for. The row is still
                # honestly marked unreachable in `unasked`; it just does not vote.
                continue
            if status == IDENTIFIER_RESOLVED:
                self._resolved_rows[side] += 1
            elif status == IDENTIFIER_MISSING:
                self._missing_rows[side] += 1
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
            missing_rows=self._missing_rows[side],
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

    if chebi_column and chebi_column == name_column:
        # Not a false all-clear — the ChEBI side rejects a compound name without a
        # request — but every row comes back name_unresolved, which reads like a
        # problem with the data rather than with the invocation.
        raise StructureCheckError(
            f"Column '{chebi_column}' was given as both --chebi-column and "
            f"--name-column. One of them is pointed at the wrong column."
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
    breaker = _OutageBreaker()

    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=output_fields)
        writer.writeheader()

        for i, row in enumerate(rows, 1):
            cas = (row.get(cas_column) or "").strip()
            chebi_id = (row.get(chebi_column) or "").strip() if chebi_column else ""
            name = (row.get(name_column) or "").strip() if name_column else ""

            key = (cas, chebi_id, name)
            from_cache = key in cache
            skipped = frozenset() if from_cache else breaker.skip()
            if not from_cache:
                # Only counted when something was actually asked, which needs the
                # pivot *and* something to compare it against: check_row returns
                # without a request when either is missing, and again when the
                # breaker has skipped every side it would have queried. Counting
                # those would overstate "N distinct lookups" in the same
                # direction the cache undercounted it — and that number is what
                # an operator uses to sanity-check runtime.
                requested = {s for s, v in (("chebi", chebi_id), ("name", name)) if v}
                if cas and requested and "cas" not in skipped:
                    attempted_keys.add(key)
                result = check_row(cas=cas, chebi_id=chebi_id, name=name, skip=skipped)
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
            tallies.add(result, cas=cas, chebi_id=chebi_id, name=name, skipped=skipped)

            # A cache hit is not evidence the upstream recovered — since
            # outages are never cached, a hit means a row that succeeded
            # *earlier*, possibly long before the upstream died. Folding one in
            # would let a sheet with repeated compounds clear the streak
            # forever. chebi_terms refuses this for the same reason.
            if not from_cache:
                present = {
                    s
                    for s, v in (("cas", cas), ("chebi", chebi_id), ("name", name))
                    if v
                }
                tripped, recovered = breaker.note(result, skipped, present=present)
                for side in tripped:
                    log.error(
                        "%s has been unreachable for %s attempts in a row, at row "
                        "%s of %s. Not asking it again for a while: those rows are "
                        "marked unreachable in the 'unasked' column, and the run "
                        "continues on the columns that still answer.",
                        side,
                        MAX_CONSECUTIVE_OUTAGE_ROWS,
                        i,
                        total,
                    )
                for side in recovered:
                    log.warning("%s is answering again at row %s.", side, i)

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
    for label, _flag, tally in summary.sides:
        if tally.outage_rows and label not in summary.outages:
            log.warning(
                "%s went unanswered on %s of the %s rows it was tried on, because "
                "the upstream could not be reached. Those rows are marked in the "
                "'unasked' column; a re-run would answer them.",
                label,
                tally.outage_rows,
                tally.outage_rows + tally.answered_rows,
            )

    for label, _flag, tally in summary.sides:
        if label in summary.outages:
            log.error(
                "%s went unanswered on %s of the %s rows it was tried on — more "
                "than the upstream answered at all (%s), so this output is not a "
                "verdict on that question. The API was unreachable, not merely "
                "unhelpful; re-run once it is back.",
                label,
                tally.outage_rows,
                tally.outage_rows + tally.answered_rows,
                tally.answered_rows,
            )
    for label, flag, tally in summary.sides:
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
