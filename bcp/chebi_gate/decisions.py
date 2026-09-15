"""Waivers and quarantine as versioned input data, keyed on CAS.

This is the module the re-run workflow turns on. A record is held back, a question
goes to the lab or to a chemist, an answer comes back weeks later, and the record
has to clear on a re-run without anyone editing code. So the decisions live in two
CSVs that a human edits and a reviewer reads:

``waivers.csv``
    ``cas, check_id, severity, reason, evidence, decided_by, decided_on`` -- a
    written decision that a finding is correct as it stands. It downgrades the
    finding to info; it never hides it.

``quarantine.csv``
    ``cas, check_id, reason, opened_on, blocked_on, resolved_on`` -- an open
    question. While ``resolved_on`` is empty the record is held, whatever the
    checks say, because "which salt does the lab hold?" is not answerable by
    looking at the file.

In the prototype both of these were Python dicts inside ``apply_text_fixes.py``,
and the waivers file was *generated* from them. That is backwards: a record left
quarantine by someone editing a script, with no date and no name against the
decision.

Two design points that cost real debugging time to learn:

**The key is the CAS number, not the NAME.** Renaming four records in the
prototype broke every name-keyed lookup and needed a rename map bolted on. One of
those four renames -- ``(S)-(+)-a-Methylhistamine dihydrobromide`` to
``...alpha-methylhistamine...`` -- was itself a quarantined record, so the rename
silently detached a record from its own open question.

**A waiver names the severity it was written against.** Without it the key is
``(cas, check_id)``, which waives every future finding of that check on that
record -- including one raised for an entirely different reason. CON-04 emits both
``medium`` for an unspecified centre and ``high`` for a name and a structure that
contradict each other, and a waiver for the first silently absorbed the second.
Severity is a stable pin, unlike the finding's prose, and a change of severity is
the signal that the check is now saying something else.

**A decision that matches nothing is reported, not ignored.** CAS is a better key
than NAME but not a permanent one: for one of the seven held records the defect
*is* the CAS number, so correcting it moves the key. A waiver or quarantine row
matching no record in the input is therefore a loud finding, which is what stops a
stale decision from silently ceasing to apply. So is a waiver whose record *is*
present but which waived nothing: 8 of the 15 shipped rows name a check the
record in question cannot emit under the default configuration, and a waiver that
applies to nothing reads exactly like one that applies.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from .checks import CHECKS

log = logging.getLogger(__name__)

DECISIONS_DIR = Path(__file__).parent / "decisions"
WAIVERS_FILE = "waivers.csv"
QUARANTINE_FILE = "quarantine.csv"

WAIVER_COLUMNS = (
    "cas",
    "check_id",
    "severity",
    "reason",
    "evidence",
    "decided_by",
    "decided_on",
)
QUARANTINE_COLUMNS = (
    "cas",
    "check_id",
    "reason",
    "opened_on",
    "blocked_on",
    "resolved_on",
)

# Minimum length for a reason to count as written down rather than gestured at.
MIN_REASON_LENGTH = 20


class DecisionsError(Exception):
    """Raised when the decisions data is malformed. Always fails the run."""


@dataclass(frozen=True)
class Waiver:
    """A recorded decision that a finding is correct as the record stands."""

    cas: str
    check_id: str
    severity: str
    reason: str
    evidence: str
    decided_by: str
    decided_on: date
    source_row: int = 0


@dataclass(frozen=True)
class Quarantine:
    """An open question against a record, and what would close it."""

    cas: str
    check_ids: tuple[str, ...]
    reason: str
    opened_on: date
    blocked_on: str
    resolved_on: date | None = None
    source_row: int = 0

    @property
    def is_open(self) -> bool:
        return self.resolved_on is None


@dataclass(frozen=True)
class Decisions:
    """Every decision recorded for a batch, indexed for lookup by CAS."""

    waivers: dict[tuple[str, str, str], Waiver]
    quarantine: dict[str, tuple[Quarantine, ...]]
    paths: tuple[Path, ...] = ()

    def waiver_for(self, cas: str, check_id: str, severity: str) -> Waiver | None:
        """The waiver for exactly this finding, or None.

        Keyed on the severity too, so a decision recorded about a medium finding
        does not absorb a high one that the same check raises later for a
        different reason.
        """
        return self.waivers.get((cas, check_id, severity))

    def open_questions(self, cas: str) -> tuple[Quarantine, ...]:
        """Unresolved quarantine rows for a CAS number."""
        return tuple(q for q in self.quarantine.get(cas, ()) if q.is_open)

    def resolved_questions(self, cas: str) -> tuple[Quarantine, ...]:
        return tuple(q for q in self.quarantine.get(cas, ()) if not q.is_open)

    def unmatched(
        self,
        present: set[str],
        applied: set[tuple[str, str, str]] | None = None,
    ) -> list[str]:
        """Decisions that did nothing, and why.

        Two ways a decision can be inert, and both read exactly like one that
        works. Its CAS may appear in no record -- a stale decision is worse than a
        missing one, because correcting a malformed CAS moves the key. Or its
        record may be right there and the finding it waives never raised, which
        ``_check_ids`` cannot catch: that check is registered, it simply cannot
        fire on this record, or fires at another severity. 8 of the 15 shipped
        waivers are in the second state under the default configuration.

        ``applied`` is the set of waiver keys that actually downgraded a finding.
        Pass None to report only the first kind, which is what a caller that did
        not apply the waivers can honestly say.
        """
        out = []
        for key in sorted(self.waivers):
            cas, check_id, severity = key
            if cas not in present:
                out.append(f"waiver {cas} {check_id} matches no record")
            elif applied is not None and key not in applied:
                out.append(
                    f"waiver {cas} {check_id}/{severity} matched a record but "
                    "waived no finding: the check did not fire on it, or fired "
                    "at another severity"
                )
        for cas in sorted(self.quarantine):
            if cas not in present:
                rows = self.quarantine[cas]
                state = "open" if any(r.is_open for r in rows) else "resolved"
                out.append(f"quarantine {cas} ({state}) matches no record")
        return out

    @property
    def counts(self) -> dict[str, int]:
        open_rows = sum(
            1 for rows in self.quarantine.values() for r in rows if r.is_open
        )
        total_rows = sum(len(rows) for rows in self.quarantine.values())
        return {
            "waivers": len(self.waivers),
            "quarantine_rows": total_rows,
            "quarantine_open": open_rows,
            "quarantine_resolved": total_rows - open_rows,
        }


def load(directory: str | Path | None = None) -> Decisions:
    """Load and validate both decision files. Raises on anything malformed."""
    base = Path(directory) if directory is not None else DECISIONS_DIR
    waivers_path = base / WAIVERS_FILE
    quarantine_path = base / QUARANTINE_FILE
    return Decisions(
        waivers=_load_waivers(waivers_path),
        quarantine=_load_quarantine(quarantine_path),
        paths=(waivers_path, quarantine_path),
    )


def _rows(path: Path, columns: tuple[str, ...]) -> list[tuple[int, dict[str, str]]]:
    """Read a decisions CSV, asserting its exact column set."""
    if not path.exists():
        raise DecisionsError(f"decisions file not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = tuple(reader.fieldnames or ())
        if header != columns:
            raise DecisionsError(
                f"{path.name} columns are {header}, expected {columns}"
            )
        # Row numbers count the header, so they match what an editor shows.
        return [(i, row) for i, row in enumerate(reader, start=2)]


def _require(path: Path, row_no: int, row: dict[str, str], field: str) -> str:
    value = (row.get(field) or "").strip()
    if not value:
        raise DecisionsError(f"{path.name} row {row_no}: {field} is empty")
    return value


def _parse_date(path: Path, row_no: int, field: str, value: str) -> date:
    try:
        return date.fromisoformat(value)
    except ValueError as exc:
        raise DecisionsError(
            f"{path.name} row {row_no}: {field} {value!r} is not an ISO date "
            "(YYYY-MM-DD)"
        ) from exc


def _check_ids(path: Path, row_no: int, raw: str) -> tuple[str, ...]:
    """Parse one or more check IDs, asserting each is registered.

    A decision naming a check that does not exist is a silent no-op: it looks like
    a considered judgement and applies to nothing. Usually it means a check was
    renamed and the decisions data was not.
    """
    ids = tuple(part for part in raw.replace("/", " ").split() if part)
    if not ids:
        raise DecisionsError(f"{path.name} row {row_no}: check_id is empty")
    unknown = [i for i in ids if i not in CHECKS]
    if unknown:
        raise DecisionsError(
            f"{path.name} row {row_no}: unknown check id(s) {unknown}; "
            f"registered checks are {sorted(CHECKS)}"
        )
    return ids


def _severity(path: Path, row_no: int, check_id: str, raw: str) -> str:
    """The severity a waiver was written against, asserted to be one the check emits.

    This is what ``_check_ids`` could not do. A registered check id says the
    decision names something real; it does not say the decision can ever apply.
    ``REL-01`` only ever emits ``info``, so a waiver for ``REL-01/high`` is a
    considered judgement about a finding that cannot exist, and it used to load
    cleanly.
    """
    allowed = CHECKS[check_id].severities
    if raw not in allowed:
        raise DecisionsError(
            f"{path.name} row {row_no}: {check_id} never emits {raw!r}; "
            f"it emits {list(allowed)}"
        )
    return raw


def _cas(path: Path, row_no: int, raw: str) -> str:
    """Validate the join key, which is useless if it is not a real CAS number."""
    from cas_registry import CAS_VALID, classify_cas

    normalised, verdict, repair = classify_cas(raw)
    if verdict != CAS_VALID or repair:
        raise DecisionsError(
            f"{path.name} row {row_no}: cas {raw!r} is not a well-formed CAS "
            f"Registry Number ({verdict}{'/' + repair if repair else ''}); "
            "it could never match a record"
        )
    return normalised


def _load_waivers(path: Path) -> dict[tuple[str, str], Waiver]:
    out: dict[tuple[str, str], Waiver] = {}
    for row_no, row in _rows(path, WAIVER_COLUMNS):
        cas = _cas(path, row_no, _require(path, row_no, row, "cas"))
        ids = _check_ids(path, row_no, _require(path, row_no, row, "check_id"))
        if len(ids) != 1:
            raise DecisionsError(
                f"{path.name} row {row_no}: a waiver covers exactly one check, "
                f"got {list(ids)}; write one row per check"
            )
        severity = _severity(
            path, row_no, ids[0], _require(path, row_no, row, "severity")
        )
        reason = _require(path, row_no, row, "reason")
        if len(reason) < MIN_REASON_LENGTH:
            raise DecisionsError(
                f"{path.name} row {row_no}: reason {reason!r} is too short to be a "
                "recorded decision"
            )
        # The handoff is explicit that a waiver with no evidence must fail the
        # build. A waiver is an assertion that the finding is wrong about this
        # record, and an assertion with nothing behind it is just a suppression.
        evidence = _require(path, row_no, row, "evidence")
        decided_by = _require(path, row_no, row, "decided_by")
        decided_on = _parse_date(
            path, row_no, "decided_on", _require(path, row_no, row, "decided_on")
        )

        key = (cas, ids[0], severity)
        if key in out:
            raise DecisionsError(
                f"{path.name} row {row_no}: duplicate waiver for {cas} {ids[0]} "
                f"{severity} (first seen at row {out[key].source_row})"
            )
        out[key] = Waiver(
            cas=cas,
            check_id=ids[0],
            severity=severity,
            reason=reason,
            evidence=evidence,
            decided_by=decided_by,
            decided_on=decided_on,
            source_row=row_no,
        )
    return out


def _load_quarantine(path: Path) -> dict[str, tuple[Quarantine, ...]]:
    grouped: dict[str, list[Quarantine]] = defaultdict(list)
    for row_no, row in _rows(path, QUARANTINE_COLUMNS):
        cas = _cas(path, row_no, _require(path, row_no, row, "cas"))
        ids = _check_ids(path, row_no, _require(path, row_no, row, "check_id"))
        reason = _require(path, row_no, row, "reason")
        if len(reason) < MIN_REASON_LENGTH:
            raise DecisionsError(
                f"{path.name} row {row_no}: reason {reason!r} is too short"
            )
        opened_on = _parse_date(
            path, row_no, "opened_on", _require(path, row_no, row, "opened_on")
        )
        # What would close this question. Without it a held record has no route
        # out, which is the state the prototype's prose left seven records in.
        blocked_on = _require(path, row_no, row, "blocked_on")

        raw_resolved = (row.get("resolved_on") or "").strip()
        resolved_on = (
            _parse_date(path, row_no, "resolved_on", raw_resolved)
            if raw_resolved
            else None
        )
        if resolved_on is not None and resolved_on < opened_on:
            raise DecisionsError(
                f"{path.name} row {row_no}: resolved_on {resolved_on} precedes "
                f"opened_on {opened_on}"
            )

        grouped[cas].append(
            Quarantine(
                cas=cas,
                check_ids=ids,
                reason=reason,
                opened_on=opened_on,
                blocked_on=blocked_on,
                resolved_on=resolved_on,
                source_row=row_no,
            )
        )
    return {cas: tuple(rows) for cas, rows in grouped.items()}
