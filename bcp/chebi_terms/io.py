"""CSV batch I/O and single-ID output for ChEBI ID verification."""

from __future__ import annotations

import csv
import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Literal, NamedTuple

from . import client
from .client import (
    CAS_NOT_IN_CHEBI,
    NAME_MISMATCH,
    NOT_CHECKED,
    OUTPUT_FIELDS_APPENDED,
    STATUS_INVALID,
    STATUS_LOOKUP_FAILED,
    STATUS_MISSING,
    STATUS_NOT_FOUND,
    STATUS_NOT_RELEASED,
    STATUS_OK,
    STATUS_SECONDARY,
    ChebiUnavailableError,
    empty_result,
    fetch_compound,
    normalize_chebi_id,
    verify_chebi_id,
    verify_payload,
)

log = logging.getLogger(__name__)

# Retrying every remaining row against a ChEBI that is plainly down costs 6s of
# backoff each — plus up to three connect timeouts if packets are being dropped
# rather than refused — and yields nothing but lookup_failed. Fail fast instead.
MAX_CONSECUTIVE_FAILURES = 5

# If every ID that reached ChEBI came back not_found, across at least this many
# *distinct IDs*, the likelier explanation is a moved endpoint than a sheet of
# compounds that all happen not to exist. EBI has already done this to the flat
# files and the SOAP service; a rename here would otherwise print a clean tally
# and exit 0. Counted per ID, not per row — see all_lookups_missed().
SUSPICIOUS_TOTAL_MISS = 5

# Statuses reported in the run summary, in the order they are logged.
SUMMARY_STATUSES = (
    STATUS_OK,
    STATUS_SECONDARY,
    STATUS_NOT_RELEASED,
    STATUS_NOT_FOUND,
    STATUS_MISSING,
    STATUS_INVALID,
    STATUS_LOOKUP_FAILED,
)


class ChebiTermsError(Exception):
    """Raised when input file or column validation fails."""


SINGLE_CHEBI_COLUMN = "chebi_id"
SINGLE_OUTPUT_FIELDS = [SINGLE_CHEBI_COLUMN, *OUTPUT_FIELDS_APPENDED]


class RunSummary(NamedTuple):
    """Outcome of a batch run: the per-status tally plus the run-level verdict."""

    status_counts: Counter[str]
    missed_ids: int
    resolved_ids: int
    invalid_values: int
    missing_rows: int
    name_mismatches: int
    cas_disagreements: int

    @property
    def suspect_endpoint(self) -> bool:
        """Every ID that reached ChEBI came back not_found. See all_lookups_missed."""
        return all_lookups_missed(self.missed_ids, self.resolved_ids)

    @property
    def suspect_column(self) -> bool:
        """
        Nothing was usable as a ChEBI ID at all, across enough rows to look like
        --chebi-column is pointed at the wrong column.

        The mirror of suspect_endpoint from the cheaper direction: it needs no
        network, so it is the easier mistake to make and the likelier one to slip
        past unnoticed. Distinct values for the invalid side — 500 rows of one
        junk string is a single mistake — and rows for the blank side, where
        there is no value to be distinct about.
        """
        if self.resolved_ids or self.missed_ids:
            return False
        return (self.invalid_values + self.missing_rows) >= SUSPICIOUS_TOTAL_MISS

    @property
    def degraded(self) -> bool:
        """True when the run as a whole cannot be trusted, whatever the rows say."""
        return bool(
            self.status_counts[STATUS_LOOKUP_FAILED]
            or self.suspect_endpoint
            or self.suspect_column
        )


def all_lookups_missed(missed_ids: int, resolved_ids: int) -> bool:
    """
    True when every *distinct* ID that reached ChEBI came back not_found, across
    enough of them to be implausible.

    Counted per ID rather than per row: five rows carrying one retired ID make a
    single request, and "five distinct compounds all happen not to exist" is the
    inference the threshold is built on — not "one ID repeated five times".

    Per-row `not_found` stays honest either way; this only decides the run-level
    verdict. It also fires on a sheet of genuinely bogus IDs — equally worth a
    human look, so the false positive is benign.
    """
    return not resolved_ids and missed_ids >= SUSPICIOUS_TOTAL_MISS


def check_output_path(output_path: Path | None) -> None:
    """
    Fail on an unwritable output before spending any requests.

    Single-ID mode writes after the fetch, so without this an unwritable path
    throws the answer away along with the request that earned it.
    """
    if output_path is None:
        return
    if output_path.is_dir():
        raise ChebiTermsError(f"Output path is a directory: {output_path}")
    parent = output_path.parent
    if not parent.is_dir():
        raise ChebiTermsError(f"Output directory does not exist: {parent}")


def build_single_row(chebi_id: str, result: dict[str, Any]) -> dict[str, Any]:
    """Build one output row for a single-ID lookup."""
    row: dict[str, Any] = {SINGLE_CHEBI_COLUMN: chebi_id}
    for field in OUTPUT_FIELDS_APPENDED:
        value = result.get(field, "")
        row[field] = "" if value is None else value
    return row


def _log_single_summary(chebi_id: str, result: dict[str, Any]) -> None:
    status = result.get("id_status", "")
    if status == STATUS_INVALID:
        log.warning("Not a ChEBI identifier: %s", chebi_id)
        return
    if status == STATUS_MISSING:
        log.warning("No ChEBI identifier given.")
        return
    if status == STATUS_LOOKUP_FAILED:
        log.error(
            "Could not reach ChEBI for %s — this says nothing about whether the "
            "compound exists.",
            chebi_id,
        )
        return
    if status != STATUS_OK:
        log.warning("ChEBI %s: %s", chebi_id, status)
    else:
        log.info("ChEBI %s: %s", result.get("chebi_accession", chebi_id), status)

    if result.get("chebi_name"):
        log.info("Name: %s", result["chebi_name"])
    for field, label in (("name_verdict", "Name check"), ("cas_verdict", "CAS check")):
        verdict = result.get(field, NOT_CHECKED)
        if verdict != NOT_CHECKED:
            log.info("%s: %s", label, verdict)


def emit_single_chebi_id(
    chebi_id: str,
    output_path: Path | None,
    *,
    fmt: Literal["json", "csv"] = "json",
    expected_name: str | None = None,
    expected_cas: str | None = None,
    max_synonyms: int | None = None,
) -> dict[str, Any]:
    """
    Resolve one ChEBI ID and write JSON or CSV to a file or stdout.

    Returns the result so callers can distinguish a verdict from a failed lookup.
    """
    chebi_id = chebi_id.strip()
    check_output_path(output_path)
    result = verify_chebi_id(
        chebi_id,
        expected_name=expected_name,
        expected_cas=expected_cas,
        max_synonyms=max_synonyms,
    )
    row = build_single_row(chebi_id, result)
    _log_single_summary(chebi_id, result)

    if fmt == "json":
        payload = json.dumps(row, indent=2, ensure_ascii=False)
        if output_path is None:
            print(payload)
        else:
            output_path.write_text(payload + "\n", encoding="utf-8")
            log.info("Output: %s", output_path)
        return result

    if output_path is None:
        out_fh = sys.stdout
    else:
        out_fh = open(output_path, "w", newline="", encoding="utf-8")

    try:
        writer = csv.DictWriter(out_fh, fieldnames=SINGLE_OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerow({k: str(v) if v is not None else "" for k, v in row.items()})
    finally:
        if output_path is not None:
            out_fh.close()
            log.info("Output: %s", output_path)

    return result


def _resolve_columns(
    fieldnames: list[str],
    chebi_column: str,
    name_column: str | None,
    cas_column: str | None,
) -> list[str]:
    """Validate requested columns exist and return the input columns to carry through."""
    if chebi_column not in fieldnames:
        raise ChebiTermsError(
            f"Column '{chebi_column}' not found. Available columns: {fieldnames}"
        )
    if chebi_column in OUTPUT_FIELDS_APPENDED:
        # Would emit two columns of the same name and corrupt the output.
        raise ChebiTermsError(
            f"Column '{chebi_column}' collides with a column this tool writes. "
            f"Rename it in the input. Reserved names: {OUTPUT_FIELDS_APPENDED}"
        )
    for label, column in (("--name-column", name_column), ("--cas-column", cas_column)):
        if column and column not in fieldnames:
            raise ChebiTermsError(
                f"Column '{column}' ({label}) not found. Available columns: {fieldnames}"
            )
        if column and column in OUTPUT_FIELDS_APPENDED:
            # The comparison itself would be correct, but the column gets
            # overwritten with ChEBI's value — so a mismatch row would show a
            # verdict with no record of what it was compared against.
            raise ChebiTermsError(
                f"Column '{column}' ({label}) collides with a column this tool "
                f"writes, which would erase the value being checked. Rename it in "
                f"the input. Reserved names: {OUTPUT_FIELDS_APPENDED}"
            )

    extra_columns = [c for c in fieldnames if c != chebi_column]
    clobbered = [c for c in extra_columns if c in OUTPUT_FIELDS_APPENDED]
    if clobbered:
        log.warning(
            "Input columns overwritten by ChEBI values: %s",
            ", ".join(clobbered),
        )
    return [c for c in extra_columns if c not in OUTPUT_FIELDS_APPENDED]


def verify_chebi_file(
    input_path: Path,
    chebi_column: str,
    output_path: Path,
    *,
    name_column: str | None = None,
    cas_column: str | None = None,
    max_synonyms: int | None = None,
) -> RunSummary:
    """
    Read a CSV of ChEBI IDs, verify each against ChEBI, and write enriched output.

    Returns a RunSummary so callers can tell a clean run from a degraded one.
    Raises ChebiTermsError if ChEBI stays unreachable for
    MAX_CONSECUTIVE_FAILURES consecutive lookups that reached the network,
    leaving the partial output in place.
    """
    # is_file(), not exists(): a directory would otherwise reach open() and
    # surface as a raw IsADirectoryError, and `--input ""` normalizes to Path('.').
    if not input_path.is_file():
        raise ChebiTermsError(f"Input file not found (or not a file): {input_path}")

    with open(input_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        extra_columns = _resolve_columns(
            fieldnames, chebi_column, name_column, cas_column
        )
        rows = list(reader)

    output_fields = [chebi_column] + extra_columns + OUTPUT_FIELDS_APPENDED
    total = len(rows)
    # One payload per unique ID: the same compound often appears on many rows,
    # and per-row expectations are evaluated against the cached payload. Only
    # definitive answers land here — see the lookup_failed branch below.
    payload_cache: dict[str, dict[str, Any] | None] = {}
    seen_ids: set[str] = set()
    invalid_values: set[str] = set()
    status_counts: Counter[str] = Counter()
    name_mismatches = 0
    cas_disagreements = 0
    consecutive_failures = 0

    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=output_fields)
        writer.writeheader()

        for i, row in enumerate(rows, 1):
            chebi_id = (row.get(chebi_column) or "").strip()
            out_row: dict[str, Any] = {
                chebi_column: chebi_id,
                **{c: row.get(c, "") for c in extra_columns},
                **{f: "" for f in OUTPUT_FIELDS_APPENDED},
            }

            if not chebi_id:
                log.warning("[%s/%s] Empty ChEBI ID.", i, total)
                out_row.update(empty_result())
                out_row["id_status"] = STATUS_MISSING
                status_counts[STATUS_MISSING] += 1
                writer.writerow(out_row)
                continue

            log.info("[%s/%s] ChEBI: %s", i, total, chebi_id)
            expected_name = (row.get(name_column) or "") if name_column else None
            expected_cas = (row.get(cas_column) or "") if cas_column else None

            parsed = normalize_chebi_id(chebi_id)
            if parsed is None:
                log.warning("  Not a ChEBI identifier: %s", chebi_id)
                out_row.update(empty_result())
                out_row["id_status"] = STATUS_INVALID
                status_counts[STATUS_INVALID] += 1
                invalid_values.add(chebi_id)
                writer.writerow(out_row)
                continue

            numeric, accession = parsed
            seen_ids.add(numeric)
            if numeric in payload_cache:
                payload = payload_cache[numeric]
            else:
                try:
                    payload = fetch_compound(numeric)
                except ChebiUnavailableError as exc:
                    # Never cached: a transient blip must not be replayed as
                    # not_found for every later row carrying this ID.
                    consecutive_failures += 1
                    log.error("  %s", exc)
                    out_row.update(empty_result())
                    out_row["id_status"] = STATUS_LOOKUP_FAILED
                    status_counts[STATUS_LOOKUP_FAILED] += 1
                    writer.writerow(out_row)
                    if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
                        raise ChebiTermsError(
                            f"ChEBI unreachable for {consecutive_failures} "
                            f"consecutive lookups — aborted at row {i} of {total}. "
                            f"Partial output: {output_path}"
                        ) from exc
                    continue
                payload_cache[numeric] = payload
                # Reset only on a lookup that actually reached ChEBI. A cache hit
                # is not evidence the outage is over, and resetting on one would
                # let a sheet with repeated IDs slip past the fail-fast guard.
                consecutive_failures = 0

            result = verify_payload(
                accession,
                payload,
                expected_name=expected_name,
                expected_cas=expected_cas,
                max_synonyms=max_synonyms,
            )
            out_row.update(result)
            writer.writerow(out_row)

            status = result["id_status"]
            status_counts[status] += 1
            if status == STATUS_OK:
                log.info("  %s — %s", result["chebi_accession"], result["chebi_name"])
            else:
                log.warning("  %s: %s", chebi_id, status)
            if result["name_verdict"] == NAME_MISMATCH:
                name_mismatches += 1
                log.warning(
                    "  Name mismatch: %r vs ChEBI %r",
                    (expected_name or "").strip(),
                    result["chebi_name"],
                )
            if result["cas_verdict"] == CAS_NOT_IN_CHEBI:
                cas_disagreements += 1
                log.warning(
                    "  CAS %s not recorded by ChEBI for %s",
                    (expected_cas or "").strip(),
                    result["chebi_accession"] or chebi_id,
                )

    log.info("=" * 55)
    log.info("Done. %s rows processed (%s unique IDs seen).", total, len(seen_ids))
    for status in SUMMARY_STATUSES:
        if status_counts[status]:
            log.info("  %-13s: %s", status, status_counts[status])
    if name_column:
        log.info("  Name mismatches : %s", name_mismatches)
    if cas_column:
        log.info("  CAS disagreements: %s", cas_disagreements)
    if status_counts[STATUS_LOOKUP_FAILED]:
        log.error(
            "ChEBI could not be reached for %s of %s rows — those rows say "
            "nothing about whether the compounds exist.",
            status_counts[STATUS_LOOKUP_FAILED],
            total,
        )
    # payload_cache holds one definitive answer per distinct ID that reached the
    # network: None for a 404, a payload otherwise. Failures were never cached.
    missed_ids = sum(1 for payload in payload_cache.values() if payload is None)
    resolved_ids = len(payload_cache) - missed_ids

    summary = RunSummary(
        status_counts=status_counts,
        missed_ids=missed_ids,
        resolved_ids=resolved_ids,
        invalid_values=len(invalid_values),
        missing_rows=status_counts[STATUS_MISSING],
        name_mismatches=name_mismatches,
        cas_disagreements=cas_disagreements,
    )

    if summary.suspect_endpoint:
        log.error(
            "All %s distinct IDs that reached ChEBI came back not_found and none "
            "resolved. Suspect a moved or renamed endpoint (%s) before trusting "
            "these rows.",
            missed_ids,
            # Read at call time so the message names the endpoint actually used.
            client.BASE,
        )
    if summary.suspect_column:
        log.error(
            "Not one of the %s rows held a usable ChEBI ID. Check that "
            "--chebi-column '%s' is the right column.",
            total,
            chebi_column,
        )
    log.info("Output: %s", output_path)
    return summary
