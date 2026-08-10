"""CSV batch I/O and single-ID output for ChEBI ID verification."""

from __future__ import annotations

import csv
import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Literal

from .client import (
    CAS_NOT_IN_CHEBI,
    NAME_MISMATCH,
    NOT_CHECKED,
    OUTPUT_FIELDS_APPENDED,
    STATUS_INVALID,
    STATUS_NOT_FOUND,
    STATUS_NOT_RELEASED,
    STATUS_OK,
    STATUS_SECONDARY,
    empty_result,
    fetch_compound,
    normalize_chebi_id,
    verify_chebi_id,
    verify_payload,
)

log = logging.getLogger(__name__)


class ChebiTermsError(Exception):
    """Raised when input file or column validation fails."""


SINGLE_CHEBI_COLUMN = "chebi_id"
SINGLE_OUTPUT_FIELDS = [SINGLE_CHEBI_COLUMN, *OUTPUT_FIELDS_APPENDED]


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
) -> None:
    """Resolve one ChEBI ID and write JSON or CSV to a file or stdout."""
    chebi_id = chebi_id.strip()
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
        return

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
) -> None:
    """Read a CSV of ChEBI IDs, verify each against ChEBI, and write enriched output."""
    if not input_path.exists():
        raise ChebiTermsError(f"Input file not found: {input_path}")

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
    # and per-row expectations are evaluated against the cached payload.
    payload_cache: dict[str, dict[str, Any] | None] = {}
    status_counts: Counter[str] = Counter()
    name_mismatches = 0
    cas_disagreements = 0

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
                log.warning("[%s/%s] Empty ChEBI ID, skipping.", i, total)
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
                writer.writerow(out_row)
                continue

            numeric, accession = parsed
            if numeric not in payload_cache:
                payload_cache[numeric] = fetch_compound(numeric)
            payload = payload_cache[numeric]

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
    log.info("Done. %s rows processed (%s unique IDs).", total, len(payload_cache))
    for status in (
        STATUS_OK,
        STATUS_SECONDARY,
        STATUS_NOT_RELEASED,
        STATUS_NOT_FOUND,
        STATUS_INVALID,
    ):
        if status_counts[status]:
            log.info("  %-13s: %s", status, status_counts[status])
    if name_column:
        log.info("  Name mismatches : %s", name_mismatches)
    if cas_column:
        log.info("  CAS disagreements: %s", cas_disagreements)
    log.info("  Output: %s", output_path)
