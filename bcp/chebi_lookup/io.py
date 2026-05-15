"""CSV batch I/O and single-CAS output for CAS → ChEBI mapping."""

from __future__ import annotations

import csv
import json
import logging
import sys
from pathlib import Path
from typing import Any, Literal

from .client import OUTPUT_FIELDS_APPENDED, lookup_cas

log = logging.getLogger(__name__)


class CasMappingError(Exception):
    """Raised when input file or column validation fails."""


SINGLE_CAS_COLUMN = "CAS"
SINGLE_OUTPUT_FIELDS = [SINGLE_CAS_COLUMN, *OUTPUT_FIELDS_APPENDED]


def build_single_cas_row(cas: str, result: dict[str, Any]) -> dict[str, Any]:
    """Build one output row for a single-CAS lookup."""
    row: dict[str, Any] = {SINGLE_CAS_COLUMN: cas}
    for field in OUTPUT_FIELDS_APPENDED:
        value = result.get(field, "")
        row[field] = "" if value is None else value
    return row


def _log_single_cas_summary(cas: str, result: dict[str, Any]) -> None:
    cid = result.get("pubchem_cid", "")
    if not cid:
        log.warning("No PubChem CID for CAS: %s", cas)
        return
    log.info("CID: %s", cid)
    chebi_id = result.get("chebi_id", "")
    if chebi_id:
        log.info("ChEBI: %s", chebi_id)
    else:
        log.info("No ChEBI xref for CID %s", cid)


def emit_single_cas(
    cas: str,
    output_path: Path | None,
    *,
    fmt: Literal["json", "csv"] = "json",
) -> None:
    """Resolve one CAS and write JSON or CSV to a file or stdout."""
    cas = cas.strip()
    result = lookup_cas(cas)
    row = build_single_cas_row(cas, result)
    _log_single_cas_summary(cas, result)

    if fmt == "json":
        payload = json.dumps(row, indent=2)
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


def map_cas_file(
    input_path: Path,
    cas_column: str,
    output_path: Path,
) -> None:
    """Read a CSV, resolve each CAS via PubChem, and write enriched output."""
    if not input_path.exists():
        raise CasMappingError(f"Input file not found: {input_path}")

    with open(input_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        if cas_column not in fieldnames:
            raise CasMappingError(
                f"Column '{cas_column}' not found. Available columns: {fieldnames}"
            )
        rows = list(reader)
        extra_columns = [c for c in fieldnames if c != cas_column]

    output_fields = [cas_column] + extra_columns + OUTPUT_FIELDS_APPENDED
    total = len(rows)
    found_cid = 0
    found_chebi = 0

    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=output_fields)
        writer.writeheader()

        for i, row in enumerate(rows, 1):
            cas = row[cas_column].strip()
            log.info("[%s/%s] CAS: %s", i, total, cas)

            out_row = {
                cas_column: cas,
                **{c: row[c] for c in extra_columns},
                **{f: "" for f in OUTPUT_FIELDS_APPENDED},
            }

            if not cas:
                log.warning("Empty CAS at row %s, skipping.", i)
                writer.writerow(out_row)
                continue

            props = lookup_cas(cas)
            cid = props.get("pubchem_cid", "")
            if not cid:
                log.warning("No PubChem CID for CAS: %s", cas)
                writer.writerow(out_row)
                continue

            found_cid += 1
            log.info("  CID: %s", cid)

            out_row.update(props)

            if props.get("chebi_id"):
                found_chebi += 1
                log.info("  ChEBI: %s", props["chebi_id"])
            else:
                log.info("No ChEBI xref for CID %s", cid)

            writer.writerow(out_row)

    log.info("=" * 55)
    log.info("Done. %s CAS numbers processed.", total)
    log.info("  Resolved to PubChem CID : %s/%s", found_cid, total)
    log.info("  Resolved to ChEBI ID    : %s/%s", found_chebi, total)
    log.info("  Output: %s", output_path)
