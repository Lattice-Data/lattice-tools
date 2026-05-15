"""CSV batch I/O for CAS → ChEBI mapping."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from .client import OUTPUT_FIELDS_APPENDED, lookup_cas

log = logging.getLogger(__name__)


class CasMappingError(Exception):
    """Raised when input file or column validation fails."""


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
