from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .io import StructureCheckError, check_file, emit_single_row

log = logging.getLogger(__name__)


def main() -> None:
    """CLI for cross-checking curation rows by structure rather than by name."""
    parser = argparse.ArgumentParser(
        description=(
            "Cross-check a curation sheet by structure: resolve each row's CAS "
            "number, ChEBI ID, and compound name to an InChIKey independently, "
            "then compare. Answers 'has ChEBI got something else under this ID?' "
            "without relying on how anyone spells the compound."
        )
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", "-i", help="Input CSV file.")
    mode.add_argument("--cas", help="Single CAS number to check (no CSV required).")

    parser.add_argument(
        "--cas-column",
        default="CAS",
        help="Batch mode: CAS number column, the pivot for both checks "
        "(default: 'CAS').",
    )
    parser.add_argument(
        "--chebi-column",
        default=None,
        help="Batch mode: ChEBI ID column, to check the ID against the CAS.",
    )
    parser.add_argument(
        "--name-column",
        default=None,
        help="Batch mode: compound name column, to check the name against the CAS.",
    )
    parser.add_argument(
        "--chebi",
        default=None,
        help="Single mode: ChEBI ID to check against --cas.",
    )
    parser.add_argument(
        "--name",
        default=None,
        help="Single mode: compound name to check against --cas.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output path. Batch default: <input>_structure_checked.csv. "
        "Single default: stdout.",
    )
    parser.add_argument(
        "--format",
        choices=["json", "csv"],
        default="json",
        help="Output format for single mode (default: json).",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging."
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.cas is not None:
        if not args.chebi and not args.name:
            log.error("Pass --chebi and/or --name to check the CAS against.")
            sys.exit(1)
        try:
            emit_single_row(
                Path(args.output) if args.output else None,
                name=args.name,
                cas=args.cas,
                chebi_id=args.chebi,
                fmt=args.format,
            )
        except StructureCheckError as exc:
            log.error("%s", exc)
            sys.exit(1)
        return

    for flag, value in (("--chebi", args.chebi), ("--name", args.name)):
        if value is not None:
            log.warning(
                "%s is ignored in batch mode; use --chebi-column/--name-column.", flag
            )
    if args.format != "json":
        log.warning("--format is ignored in batch mode; output is always CSV.")

    input_path = Path(args.input)
    output_path = (
        Path(args.output)
        if args.output
        else input_path.parent / f"{input_path.stem}_structure_checked.csv"
    )

    try:
        summary = check_file(
            input_path,
            output_path,
            cas_column=args.cas_column,
            chebi_column=args.chebi_column,
            name_column=args.name_column,
        )
    except StructureCheckError as exc:
        log.error("%s", exc)
        sys.exit(1)

    # A run where nothing was compared is a broken run, not a clean sheet, and must
    # not be mistaken for one by a wrapper script.
    if summary.degraded:
        sys.exit(1)

    # Findings themselves are the product, not a failure: exit 0 and let the review
    # ranking drive what happens next.
    if summary.needs_attention:
        log.info(
            "%s row(s) need attention — sort the output by 'review_rank'.",
            summary.needs_attention,
        )


if __name__ == "__main__":
    main()
