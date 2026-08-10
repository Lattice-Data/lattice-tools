from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .io import ChebiTermsError, emit_single_chebi_id, verify_chebi_file

log = logging.getLogger(__name__)

DEFAULT_CHEBI_COLUMN = "chebi_id"


def main() -> None:
    """CLI for resolving and verifying ChEBI IDs via the ChEBI REST API."""
    parser = argparse.ArgumentParser(
        description=(
            "Resolve ChEBI identifiers to their authoritative name and synonyms, "
            "and check them for correctness (current, released, and consistent "
            "with a recorded name or CAS number)."
        )
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--input",
        "-i",
        help="Input CSV file with a ChEBI ID column.",
    )
    mode.add_argument(
        "--chebi",
        help="Single ChEBI ID to look up, e.g. CHEBI:16236 or 16236 (no CSV required).",
    )
    parser.add_argument(
        "--chebi-column",
        default=DEFAULT_CHEBI_COLUMN,
        help=(
            "Name of the ChEBI ID column for batch mode "
            f"(default: '{DEFAULT_CHEBI_COLUMN}')."
        ),
    )
    parser.add_argument(
        "--name-column",
        default=None,
        help="Batch mode: column holding the recorded compound name to check.",
    )
    parser.add_argument(
        "--cas-column",
        default=None,
        help="Batch mode: column holding the recorded CAS number to check.",
    )
    parser.add_argument(
        "--expect-name",
        default=None,
        help="Single-ID mode: name to check against ChEBI's name and synonyms.",
    )
    parser.add_argument(
        "--expect-cas",
        default=None,
        help="Single-ID mode: CAS number to check against ChEBI's own CAS xrefs.",
    )
    parser.add_argument(
        "--max-synonyms",
        type=int,
        default=None,
        help="Cap the number of synonyms reported (default: no cap).",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help=(
            "Output path. Batch default: <input>_chebi_checked.csv. "
            "Single default: stdout."
        ),
    )
    parser.add_argument(
        "--format",
        choices=["json", "csv"],
        default="json",
        help="Output format for single-ID mode (default: json).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.max_synonyms is not None and args.max_synonyms < 1:
        log.error("--max-synonyms must be at least 1 (got %s).", args.max_synonyms)
        sys.exit(1)

    if args.chebi:
        for flag, value in (
            ("--chebi-column", args.chebi_column != DEFAULT_CHEBI_COLUMN),
            ("--name-column", args.name_column is not None),
            ("--cas-column", args.cas_column is not None),
        ):
            if value:
                log.warning(
                    "%s is ignored in single-ID mode; use --expect-name/--expect-cas.",
                    flag,
                )
        emit_single_chebi_id(
            args.chebi,
            Path(args.output) if args.output else None,
            fmt=args.format,
            expected_name=args.expect_name,
            expected_cas=args.expect_cas,
            max_synonyms=args.max_synonyms,
        )
        return

    for flag, value in (
        ("--expect-name", args.expect_name),
        ("--expect-cas", args.expect_cas),
    ):
        if value is not None:
            log.warning(
                "%s is ignored in batch mode; use --name-column/--cas-column.",
                flag,
            )

    input_path = Path(args.input)
    if args.output is None:
        output_path = input_path.parent / f"{input_path.stem}_chebi_checked.csv"
    else:
        output_path = Path(args.output)

    try:
        verify_chebi_file(
            input_path,
            args.chebi_column,
            output_path,
            name_column=args.name_column,
            cas_column=args.cas_column,
            max_synonyms=args.max_synonyms,
        )
    except ChebiTermsError as exc:
        log.error("%s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
