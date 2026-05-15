from __future__ import annotations

import argparse
import sys


def main() -> None:
    """CLI for looking up ChEBI records from a file of identifiers."""
    parser = argparse.ArgumentParser(
        description="Look up ChEBI records for identifiers listed in an input file."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to input file containing identifiers (one per line or CSV/TSV)",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Path to write results (default: stdout)",
    )
    parser.add_argument(
        "--id-column",
        help="Column name when input is CSV/TSV (default: first column)",
    )
    args = parser.parse_args()

    # Implementation to be added.
    print(
        f"chebi_lookup: not implemented yet (input={args.input!r})",
        file=sys.stderr,
    )
    sys.exit(2)


if __name__ == "__main__":
    main()
