from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .io import CasMappingError, map_cas_file

log = logging.getLogger(__name__)


def main() -> None:
    """CLI for mapping CAS numbers to ChEBI + PubChem properties via PubChem PUG REST."""
    parser = argparse.ArgumentParser(
        description=(
            "Map CAS Registry Numbers to ChEBI identifiers and PubChem properties "
            "via the PubChem PUG REST API."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input CSV file with a CAS column.",
    )
    parser.add_argument(
        "--cas-column",
        "-c",
        default="CAS",
        help="Name of the CAS column (default: 'CAS').",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output CSV (default: <input>_chebi_mapped.csv).",
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

    input_path = Path(args.input)
    if args.output is None:
        output_path = input_path.parent / f"{input_path.stem}_chebi_mapped.csv"
    else:
        output_path = Path(args.output)

    try:
        map_cas_file(input_path, args.cas_column, output_path)
    except CasMappingError as exc:
        log.error("%s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
