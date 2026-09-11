"""
Dash + Cytoscape explorer for the DB2 reference graph.

Seeded from a single object, the graph grows only where the user clicks: each
tap resolves that node's neighbors and folds them in. Nothing is precomputed,
so the seed can sit anywhere in the database.

The seed can be an object path, an alias or a bare uuid; whichever it is, it is
resolved against the instance and the graph is keyed by the canonical
/object_type/uuid/ that comes back.

Run from the bcp directory in an environment with
    db2_flattener installed
    DB2_<MODE>_KEY / _SECRET / _SERVER set:

    python -m graph_db2
    python -m graph_db2 --mode db2_demo --seed /matrix_file_sets/<uuid>/
    python -m graph_db2 --mode db2_demo --seed some-lab:my_matrixfileset

A seed only has meaning on the server it came from, so --mode without --seed
starts on an empty canvas rather than guessing a path from another deployment.
"""

from __future__ import annotations

import argparse

from .constants import DEFAULT_MODE, DEFAULT_PORT, FETCH_NEW
from .explorer import SAMPLE_SEED, build_app


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="graph_db2",
        description=__doc__,
    )
    parser.add_argument(
        "--seed",
        help="object path, alias or uuid to start from; defaults to a sample on "
        f"{DEFAULT_MODE} and to an empty canvas elsewhere",
    )
    parser.add_argument(
        "--mode",
        default=DEFAULT_MODE,
        help=f"instance to run on, defaults to {DEFAULT_MODE}",
    )
    parser.add_argument(
        "--fetch-new",
        action="store_true",
        default=FETCH_NEW,
        help="Fetch new profile schemas from DB2 instance instead of loading the constants.yaml",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Set port, default is {DEFAULT_PORT}",
    )
    parser.add_argument("--debug", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    # SAMPLE_SEED only exists on DEFAULT_MODE's server; carrying it over to
    # another deployment just 404s into an empty canvas
    seed = args.seed or (SAMPLE_SEED if args.mode == DEFAULT_MODE else "")
    build_app(seed, args.mode, args.fetch_new).run(debug=args.debug, port=args.port)


if __name__ == "__main__":
    main()
