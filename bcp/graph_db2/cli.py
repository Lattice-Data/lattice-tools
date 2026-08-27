from __future__ import annotations

import argparse

from .constants import DEFAULT_MODE, FETCH_NEW
from .explorer import SAMPLE_SEED, build_app


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="graph_db2",
        description=__doc__,
    )
    parser.add_argument(
        "--seed",
        help=f"object path to start from; defaults to a sample on {DEFAULT_MODE} "
        "and to an empty canvas elsewhere",
    )
    parser.add_argument("--mode", default=DEFAULT_MODE)
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--fetch-new", action="store_true", default=FETCH_NEW)

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
