"""CLI for protospacer set signatures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .core import DEFAULT_COLUMN, GuideSigError, protospacer_set, signature


def _format_examples(sequences: set[str], limit: int = 10) -> str:
    examples = sorted(sequences)[:limit]
    return ", ".join(examples) if examples else "(none)"


def _compare(file_a: Path, file_b: Path, column: str) -> int:
    set_a = protospacer_set(file_a, column=column)
    set_b = protospacer_set(file_b, column=column)

    if set_a == set_b:
        print("match")
        return 0

    only_a = set_a - set_b
    only_b = set_b - set_a
    intersection = set_a & set_b

    print("mismatch")
    print(
        f"|A|={len(set_a)} |B|={len(set_b)} "
        f"|A n B|={len(intersection)} "
        f"|A - B|={len(only_a)} |B - A|={len(only_b)}"
    )
    print(f"A - B examples: {_format_examples(only_a)}")
    print(f"B - A examples: {_format_examples(only_b)}")
    return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="guidesig",
        description=(
            "Compute deterministic protospacer-set signatures for "
            "guide-template TSV files."
        ),
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="One or more guide-template TSV files to signature",
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("FILE_A", "FILE_B"),
        type=Path,
        help="Compare two files and report match or set-difference summary",
    )
    parser.add_argument(
        "--column",
        default=DEFAULT_COLUMN,
        metavar="NAME",
        help=f"Protospacer column name (default: {DEFAULT_COLUMN})",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.compare is not None and args.files:
        parser.error("do not pass FILE arguments with --compare")
    if args.compare is None and not args.files:
        parser.error("provide FILE arguments, or use --compare FILE_A FILE_B")

    try:
        if args.compare is not None:
            file_a, file_b = args.compare
            raise SystemExit(_compare(file_a, file_b, args.column))

        for path in args.files:
            print(f"{signature(path, column=args.column)}\t{path}")
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
