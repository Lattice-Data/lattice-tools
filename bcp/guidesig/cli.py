"""CLI for protospacer set signatures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .core import (
    DEFAULT_COLUMN,
    FILE_FORMATS,
    GuideSigError,
    protospacer_set,
    protospacers_and_guide_ids,
    signature,
    signature_of,
)

# Matches file_extract's DEFAULT_H5_TARGET_FILENAME; spelled out here only for
# the help text, so building the parser never imports the S3/h5 stack.
H5_FILENAME_HELP_DEFAULT = "sample_filtered_feature_bc_matrix.h5"
DEFAULT_RETRIES = 5


def _positive_int(value: str) -> int:
    """argparse type for a count of at least 1."""
    try:
        parsed = int(value)
    except ValueError:
        parsed = 0
    if parsed < 1:
        raise argparse.ArgumentTypeError(
            f"invalid value {value!r}: expected a positive integer"
        )
    return parsed


def _format_examples(sequences: set[str], limit: int = 10) -> str:
    examples = sorted(sequences)[:limit]
    return ", ".join(examples) if examples else "(none)"


def _compare(file_a: Path, file_b: Path, column: str, file_format: str) -> int:
    set_a = protospacer_set(file_a, file_format=file_format, column=column)
    set_b = protospacer_set(file_b, file_format=file_format, column=column)

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
            "guide-template TSV or CSV files."
        ),
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="One or more guide-template files to signature",
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("FILE_A", "FILE_B"),
        type=Path,
        help=(
            "Compare two files and report match or set-difference summary; "
            "both files must be in the format given by --format"
        ),
    )
    parser.add_argument(
        "--format",
        dest="file_format",
        required=True,
        choices=sorted(FILE_FORMATS),
        help="Input file format; required, never inferred from the file name",
    )
    parser.add_argument(
        "--column",
        default=DEFAULT_COLUMN,
        metavar="NAME",
        help=f"Protospacer column name (default: {DEFAULT_COLUMN})",
    )
    parser.add_argument(
        "--h5-uri",
        metavar="S3_URI",
        help=(
            "Also check that every guide_id in the one FILE appears in each "
            "sample's h5 under this S3 order prefix (s3://bucket/proj/order); "
            "not allowed with --compare"
        ),
    )
    parser.add_argument(
        "--h5-filename",
        metavar="NAME",
        help=f"h5 basename to check (default: {H5_FILENAME_HELP_DEFAULT})",
    )
    parser.add_argument(
        "--workers",
        type=_positive_int,
        help="Threads for reading h5 files (default: 8)",
    )
    parser.add_argument(
        "--retries",
        type=_positive_int,
        help=f"Max attempts per transient S3 error (default: {DEFAULT_RETRIES})",
    )
    return parser


def _run_h5_check(path: Path, args: argparse.Namespace) -> None:
    """Signature one file, then check its guide_ids against the order's h5s.

    The sheet's own checks run first; if they fail there is no trustworthy
    guide_id list, so S3 is never touched. The signature is printed only if
    every h5 passes.
    """
    try:
        sequences, guide_ids = protospacers_and_guide_ids(
            path, file_format=args.file_format, column=args.column
        )
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    try:
        from . import h5check
    except ImportError as exc:
        print(
            f"error: --h5-uri needs the file_extract package and its "
            f"dependencies (missing: {exc.name}); run from bcp/ with "
            "requirements.txt installed",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    options: dict[str, object] = {"workers": args.workers}
    if args.h5_filename is not None:
        options["filename"] = args.h5_filename
    if args.retries is not None:
        options["retries"] = args.retries
    try:
        report = h5check.check_guide_ids_against_h5(guide_ids, args.h5_uri, **options)
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    for note in report.notes:
        print(f"note: {note}", file=sys.stderr)
    if not report.ok:
        for failure in report.failures:
            print(f"error: {failure}", file=sys.stderr)
        print(
            f"error: guide_id check failed: {len(report.failures)} problem(s); "
            f"{report.h5_failed} of {report.h5_total} h5 file(s) failed",
            file=sys.stderr,
        )
        raise SystemExit(1)

    print(f"{signature_of(sequences)}\t{path}")


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.compare is not None and args.files:
        parser.error("do not pass FILE arguments with --compare")
    if args.compare is None and not args.files:
        parser.error("provide FILE arguments, or use --compare FILE_A FILE_B")

    h5_only = [
        flag
        for flag, value in (
            ("--h5-filename", args.h5_filename),
            ("--workers", args.workers),
            ("--retries", args.retries),
        )
        if value is not None
    ]
    if args.h5_uri is None and h5_only:
        parser.error(f"{', '.join(h5_only)} require --h5-uri")
    if args.h5_uri is not None:
        if args.compare is not None:
            parser.error("--h5-uri is not allowed with --compare")
        if len(args.files) != 1:
            parser.error("--h5-uri checks exactly one FILE")
        _run_h5_check(args.files[0], args)
        return

    try:
        if args.compare is not None:
            file_a, file_b = args.compare
            raise SystemExit(_compare(file_a, file_b, args.column, args.file_format))

        # Compute every signature before printing so a failure on a later file
        # does not leave partial output on stdout for a piped consumer.
        lines = []
        for path in args.files:
            sig = signature(path, file_format=args.file_format, column=args.column)
            lines.append(f"{sig}\t{path}")
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    print("\n".join(lines))


if __name__ == "__main__":
    main()
