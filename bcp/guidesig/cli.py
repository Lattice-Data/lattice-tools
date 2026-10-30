"""CLI for protospacer set signatures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .core import (
    DEFAULT_COLUMN,
    FILE_FORMATS,
    GuideSigError,
    header_warnings,
    protospacer_set,
    protospacers_and_guide_ids,
    signature,
    signature_of,
)

# Mirrors h5check's scope names and default basenames; spelled out here only
# for the help text and argument validation, so building the parser never
# imports the S3/h5 stack.
H5_SCOPE_MULTI = "multi"
H5_SCOPE_PER_SAMPLE = "per-sample"
H5_SCOPES = (H5_SCOPE_MULTI, H5_SCOPE_PER_SAMPLE)
H5_MULTI_FILENAME_HELP_DEFAULT = "raw_feature_bc_matrix.h5"
H5_PER_SAMPLE_FILENAME_HELP_DEFAULT = "sample_filtered_feature_bc_matrix.h5"
# Mirror h5check's per-scope worker defaults, for the help text only.
H5_MULTI_WORKERS_HELP_DEFAULT = 24
H5_PER_SAMPLE_WORKERS_HELP_DEFAULT = 8
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
            "Also check that every guide_id in the one FILE appears in the h5 "
            "files under this S3 order prefix (s3://bucket/proj/order); "
            "not allowed with --compare"
        ),
    )
    parser.add_argument(
        "--h5-scope",
        choices=list(H5_SCOPES),
        # Left unset so main() can tell an explicit --h5-scope from the default
        # and still require --h5-uri alongside it.
        default=None,
        help=(
            f"Which h5 files to check (default: {H5_SCOPE_MULTI}). "
            f"{H5_SCOPE_MULTI!r} reads the run-level "
            f"{H5_MULTI_FILENAME_HELP_DEFAULT}, one per library, from either "
            "outs/ (Cell Ranger 10.x) or outs/multi/count/ (9.0.1); that is "
            "the whole guide_id check because a cellranger multi run's feature "
            "reference is per-run, not per-sample; a library with none fails, "
            f"and one with several fails as an ambiguity. "
            f"{H5_SCOPE_PER_SAMPLE!r} instead sweeps every "
            f"per_sample_outs/<sample>/count/"
            f"{H5_PER_SAMPLE_FILENAME_HELP_DEFAULT}: since those all carry the "
            "same features table, it is a file integrity check (missing, "
            "truncated or unreadable per-sample files), not a guide check"
        ),
    )
    parser.add_argument(
        "--h5-filename",
        metavar="NAME",
        help=(
            "h5 basename to check, overriding the scope default "
            f"({H5_MULTI_FILENAME_HELP_DEFAULT} for --h5-scope "
            f"{H5_SCOPE_MULTI}, {H5_PER_SAMPLE_FILENAME_HELP_DEFAULT} for "
            f"--h5-scope {H5_SCOPE_PER_SAMPLE})"
        ),
    )
    parser.add_argument(
        "--h5-samples-file",
        type=Path,
        metavar="PATH",
        help=(
            "Only check the names in PATH (one exact name per line, no "
            "header). Use it when the sheet "
            "describes only some of the order's libraries, so the rest are "
            "excluded instead of failing for guides they were never built "
            f"with. Under --h5-scope {H5_SCOPE_MULTI} the names are library "
            "folders (the folder before processed/); under "
            f"{H5_SCOPE_PER_SAMPLE} they match per_sample_outs/<sample>/, "
            "falling back to the library folder. A listed name with no h5 fails"
        ),
    )
    parser.add_argument(
        "--workers",
        type=_positive_int,
        help=(
            "Threads for reading h5 files, and the S3 connection pool they "
            f"share (default: {H5_MULTI_WORKERS_HELP_DEFAULT} for --h5-scope "
            f"{H5_SCOPE_MULTI}, which opens one h5 per library; "
            f"{H5_PER_SAMPLE_WORKERS_HELP_DEFAULT} for "
            f"{H5_SCOPE_PER_SAMPLE}, which opens one per sample)"
        ),
    )
    parser.add_argument(
        "--retries",
        type=_positive_int,
        help=f"Max attempts per transient S3 error (default: {DEFAULT_RETRIES})",
    )
    parser.add_argument(
        "--guidescan",
        action="store_true",
        help=(
            "Also run guidescan_pipeline on the one FILE's targeting guides and "
            "report the exact-match rate (warning below 90%%); runs after "
            "--h5-uri; needs --species; not allowed with --compare"
        ),
    )
    parser.add_argument(
        "--species",
        choices=["human", "mouse"],
        help="Genome for --guidescan; cross-checked against the sheet's gene IDs",
    )
    return parser


def _print_header_warnings(paths: list[Path], file_format: str) -> None:
    """Print a warning per recommended column a file lacks; never fails."""
    for path in paths:
        for warning in header_warnings(path, file_format=file_format):
            print(f"warning: {warning}", file=sys.stderr)


def _run_single_file_checks(path: Path, args: argparse.Namespace) -> None:
    """Signature one file after its --h5-uri and/or --guidescan checks.

    The sheet's own checks run first; if they fail there is no trustworthy
    guide list, so S3 is never touched. A failed h5 check stops the run before
    guidescan. The guidescan rate is reported but never fails the run. The
    guidescan species cross-check and setup errors do fail the run, as do
    sheet-content errors from read_targeting_guides: a missing guide_PAM, an
    unparseable or 5' PAM, mixed PAM patterns, and a targeting guide with no
    guide_id or one guide_id carrying two different protospacers.
    """
    try:
        sequences, guide_ids = protospacers_and_guide_ids(
            path, file_format=args.file_format, column=args.column
        )
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    if args.h5_uri is not None:
        _run_h5_check(guide_ids, args)
    if args.guidescan:
        _run_guidescan_check(path, args)

    print(f"{signature_of(sequences)}\t{path}")


def _run_guidescan_check(path: Path, args: argparse.Namespace) -> None:
    """Report the sheet's exact-match rate; exit 1 only if it cannot run."""
    from . import guidescancheck

    try:
        guides, target_ids = guidescancheck.read_targeting_guides(
            path, file_format=args.file_format, column=args.column
        )
        species_note = guidescancheck.check_species(target_ids, args.species)
        if species_note:
            print(f"note: {species_note}", file=sys.stderr)
        index_prefix = None
        if guides:
            index_prefix = guidescancheck.ensure_references(
                args.species, log=guidescancheck.stderr_log
            )
        report = guidescancheck.run_guidescan_check(guides, index_prefix)
    except GuideSigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    for line in guidescancheck.format_report(report):
        print(line, file=sys.stderr)


def _run_h5_check(guide_ids: set[str], args: argparse.Namespace) -> None:
    """Check the sheet's guide_ids against the order's h5s; exit 1 on failure.

    The scope decides which h5 files those are; see ``--h5-scope``. An unset
    ``--h5-scope`` means the default, resolved here rather than by argparse so
    main() can still require ``--h5-uri`` alongside an explicit scope.
    """
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

    options: dict[str, object] = {
        "workers": args.workers,
        "scope": args.h5_scope or H5_SCOPE_MULTI,
    }
    if args.h5_filename is not None:
        options["filename"] = args.h5_filename
    if args.retries is not None:
        options["retries"] = args.retries
    try:
        if args.h5_samples_file is not None:
            options["samples"] = h5check.read_sample_names(args.h5_samples_file)
            options["samples_label"] = str(args.h5_samples_file)
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
            ("--h5-scope", args.h5_scope),
            ("--h5-filename", args.h5_filename),
            ("--h5-samples-file", args.h5_samples_file),
            ("--workers", args.workers),
            ("--retries", args.retries),
        )
        if value is not None
    ]
    if args.h5_uri is None and h5_only:
        parser.error(f"{', '.join(h5_only)} require --h5-uri")

    if args.guidescan and args.species is None:
        parser.error("--guidescan requires --species")
    if args.species is not None and not args.guidescan:
        parser.error("--species requires --guidescan")
    for flag, used in (
        ("--h5-uri", args.h5_uri is not None),
        ("--guidescan", args.guidescan),
    ):
        if not used:
            continue
        if args.compare is not None:
            parser.error(f"{flag} is not allowed with --compare")
        if len(args.files) != 1:
            parser.error(f"{flag} checks exactly one FILE")
    _print_header_warnings(
        list(args.compare) if args.compare is not None else args.files,
        args.file_format,
    )
    if args.h5_uri is not None or args.guidescan:
        _run_single_file_checks(args.files[0], args)
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
