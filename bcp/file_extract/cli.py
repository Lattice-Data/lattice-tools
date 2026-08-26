from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import NoReturn

import boto3

from .constants import CRAM_SLOT_COLUMNS, DEFAULT_H5_TARGET_FILENAME
from .cram import (
    default_cram_output_name,
    extract_cram,
    only_unmatched_cram_warning,
    scan_cram_listing,
    ucram_found_warning,
)
from .fastq import (
    default_fastq_output_name,
    extract_fastq,
    r1_r2_mismatch_warning,
)
from .h5 import default_h5_output_name, extract_h5
from .h5_introspect import check_introspection_deps
from .s3_utils import parse_s3_uri
from .sheets import (
    LabIdentity,
    SheetBuildError,
    SheetOptions,
    cro_order_mismatch_warning,
    default_sheet_output_names,
    validate_cro_order,
)

log = logging.getLogger(__name__)


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def _print_failures(failures: list[tuple[str, str, str]], limit: int = 10) -> None:
    if not failures:
        return
    print(f"\nFailures (first {limit}):")
    for key, err1, err2 in failures[:limit]:
        print(f"  - {key}")
        if err1:
            print(f"      crc: {err1}")
        if err2:
            print(f"      detail: {err2}")
    if len(failures) > limit:
        print(f"  ... and {len(failures) - limit} more")


def _print_warnings(warnings: list[str], limit: int = 10) -> None:
    for warning in warnings[:limit]:
        print(f"  WARNING: {warning}")
    if len(warnings) > limit:
        print(f"  ... and {len(warnings) - limit} more warning(s)")


def _pilot_flag(value: str) -> bool:
    lowered = value.strip().lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    raise argparse.ArgumentTypeError(
        f"invalid value {value!r}: expected 'true' or 'false'"
    )


def _sheet_options(args: argparse.Namespace, *, output: str) -> SheetOptions:
    """Build the submission-sheet options from the required CLI inputs."""
    lab = LabIdentity.parse(args.lab, namespace=args.alias_namespace)
    cro_order = validate_cro_order(args.cro_order)
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = Path(output).parent
    sequence_file_path, sequence_file_set_path = default_sheet_output_names(
        cro_order, out_dir=out_dir
    )
    return SheetOptions(
        lab=lab,
        cro_order=cro_order,
        is_pilot_order=args.is_pilot_order,
        sequencing_platform=args.sequencing_platform,
        sequence_file_path=sequence_file_path,
        sequence_file_set_path=sequence_file_set_path,
        cram_slot=getattr(args, "cram_slot", None),
    )


def _print_sheet_outputs(options: SheetOptions) -> None:
    print(f"  SequenceFile sheet: {options.sequence_file_path}")
    print(f"  SequenceFileSet sheet: {options.sequence_file_set_path}")


def _run_fastq(args: argparse.Namespace) -> int:
    location = parse_s3_uri(args.s3_uri)
    output = args.output or default_fastq_output_name(location.prefix)
    require_raw = not args.no_require_raw
    sheet_options = _sheet_options(args, output=output)

    print(f"Bucket: {location.bucket}")
    print(f"Prefix: {location.prefix}")
    print(f"Require /raw/: {require_raw}")
    print(
        f"Lab: {sheet_options.lab.path} | alias namespace: "
        f"{sheet_options.lab.namespace}"
    )
    print(
        f"CRO order: {sheet_options.cro_order} | pilot: "
        f"{'TRUE' if sheet_options.is_pilot_order else 'FALSE'} | platform: "
        f"{sheet_options.sequencing_platform}"
    )
    order_warning = cro_order_mismatch_warning(sheet_options.cro_order, location.prefix)
    if order_warning:
        print(f"  WARNING: {order_warning}")
    print("Listing fastq.gz files ...")

    s3_client = boto3.client("s3")
    summary = extract_fastq(
        s3_client,
        location.bucket,
        location.prefix,
        output,
        require_raw=require_raw,
        workers=args.workers,
        retries=args.retries,
        show_progress=not args.quiet,
        sheets=sheet_options,
    )

    print(f"Found {summary.total} matching files")
    if summary.total == 0:
        print("Nothing to do. (If this order doesn't use /raw/, try --no-require-raw.)")
        return 0

    print(f"\nDone. Total: {summary.total}")
    print(
        f"  CRC64NVME retrieved: {summary.crc_ok} | failed: {summary.total - summary.crc_ok}"
    )
    print(
        f"  read_count retrieved: {summary.enrichment_ok} | "
        f"failed: {summary.total - summary.enrichment_ok}"
    )
    if summary.read_tally:
        tally = ", ".join(f"{k}={v}" for k, v in sorted(summary.read_tally.items()))
        print(f"  By read: {tally}")
    warning = r1_r2_mismatch_warning(summary.read_tally)
    if warning:
        print(f"  WARNING: {warning}")
    print(f"  SequenceFileSets: {summary.set_count}")
    print(f"  Output: {output}")
    _print_sheet_outputs(sheet_options)

    _print_warnings(summary.warnings)
    _print_failures(summary.failures)
    if args.strict and summary.has_failures:
        return 1
    return 0


def _run_cram(args: argparse.Namespace) -> int:
    location = parse_s3_uri(args.s3_uri)
    output = args.output or default_cram_output_name(location.prefix)
    require_raw = not args.no_require_raw
    sheet_options = _sheet_options(args, output=output)

    print(f"Bucket: {location.bucket}")
    print(f"Prefix: {location.prefix}")
    print(f"Require /raw/: {require_raw}")
    print(
        f"Lab: {sheet_options.lab.path} | alias namespace: "
        f"{sheet_options.lab.namespace}"
    )
    print(
        f"CRO order: {sheet_options.cro_order} | pilot: "
        f"{'TRUE' if sheet_options.is_pilot_order else 'FALSE'} | platform: "
        f"{sheet_options.sequencing_platform}"
    )
    print(f"CRAM slot: {CRAM_SLOT_COLUMNS[sheet_options.cram_slot]}")
    order_warning = cro_order_mismatch_warning(sheet_options.cro_order, location.prefix)
    if order_warning:
        print(f"  WARNING: {order_warning}")
    print("Listing cram files ...")

    s3_client = boto3.client("s3")
    # One traversal feeds both the guardrail warnings printed here and the work
    # extract_cram does, so the two cannot disagree about what the prefix holds.
    listing = scan_cram_listing(
        s3_client, location.bucket, location.prefix, require_raw=require_raw
    )

    ucram_warning = ucram_found_warning(listing.warnings.ucram_count)
    if ucram_warning:
        print(f"  WARNING: {ucram_warning}")

    print(f"Found {len(listing.targets)} matching files")
    if not listing.targets:
        print("Nothing to do. (If this order doesn't use /raw/, try --no-require-raw.)")
        unmatched_warning = only_unmatched_cram_warning(
            listing.warnings.unmatched_cram_count,
            require_raw=require_raw,
        )
        if unmatched_warning:
            print(f"  WARNING: {unmatched_warning}")
        return 0

    summary = extract_cram(
        s3_client,
        location.bucket,
        location.prefix,
        output,
        require_raw=require_raw,
        workers=args.workers,
        retries=args.retries,
        show_progress=not args.quiet,
        sheets=sheet_options,
        listing=listing,
    )

    print(f"\nDone. Total: {summary.total}")
    print(
        f"  CRC64NVME retrieved: {summary.crc_ok} | failed: {summary.total - summary.crc_ok}"
    )
    print(
        f"  read_count retrieved: {summary.enrichment_ok} | "
        f"failed: {summary.total - summary.enrichment_ok}"
    )
    print(f"  SequenceFileSets: {summary.set_count}")
    print(f"  Output: {output}")
    _print_sheet_outputs(sheet_options)

    _print_warnings(summary.warnings)
    _print_failures(summary.failures)
    if args.strict and summary.has_failures:
        return 1
    return 0


def _run_h5(args: argparse.Namespace) -> int:
    location = parse_s3_uri(args.s3_uri)
    output = args.output or default_h5_output_name(location.prefix)
    do_introspect = not args.no_introspect

    if do_introspect:
        check_introspection_deps()

    print(f"Bucket: {location.bucket}")
    print(f"Prefix: {location.prefix}")
    print(f"Target filename: {args.target_filename}")
    print(
        f"Introspect: {do_introspect} | genome: {args.genome} | metrics: {args.metrics}"
    )
    print("Listing matching h5 files ...")

    s3_client = boto3.client("s3")
    summary = extract_h5(
        s3_client,
        location.bucket,
        location.prefix,
        output,
        target_filename=args.target_filename,
        do_introspect=do_introspect,
        do_genome=args.genome,
        do_metrics=args.metrics,
        workers=args.workers,
        retries=args.retries,
        show_progress=not args.quiet,
    )

    print(f"Found {summary.total} matching files")
    if summary.total == 0:
        print("Nothing to do.")
        return 0

    print(f"\nDone. Total: {summary.total} | checksum OK: {summary.crc_ok}", end="")
    if do_introspect:
        print(f" | introspect OK: {summary.enrichment_ok}", end="")
    print(f"\nOutput: {output}")

    _print_failures(summary.failures)
    if args.strict and summary.has_failures:
        return 1
    return 0


def _invalid_uri_exit(message: str) -> NoReturn:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def _add_sheet_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the inputs the Lattice submission sheets cannot derive from S3."""
    parser.add_argument(
        "--lab",
        required=True,
        help="Lab name or path, e.g. example-lab or /labs/example-lab/",
    )
    parser.add_argument(
        "--cro-order",
        required=True,
        help="CRO order identifier, e.g. NVUS0000000000-15",
    )
    parser.add_argument(
        "--is-pilot-order",
        required=True,
        type=_pilot_flag,
        metavar="{true,false}",
        help="Whether this order is a pilot (written to SequenceFileSet)",
    )
    parser.add_argument(
        "--sequencing-platform",
        required=True,
        help='Sequencing platform, e.g. "Ultima Genomics UG 100"',
    )
    parser.add_argument(
        "--alias-namespace",
        default=None,
        help="Alias namespace (default: the lab name)",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Directory for the sheet TSVs (default: alongside --output)",
    )


def build_parser() -> argparse.ArgumentParser:
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parent.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Disable progress bars.",
    )

    parser = argparse.ArgumentParser(
        description=(
            "Extract S3 metadata for FASTQ.gz, CRAM, and Cell Ranger h5 deliverables."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[parent],
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    fastq = subparsers.add_parser(
        "fastq",
        help="Extract FASTQ.gz metadata under an S3 order prefix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[parent],
    )
    fastq.add_argument("s3_uri", help="s3://bucket/path/order")
    fastq.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV (default: <order>_fastq_info.tsv)",
    )
    fastq.add_argument(
        "--no-require-raw",
        action="store_true",
        help="Don't require a /raw/ subfolder in the key",
    )
    fastq.add_argument("--workers", type=int, default=None, help="Process count")
    fastq.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Max attempts per transient operation (default: 5)",
    )
    fastq.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any per-file enrichment fails",
    )
    _add_sheet_arguments(fastq)

    cram = subparsers.add_parser(
        "cram",
        help="Extract CRAM metadata under an S3 order prefix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[parent],
    )
    cram.add_argument("s3_uri", help="s3://bucket/path/order")
    cram.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV (default: <order>_cram_info.tsv)",
    )
    cram.add_argument(
        "--no-require-raw",
        action="store_true",
        help="Don't require a /raw/ subfolder in the key",
    )
    cram.add_argument("--workers", type=int, default=None, help="Process count")
    cram.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Max attempts per transient operation (default: 5)",
    )
    cram.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any per-file enrichment fails",
    )
    _add_sheet_arguments(cram)
    cram.add_argument(
        "--cram-slot",
        required=True,
        choices=sorted(CRAM_SLOT_COLUMNS),
        help="SequenceFileSet slot the deliverable CRAM fills",
    )

    h5 = subparsers.add_parser(
        "h5",
        help="Extract Cell Ranger h5 matrix metadata.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[parent],
    )
    h5.add_argument("s3_uri", help="s3://bucket/.../outs/per_sample_outs")
    h5.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV (default: <run-or-dir>_h5_info.tsv)",
    )
    h5.add_argument(
        "--target-filename",
        default=DEFAULT_H5_TARGET_FILENAME,
        help=f"h5 basename to match (default: {DEFAULT_H5_TARGET_FILENAME})",
    )
    h5.add_argument(
        "--no-introspect",
        action="store_true",
        help="Skip opening the h5; emit checksums/size only",
    )
    h5.add_argument(
        "--genome",
        action="store_true",
        help="Add per-genome Gene Expression counts",
    )
    h5.add_argument(
        "--metrics",
        action="store_true",
        help="Cross-check cell count against sibling metrics_summary.csv",
    )
    h5.add_argument("--workers", type=int, default=None, help="Thread count")
    h5.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Max attempts per transient operation (default: 5)",
    )
    h5.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any per-file enrichment fails",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)

    try:
        parse_s3_uri(args.s3_uri)
    except ValueError as exc:
        _invalid_uri_exit(str(exc))

    try:
        if args.command == "fastq":
            return _run_fastq(args)
        if args.command == "cram":
            return _run_cram(args)
        if args.command == "h5":
            return _run_h5(args)
    except SheetBuildError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    parser.error(f"Unknown command: {args.command}")
    return 2
