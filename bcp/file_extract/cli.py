from __future__ import annotations

import argparse
import logging
import sys
from typing import NoReturn

import boto3

from .constants import DEFAULT_H5_TARGET_FILENAME
from .fastq import (
    default_fastq_output_name,
    extract_fastq,
    r1_r2_mismatch_warning,
)
from .h5 import default_h5_output_name, extract_h5
from .h5_introspect import check_introspection_deps
from .s3_utils import parse_s3_uri

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


def _run_fastq(args: argparse.Namespace) -> int:
    location = parse_s3_uri(args.s3_uri)
    output = args.output or default_fastq_output_name(location.prefix)
    require_raw = not args.no_require_raw

    print(f"Bucket: {location.bucket}")
    print(f"Prefix: {location.prefix}")
    print(f"Require /raw/: {require_raw}")
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
    print(f"  Output: {output}")

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
        description="Extract S3 metadata for FASTQ.gz and Cell Ranger h5 deliverables.",
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

    if args.command == "fastq":
        return _run_fastq(args)
    if args.command == "h5":
        return _run_h5(args)
    parser.error(f"Unknown command: {args.command}")
    return 2
