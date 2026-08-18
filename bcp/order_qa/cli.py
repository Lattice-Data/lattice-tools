"""Command line for order QA.

stdout carries exactly one line: the JSON status line (or, under ``--probe``, the
probe result). Every human-readable message goes to stderr through logging, so
``order-qa ... | jq .verdict`` works and a tracker never has to strip log noise
out of what it parses.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from . import exits
from .integrity import probe_capabilities
from .report import DEFAULT_OUTPUT_ROOT, status_line
from .runner import RunOptions, run
from .spec import SpecError, resolve_spec
from .verify import DEFAULT_QUIESCENCE_MINUTES, list_order_objects

log = logging.getLogger("order_qa")

EPILOG = """\
exit codes:
  0  verified, quiescent, every applicable check passed
  2  bad arguments
  3  still uploading (objects newer than the quiescence window)
  4  verification failed (nothing there, or manifest mismatch)
  5  QA found defects
  6  degraded: an applicable check could not run, so there is no verdict
  7  the spec named nothing QA-able

examples:
  python -m order_qa czi-psomagen/marson-crispra/AN00028352 --assay 10x
  python -m order_qa czi-psomagen/marson-crispra/AN00028352 --assay 10x --dry-run
  python -m order_qa s3://czi-novogene/proj/NVUS2024101701-36 --assay 10x --force
  python -m order_qa czi-psomagen/marson-crispra/AN00028352 --probe
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="order-qa",
        description=(
            "QA one sequencing order without opening the notebook: verify the "
            "upload against S3, run every check that applies to the assay, and "
            "write a report. Every declared check appears in the output with a "
            "status, so a check that did not run is visible rather than absent."
        ),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "spec",
        help=(
            "Order to QA: bucket/project/order, s3://bucket/project/order, or a "
            "bare order ID resolved against the watched projects."
        ),
    )
    parser.add_argument(
        "--assay",
        default="",
        help=(
            "Raw assay: 10x, 10x_cram, 10x_viral_ORF, scale, sci_jumbo, sci_plex, "
            "seahub_sci. Required for QA -- it is what determines which files an "
            "order is expected to contain, and it varies between orders from the "
            "same vendor, so it cannot be derived from the path."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Verify the upload and stop; do not gather or run QA checks.",
    )
    parser.add_argument(
        "--probe",
        action="store_true",
        help=(
            "Probe which S3 integrity calls this identity may make against this "
            "order, print the result, and exit. Run this once per bucket to find "
            "out what the reports can claim."
        ),
    )
    parser.add_argument(
        "--quiescence-minutes",
        type=int,
        default=DEFAULT_QUIESCENCE_MINUTES,
        help=(
            "How old the newest object must be before the upload counts as "
            f"finished (default: {DEFAULT_QUIESCENCE_MINUTES})."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "QA the order even if it is still receiving objects. Recorded in the "
            "report, which will say the upload may have been in flight."
        ),
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help=(
            "CSV/TSV of expected s3:// URIs to compare the listing against. "
            "Without one, expected contents come from the assay's naming "
            "convention and the report says the manifest check did not run."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=f"Report root (default: {DEFAULT_OUTPUT_ROOT}).",
    )
    raw_group = parser.add_mutually_exclusive_group()
    raw_group.add_argument(
        "--raw",
        dest="validate_raw",
        action="store_true",
        default=None,
        help="Force raw validation on (default: on if raw/ is found).",
    )
    raw_group.add_argument(
        "--no-raw",
        dest="validate_raw",
        action="store_false",
        help="Force raw validation off.",
    )
    processed_group = parser.add_mutually_exclusive_group()
    processed_group.add_argument(
        "--processed",
        dest="validate_processed",
        action="store_true",
        default=None,
        help="Force processed validation on (default: on if processed/ is found).",
    )
    processed_group.add_argument(
        "--no-processed",
        dest="validate_processed",
        action="store_false",
        help="Force processed validation off.",
    )
    parser.add_argument(
        "--allow-truncated-stats-name",
        action="store_true",
        help=(
            "Accept *_stats.csv as an alias of *_trimmer-stats.csv. Only for "
            "collaborator uploads with the truncated naming."
        ),
    )
    parser.add_argument(
        "--head-limit",
        type=int,
        default=500,
        help=(
            "Cap on per-object HeadObject calls when ListObjectVersions is "
            "denied (default: 500). A run that hits the cap says so."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging."
    )
    return parser


def _configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


def _make_client():
    """Build the S3 client.

    No profile, no keys, no region argument: in the JupyterHub pod the credentials
    come from IRSA via the default chain, and anything explicit here would either
    override that or pick up the stale static keys sitting in ~/.aws/credentials.
    """
    import boto3

    return boto3.client("s3")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)

    try:
        spec = resolve_spec(args.spec)
    except SpecError as exc:
        log.error("%s", exc)
        return exits.RESOLUTION_FAILED

    log.info("Order %s -> s3://%s/%s", spec.key, spec.bucket, spec.prefix)

    manifest_path = Path(args.manifest) if args.manifest else None
    if manifest_path is not None and not manifest_path.exists():
        log.error("Manifest not found: %s", manifest_path)
        return exits.USAGE

    if args.probe:
        return _run_probe(spec, _make_client())

    if not args.assay.strip() and not args.dry_run:
        # Not defaulted: the assay decides which files the order is expected to
        # contain, so a wrong default does not fail loudly -- it silently checks
        # the order against the wrong expectations and reports a clean pass.
        # --dry-run needs no assay, because verification is assay-independent.
        log.error(
            "--assay is required (10x, 10x_cram, 10x_viral_ORF, scale, "
            "sci_jumbo, sci_plex, seahub_sci). Pass --dry-run to verify the "
            "upload without it."
        )
        return exits.USAGE

    output_root = (
        Path(args.output_dir).expanduser() if args.output_dir else DEFAULT_OUTPUT_ROOT
    )
    options = RunOptions(
        spec=spec,
        assay=args.assay.strip(),
        output_root=output_root,
        quiescence_minutes=args.quiescence_minutes,
        force=args.force,
        dry_run=args.dry_run,
        manifest=manifest_path,
        validate_raw=args.validate_raw,
        validate_processed=args.validate_processed,
        allow_truncated_stats_name=args.allow_truncated_stats_name,
        head_limit=args.head_limit,
    )

    try:
        report = run(options, _make_client())
    except ValueError as exc:
        # resolve_qa_run_context rejects an unknown assay this way, which is a
        # usage error rather than a fact about the order.
        log.error("%s", exc)
        return exits.USAGE
    except OSError as exc:
        # Writing the report to EFS is the one step whose failure says nothing
        # about the order. Reported as an internal error so it cannot be mistaken
        # for a verdict, and named so a wrapper can tell it apart from one.
        log.error("Could not complete or record the run: %s", exc)
        return exits.INTERNAL

    print(status_line(report))
    _log_verdict(report)
    return report.exit_code


def _log_verdict(report) -> None:
    counts = report.counts()
    log.info(
        "Verdict %s: %d passed, %d with findings, %d not applicable, %d unanswered",
        report.verdict.upper(),
        counts["pass"],
        counts["findings"],
        counts["skipped"],
        counts["not_run"] + counts["error"],
    )
    for result in report.unanswered:
        log.warning("unanswered: %s -- %s", result.name, result.reason)


def _run_probe(spec, s3_client) -> int:
    """Find out what this identity may actually do, and print it.

    Exists because the answer cannot be derived: GetObjectAttributes has never
    been exercised against these buckets, ListObjectVersions needs a permission
    separate from the listing we know we have, and the objects are encrypted with
    a key in a third account, so a call can fail on the key policy rather than the
    bucket policy. Probing beats assuming in both directions -- claiming a
    checksum check that was denied is a false assurance, and skipping one that was
    available is a wasted capability.
    """
    listing = list_order_objects(s3_client, spec.bucket, spec.prefix)
    if not listing.complete:
        log.error(
            "Could not list s3://%s/%s: %s", spec.bucket, spec.prefix, listing.error
        )
        return exits.DEGRADED
    files = listing.files
    if not files:
        log.error(
            "No objects under s3://%s/%s to probe against.", spec.bucket, spec.prefix
        )
        return exits.VERIFICATION_FAILED

    multipart = next((o for o in files if o.multipart_parts), None)
    probe_key = (multipart or files[0]).key
    log.info(
        "Probing against %s (%s)",
        probe_key,
        f"multipart, {multipart.multipart_parts} parts" if multipart else "single-part",
    )
    report = probe_capabilities(
        s3_client, spec.bucket, probe_key=probe_key, prefix=spec.prefix
    )
    payload = report.as_dict()
    payload["listing"] = {
        "file_count": len(files),
        "multipart_files": sum(1 for o in files if o.multipart_parts),
        "single_part_files": sum(1 for o in files if o.etag_is_md5),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    for capability in report.capabilities:
        log.info("%-26s %s", capability.name, capability.status)
    return exits.OK


if __name__ == "__main__":
    raise SystemExit(main())
