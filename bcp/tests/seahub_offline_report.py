"""Re-run a SeaHub QA report from a saved S3 listing, with no network access.

Why this exists
---------------
The numbers quoted for a SeaHub upload -- SOP rows, per-well verdicts, rename
statuses -- come from running the notebook against a bucket. Anyone checking
those numbers afterwards would have to have the same S3 access on the same day.
This drives the same production code from a *listing file* instead, so a claim
about an upload can be re-derived by whoever has the listing, and a change can be
checked against real shapes before it goes near a bucket.

It is a reporting tool, not a test. It lives in ``tests/`` because that is where
``MockS3Client`` lives and duplicating an S3 stub to avoid the import would be
worse. ``pytest`` does not collect it (``python_files = test_*.py``).

Sanitization
------------
No bucket, lab, project, order or experiment identifier appears here. Listings
are passed as paths, and every identity is read out of the listing at run time --
including the ExperimentID, which ``resolve_qa_run_context`` recovers from the
keys. Keep it that way: this file is inside the tree
``test_sanitized_identifiers.py`` guards.

Input
-----
Any CSV/TSV with a column of ``s3://bucket/key`` URIs -- the column is found by
looking for one, so an ``aws s3 ls`` dump or a console export both work. A
``Size_Bytes`` column is used when present, since the trimmed-vs-vendor size
check reads it.

Usage
-----
    python -m tests.seahub_offline_report LISTING [LISTING ...]
    python -m tests.seahub_offline_report LISTING --vendor VENDOR [VENDOR ...]
    python -m tests.seahub_offline_report LISTING --mode both --json > after.json

``--mode both`` runs the same keys through the S3 walk and through manifest
mode, which must agree; ``--json`` emits a diffable snapshot, so the effect of a
change is ``diff before.json after.json``.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from qa_gather import gather_qa_data
from qa_mods import QARunContext, is_s3_folder_marker, resolve_qa_run_context
from qa_seahub_rename import build_rename_mapping, roll_up_wells, rollup_summary
from qa_seahub_sop import sop_violation_summary, validate_seahub_stems
from qa_seahub_source import (
    discover_untrimmed_sources,
    index_trimmed_upload,
    seahub_wafer_seeds,
)

from tests.test_qa_gather import MockS3Client

# Enough of a fail CSV and a sidecar to exercise the parsers the gatherer runs on
# every well. The listing records that these objects exist but not what is in
# them, and the report is about names and inventory, not read counts.
_FAIL_CSV = (
    "format,reason,failed read count,total read count\n"
    "JumboSciGEX,barcode,100,10000\n"
    "JumboSciHash,barcode,50,8000\n"
)
_CRAM_METADATA = '{"read_count": 10000}'


def load_listing(path: str | Path) -> tuple[str, list[str], dict[str, int]]:
    """Read ``(bucket, keys, sizes)`` out of a listing file.

    The S3 column is located by content rather than by name or position, because
    these files come from several tools and agree only on holding URIs somewhere.
    """
    rows = _read_rows(Path(path))
    bucket = ""
    keys: list[str] = []
    sizes: dict[str, int] = {}
    for row in rows:
        uri = next((cell for cell in row if str(cell).startswith("s3://")), "")
        if not uri:
            continue
        found, _, key = uri[len("s3://") :].partition("/")
        if not key:
            continue
        if is_s3_folder_marker(key):
            continue
        bucket = bucket or found
        keys.append(key)
        size = _first_int(row)
        # Only when the row actually carries one. A key absent from the mapping
        # means no size was collected; a key present with 0 means the object is
        # genuinely empty, and a 0-byte CRAM is a DATA_GAP. Defaulting to 0 here
        # conflated the two, so any listing without a size column -- a bare
        # manifest of URIs, for instance -- reported every well in the upload as a
        # data gap. The tool exists to check numbers, so it must not invent them.
        if size is not None:
            sizes[key] = size
    if not keys:
        raise SystemExit(f"{path}: no s3:// URIs found in any column")
    return bucket, keys, sizes


def _read_rows(path: Path) -> list[list[str]]:
    delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    with path.open(newline="") as fh:
        sample = fh.readline()
        fh.seek(0)
        if "\t" in sample and delimiter == ",":
            delimiter = "\t"
        return [row for row in csv.reader(fh, delimiter=delimiter) if row]


def _first_int(row: list[str]) -> int | None:
    """The first plain integer in the row -- the size column -- or None.

    None rather than 0, because the two mean different things downstream and only
    one of them is a finding.
    """
    for cell in row:
        text = str(cell).strip()
        if text.isdigit():
            return int(text)
    return None


def _synthetic_contents(keys: list[str]) -> dict[str, str]:
    """Bodies for the objects the gatherer downloads, keyed by S3 key."""
    contents: dict[str, str] = {}
    for key in keys:
        if key.endswith((".trim_fail.csv", "_fail.csv", "failure_codes.csv")):
            contents[key] = _FAIL_CSV
        elif key.endswith("-metadata.json"):
            contents[key] = _CRAM_METADATA
    return contents


def _manifest_context(
    listing_keys: list[str], bucket: str, label: str, scratch: Path
) -> QARunContext:
    """Resolve a manifest context, which also recovers the ExperimentID."""
    manifest = scratch / f"{label}.manifest.tsv"
    manifest.write_text("".join(f"s3://{bucket}/{key}\n" for key in listing_keys))
    return resolve_qa_run_context(
        data_source="manifest",
        raw_assay="seahub_sci",
        manifest_path=str(manifest),
        manifest_delimiter="\t",
        manifest_s3_column=0,
        manifest_has_header=False,
        run_label=label,
    )


def _s3_context(manifest_ctx: QARunContext, proj: str) -> QARunContext:
    """The s3-mode twin of a manifest context, over the same keys."""
    return QARunContext(
        data_source="s3",
        raw_assay="seahub_sci",
        bucket=manifest_ctx.bucket,
        provider=manifest_ctx.provider,
        proj=proj,
        order=manifest_ctx.order,
        output_label=manifest_ctx.output_label,
        listing_prefix=f"{proj}/{manifest_ctx.order}/",
    )


def report(
    listing: str | Path,
    vendors: list[str] | None = None,
    mode: str = "s3",
    scratch: Path | None = None,
) -> dict[str, Any]:
    """Run one listing through the gatherer and the notebook's SeaHub cells."""
    if scratch is None:
        scratch = Path(tempfile.mkdtemp(prefix="seahub-offline-"))
    bucket, keys, sizes = load_listing(listing)
    proj = keys[0].split("/")[0]
    label = Path(listing).stem

    manifest_ctx = _manifest_context(keys, bucket, label, scratch)
    if mode == "manifest":
        ctx = manifest_ctx
        client = MockS3Client(keys=[], file_contents=_synthetic_contents(keys))
    else:
        ctx = _s3_context(manifest_ctx, proj)
        client = MockS3Client(
            keys=keys, sizes=sizes, file_contents=_synthetic_contents(keys)
        )

    data = gather_qa_data(ctx, client)

    untrimmed_index: dict = {}
    assay_by_identity: dict = {}
    findings: list = []
    search = None
    if vendors:
        vendor_keys: list[str] = []
        vendor_bucket = ""
        for path in vendors:
            found, more, _ = load_listing(path)
            vendor_bucket = vendor_bucket or found
            vendor_keys.extend(more)
        # Search roots, not order prefixes: the project alone, so the run
        # exercises the same discovery the notebook does rather than being handed
        # the answer. Deriving {project}/{order} from the vendor keys, as this
        # used to, told the indexer where to look and so could not have caught a
        # descent that failed to find it.
        roots = sorted({k.split("/")[0] for k in vendor_keys})
        seeds = seahub_wafer_seeds(
            trimmed_keys=data.all_raw_files,
            trimmed_index=index_trimmed_upload(data.all_raw_files),
            discovered_wafers=data.discovered_wafers,
        )
        sources = discover_untrimmed_sources(
            MockS3Client(
                buckets={vendor_bucket: vendor_keys},
                file_contents=_synthetic_contents(vendor_keys),
            ),
            [f"s3://{vendor_bucket}/{root}" for root in roots],
            seeds,
        )
        untrimmed_index = sources.index
        findings = list(sources.findings)
        search = sources.search
        assay_by_identity = {k: e.assay for k, e in untrimmed_index.items() if e.assay}

    violations = validate_seahub_stems(
        ctx.bucket, data.all_raw_files, assay_by_identity=assay_by_identity
    )
    mapping = build_rename_mapping(ctx.bucket, data.all_raw_files, untrimmed_index)
    rollup = roll_up_wells(
        ctx.bucket,
        data.all_raw_files,
        untrimmed_index,
        # None, not an all-zero mapping, when the listing carried no size column.
        # A key present with 0 means the CRAM is genuinely empty, which is a
        # DATA_GAP; a key absent means nobody measured it. MockS3Client fills a
        # missing Size with 0 in its listing output -- reasonable for a stub of an
        # API that always returns one -- so handing it no sizes produced a
        # complete set of zeros and reported every well of a clean upload as a
        # data gap. The distinction is only knowable here, where whether the
        # listing had the column is still in scope.
        sizes=data.raw_file_sizes if sizes else None,
    )

    return {
        "experiment_id": ctx.order,
        "listing_objects": len(keys),
        "all_raw_files": len(data.all_raw_files),
        "raw_file_sizes": len(data.raw_file_sizes),
        # Populated by the folder walk, so s3 mode only -- see the parity note in
        # SEAHUB_QA.md. Not part of what the two modes are asserted to agree on.
        "discovered_wafers": len(data.discovered_wafers),
        "sop_rows": len(violations),
        "sop_rules": dict(sop_violation_summary(violations)),
        "wells": len(rollup.rows),
        "well_verdicts": rollup_summary(rollup.rows),
        "unaccounted": rollup.unaccounted,
        "rename_total": mapping.total_objects,
        "rename_compliant": mapping.compliant_objects,
        "rename_counts": dict(mapping.counts),
        "vendor_indexed": len(untrimmed_index),
        "source_findings": dict(Counter(f["category"] for f in findings)),
        # Discovery accounting, so a differential shows whether the search saw
        # everything and not only what it did with what it saw.
        "wafers_sought": len(search.seeds.wafers) if search else 0,
        "wafers_located": len(search.located) if search else 0,
        "wafers_not_located": sorted(search.not_located) if search else [],
        "search_listings": search.listings if search else 0,
        "search_complete": bool(search.complete) if search else None,
        "warnings": len(data.gathering_warnings),
        "errors": len(data.gathering_errors),
        "warning_kinds": dict(
            Counter(w.split(":")[0][:44] for w in data.gathering_warnings)
        ),
        "error_kinds": dict(
            Counter(e.split(":")[0][:44] for e in data.gathering_errors)
        ),
    }


# The keys whose agreement between the two modes is a property of the gatherer.
# discovered_wafers is deliberately absent: only the s3 walk enumerates folders.
PARITY_FIELDS = (
    "all_raw_files",
    "sop_rows",
    "sop_rules",
    "wells",
    "well_verdicts",
    "rename_counts",
    "unaccounted",
    "errors",
)


def parity_diff(from_s3: dict, from_manifest: dict) -> dict[str, tuple]:
    """Fields on which the two modes disagree; empty means they agree."""
    return {
        field: (from_s3[field], from_manifest[field])
        for field in PARITY_FIELDS
        if from_s3[field] != from_manifest[field]
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Re-run a SeaHub QA report from a saved S3 listing.",
    )
    parser.add_argument("listings", nargs="+", help="S3 listing files (CSV/TSV)")
    parser.add_argument(
        "--vendor",
        nargs="+",
        default=[],
        metavar="LISTING",
        help="untrimmed vendor listing(s) to reconcile against",
    )
    parser.add_argument("--mode", default="s3", choices=("s3", "manifest", "both"))
    parser.add_argument("--json", action="store_true", help="diffable snapshot")
    parser.add_argument(
        "--scratch",
        default="",
        help="where to write intermediate manifests (default: a temp directory)",
    )
    args = parser.parse_args(argv)

    modes = ("s3", "manifest") if args.mode == "both" else (args.mode,)
    scratch = Path(args.scratch) if args.scratch else None
    results: dict[str, dict] = {}

    # gather_qa_data prints download progress; keep it off a --json snapshot.
    with contextlib.redirect_stdout(sys.stderr):
        for listing in args.listings:
            for mode in modes:
                label = f"{Path(listing).stem}/{mode}"
                results[label] = report(
                    listing, vendors=args.vendor, mode=mode, scratch=scratch
                )

    if args.json:
        print(json.dumps(results, indent=2, sort_keys=True))
        return 0

    for label, result in results.items():
        print(f"=== {label}")
        for key, value in result.items():
            print(f"    {key:20} {value}")
        print()

    if args.mode == "both":
        for listing in args.listings:
            stem = Path(listing).stem
            disagreement = parity_diff(
                results[f"{stem}/s3"], results[f"{stem}/manifest"]
            )
            verdict = disagreement or "s3 and manifest agree"
            print(f"parity {stem}: {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
