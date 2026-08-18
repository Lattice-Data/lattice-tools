"""Order of operations for one QA run.

Verification comes before QA, and quiescence before both. Running QA against an
order that is still uploading burns a full gather and reports absences as defects,
so the cheap listing-based questions are asked first and a still-arriving upload
stops there.

Nothing in here decides what a check means; it decides what gets asked, in what
order, and what the combination adds up to.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from qa_gather import gather_qa_data
from qa_mods import resolve_qa_run_context

from . import exits
from .integrity import collect_version_state, probe_capabilities
from .ledger import Ledger, diff_against_ledger, load_ledger, save_ledger
from .registry import CheckContext, Scope, Status, run_registry
from .report import (
    RunReport,
    make_run_dir,
    order_root,
    utc_now,
    write_report,
)
from .spec import OrderSpec
from .verify import (
    assess_quiescence,
    compare_to_manifest,
    list_order_objects,
)

log = logging.getLogger(__name__)

# Gathering warnings that mean a check did not see everything it should have.
# ``qa_gather`` emits these alongside genuinely advisory warnings; the notebook
# writes all of them to the error log precisely because these four say coverage
# was incomplete. An incomplete check is not a passed check, so these deny the
# verdict rather than sitting in a notes section.
COVERAGE_WARNING_MARKERS = (
    "LISTING TRUNCATED",
    "METADATA UNREADABLE",
    "TRIM FAIL UNREADABLE",
    "EXPERIMENT UNKNOWN",
)


@dataclass
class RunOptions:
    """Everything the CLI decided, resolved."""

    spec: OrderSpec
    assay: str
    output_root: Path
    quiescence_minutes: int
    force: bool = False
    dry_run: bool = False
    manifest: Path | None = None
    validate_raw: bool | None = None
    validate_processed: bool | None = None
    allow_truncated_stats_name: bool = False
    head_limit: int = 500
    now: datetime | None = None


def run(options: RunOptions, s3_client: Any) -> RunReport:
    """Verify, then QA, then report. Returns the report with its exit code set."""
    started = options.now or utc_now()
    spec = options.spec
    report = RunReport(
        order_key=spec.key,
        bucket=spec.bucket,
        project=spec.project,
        order=spec.order,
        order_shape=spec.order_shape,
        assay=options.assay,
        started_at=started,
    )

    log.info("Listing s3://%s/%s", spec.bucket, spec.prefix)
    listing = list_order_objects(s3_client, spec.bucket, spec.prefix)
    files = listing.files
    report.listing_summary = {
        "file_count": len(files),
        "object_count": len(listing.objects),
        "folder_markers": len(listing.objects) - len(files),
        "total_bytes": listing.total_bytes,
        "complete": listing.complete,
        "error": listing.error or "",
        "multipart_files": sum(1 for o in files if o.multipart_parts),
        "etag_is_md5_files": sum(1 for o in files if o.etag_is_md5),
    }
    log.info(
        "%d file(s), %d byte(s); %d multipart (ETag is not MD5 for those)",
        len(files),
        listing.total_bytes,
        report.listing_summary["multipart_files"],
    )

    quiescence = assess_quiescence(
        listing,
        threshold_minutes=options.quiescence_minutes,
        force=options.force,
        now=started,
    )
    report.quiescence = {
        "checked": quiescence.checked,
        "quiet": quiescence.quiet,
        "forced": quiescence.forced,
        "newest": quiescence.newest.isoformat() if quiescence.newest else "",
        "age_seconds": quiescence.age_seconds,
        "threshold_seconds": quiescence.threshold_seconds,
        "summary": quiescence.summary,
    }
    log.info("Quiescence: %s", quiescence.summary)

    report.manifest = _manifest_dict(compare_to_manifest(listing, options.manifest))

    probe_key = _pick_probe_key(files)
    capabilities = probe_capabilities(
        s3_client, spec.bucket, probe_key=probe_key, prefix=spec.prefix
    )
    report.capabilities = capabilities.as_dict()
    for capability in capabilities.capabilities:
        log.info("capability %-26s %s", capability.name, capability.status)

    version_state = collect_version_state(
        s3_client,
        spec.bucket,
        spec.prefix,
        [o.key for o in files],
        capabilities=capabilities,
        head_limit=options.head_limit,
    )
    report.versions = {
        "checked": version_state.checked,
        "method": version_state.method,
        "covered": version_state.covered,
        "total": version_state.total,
        "partial": version_state.partial,
        "delete_markers": version_state.delete_markers,
        "reason": version_state.reason,
    }
    if version_state.partial:
        log.warning(
            "Version identities cover %d of %d key(s); a re-upload of an "
            "unrecorded key would go unnoticed.",
            version_state.covered,
            version_state.total,
        )

    root = order_root(options.output_root, spec.bucket, spec.project, spec.order)
    previous = load_ledger(root)
    diff = diff_against_ledger(
        previous,
        current_keys=listing.keys,
        current_versions=version_state.version_by_key,
    )
    report.rerun = {
        "compared": diff.compared,
        "reason": diff.reason,
        "previous_run_at": diff.previous_run_at,
        "previous_verdict": diff.previous_verdict,
        "new_keys": diff.new_keys,
        "removed_keys": diff.removed_keys,
        "reuploaded_keys": diff.reuploaded_keys,
        "unchanged": diff.unchanged,
        "uncomparable_keys": diff.uncomparable_keys,
        "any_change": diff.any_change,
        "summary": diff.summary,
    }
    log.info("Re-run: %s", diff.summary)

    # Only a complete listing of a settled order is a baseline worth keeping.
    # Recording a partial or in-flight observation would overwrite the last good
    # one, and the next run would then diff against it: every object the failed
    # run never saw comes back as "new", and -- worse -- the version identities
    # that would have caught a silent re-upload are gone.
    trustworthy_baseline = listing.complete and quiescence.quiet

    stop_code = _verification_verdict(report, quiescence, listing)
    if stop_code is not None:
        report.exit_code = stop_code
        report.qa_skipped_reason = (
            f"Verification stopped the run before QA: {exits.name_for(stop_code)}."
        )
        _finish(
            report, options, root, version_state, update_ledger=trustworthy_baseline
        )
        return report

    if options.dry_run:
        log.info("--dry-run: verification only, QA checks not run.")
        report.exit_code = exits.OK
        report.qa_skipped_reason = (
            "--dry-run was passed, so the upload was verified but no QA check "
            "was run. This is not a statement about the data."
        )
        _finish(
            report, options, root, version_state, update_ledger=trustworthy_baseline
        )
        return report

    _run_qa(report, options, s3_client)
    report.exit_code = _qa_verdict(report)
    _finish(report, options, root, version_state, update_ledger=trustworthy_baseline)
    return report


def _manifest_dict(comparison) -> dict[str, Any]:
    return {
        "checked": comparison.checked,
        "reason": comparison.reason,
        "expected": comparison.expected,
        "present": comparison.present,
        "missing": comparison.missing,
        "extra": comparison.extra,
        "size_mismatches": [
            {"key": key, "expected": expected, "actual": actual}
            for key, expected, actual in comparison.size_mismatches
        ],
        "ok": comparison.ok,
    }


def _pick_probe_key(files) -> str:
    """Choose an object to probe capabilities against.

    Prefers a multipart object, because the part-level calls are the ones whose
    availability is in question and a single-part object cannot exercise them.
    """
    multipart = next((o for o in files if o.multipart_parts), None)
    if multipart is not None:
        return multipart.key
    return files[0].key if files else ""


def _verification_verdict(report: RunReport, quiescence, listing) -> int | None:
    """Verification outcomes that stop the run before QA, or None to continue."""
    if not listing.complete:
        log.error(
            "The listing did not complete (%s), so neither completeness nor "
            "quiescence means anything for this order.",
            listing.error,
        )
        return exits.DEGRADED
    if not listing.files:
        log.error("No objects under s3://%s/%s.", listing.bucket, listing.prefix)
        return exits.VERIFICATION_FAILED
    if not quiescence.quiet:
        log.error(
            "Objects were still arriving %s ago; the upload may be in flight. "
            "Re-run later, or pass --force to QA it anyway.",
            quiescence.age_human,
        )
        return exits.STILL_UPLOADING
    manifest = report.manifest
    if manifest["checked"] and (manifest["missing"] or manifest["size_mismatches"]):
        log.error(
            "Manifest comparison failed: %d missing, %d size mismatch(es).",
            len(manifest["missing"]),
            len(manifest["size_mismatches"]),
        )
        return exits.VERIFICATION_FAILED
    return None


def _run_qa(report: RunReport, options: RunOptions, s3_client: Any) -> None:
    """Gather with the notebook's own gatherer, then run every declared check."""
    ctx = resolve_qa_run_context(
        data_source="s3",
        raw_assay=options.assay,
        s3_path=options.spec.s3_path,
        allow_truncated_stats_name=options.allow_truncated_stats_name,
    )
    log.info("Gathering QA data for %s/%s", ctx.bucket, ctx.listing_prefix)
    data = gather_qa_data(ctx, s3_client)
    report.gathering_errors = list(data.gathering_errors)
    report.gathering_warnings = list(data.gathering_warnings)
    for warning in report.gathering_warnings:
        log.warning("%s", warning)
    for error in report.gathering_errors:
        log.error("%s", error)

    scope = Scope.resolve(
        raw_assay=ctx.raw_assay,
        validate_raw=(
            data.has_raw if options.validate_raw is None else options.validate_raw
        ),
        validate_processed=(
            data.has_processed
            if options.validate_processed is None
            else options.validate_processed
        ),
    )
    report.scope = scope.describe()
    log.info("Scope: %s", scope.describe())

    check_ctx = CheckContext(
        data=data,
        scope=scope,
        bucket=ctx.bucket,
        order_label=ctx.output_label,
        s3_client=s3_client,
        output_dir="",
        allow_truncated_stats_name=ctx.allow_truncated_stats_name,
    )
    report.checks = run_registry(check_ctx)


def _qa_verdict(report: RunReport) -> int:
    """Combine check outcomes into one exit code.

    Precedence follows what the reader should do first: an order with unanswered
    questions needs those answered before its findings are worth triaging, because
    the findings may not be the whole story.
    """
    coverage_warnings = [
        warning
        for warning in report.gathering_warnings
        if any(marker in warning for marker in COVERAGE_WARNING_MARKERS)
    ]
    if report.unanswered or report.gathering_errors or coverage_warnings:
        for warning in coverage_warnings:
            log.error("Incomplete coverage: %s", warning)
        return exits.DEGRADED

    # An order where every check was inapplicable was not validated at all, and
    # must not report as validated. This happens for real: the gather finds
    # neither raw/ nor processed/ -- a vendor uploading to an unexpected layout,
    # or the wrong assay passed -- and every check skips legitimately, which
    # summed to a clean pass on an order nothing had looked at.
    answered = [r for r in report.checks if r.status in (Status.PASS, Status.FINDINGS)]
    if not answered:
        log.error(
            "No check applied to this order: %d of %d were skipped. Nothing was "
            "validated, so there is no verdict. Check --assay, and whether the "
            "upload has the raw/ and processed/ layout QA expects.",
            sum(1 for r in report.checks if r.status is Status.SKIPPED),
            len(report.checks),
        )
        return exits.DEGRADED

    if any(r.status is Status.FINDINGS for r in report.checks):
        return exits.QA_FINDINGS
    return exits.OK


def _finish(
    report: RunReport,
    options: RunOptions,
    root: Path,
    version_state,
    *,
    update_ledger: bool,
) -> None:
    """Write the run directory, and replace the ledger only if this run earned it.

    The report is always written -- a failed run is exactly the one someone needs
    to read. The ledger is not, because it is the baseline the *next* run compares
    against; see the caller for why a partial observation must not become one.
    """
    run_dir = make_run_dir(root, report.started_at)
    report.report_dir = run_dir
    written = write_report(report, run_dir)
    log.info("Report: %s", written["report"])

    if not update_ledger:
        log.warning(
            "Not updating the re-run baseline: this run did not observe the whole "
            "order (listing complete=%s, quiescent=%s). The previous baseline is "
            "kept so the next run still has something to compare against.",
            report.listing_summary.get("complete"),
            report.quiescence.get("quiet"),
        )
        return

    save_ledger(
        root,
        Ledger(
            order_key=report.order_key,
            run_at=report.started_at.isoformat(),
            verdict=report.verdict,
            file_count=report.listing_summary.get("file_count", 0),
            total_bytes=report.listing_summary.get("total_bytes", 0),
            newest_object=report.quiescence.get("newest", ""),
            version_by_key=version_state.version_by_key,
            version_method=version_state.method,
            version_coverage=f"{version_state.covered}/{version_state.total}",
            report_dir=str(run_dir),
        ),
    )
