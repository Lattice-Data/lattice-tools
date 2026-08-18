"""Writing the report to EFS, and the status line to stdout.

Two audiences, deliberately separated. Everything a person needs to act is in the
report on disk; stdout carries exactly one line of JSON so a tracker can consume a
run without parsing prose. Human-readable progress goes to stderr via logging, so
``order-qa ... | jq`` works and nothing has to be scraped.

Run directories are timestamped and never rewritten, because a re-run after a
resequence must not destroy the report that justified the first verdict.

The report separates what was *verified* from what was *assumed*. A denied
capability, a bounded sweep, or a check that could not run each appear as such --
an order with nine passes and three unanswered questions is not the same thing as
an order with twelve passes, and the report must not let those look alike.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import exits
from .registry import CheckResult, Status

log = logging.getLogger(__name__)

DEFAULT_OUTPUT_ROOT = Path.home() / "qa-reports"
TIMESTAMP_FORMAT = "%Y%m%dT%H%M%SZ"


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def run_stamp(moment: datetime) -> str:
    return moment.strftime(TIMESTAMP_FORMAT)


@dataclass
class RunReport:
    """Everything one run produced."""

    order_key: str
    bucket: str
    project: str
    order: str
    order_shape: str
    assay: str
    started_at: datetime
    scope: str = ""
    exit_code: int = exits.OK
    listing_summary: dict[str, Any] = field(default_factory=dict)
    quiescence: dict[str, Any] = field(default_factory=dict)
    manifest: dict[str, Any] = field(default_factory=dict)
    capabilities: dict[str, Any] = field(default_factory=dict)
    versions: dict[str, Any] = field(default_factory=dict)
    rerun: dict[str, Any] = field(default_factory=dict)
    gathering_errors: list[str] = field(default_factory=list)
    gathering_warnings: list[str] = field(default_factory=list)
    checks: list[CheckResult] = field(default_factory=list)
    report_dir: Path | None = None
    # Why the check table is empty, when it is. An empty table with no
    # explanation reads as "there was nothing to check", which is never true.
    qa_skipped_reason: str = ""

    @property
    def verdict(self) -> str:
        return exits.name_for(self.exit_code)

    def counts(self) -> dict[str, int]:
        tally = {status.value: 0 for status in Status}
        for result in self.checks:
            tally[result.status.value] += 1
        return tally

    @property
    def total_findings(self) -> int:
        return sum(len(r.findings) for r in self.checks)

    @property
    def unanswered(self) -> list[CheckResult]:
        """Applicable checks with no answer. The reason the verdict can be denied."""
        return [r for r in self.checks if r.status.denies_verdict]

    def as_dict(self) -> dict[str, Any]:
        return {
            "order": self.order_key,
            "bucket": self.bucket,
            "project": self.project,
            "order_id": self.order,
            "order_shape": self.order_shape,
            "assay": self.assay,
            "scope": self.scope,
            "started_at": self.started_at.isoformat(),
            "verdict": self.verdict,
            "exit_code": self.exit_code,
            "counts": self.counts(),
            "total_findings": self.total_findings,
            "listing": self.listing_summary,
            "quiescence": self.quiescence,
            "manifest": self.manifest,
            "capabilities": self.capabilities,
            "versions": self.versions,
            "rerun": self.rerun,
            "gathering_errors": self.gathering_errors,
            "gathering_warnings": self.gathering_warnings,
            "qa_skipped_reason": self.qa_skipped_reason,
            "checks": [r.as_dict() for r in self.checks],
            "report_dir": str(self.report_dir) if self.report_dir else "",
        }


def order_root(output_root: Path, bucket: str, project: str, order: str) -> Path:
    """``{root}/{bucket}/{project}/{order}`` -- the per-order directory."""
    return output_root / bucket / project / order


def make_run_dir(root: Path, moment: datetime) -> Path:
    """Create a fresh timestamped run directory under the order root.

    Two runs inside the same second get ``-2``, ``-3`` suffixes rather than
    sharing a directory. Same-second collisions are unlikely by hand but entirely
    normal for a retry loop or a future CronJob, and silently merging two runs'
    outputs would leave a report describing neither.
    """
    stamp = run_stamp(moment)
    candidate = root / stamp
    suffix = 2
    while candidate.exists():
        candidate = root / f"{stamp}-{suffix}"
        suffix += 1
    candidate.mkdir(parents=True)
    return candidate


def status_line(report: RunReport) -> str:
    """The single line of JSON a downstream tracker parses."""
    counts = report.counts()
    return json.dumps(
        {
            "order": report.order_key,
            "assay": report.assay,
            "verdict": report.verdict,
            "exit_code": report.exit_code,
            "checks_passed": counts[Status.PASS.value],
            "checks_with_findings": counts[Status.FINDINGS.value],
            "checks_skipped": counts[Status.SKIPPED.value],
            "checks_unanswered": counts[Status.NOT_RUN.value]
            + counts[Status.ERROR.value],
            "findings": report.total_findings,
            "files": report.listing_summary.get("file_count", 0),
            "bytes": report.listing_summary.get("total_bytes", 0),
            "quiet": report.quiescence.get("quiet"),
            "changed_since_last_run": report.rerun.get("any_change"),
            "report_dir": str(report.report_dir) if report.report_dir else "",
            "timestamp": report.started_at.isoformat(),
        },
        sort_keys=True,
    )


_STATUS_MARK = {
    Status.PASS: "PASS",
    Status.FINDINGS: "FINDINGS",
    Status.SKIPPED: "skipped",
    Status.NOT_RUN: "NOT RUN",
    Status.ERROR: "ERROR",
}


def render_markdown(report: RunReport) -> str:
    """The report a person reads."""
    lines: list[str] = []
    add = lines.append

    add(f"# QA report: {report.order_key}")
    add("")
    add(f"**Verdict: {report.verdict.upper()}** (exit {report.exit_code})")
    add("")
    add(f"- Assay: `{report.assay}`")
    add(f"- Scope: {report.scope}")
    add(f"- Order ID shape: {report.order_shape}")
    add(f"- Started: {report.started_at.isoformat()}")
    add("")

    if report.unanswered:
        add(
            "> **This order has no clean verdict available.** "
            f"{len(report.unanswered)} applicable check(s) produced no answer; "
            "they are listed as NOT RUN or ERROR below. A check that could not "
            "run is not a check that passed."
        )
        add("")

    add("## Upload verification")
    add("")
    listing = report.listing_summary
    add(
        f"- Files: {listing.get('file_count', 0)} "
        f"({listing.get('folder_markers', 0)} folder marker(s) excluded)"
    )
    add(f"- Bytes: {listing.get('total_bytes', 0):,}")
    if not listing.get("complete", True):
        add(
            f"- **Listing incomplete**: {listing.get('error', 'unknown error')} "
            "— every count above is a floor, not a total."
        )
    add(f"- Quiescence: {report.quiescence.get('summary', 'not assessed')}")
    manifest = report.manifest
    if manifest.get("checked"):
        add(
            f"- Manifest: {manifest.get('present', 0)} present of "
            f"{manifest.get('expected', 0)} expected, "
            f"{len(manifest.get('missing', []))} missing, "
            f"{len(manifest.get('extra', []))} extra"
        )
    else:
        add(f"- Manifest: not checked — {manifest.get('reason', 'no manifest')}")
    add(f"- Re-run: {report.rerun.get('summary', 'not compared')}")
    add("")

    add("## What could be checked")
    add("")
    add(
        "Capabilities probed against a real object in this order, so the limits "
        "below are measured rather than assumed."
    )
    add("")
    caps = report.capabilities.get("capabilities") or {}
    if caps:
        add("| capability | status | detail |")
        add("| --- | --- | --- |")
        for name, info in sorted(caps.items()):
            add(f"| `{name}` | {info.get('status')} | {info.get('detail', '')} |")
    else:
        add(f"_No capabilities probed: {report.capabilities.get('reason', 'unknown')}_")
    add("")
    versions = report.versions
    if versions.get("checked"):
        coverage = f"{versions.get('covered', 0)}/{versions.get('total', 0)}"
        add(
            f"- Version identities: {coverage} key(s) via "
            f"`{versions.get('method', '?')}`"
            + (
                " — **partial**, so an unrecorded key may have been re-uploaded "
                "without being noticed."
                if versions.get("partial")
                else ""
            )
        )
    else:
        add(
            f"- Version identities: not collected — {versions.get('reason', 'unknown')}"
        )
    add("")

    add("## Checks")
    add("")
    if not report.checks:
        add(
            f"**No checks were run.** {report.qa_skipped_reason or 'Reason unrecorded.'}"
        )
        add("")
        return "\n".join(lines) + "\n"
    counts = report.counts()
    add(
        f"{len(report.checks)} declared · {counts[Status.PASS.value]} passed · "
        f"{counts[Status.FINDINGS.value]} with findings · "
        f"{counts[Status.SKIPPED.value]} not applicable · "
        f"{counts[Status.NOT_RUN.value]} not run · "
        f"{counts[Status.ERROR.value]} errored"
    )
    add("")
    add("| check | status | detail |")
    add("| --- | --- | --- |")
    for result in report.checks:
        detail = result.summary or result.reason
        add(f"| `{result.name}` | {_STATUS_MARK[result.status]} | {detail} |")
    add("")

    findings_present = [r for r in report.checks if r.findings]
    if findings_present:
        add("## Findings")
        add("")
        for result in findings_present:
            add(f"### {result.name} — {result.title}")
            add("")
            for finding in result.findings:
                add(f"- {finding}")
            add("")

    noted = [r for r in report.checks if r.notes]
    if noted:
        add("## Notes")
        add("")
        for result in noted:
            add(f"### {result.name}")
            add("")
            for note in result.notes:
                add(f"- {note}")
            add("")

    if report.gathering_errors or report.gathering_warnings:
        add("## Gathering")
        add("")
        for warning in report.gathering_warnings:
            add(f"- WARNING {warning}")
        for error in report.gathering_errors:
            add(f"- {error}")
        add("")

    return "\n".join(lines) + "\n"


def write_report(report: RunReport, run_dir: Path) -> dict[str, Path]:
    """Write the markdown report, the JSON summary, and any per-check tables."""
    written: dict[str, Path] = {}

    summary_path = run_dir / "summary.json"
    summary_path.write_text(
        json.dumps(report.as_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    written["summary"] = summary_path

    markdown_path = run_dir / "report.md"
    markdown_path.write_text(render_markdown(report), encoding="utf-8")
    written["report"] = markdown_path

    # The notebook's per-check CSVs, kept because people already read them:
    # {order}_raw_missing.csv, {order}_process_alerts.csv and friends.
    for result in report.checks:
        for table_name, rows in result.tables.items():
            if not rows:
                continue
            path = run_dir / f"{table_name}.csv"
            fieldnames = sorted({key for row in rows for key in row})
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
            written[table_name] = path

    # One flat file of every finding, matching the notebook's {order}_errors.txt,
    # which is what people grep.
    findings = [
        f"[{result.name}] {finding}"
        for result in report.checks
        for finding in result.findings
    ]
    errors_path = run_dir / "errors.txt"
    errors_path.write_text(
        "\n".join(findings) + ("\n" if findings else ""), encoding="utf-8"
    )
    written["errors"] = errors_path

    return written
