"""The QA checks, declared as data.

``qa.ipynb`` decides what applies inside each cell: ``if validate_raw and
raw_assay == 'scale'``, and so on, twenty-odd times. Every branch prints why it
skipped, which is correct and also invisible -- a skipped cell scrolls past
exactly like a passing one, and there is no point at which the notebook can tell
you that eleven of its checks never ran.

Here each check declares the conditions it applies under, so the run can be
enumerated before it starts and every declared check appears in the output with a
status. "Skipped" is a reported outcome. The distinction that matters is between
SKIPPED, which means the check does not apply to this order and implies nothing,
and NOT_RUN, which means it does apply and could not be performed -- the second
one denies the order a verdict.

What a check *means* is not defined here. Every one delegates to the same
``qa_checks`` function the notebook calls.
"""

from __future__ import annotations

import logging
import os
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from qa_checks import (
    MIN_PCT_Q30_THRESHOLD,
    check_expected_raw_files,
    check_extra_raw_files,
    summarize_fastq_count_validation,
    validate_fastq_counts,
    validate_pct_q30,
    validate_processed_group,
    validate_read_metadata,
    validate_scale_cb_tag,
    validate_scale_processed_files,
    validate_scale_samples_csv,
    validate_scale_workflow_info,
)
from qa_mods import parse_met_summ, parse_web_summ

log = logging.getLogger(__name__)

SCALE = "scale"
SEAHUB = "seahub_sci"


class Status(str, Enum):
    """Outcome of one check.

    ``str`` mixin so these serialize straight into the JSON summary.
    """

    PASS = "pass"
    FINDINGS = "findings"
    SKIPPED = "skipped"
    NOT_RUN = "not_run"
    ERROR = "error"

    @property
    def denies_verdict(self) -> bool:
        """True when this outcome means the order cannot be called clean.

        SKIPPED does not: a Scale check on a 10x order is not a gap. NOT_RUN and
        ERROR do -- an applicable check that did not happen leaves a question
        open, and the whole point of this tool is that open questions are loud.
        """
        return self in (Status.NOT_RUN, Status.ERROR)


class CheckNotRun(Exception):
    """Raised by a check that applies but could not be performed."""


@dataclass(frozen=True)
class Scope:
    """The gating axes every notebook cell branches on."""

    raw_assay: str
    validate_raw: bool
    validate_processed: bool
    has_untrimmed_roots: bool = False

    @classmethod
    def resolve(
        cls,
        *,
        raw_assay: str,
        validate_raw: bool,
        validate_processed: bool,
        has_untrimmed_roots: bool = False,
    ) -> Scope:
        """Build a scope, applying the assay-level rules the notebook applies.

        SeaHub lab uploads are raw-only QA: the notebook forces
        ``validate_processed = False`` for ``seahub_sci`` in its gather cell.
        Enforced here rather than as an ``except_assays`` on each processed check,
        because it is a fact about the assay and not about any one check -- listing
        it per-check means the next processed check added inherits the wrong
        default and silently runs against a SeaHub upload.
        """
        if raw_assay == SEAHUB:
            validate_processed = False
        return cls(
            raw_assay=raw_assay,
            validate_raw=validate_raw,
            validate_processed=validate_processed,
            has_untrimmed_roots=has_untrimmed_roots,
        )

    def describe(self) -> str:
        return (
            f"assay={self.raw_assay} raw={self.validate_raw} "
            f"processed={self.validate_processed}"
        )


@dataclass
class CheckContext:
    """Everything the checks read, plus the state they hand each other."""

    data: Any
    scope: Scope
    bucket: str
    order_label: str
    s3_client: Any
    output_dir: str
    allow_truncated_stats_name: bool = False
    # Passed between checks: raw_expected_files computes the recognised set that
    # raw_extra_files needs, and read_metadata computes the per-group read counts
    # the CellRanger check cross-references. Threaded explicitly rather than left
    # in a shared namespace, which is what made the notebook's cell order load
    # bearing.
    shared: dict[str, Any] = field(default_factory=dict)


@dataclass
class CheckOutcome:
    """What a check produces."""

    findings: list[str] = field(default_factory=list)
    summary: str = ""
    notes: list[str] = field(default_factory=list)
    tables: dict[str, list[dict]] = field(default_factory=dict)


@dataclass
class CheckResult:
    """A declared check plus what became of it."""

    name: str
    title: str
    status: Status
    summary: str = ""
    reason: str = ""
    findings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    tables: dict[str, list[dict]] = field(default_factory=dict)

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "title": self.title,
            "status": self.status.value,
            "summary": self.summary,
            "reason": self.reason,
            "finding_count": len(self.findings),
            "findings": self.findings,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Check:
    """One declared check and the conditions under which it applies."""

    name: str
    title: str
    run: Callable[[CheckContext], CheckOutcome]
    requires_raw: bool = False
    requires_processed: bool = False
    only_assays: frozenset[str] | None = None
    except_assays: frozenset[str] = frozenset()
    implemented: bool = True
    unimplemented_reason: str = ""

    def applicability(self, scope: Scope) -> tuple[bool, str]:
        """Does this check apply to this order, and if not, why not?"""
        if self.only_assays is not None and scope.raw_assay not in self.only_assays:
            wanted = ", ".join(sorted(self.only_assays))
            return False, f"applies to {wanted} only; this order is {scope.raw_assay}"
        if scope.raw_assay in self.except_assays:
            return False, f"not applicable to {scope.raw_assay}"
        if self.requires_raw and not scope.validate_raw:
            return False, "raw validation is off (no raw/ found, or --no-raw)"
        if self.requires_processed and not scope.validate_processed:
            return False, (
                "processed validation is off (no processed/ found, or --no-processed)"
            )
        return True, ""


# --------------------------------------------------------------------------
# The checks. One per notebook cell that produces a verdict; each delegates to
# the same qa_checks function the notebook calls.
# --------------------------------------------------------------------------


def _check_pct_q30(ctx: CheckContext) -> CheckOutcome:
    values = ctx.data.pct_q30_values
    if not values:
        raise CheckNotRun(
            "no PCT_PF_Q30_bases values were gathered; the aggregate raw stats "
            "CSVs were absent or unreadable, so sublibrary Q30 is unknown"
        )
    findings = validate_pct_q30(values, MIN_PCT_Q30_THRESHOLD)
    return CheckOutcome(
        findings=findings,
        summary=(
            f"{len(values) - len(findings)}/{len(values)} sublibraries at or above "
            f"{MIN_PCT_Q30_THRESHOLD}% Q30"
        ),
    )


def _check_fastq_counts(ctx: CheckContext) -> CheckOutcome:
    if not ctx.data.fastq_log:
        raise CheckNotRun(
            "no FASTQ groups were gathered, so there was nothing to compare "
            "counts across -- the validator itself reports 'nothing to compare' "
            "here, which is not the same as counts that agree"
        )
    findings = validate_fastq_counts(ctx.data.fastq_log, ctx.scope.raw_assay)
    summary = summarize_fastq_count_validation(
        ctx.data.fastq_log, ctx.scope.raw_assay, findings
    )
    return CheckOutcome(findings=findings, summary=summary)


def _check_read_metadata(ctx: CheckContext) -> CheckOutcome:
    if not ctx.data.read_metadata:
        raise CheckNotRun(
            "no metadata.json objects were gathered, so R1/R2 read counts could "
            "not be compared for any group"
        )
    group_read_counts, findings, pairing = validate_read_metadata(
        ctx.data.read_metadata,
        ctx.scope.raw_assay,
        print_success=False,
    )
    # The CellRanger check cross-references these, so publish them even when this
    # check found problems -- a partial mapping still constrains more than none.
    ctx.shared["group_read_counts"] = group_read_counts
    return CheckOutcome(
        findings=findings,
        summary=(
            f"{len(group_read_counts)} of {len(ctx.data.read_metadata)} "
            "metadata entries yielded usable read counts"
        ),
        tables={"read_pairing": _as_rows(pairing)},
    )


def _as_rows(value: Any) -> list[dict]:
    """Coerce a validator's auxiliary return into report rows, or drop it.

    These shapes vary by assay and are advisory, so an unexpected one is not worth
    failing a run over.
    """
    if isinstance(value, dict):
        return [{"key": str(k), "value": str(v)} for k, v in value.items()]
    if isinstance(value, list):
        return [v if isinstance(v, dict) else {"value": str(v)} for v in value]
    return []


def _check_expected_raw_files(ctx: CheckContext) -> CheckOutcome:
    if not ctx.data.all_raw_files:
        raise CheckNotRun("no raw objects were listed, so nothing could be compared")
    all_good, raw_lost, raw_found = check_expected_raw_files(
        ctx.data.all_raw_files,
        ctx.scope.raw_assay,
        allow_truncated_stats_name=ctx.allow_truncated_stats_name,
    )
    total = all_good + len(raw_lost)
    if total == 0:
        # Zero recognised groups is not zero missing files. The convention-driven
        # check found nothing it knew how to have an opinion about -- a vendor
        # using unexpected naming, or the wrong --assay -- and reporting that as
        # "complete" turned an unreadable order into four green checks.
        raise CheckNotRun(
            f"no {ctx.scope.raw_assay} file groups were recognised among "
            f"{len(ctx.data.all_raw_files)} raw object(s), so completeness could "
            "not be assessed -- check --assay and the upload's naming"
        )
    ctx.shared["raw_found"] = raw_found
    findings = [f"MISSING RAW: {row}" for row in raw_lost]
    return CheckOutcome(
        findings=findings,
        summary=f"{all_good}/{total} file groups complete",
        tables={
            "raw_missing": [
                r if isinstance(r, dict) else {"item": str(r)} for r in raw_lost
            ]
        },
    )


def _check_extra_raw_files(ctx: CheckContext) -> CheckOutcome:
    if "raw_found" not in ctx.shared:
        raise CheckNotRun(
            "the expected-raw-files check did not produce a recognised-file set, "
            "so an unexpected file cannot be told from an expected one"
        )
    extra = check_extra_raw_files(
        ctx.data.all_raw_files,
        ctx.shared["raw_found"],
        ctx.scope.raw_assay,
        allow_truncated_stats_name=ctx.allow_truncated_stats_name,
    )
    return CheckOutcome(
        # Extras are reported, not failed on: an order routinely carries vendor
        # extras that are harmless. They still belong in the report, because an
        # unexpected file is sometimes a misnamed expected one.
        findings=[],
        summary=f"{len(extra)} unexpected raw object(s)",
        notes=[f"EXTRA RAW: {name}" for name in extra],
        tables={"raw_extra": [{"key": name} for name in extra]},
    )


def _check_cellranger_processed(ctx: CheckContext) -> CheckOutcome:
    groups = ctx.data.all_proc_files
    if not groups:
        raise CheckNotRun(
            "processed/ was detected but no per-group files were gathered"
        )
    group_read_counts = ctx.shared.get("group_read_counts", {})
    findings: list[str] = []
    notes: list[str] = []
    missing_rows: list[dict] = []
    alert_rows: list[dict] = []
    extra_rows: list[dict] = []

    if not group_read_counts:
        notes.append(
            "no per-group read counts available (raw validation off or metadata "
            "unreadable), so the reads-vs-CellRanger cross-check is weaker"
        )

    # The notebook downloads into the working directory under fixed filenames.
    # A temp dir instead: two concurrent runs would otherwise overwrite each
    # other's web_summary.html and QA the wrong order's metrics.
    with tempfile.TemporaryDirectory(prefix="order-qa-cellranger-") as scratch:
        for group, proc_files in groups.items():
            report: dict[str, Any] = {}
            for basename, parser in (
                ("web_summary.html", parse_web_summ),
                ("metrics_summary.csv", parse_met_summ),
            ):
                key = next(
                    (f for f in proc_files if f.rsplit("/", 1)[-1] == basename), None
                )
                if key is None:
                    continue
                local = os.path.join(scratch, f"{group}_{basename}")
                try:
                    ctx.s3_client.download_file(ctx.bucket, key, local)
                    report.update(parser(local))
                except Exception as exc:  # noqa: BLE001 - one group, not the run
                    notes.append(f"{group}: could not read {basename}: {exc}")

            result = validate_processed_group(
                group, proc_files, report, group_read_counts
            )
            findings.extend(result["errors"])
            alert_rows.extend(
                a if isinstance(a, dict) else {"alert": str(a)}
                for a in result["alerts"]
            )
            missing_rows.extend(
                m if isinstance(m, dict) else {"item": str(m)}
                for m in result["proc_missing"]
            )
            if result["process_extra"]:
                extra_rows.append(
                    {"group": group, "extra": ",".join(result["process_extra"])}
                )
            if "Probe Set Name" not in report:
                notes.append(
                    f"TRANSCRIPTOME: {group} used {report.get('Transcriptome', 'N/A')}"
                )

    return CheckOutcome(
        findings=findings,
        summary=f"{len(groups)} processed group(s) checked",
        notes=notes,
        tables={
            "process_missing": missing_rows,
            "process_alerts": alert_rows,
            "process_extra": extra_rows,
        },
    )


def _check_scale_cb_tag(ctx: CheckContext) -> CheckOutcome:
    if not ctx.data.read_metadata:
        raise CheckNotRun("no CRAM metadata entries were gathered")
    findings = validate_scale_cb_tag(ctx.data.read_metadata)
    return CheckOutcome(
        findings=findings,
        summary=(
            f"{len(ctx.data.read_metadata) - len(findings)}/"
            f"{len(ctx.data.read_metadata)} CRAM metadata entries have cb_tag=True"
        ),
    )


def _check_scale_workflow(ctx: CheckContext) -> CheckOutcome:
    infos = ctx.data.scale_workflow_infos
    samples = ctx.data.scale_samples_info
    if not infos and not samples:
        raise CheckNotRun(
            "neither workflow_info.json nor samples.csv was gathered for any group"
        )
    findings: list[str] = []
    notes: list[str] = []
    for group, params in infos.items():
        group_findings, info = validate_scale_workflow_info(params)
        findings.extend(group_findings)
        notes.extend(f"{group}: {msg}" for msg in info)
    for group, sinfo in samples.items():
        csv_findings = validate_scale_samples_csv(sinfo["columns"])
        findings.extend(csv_findings)
        if not csv_findings:
            notes.append(f"{group}: samples.csv OK ({len(sinfo['samples'])} samples)")
    return CheckOutcome(
        findings=findings,
        summary=f"{len(infos)} workflow_info, {len(samples)} samples.csv checked",
        notes=notes,
    )


def _check_scale_processed_files(ctx: CheckContext) -> CheckOutcome:
    groups = ctx.data.scale_proc_files
    if not groups:
        raise CheckNotRun("no Scale processed groups were gathered")
    findings: list[str] = []
    missing_rows: list[dict] = []
    for group, proc_files in groups.items():
        sinfo = ctx.data.scale_samples_info.get(group, {})
        result = validate_scale_processed_files(proc_files, sinfo)
        findings.extend(result["errors"])
        missing_rows.extend(
            {"group": group, "missing": str(m)} for m in result["missing_files"]
        )
    return CheckOutcome(
        findings=findings,
        summary=f"{len(groups)} Scale processed group(s) checked",
        tables={"scale_processed_missing": missing_rows},
    )


def _check_seahub_sop(ctx: CheckContext) -> CheckOutcome:
    violations = ctx.data.sop_violations
    return CheckOutcome(
        findings=[
            f"SOP {v.get('type', 'violation')}: {v.get('detail', v)}"
            for v in violations
        ],
        summary=f"{len(violations)} SOP naming violation(s)",
        tables={"sop_violations": [dict(v) for v in violations]},
    )


def _unimplemented(ctx: CheckContext) -> CheckOutcome:  # pragma: no cover
    raise CheckNotRun("not implemented in this CLI")


REGISTRY: tuple[Check, ...] = (
    Check(
        name="raw_expected_files",
        title="Expected raw files present",
        run=_check_expected_raw_files,
        requires_raw=True,
    ),
    Check(
        name="raw_extra_files",
        title="Unexpected raw files",
        run=_check_extra_raw_files,
        requires_raw=True,
    ),
    Check(
        name="raw_fastq_counts",
        title="FASTQ counts consistent across sublibrary types",
        run=_check_fastq_counts,
        requires_raw=True,
    ),
    Check(
        name="raw_read_metadata",
        title="R1/R2 read counts agree within each group",
        run=_check_read_metadata,
        requires_raw=True,
    ),
    Check(
        name="raw_pct_q30",
        title="Sublibrary Q30 at or above threshold",
        run=_check_pct_q30,
        requires_raw=True,
    ),
    Check(
        name="scale_cb_tag",
        title="Scale CRAMs carry cb_tag",
        run=_check_scale_cb_tag,
        requires_raw=True,
        only_assays=frozenset({SCALE}),
    ),
    Check(
        name="processed_cellranger",
        title="CellRanger run metadata and outputs",
        run=_check_cellranger_processed,
        requires_processed=True,
        # Scale processed data is validated by the three Scale checks instead.
        except_assays=frozenset({SCALE}),
    ),
    Check(
        name="scale_workflow_info",
        title="Scale workflow_info.json and samples.csv",
        run=_check_scale_workflow,
        requires_processed=True,
        only_assays=frozenset({SCALE}),
    ),
    Check(
        name="scale_processed_files",
        title="Scale per-sample and per-sublibrary outputs",
        run=_check_scale_processed_files,
        requires_processed=True,
        only_assays=frozenset({SCALE}),
    ),
    Check(
        name="seahub_sop_naming",
        title="SeaHub SOP path and filename conventions",
        run=_check_seahub_sop,
        requires_raw=True,
        only_assays=frozenset({SEAHUB}),
    ),
    # Declared but not built. A seahub_sci order therefore reports these as
    # NOT_RUN and exits degraded, rather than passing on nine checks while
    # silently omitting the cross-bucket reconciliation that is the entire point
    # of QA-ing a trimmed upload. Declaring it is the difference between a known
    # gap and an invisible one.
    Check(
        name="seahub_trimming_completeness",
        title="Trimmed upload reconciled against untrimmed vendor sources",
        run=_unimplemented,
        requires_raw=True,
        only_assays=frozenset({SEAHUB}),
        implemented=False,
        unimplemented_reason=(
            "cross-bucket reconciliation (index_trimmed_upload / "
            "discover_untrimmed_sources / reconcile_trimming) is not wired into "
            "the CLI yet -- run qa.ipynb for seahub_sci orders"
        ),
    ),
    Check(
        name="seahub_well_status",
        title="Per-well verdict and corrected names",
        run=_unimplemented,
        requires_raw=True,
        only_assays=frozenset({SEAHUB}),
        implemented=False,
        unimplemented_reason=(
            "depends on the reconciliation above -- run qa.ipynb for seahub_sci orders"
        ),
    ),
)


def run_registry(ctx: CheckContext) -> list[CheckResult]:
    """Run every declared check, in order, and return one result per declaration.

    The returned list always has one entry per registry entry. That invariant is
    the product: a caller can compare declared against reported and get zero, so
    there is no way for a check to go missing between here and the report.
    """
    results: list[CheckResult] = []
    for check in REGISTRY:
        applicable, reason = check.applicability(ctx.scope)
        if not applicable:
            log.info("skip  %-32s %s", check.name, reason)
            results.append(
                CheckResult(
                    name=check.name,
                    title=check.title,
                    status=Status.SKIPPED,
                    reason=reason,
                )
            )
            continue
        if not check.implemented:
            log.warning("NOT RUN %-30s %s", check.name, check.unimplemented_reason)
            results.append(
                CheckResult(
                    name=check.name,
                    title=check.title,
                    status=Status.NOT_RUN,
                    reason=check.unimplemented_reason,
                )
            )
            continue
        results.append(_run_one(check, ctx))
    return results


def _run_one(check: Check, ctx: CheckContext) -> CheckResult:
    try:
        outcome = check.run(ctx)
    except CheckNotRun as exc:
        log.warning("NOT RUN %-30s %s", check.name, exc)
        return CheckResult(
            name=check.name,
            title=check.title,
            status=Status.NOT_RUN,
            reason=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 - one check must not end the run
        log.exception("ERROR %-32s %s", check.name, exc)
        return CheckResult(
            name=check.name,
            title=check.title,
            status=Status.ERROR,
            reason=f"{type(exc).__name__}: {exc}",
        )
    status = Status.FINDINGS if outcome.findings else Status.PASS
    log.info(
        "%-5s %-32s %s",
        "FIND" if outcome.findings else "ok",
        check.name,
        outcome.summary,
    )
    return CheckResult(
        name=check.name,
        title=check.title,
        status=status,
        summary=outcome.summary,
        findings=outcome.findings,
        notes=outcome.notes,
        tables=outcome.tables,
    )
