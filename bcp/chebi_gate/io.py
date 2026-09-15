"""Writing a run's outputs: two SDFs, a findings table, open questions, a manifest.

Every file is always written, even when empty, so a downstream step can rely on it
existing rather than branching on whether the gate found anything.

One output exists purely for the human in the re-run workflow. The decisions data
is keyed on CAS, which is right for joining and unhelpful for reading: nobody
remembers which compound 23256-33-9 is. So the gate writes ``open_questions.csv``,
which carries the compound name alongside the key. Display text belongs in
generated artifacts; input data keeps the key.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from . import sdf as sdf_mod
from .checks import CHECKS, SEVERITY_RANK
from .client import STATUS_HELD, GateRun
from .manifest import MANIFEST_FILENAME

log = logging.getLogger(__name__)

CLEARED_SUFFIX = "_cleared.sdf"
HELD_SUFFIX = "_held.sdf"
FINDINGS_FILENAME = "findings.csv"
OPEN_QUESTIONS_FILENAME = "open_questions.csv"

FINDINGS_COLUMNS = (
    "record",
    "cas",
    "name",
    "check",
    "severity",
    "independent",
    "holds",
    "waived",
    "detail",
    "evidence",
)
OPEN_QUESTIONS_COLUMNS = ("cas", "name", "checks", "detail")


@dataclass(frozen=True)
class Outputs:
    """Where a run's outputs were written."""

    cleared: Path
    held: Path
    findings: Path
    open_questions: Path
    manifest: Path


def write(gate_run: GateRun, out_dir: str | Path, *, stem: str) -> Outputs:
    """Write every output for one run and return their paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cleared_path = out_dir / f"{stem}{CLEARED_SUFFIX}"
    held_path = out_dir / f"{stem}{HELD_SUFFIX}"

    sdf_mod.write_records([r.record for r in gate_run.cleared], cleared_path)
    sdf_mod.write_annotated(
        [(r.record, _annotation(gate_run, r)) for r in gate_run.held], held_path
    )

    findings_path = _write_findings(gate_run, out_dir / FINDINGS_FILENAME)
    questions_path = _write_open_questions(gate_run, out_dir / OPEN_QUESTIONS_FILENAME)
    manifest_path = gate_run.manifest.write(out_dir / MANIFEST_FILENAME)

    return Outputs(
        cleared=cleared_path,
        held=held_path,
        findings=findings_path,
        open_questions=questions_path,
        manifest=manifest_path,
    )


def _annotation(gate_run: GateRun, result) -> dict[str, str]:
    """The GATE_* fields for one held record.

    ``GATE_RUN`` carries the run identifier, which is derived from the inputs
    rather than the clock -- so two runs over the same material annotate
    identically and the byte-identity invariant is testable.
    """
    findings = sorted(
        result.findings + gate_run.file_findings,
        key=lambda f: (SEVERITY_RANK[f.severity], f.check),
    )
    reasons = "\n".join(
        f"{f.check} [{f.severity}]"
        f"{' HOLDS' if f in result.blocking else ''}"
        f"{' independent' if f.is_independent else ' circular'} {f.detail}"
        for f in findings
    )
    evidence = "\n".join(sorted({f.evidence for f in findings if f.evidence}))
    held_by = sorted({f.check for f in result.blocking})
    return {
        "GATE_STATUS": STATUS_HELD,
        "GATE_CHECKS_FAILED": "; ".join(held_by) or "file-level",
        "GATE_REASONS": reasons,
        "GATE_EVIDENCE": evidence or "none",
        "GATE_RUN": f"{gate_run.manifest.run_id} role={gate_run.manifest.checks['role']}"
        f" medium={'holds' if not gate_run.allow_medium else 'allowed'}",
    }


def _write_findings(gate_run: GateRun, path: Path) -> Path:
    """Every finding, including the ones that never hold a record.

    Sorted by severity then check then record, so two runs over the same inputs
    write the same bytes: dictionary iteration order was one of the two things
    that made the prototype's outputs differ between identical runs.
    """
    rows = []
    for result in gate_run.results:
        for finding in result.findings:
            rows.append(
                {
                    "record": finding.record_index,
                    "cas": finding.cas,
                    "name": finding.name,
                    "check": finding.check,
                    "severity": finding.severity,
                    "independent": "yes" if finding.is_independent else "no",
                    "holds": "yes" if finding in result.blocking else "no",
                    "waived": finding.waiver or "",
                    "detail": finding.detail,
                    "evidence": finding.evidence,
                }
            )
    for finding in gate_run.file_findings:
        rows.append(
            {
                "record": 0,
                "cas": "",
                "name": "(whole file)",
                "check": finding.check,
                "severity": finding.severity,
                "independent": "yes" if finding.is_independent else "no",
                "holds": "yes",
                "waived": "",
                "detail": finding.detail,
                "evidence": finding.evidence,
            }
        )
    rows.sort(
        key=lambda r: (
            SEVERITY_RANK[r["severity"]],
            r["check"],
            r["record"],
            r["detail"],
        )
    )
    # backslashreplace, not strict: a finding's detail can echo a non-ASCII byte
    # from the record it is about, carried as a lone surrogate. Encoding strictly
    # raised UnicodeEncodeError part-way through writerows, leaving a truncated
    # findings.csv that still looked like valid CSV, no run manifest at all, and
    # an exit status of 1 that is indistinguishable from "records were held".
    with path.open(
        "w", newline="", encoding="utf-8", errors="backslashreplace"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FINDINGS_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_open_questions(gate_run: GateRun, path: Path) -> Path:
    """Held records with an open question, with the compound name attached.

    The decisions data is keyed on CAS because a name changes and a key must not.
    This file is the other half of that trade: the human resolving a question
    needs to know which compound it is.
    """
    rows = []
    for result in gate_run.results:
        questions = [f for f in result.findings if f.check == "QUAR-01"]
        if not questions:
            continue
        rows.append(
            {
                "cas": result.context.cas,
                "name": result.context.name,
                "checks": "; ".join(
                    sorted({q.evidence for q in questions if q.evidence})
                ),
                "detail": " | ".join(q.detail for q in questions),
            }
        )
    rows.sort(key=lambda r: r["cas"])
    with path.open(
        "w", newline="", encoding="utf-8", errors="backslashreplace"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OPEN_QUESTIONS_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)
    return path


def summary(gate_run: GateRun, outputs: Outputs | None = None) -> str:
    """A human-readable report, with every figure computed from the run.

    No number here is written by hand. A hard-coded "122 structures" survived into
    a cover note when the real figure was 131, so the totals are derived and the
    independent and circular counts are printed separately -- adding them together
    is the mistake this whole design exists to prevent.
    """
    counts = gate_run.counts()
    independence = gate_run.independence()
    reference = gate_run.manifest.reference
    lines = [
        f"run {gate_run.manifest.run_id} at {gate_run.manifest.timestamp}",
        f"  records  {counts['records']:>5}",
        f"  cleared  {counts['cleared']:>5}",
        f"  held     {counts['held']:>5}",
    ]
    if outputs is not None:
        lines += [
            f"    cleared -> {outputs.cleared}",
            f"    held    -> {outputs.held}",
        ]

    available = reference.get("sources_available") or []
    independent = reference.get("independent_sources") or []
    lines.append("")
    lines.append(f"  reference sources: {', '.join(available) or 'none'}")
    lines.append(f"  of which independent: {', '.join(independent) or 'none'}")
    if not independent:
        lines.append(
            "  NOTE: with no independent source, clearance means internally "
            "consistent, not confirmed."
        )
    registry = reference.get("cas_registry")
    if registry:
        coverage = registry["coverage"]
        lines.append(
            f"  CAS registry rows for {coverage['with_registry_row']}"
            f"/{coverage['records']} records"
            f" ({coverage['without_registry_row']} unchecked on stoichiometry)"
        )

    lines.append("")
    lines.append(
        f"  findings {counts['findings']:>5}"
        f"   independent {independence['independent']}"
        f", circular {independence['circular']}"
    )
    lines.append(
        f"  of those, holding a record: independent "
        f"{independence['independent_holding']}, circular "
        f"{independence['circular_holding']}"
    )

    by_check = gate_run.by_check()
    if by_check:
        lines.append("")
        lines.append("  findings by check:")
        for key, count in sorted(
            by_check.items(),
            key=lambda kv: (SEVERITY_RANK[kv[0].split("/")[1]], kv[0]),
        ):
            check_id, severity = key.split("/")
            title = CHECKS[check_id].title
            lines.append(f"    {check_id:8s} {severity:7s} {count:>4}  {title}")

    if gate_run.unmatched_decisions:
        lines.append("")
        lines.append("  decisions matching no record in this input:")
        for message in gate_run.unmatched_decisions:
            lines.append(f"    {message}")

    return "\n".join(lines)
