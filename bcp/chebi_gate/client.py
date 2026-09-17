"""The gate: one SDF in, a cleared file and a held file out, with reasons.

    cleared   records that passed every check. Byte-identical to the input
              records, no fields added, submittable as they stand.
    held      records that did not. Each carries GATE_* fields saying which
              checks failed, why, and what the evidence was.

The held file is an SDF rather than a report so that a fixed record can be fed
straight back in; the cleared file is left untouched precisely so it stays
submittable, since annotating it would break the field template ChEBI loads it
with.

Two invariants are asserted rather than hoped for. Cleared plus held always equals
the input count, and a cleared record is byte-identical to its input. Running the
gate twice over the same inputs produces byte-identical outputs, which is possible
only because the run identifier is derived from the inputs instead of the clock.

Severity decides what holds a record: ``high`` always, ``medium`` unless the caller
allows it, ``low`` and ``info`` never. On top of that, an unresolved question in
the decisions data holds a record whatever the checks say -- "which salt does the
lab hold?" is not answerable by looking at the file.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path

from . import checks as checks_mod
from . import external as external_mod
from . import manifest as manifest_mod
from . import sdf as sdf_mod
from . import structure as structure_mod
from .checks import HIGH, INFO, LOW, MEDIUM, Finding, RecordContext
from .decisions import Decisions

log = logging.getLogger(__name__)

GATE_FIELDS = (
    "GATE_STATUS",
    "GATE_CHECKS_FAILED",
    "GATE_REASONS",
    "GATE_EVIDENCE",
    "GATE_RUN",
)

STATUS_HELD = "HELD"


class GateError(Exception):
    """Raised when the gate cannot produce a trustworthy verdict."""


@dataclass
class RecordResult:
    """One record's verdict and everything that produced it."""

    context: RecordContext
    findings: list[Finding] = field(default_factory=list)
    cleared: bool = True
    # Findings that held this record. Stored rather than recomputed, because
    # whether a medium holds depends on how the run was invoked, and a property
    # that silently guessed would disagree with the file that was written.
    blocking: list[Finding] = field(default_factory=list)

    @property
    def record(self) -> sdf_mod.SdfRecord:
        return self.context.record


@dataclass
class GateRun:
    """The outcome of one run over one input file."""

    results: list[RecordResult]
    file_findings: list[Finding]
    manifest: manifest_mod.Manifest
    unmatched_decisions: list[str]
    allow_medium: bool = False

    @property
    def cleared(self) -> list[RecordResult]:
        return [r for r in self.results if r.cleared]

    @property
    def held(self) -> list[RecordResult]:
        return [r for r in self.results if not r.cleared]

    @property
    def all_findings(self) -> list[Finding]:
        return [f for r in self.results for f in r.findings] + self.file_findings

    def counts(self) -> dict[str, int]:
        return {
            "records": len(self.results),
            "cleared": len(self.cleared),
            "held": len(self.held),
            "findings": len(self.all_findings),
        }

    def by_check(self) -> dict[str, int]:
        return dict(
            sorted(
                Counter(f"{f.check}/{f.severity}" for f in self.all_findings).items()
            )
        )

    def independence(self) -> dict[str, int]:
        """Findings split by whether their evidence had an independent origin.

        Reported separately and never added together. A circular check passed 284
        of 290 records on the reference batch, and the number mostly says the
        generation step ran.
        """
        # "Holding" has to mean what this run actually did, not what the
        # severities would do by default: with --allow-medium a medium finding is
        # reported and does not hold, and a summary saying 6 findings held a
        # record while 0 records were held is worse than no summary.
        holding = [
            f
            for f in self.all_findings
            if not severity_clears(f, allow_medium=self.allow_medium)
        ]
        return {
            "independent": sum(1 for f in self.all_findings if f.is_independent),
            "circular": sum(1 for f in self.all_findings if not f.is_independent),
            "independent_holding": sum(1 for f in holding if f.is_independent),
            "circular_holding": sum(1 for f in holding if not f.is_independent),
        }


def severity_clears(finding: Finding, *, allow_medium: bool) -> bool:
    if finding.severity in (LOW, INFO):
        return True
    if finding.severity == MEDIUM and allow_medium:
        return True
    return False


def run(
    input_path: str | Path,
    *,
    evidence: external_mod.Evidence | None = None,
    decisions: Decisions | None = None,
    role: str = checks_mod.ROLE_AUTO,
    allow_medium: bool = False,
    repo: str | Path | None = None,
    timestamp: str | None = None,
) -> GateRun:
    """Judge every record in one SDF. Does not write anything; see :mod:`.io`."""
    if role not in checks_mod.ROLES:
        raise GateError(f"unknown role {role!r}; expected one of {checks_mod.ROLES}")

    evidence = evidence or external_mod.Evidence()
    input_path = Path(input_path)
    parsed = sdf_mod.parse_file(input_path)
    if not parsed.records:
        raise GateError(f"no records found in {input_path}")

    if parsed.malformed:
        # Refused, not worked around. If a record boundary is in doubt then so is
        # every finding attributed to a record, and re-emitting invents bytes: a
        # truncated input used to clear silently, with a cleared file five bytes
        # longer than its input because the gate supplied the missing terminator.
        raise GateError(
            f"{input_path} is not a well-formed SDF: "
            + "; ".join(parsed.malformed)
            + ". Fix the file before gating it."
        )

    # A held record carries GATE_* fields, and the whole point of writing the held
    # file as an SDF is that a fixed record goes straight back in. It could not:
    # none of the GATE_* tags is in EXPECTED_TAGS, so INT-02 reported "unexpected
    # tags" at medium and held every record of the held file -- measured, 193 of
    # 193 -- and with --allow-medium the record cleared instead and tripped the
    # byte-identity invariant, which aborts the whole run with no outputs written
    # at all. The documented five-step recipe could not complete under either
    # policy.
    #
    # Stripping them on load is the resolution that keeps invariant 2 rather than
    # weakening it: the annotations are removed, not tolerated, so a cleared record
    # still carries none and an annotation still cannot be laundered through a
    # re-run. What changes is the baseline for "byte-identical" -- it is the input
    # with any previous run's annotation taken off, which is what the record was
    # before the gate touched it.
    annotated = [r for r in parsed.records if any(t in r.data for t in GATE_FIELDS)]
    if annotated:
        log.info(
            "%d of %d records arrived carrying a previous run's GATE_* annotation; "
            "stripping it before judging them",
            len(annotated),
            len(parsed.records),
        )
    records = (
        [_strip_annotation(r) for r in parsed.records]
        if annotated
        else (parsed.records)
    )

    contexts = [
        RecordContext(record=record, structure=structure_mod.analyse(record), role=role)
        for record in records
    ]

    results = [RecordResult(context=ctx) for ctx in contexts]
    for result in results:
        local = checks_mod.run_record_checks(result.context)
        remote = external_mod.run_external_checks(evidence, result.context)
        result.findings = local + remote

    file_findings = checks_mod.run_file_checks(
        checks_mod.FileContext(contexts=contexts, raw=parsed.raw)
    )
    by_index: dict[int, list[Finding]] = defaultdict(list)
    whole_file: list[Finding] = []
    for finding in file_findings:
        if finding.record_index:
            by_index[finding.record_index].append(finding)
        else:
            whole_file.append(finding)
    for result in results:
        result.findings.extend(by_index.get(result.record.index, []))

    unmatched: list[str] = []
    if decisions is not None:
        applied = _apply_decisions(results, decisions)
        unmatched = decisions.unmatched(
            {c.cas_key for c in contexts if c.cas_key}, applied
        )
        for message in unmatched:
            log.warning("%s", message)

    # A whole-file problem holds every record, because the file as a whole is not
    # submittable -- but it is counted once, not once per record.
    file_blocks = any(
        not severity_clears(f, allow_medium=allow_medium) for f in whole_file
    )
    for result in results:
        blocking = [
            f
            for f in result.findings
            if not severity_clears(f, allow_medium=allow_medium)
        ]
        result.blocking = blocking
        result.cleared = not blocking and not file_blocks

    run_manifest = _build_manifest(
        input_path=input_path,
        parsed=parsed,
        evidence=evidence,
        decisions=decisions,
        role=role,
        allow_medium=allow_medium,
        repo=repo,
        timestamp=timestamp,
    )
    gate_run = GateRun(
        results=results,
        file_findings=whole_file,
        manifest=run_manifest,
        unmatched_decisions=unmatched,
        allow_medium=allow_medium,
    )
    run_manifest.results = {
        **gate_run.counts(),
        "by_check": gate_run.by_check(),
        "independence": gate_run.independence(),
        "unmatched_decisions": unmatched,
    }
    _assert_invariants(gate_run, parsed)
    return gate_run


def _strip_annotation(record: sdf_mod.SdfRecord) -> sdf_mod.SdfRecord:
    """One record with any GATE_* field removed, re-parsed from the stripped bytes.

    Re-parsed rather than edited in place so every span, tag order and value comes
    from the same reader the rest of the gate uses -- the alternative is a second
    notion of where a field ends, which is the defect the byte-span removal was
    introduced to fix.
    """
    if not any(tag in record.data for tag in GATE_FIELDS):
        return record
    stripped = record.without_fields(GATE_FIELDS)
    reparsed = sdf_mod.parse_bytes(stripped + sdf_mod.RECORD_TERMINATOR)
    if not reparsed.records:  # pragma: no cover - a record cannot strip to nothing
        return record
    return replace(
        reparsed.records[0], index=record.index, terminator=record.terminator
    )


def _apply_decisions(
    results: list[RecordResult], decisions: Decisions
) -> set[tuple[str, str, str]]:
    """Downgrade waived findings, hold records with an open question.

    Returns the waiver keys that actually downgraded something, so the caller can
    report the ones that did not. Two reasons that used to be silent:

    The severity filter. Only findings that would have held the record were waived
    -- high and medium -- so a waiver written against a ``low`` or ``info`` finding
    did nothing at all and nothing said so, and five checks emit ``low`` as their
    only severity. A
    waiver is a recorded judgement about a finding; whether that finding would
    have held the record is a separate question, and not one the waiver asked.

    And a waiver that applies to nothing looks exactly like one that applies. The
    key includes the severity now, so a decision about a medium CON-04 no longer
    absorbs a high one; the cost of that narrowing is that a waiver can miss, and
    a miss has to be reported rather than inferred.
    """
    applied: set[tuple[str, str, str]] = set()
    for result in results:
        # The normalised key, so a decision keeps applying while a malformed CAS
        # spelling is being fixed. INT-04 reports the spelling independently.
        cas = result.context.cas_key
        waived: list[Finding] = []
        for finding in result.findings:
            waiver = (
                decisions.waiver_for(cas, finding.check, finding.severity)
                if cas
                else None
            )
            if waiver is not None:
                applied.add((cas, finding.check, finding.severity))
                waived.append(checks_mod.waive(finding, waiver.reason, waiver.evidence))
            else:
                waived.append(finding)
        result.findings = waived

        for question in decisions.open_questions(cas) if cas else ():
            result.findings.append(
                result.context.finding(
                    "QUAR-01",
                    HIGH,
                    f"open question since {question.opened_on.isoformat()}: "
                    f"{question.reason} Unblocked by: {question.blocked_on}",
                    evidence=", ".join(question.check_ids),
                )
            )
    return applied


def _build_manifest(
    *,
    input_path: Path,
    parsed: sdf_mod.SdfFile,
    evidence: external_mod.Evidence,
    decisions: Decisions | None,
    role: str,
    allow_medium: bool,
    repo: str | Path | None,
    timestamp: str | None,
) -> manifest_mod.Manifest:
    run_manifest = manifest_mod.Manifest()
    run_manifest.add_input(input_path, len(parsed.records))
    run_manifest.record_environment(repo)
    run_manifest.timestamp = timestamp or datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )

    cas_values = {
        r.cas for r in parsed.records if r.cas
    }  # raw, as the file spells them
    reference: dict = {
        "sources_available": list(evidence.available),
        "independent_sources": list(evidence.independent_sources),
    }
    if len(evidence.registry):
        reference["cas_registry"] = {
            "path": str(evidence.registry.path),
            "sha256": evidence.registry.sha256,
            "rows": len(evidence.registry),
            "source_pdfs": list(evidence.registry.source_pdfs),
            "coverage": evidence.registry.coverage(cas_values),
        }
    else:
        reference["cas_registry"] = None
    if len(evidence.chebi):
        reference["chebi_release"] = evidence.chebi.manifest
    else:
        reference["chebi_release"] = None
    for name, directory in (
        ("pubchem_cache", evidence.pubchem_dir),
        ("cas_common_chemistry_cache", evidence.common_chemistry_dir),
    ):
        if directory and Path(directory).is_dir():
            reference[name] = manifest_mod.directory_digest(directory, "*.json")
        else:
            reference[name] = None
    run_manifest.reference = reference

    if decisions is not None:
        run_manifest.decisions = {
            "files": [
                {
                    "path": str(path),
                    "sha256": manifest_mod.sha256_file(path),
                }
                for path in decisions.paths
                if Path(path).exists()
            ],
            "counts": decisions.counts,
        }
    else:
        run_manifest.decisions = {"files": [], "counts": {}}

    run_manifest.checks = {
        "role": role,
        "medium_holds": not allow_medium,
        "registered": {
            check_id: {
                "severities": list(check.severities),
                "independent": check.independent,
                "title": check.title,
            }
            for check_id, check in sorted(checks_mod.CHECKS.items())
        },
    }
    run_manifest.compute_run_id()
    return run_manifest


def _assert_invariants(gate_run: GateRun, parsed: sdf_mod.SdfFile) -> None:
    """Handoff invariants 19 and 20, checked at run time and not only in tests."""
    total = len(gate_run.cleared) + len(gate_run.held)
    if total != len(parsed.records):
        raise GateError(
            f"cleared {len(gate_run.cleared)} plus held {len(gate_run.held)} is "
            f"{total}, but the input had {len(parsed.records)} records"
        )
    for result in gate_run.cleared:
        # Invariant 19 is about annotation, so test for annotation. An earlier
        # version compared strip_fields(raw) against raw, which also normalises
        # trailing newlines -- so it fired on any input whose records end with a
        # single newline rather than a blank line, including the prototype's own
        # output. Checking the parsed field names says what was meant.
        present = [tag for tag in GATE_FIELDS if tag in result.record.data]
        if present:
            raise GateError(
                f"record {result.record.index} ({result.record.name}) cleared but "
                f"carries {present}; annotate only held records, and re-run the "
                "gate on the original input rather than on its own output"
            )
