"""The gate end to end: the two invariants, the two workflows, and the manifest.

Handoff invariants 19 and 20 are asserted here and also inside the gate itself,
because a wrong answer about which records cleared is not something to discover
from a test suite months later.

The two workflows the whole design exists for:

new batch
    a file arrives, the gate splits it into submittable and held, and both
    outcomes are explained and reproducible later.

re-run
    a question is answered, one row of the decisions data gets a resolution date,
    the gate runs again and the record clears. No code is edited.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path

import pytest

from chebi_gate import casreg, chebi_release, checks, client, decisions, external
from chebi_gate import io as gate_io
from chebi_gate import manifest as manifest_mod
from chebi_gate import sdf
from chebi_gate.cli import EXIT_HELD, EXIT_OK, EXIT_USAGE, main
from tests import test_chebi_gate_release as release_fixtures
from tests.chebi_gate_helpers import (
    ISA_MALEATE,
    neutral_record,
    salt_record,
    sdf_bytes,
)

CLEAN_CAS = "557-66-4"
OTHER_CAS = "64-17-5"


@pytest.fixture
def clean_input(tmp_path: Path) -> Path:
    """Two records with no holding finding at all."""
    path = tmp_path / "batch.sdf"
    path.write_bytes(
        sdf_bytes(
            salt_record(iupac="ethanamine hydrochloride"),
            neutral_record(cas=OTHER_CAS, iupac="ethanol"),
        )
    )
    return path


@pytest.fixture
def mixed_input(tmp_path: Path) -> Path:
    """One clean record and one with a high finding."""
    path = tmp_path / "mixed.sdf"
    path.write_bytes(
        sdf_bytes(
            salt_record(iupac="ethanamine;hydrochloride"),
            # A different molfile on purpose: the same one would make these two
            # records duplicate structures and INT-05 would hold both, which is
            # correct behaviour and not what this fixture is for.
            salt_record(
                name="ethylamine maleate",
                mol="amine_2hcl",
                cas=OTHER_CAS,
                relationship=ISA_MALEATE,
                iupac="ethanamine maleate",
            ),
        )
    )
    return path


def empty_decisions(tmp_path: Path) -> decisions.Decisions:
    directory = tmp_path / "decisions"
    directory.mkdir(exist_ok=True)
    for name, columns in (
        (decisions.WAIVERS_FILE, decisions.WAIVER_COLUMNS),
        (decisions.QUARANTINE_FILE, decisions.QUARANTINE_COLUMNS),
    ):
        with (directory / name).open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle).writerow(columns)
    return directory


def write_decisions(tmp_path: Path, *, waivers=(), quarantine=()) -> Path:
    directory = empty_decisions(tmp_path)
    for name, columns, rows in (
        (decisions.WAIVERS_FILE, decisions.WAIVER_COLUMNS, waivers),
        (decisions.QUARANTINE_FILE, decisions.QUARANTINE_COLUMNS, quarantine),
    ):
        with (directory / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(columns))
            writer.writeheader()
            writer.writerows(rows)
    return directory


# ------------------------------------------------------------- the invariants


def test_cleared_plus_held_always_equals_the_input_count(mixed_input):
    """Handoff invariant 20."""
    run = client.run(mixed_input, allow_medium=True)
    assert len(run.cleared) + len(run.held) == 2
    assert run.counts()["records"] == 2


def test_a_cleared_record_is_byte_identical_to_its_input(tmp_path, mixed_input):
    """Handoff invariant 19: no annotation fields, no drift of any kind."""
    run = client.run(mixed_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="mixed")

    source = {
        r.raw.split(b"\n", 1)[0]: r.raw for r in sdf.parse_file(mixed_input).records
    }
    cleared = sdf.parse_file(outputs.cleared).records
    assert cleared
    for record in cleared:
        assert source[record.raw.split(b"\n", 1)[0]] == record.raw
    assert b"GATE_" not in outputs.cleared.read_bytes()


def test_the_gate_refuses_to_clear_a_record_that_is_already_annotated(tmp_path):
    """Re-running on the gate's own held file must not launder an annotation."""
    record = sdf.parse_bytes(sdf_bytes(salt_record(iupac="x"))).records[0]
    annotated = record.with_fields({"GATE_STATUS": "HELD"}) + sdf.RECORD_TERMINATOR
    path = tmp_path / "already.sdf"
    path.write_bytes(annotated)

    with pytest.raises(client.GateError, match="GATE_STATUS"):
        client.run(path, allow_medium=True)


def test_every_output_file_is_written_even_when_empty(tmp_path, clean_input):
    run = client.run(clean_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="batch")
    assert not run.held
    for path in (
        outputs.cleared,
        outputs.held,
        outputs.findings,
        outputs.open_questions,
        outputs.manifest,
    ):
        assert path.exists(), path
    assert outputs.held.read_bytes() == b""


def test_two_runs_over_the_same_inputs_write_identical_bytes(tmp_path, mixed_input):
    """Handoff 4.4. Possible only because the run id comes from the inputs.

    The prototype wrote a wall-clock timestamp into every annotated record, so no
    two runs could agree and the invariant could not even be tested.
    """
    first = gate_io.write(
        client.run(mixed_input, allow_medium=True), tmp_path / "a", stem="m"
    )
    second = gate_io.write(
        client.run(mixed_input, allow_medium=True), tmp_path / "b", stem="m"
    )
    assert first.cleared.read_bytes() == second.cleared.read_bytes()
    assert first.held.read_bytes() == second.held.read_bytes()
    assert first.findings.read_text() == second.findings.read_text()

    left = json.loads(first.manifest.read_text())
    right = json.loads(second.manifest.read_text())
    assert left["run_id"] == right["run_id"]
    assert {k for k in left if left[k] != right[k]} <= {"timestamp"}


# ------------------------------------------------------------------ severity


def test_a_high_finding_always_holds_a_record(mixed_input):
    run = client.run(mixed_input, allow_medium=True)
    held = run.held
    assert len(held) == 1
    assert any(f.severity == checks.HIGH for f in held[0].blocking)


def test_a_medium_finding_holds_by_default_and_clears_when_allowed(tmp_path):
    path = tmp_path / "medium.sdf"
    path.write_bytes(
        sdf_bytes(salt_record(synonym="ethylamine hcl;CHEMBL1234", iupac="x"))
    )

    strict = client.run(path)
    assert len(strict.held) == 1
    assert any(f.check == "SYN-01" for f in strict.held[0].blocking)

    lenient = client.run(path, allow_medium=True)
    assert not lenient.held


def test_low_and_info_findings_never_hold_a_record(tmp_path):
    path = tmp_path / "low.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="ethanamine;hydrochloride")))
    run = client.run(path)
    assert not run.held
    assert any(f.check == "IUP-01" for f in run.all_findings)


def test_a_whole_file_problem_holds_every_record_but_is_counted_once(tmp_path):
    """The prototype appended it per record, so one stray byte held all 290."""
    raw = sdf_bytes(
        salt_record(iupac="x"),
        neutral_record(cas=OTHER_CAS, synonym="cafe_", iupac="y"),
    ).replace(b"cafe_", b"caf\xe9")
    path = tmp_path / "nonascii.sdf"
    path.write_bytes(raw)

    run = client.run(path, allow_medium=True)
    assert len(run.held) == 2
    assert len(run.file_findings) == 1
    assert run.file_findings[0].check == "INT-06"


# ------------------------------------------------------------- the two workflows


def test_new_batch_workflow_explains_every_held_record(tmp_path, mixed_input):
    """Each held record says which checks held it, why, and on what evidence."""
    run = client.run(mixed_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="mixed")

    held = sdf.parse_file(outputs.held).records
    assert len(held) == 1
    record = held[0]
    assert record.data["GATE_STATUS"] == "HELD"
    assert "CON-01" in record.data["GATE_CHECKS_FAILED"]
    assert "HOLDS" in record.data["GATE_REASONS"]
    assert "independent" in record.data["GATE_REASONS"]
    assert run.manifest.run_id in record.data["GATE_RUN"]


def test_a_held_record_can_be_fixed_and_fed_straight_back_in(tmp_path, mixed_input):
    """Which is why the held file is an SDF and not a report."""
    run = client.run(mixed_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="mixed")

    held = sdf.parse_file(outputs.held).records[0]
    repaired = sdf.strip_fields(held.raw, client.GATE_FIELDS)
    fixed_path = tmp_path / "fixed.sdf"
    fixed_path.write_bytes(repaired + sdf.RECORD_TERMINATOR)

    reparsed = sdf.parse_file(fixed_path).records[0]
    assert not [t for t in client.GATE_FIELDS if t in reparsed.data]
    assert reparsed.data["NAME"] == held.data["NAME"]


def test_rerun_workflow_a_waiver_clears_a_record_by_a_data_edit_alone(tmp_path):
    """No code changes: one row, with evidence, a decider and a date."""
    path = tmp_path / "waivable.sdf"
    path.write_bytes(
        sdf_bytes(salt_record(synonym="ethylamine hcl;CHEMBL1234", iupac="x"))
    )

    without = client.run(path, decisions=decisions.load(empty_decisions(tmp_path)))
    assert len(without.held) == 1

    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "SYN-01",
                "severity": "medium",
                "reason": "The catalogue identifier is a deliberate cross-reference here.",
                "evidence": "reviewed against the vendor page",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    with_waiver = client.run(path, decisions=decisions.load(directory))
    assert not with_waiver.held

    waived = [f for f in with_waiver.all_findings if f.check == "SYN-01"]
    assert len(waived) == 1
    assert waived[0].severity == checks.INFO
    assert "WAIVED" in waived[0].detail, "a waiver downgrades, it never hides"


def test_rerun_workflow_an_open_question_holds_a_record_the_checks_like(tmp_path):
    """ "Which salt does the lab hold?" is not answerable by looking at the file."""
    path = tmp_path / "clean.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="ethanamine hydrochloride")))

    directory = write_decisions(
        tmp_path,
        quarantine=[
            {
                "cas": CLEAN_CAS,
                "check_id": "CON-02",
                "reason": "Which salt form does the lab actually hold? Name and label disagree.",
                "opened_on": "2026-09-15",
                "blocked_on": "Lab answer naming the salt form.",
                "resolved_on": "",
            }
        ],
    )
    held_run = client.run(path, decisions=decisions.load(directory))
    assert len(held_run.held) == 1
    assert any(f.check == "QUAR-01" for f in held_run.held[0].blocking)


def test_rerun_workflow_a_resolution_date_releases_the_record(tmp_path):
    """The whole point: a record leaves quarantine by one edit to a data file."""
    path = tmp_path / "clean.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="ethanamine hydrochloride")))
    row = {
        "cas": CLEAN_CAS,
        "check_id": "CON-02",
        "reason": "Which salt form does the lab actually hold? Name and label disagree.",
        "opened_on": "2026-09-15",
        "blocked_on": "Lab answer naming the salt form.",
        "resolved_on": "",
    }
    held = client.run(
        path, decisions=decisions.load(write_decisions(tmp_path, quarantine=[row]))
    )
    assert len(held.held) == 1

    resolved = client.run(
        path,
        decisions=decisions.load(
            write_decisions(tmp_path, quarantine=[{**row, "resolved_on": "2026-10-01"}])
        ),
    )
    assert not resolved.held
    assert not [f for f in resolved.all_findings if f.check == "QUAR-01"]


def test_open_questions_are_written_with_the_compound_name(tmp_path):
    """The decisions data is keyed on CAS; nobody remembers which compound that is."""
    path = tmp_path / "clean.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="ethanamine hydrochloride")))
    directory = write_decisions(
        tmp_path,
        quarantine=[
            {
                "cas": CLEAN_CAS,
                "check_id": "CON-02",
                "reason": "Which salt form does the lab actually hold? Name and label disagree.",
                "opened_on": "2026-09-15",
                "blocked_on": "Lab answer naming the salt form.",
                "resolved_on": "",
            }
        ],
    )
    run = client.run(path, decisions=decisions.load(directory))
    outputs = gate_io.write(run, tmp_path / "out", stem="clean")

    rows = list(csv.DictReader(outputs.open_questions.open()))
    assert len(rows) == 1
    assert rows[0]["cas"] == CLEAN_CAS
    assert rows[0]["name"] == "ethylamine hydrochloride"
    assert "Lab answer" in rows[0]["detail"]


def test_a_waiver_on_a_low_finding_applies_instead_of_doing_nothing(tmp_path):
    """Only high and medium were waived, so a low waiver was a silent no-op.

    Five checks emit `low` as their only severity, and a waiver is a recorded
    judgement about a finding -- whether that finding would have held the record
    is a separate question, and not the one the waiver asked.
    """
    path = tmp_path / "machine_name.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="ethanamine;hydrochloride")))
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "IUP-01",
                "severity": "low",
                "reason": "The machine name is what the depositor supplied; kept.",
                "evidence": "reviewed against the submission sheet",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    run = client.run(path, decisions=decisions.load(directory), allow_medium=True)
    waived = [f for f in run.all_findings if f.check == "IUP-01" and f.waiver]
    assert waived, [f.check for f in run.all_findings]
    assert "WAIVED" in waived[0].detail


def test_a_waiver_that_matched_a_record_but_no_finding_is_reported(
    tmp_path, clean_input
):
    """The state `_check_ids` cannot catch: the check is real, it just never fires.

    A waiver that applies to nothing reads exactly like one that applies, and
    "matches no record" did not cover it -- the record is right there.
    """
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "CON-05",
                "severity": "high",
                "reason": "A decision about a finding this record does not raise.",
                "evidence": "somewhere",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    run = client.run(
        clean_input, decisions=decisions.load(directory), allow_medium=True
    )
    assert any(
        "waived no finding" in message and CLEAN_CAS in message
        for message in run.unmatched_decisions
    ), run.unmatched_decisions


def test_a_waiver_that_did_its_job_is_not_reported_as_inert(tmp_path):
    path = tmp_path / "clean.sdf"
    path.write_bytes(sdf_bytes(salt_record(synonym="CHEMBL25", iupac="x")))
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "SYN-01",
                "severity": "medium",
                "reason": "The catalogue identifier is a deliberate cross-reference.",
                "evidence": "reviewed against the vendor page",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    run = client.run(path, decisions=decisions.load(directory), allow_medium=True)
    assert not [m for m in run.unmatched_decisions if "SYN-01" in m]


def test_a_decision_matching_no_record_is_reported_on_the_run(tmp_path, clean_input):
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": "50-00-0",
                "check_id": "SYN-01",
                "severity": "medium",
                "reason": "A decision about a record that is not in this batch at all.",
                "evidence": "somewhere",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    run = client.run(
        clean_input, decisions=decisions.load(directory), allow_medium=True
    )
    assert run.unmatched_decisions
    assert "50-00-0" in run.unmatched_decisions[0]
    assert run.manifest.results["unmatched_decisions"] == run.unmatched_decisions


def test_a_decision_still_applies_when_the_cas_spelling_needs_repair(tmp_path):
    """The decision is about the compound; INT-04 is about the string.

    Joining decisions on the raw CAS_NO meant a record spelled ``0557-66-4``
    could not match its own waiver, and the run then reported "matches no record"
    while the record sat right there in the file. INT-04 still holds it for the
    spelling, which is the separate and correct complaint.
    """
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "SYN-01",
                "severity": "medium",
                "reason": "The catalogue identifier is a deliberate cross-reference here.",
                "evidence": "reviewed against the vendor page",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    for spelling in ("557-66-4", "0557-66-4", "557-66-4 00:00:00"):
        path = tmp_path / "spelled.sdf"
        path.write_bytes(
            sdf_bytes(
                salt_record(
                    cas=spelling, synonym="ethylamine hcl;CHEMBL1234", iupac="x"
                )
            )
        )
        run = client.run(path, decisions=decisions.load(directory))

        waived = [f for f in run.all_findings if f.check == "SYN-01"]
        assert len(waived) == 1, spelling
        assert waived[0].waiver, f"the waiver must apply to {spelling}"
        assert not run.unmatched_decisions, spelling

        malformed = [f for f in run.all_findings if f.check == "INT-04"]
        if spelling != "557-66-4":
            assert malformed, f"INT-04 must still report {spelling}"
            assert run.held


def test_the_findings_table_reports_the_cas_as_the_file_spells_it(tmp_path):
    """Joining normalises; reporting must not, or the fix instruction is wrong."""
    path = tmp_path / "spelled.sdf"
    path.write_bytes(sdf_bytes(salt_record(cas="0557-66-4", iupac="x")))
    run = client.run(path, allow_medium=True)
    assert all(f.cas == "0557-66-4" for f in run.all_findings if f.cas)


# --------------------------------------------------------------- independence


def test_the_summary_reports_independent_and_circular_separately(tmp_path, mixed_input):
    """Adding them together is the mistake the whole design exists to prevent."""
    run = client.run(mixed_input, allow_medium=True)
    text = gate_io.summary(run)
    assert "independent" in text
    assert "circular" in text
    counts = run.independence()
    assert set(counts) == {
        "independent",
        "circular",
        "independent_holding",
        "circular_holding",
    }
    assert counts["independent"] + counts["circular"] == len(run.all_findings)


def test_the_holding_counts_match_what_the_run_actually_did(tmp_path):
    """A summary claiming findings held a record while none was held is worse than none."""
    path = tmp_path / "medium.sdf"
    path.write_bytes(
        sdf_bytes(salt_record(synonym="ethylamine hcl;CHEMBL1234", iupac="x"))
    )

    lenient = client.run(path, allow_medium=True)
    assert not lenient.held
    counts = lenient.independence()
    assert counts["independent_holding"] == 0
    assert counts["circular_holding"] == 0

    strict = client.run(path)
    assert strict.held
    assert strict.independence()["independent_holding"] > 0


def test_with_no_reference_source_the_summary_says_clearance_is_weaker(clean_input):
    run = client.run(clean_input, allow_medium=True)
    text = gate_io.summary(run)
    assert "internally consistent, not confirmed" in text


def test_a_pubchem_only_finding_is_reported_as_circular(tmp_path):
    """The SDFs were generated from PubChem, so its agreement confirms nothing."""
    cache = tmp_path / "pubchem"
    cache.mkdir()
    (cache / f"{CLEAN_CAS}.json").write_text(
        json.dumps(
            {
                "cas": CLEAN_CAS,
                "cids_status": 200,
                "cids": [1],
                "cid": 1,
                "property_status": 200,
                "properties": {
                    "InChIKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "MolecularFormula": "C9H9NO",
                },
                "status": "ok",
            }
        )
    )
    path = tmp_path / "one.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="x")))

    run = client.run(
        path,
        evidence=external.Evidence(pubchem_dir=cache),
        allow_medium=True,
    )
    ext01 = [f for f in run.all_findings if f.check == "EXT-01"]
    assert len(ext01) == 1
    assert ext01[0].is_independent is False
    assert run.independence()["circular"] >= 1


# ------------------------------------------------------------------ manifest


def test_the_manifest_pins_the_input_by_hash(tmp_path, clean_input):
    run = client.run(clean_input, allow_medium=True)
    entry = run.manifest.inputs[0]
    assert entry["sha256"] == manifest_mod.sha256_file(clean_input)
    assert entry["records"] == 2
    assert entry["bytes"] == clean_input.stat().st_size


def test_the_manifest_records_the_environment_and_tool_version(clean_input):
    run = client.run(clean_input, allow_medium=True)
    environment = run.manifest.environment
    assert environment["tool"] == manifest_mod.TOOL
    assert environment["tool_version"] == manifest_mod.TOOL_VERSION
    assert environment["rdkit"]
    assert environment["python"]


def test_the_manifest_records_which_checks_ran_and_how(clean_input):
    run = client.run(clean_input, allow_medium=True, role=checks.ROLE_SALT)
    assert run.manifest.checks["role"] == checks.ROLE_SALT
    assert run.manifest.checks["medium_holds"] is False
    registered = run.manifest.checks["registered"]
    assert set(registered) == set(checks.CHECKS)
    assert registered["EXT-01"]["independent"] is False
    assert registered["EXT-04"]["independent"] is True


def test_the_manifest_says_which_sources_were_unavailable(clean_input):
    run = client.run(clean_input, allow_medium=True)
    reference = run.manifest.reference
    assert reference["cas_registry"] is None
    assert reference["chebi_release"] is None
    assert reference["pubchem_cache"] is None
    assert reference["sources_available"] == []


def test_the_manifest_pins_the_decisions_data(tmp_path, clean_input):
    directory = write_decisions(tmp_path)
    run = client.run(
        clean_input, decisions=decisions.load(directory), allow_medium=True
    )
    files = run.manifest.decisions["files"]
    assert len(files) == 2
    for entry in files:
        assert len(entry["sha256"]) == 64
    assert run.manifest.decisions["counts"]["waivers"] == 0


def test_the_run_id_changes_when_the_input_changes(tmp_path, clean_input, mixed_input):
    first = client.run(clean_input, allow_medium=True).manifest.run_id
    second = client.run(mixed_input, allow_medium=True).manifest.run_id
    assert first != second
    assert len(first) == 12


def test_the_run_id_changes_when_the_decisions_change(tmp_path, clean_input):
    plain = client.run(
        clean_input,
        decisions=decisions.load(empty_decisions(tmp_path)),
        allow_medium=True,
    )
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": CLEAN_CAS,
                "check_id": "SYN-01",
                "severity": "medium",
                "reason": "A decision recorded so the run identifier has to change.",
                "evidence": "somewhere",
                "decided_by": "A Curator",
                "decided_on": "2026-09-15",
            }
        ],
    )
    changed = client.run(
        clean_input, decisions=decisions.load(directory), allow_medium=True
    )
    assert plain.manifest.run_id != changed.manifest.run_id


def _mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache(directory, payload):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "557-66-4.json").write_text(json.dumps(payload))
    return directory


def test_the_manifest_pins_a_cache_by_its_contents_not_its_file_count(
    tmp_path, clean_input
):
    """Re-fetching a cache in place left the manifest byte for byte unchanged.

    The entry was a path and a `*.json` count, so two runs over completely
    different evidence claimed to be the same run -- including the same run_id,
    which is what the held records are stamped with.
    """
    cache = _cache(tmp_path / "pubchem", {"cid": 1})
    first = client.run(
        clean_input, evidence=external.Evidence(pubchem_dir=cache), allow_medium=True
    )
    _cache(cache, {"cid": 99999})
    second = client.run(
        clean_input, evidence=external.Evidence(pubchem_dir=cache), allow_medium=True
    )

    assert first.manifest.reference["pubchem_cache"]["records"] == 1
    assert (
        first.manifest.reference["pubchem_cache"]["sha256"]
        != second.manifest.reference["pubchem_cache"]["sha256"]
    )
    assert first.manifest.run_id != second.manifest.run_id


def test_the_same_cache_contents_pin_the_same_way(tmp_path, clean_input):
    here = _cache(tmp_path / "a" / "pubchem", {"cid": 1})
    there = _cache(tmp_path / "b" / "pubchem", {"cid": 1})
    runs = [
        client.run(
            clean_input, evidence=external.Evidence(pubchem_dir=d), allow_medium=True
        )
        for d in (here, there)
    ]
    assert (
        runs[0].manifest.reference["pubchem_cache"]["sha256"]
        == runs[1].manifest.reference["pubchem_cache"]["sha256"]
    )


def test_the_run_id_does_not_change_when_the_evidence_moves(tmp_path, clean_input):
    """The manifest claims a run is replayable; an absolute path is not portable."""
    here = _cache(tmp_path / "a" / "pubchem", {"cid": 1})
    there = _cache(tmp_path / "b" / "pubchem", {"cid": 1})
    first = client.run(
        clean_input, evidence=external.Evidence(pubchem_dir=here), allow_medium=True
    )
    second = client.run(
        clean_input, evidence=external.Evidence(pubchem_dir=there), allow_medium=True
    )

    assert first.manifest.run_id == second.manifest.run_id
    assert (
        first.manifest.reference["pubchem_cache"]["path"]
        != second.manifest.reference["pubchem_cache"]["path"]
    )


def test_the_run_id_does_not_change_when_the_decisions_data_moves(
    tmp_path, clean_input
):
    runs = [
        client.run(
            clean_input,
            decisions=decisions.load(empty_decisions(_mkdir(tmp_path / where))),
            allow_medium=True,
        )
        for where in ("first", "second")
    ]
    assert runs[0].manifest.run_id == runs[1].manifest.run_id
    assert (
        runs[0].manifest.decisions["files"][0]["path"]
        != runs[1].manifest.decisions["files"][0]["path"]
    )


def test_the_chebi_index_build_time_does_not_reach_the_run_id():
    """The docstring said no clock value reaches it; one did, embedded whole.

    `reference["chebi_release"]` is the index manifest verbatim, and that carries
    `generated` -- the wall-clock string from when the index was distilled.
    """
    ids = []
    for when in ("2026-08-01T00:00:00Z", "2026-09-15T12:00:00Z"):
        one = manifest_mod.Manifest()
        one.reference = {
            "chebi_release": {
                "generated": when,
                "release_dir": f"/somewhere/{when}",
                "index": {"cas": {"sha256": "abc"}},
            }
        }
        ids.append(one.compute_run_id())
    assert ids[0] == ids[1]


def test_the_run_id_still_changes_when_the_evidence_content_changes():
    ids = []
    for digest in ("abc", "def"):
        one = manifest_mod.Manifest()
        one.reference = {"chebi_release": {"index": {"cas": {"sha256": digest}}}}
        ids.append(one.compute_run_id())
    assert ids[0] != ids[1]


def test_git_commit_never_raises_outside_a_work_tree(tmp_path):
    assert (
        manifest_mod.git_commit(tmp_path) in ("unknown",)
        or len(manifest_mod.git_commit(tmp_path)) == 40
    )


def _repo(tmp_path: Path) -> Path:
    """A one-commit git work tree, for the manifest's environment block."""
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "code.py").write_text("SEVERITY = 'high'\n")
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "T",
        "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": "T",
        "GIT_COMMITTER_EMAIL": "t@example.invalid",
    }
    for args in (["init", "-q"], ["add", "-A"], ["commit", "-qm", "one"]):
        subprocess.run(["git", *args], cwd=repo, env=env, check=True)
    return repo


def test_the_manifest_says_whether_the_work_tree_matched_the_commit(tmp_path):
    """A commit alone does not say what ran.

    An uncommitted edit to a check is the normal state of a working session, and
    it left the manifest naming a commit whose code produces different findings
    from the ones recorded beside it, with nothing to mark the discrepancy.
    """
    repo = _repo(tmp_path)
    clean = manifest_mod.git_worktree_state(repo)
    assert clean["git_dirty"] is False
    assert clean["git_diff_sha256"] == ""

    (repo / "code.py").write_text("SEVERITY = 'low'\n")
    dirty = manifest_mod.git_worktree_state(repo)
    assert dirty["git_dirty"] is True
    assert len(dirty["git_diff_sha256"]) == 64


def test_the_same_uncommitted_edit_hashes_the_same_and_a_different_one_does_not(
    tmp_path,
):
    repo = _repo(tmp_path)
    (repo / "code.py").write_text("SEVERITY = 'low'\n")
    first = manifest_mod.git_worktree_state(repo)["git_diff_sha256"]
    (repo / "code.py").write_text("SEVERITY = 'low'\n")
    assert manifest_mod.git_worktree_state(repo)["git_diff_sha256"] == first
    (repo / "code.py").write_text("SEVERITY = 'info'\n")
    assert manifest_mod.git_worktree_state(repo)["git_diff_sha256"] != first


def test_an_untracked_file_does_not_make_the_work_tree_dirty(tmp_path):
    """Otherwise a gitignored run directory makes every manifest dirty."""
    repo = _repo(tmp_path)
    (repo / "chebi_run_2026_09a").mkdir()
    (repo / "chebi_run_2026_09a" / "batch.sdf").write_text("x\n")
    assert manifest_mod.git_worktree_state(repo)["git_dirty"] is False


def test_the_work_tree_state_is_unknown_rather_than_clean_outside_a_repo(tmp_path):
    """None and False have to read differently: one is a fact, one is a silence."""
    state = manifest_mod.git_worktree_state(tmp_path)
    assert state["git_dirty"] is None
    assert state["git_diff_sha256"] == "unknown"


def test_the_environment_block_carries_the_work_tree_state(tmp_path, clean_input):
    run = client.run(clean_input, allow_medium=True, repo=_repo(tmp_path))
    assert run.manifest.environment["git_dirty"] is False
    assert len(run.manifest.environment["git_commit"]) == 40


def test_two_inputs_gated_into_one_directory_keep_both_sets_of_tables(
    tmp_path, clean_input, mixed_input
):
    """The SDFs were stem-prefixed and the tables were not.

    So the second run replaced the first run's findings, open questions and
    manifest, silently, while leaving both SDF pairs in place -- a held record
    stamped with its own GATE_RUN sitting beside a manifest for a different run.
    The thing that explains the verdict was the thing that got overwritten.
    """
    out = tmp_path / "out"
    first = gate_io.write(
        client.run(clean_input, allow_medium=True), out, stem=clean_input.stem
    )
    before = first.findings.read_bytes()
    second = gate_io.write(
        client.run(mixed_input, allow_medium=True), out, stem=mixed_input.stem
    )

    assert first.findings != second.findings
    assert first.manifest != second.manifest
    assert first.open_questions != second.open_questions
    assert first.findings.read_bytes() == before
    for path in (first.findings, first.manifest, second.findings, second.manifest):
        assert path.exists()


def test_every_output_carries_the_input_stem(tmp_path, clean_input):
    outputs = gate_io.write(
        client.run(clean_input, allow_medium=True), tmp_path / "out", stem="batch7"
    )
    written = sorted(p.name for p in (tmp_path / "out").iterdir())
    assert written == [
        "batch7_cleared.sdf",
        "batch7_findings.csv",
        "batch7_held.sdf",
        "batch7_open_questions.csv",
        "batch7_run_manifest.json",
    ]
    assert outputs.manifest.name == "batch7_run_manifest.json"


def test_a_crlf_file_clears_every_record(tmp_path):
    """What CHEBI_GATE.md claims, asserted rather than assumed.

    INT-06 rates CRLF `low` so it holds nothing -- but every record was held
    anyway, on INT-02 "mol title differs from NAME", a defect none of them had.
    `allow_medium=True` in the older CLI test is why this stayed invisible.
    """
    path = tmp_path / "crlf.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="x")).replace(b"\n", b"\r\n"))
    run = client.run(path)

    assert not run.held, [(f.check, f.detail) for r in run.held for f in r.blocking]
    assert len(run.cleared) == 1
    assert [f.check for f in run.file_findings] == ["INT-06"]


def test_a_holding_file_finding_is_marked_as_holding_in_the_annotation(tmp_path):
    """`result.blocking` holds only the record's own findings.

    So a high INT-06 -- a non-ASCII byte, which holds every record through
    `file_blocks` -- was rendered without HOLDS and dropped out of
    GATE_CHECKS_FAILED whenever the record also had a blocking finding of its own.
    The annotated SDF is the artifact a chemist reads; it was the last output
    still guessing.
    """
    path = tmp_path / "nonascii.sdf"
    record = salt_record(name="cafe hydrochloride", iupac="x")
    path.write_bytes(record.replace("cafe", "caf\u00e9").encode("utf-8") + b"$$$$\n")
    run = client.run(path, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="nonascii")

    held = outputs.held.read_text(errors="replace")
    reasons = held.split("> <GATE_REASONS>")[1]
    assert "INT-06 [high] HOLDS" in reasons
    failed = held.split("> <GATE_CHECKS_FAILED>")[1].splitlines()[1]
    assert "INT-06" in failed


def test_a_file_finding_that_holds_nothing_is_not_written_as_holding(tmp_path):
    """The column was the constant "yes" for every whole-file finding.

    INT-06 reports CRLF line endings at low, which holds nothing, so the table
    said a finding held a record when the gate had cleared it.
    """
    path = tmp_path / "crlf.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="x")).replace(b"\n", b"\r\n"))
    run = client.run(path, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="crlf")

    rows = {r["check"]: r for r in csv.DictReader(outputs.findings.open())}
    assert rows["INT-06"]["severity"] == checks.LOW
    assert rows["INT-06"]["holds"] == "no"


def test_a_file_finding_that_does_hold_is_still_written_as_holding(tmp_path):
    """The narrowing must not turn the column off for a real whole-file block."""
    path = tmp_path / "nonascii.sdf"
    record = salt_record(name="cafe hydrochloride", iupac="x")
    path.write_bytes(record.replace("cafe", "caf\u00e9").encode("utf-8") + b"$$$$\n")
    run = client.run(path, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="nonascii")

    rows = [r for r in csv.DictReader(outputs.findings.open()) if r["record"] == "0"]
    assert rows, "expected a whole-file finding"
    assert all(r["holds"] == "yes" for r in rows if r["severity"] in ("high", "medium"))


# ------------------------------------------------------------- findings table


def test_the_findings_table_carries_every_finding_including_the_quiet_ones(
    tmp_path, mixed_input
):
    run = client.run(mixed_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="mixed")
    rows = list(csv.DictReader(outputs.findings.open()))
    assert len(rows) == len(run.all_findings)
    assert set(rows[0]) == set(gate_io.FINDINGS_COLUMNS)
    assert any(r["severity"] == "low" for r in rows)
    assert all(r["independent"] in ("yes", "no") for r in rows)
    assert all(r["holds"] in ("yes", "no") for r in rows)


def test_the_findings_table_is_keyed_on_cas(tmp_path, mixed_input):
    run = client.run(mixed_input, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="mixed")
    rows = list(csv.DictReader(outputs.findings.open()))
    assert all(r["cas"] for r in rows if r["record"] != "0")


# ------------------------------------------------------------------------ cli


def test_the_cli_exits_zero_when_everything_clears(tmp_path, clean_input, capsys):
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-decisions",
            "--allow-medium",
        ]
    )
    assert code == EXIT_OK
    assert "cleared" in capsys.readouterr().out


def test_the_cli_exits_one_when_a_record_is_held(tmp_path, mixed_input):
    code = main(
        [
            str(mixed_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-decisions",
            "--allow-medium",
        ]
    )
    assert code == EXIT_HELD


def test_the_cli_exits_two_on_a_missing_input(tmp_path, capsys):
    code = main([str(tmp_path / "nope.sdf")])
    assert code == EXIT_USAGE
    assert "no such file" in capsys.readouterr().err


def test_the_cli_exits_two_on_malformed_decisions(tmp_path, clean_input, capsys):
    directory = tmp_path / "bad"
    directory.mkdir()
    (directory / decisions.WAIVERS_FILE).write_text("cas,check_id\n64-17-5,SYN-01\n")
    (directory / decisions.QUARANTINE_FILE).write_text(
        ",".join(decisions.QUARANTINE_COLUMNS) + "\n"
    )
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--decisions",
            str(directory),
        ]
    )
    assert code == EXIT_USAGE
    assert "error:" in capsys.readouterr().err


@pytest.mark.parametrize("option", ["--pubchem-cache", "--cas-common-chemistry"])
def test_the_cli_exits_two_on_a_cache_path_that_is_not_there(
    tmp_path, clean_input, capsys, option
):
    """A typo lost a whole source and the run reported the loss as the truth.

    Both caches were wrapped in Path() and handed to Evidence, which tests
    is_dir() and reads a false as "not configured". So the gate ran without the
    source, cleared what it could and exited 0 -- printing "with no independent
    source, clearance means internally consistent, not confirmed" on a run where
    the caller had asked for one.
    """
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-decisions",
            "--allow-medium",
            option,
            str(tmp_path / "typo"),
        ]
    )
    assert code == EXIT_USAGE
    error = capsys.readouterr().err
    assert option in error
    assert "does not exist" in error


def test_the_cli_says_when_a_cache_path_is_a_file_rather_than_a_directory(
    tmp_path, clean_input, capsys
):
    path = tmp_path / "pubchem.json"
    path.write_text("{}")
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-decisions",
            "--allow-medium",
            "--pubchem-cache",
            str(path),
        ]
    )
    assert code == EXIT_USAGE
    assert "is not a directory" in capsys.readouterr().err


def test_the_cli_accepts_a_cache_directory_that_is_there(tmp_path, clean_input):
    cache = tmp_path / "pubchem"
    cache.mkdir()
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-decisions",
            "--allow-medium",
            "--pubchem-cache",
            str(cache),
        ]
    )
    assert code in (EXIT_OK, EXIT_HELD)


def test_the_cli_exits_two_on_a_bad_release_directory(tmp_path, capsys):
    """--distil ran outside the error handler, so a typo was a traceback.

    The documented contract, which the parser's own epilog repeats, is "2 on a
    usage or input error".
    """
    code = main(["--distil", str(tmp_path / "nope"), str(tmp_path / "index")])
    assert code == EXIT_USAGE
    assert "error:" in capsys.readouterr().err


def test_the_cli_exits_two_when_the_output_directory_cannot_be_written(
    tmp_path, clean_input, capsys
):
    """gate_io.write ran after the handler, so an unwritable --out-dir did too."""
    blocker = tmp_path / "blocked"
    blocker.write_text("a file where a directory would go")
    code = main(
        [
            str(clean_input),
            "--out-dir",
            str(blocker / "out"),
            "--no-decisions",
            "--allow-medium",
        ]
    )
    assert code == EXIT_USAGE
    assert "error:" in capsys.readouterr().err


def test_distilling_from_the_cli_records_when_the_index_was_built(tmp_path):
    """Only tests ever passed `generated`, so every real manifest recorded "".

    CHEBI_GATE.md describes that field as the reason index_manifest.json is not
    byte-identical between builds, which was not true of any index the CLI made.
    """
    release = tmp_path / "release"
    release.mkdir()
    for name, body in (
        ("structures.tsv", release_fixtures.STRUCTURES),
        ("compounds.tsv", release_fixtures.COMPOUNDS),
        ("database_accession.tsv", release_fixtures.ACCESSIONS),
        ("secondary_ids.tsv", release_fixtures.SECONDARY),
        ("status.tsv", release_fixtures.STATUS),
    ):
        (release / name).write_text(body)
    index_dir = tmp_path / "index"

    assert main(["--distil", str(release), str(index_dir)]) == EXIT_OK
    manifest = json.loads((index_dir / chebi_release.INDEX_MANIFEST).read_text())
    assert manifest["generated"].endswith("Z")
    assert manifest["generated"][:4].isdigit()


def test_the_cli_rejects_an_unknown_role(clean_input):
    with pytest.raises(SystemExit):
        main([str(clean_input), "--role", "nonsense"])


def test_the_gate_rejects_an_unknown_role_in_the_api(clean_input):
    with pytest.raises(client.GateError, match="unknown role"):
        client.run(clean_input, role="nonsense")


def test_an_empty_sdf_is_an_error_not_an_empty_success(tmp_path):
    path = tmp_path / "empty.sdf"
    path.write_bytes(b"")
    with pytest.raises(client.GateError, match="no records"):
        client.run(path)


def test_the_cli_can_distil_a_release(tmp_path, capsys):
    release = tmp_path / "release"
    release.mkdir()
    (release / "structures.tsv").write_text(
        "compound_id\tstandard_inchi_key\n16236\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N\n"
    )
    (release / "compounds.tsv").write_text("id\tascii_name\tstars\n16236\tethanol\t3\n")
    (release / "database_accession.tsv").write_text(
        "compound_id\taccession_number\ttype\n16236\t64-17-5\tCAS\n"
    )
    code = main(["unused.sdf", "--distil", str(release), str(tmp_path / "index")])
    assert code == EXIT_OK
    assert "indexed 1 structures" in capsys.readouterr().out
    assert chebi_release.load_index(tmp_path / "index").cas("64-17-5") == "CHEBI:16236"


def test_the_shipped_registry_and_index_defaults_are_absent_not_guessed():
    """No evidence is a valid run that reports reduced coverage, not a crash."""
    assert len(casreg.EMPTY) == 0
    assert len(chebi_release.EMPTY) == 0
    assert external.Evidence().available == ()


# --------------------------------------------------- malformed and non-ASCII input


@pytest.mark.parametrize(
    "label,build",
    [
        ("truncated", lambda r: sdf_bytes(r, terminate_last=False)),
        ("dollars with no newline", lambda r: sdf_bytes(r)[:-1]),
        ("trailing content", lambda r: sdf_bytes(r) + b"leftover\n"),
    ],
)
def test_the_gate_refuses_a_malformed_file_rather_than_working_around_it(
    tmp_path, label, build
):
    """If a record boundary is in doubt, so is every finding attributed to a record.

    This used to clear silently: exit 0, a header-only findings table, a manifest
    claiming 2 records and no mention of truncation, and a cleared file five bytes
    longer than its input.
    """
    path = tmp_path / "malformed.sdf"
    path.write_bytes(build(salt_record(iupac="ethanamine hydrochloride")))
    with pytest.raises(client.GateError, match="not a well-formed SDF"):
        client.run(path, allow_medium=True)


def test_the_cli_exits_two_on_a_malformed_file(tmp_path, capsys):
    path = tmp_path / "malformed.sdf"
    path.write_bytes(sdf_bytes(salt_record(iupac="x"), terminate_last=False))
    code = main([str(path), "--out-dir", str(tmp_path / "out"), "--no-decisions"])
    assert code == EXIT_USAGE
    assert "not a well-formed SDF" in capsys.readouterr().err


def test_a_non_ascii_byte_does_not_truncate_the_findings_table(tmp_path):
    """The CSV writers opened with strict UTF-8, so a detail echoing a non-ASCII
    byte raised part-way through writerows: a truncated findings.csv that still
    looked like valid CSV, no run manifest at all, and an exit status of 1 that is
    indistinguishable from "records were held".
    """
    raw = sdf_bytes(salt_record(iupac="x"), neutral_record(cas=OTHER_CAS, iupac="y"))
    path = tmp_path / "nonascii.sdf"
    path.write_bytes(raw.replace(b"ISA36807", b"\xe9SA36807"))

    run = client.run(path, allow_medium=True)
    outputs = gate_io.write(run, tmp_path / "out", stem="nonascii")

    assert outputs.manifest.exists(), "the manifest must still be written"
    rows = list(csv.DictReader(outputs.findings.open(encoding="utf-8")))
    assert len(rows) == len(run.all_findings)
    assert outputs.held.read_bytes(), "the non-ASCII record is held, not lost"
