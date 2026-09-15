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
from pathlib import Path

import pytest

from chebi_gate import casreg, chebi_release, checks, client, decisions, external
from chebi_gate import io as gate_io
from chebi_gate import manifest as manifest_mod
from chebi_gate import sdf
from chebi_gate.cli import EXIT_HELD, EXIT_OK, EXIT_USAGE, main
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


def test_a_decision_matching_no_record_is_reported_on_the_run(tmp_path, clean_input):
    directory = write_decisions(
        tmp_path,
        waivers=[
            {
                "cas": "50-00-0",
                "check_id": "SYN-01",
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


def test_git_commit_never_raises_outside_a_work_tree(tmp_path):
    assert (
        manifest_mod.git_commit(tmp_path) in ("unknown",)
        or len(manifest_mod.git_commit(tmp_path)) == 40
    )


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
