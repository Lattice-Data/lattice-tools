"""Waivers and quarantine as input data: what must load, and what must fail the run.

The re-run workflow is the reason this module exists. A record is held, a question
goes out, an answer comes back weeks later, and the record has to clear without
anyone touching code. Every rule below exists so that a decision cannot quietly
stop applying: a malformed date, a check ID that no longer exists, a waiver with
nothing behind it, or a key that matches no record are all loud.
"""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path

import pytest

from chebi_gate import decisions

SHIPPED = Path(decisions.DECISIONS_DIR)

GOOD_WAIVER = {
    "cas": "64-17-5",
    "check_id": "CON-04",
    "reason": "Reviewed against the registry entry; the drawing is correct as it stands.",
    "evidence": "cas_registry.csv row for 64-17-5",
    "decided_by": "A Curator",
    "decided_on": "2026-09-15",
}
GOOD_QUARANTINE = {
    "cas": "64-17-5",
    "check_id": "CON-02",
    "reason": "Which salt form does the lab actually hold? The name and structure disagree.",
    "opened_on": "2026-09-15",
    "blocked_on": "Lab answer confirming the salt form.",
    "resolved_on": "",
}


def write(tmp_path: Path, *, waivers=(), quarantine=()) -> Path:
    """Build a decisions directory from row dicts."""
    for name, columns, rows in (
        (decisions.WAIVERS_FILE, decisions.WAIVER_COLUMNS, waivers),
        (decisions.QUARANTINE_FILE, decisions.QUARANTINE_COLUMNS, quarantine),
    ):
        with (tmp_path / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(columns))
            writer.writeheader()
            writer.writerows(rows)
    return tmp_path


def rejects(tmp_path, match, *, waivers=(), quarantine=()):
    write(tmp_path, waivers=waivers, quarantine=quarantine)
    with pytest.raises(decisions.DecisionsError, match=match):
        decisions.load(tmp_path)


# ------------------------------------------------------- the data that ships


def test_the_shipped_decisions_load():
    loaded = decisions.load()
    assert loaded.counts["waivers"] == 15
    assert loaded.counts["quarantine_rows"] == 7


def test_the_shipped_quarantine_is_the_seven_held_records():
    """All seven are still open; each names what would close it."""
    loaded = decisions.load()
    assert loaded.counts["quarantine_open"] == 7
    assert loaded.counts["quarantine_resolved"] == 0
    for rows in loaded.quarantine.values():
        for row in rows:
            assert row.blocked_on
            assert row.is_open


def test_every_shipped_waiver_names_evidence_and_a_decider():
    loaded = decisions.load()
    for waiver in loaded.waivers.values():
        assert waiver.evidence
        assert waiver.decided_by
        assert isinstance(waiver.decided_on, date)


def test_the_shipped_decisions_are_keyed_on_cas_not_on_name():
    """Handoff 4.1. One of the seven held records was itself renamed mid-run.

    ``(S)-(+)-a-Methylhistamine dihydrobromide`` became
    ``(S)-(+)-alpha-methylhistamine dihydrobromide``, which detached the record
    from its own open question under a name-keyed scheme.
    """
    loaded = decisions.load()
    for cas, _check in loaded.waivers:
        assert cas[0].isdigit()
    assert "75614-93-6" in loaded.quarantine


# --------------------------------------------------------------- happy paths


def test_a_waiver_is_found_by_cas_and_check(tmp_path):
    loaded = decisions.load(write(tmp_path, waivers=[GOOD_WAIVER]))
    waiver = loaded.waiver_for("64-17-5", "CON-04")
    assert waiver is not None
    assert waiver.decided_on == date(2026, 9, 15)
    assert loaded.waiver_for("64-17-5", "CON-02") is None
    assert loaded.waiver_for("557-66-4", "CON-04") is None


def test_one_record_may_carry_waivers_for_several_checks(tmp_path):
    second = {**GOOD_WAIVER, "check_id": "EXT-01"}
    loaded = decisions.load(write(tmp_path, waivers=[GOOD_WAIVER, second]))
    assert loaded.waiver_for("64-17-5", "CON-04") is not None
    assert loaded.waiver_for("64-17-5", "EXT-01") is not None


def test_a_quarantine_row_may_name_several_checks(tmp_path):
    row = {**GOOD_QUARANTINE, "check_id": "CON-02/EXT-04"}
    loaded = decisions.load(write(tmp_path, quarantine=[row]))
    (question,) = loaded.open_questions("64-17-5")
    assert question.check_ids == ("CON-02", "EXT-04")


def test_a_resolved_question_stops_being_open(tmp_path):
    """This is the whole re-run workflow: one data edit and the record moves on."""
    resolved = {**GOOD_QUARANTINE, "resolved_on": "2026-10-01"}
    loaded = decisions.load(write(tmp_path, quarantine=[resolved]))
    assert loaded.open_questions("64-17-5") == ()
    (done,) = loaded.resolved_questions("64-17-5")
    assert done.resolved_on == date(2026, 10, 1)
    assert done.is_open is False


def test_counts_separate_open_from_resolved(tmp_path):
    rows = [
        GOOD_QUARANTINE,
        {**GOOD_QUARANTINE, "cas": "557-66-4", "resolved_on": "2026-10-01"},
    ]
    loaded = decisions.load(write(tmp_path, quarantine=rows))
    assert loaded.counts == {
        "waivers": 0,
        "quarantine_rows": 2,
        "quarantine_open": 1,
        "quarantine_resolved": 1,
    }


# -------------------------------------------------------- stale decisions


def test_a_decision_matching_no_record_is_reported(tmp_path):
    """CAS beats NAME as a key but is not permanent: sometimes the CAS is the defect.

    For one of the seven held records the CAS itself is wrong, so correcting it
    moves the key. A decision that silently stops applying is worse than a missing
    one, because it still reads as applying.
    """
    loaded = decisions.load(
        write(tmp_path, waivers=[GOOD_WAIVER], quarantine=[GOOD_QUARANTINE])
    )
    assert loaded.unmatched({"64-17-5"}) == []

    stale = loaded.unmatched({"557-66-4"})
    assert len(stale) == 2
    assert any("waiver 64-17-5 CON-04" in s for s in stale)
    assert any("quarantine 64-17-5 (open)" in s for s in stale)


def test_a_stale_resolved_quarantine_says_it_is_resolved(tmp_path):
    resolved = {**GOOD_QUARANTINE, "resolved_on": "2026-10-01"}
    loaded = decisions.load(write(tmp_path, quarantine=[resolved]))
    (message,) = loaded.unmatched(set())
    assert "(resolved)" in message


# --------------------------------------------------------- what fails the run


def test_a_waiver_with_no_evidence_fails_the_run(tmp_path):
    """The handoff is explicit. A waiver with nothing behind it is a suppression."""
    rejects(
        tmp_path,
        "evidence is empty",
        waivers=[{**GOOD_WAIVER, "evidence": ""}],
    )


@pytest.mark.parametrize(
    "field", ["cas", "check_id", "reason", "decided_by", "decided_on"]
)
def test_a_waiver_missing_any_required_field_fails_the_run(tmp_path, field):
    rejects(tmp_path, f"{field} is empty", waivers=[{**GOOD_WAIVER, field: ""}])


@pytest.mark.parametrize(
    "field", ["cas", "check_id", "reason", "opened_on", "blocked_on"]
)
def test_a_quarantine_row_missing_any_required_field_fails_the_run(tmp_path, field):
    rejects(tmp_path, f"{field} is empty", quarantine=[{**GOOD_QUARANTINE, field: ""}])


def test_a_quarantine_row_with_no_route_out_fails_the_run(tmp_path):
    """A held record with no stated unblock condition is where the seven sat."""
    rejects(
        tmp_path,
        "blocked_on is empty",
        quarantine=[{**GOOD_QUARANTINE, "blocked_on": ""}],
    )


def test_a_gestured_reason_fails_the_run(tmp_path):
    rejects(tmp_path, "too short", waivers=[{**GOOD_WAIVER, "reason": "ok"}])


def test_an_unknown_check_id_fails_the_run(tmp_path):
    """A decision naming a check that does not exist applies to nothing, silently.

    Usually it means a check was renamed and the decisions data was not.
    """
    rejects(
        tmp_path,
        "unknown check id",
        waivers=[{**GOOD_WAIVER, "check_id": "CON-99"}],
    )


def test_a_waiver_covering_several_checks_fails_the_run(tmp_path):
    """One row per check, so each decision has its own reason and evidence."""
    rejects(
        tmp_path,
        "exactly one check",
        waivers=[{**GOOD_WAIVER, "check_id": "CON-04 EXT-01"}],
    )


@pytest.mark.parametrize("cas", ["not-a-cas", "64-17-4", "42971-09-05", "064-17-5"])
def test_a_key_that_could_never_match_a_record_fails_the_run(tmp_path, cas):
    """Including values that only verify after repair: the record holds the raw string."""
    rejects(tmp_path, "well-formed CAS", waivers=[{**GOOD_WAIVER, "cas": cas}])


@pytest.mark.parametrize(
    "field,value",
    [("decided_on", "15/09/2026"), ("decided_on", "2026-13-01")],
)
def test_a_non_iso_date_fails_the_run(tmp_path, field, value):
    rejects(tmp_path, "not an ISO date", waivers=[{**GOOD_WAIVER, field: value}])


def test_a_resolution_before_the_question_was_opened_fails_the_run(tmp_path):
    rejects(
        tmp_path,
        "precedes",
        quarantine=[{**GOOD_QUARANTINE, "resolved_on": "2026-09-01"}],
    )


def test_a_duplicate_waiver_fails_the_run_and_names_the_first_row(tmp_path):
    rejects(tmp_path, "duplicate waiver", waivers=[GOOD_WAIVER, GOOD_WAIVER])


def test_unexpected_columns_fail_the_run(tmp_path):
    write(tmp_path)
    path = tmp_path / decisions.WAIVERS_FILE
    path.write_text("cas,check_id,reason\n64-17-5,CON-04,because\n")
    with pytest.raises(decisions.DecisionsError, match="columns are"):
        decisions.load(tmp_path)


def test_a_missing_decisions_file_fails_the_run(tmp_path):
    with pytest.raises(decisions.DecisionsError, match="not found"):
        decisions.load(tmp_path)


def test_the_error_names_the_row_number_an_editor_would_show(tmp_path):
    """Row 2 is the first data row, because row 1 is the header."""
    write(tmp_path, waivers=[{**GOOD_WAIVER, "evidence": ""}])
    with pytest.raises(decisions.DecisionsError, match="row 2"):
        decisions.load(tmp_path)
