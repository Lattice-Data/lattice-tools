"""Unit tests for the cumulative ChEBI answer ledger (no network)."""

from __future__ import annotations

from pathlib import Path

import pytest

from chebi_ledger import (
    FIELDS,
    LEDGER_PATH,
    UNVERIFIED,
    VERIFIED,
    answer_key,
    load,
    lookup,
    merge,
    save,
)


def _entry(name: str, cas: str, **over: str) -> dict[str, str]:
    row = {
        "experimental_perturbation": name,
        "cas": cas,
        "verdict": VERIFIED,
        "term_id": "CHEBI:1",
        "chebi_name": "thing",
        "reason_code": "",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "identity_evidence": "name and CAS resolved independently to the same structure",
    }
    row.update(over)
    return row


# --------------------------------------------------------------------------
# The key
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    ["PHA 767491", "PHA-767491", "PHA_767491", "pha767491", "PHA–767491"],
)
def test_spelling_variants_fold_to_one_key(name: str) -> None:
    """Three spellings of PHA-767491 really did appear in one batch."""
    assert answer_key(name, "845714-00-3") == answer_key("PHA 767491", "845714-00-3")


@pytest.mark.parametrize(
    "name", ["(±)-Blebbistatin", "(R)-Blebbistatin", "(?)-Blebbistatin"]
)
def test_stereo_and_mojibake_markers_do_not_fold_away(name: str) -> None:
    """The one folding this must never do.

    Collapsing `(±)-Blebbistatin` onto `Blebbistatin` would hand a racemate the
    answer found for an unspecified-stereo record, and letting a `?`-damaged name
    inherit its repaired twin's identifier would launder a guess into a verdict.
    Both must pay for their own lookup.
    """
    assert answer_key(name, "1") != answer_key("Blebbistatin", "1")


def test_the_same_name_with_two_registry_numbers_is_two_answers() -> None:
    """A free base and its salt share a name and are not the same compound."""
    assert answer_key("Doxorubicin hydrochloride", "25316-40-9") != answer_key(
        "Doxorubicin hydrochloride", "23214-92-8"
    )


# --------------------------------------------------------------------------
# Load, save, lookup
# --------------------------------------------------------------------------


def test_a_missing_ledger_reads_as_empty_not_an_error() -> None:
    assert load(Path("/nonexistent/ledger.tsv")) == {}


def test_the_ledger_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "ledger.tsv"
    ledger: dict = {}
    merge(ledger, [_entry("Aspirin", "50-78-2"), _entry("Ethanol", "64-17-5")], "1")
    assert save(ledger, path) == 2
    assert set(load(path)) == set(ledger)


def test_lookup_reports_a_full_match() -> None:
    ledger: dict = {}
    merge(ledger, [_entry("Aspirin", "50-78-2")], "1")
    entry, how = lookup(ledger, "aspirin", "50-78-2")
    assert how == "name+cas"
    assert entry is not None and entry["term_id"] == "CHEBI:1"


def test_lookup_distinguishes_a_half_match_from_an_answer() -> None:
    """A half match must be reported as one, never adopted.

    Same name with a different CAS is exactly the signal the method interrogates:
    nine such rows in the second batch produced two real findings, one a doxorubicin
    row that had swapped to the free base's registry number.
    """
    ledger: dict = {}
    merge(ledger, [_entry("Doxorubicin hydrochloride", "25316-40-9")], "1")
    _, how = lookup(ledger, "Doxorubicin hydrochloride", "23214-92-8")
    assert how == "name_only"
    _, how = lookup(ledger, "Something else entirely", "25316-40-9")
    assert how == "cas_only"


def test_lookup_reports_nothing_when_nothing_matches() -> None:
    entry, how = lookup({}, "Aspirin", "50-78-2")
    assert (entry, how) == (None, "")


# --------------------------------------------------------------------------
# Merging a later batch
# --------------------------------------------------------------------------


def test_a_repeat_keeps_its_first_seen_batch() -> None:
    ledger: dict = {}
    merge(ledger, [_entry("Aspirin", "50-78-2")], "1")
    merge(ledger, [_entry("Aspirin", "50-78-2")], "2")
    entry = ledger[answer_key("Aspirin", "50-78-2")]
    assert (entry["first_seen_batch"], entry["last_updated_batch"]) == ("1", "2")


def test_an_unchanged_answer_is_not_reported_as_a_change() -> None:
    ledger: dict = {}
    merge(ledger, [_entry("Aspirin", "50-78-2")], "1")
    assert merge(ledger, [_entry("Aspirin", "50-78-2")], "2") == []


def test_a_withdrawn_identifier_is_returned_not_swallowed() -> None:
    """The case this return value exists for.

    A stereochemistry correction withdrew eight identifiers that a previous batch had
    published. The later batch wins, but an answer moving is a finding about what was
    already deposited, so the caller has to be told.
    """
    ledger: dict = {}
    merge(ledger, [_entry("Formoterol hemifumarate", "43229-80-7")], "1")
    changed = merge(
        ledger,
        [
            _entry(
                "Formoterol hemifumarate",
                "43229-80-7",
                verdict=UNVERIFIED,
                term_id="",
                reason_code="cas_and_name_are_different_stereoisomers",
            )
        ],
        "2",
    )
    assert len(changed) == 1
    _, before, after = changed[0]
    assert (before["verdict"], before["term_id"]) == (VERIFIED, "CHEBI:1")
    assert (after["verdict"], after["term_id"]) == (UNVERIFIED, "")


# --------------------------------------------------------------------------
# The committed ledger itself
# --------------------------------------------------------------------------


def test_the_committed_ledger_is_internally_consistent() -> None:
    """Integrity of the real file, not a fixture.

    Every verified row must carry both an identifier and the structure it was matched
    on, or a later batch would inherit an answer it cannot re-check. Every unverified
    row must carry a reason code, or the dead end is undiagnosed and the next batch
    cannot tell "already asked" from "never looked".
    """
    ledger = load(LEDGER_PATH)
    assert ledger, "the committed ledger is missing or empty"
    verified = [r for r in ledger.values() if r["verdict"] == VERIFIED]
    unverified = [r for r in ledger.values() if r["verdict"] != VERIFIED]
    assert verified and unverified
    assert all(r["term_id"].startswith("CHEBI:") for r in verified)
    assert all(r["inchikey"] for r in verified)
    assert all(r["reason_code"] for r in unverified)
    assert all(not r["term_id"] for r in unverified)
    assert all(set(r) == set(FIELDS) for r in ledger.values())
