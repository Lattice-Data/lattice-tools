"""Unit tests for chebi_terms client and io (mocked HTTP)."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import requests

from chebi_terms.client import (
    CAS_CONFIRMED,
    CAS_NOT_IN_CHEBI,
    NAME_MATCH,
    NAME_MISMATCH,
    NAME_SYNONYM_MATCH,
    NOT_CHECKED,
    OUTPUT_FIELDS_APPENDED,
    STATUS_INVALID,
    STATUS_LOOKUP_FAILED,
    STATUS_MISSING,
    STATUS_NOT_FOUND,
    STATUS_NOT_RELEASED,
    MAX_RETRIES,
    RETRY_BACKOFF,
    STATUS_OK,
    STATUS_SECONDARY,
    ChebiUnavailableError,
    clean_name,
    describe_chebi_id,
    extract_cas_numbers,
    extract_synonyms,
    fetch_compound,
    get_with_retry,
    match_key,
    normalize_chebi_id,
    verify_chebi_id,
    verify_payload,
)
from chebi_terms.io import (
    MAX_CONSECUTIVE_FAILURES,
    SUMMARY_STATUSES,
    ChebiTermsError,
    build_single_row,
    emit_single_chebi_id,
    verify_chebi_file,
)
from tests.chebi_terms_helpers import (
    ETHANOL,
    EPICATECHIN_SECONDARY,
    MARKUP_NAME,
    load_payload,
    mock_response,
    payload_copy,
    route_chebi_get,
)


# --------------------------------------------------------------------------
# identifier parsing
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    ["16236", "CHEBI:16236", "chebi:16236", " CHEBI:16236 ", "CHEBI:016236"],
)
def test_normalize_chebi_id_accepts_common_forms(raw: str) -> None:
    assert normalize_chebi_id(raw) == ("16236", "CHEBI:16236")


@pytest.mark.parametrize(
    "raw",
    ["", "   ", "abc", "CHEBI:", "CHEBI:abc", "16236a", "CHEBI:16236:1", None, "0"],
)
def test_normalize_chebi_id_rejects_junk(raw: str | None) -> None:
    assert normalize_chebi_id(raw) is None


# --------------------------------------------------------------------------
# name cleaning and comparison
# --------------------------------------------------------------------------


def test_clean_name_strips_nested_markup() -> None:
    raw = "[CH<small><sub>2</sub></small>Me(OH)]"
    assert clean_name(raw) == "[CH2Me(OH)]"


def test_clean_name_strips_stereo_markup_and_collapses_whitespace() -> None:
    raw = "(3<i>R</i>,4<i>S</i>)-diol   1,1-dioxide"
    assert clean_name(raw) == "(3R,4S)-diol 1,1-dioxide"


def test_clean_name_unescapes_entities() -> None:
    assert clean_name("alpha&ndash;beta") == "alpha–beta"


def test_clean_name_handles_none() -> None:
    assert clean_name(None) == ""


def test_match_key_normalizes_unicode_minus_to_hyphen() -> None:
    # ChEBI stores U+2212; curators type an ASCII hyphen.
    assert match_key("(−)-epicatechin") == match_key("(-)-epicatechin")


def test_match_key_is_case_insensitive() -> None:
    assert match_key("ETHANOL") == match_key("ethanol")


# --------------------------------------------------------------------------
# payload extraction
# --------------------------------------------------------------------------


def test_extract_synonyms_excludes_brand_names() -> None:
    synonyms = extract_synonyms(load_payload(ETHANOL))
    # 'Dehydrated ethanol' is a BRAND NAME in ChEBI.
    assert "Dehydrated ethanol" not in synonyms


def test_extract_synonyms_excludes_non_english() -> None:
    synonyms = extract_synonyms(load_payload(ETHANOL))
    for foreign in ("alcohol etílico", "Alkohol", "spiritus vini", "éthanol"):
        assert foreign not in synonyms


def test_extract_synonyms_omits_official_name_and_case_duplicates() -> None:
    synonyms = extract_synonyms(load_payload(ETHANOL))
    assert [s for s in synonyms if s.lower() == "ethanol"] == []


def test_extract_synonyms_strips_markup() -> None:
    synonyms = extract_synonyms(load_payload(ETHANOL))
    assert "C2H5OH" in synonyms
    assert not any("<" in s for s in synonyms)


def test_extract_synonyms_respects_max() -> None:
    assert len(extract_synonyms(load_payload(ETHANOL), max_synonyms=3)) == 3


def test_extract_synonyms_handles_missing_names_key() -> None:
    assert extract_synonyms({"name": "x"}) == []


def test_extract_synonyms_tolerates_malformed_names_block() -> None:
    assert extract_synonyms({"name": "x", "names": {"SYNONYM": "not-a-list"}}) == []


def test_extract_cas_numbers_finds_chebi_own_xrefs() -> None:
    assert "64-17-5" in extract_cas_numbers(load_payload(ETHANOL))


def test_extract_cas_numbers_handles_missing_accessions() -> None:
    assert extract_cas_numbers({"name": "x"}) == []


# --------------------------------------------------------------------------
# verify_payload — id_status
# --------------------------------------------------------------------------


def test_verify_payload_ok() -> None:
    result = verify_payload("CHEBI:16236", load_payload(ETHANOL))
    assert result["id_status"] == STATUS_OK
    assert result["chebi_accession"] == "CHEBI:16236"
    assert result["chebi_name"] == "ethanol"
    assert result["chebi_stars"] == 3
    assert "Ethyl alcohol" in result["chebi_synonyms"]
    assert result["chebi_definition"].startswith("A primary alcohol")


def test_verify_payload_detects_secondary_id_and_reports_primary() -> None:
    result = verify_payload("CHEBI:18484", load_payload(EPICATECHIN_SECONDARY))
    assert result["id_status"] == STATUS_SECONDARY
    assert result["chebi_accession"] == "CHEBI:90"


def test_verify_payload_not_found() -> None:
    result = verify_payload("CHEBI:99999999", None)
    assert result["id_status"] == STATUS_NOT_FOUND
    assert result["chebi_name"] == ""
    assert result["name_verdict"] == NOT_CHECKED
    assert result["cas_verdict"] == NOT_CHECKED


def test_verify_payload_not_released_takes_precedence_over_secondary() -> None:
    payload = payload_copy(EPICATECHIN_SECONDARY)
    payload["is_released"] = False
    result = verify_payload("CHEBI:18484", payload)
    assert result["id_status"] == STATUS_NOT_RELEASED
    # Facts are still reported so the row remains useful.
    assert result["chebi_name"] == "(−)-epicatechin"


def test_verify_payload_strips_markup_from_name() -> None:
    result = verify_payload("CHEBI:44567", load_payload(MARKUP_NAME))
    assert "<" not in result["chebi_name"]
    assert result["chebi_name"].startswith("(3R,4S,5S,6R)-")


def test_verify_payload_blank_stars_when_absent() -> None:
    payload = payload_copy(ETHANOL)
    payload["stars"] = None
    assert verify_payload("CHEBI:16236", payload)["chebi_stars"] == ""


# --------------------------------------------------------------------------
# verify_payload — name and CAS verdicts
# --------------------------------------------------------------------------


def test_name_verdict_exact_match() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_name="Ethanol"
    )
    assert result["name_verdict"] == NAME_MATCH


def test_name_verdict_matches_across_unicode_dash() -> None:
    result = verify_payload(
        "CHEBI:18484",
        load_payload(EPICATECHIN_SECONDARY),
        expected_name="(-)-epicatechin",
    )
    assert result["name_verdict"] == NAME_MATCH


def test_name_verdict_synonym_match() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_name="ethyl alcohol"
    )
    assert result["name_verdict"] == NAME_SYNONYM_MATCH


def test_name_verdict_synonym_match_accepts_brand_name() -> None:
    # Excluded from output, still accepted as evidence the ID is right.
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_name="Dehydrated ethanol"
    )
    assert result["name_verdict"] == NAME_SYNONYM_MATCH


def test_name_verdict_synonym_match_accepts_foreign_spelling() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_name="alcohol etilico"
    )
    assert result["name_verdict"] == NAME_SYNONYM_MATCH


def test_name_verdict_mismatch_flags_wrong_compound() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_name="caffeine"
    )
    assert result["name_verdict"] == NAME_MISMATCH


def test_name_verdict_not_checked_when_expected_blank() -> None:
    result = verify_payload("CHEBI:16236", load_payload(ETHANOL), expected_name="   ")
    assert result["name_verdict"] == NOT_CHECKED


def test_cas_verdict_confirmed() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_cas=" 64-17-5 "
    )
    assert result["cas_verdict"] == CAS_CONFIRMED


def test_cas_verdict_not_in_chebi() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_cas="58-08-2"
    )
    assert result["cas_verdict"] == CAS_NOT_IN_CHEBI


def test_cas_verdict_not_checked_when_expected_blank() -> None:
    result = verify_payload("CHEBI:16236", load_payload(ETHANOL), expected_cas="")
    assert result["cas_verdict"] == NOT_CHECKED


def test_cas_verdict_confirmed_across_unicode_hyphen() -> None:
    # U+2011 non-breaking hyphen, the same class of paste artefact _DASH_MAP
    # absorbs for names.
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_cas="64‑17‑5"
    )
    assert result["cas_verdict"] == CAS_CONFIRMED


def test_cas_verdict_confirmed_ignoring_internal_whitespace() -> None:
    result = verify_payload(
        "CHEBI:16236", load_payload(ETHANOL), expected_cas="64 - 17 - 5"
    )
    assert result["cas_verdict"] == CAS_CONFIRMED


def test_cas_verdict_does_not_confirm_dash_stripped_variant() -> None:
    # 64175 must NOT match 64-17-5: a false 'confirmed' is the dangerous direction.
    result = verify_payload("CHEBI:16236", load_payload(ETHANOL), expected_cas="64175")
    assert result["cas_verdict"] == CAS_NOT_IN_CHEBI


# --------------------------------------------------------------------------
# an unreachable ChEBI must never read as "no such compound"
# --------------------------------------------------------------------------


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_get_with_retry_returns_none_on_404(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.return_value = mock_response(404)
    assert get_with_retry("https://example.invalid/x") is None
    # A 404 is an answer, so it is not retried.
    assert mock_get.call_count == 1


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_get_with_retry_raises_on_connection_error(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.side_effect = requests.exceptions.ConnectionError("dns")
    with pytest.raises(ChebiUnavailableError):
        get_with_retry("https://example.invalid/x")
    assert mock_get.call_count == MAX_RETRIES


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_get_with_retry_raises_on_persistent_500(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.return_value = mock_response(500)
    with pytest.raises(ChebiUnavailableError):
        get_with_retry("https://example.invalid/x")


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_get_with_retry_does_not_sleep_after_the_final_attempt(
    mock_get: MagicMock, mock_sleep: MagicMock
) -> None:
    """Sleeping after the last try only delays the failure."""
    mock_get.side_effect = requests.exceptions.ConnectionError("down")
    with pytest.raises(ChebiUnavailableError):
        get_with_retry("https://example.invalid/x")

    assert mock_get.call_count == MAX_RETRIES
    assert mock_sleep.call_count == MAX_RETRIES - 1
    slept = [call.args[0] for call in mock_sleep.call_args_list]
    assert slept == [RETRY_BACKOFF, RETRY_BACKOFF * 2]
    # The figure MAX_CONSECUTIVE_FAILURES is justified by.
    assert sum(slept) == RETRY_BACKOFF * (2 ** (MAX_RETRIES - 1) - 1)


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_fetch_compound_returns_none_only_for_404(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.return_value = mock_response(404)
    assert fetch_compound("99999999") is None


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_fetch_compound_raises_on_non_dict_payload(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = ["not", "a", "dict"]
    mock_get.return_value = resp
    with pytest.raises(ChebiUnavailableError):
        fetch_compound("16236")


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_id_reports_lookup_failed_not_not_found(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.side_effect = requests.exceptions.Timeout("slow")
    result = verify_chebi_id(
        "CHEBI:16236", expected_name="ethanol", expected_cas="64-17-5"
    )
    assert result["id_status"] == STATUS_LOOKUP_FAILED
    assert result["id_status"] != STATUS_NOT_FOUND
    assert result["chebi_name"] == ""
    # Nothing was checked, so nothing may be asserted about the name or CAS.
    assert result["name_verdict"] == NOT_CHECKED
    assert result["cas_verdict"] == NOT_CHECKED


@pytest.mark.parametrize("blank", ["", "   "])
def test_verify_chebi_id_blank_is_missing(blank: str) -> None:
    assert verify_chebi_id(blank)["id_status"] == STATUS_MISSING


def test_verify_chebi_id_junk_is_invalid_not_missing() -> None:
    assert verify_chebi_id("not-an-id")["id_status"] == STATUS_INVALID


# --------------------------------------------------------------------------
# verify_chebi_id / describe_chebi_id
# --------------------------------------------------------------------------


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_id_end_to_end(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.side_effect = route_chebi_get
    result = verify_chebi_id(
        "CHEBI:16236", expected_name="ethanol", expected_cas="64-17-5"
    )
    assert result["id_status"] == STATUS_OK
    assert result["name_verdict"] == NAME_MATCH
    assert result["cas_verdict"] == CAS_CONFIRMED


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_id_invalid_makes_no_request(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    result = verify_chebi_id("not-an-id")
    assert result["id_status"] == STATUS_INVALID
    mock_get.assert_not_called()


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_id_unknown_id(mock_get: MagicMock, _sleep: MagicMock) -> None:
    mock_get.return_value = mock_response(404)
    assert verify_chebi_id("CHEBI:99999999")["id_status"] == STATUS_NOT_FOUND


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_id_non_json_response(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    resp = MagicMock()
    resp.status_code = 200
    resp.json.side_effect = ValueError("not json")
    mock_get.return_value = resp
    # A 200 with a garbled body is a server fault, not evidence about the compound.
    assert verify_chebi_id("CHEBI:16236")["id_status"] == STATUS_LOOKUP_FAILED


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_describe_chebi_id_leaves_verdicts_unchecked(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    mock_get.side_effect = route_chebi_get
    result = describe_chebi_id("16236")
    assert result["chebi_name"] == "ethanol"
    assert result["name_verdict"] == NOT_CHECKED
    assert result["cas_verdict"] == NOT_CHECKED


# --------------------------------------------------------------------------
# io — single ID
# --------------------------------------------------------------------------


def test_build_single_row_fills_all_fields() -> None:
    row = build_single_row("CHEBI:16236", {"chebi_name": "ethanol"})
    assert row["chebi_id"] == "CHEBI:16236"
    assert set(row) == {"chebi_id", *OUTPUT_FIELDS_APPENDED}
    assert row["chebi_accession"] == ""


def test_build_single_row_coerces_none() -> None:
    row = build_single_row("CHEBI:16236", {"chebi_name": None})
    assert row["chebi_name"] == ""


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_emit_single_json_to_file(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    out = tmp_path / "result.json"
    emit_single_chebi_id("CHEBI:16236", out)
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["chebi_name"] == "ethanol"
    assert payload["id_status"] == STATUS_OK


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_emit_single_csv_to_file(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    out = tmp_path / "result.csv"
    emit_single_chebi_id("CHEBI:16236", out, fmt="csv")
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert len(rows) == 1
    assert rows[0]["chebi_name"] == "ethanol"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_emit_single_json_to_stdout(
    mock_get: MagicMock, _sleep: MagicMock, capsys: pytest.CaptureFixture
) -> None:
    mock_get.side_effect = route_chebi_get
    emit_single_chebi_id("CHEBI:16236", None)
    payload = json.loads(capsys.readouterr().out)
    assert payload["chebi_accession"] == "CHEBI:16236"


# --------------------------------------------------------------------------
# io — batch
# --------------------------------------------------------------------------


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_happy_path(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id", "note"], [["CHEBI:16236", "keep me"]])
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out)

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["note"] == "keep me"
    assert rows[0]["chebi_name"] == "ethanol"
    assert rows[0]["id_status"] == STATUS_OK
    assert rows[0]["name_verdict"] == NOT_CHECKED


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_checks_name_and_cas_columns(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["chebi_id", "name", "CAS"],
        [
            ["CHEBI:16236", "ethanol", "64-17-5"],
            ["CHEBI:16236", "caffeine", "58-08-2"],
        ],
    )
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out, name_column="name", cas_column="CAS")

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["name_verdict"] == NAME_MATCH
    assert rows[0]["cas_verdict"] == CAS_CONFIRMED
    assert rows[1]["name_verdict"] == NAME_MISMATCH
    assert rows[1]["cas_verdict"] == CAS_NOT_IN_CHEBI


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_fetches_each_unique_id_once(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["chebi_id", "name"],
        [
            ["CHEBI:16236", "ethanol"],
            ["16236", "ethyl alcohol"],
            ["CHEBI:16236", "caffeine"],
        ],
    )
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out, name_column="name")

    assert mock_get.call_count == 1
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    # Same cached payload, three different per-row verdicts.
    assert [r["name_verdict"] for r in rows] == [
        NAME_MATCH,
        NAME_SYNONYM_MATCH,
        NAME_MISMATCH,
    ]


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_blank_and_invalid_ids(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id"], [[""], ["not-an-id"]])
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out)

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["id_status"] == STATUS_MISSING
    assert rows[1]["id_status"] == STATUS_INVALID
    # Blank rows still carry documented verdicts, not empty strings.
    for row in rows:
        assert row["name_verdict"] == NOT_CHECKED
        assert row["cas_verdict"] == NOT_CHECKED
    mock_get.assert_not_called()


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_overwrites_colliding_input_column(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id", "chebi_name"], [["CHEBI:16236", "stale value"]])
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out)

    with out.open(encoding="utf-8") as fh:
        header = next(csv.reader(fh))
    assert header.count("chebi_name") == 1
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["chebi_name"] == "ethanol"


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_does_not_cache_a_transient_failure(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    """One blip must not be replayed as a verdict for every later row."""
    calls = {"n": 0}

    def flaky(url: str, *args, **kwargs):
        calls["n"] += 1
        # Exhaust the first row's retries, then recover.
        if calls["n"] <= MAX_RETRIES:
            raise requests.exceptions.ConnectionError("blip")
        return route_chebi_get(url)

    mock_get.side_effect = flaky
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id"], [["CHEBI:16236"], ["CHEBI:16236"]])
    out = tmp_path / "out.csv"

    counts = verify_chebi_file(src, "chebi_id", out)

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["id_status"] == STATUS_LOOKUP_FAILED
    # Retried rather than served a cached failure.
    assert rows[1]["id_status"] == STATUS_OK
    assert rows[1]["chebi_name"] == "ethanol"
    assert counts[STATUS_LOOKUP_FAILED] == 1


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_aborts_after_consecutive_failures(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = requests.exceptions.ConnectionError("down")
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["chebi_id"],
        [[f"CHEBI:{16236 + i}"] for i in range(MAX_CONSECUTIVE_FAILURES + 3)],
    )
    out = tmp_path / "out.csv"

    with pytest.raises(ChebiTermsError, match="consecutive"):
        verify_chebi_file(src, "chebi_id", out)

    # Partial output survives, and stops at the abort point rather than
    # grinding through the rest of the sheet.
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert len(rows) == MAX_CONSECUTIVE_FAILURES
    assert all(row["id_status"] == STATUS_LOOKUP_FAILED for row in rows)


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_status_counts_sum_to_row_total(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    """The per-status summary lines must add up to the reported row total."""
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["chebi_id"],
        [["CHEBI:16236"], ["CHEBI:18484"], [""], ["not-an-id"], ["CHEBI:99999999"]],
    )
    out = tmp_path / "out.csv"

    counts = verify_chebi_file(src, "chebi_id", out)

    assert sum(counts.values()) == 5
    # Every status produced is one the summary actually prints.
    assert sum(counts[status] for status in SUMMARY_STATUSES) == 5
    assert counts[STATUS_MISSING] == 1
    assert counts[STATUS_NOT_FOUND] == 1


def test_verify_chebi_file_missing_input(tmp_path: Path) -> None:
    with pytest.raises(ChebiTermsError, match="Input file not found"):
        verify_chebi_file(tmp_path / "nope.csv", "chebi_id", tmp_path / "out.csv")


def test_verify_chebi_file_missing_chebi_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["other"], [["x"]])
    with pytest.raises(ChebiTermsError, match="Column 'chebi_id' not found"):
        verify_chebi_file(src, "chebi_id", tmp_path / "out.csv")


def test_verify_chebi_file_rejects_reserved_chebi_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_accession"], [["CHEBI:16236"]])
    with pytest.raises(ChebiTermsError, match="collides"):
        verify_chebi_file(src, "chebi_accession", tmp_path / "out.csv")


@pytest.mark.parametrize("kwarg", ["name_column", "cas_column"])
def test_verify_chebi_file_rejects_reserved_check_column(
    tmp_path: Path, kwarg: str
) -> None:
    """Checking against a column we overwrite would erase the evidence."""
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id", "chebi_name"], [["CHEBI:16236", "ethanol"]])
    with pytest.raises(ChebiTermsError, match="erase the value being checked"):
        verify_chebi_file(
            src, "chebi_id", tmp_path / "out.csv", **{kwarg: "chebi_name"}
        )


def test_verify_chebi_file_missing_name_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id"], [["CHEBI:16236"]])
    with pytest.raises(ChebiTermsError, match="--name-column"):
        verify_chebi_file(
            src, "chebi_id", tmp_path / "out.csv", name_column="missing_name"
        )


def test_verify_chebi_file_missing_cas_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["chebi_id"], [["CHEBI:16236"]])
    with pytest.raises(ChebiTermsError, match="--cas-column"):
        verify_chebi_file(
            src, "chebi_id", tmp_path / "out.csv", name_column=None, cas_column="nope"
        )


@patch("chebi_terms.client.time.sleep")
@patch("chebi_terms.client.requests.get")
def test_verify_chebi_file_reads_utf8_bom(
    mock_get: MagicMock, _sleep: MagicMock, tmp_path: Path
) -> None:
    mock_get.side_effect = route_chebi_get
    src = tmp_path / "in.csv"
    src.write_text("chebi_id\nCHEBI:16236\n", encoding="utf-8-sig")
    out = tmp_path / "out.csv"

    verify_chebi_file(src, "chebi_id", out)

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["chebi_name"] == "ethanol"
