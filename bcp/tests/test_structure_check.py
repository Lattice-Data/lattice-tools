"""Unit tests for structure_check client and io (mocked HTTP)."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from structure_check.client import (
    CAS_UNRESOLVED,
    CHEBI_NO_STRUCTURE,
    CHEBI_UNREACHABLE,
    CHEBI_UNRESOLVED,
    MATCH,
    NAME_UNRESOLVED,
    NOT_CHECKED,
    OUTPUT_FIELDS_APPENDED,
    REVIEW_CHECK,
    REVIEW_INVESTIGATE,
    REVIEW_OK,
    REVIEW_UNVERIFIED,
    SALT_DIFFERS,
    SKELETON_DIFFERS,
    STEREO_DIFFERS,
    check_row,
    compare_structures,
    empty_result,
    name_candidates,
    review_level,
    skeleton,
)
from structure_check.io import (
    RunSummary,
    SideTally,
    StructureCheckError,
    check_file,
    check_output_path,
    emit_single_row,
)

# Real InChIKeys, so the skeleton/stereo split is exercised against genuine data.
ETHANOL = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
UCN01 = "PBCZSGKMGDDXIJ-KRUBCLEUSA-N"
UCN01_OTHER_EPIMER = "PBCZSGKMGDDXIJ-HQCWYSJUSA-N"  # same skeleton, different stereo
ALEXIDINE = "BRJJFBHTDVWTCJ-UHFFFAOYSA-N"


# --------------------------------------------------------------------------
# InChIKey comparison
# --------------------------------------------------------------------------


def test_skeleton_is_the_connectivity_block() -> None:
    assert skeleton(ETHANOL) == "LFQSCWFLJHTTHZ"


def test_skeleton_of_nothing_is_empty() -> None:
    assert skeleton("") == ""


def test_compare_identical_keys_is_a_match() -> None:
    assert compare_structures(ETHANOL, [ETHANOL]) == MATCH


def test_compare_same_skeleton_different_stereo() -> None:
    """The distinction the whole tool rests on: isomer, not different molecule."""
    assert compare_structures(UCN01, [UCN01_OTHER_EPIMER]) == STEREO_DIFFERS


def test_compare_different_skeleton() -> None:
    assert compare_structures(ETHANOL, [ALEXIDINE]) == SKELETON_DIFFERS


def test_compare_matches_against_any_candidate() -> None:
    # A name resolving to several entries is not evidence of an error.
    assert compare_structures(ETHANOL, [ALEXIDINE, ETHANOL]) == MATCH


def test_compare_prefers_exact_match_over_stereo_sibling() -> None:
    assert compare_structures(UCN01, [UCN01_OTHER_EPIMER, UCN01]) == MATCH


# --------------------------------------------------------------------------
# name candidate ladder
# --------------------------------------------------------------------------


def test_name_candidates_tries_the_raw_string_first() -> None:
    whole, tokens = name_candidates("Scriptaid")
    assert whole == ["Scriptaid"]
    assert tokens == []


def test_name_candidates_adds_an_underscore_free_form() -> None:
    whole, _tokens = name_candidates("Alexidine_dihydrochloride")
    assert whole == ["Alexidine_dihydrochloride", "Alexidine dihydrochloride"]


@pytest.mark.parametrize(
    "raw",
    [
        "PIK-75_HCl",
        "Alexidine_dihydrochloride",
        "Omeprazole_Sodium",
        "Ara-G_hydrate",
        "Olanexidine_Hydrochloride_semihydrate",
        "PARAROSANILINE_PAMOATE",
    ],
)
def test_name_candidates_never_queries_a_salt_free_base(raw: str) -> None:
    """
    'PIK-75 HCl' and 'PIK-75' are different molecules.

    Falling back to tokens here would compare the free base against a salt's CAS
    and report a difference of our own manufacture.
    """
    whole, tokens = name_candidates(raw)
    assert tokens == []
    # Every form still carries the salt word.
    salt_word = raw.split("_")[-1].lower()
    assert whole and all(salt_word in form.lower() for form in whole)


def test_name_candidates_offers_tokens_for_dual_alias_cells() -> None:
    _whole, tokens = name_candidates("ABT_263_Navitoclax")
    assert tokens == ["ABT", "Navitoclax"]  # '263' has no letters


def test_name_candidates_skips_short_and_numeric_tokens() -> None:
    _whole, tokens = name_candidates("Epirubicin_IMI_28")
    assert tokens == ["Epirubicin", "IMI"]


def test_name_candidates_of_blank_is_empty() -> None:
    assert name_candidates("   ") == ([], [])


# --------------------------------------------------------------------------
# review level
# --------------------------------------------------------------------------


def test_skeleton_difference_outranks_everything() -> None:
    assert review_level(MATCH, SKELETON_DIFFERS) == REVIEW_INVESTIGATE
    assert review_level(SKELETON_DIFFERS, STEREO_DIFFERS) == REVIEW_INVESTIGATE


def test_stereo_difference_is_a_softer_flag() -> None:
    assert review_level(MATCH, STEREO_DIFFERS) == REVIEW_CHECK


def test_salt_difference_is_a_softer_flag() -> None:
    # 'Doxorubicin' with the hydrochloride's CAS is looseness, not a wrong row.
    assert review_level(MATCH, SALT_DIFFERS) == REVIEW_CHECK


def test_salt_difference_never_outranks_a_real_one() -> None:
    assert review_level(SKELETON_DIFFERS, SALT_DIFFERS) == REVIEW_INVESTIGATE


def test_one_successful_comparison_is_enough_to_be_ok() -> None:
    # The ID/CAS pairing is what the curator feared; an unresolvable name does not
    # undo that, and its own column still records why it was not compared.
    assert review_level(MATCH, NAME_UNRESOLVED) == REVIEW_OK


def test_nothing_comparable_is_unverified_not_ok() -> None:
    assert review_level(NOT_CHECKED, NOT_CHECKED) == REVIEW_UNVERIFIED
    assert review_level(CHEBI_UNRESOLVED, NAME_UNRESOLVED) == REVIEW_UNVERIFIED


def test_empty_result_starts_unverified() -> None:
    from structure_check.client import IDENTIFIER_NOT_CHECKED, STATUS_FIELDS

    result = empty_result()
    assert result["review"] == REVIEW_UNVERIFIED
    assert set(result) == set(OUTPUT_FIELDS_APPENDED) | set(STATUS_FIELDS)
    # Nothing was asked, so nothing failed.
    assert all(result[f] == IDENTIFIER_NOT_CHECKED for f in STATUS_FIELDS)
    assert result["unasked"] == ""


# --------------------------------------------------------------------------
# salt-form refinement
# --------------------------------------------------------------------------

DOXORUBICIN_FREE_BASE = "AOJJSUZBOXZQNB-TZSSRYMLSA-N"
DOXORUBICIN_HCL = "MWWSFMDVAYGXBV-RUELKSSGSA-N"


@patch("structure_check.client.parent_inchikey")
def test_refine_downgrades_a_shared_parent_to_salt_differs(
    mock_parent: MagicMock,
) -> None:
    """A counterion is part of the skeleton block; a shared parent collapses it."""
    from structure_check.client import refine_skeleton_difference

    mock_parent.return_value = DOXORUBICIN_FREE_BASE
    assert (
        refine_skeleton_difference(DOXORUBICIN_HCL, [DOXORUBICIN_FREE_BASE])
        == SALT_DIFFERS
    )


@patch("structure_check.client.parent_inchikey")
def test_refine_keeps_a_genuine_difference(mock_parent: MagicMock) -> None:
    from structure_check.client import refine_skeleton_difference

    mock_parent.side_effect = lambda key: {ETHANOL: ETHANOL, ALEXIDINE: ALEXIDINE}[key]
    assert refine_skeleton_difference(ETHANOL, [ALEXIDINE]) == SKELETON_DIFFERS


@patch("structure_check.client.parent_inchikey", return_value="")
def test_refine_leaves_the_row_flagged_when_parents_are_unavailable(
    _mock_parent: MagicMock,
) -> None:
    """
    PubChem's parent is wrong for some multi-component salts (pyrvinium pamoate
    resolves to the pamoic acid counterion), so an unusable answer must only ever
    leave a row flagged, never wave one through.
    """
    from structure_check.client import refine_skeleton_difference

    assert refine_skeleton_difference(ETHANOL, [ALEXIDINE]) == SKELETON_DIFFERS


# --------------------------------------------------------------------------
# check_row, with both upstream lookups mocked
# --------------------------------------------------------------------------


def _patch_lookups(cas_key=ETHANOL, cas_name="Ethanol", chebi=(ETHANOL, ""), name=None):
    """Patch the three resolvers structure_check composes."""
    if name is None:
        name = ("Ethanol", [ETHANOL], 1)
    return (
        patch("structure_check.client.cas_structure", return_value=(cas_key, cas_name)),
        patch("structure_check.client.chebi_structure", return_value=chebi),
        patch("structure_check.client.name_structure", return_value=name),
        # No network in unit tests: a skeleton difference stays one.
        patch(
            "structure_check.client.refine_skeleton_difference",
            side_effect=lambda ref, cands: SKELETON_DIFFERS,
        ),
    )


def test_check_row_all_three_agree() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups()
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol")
    assert result["id_cas_verdict"] == MATCH
    assert result["name_cas_verdict"] == MATCH
    assert result["review"] == REVIEW_OK
    assert result["cas_inchikey"] == ETHANOL


def test_check_row_isolates_a_wrong_name_from_a_right_id() -> None:
    """The Alexidine shape: ChEBI matches the CAS, the row's own name does not."""
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        name=("Alexidine", [ALEXIDINE], 1)
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="22573-88-2", chebi_id="CHEBI:27391", name="Alexidine")
    assert result["id_cas_verdict"] == MATCH
    assert result["name_cas_verdict"] == SKELETON_DIFFERS
    assert result["review"] == REVIEW_INVESTIGATE


def test_check_row_flags_a_chebi_id_holding_another_molecule() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(chebi=(ALEXIDINE, ""))
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", chebi_id="CHEBI:27391", name="Ethanol")
    assert result["id_cas_verdict"] == SKELETON_DIFFERS
    assert result["review"] == REVIEW_INVESTIGATE


def test_check_row_reports_stereo_difference_as_check_not_investigate() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        cas_key=UCN01_OTHER_EPIMER,
        chebi=(UCN01_OTHER_EPIMER, ""),
        name=("UCN-01", [UCN01], 1),
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="112953-11-4", chebi_id="CHEBI:221840", name="UCN-01")
    assert result["name_cas_verdict"] == STEREO_DIFFERS
    assert result["review"] == REVIEW_CHECK


def test_check_row_unresolved_name_is_not_a_finding() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(name=("", [], 0))
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", chebi_id="CHEBI:16236", name="Widget_9000")
    assert result["name_cas_verdict"] == NAME_UNRESOLVED
    assert result["name_cas_verdict"] != SKELETON_DIFFERS
    assert result["review"] == REVIEW_OK


def test_check_row_unresolved_cas_blocks_both_comparisons() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(cas_key="", cas_name="")
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="00-00-0", chebi_id="CHEBI:16236", name="Ethanol")
    assert result["id_cas_verdict"] == CAS_UNRESOLVED
    assert result["name_cas_verdict"] == CAS_UNRESOLVED
    assert result["review"] == REVIEW_UNVERIFIED


def test_check_row_blank_cas_is_not_checked_rather_than_unresolved() -> None:
    """A blank CAS cell was never asked about; that is not the same failure as
    PubChem drawing a blank on a CAS it was given, so it must not be counted the
    same way (dead_checks divides unresolved rows by attempted rows)."""
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(cas_key="", cas_name="")
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="", chebi_id="CHEBI:16236", name="Ethanol")
    assert result["id_cas_verdict"] == NOT_CHECKED
    assert result["name_cas_verdict"] == NOT_CHECKED
    assert result["review"] == REVIEW_UNVERIFIED


def test_check_row_chebi_without_a_structure() -> None:
    """Class terms and R-group entries carry no InChIKey."""
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(chebi=("", CHEBI_NO_STRUCTURE))
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", chebi_id="CHEBI:33697", name="Ethanol")
    assert result["id_cas_verdict"] == CHEBI_NO_STRUCTURE
    assert result["review"] == REVIEW_OK  # the name check still succeeded


def test_check_row_skips_a_check_it_was_given_nothing_for() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups()
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", chebi_id="CHEBI:16236")
    assert result["name_cas_verdict"] == NOT_CHECKED
    assert result["id_cas_verdict"] == MATCH


def test_check_row_records_the_query_it_actually_used() -> None:
    """A difference must never be traceable to a query the caller cannot see."""
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        name=("tokens: Vorinostat|SAHA", [ETHANOL], 1)
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", name="Vorinostat_SAHA")
    assert result["name_query"] == "tokens: Vorinostat|SAHA"


# --------------------------------------------------------------------------
# resolvers against mocked HTTP
# --------------------------------------------------------------------------


def _pubchem_name_response(keys: list[str]) -> MagicMock:
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = {
        "PropertyTable": {"Properties": [{"InChIKey": k} for k in keys]}
    }
    return resp


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_inchikeys_for_name_collects_every_key(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    from structure_check.client import inchikeys_for_name

    mock_get.return_value = _pubchem_name_response([ETHANOL, ALEXIDINE, ETHANOL])
    # Deduplicated, order preserved.
    assert inchikeys_for_name("ethanol") == [ETHANOL, ALEXIDINE]


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_inchikeys_for_name_url_encodes_awkward_names(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """Chemical names carry slashes, parens, and commas."""
    from structure_check.client import inchikeys_for_name

    mock_get.return_value = _pubchem_name_response([ETHANOL])
    inchikeys_for_name("4-(4-((4'-Chloro(1,1'-biphenyl)-2-yl)methyl)/x)")
    url = mock_get.call_args[0][0]
    assert "/" not in url.split("/compound/name/")[1].split("/property")[0]
    assert "(" not in url.split("/compound/name/")[1].split("/property")[0]


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_name_structure_stops_at_the_first_whole_string_hit(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """Salt fidelity: if the raw form resolves, tokens are never tried."""
    from structure_check.client import name_structure

    mock_get.return_value = _pubchem_name_response([ETHANOL])
    query, keys, total = name_structure("PIK-75_HCl")
    assert query == "PIK-75_HCl"
    assert keys == [ETHANOL]
    assert total == 1
    assert mock_get.call_count == 1


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_name_structure_falls_back_to_a_token_union(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    from structure_check.client import name_structure

    empty = _pubchem_name_response([])
    empty.status_code = 404
    responses = {
        "ABT_263_Navitoclax": empty,
        "ABT%20263%20Navitoclax": empty,
        "ABT": empty,
        "Navitoclax": _pubchem_name_response([ALEXIDINE]),
    }

    def route(url: str, *args, **kwargs):
        name = url.split("/compound/name/")[1].split("/property")[0]
        return responses.get(name, empty)

    mock_get.side_effect = route
    query, keys, _total = name_structure("ABT_263_Navitoclax")
    assert query == "tokens: Navitoclax"
    assert keys == [ALEXIDINE]


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_name_structure_unresolved_returns_no_keys(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    from structure_check.client import name_structure

    resp = MagicMock()
    resp.status_code = 404
    mock_get.return_value = resp
    assert name_structure("Widget_9000") == ("", [], 0)


@patch("structure_check.client.fetch_compound")
def test_chebi_structure_reads_the_default_structure(mock_fetch: MagicMock) -> None:
    from structure_check.client import chebi_structure

    mock_fetch.return_value = {"default_structure": {"standard_inchi_key": ETHANOL}}
    assert chebi_structure("CHEBI:16236") == (ETHANOL, "")


@patch("structure_check.client.fetch_compound")
def test_chebi_structure_without_a_structure_block(mock_fetch: MagicMock) -> None:
    from structure_check.client import chebi_structure

    mock_fetch.return_value = {"default_structure": None}
    assert chebi_structure("CHEBI:33697") == ("", CHEBI_NO_STRUCTURE)


@patch("structure_check.client.fetch_compound")
def test_chebi_structure_missing_record(mock_fetch: MagicMock) -> None:
    from structure_check.client import chebi_structure

    mock_fetch.return_value = None
    assert chebi_structure("CHEBI:99999999") == ("", CHEBI_UNRESOLVED)


@patch("structure_check.client.fetch_compound")
def test_chebi_structure_unreachable_is_not_a_finding(mock_fetch: MagicMock) -> None:
    from chebi_terms.client import ChebiUnavailableError
    from structure_check.client import chebi_structure

    mock_fetch.side_effect = ChebiUnavailableError("down")
    # Distinct from "no such record": only a universal outage implies degradation.
    assert chebi_structure("CHEBI:16236") == ("", CHEBI_UNREACHABLE)


def test_chebi_structure_rejects_junk_without_fetching() -> None:
    from structure_check.client import chebi_structure

    with patch("structure_check.client.fetch_compound") as mock_fetch:
        assert chebi_structure("not-an-id") == ("", CHEBI_UNRESOLVED)
        mock_fetch.assert_not_called()


# --------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


@patch("structure_check.io.check_row")
def test_check_file_appends_verdicts_and_keeps_input_columns(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    mock_check.return_value = {**empty_result(), "review": REVIEW_OK}
    src = tmp_path / "in.csv"
    _write_csv(src, ["Name", "CAS", "note"], [["Ethanol", "64-17-5", "keep"]])
    out = tmp_path / "out.csv"

    check_file(src, out, cas_column="CAS", name_column="Name")

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["note"] == "keep"
    assert rows[0]["review"] == REVIEW_OK


@patch("structure_check.io.check_row")
def test_check_file_caches_identical_triples(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    mock_check.return_value = empty_result()
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["Name", "CAS"],
        [["Ethanol", "64-17-5"], ["Ethanol", "64-17-5"], ["Other", "64-17-5"]],
    )
    check_file(src, tmp_path / "out.csv", cas_column="CAS", name_column="Name")
    assert mock_check.call_count == 2


@patch("structure_check.io.check_row")
def test_check_file_summary_counts_reviews(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    mock_check.side_effect = [
        {**empty_result(), "review": REVIEW_INVESTIGATE},
        {**empty_result(), "review": REVIEW_CHECK},
        {**empty_result(), "review": REVIEW_OK},
    ]
    src = tmp_path / "in.csv"
    _write_csv(src, ["Name", "CAS"], [["a", "1"], ["b", "2"], ["c", "3"]])

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", name_column="Name"
    )

    assert isinstance(summary, RunSummary)
    assert summary.rows == 3
    assert summary.review_counts[REVIEW_INVESTIGATE] == 1
    assert summary.needs_attention == 2


def test_check_file_requires_something_to_check_against(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["CAS"], [["64-17-5"]])
    with pytest.raises(StructureCheckError, match="Nothing to check"):
        check_file(src, tmp_path / "out.csv", cas_column="CAS")


def test_check_file_missing_cas_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["Name"], [["Ethanol"]])
    with pytest.raises(StructureCheckError, match="--cas-column"):
        check_file(src, tmp_path / "out.csv", cas_column="CAS", name_column="Name")


def test_check_file_missing_input(tmp_path: Path) -> None:
    with pytest.raises(StructureCheckError, match="Input file not found"):
        check_file(tmp_path / "nope.csv", tmp_path / "out.csv", cas_column="CAS")


def test_check_file_rejects_reserved_column(tmp_path: Path) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["review", "CAS"], [["x", "64-17-5"]])
    with pytest.raises(StructureCheckError, match="collides"):
        check_file(src, tmp_path / "out.csv", cas_column="CAS", name_column="review")


def test_check_file_refuses_to_overwrite_its_input(tmp_path: Path) -> None:
    src = tmp_path / "sheet.csv"
    _write_csv(src, ["Name", "CAS"], [["Ethanol", "64-17-5"]])
    before = src.read_text(encoding="utf-8")
    with pytest.raises(StructureCheckError, match="input file"):
        check_file(src, src, cas_column="CAS", name_column="Name")
    assert src.read_text(encoding="utf-8") == before


def test_check_output_path_rejects_a_directory(tmp_path: Path) -> None:
    with pytest.raises(StructureCheckError, match="is a directory"):
        check_output_path(tmp_path)


def test_check_output_path_rejects_a_missing_parent(tmp_path: Path) -> None:
    with pytest.raises(StructureCheckError, match="does not exist"):
        check_output_path(tmp_path / "nope" / "out.csv")


@patch("structure_check.io.check_row")
def test_emit_single_row_json_to_stdout(
    mock_check: MagicMock, capsys: pytest.CaptureFixture
) -> None:
    mock_check.return_value = {**empty_result(), "review": REVIEW_INVESTIGATE}
    emit_single_row(None, name="Alexidine", cas="22573-88-2", chebi_id="CHEBI:27391")
    payload = json.loads(capsys.readouterr().out)
    assert payload["review"] == REVIEW_INVESTIGATE
    assert payload["cas"] == "22573-88-2"


@patch("structure_check.io.check_row")
def test_emit_single_row_csv_to_file(mock_check: MagicMock, tmp_path: Path) -> None:
    mock_check.return_value = {**empty_result(), "review": REVIEW_OK}
    out = tmp_path / "one.csv"
    emit_single_row(out, name="Ethanol", cas="64-17-5", fmt="csv")
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["review"] == REVIEW_OK


@patch("structure_check.io.check_row")
def test_emit_single_row_checks_output_before_looking_anything_up(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    with pytest.raises(StructureCheckError):
        emit_single_row(tmp_path / "nope" / "out.json", cas="64-17-5", name="Ethanol")
    mock_check.assert_not_called()


# --------------------------------------------------------------------------
# resolvers with no mocking above them: the request chains and JSON shapes
# --------------------------------------------------------------------------


def _json_response(payload: dict) -> MagicMock:
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = payload
    return resp


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_structure_uses_two_requests_not_four(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """lookup_cas costs four; only InChIKey and Title are ever used."""
    from structure_check.client import cas_structure

    mock_get.side_effect = [
        _json_response({"IdentifierList": {"CID": [702]}}),
        _json_response(
            {
                "PropertyTable": {
                    "Properties": [{"InChIKey": ETHANOL, "Title": "Ethanol"}]
                }
            }
        ),
    ]
    assert cas_structure("64-17-5") == (ETHANOL, "Ethanol")
    assert mock_get.call_count == 2


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_structure_unresolved_cas(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    from structure_check.client import cas_structure

    missing = MagicMock()
    missing.status_code = 404
    mock_get.return_value = missing
    assert cas_structure("00-00-0") == ("", "")


def test_cas_structure_blank_costs_nothing() -> None:
    from structure_check.client import cas_structure

    with patch("chebi_lookup.client.requests.get") as mock_get:
        assert cas_structure("   ") == ("", "")
        mock_get.assert_not_called()


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_walks_the_three_request_chain(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """Three distinct JSON shapes, none of them exercised by the mocked tests above."""
    import structure_check.client as sc

    sc._parent_cache.clear()
    mock_get.side_effect = [
        _json_response({"IdentifierList": {"CID": [443939]}}),  # inchikey -> cid
        _json_response({"IdentifierList": {"CID": [31703]}}),  # cid -> parent cid
        _json_response(
            {"PropertyTable": {"Properties": [{"InChIKey": DOXORUBICIN_FREE_BASE}]}}
        ),
    ]
    assert sc.parent_inchikey(DOXORUBICIN_HCL) == DOXORUBICIN_FREE_BASE
    assert mock_get.call_count == 3
    # The parent request must ask for the parent, not just any related CID.
    assert "cids_type=parent" in mock_get.call_args_list[1][0][0]


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_caches_successes(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    import structure_check.client as sc

    sc._parent_cache.clear()
    mock_get.side_effect = [
        _json_response({"IdentifierList": {"CID": [443939]}}),
        _json_response({"IdentifierList": {"CID": [31703]}}),
        _json_response(
            {"PropertyTable": {"Properties": [{"InChIKey": DOXORUBICIN_FREE_BASE}]}}
        ),
    ]
    sc.parent_inchikey(DOXORUBICIN_HCL)
    sc.parent_inchikey(DOXORUBICIN_HCL)
    assert mock_get.call_count == 3  # second call served from cache


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_does_not_cache_an_outage(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """One transient failure must not disable salt demotion for the whole run."""
    import structure_check.client as sc

    sc._parent_cache.clear()
    down = MagicMock()
    down.status_code = 503
    mock_get.return_value = down

    assert sc.parent_inchikey(DOXORUBICIN_HCL) == ""
    assert DOXORUBICIN_HCL not in sc._parent_cache


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_caches_a_definitive_no_parent(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """
    "PubChem has no parent for this structure" is an answer, not a silence.

    It will not change during the run, so re-walking three requests for every
    repeat of the same structure buys nothing.
    """
    import structure_check.client as sc

    sc._parent_cache.clear()
    missing = MagicMock()
    missing.status_code = 404
    mock_get.return_value = missing

    assert sc.parent_inchikey(DOXORUBICIN_HCL) == ""
    assert sc._parent_cache[DOXORUBICIN_HCL] == ""
    calls = mock_get.call_count
    assert sc.parent_inchikey(DOXORUBICIN_HCL) == ""
    assert mock_get.call_count == calls  # served from cache, no new requests


@pytest.mark.parametrize(
    "payload",
    [
        {"IdentifierList": {"CID": []}},
        {"IdentifierList": {}},
        {"unexpected": "shape"},
    ],
)
@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_survives_unexpected_json(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock, payload: dict
) -> None:
    import structure_check.client as sc

    sc._parent_cache.clear()
    mock_get.return_value = _json_response(payload)
    assert sc.parent_inchikey(DOXORUBICIN_HCL) == ""


def test_parent_inchikey_of_nothing_costs_nothing() -> None:
    import structure_check.client as sc

    with patch("chebi_lookup.client.requests.get") as mock_get:
        assert sc.parent_inchikey("") == ""
        mock_get.assert_not_called()


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_name_structure_caps_a_runaway_name_and_reports_the_total(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """Each extra candidate costs three more requests during refinement."""
    from structure_check.client import MAX_NAME_CANDIDATES, name_structure

    many = [f"AAAAAAAAAAAAA{n}-UHFFFAOYSA-N" for n in range(MAX_NAME_CANDIDATES + 15)]
    mock_get.return_value = _json_response(
        {"PropertyTable": {"Properties": [{"InChIKey": k} for k in many]}}
    )
    _query, keys, total = name_structure("acid")
    assert len(keys) == MAX_NAME_CANDIDATES
    # The caller needs the real total to know the comparison was incomplete.
    assert total == MAX_NAME_CANDIDATES + 15


# --------------------------------------------------------------------------
# a run that compared nothing is a broken run, not a clean sheet
# --------------------------------------------------------------------------


def test_degraded_when_no_row_compared_anything() -> None:
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 8}),
        verdict_counts=Counter(),
        rows=8,
    )
    assert summary.nothing_compared
    assert summary.degraded


def test_not_degraded_when_anything_resolved() -> None:
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 7, REVIEW_OK: 1}),
        verdict_counts=Counter(),
        rows=8,
    )
    assert not summary.degraded


def test_not_degraded_on_a_tiny_spot_check() -> None:
    """Two unverified rows is as likely a spot check as a broken run."""
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 2}), verdict_counts=Counter(), rows=2
    )
    assert not summary.degraded


def test_nothing_compared_survives_a_column_that_resolved_perfectly() -> None:
    """
    The case the identifier-level signals cannot see.

    A ChEBI column of class terms resolves every single value — ChEBI holds the
    records — and still compares nothing, because a class term carries no
    structure. Both other signals are silent, so this term has to exist.
    """
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 117}),
        verdict_counts=Counter(),
        rows=117,
        cas=SideTally(resolved_values=117, resolved_rows=117),
        chebi=SideTally(resolved_values=117, resolved_rows=117),
    )
    assert not summary.outages
    assert not summary.suspect_columns
    assert summary.nothing_compared
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_check_file_flags_a_wrong_cas_column_as_degraded(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """--cas-column pointed at a note column: every lookup 404s, nothing compared."""
    mock_check.return_value = empty_result()
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["Name", "note"],
        [[f"Compound {n}", f"internal note {n}"] for n in range(6)],
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="note", name_column="Name"
    )

    assert summary.review_counts[REVIEW_UNVERIFIED] == 6
    assert summary.degraded


@pytest.mark.parametrize("other", ["chebi_column", "name_column"])
def test_check_file_rejects_a_column_that_is_also_the_cas_column(
    tmp_path: Path, other: str
) -> None:
    src = tmp_path / "in.csv"
    _write_csv(src, ["CAS"], [["64-17-5"]])
    with pytest.raises(StructureCheckError, match="same as --cas-column"):
        check_file(src, tmp_path / "out.csv", cas_column="CAS", **{other: "CAS"})


# --------------------------------------------------------------------------
# review_rank
# --------------------------------------------------------------------------


def test_review_rank_sorts_worst_first() -> None:
    """Alphabetically 'check' < 'investigate'; the rank exists to fix that."""
    from structure_check.client import REVIEW_RANK

    ordered = sorted(REVIEW_RANK, key=lambda level: REVIEW_RANK[level])
    assert ordered == [
        REVIEW_INVESTIGATE,
        REVIEW_CHECK,
        REVIEW_UNVERIFIED,
        REVIEW_OK,
    ]


def test_check_row_sets_a_rank_matching_its_review() -> None:
    from structure_check.client import REVIEW_RANK

    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        name=("Alexidine", [ALEXIDINE], 1)
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="22573-88-2", chebi_id="CHEBI:27391", name="Alexidine")
    assert result["review"] == REVIEW_INVESTIGATE
    assert result["review_rank"] == REVIEW_RANK[REVIEW_INVESTIGATE]


def test_check_row_with_nothing_to_compare_spends_no_requests() -> None:
    with patch("structure_check.client.cas_structure") as mock_cas:
        result = check_row(cas="64-17-5")
        mock_cas.assert_not_called()
    assert result["review"] == REVIEW_UNVERIFIED


# --------------------------------------------------------------------------
# compare_structures is public, so "could not compare" must not read as a difference
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "reference, candidates",
    [("", [ETHANOL]), (ETHANOL, []), (ETHANOL, [""]), ("", [])],
)
def test_compare_reports_incomparable_rather_than_a_difference(
    reference: str, candidates: list[str]
) -> None:
    from structure_check.client import NOT_COMPARABLE

    assert compare_structures(reference, candidates) == NOT_COMPARABLE


# --------------------------------------------------------------------------
# the salt guard's non-obvious branches
#
# The six cases parametrized above all resolve through exact SALT_TOKENS
# membership, so none of them touched the machinery added to close the escapes.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("raw", ["PIK-75_2HCl", "Compound_3HCl", "Thing_2hydrate"])
def test_salt_guard_strips_a_count_prefix_digit(raw: str) -> None:
    """'2HCl' is not a SALT_TOKENS entry; the digit strip is what catches it."""
    _whole, tokens = name_candidates(raw)
    assert tokens == []


@pytest.mark.parametrize("fragment", ["trisodium", "dipotassium", "dinitrate"])
def test_salt_guard_strips_a_multiplier_word(fragment: str) -> None:
    """
    These match only via _COUNT_PREFIXES.

    Not redundant with the suffix list: 'sodium', 'potassium', and 'nitrate' are
    SALT_TOKENS entries but deliberately not _SALT_SUFFIXES, so nothing else here
    would catch them.
    """
    from structure_check.client import SALT_TOKENS, _is_salt_word, _SALT_SUFFIXES

    assert fragment not in SALT_TOKENS
    assert not fragment.endswith(_SALT_SUFFIXES)
    assert _is_salt_word(fragment)
    assert name_candidates(f"Compound_{fragment}")[1] == []


@pytest.mark.parametrize(
    "fragment", ["Doxorubicinhcl", "compoundhydrogensulfate", "thingpamoate"]
)
def test_salt_guard_matches_a_counterion_suffix(fragment: str) -> None:
    """Run-together spellings that are neither a token nor prefix+token."""
    from structure_check.client import SALT_TOKENS, _is_salt_word

    assert fragment.casefold() not in SALT_TOKENS
    assert _is_salt_word(fragment)


def test_salt_guard_sees_a_counterion_behind_a_space() -> None:
    """
    'Doxorubicin HCl_Adriamycin' hides its counterion behind a space.

    Splitting the guard on '_' alone would leave the fallback enabled and compare
    the free base 'Adriamycin' against the hydrochloride's CAS.
    """
    whole, tokens = name_candidates("Doxorubicin HCl_Adriamycin")
    assert tokens == []
    assert whole[0] == "Doxorubicin HCl_Adriamycin"


def test_salt_guard_leaves_genuine_alias_pairs_alone() -> None:
    """The guard must not swallow the case the fallback exists for."""
    for raw in ("Vorinostat_SAHA", "ABT_263_Navitoclax", "Actinomycin_D_Dactinomycin"):
        assert name_candidates(raw)[1], raw


# --------------------------------------------------------------------------
# truncation must not manufacture a finding
# --------------------------------------------------------------------------


def test_truncated_name_with_no_match_is_ambiguous_not_a_difference() -> None:
    """
    PubChem returns structures in CID order, not relevance order, so the true match
    may sit past the cap. Calling that a different molecule would be a finding of
    the tool's own making.
    """
    from structure_check.client import MAX_NAME_CANDIDATES, NAME_AMBIGUOUS

    capped = [f"BBBBBBBBBBBBB{n}-UHFFFAOYSA-N" for n in range(MAX_NAME_CANDIDATES)]
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        name=("acid", capped, MAX_NAME_CANDIDATES + 30)
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", name="acid")

    assert result["name_cas_verdict"] == NAME_AMBIGUOUS
    assert result["name_cas_verdict"] != SKELETON_DIFFERS


def test_truncated_name_records_the_truncation_in_the_row() -> None:
    """A flagged row has to be auditable from the CSV alone."""
    from structure_check.client import MAX_NAME_CANDIDATES

    capped = [f"BBBBBBBBBBBBB{n}-UHFFFAOYSA-N" for n in range(MAX_NAME_CANDIDATES)]
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(name=("acid", capped, 40))
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", name="acid")
    assert "truncated: 40" in result["name_query"]


def test_a_match_inside_a_truncated_set_is_still_a_match() -> None:
    """Matching any candidate is sound; only a non-match is unsafe when truncated."""
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(
        name=("acid", [ALEXIDINE, ETHANOL], 40)
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="64-17-5", name="acid")
    assert result["name_cas_verdict"] == MATCH


def test_truncated_skeleton_difference_skips_refinement() -> None:
    """
    A truncated non-match is name_ambiguous no matter what refinement finds, so
    calling it would spend up to 3 requests per candidate on a result that is
    discarded either way.
    """
    from structure_check.client import MAX_NAME_CANDIDATES, NAME_AMBIGUOUS

    capped = [f"BBBBBBBBBBBBB{n}-UHFFFAOYSA-N" for n in range(MAX_NAME_CANDIDATES)]
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(name=("acid", capped, 40))
    with cas_p, chebi_p, name_p, refine_p as mock_refine:
        result = check_row(cas="64-17-5", name="acid")
    assert result["name_cas_verdict"] == NAME_AMBIGUOUS
    mock_refine.assert_not_called()


# --------------------------------------------------------------------------
# untrusted sheet cells become URL path segments
# --------------------------------------------------------------------------


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_structure_quotes_the_sheet_cell(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """A stray '/' or '?' in a CAS cell must not rewrite the request path."""
    from structure_check.client import cas_structure

    missing = MagicMock()
    missing.status_code = 404
    mock_get.return_value = missing

    cas_structure("64-17-5/../../evil?x=1")

    path = mock_get.call_args_list[0][0][0].split("/compound/name/")[1]
    segment = path.split("/cids")[0]
    assert "/" not in segment
    assert "?" not in segment


# --------------------------------------------------------------------------
# is this run trustworthy: outage vs wrong column vs nothing compared
#
# One boolean answering two questions was wrong four rounds running, each fix
# correct for its own input and blind to the next. These cases are the input
# shapes that broke it, kept together so the next change has to face all of them
# at once rather than one at a time.
# --------------------------------------------------------------------------


def _summary(**kwargs) -> RunSummary:
    defaults = dict(review_counts=Counter(), verdict_counts=Counter(), rows=117)
    return RunSummary(**{**defaults, **kwargs})


# --- outage: the majority rule --------------------------------------------


def test_chebi_down_run_wide_while_names_matched_is_an_outage() -> None:
    """
    The hole round 2 closed, restated.

    ChEBI unreachable, names fine: every row gets a name match, review_level says
    ok, review_rank is 4, and the ID question a curator actually asked was never
    answered once. The review column alone cannot show this.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 117}),
        chebi=SideTally(outage_rows=117),
        name=SideTally(resolved_values=117, resolved_rows=117),
    )
    assert summary.review_counts[REVIEW_UNVERIFIED] == 0
    assert summary.outages == ("ChEBI ID",)
    assert summary.degraded


def test_a_few_transient_failures_do_not_condemn_a_healthy_run() -> None:
    """
    The false-alarm direction, and the reason outages are not "fires on >= 1".

    This tool makes 3-44 requests per row against an API that rate-limits, and a
    throttle that outlasts its retries is indistinguishable from an outage. Three
    unlucky rows in 117 is a normal large run; exiting 1 on it teaches an operator
    to stop reading the exit code, which costs more than the case it would catch.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 114, REVIEW_UNVERIFIED: 3}),
        cas=SideTally(
            resolved_values=114, missing_values=0, resolved_rows=114, outage_rows=3
        ),
    )
    assert summary.outages == ()
    assert not summary.degraded


def test_an_outage_that_outweighs_what_it_answered_is_degraded() -> None:
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 57, REVIEW_UNVERIFIED: 60}),
        cas=SideTally(resolved_values=57, resolved_rows=57, outage_rows=60),
    )
    assert summary.outages == ("CAS",)
    assert summary.degraded


def test_an_outage_exactly_matching_what_resolved_is_not_a_majority() -> None:
    """The boundary is strict: a tie still has as much answered as unanswered."""
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 50}),
        cas=SideTally(resolved_values=50, resolved_rows=50, outage_rows=50),
    )
    assert summary.outages == ()


def test_the_majority_rule_needs_a_floor_because_it_degenerates_against_zero() -> None:
    """
    Against zero answers a single outage is a "majority", which is the normal
    shape of a sparse column, not of a broken run.

    117 rows of which only 2 carry a ChEBI ID, and both are unlucky. Without a
    floor on the outage count this condemns the run and throws away 115 good rows
    — the same shape as the retired-ID bug, with `unreachable` swapped for
    `missing`.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 115, REVIEW_UNVERIFIED: 2}),
        cas=SideTally(resolved_values=117, resolved_rows=117),
        chebi=SideTally(outage_rows=2),
    )
    assert summary.outages == ()
    assert not summary.degraded


def test_a_substantial_outage_still_fires_on_a_side_that_resolved_nothing() -> None:
    summary = _summary(chebi=SideTally(outage_rows=5))
    assert summary.outages == ("ChEBI ID",)


def test_one_blip_on_a_wrong_column_does_not_rename_the_diagnosis() -> None:
    """
    A wrong column resolves nothing by definition, so any single throttled
    request on it would otherwise read as an outage and suppress the real
    finding — sending the operator to re-run a sheet that will fail again for
    the reason nobody named.
    """
    summary = _summary(
        name=SideTally(resolved_values=0, missing_values=116, outage_rows=1),
    )
    assert summary.outages == ()
    assert summary.suspect_columns == ("name",)
    assert summary.degraded


def test_a_genuine_outage_still_outranks_a_wrong_column_diagnosis() -> None:
    """When the upstream really is down, it is a competing explanation again."""
    summary = _summary(
        chebi=SideTally(resolved_values=0, missing_values=6, outage_rows=40),
    )
    assert summary.outages == ("ChEBI ID",)
    assert summary.suspect_columns == ()


# --- suspect column: distinct values, per-side guard -----------------------


def test_a_column_of_junk_is_suspect() -> None:
    summary = _summary(
        review_counts=Counter({REVIEW_UNVERIFIED: 117}),
        chebi=SideTally(missing_values=117),
        name=SideTally(resolved_values=117, resolved_rows=117),
    )
    assert summary.suspect_columns == ("ChEBI ID",)
    assert summary.degraded


def test_two_retired_ids_on_a_big_sheet_are_a_finding_not_a_broken_run() -> None:
    """
    The live bug this redesign was called for.

    117 rows, only 2 carry a ChEBI ID, and both are genuinely obsolete. That is
    exactly the finding this tool exists to surface, and the old floor reported it
    as a broken run and exited 1, throwing the other 115 good rows away.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 115, REVIEW_UNVERIFIED: 2}),
        cas=SideTally(resolved_values=117, resolved_rows=117),
        chebi=SideTally(missing_values=2),
    )
    assert summary.suspect_columns == ()
    assert summary.outages == ()
    assert not summary.degraded


def test_repeating_one_junk_value_is_one_mistake_not_five_hundred() -> None:
    """Distinct values, not rows: a dose-response sheet repeats every compound."""
    summary = _summary(
        rows=500,
        review_counts=Counter({REVIEW_OK: 500}),
        chebi=SideTally(missing_values=1),
    )
    assert summary.suspect_columns == ()


def test_a_column_of_class_terms_is_a_real_chebi_column() -> None:
    """
    chebi_no_structure means ChEBI held the record and it carries no structure, as
    class terms and R-group entries do. That is a correct ChEBI ID, so it must not
    read as evidence that the column holds something else.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_UNVERIFIED: 117}),
        chebi=SideTally(resolved_values=117, resolved_rows=117),
    )
    assert summary.suspect_columns == ()


def test_an_outage_on_one_upstream_does_not_excuse_another_columns_junk() -> None:
    """
    The guard is per side, never global.

    A ChEBI column full of free text is rejected by a regex without one request,
    so a PubChem outage is logically incapable of explaining it. Blanket
    disqualification would hide a real misconfiguration behind an unrelated blip.
    """
    summary = _summary(
        chebi=SideTally(missing_values=117),
        cas=SideTally(resolved_values=50, resolved_rows=50, outage_rows=67),
    )
    assert "ChEBI ID" in summary.suspect_columns
    assert summary.outages == ("CAS",)


def test_an_all_blank_column_is_not_suspect() -> None:
    """A partly-filled sheet is normal input; blanks never enter these counts."""
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 117}),
        cas=SideTally(resolved_values=117, resolved_rows=117),
        name=SideTally(resolved_values=117, resolved_rows=117),
    )
    assert summary.suspect_columns == ()
    assert not summary.degraded


def test_an_empty_sheet_has_not_failed_anything() -> None:
    assert not _summary(rows=0).degraded


# --------------------------------------------------------------------------
# end to end through check_file: the counters must agree with what happened
#
# Every previous version of this verdict divided "rows whose cell was non-empty"
# by "rows whose verdict looked unresolved" — two numbers derived from different
# things, which is why one blank cell, and later one not_checked verdict, could
# defeat it. check_row now reports what actually happened per identifier, and
# check_file counts only that. These tests pin that agreement.
# --------------------------------------------------------------------------


def _row_result(**over) -> dict:
    return {**empty_result(), **over}


@patch("structure_check.io.check_row")
def test_check_file_counts_an_outage_against_rows_not_cells(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """One blank ChEBI cell among six rows used to hide a total ChEBI outage."""
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
        NOT_CHECKED,
    )

    def per_row(*, cas, chebi_id, name):
        return _row_result(
            review=REVIEW_OK,
            id_cas_verdict=CHEBI_UNREACHABLE if chebi_id else NOT_CHECKED,
            name_cas_verdict=MATCH,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=IDENTIFIER_UNREACHABLE if chebi_id else "not_checked",
            name_status=IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    rows = [[f"C{n}", f"{n}-00-0", f"CHEBI:{n}"] for n in range(5)]
    rows.append(["C5", "5-00-0", ""])  # the blank that used to hide everything
    _write_csv(src, ["Name", "CAS", "ChEBI ID"], rows)

    summary = check_file(
        src,
        tmp_path / "out.csv",
        cas_column="CAS",
        chebi_column="ChEBI ID",
        name_column="Name",
    )

    assert summary.chebi.outage_rows == 5
    assert summary.chebi.resolved_rows == 0
    assert summary.outages == ("ChEBI ID",)
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_check_file_does_not_count_a_cas_it_never_asked_about(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A filled CAS column with nothing to compare it against.

    check_row short-circuits before making a request, so those CAS values were
    never asked about. Counting them as failures would accuse a perfectly good
    --cas-column of being the wrong column.
    """
    mock_check.return_value = empty_result()  # all statuses not_checked
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["CAS", "ChEBI ID"],
        [[f"{n}-00-0", ""] for n in range(6)],
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    assert summary.cas.missing_values == 0
    assert summary.suspect_columns == ()
    # It is still a run that compared nothing, and that is what it says.
    assert summary.nothing_compared
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_check_file_counts_distinct_values_not_rows_for_a_column(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """A dose-response sheet repeats each compound; that is one value, not thirty."""
    from structure_check.client import CHEBI_UNRESOLVED, IDENTIFIER_MISSING

    mock_check.return_value = _row_result(
        review=REVIEW_OK,
        id_cas_verdict=CHEBI_UNRESOLVED,
        chebi_status=IDENTIFIER_MISSING,
    )
    src = tmp_path / "in.csv"
    _write_csv(src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", "junk"] for n in range(30)])

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    assert summary.chebi.missing_values == 1
    assert summary.suspect_columns == ()


@patch("structure_check.io.check_row")
def test_check_file_flags_a_wrong_column_across_enough_distinct_values(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    from structure_check.client import CHEBI_UNRESOLVED, IDENTIFIER_MISSING

    mock_check.return_value = _row_result(
        review=REVIEW_UNVERIFIED,
        id_cas_verdict=CHEBI_UNRESOLVED,
        chebi_status=IDENTIFIER_MISSING,
    )
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["CAS", "note"],
        [[f"{n}-00-0", f"internal note {n}"] for n in range(6)],
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="note"
    )

    assert summary.chebi.missing_values == 6
    assert summary.suspect_columns == ("ChEBI ID",)
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_check_file_does_not_cache_an_outage(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A transient failure on one triple must not be replayed as a settled answer.

    The same rule parent_inchikey already applies to its own cache: caching a
    failure would let one unlucky moment decide every row that repeats it.
    """
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )

    calls = []

    def per_row(*, cas, chebi_id, name):
        calls.append((cas, chebi_id, name))
        # Fails the first time, succeeds after — the shape of a transient blip.
        if len(calls) == 1:
            return _row_result(
                review=REVIEW_UNVERIFIED,
                id_cas_verdict=CHEBI_UNREACHABLE,
                chebi_status=IDENTIFIER_UNREACHABLE,
                cas_status=IDENTIFIER_RESOLVED,
            )
        return _row_result(
            review=REVIEW_OK,
            id_cas_verdict=MATCH,
            chebi_status=IDENTIFIER_RESOLVED,
            cas_status=IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    _write_csv(src, ["CAS", "ChEBI ID"], [["64-17-5", "CHEBI:16236"]] * 3)

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    # Row 1 outaged and was retried rather than replayed; rows 2-3 came from cache.
    assert len(calls) == 2
    assert summary.chebi.outage_rows == 1
    assert summary.chebi.resolved_rows == 2
    assert summary.outages == ()


@patch("structure_check.io.check_row")
def test_check_file_stops_rather_than_grinding_through_a_dead_upstream(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A wholly unreachable upstream is decided by the first few rows.

    Every remaining row would cost seconds of backoff to reach a verdict already
    settled, so the run stops and says so, as chebi_terms does.
    """
    from structure_check.client import CAS_UNREACHABLE, IDENTIFIER_UNREACHABLE

    mock_check.return_value = _row_result(
        review=REVIEW_UNVERIFIED,
        id_cas_verdict=CAS_UNREACHABLE,
        cas_status=IDENTIFIER_UNREACHABLE,
    )
    src = tmp_path / "in.csv"
    _write_csv(
        src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", f"CHEBI:{n}"] for n in range(40)]
    )
    out = tmp_path / "out.csv"

    with pytest.raises(StructureCheckError, match="could not be reached"):
        check_file(src, out, cas_column="CAS", chebi_column="ChEBI ID")

    # It stopped early rather than walking all 40 rows.
    assert mock_check.call_count < 40
    # The partial output is still there for inspection.
    assert out.exists()


# --------------------------------------------------------------------------
# an answer and a silence are different things, on the PubChem side too
#
# chebi_lookup.get_with_retry returns None for a 404 and for exhausted retries
# alike, which is why cas_unresolved could not be told from "PubChem was down".
# --------------------------------------------------------------------------


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_structure_raises_when_pubchem_is_unreachable(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    from structure_check.client import PubChemUnavailableError, cas_structure

    down = MagicMock()
    down.status_code = 503
    mock_get.return_value = down

    with pytest.raises(PubChemUnavailableError):
        cas_structure("64-17-5")


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_structure_still_answers_no_for_a_real_404(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """A 404 is PubChem answering. It must not be mistaken for an outage."""
    from structure_check.client import cas_structure

    missing = MagicMock()
    missing.status_code = 404
    mock_get.return_value = missing

    assert cas_structure("00-00-0") == ("", "")


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_two_hundred_that_is_not_json_is_an_outage_not_an_answer(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """
    A maintenance page served with HTTP 200 is something other than this API
    answering. Read as 'no such compound' it would look like a wrong column.
    """
    from structure_check.client import PubChemUnavailableError, cas_structure

    html = MagicMock()
    html.status_code = 200
    html.json.side_effect = ValueError("not json")
    mock_get.return_value = html

    with pytest.raises(PubChemUnavailableError):
        cas_structure("64-17-5")


@patch("structure_check.client.cas_structure")
@patch("structure_check.client.chebi_structure")
def test_check_row_never_raises_on_a_pubchem_outage(
    mock_chebi: MagicMock, mock_cas: MagicMock
) -> None:
    """check_row's contract is that nothing upstream can make it throw."""
    from structure_check.client import (
        CAS_UNREACHABLE,
        IDENTIFIER_UNREACHABLE,
        PubChemUnavailableError,
    )

    mock_cas.side_effect = PubChemUnavailableError("down")
    mock_chebi.return_value = (ETHANOL, "")

    result = check_row(cas="64-17-5", chebi_id="CHEBI:16236")

    assert result["id_cas_verdict"] == CAS_UNREACHABLE
    assert result["cas_status"] == IDENTIFIER_UNREACHABLE
    # Only the ID check was requested, so "id" is the only correct answer.
    assert result["unasked"] == "id"


@patch("structure_check.client.time.sleep")
@patch("structure_check.client.inchikeys_for_name")
def test_an_unreachable_whole_name_does_not_fall_through_to_tokens(
    mock_keys: MagicMock, _sleep: MagicMock
) -> None:
    """
    The token fallback's premise is that no whole-string form resolved.

    An unreachable whole-string query never established that, so falling through
    would compare the row against a free base the sheet never named — and the
    audit trail would read like an ordinary deliberate fallback.
    """
    from structure_check.client import PubChemUnavailableError, name_structure

    mock_keys.side_effect = PubChemUnavailableError("down")

    with pytest.raises(PubChemUnavailableError):
        name_structure("ABT_263_Navitoclax")

    # It stopped at the first whole-string form rather than trying the tokens.
    assert mock_keys.call_count == 1


def test_get_with_retry_still_collapses_both_silences_for_old_callers() -> None:
    """The wrapper must stay exactly as it was for every existing caller."""
    from chebi_lookup.client import (
        OUTCOME_NOT_FOUND,
        OUTCOME_UNREACHABLE,
        get_with_retry,
        get_with_retry_status,
    )

    with (
        patch("chebi_lookup.client.requests.get") as mock_get,
        patch("chebi_lookup.client.time.sleep"),
    ):
        missing = MagicMock()
        missing.status_code = 404
        mock_get.return_value = missing
        assert get_with_retry("https://example.invalid/x") is None
        assert (
            get_with_retry_status("https://example.invalid/x")[1] == OUTCOME_NOT_FOUND
        )

        down = MagicMock()
        down.status_code = 503
        mock_get.return_value = down
        assert get_with_retry("https://example.invalid/x") is None
        assert (
            get_with_retry_status("https://example.invalid/x")[1] == OUTCOME_UNREACHABLE
        )


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_parent_lookup_outage_never_degrades_the_run(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock, caplog
) -> None:
    """
    The one path where an outage cannot produce a false clean sheet.

    Refinement only ever downgrades a difference already found, so a failure there
    leaves a row flagged for a human — the safe direction. It is also by far the
    highest-volume caller, so counting it would exit 1 on a run whose only content
    is a genuine finding.
    """
    from structure_check.client import _parent_cache, parent_inchikey

    _parent_cache.clear()
    down = MagicMock()
    down.status_code = 503
    mock_get.return_value = down

    # Returns benignly rather than raising, and says why in the log.
    assert parent_inchikey(ETHANOL) == ""
    assert "unreachable" in caplog.text.lower()


def test_a_small_batch_is_not_quieter_than_the_same_row_in_single_mode() -> None:
    """
    Four rows against a dead ChEBI sit under every threshold: too few outage rows
    for `outages`, no `missing` values for `suspect_columns` (unreachable is not
    an answer), too few rows for the `nothing_compared` floor. Meanwhile single
    mode exits 1 on one such row. An outage is positive evidence that we failed
    to ask, so it does not need volume to be believed.
    """
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 4}),
        verdict_counts=Counter(),
        rows=4,
        cas=SideTally(resolved_values=4, resolved_rows=4),
        chebi=SideTally(outage_rows=4),
    )
    assert summary.outages == ()
    assert summary.suspect_columns == ()
    assert summary.nothing_compared
    assert summary.degraded


def test_lifting_that_floor_does_not_reopen_the_sparse_column_hole() -> None:
    """The 115 rows that compared fine keep nothing_compared false on the count."""
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 115, REVIEW_UNVERIFIED: 2}),
        cas=SideTally(resolved_values=117, resolved_rows=117),
        chebi=SideTally(outage_rows=2),
    )
    assert not summary.nothing_compared
    assert not summary.degraded


def test_a_tiny_spot_check_with_no_outage_still_stays_quiet() -> None:
    """Two compounds PubChem happens not to know is a spot check, not a breakage."""
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 2}),
        verdict_counts=Counter(),
        rows=2,
        cas=SideTally(missing_values=2),
    )
    assert not summary.nothing_compared
    assert not summary.degraded


@patch("structure_check.io.check_row")
def test_a_cache_hit_does_not_reset_the_fail_fast_counter(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A cache hit is not evidence the outage is over.

    Outages are never cached, so a hit means a row that succeeded *earlier* —
    possibly long before the upstream died. Resetting on one lets a sheet with
    repeated compounds alternate free hits with live rows and never trip the
    guard, spending full backoff on each one to reach a settled verdict.
    """
    from structure_check.client import (
        CAS_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )

    good = _row_result(
        review=REVIEW_OK,
        id_cas_verdict=MATCH,
        cas_status=IDENTIFIER_RESOLVED,
        chebi_status=IDENTIFIER_RESOLVED,
    )
    dead = _row_result(
        review=REVIEW_UNVERIFIED,
        id_cas_verdict=CAS_UNREACHABLE,
        cas_status=IDENTIFIER_UNREACHABLE,
    )

    def per_row(*, cas, chebi_id, name):
        return good if cas == "repeat" else dead

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    # The repeated compound resolves once in row 1, then alternates with dead
    # rows: every other row is a free cache hit.
    rows = [["repeat", "CHEBI:1"]]
    for n in range(12):
        rows.append([f"{n}-00-0", f"CHEBI:{n}"])
        rows.append(["repeat", "CHEBI:1"])
    _write_csv(src, ["CAS", "ChEBI ID"], rows)

    with pytest.raises(StructureCheckError, match="could not be reached"):
        check_file(src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID")
