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
    MAX_CONSECUTIVE_OUTAGE_ROWS,
    QuestionTally,
    RunSummary,
    SideTally,
    StructureCheckError,
    check_file,
    check_output_path,
    emit_single_row,
)


@pytest.fixture(autouse=True)
def _clear_parent_cache():
    """
    parent_inchikey memoises across a run, and the cache is a module global.

    Cleared around every test rather than by hand in the five that remember to:
    a test added later that exercises the real resolver would otherwise be
    served another test's entries, and pass or fail depending on ordering.
    """
    import structure_check.client as sc

    sc._parent_cache.clear()
    yield
    sc._parent_cache.clear()


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
    same way. cas_status stays not_checked and enters none of the run tallies."""
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
    # The column both the docs and the runtime WARNING send operators to has to
    # actually reach the file.
    assert "unasked" in rows[0]
    assert rows[0]["unasked"] == ""


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
    """
    A blank ChEBI cell among the rows used to hide a mostly-dead ChEBI.

    One row resolves in the middle, so the per-side fail-fast streak never
    reaches its limit and the run completes — the tally is what is under test
    here, not the abort.
    """
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_NOT_CHECKED,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
        NOT_CHECKED,
    )

    def per_row(*, cas, chebi_id, name, skip=()):
        if not chebi_id:
            status, verdict = IDENTIFIER_NOT_CHECKED, NOT_CHECKED
        elif chebi_id == "CHEBI:live":
            status, verdict = IDENTIFIER_RESOLVED, MATCH
        else:
            status, verdict = IDENTIFIER_UNREACHABLE, CHEBI_UNREACHABLE
        return _row_result(
            review=REVIEW_OK,
            id_cas_verdict=verdict,
            name_cas_verdict=MATCH,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=status,
            name_status=IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    rows = [[f"C{n}", f"{n}-00-0", f"CHEBI:{n}"] for n in range(3)]
    rows.append(["Clive", "9-00-0", "CHEBI:live"])  # breaks the streak
    rows += [[f"D{n}", f"1{n}-00-0", f"CHEBI:1{n}"] for n in range(3)]
    rows.append(["Cblank", "5-00-0", ""])  # the blank that used to hide everything
    _write_csv(src, ["Name", "CAS", "ChEBI ID"], rows)

    summary = check_file(
        src,
        tmp_path / "out.csv",
        cas_column="CAS",
        chebi_column="ChEBI ID",
        name_column="Name",
    )

    # The blank contributed to neither count, and the one live row is not enough
    # to outweigh six unanswered ones.
    assert summary.chebi.outage_rows == 6
    assert summary.chebi.resolved_rows == 1
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

    def per_row(*, cas, chebi_id, name, skip=()):
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
def test_check_file_stops_querying_a_dead_upstream_but_finishes_the_run(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A dead upstream stops being asked; the run still completes.

    Aborting would throw away every remaining row, including the columns that
    are still answering. Skipped rows are marked unreachable — which is what
    they are — so the verdict is unchanged and the CSV is complete.
    """
    from structure_check.client import CAS_UNREACHABLE, IDENTIFIER_UNREACHABLE

    seen_skips = []

    def per_row(*, cas, chebi_id, name, skip=()):
        seen_skips.append(frozenset(skip))
        return _row_result(
            review=REVIEW_UNVERIFIED,
            id_cas_verdict=CAS_UNREACHABLE,
            cas_status=IDENTIFIER_UNREACHABLE,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    _write_csv(
        src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", f"CHEBI:{n}"] for n in range(40)]
    )
    out = tmp_path / "out.csv"

    summary = check_file(src, out, cas_column="CAS", chebi_column="ChEBI ID")

    # Every row is present and accounted for, not just the ones before the trip.
    assert summary.rows == 40
    assert len(list(csv.DictReader(out.open(encoding="utf-8")))) == 40
    # The CAS side stopped being queried once it had failed five times running.
    assert any("cas" in sk for sk in seen_skips)
    # And the run is still reported as untrustworthy.
    assert summary.outages == ("CAS",)
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_a_tripped_side_is_retried_rather_than_written_off(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    Being tripped is a guess that the upstream is down, so it gets re-tested.

    Without a probe, five scattered blips would poison every remaining row with
    an `unreachable` no request was ever made for — an answer manufactured from
    a silence, which is the one thing this module refuses to do.
    """
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )
    from structure_check.io import OUTAGE_PROBE_INTERVAL

    attempts = []

    def per_row(*, cas, chebi_id, name, skip=()):
        if "chebi" in skip:
            return _row_result(
                review=REVIEW_UNVERIFIED,
                id_cas_verdict=CHEBI_UNREACHABLE,
                cas_status=IDENTIFIER_RESOLVED,
                chebi_status=IDENTIFIER_UNREACHABLE,
            )
        attempts.append(chebi_id)
        # ChEBI is down for the first five attempts, then comes back.
        down = len(attempts) <= 5
        return _row_result(
            review=REVIEW_UNVERIFIED if down else REVIEW_OK,
            id_cas_verdict=CHEBI_UNREACHABLE if down else MATCH,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=IDENTIFIER_UNREACHABLE if down else IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    # A long tail after recovery, so the rows that were answered outnumber the
    # ones lost while it was down and the run is not condemned for the outage.
    rows = 5 + OUTAGE_PROBE_INTERVAL + 40
    _write_csv(
        src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", f"CHEBI:{n}"] for n in range(rows)]
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    # It tripped after 5, waited, probed, found ChEBI back, and resumed.
    assert len(attempts) > 5
    assert summary.chebi.resolved_rows > 0
    # Having recovered, the run is not condemned.
    assert summary.outages == ()


@patch("structure_check.io.check_row")
def test_a_cache_hit_does_not_clear_an_outage_streak(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    A cache hit is not evidence the upstream recovered.

    Outages are never cached, so a hit means a row that succeeded *earlier*,
    possibly long before the upstream died. Folding one in would let a sheet with
    repeated compounds alternate free hits with live rows and never trip at all.
    """
    from structure_check.client import (
        CAS_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )

    saw_cas_skipped = []

    def per_row(*, cas, chebi_id, name, skip=()):
        saw_cas_skipped.append("cas" in skip)
        if cas == "repeat":
            return _row_result(
                review=REVIEW_OK,
                id_cas_verdict=MATCH,
                cas_status=IDENTIFIER_RESOLVED,
                chebi_status=IDENTIFIER_RESOLVED,
            )
        return _row_result(
            review=REVIEW_UNVERIFIED,
            id_cas_verdict=CAS_UNREACHABLE,
            cas_status=IDENTIFIER_UNREACHABLE,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    # The repeated compound resolves once in row 1, then alternates with dead
    # rows: every other row is a free cache hit.
    rows = [["repeat", "CHEBI:1"]]
    for n in range(12):
        rows.append([f"{n}-00-0", f"CHEBI:{n}"])
        rows.append(["repeat", "CHEBI:1"])
    _write_csv(src, ["CAS", "ChEBI ID"], rows)

    check_file(src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID")

    # The streak survived the interleaved cache hits and the side was given up
    # on. (Note `outages` does not fire here: the cached rows count as answered,
    # and they outnumber the dead ones — the breaker and the run verdict ask
    # different questions, and this test is about the breaker.)
    assert any(saw_cas_skipped), "the CAS side was never given up on"


@patch("structure_check.io.check_row")
def test_a_dead_side_is_dropped_while_the_live_one_finishes_the_sheet(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    The case a whole-row rule cannot see, and the reason this is not an abort.

    With ChEBI dead run-wide and a live name column whose names match, every row
    reads `ok`, so a row-level streak resets on every row. Tracking it per side
    catches it — and dropping only the dead side means the name check still
    completes for the whole sheet instead of being thrown away with it.
    """
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )

    chebi_attempts = []

    def per_row(*, cas, chebi_id, name, skip=()):
        if "chebi" not in skip:
            chebi_attempts.append(chebi_id)
        return _row_result(
            review=REVIEW_OK,  # the name side is fine, so the row reads as ok
            id_cas_verdict=CHEBI_UNREACHABLE,
            name_cas_verdict=MATCH,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=IDENTIFIER_UNREACHABLE,
            name_status=IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["Name", "CAS", "ChEBI ID"],
        [[f"C{n}", f"{n}-00-0", f"CHEBI:{n}"] for n in range(40)],
    )

    summary = check_file(
        src,
        tmp_path / "out.csv",
        cas_column="CAS",
        chebi_column="ChEBI ID",
        name_column="Name",
    )

    # ChEBI stopped being asked, far short of once per row...
    assert len(chebi_attempts) < 40
    # ...while the name column was checked on every single row.
    assert summary.name.resolved_rows == 40
    assert summary.rows == 40
    # And the ID question is still correctly reported as unanswered.
    assert summary.outages == ("ChEBI ID",)
    assert summary.degraded


@patch("structure_check.io.check_row")
def test_a_blank_cell_neither_extends_nor_breaks_an_outage_streak(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """A sparsely mapped column says nothing about whether the upstream is up."""
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_NOT_CHECKED,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        NOT_CHECKED,
    )

    def per_row(*, cas, chebi_id, name, skip=()):
        unreachable = bool(chebi_id)
        return _row_result(
            review=REVIEW_UNVERIFIED,
            id_cas_verdict=CHEBI_UNREACHABLE if unreachable else NOT_CHECKED,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=(
                IDENTIFIER_UNREACHABLE if unreachable else IDENTIFIER_NOT_CHECKED
            ),
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    # Only 4 ChEBI IDs, spread across 20 rows of blanks. Four unreachable rows is
    # under the limit, and the 16 blanks must not be read as recoveries either.
    rows = [[f"{n}-00-0", f"CHEBI:{n}" if n % 5 == 0 else ""] for n in range(20)]
    _write_csv(src, ["CAS", "ChEBI ID"], rows)

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )
    assert summary.chebi.outage_rows == 4


@patch("structure_check.client.cas_structure")
@patch("structure_check.client.chebi_structure")
@patch("structure_check.client.name_structure")
def test_an_unreachable_cas_leaves_both_questions_unasked(
    mock_name: MagicMock, mock_chebi: MagicMock, mock_cas: MagicMock
) -> None:
    """
    The CAS is the pivot, so losing it costs both comparisons at once.

    `unasked` is the whole record of which question went unanswered, so the
    branch that says "neither of them" is the one most worth pinning.
    """
    from structure_check.client import (
        CAS_UNREACHABLE,
        IDENTIFIER_UNREACHABLE,
        PubChemUnavailableError,
    )

    mock_cas.side_effect = PubChemUnavailableError("down")
    mock_chebi.return_value = (ETHANOL, "")
    mock_name.return_value = ("Ethanol", [ETHANOL], 1)

    result = check_row(cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol")

    assert result["id_cas_verdict"] == CAS_UNREACHABLE
    assert result["name_cas_verdict"] == CAS_UNREACHABLE
    assert result["unasked"] == "both"
    assert result["cas_status"] == IDENTIFIER_UNREACHABLE
    assert result["review"] == REVIEW_UNVERIFIED


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_client_error_is_an_answer_about_the_request_not_the_network(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """
    PUG REST answers 400 for input it cannot parse — which is what a column of
    free-text notes produces. Reported as unreachable it would blame the network
    for a bad cell, and suppress the wrong-column diagnosis, since an outage
    disqualifies suspect_columns.
    """
    from chebi_lookup.client import (
        MALFORMED_REQUEST_CODES,
        OUTCOME_NOT_FOUND,
        get_with_retry_status,
    )

    for code in sorted(MALFORMED_REQUEST_CODES):
        mock_get.reset_mock()
        mock_get.return_value = MagicMock(status_code=code)
        _resp, outcome = get_with_retry_status("https://example.invalid/x")
        assert outcome == OUTCOME_NOT_FOUND, code
        assert mock_get.call_count == 1, f"{code} must not be retried"


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_blocked_client_is_an_outage_not_a_missing_compound(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """
    401/403 say nothing about the compound — they say we were not let in.

    Filed with the malformed-request codes, a run-wide block would report every
    good CAS as "no such compound", leave `outages` empty, and trip
    `suspect_columns` — telling the operator to check their column while the
    real cause was access. That is the misdiagnosis the three signals exist to
    keep apart, and the one direction a re-run can fix.
    """
    from chebi_lookup.client import OUTCOME_UNREACHABLE, get_with_retry_status

    for code in (401, 403, 407):
        mock_get.reset_mock()
        mock_get.return_value = MagicMock(status_code=code)
        _resp, outcome = get_with_retry_status("https://example.invalid/x")
        assert outcome == OUTCOME_UNREACHABLE, code


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_server_error_is_still_retried_and_still_an_outage(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    from chebi_lookup.client import OUTCOME_UNREACHABLE, get_with_retry_status

    mock_get.return_value = MagicMock(status_code=500)
    _resp, outcome = get_with_retry_status("https://example.invalid/x")
    assert outcome == OUTCOME_UNREACHABLE
    assert mock_get.call_count > 1


def test_check_file_rejects_the_same_column_for_chebi_and_name(tmp_path: Path) -> None:
    """Every row would come back name_unresolved, reading like a data problem."""
    src = tmp_path / "in.csv"
    _write_csv(src, ["CAS", "ChEBI ID"], [["64-17-5", "CHEBI:16236"]])
    with pytest.raises(StructureCheckError, match="both --chebi-column"):
        check_file(
            src,
            tmp_path / "out.csv",
            cas_column="CAS",
            chebi_column="ChEBI ID",
            name_column="ChEBI ID",
        )


def test_an_answered_404_counts_as_an_answer_against_an_outage() -> None:
    """
    "No such compound" is an answer, so it belongs on the answered side of the
    majority. Counting only rows that *resolved* let five throttled rows outweigh
    a hundred honest 404s — a column of internal codes PubChem does not carry.
    """
    summary = RunSummary(
        review_counts=Counter({REVIEW_OK: 3, REVIEW_UNVERIFIED: 117}),
        verdict_counts=Counter(),
        rows=120,
        name=SideTally(
            resolved_values=3,
            missing_values=112,
            resolved_rows=3,
            missing_rows=112,
            outage_rows=5,
        ),
    )
    assert summary.name.answered_rows == 115
    assert summary.outages == ()
    assert not summary.degraded


def test_a_wrong_column_is_still_named_when_a_few_rows_were_throttled() -> None:
    """
    The costly version of the same bug: against a wrong column nothing resolves
    by construction, so without counting answered-404s the outage branch fired
    and suppressed the diagnosis that would have named the real problem.
    """
    summary = RunSummary(
        review_counts=Counter({REVIEW_UNVERIFIED: 120}),
        verdict_counts=Counter(),
        rows=120,
        name=SideTally(missing_values=116, missing_rows=116, outage_rows=5),
    )
    assert summary.outages == ()
    assert summary.suspect_columns == ("name",)
    assert summary.degraded


def test_a_real_outage_still_outweighs_the_rows_it_answered() -> None:
    summary = _summary(
        chebi=SideTally(resolved_rows=3, missing_rows=4, outage_rows=40),
    )
    assert summary.outages == ("ChEBI ID",)


@patch("structure_check.io.check_row")
def test_check_file_counts_answered_and_unanswered_rows_apart(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """End to end: the three row counters must add up to the rows attempted."""
    from structure_check.client import (
        CHEBI_UNRESOLVED,
        IDENTIFIER_MISSING,
        IDENTIFIER_RESOLVED,
        MATCH,
    )

    def per_row(*, cas, chebi_id, name, skip=()):
        hit = chebi_id.endswith(("0", "1", "2"))
        return _row_result(
            review=REVIEW_OK if hit else REVIEW_UNVERIFIED,
            id_cas_verdict=MATCH if hit else CHEBI_UNRESOLVED,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=IDENTIFIER_RESOLVED if hit else IDENTIFIER_MISSING,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    _write_csv(
        src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", f"CHEBI:{n}"] for n in range(10)]
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    assert summary.chebi.resolved_rows == 3
    assert summary.chebi.missing_rows == 7
    assert summary.chebi.answered_rows == 10
    assert summary.chebi.outage_rows == 0


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_a_200_carrying_unrelated_json_is_an_outage_not_a_missing_cas(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """
    A proxy or maintenance page answering `{"Fault": ...}` parses fine as JSON.

    Defaulting the missing IdentifierList to [] would report NOT_FOUND, and five
    such values would then blame the CAS column for what the network did.
    """
    from chebi_lookup.client import OUTCOME_UNREACHABLE, cas_to_cid_status

    fault = MagicMock(status_code=200)
    fault.json.return_value = {"Fault": {"Code": "PUGREST.ServerBusy"}}
    mock_get.return_value = fault

    assert cas_to_cid_status("64-17-5") == (None, OUTCOME_UNREACHABLE)


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_an_empty_cid_list_is_still_a_real_answer(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    from chebi_lookup.client import OUTCOME_NOT_FOUND, cas_to_cid_status

    empty = MagicMock(status_code=200)
    empty.json.return_value = {"IdentifierList": {"CID": []}}
    mock_get.return_value = empty

    assert cas_to_cid_status("00-00-0") == (None, OUTCOME_NOT_FOUND)


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_cas_to_cid_status_quotes_the_sheet_cell_for_every_caller(
    mock_get: MagicMock, _sleep: MagicMock
) -> None:
    """Quoted inside, so lookup_cas is covered by the same guard as this module."""
    from chebi_lookup.client import cas_to_cid_status

    mock_get.return_value = MagicMock(status_code=404)
    cas_to_cid_status("64-17-5/../../evil?x=1")

    segment = mock_get.call_args[0][0].split("/compound/name/")[1].split("/cids")[0]
    assert "/" not in segment
    assert "?" not in segment


@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_no_backoff_sleep_after_the_final_attempt(
    mock_get: MagicMock, mock_sleep: MagicMock
) -> None:
    """
    The backoff spaces out the *next* try, and after the last one there isn't one.

    It is 8 of the ~14s an exhausted request costs, paid on every row of a run
    against a dead upstream — which is what MAX_CONSECUTIVE_OUTAGE_ROWS is sized
    against. chebi_terms already guards this.
    """
    from chebi_lookup.client import MAX_RETRIES, get_with_retry_status

    mock_get.return_value = MagicMock(status_code=503)
    get_with_retry_status("https://example.invalid/x")

    assert mock_get.call_count == MAX_RETRIES
    assert mock_sleep.call_count == MAX_RETRIES - 1


@pytest.mark.parametrize(
    "raw", ["Diclofenac_2Na", "Compound_Na2", "Foo_NaCl", "Bar_2K", "Baz_TsOH"]
)
def test_a_counterion_written_as_a_formula_also_disables_the_fallback(
    raw: str,
) -> None:
    """
    The guard covered counterions spelled as words but not as formulae.

    'Diclofenac_2Na' split into tokens, so a cell that failed to resolve whole
    would have queried the free base against the sodium salt's CAS — the
    difference of the tool's own making this guard exists to prevent, and which
    the docs claim it prevents absolutely.
    """
    assert name_candidates(raw)[1] == []


def test_the_formula_guard_leaves_genuine_alias_pairs_alone() -> None:
    """The guard must not swallow the case the fallback exists for."""
    for raw in ("Vorinostat_SAHA", "ABT_263_Navitoclax", "Actinomycin_D_Dactinomycin"):
        assert name_candidates(raw)[1], raw


@pytest.mark.parametrize(
    ("chebi_down", "name_down", "expected"),
    [(True, False, "id"), (False, True, "name"), (True, True, "both")],
)
def test_unasked_names_whichever_side_went_unanswered(
    chebi_down: bool, name_down: bool, expected: str
) -> None:
    """`unasked` is a documented output column; all three of its values are real."""
    from structure_check.client import CHEBI_UNREACHABLE, PubChemUnavailableError

    chebi = ("", CHEBI_UNREACHABLE) if chebi_down else (ETHANOL, "")
    name_p = patch("structure_check.client.name_structure")
    with (
        patch("structure_check.client.cas_structure", return_value=(ETHANOL, "Eth")),
        patch("structure_check.client.chebi_structure", return_value=chebi),
        name_p as mock_name,
    ):
        if name_down:
            mock_name.side_effect = PubChemUnavailableError("down")
        else:
            mock_name.return_value = ("Ethanol", [ETHANOL], 1)
        result = check_row(cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol")

    assert result["unasked"] == expected


# --------------------------------------------------------------------------
# check_row's skip= contract, exercised directly
#
# Every breaker test patches structure_check.io.check_row with a fake, so the
# real branches were never executed. skip= is a public keyword argument whose
# behaviour the docs describe ("marked unreachable without a request").
# --------------------------------------------------------------------------


@patch("structure_check.client.name_structure")
@patch("structure_check.client.chebi_structure")
@patch("structure_check.client.cas_structure")
def test_skipping_the_cas_side_costs_no_request_at_all(
    mock_cas: MagicMock, mock_chebi: MagicMock, mock_name: MagicMock
) -> None:
    """The CAS is the pivot, so skipping it settles both verdicts."""
    from structure_check.client import CAS_UNREACHABLE, IDENTIFIER_UNREACHABLE

    result = check_row(
        cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol", skip={"cas"}
    )

    mock_cas.assert_not_called()
    mock_chebi.assert_not_called()
    mock_name.assert_not_called()
    assert result["cas_status"] == IDENTIFIER_UNREACHABLE
    assert result["id_cas_verdict"] == CAS_UNREACHABLE
    assert result["name_cas_verdict"] == CAS_UNREACHABLE
    assert result["unasked"] == "both"


@patch("structure_check.client.name_structure")
@patch("structure_check.client.chebi_structure")
@patch("structure_check.client.cas_structure")
def test_skipping_the_chebi_side_leaves_the_name_check_running(
    mock_cas: MagicMock, mock_chebi: MagicMock, mock_name: MagicMock
) -> None:
    from structure_check.client import CHEBI_UNREACHABLE, IDENTIFIER_UNREACHABLE

    mock_cas.return_value = (ETHANOL, "Ethanol")
    mock_name.return_value = ("Ethanol", [ETHANOL], 1)

    result = check_row(
        cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol", skip={"chebi"}
    )

    mock_chebi.assert_not_called()
    mock_name.assert_called_once()
    assert result["chebi_status"] == IDENTIFIER_UNREACHABLE
    assert result["id_cas_verdict"] == CHEBI_UNREACHABLE
    assert result["name_cas_verdict"] == MATCH
    assert result["unasked"] == "id"


@patch("structure_check.client.name_structure")
@patch("structure_check.client.chebi_structure")
@patch("structure_check.client.cas_structure")
def test_skipping_the_name_side_leaves_the_id_check_running(
    mock_cas: MagicMock, mock_chebi: MagicMock, mock_name: MagicMock
) -> None:
    from structure_check.client import IDENTIFIER_UNREACHABLE, NAME_UNREACHABLE

    mock_cas.return_value = (ETHANOL, "Ethanol")
    mock_chebi.return_value = (ETHANOL, "")

    result = check_row(
        cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol", skip={"name"}
    )

    mock_name.assert_not_called()
    mock_chebi.assert_called_once()
    assert result["name_status"] == IDENTIFIER_UNREACHABLE
    assert result["name_cas_verdict"] == NAME_UNREACHABLE
    assert result["id_cas_verdict"] == MATCH
    assert result["unasked"] == "name"


@patch("structure_check.io.check_row")
def test_breaker_skipped_rows_do_not_vote_in_the_outage_majority(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """
    The breaker must not manufacture the evidence that justifies it.

    Five real failures trip it; the twenty rows it then skips are inferred from
    those same five. Counting them turned 5-vs-25 into 25-vs-5 and flipped a run
    that recovered into exit 1.
    """
    from structure_check.client import (
        CHEBI_UNREACHABLE,
        IDENTIFIER_RESOLVED,
        IDENTIFIER_UNREACHABLE,
        MATCH,
    )

    live = [0]

    def per_row(*, cas, chebi_id, name, skip=()):
        if "chebi" in skip:
            return _row_result(
                review=REVIEW_UNVERIFIED,
                id_cas_verdict=CHEBI_UNREACHABLE,
                cas_status=IDENTIFIER_RESOLVED,
                chebi_status=IDENTIFIER_UNREACHABLE,
            )
        live[0] += 1
        down = live[0] <= MAX_CONSECUTIVE_OUTAGE_ROWS
        return _row_result(
            review=REVIEW_UNVERIFIED if down else REVIEW_OK,
            id_cas_verdict=CHEBI_UNREACHABLE if down else MATCH,
            cas_status=IDENTIFIER_RESOLVED,
            chebi_status=IDENTIFIER_UNREACHABLE if down else IDENTIFIER_RESOLVED,
        )

    mock_check.side_effect = per_row
    src = tmp_path / "in.csv"
    _write_csv(
        src, ["CAS", "ChEBI ID"], [[f"{n}-00-0", f"CHEBI:{n}"] for n in range(30)]
    )

    summary = check_file(
        src, tmp_path / "out.csv", cas_column="CAS", chebi_column="ChEBI ID"
    )

    # Only the rows actually asked about count as evidence.
    assert summary.chebi.outage_rows == MAX_CONSECUTIVE_OUTAGE_ROWS
    assert summary.outages == ()
    assert not summary.degraded


# --------------------------------------------------------------------------
# a question asked on every row and compared on none
# --------------------------------------------------------------------------


def test_a_question_answered_on_no_row_is_degraded_even_when_rows_read_ok() -> None:
    """
    The shape none of the other three signals can see.

    A ChEBI column of class terms: ChEBI holds every record, so nothing is
    unresolved and nothing is unreachable; the healthy name column makes every
    row `ok`, so nothing_compared is false. The ID question — the headline
    column — was answered on no row at all.
    """
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 117}),
        chebi=SideTally(resolved_values=117, resolved_rows=117),
        name=SideTally(resolved_values=117, resolved_rows=117),
        id_question=QuestionTally(asked=117, compared=0),
        name_question=QuestionTally(asked=117, compared=117),
    )
    assert summary.outages == ()
    assert summary.suspect_columns == ()
    assert not summary.nothing_compared
    assert summary.dead_questions == ("ChEBI ID vs CAS",)
    assert summary.degraded


def test_a_question_asked_on_two_rows_is_a_finding_not_a_dead_question() -> None:
    """The retired-ID case again: two bad IDs are the product, not a broken run."""
    summary = _summary(
        review_counts=Counter({REVIEW_OK: 115, REVIEW_UNVERIFIED: 2}),
        id_question=QuestionTally(asked=2, compared=0),
        name_question=QuestionTally(asked=117, compared=117),
    )
    assert summary.dead_questions == ()
    assert not summary.degraded


def test_one_successful_comparison_keeps_a_question_alive() -> None:
    summary = _summary(
        id_question=QuestionTally(asked=117, compared=1),
    )
    assert summary.dead_questions == ()


@patch("structure_check.io.check_row")
def test_check_file_flags_a_chebi_column_that_never_compares(
    mock_check: MagicMock, tmp_path: Path
) -> None:
    """End to end: ChEBI answers every row, and no comparison ever happens."""
    from structure_check.client import CHEBI_NO_STRUCTURE, IDENTIFIER_RESOLVED, MATCH

    mock_check.return_value = _row_result(
        review=REVIEW_OK,
        id_cas_verdict=CHEBI_NO_STRUCTURE,
        name_cas_verdict=MATCH,
        cas_status=IDENTIFIER_RESOLVED,
        chebi_status=IDENTIFIER_RESOLVED,
        name_status=IDENTIFIER_RESOLVED,
    )
    src = tmp_path / "in.csv"
    _write_csv(
        src,
        ["Name", "CAS", "ChEBI ID"],
        [[f"C{n}", f"{n}-00-0", f"CHEBI:{n}"] for n in range(8)],
    )

    summary = check_file(
        src,
        tmp_path / "out.csv",
        cas_column="CAS",
        chebi_column="ChEBI ID",
        name_column="Name",
    )

    assert summary.review_counts[REVIEW_OK] == 8
    assert summary.id_question == QuestionTally(asked=8, compared=0)
    assert summary.name_question.compared == 8
    assert summary.dead_questions == ("ChEBI ID vs CAS",)
    assert summary.degraded


def test_an_upstream_that_dies_midway_still_degrades_the_run() -> None:
    """
    Skipped rows count when the side never came back.

    outage_rows stops growing once the breaker trips, so without this a side
    that died at row 100 of 500 shows ~24 unanswered against 99 answered and
    exits 0 — the majority could only ever fire on a side dead from the start.
    """
    summary = _summary(
        rows=500,
        review_counts=Counter({REVIEW_OK: 500}),
        chebi=SideTally(
            resolved_values=100, resolved_rows=100, outage_rows=24, skipped_rows=377
        ),
        abandoned=frozenset({"chebi"}),
    )
    assert summary.outages == ("ChEBI ID",)
    assert summary.degraded


def test_a_side_that_came_back_is_not_judged_on_its_skipped_rows() -> None:
    """The same tally, minus the abandonment, must stay quiet."""
    summary = _summary(
        rows=500,
        review_counts=Counter({REVIEW_OK: 500}),
        chebi=SideTally(
            resolved_values=100, resolved_rows=100, outage_rows=24, skipped_rows=377
        ),
        abandoned=frozenset(),
    )
    assert summary.outages == ()
    assert not summary.degraded


@patch("structure_check.client.time.sleep")
@patch("chebi_lookup.client.time.sleep")
@patch("chebi_lookup.client.requests.get")
def test_parent_inchikey_does_not_cache_a_garbled_two_hundred(
    mock_get: MagicMock, _s1: MagicMock, _s2: MagicMock
) -> None:
    """
    A maintenance page is a silence, not "this structure has no parent".

    _first_cid swallowed an unparseable body into None, so the empty result was
    memoised for the rest of the run — against the cache's own stated invariant,
    and everywhere else in the module a 200 that is not this endpoint's JSON is
    treated as unreachable.
    """
    import structure_check.client as sc

    html = MagicMock(status_code=200)
    html.json.side_effect = ValueError("not json")
    mock_get.return_value = html

    assert sc.parent_inchikey(DOXORUBICIN_HCL) == ""
    assert DOXORUBICIN_HCL not in sc._parent_cache
