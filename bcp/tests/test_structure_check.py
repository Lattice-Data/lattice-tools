"""Unit tests for structure_check client and io (mocked HTTP)."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from structure_check.client import (
    CAS_UNRESOLVED,
    CHEBI_NO_STRUCTURE,
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


@pytest.mark.parametrize(
    "reference, candidates",
    [("", [ETHANOL]), (ETHANOL, []), (ETHANOL, [""])],
)
def test_compare_without_both_sides_is_not_a_match(
    reference: str, candidates: list[str]
) -> None:
    assert compare_structures(reference, candidates) != MATCH


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
    result = empty_result()
    assert result["review"] == REVIEW_UNVERIFIED
    assert set(result) == set(OUTPUT_FIELDS_APPENDED)


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
        name = ("Ethanol", [ETHANOL])
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
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(name=("Alexidine", [ALEXIDINE]))
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
        name=("UCN-01", [UCN01]),
    )
    with cas_p, chebi_p, name_p, refine_p:
        result = check_row(cas="112953-11-4", chebi_id="CHEBI:221840", name="UCN-01")
    assert result["name_cas_verdict"] == STEREO_DIFFERS
    assert result["review"] == REVIEW_CHECK


def test_check_row_unresolved_name_is_not_a_finding() -> None:
    cas_p, chebi_p, name_p, refine_p = _patch_lookups(name=("", []))
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
        name=("tokens: Vorinostat|SAHA", [ETHANOL])
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
    query, keys = name_structure("PIK-75_HCl")
    assert query == "PIK-75_HCl"
    assert keys == [ETHANOL]
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
    query, keys = name_structure("ABT_263_Navitoclax")
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
    assert name_structure("Widget_9000") == ("", [])


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
    assert chebi_structure("CHEBI:16236") == ("", CHEBI_UNRESOLVED)


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
