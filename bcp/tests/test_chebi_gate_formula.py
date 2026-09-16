"""Registry formula against drawn stoichiometry, and the CAS table it comes from.

EXT-04 is the only stoichiometry check whose evidence does not come from the step
that generated the SDFs, so it is the only one that can see a missing counterion.
Most of these tests are a defect that produced a clean report: a status of "agree"
where nothing had been confirmed is the exact failure mode to guard.
"""

from __future__ import annotations

import csv
from fractions import Fraction
from pathlib import Path

import pytest

from chebi_gate import casreg, formula

REGISTRY_ROW = {name: "" for name in casreg.REGISTRY_COLUMNS}


def registry(tmp_path: Path, *rows: dict) -> Path:
    path = tmp_path / "cas_registry.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(casreg.REGISTRY_COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow({**REGISTRY_ROW, **row})
    return path


# ----------------------------------------------------------------- parsing


def test_a_simple_formula_parses_to_element_counts():
    assert formula.parse("C2H6O") == {"C": 2, "H": 6, "O": 1}


def test_an_element_with_no_count_means_one():
    assert formula.parse("HCl") == {"H": 1, "Cl": 1}


def test_a_dotted_formula_sums_its_components():
    assert formula.parse("C2H7N.2HCl") == {"C": 2, "H": 9, "N": 1, "Cl": 2}


def test_html_markup_from_a_scifinder_export_is_stripped():
    assert formula.parse("C<sub>2</sub>H<sub>6</sub>O") == {"C": 2, "H": 6, "O": 1}


def test_a_fractional_multiplier_stays_exact_until_reduction():
    parsed = formula.parse("C12H18N2O.3/2C4H4O4")
    assert parsed["C"] == Fraction(18)
    assert parsed["O"] == Fraction(7)


def test_a_decimal_multiplier_means_the_same_as_the_fraction_it_spells():
    """`split(".")` cut the decimal point and made the tail its own multiplier.

    "C19H23N.1.5C4H4O4" split to ["C19H23N", "1", "5C4H4O4"]: the "1" was a
    multiplier with no body and was skipped, and the "5" multiplied the fumarate,
    so a sesquifumarate parsed as C39H43NO20 -- silently, and confidently enough
    to be compared against a correct drawing and called a disagreement.
    """
    assert formula.parse("C19H23N.1.5C4H4O4") == formula.parse("C19H23N.3/2C4H4O4")
    assert formula.reduce_ratio(formula.parse("C19H23N.1.5C4H4O4")) == {
        "C": 25,
        "H": 29,
        "N": 1,
        "O": 6,
    }


def test_a_dot_before_a_hydrate_multiplier_is_still_a_separator():
    """The half of the fix that can be got wrong: C4H10O2.2H2O is not a decimal."""
    assert formula.reduce_ratio(formula.parse("C4H10O2.2H2O")) == {
        "C": 2,
        "H": 7,
        "O": 2,
    }


def test_a_zero_denominator_multiplier_is_unparseable_not_an_exception():
    """Fraction raises ZeroDivisionError, which left parse and aborted the run.

    One malformed cell in a hand-parsed table took down compare, ext04 and
    run_external_checks with it, so every other record went unjudged too.
    """
    assert formula.parse("C19H23N.1/0C4H4O4") is None
    assert formula.compare("C19H23N.1/0C4H4O4", "C23H27NO4").status == (
        formula.UNPARSEABLE
    )


def test_deuterium_and_tritium_are_counted_as_hydrogen():
    """Decided, not defaulted: RDKit's CalcMolFormula writes CH4O for CD3OH.

    The drawn side can never say D, so counting it as its own element would make
    every deuterated registry row DISAGREE with no edit that could clear it. The
    alternative -- leaving it unparseable -- loses a stoichiometry check over a
    formula that is perfectly well formed.
    """
    assert formula.parse("C2D6O") == formula.parse("C2H6O")
    assert formula.parse("C6H5T") == formula.parse("C6H6")


def test_an_element_whose_symbol_starts_with_d_or_t_is_still_that_element():
    assert formula.parse("Dy2O3") == {"Dy": Fraction(2), "O": Fraction(3)}
    assert formula.parse("TeO2") == {"Te": Fraction(1), "O": Fraction(2)}


def test_an_isotope_agreement_says_it_did_not_check_the_labelling():
    """Otherwise "formula agrees" reads as confirming a label nobody compared."""
    verdict = formula.compare("C2D6O", "C2H6O")
    assert verdict.status == formula.AGREE
    assert "does not check it" in verdict.note

    plain = formula.compare("C2H6O", "C2H6O")
    assert "does not check it" not in plain.note


def test_one_trailing_charge_is_dropped():
    """RDKit writes the net charge onto a whole-molecule formula."""
    assert formula.parse("C5H11NO5P-") == formula.parse("C5H11NO5P")
    assert formula.parse("C6H5O7-3") == formula.parse("C6H5O7")


def test_an_empty_formula_is_none():
    assert formula.parse("") is None
    assert formula.parse(None) is None


@pytest.mark.parametrize("text", ["WAG994", "Xx12", "C2H6Q"])
def test_a_token_that_is_not_an_element_makes_the_formula_unparseable(text):
    """Handoff case 3: a vendor code that lost its space parses as a formula.

    "WAG 994" de-spaces to "WAG994", and the element pattern [A-Z][a-z]? happily
    reads W, A, G994 -- a confident, meaningless formula that then overwrote the
    real one. Symbols are checked against the periodic table instead.
    """
    assert formula.parse(text) is None


def test_leftover_characters_make_the_formula_unparseable():
    """Partly read is worse than not read: it compares as though it were complete."""
    assert formula.parse("C2H6O(ish)") is None


# ---------------------------------------------------------------- reduction


def test_reduction_gives_the_smallest_whole_number_ratio():
    assert formula.reduce_ratio(formula.parse("C36H48N4O14")) == {
        "C": 18,
        "H": 24,
        "N": 2,
        "O": 7,
    }


def test_a_sesquifumarate_written_one_to_one_point_five_equals_one_drawn_two_to_three():
    """Handoff case 6, the whole reason formulae are compared as ratios."""
    registered = formula.reduce_ratio(formula.parse("C12H18N2O.3/2C4H4O4"))
    drawn = formula.reduce_ratio(formula.parse("C36H48N4O14"))
    assert registered == drawn


def test_reduction_of_nothing_is_none():
    assert formula.reduce_ratio(None) is None
    assert formula.reduce_ratio(formula.parse("")) is None


# --------------------------------------------------------- the five statuses


def test_identical_formulae_agree():
    result = formula.compare("C2H8ClN", "C2H8ClN")
    assert result.status == formula.AGREE
    assert result.confirms is True
    assert result.checked is True


def test_a_reduced_ratio_match_agrees():
    result = formula.compare("C12H18N2O.3/2C4H4O4", "C36H48N4O14")
    assert result.status == formula.AGREE


def test_a_hydrogen_only_difference_is_a_drawing_convention():
    """Handoff case 7: a neutral acid registered against an ionic drawing.

    10 of 290 records on the reference batch. Reporting these as defects would
    bury the three real disagreements.
    """
    result = formula.compare("C2H8ClN", "C2H7ClN")
    assert result.status == formula.CONVENTION
    assert result.diff == {"H": 1}
    assert result.confirms is True
    assert "neutral-acid versus ionic" in result.note


def test_a_formula_with_no_components_is_compared_as_it_stands():
    """Dividing out the gcd compared proportions, not composition.

    C6H12O6 and C2H4O2 both reduce to {C:1, H:2, O:1}, so glucose against acetic
    acid came back "agrees ... (same ratio)". Nothing else would catch it: a CAS
    source cannot refute, and the only source that can is PubChem, which is
    circular -- so the one independent stoichiometry check went quiet on a
    molecule three times the size it should be.
    """
    assert formula.compare("C6H12O6", "C2H4O2").status == formula.DISAGREE
    assert formula.compare("C4H4O4", "C2H2O2").status == formula.DISAGREE
    assert formula.compare("C2H6O", "C2H6O").status == formula.AGREE


def test_a_formula_that_declares_components_still_reduces():
    """The sesquifumarate case the reduction exists for."""
    verdict = formula.compare("C12H18N2O.3/2C4H4O4", "C36H48N4O14")
    assert verdict.status == formula.AGREE
    assert formula.declares_components("C12H18N2O.3/2C4H4O4")
    assert formula.declares_components("2C4H4O4")
    assert not formula.declares_components("C6H12O6")


def test_a_difference_in_anything_but_hydrogen_disagrees():
    result = formula.compare("C6H15N3S.2ClH", "C6H16ClN3S")
    assert result.status == formula.DISAGREE
    assert result.confirms is False
    assert "Cl+1" in result.note


def test_the_difference_sign_says_which_side_has_more():
    result = formula.compare("C6H15N3S.2ClH", "C6H16ClN3S")
    assert result.diff["Cl"] == 1, "positive means the registry has more"


def test_a_missing_registry_row_is_reported_as_unchecked_not_as_agreement():
    """Handoff section 1: a record with no row must say it could not be checked."""
    result = formula.compare(None, "C2H8ClN")
    assert result.status == formula.NO_RECORD
    assert result.checked is False
    assert result.confirms is False
    assert "unchecked" in result.note


def test_an_unparseable_formula_is_reported_as_unparseable():
    result = formula.compare("WAG994", "C2H8ClN")
    assert result.status == formula.UNPARSEABLE
    assert result.checked is False
    assert "registry formula could not be parsed" in result.note


# ------------------------------------------------- the indefinite multiplier


@pytest.mark.parametrize(
    "registry_formula,registry_name",
    [
        ("C19H23N.xC4H4O4", "some amine, maleate (1:?)"),
        ("C20H23N.xC4H4O4", "spiro compound, maleate"),
        ("C4H6N2O2S.xH2O4S", "propenoic acid, sulfate (9 CI)"),
        ("C12H18N2O", "something, maleate (1:?)"),
        ("C12H18N2O", "something, maleate (?:1)"),
    ],
)
def test_an_indefinite_multiplier_is_never_read_as_one(registry_formula, registry_name):
    """Handoff case 1, the comparator half. Three records reported "formula agrees".

    The registry writes ``.xC4H4O4`` and ``(1:?)`` when it declines to fix the
    component ratio. Reading the ``x`` as 1 turned three unconfirmed records into
    confirmations.
    """
    result = formula.compare(
        registry_formula, "C19H23N.C4H4O4", registry_name=registry_name
    )
    assert result.status == formula.RATIO_UNKNOWN
    assert result.confirms is False
    assert result.checked is False


def test_the_indefinite_gate_runs_before_parsing():
    """By the time a multiplier is parsed, the fact that it was indefinite is gone."""
    indefinite = "C19H23N.xC4H4O4"
    assert formula.parse(indefinite) is None
    result = formula.compare(indefinite, "C19H23N.C4H4O4")
    assert result.status == formula.RATIO_UNKNOWN, (
        "must be ratio-unknown, not unparseable: the registry said something "
        "specific, which is that it does not know"
    )


def test_the_registrys_own_ratio_marker_is_honoured():
    result = formula.compare(
        "C4H6N2O2S.H2O4S", "C4H8N2O6S2", registry_ratio="unknown-to-CAS"
    )
    assert result.status == formula.RATIO_UNKNOWN


def test_a_definite_ratio_is_compared_normally():
    result = formula.compare(
        "C19H23N.C4H4O4", "C23H27NO4", registry_name="some amine, maleate (1:1)"
    )
    assert result.status == formula.AGREE


@pytest.mark.parametrize(
    "text,indefinite",
    [
        ("C19H23N.xC4H4O4", True),
        ("C19H23N.nC4H4O4", True),
        ("maleate (1:?)", True),
        ("maleate (?:1)", True),
        ("C19H23N.C4H4O4", False),
        ("C19H23N.2C4H4O4", False),
        ("Xenon compound", False),
    ],
)
def test_indefinite_detection(text, indefinite):
    assert formula.is_indefinite(text) is indefinite


# ------------------------------------------------------------- the CAS table


def test_the_registry_loads_and_indexes_by_cas(tmp_path):
    path = registry(
        tmp_path,
        {"cas": "64-17-5", "molecular_formula": "C2H6O", "source_pdf": "a.pdf"},
        {"cas": "557-66-4", "molecular_formula": "C2H8ClN", "source_pdf": "b.pdf"},
    )
    table = casreg.load(path)
    assert len(table) == 2
    assert table.get("64-17-5").molecular_formula == "C2H6O"
    assert table.get(" 64-17-5 ").molecular_formula == "C2H6O"
    assert table.get("1-1-1") is None
    assert table.source_pdfs == ("a.pdf", "b.pdf")
    assert len(table.sha256) == 64


def test_a_row_whose_cas_spelling_needs_repair_is_still_found(tmp_path):
    """Every lookup passes the normalised `ctx.cas_key`; the table used the spelling.

    So a row written `0557-66-4` could never be found, and EXT-01 and EXT-04
    reported "no registry row" for a record whose row was in the table. The
    registry is assembled by hand from PDF exports, so a leading zero is the
    expected kind of defect.
    """
    loaded = casreg.load(
        registry(tmp_path, {"cas": "0557-66-4", "molecular_formula": "C2H7N.ClH"})
    )
    row = loaded.get("557-66-4")
    assert row is not None
    assert row.cas == "0557-66-4", "the row keeps the spelling the table used"
    assert loaded.get("0557-66-4") is not None
    assert loaded.coverage({"557-66-4"})["with_registry_row"] == 1


def test_coverage_says_how_much_of_a_batch_the_table_can_speak_to(tmp_path):
    """Clearance on a record with no row is "internally consistent", not "confirmed"."""
    table = casreg.load(
        registry(tmp_path, {"cas": "64-17-5", "molecular_formula": "C2H6O"})
    )
    assert table.coverage({"64-17-5", "557-66-4"}) == {
        "records": 2,
        "with_registry_row": 1,
        "without_registry_row": 1,
        "registry_rows": 1,
    }


def test_stoichiometry_never_comes_from_the_structure_string(tmp_path):
    """Handoff case 4, reproduced from the real table.

    Both ratio-unknown records in the batch have a formula saying the ratio is
    indefinite and an InChI carrying a definite 1:1. The formula is authoritative
    about composition; the structure string is not. Reading the ratio off the
    InChI would turn "the registry does not know" into "confirmed 1:1".
    """
    path = registry(
        tmp_path,
        {
            "cas": "207455-21-8",
            "molecular_formula": "C20H23N.xC4H4O4",
            "name_ratio": "unknown-to-CAS",
            "inchi": "InChI=1S/C20H23N.C4H4O4/c1-2-6-17...",
        },
    )
    row = casreg.load(path).get("207455-21-8")

    assert "x" in row.molecular_formula
    assert "C20H23N.C4H4O4" in row.structure_strings["inchi"]

    from_formula = formula.compare(
        row.molecular_formula,
        "C24H27NO4",
        registry_name=row.registry_name,
        registry_ratio=row.name_ratio,
    )
    assert from_formula.status == formula.RATIO_UNKNOWN

    from_structure = formula.compare("C20H23N.C4H4O4", "C24H27NO4")
    assert from_structure.status == formula.AGREE, (
        "this is the wrong answer the structure string would have given"
    )


def test_structure_columns_are_named_as_structure_columns():
    """Using one for stoichiometry should be a visible decision, not a reach into raw."""
    assert set(casreg.STRUCTURE_COLUMNS).isdisjoint(casreg.STOICHIOMETRY_COLUMNS)
    assert "molecular_formula" in casreg.STOICHIOMETRY_COLUMNS
    assert "inchi" in casreg.STRUCTURE_COLUMNS


def test_a_missing_table_says_where_it_comes_from(tmp_path):
    with pytest.raises(casreg.RegistryError, match="never fetches"):
        casreg.load(tmp_path / "nope.csv")


def test_a_missing_column_fails_the_load(tmp_path):
    path = tmp_path / "cas_registry.csv"
    path.write_text("cas,molecular_formula\n64-17-5,C2H6O\n")
    with pytest.raises(casreg.RegistryError, match="missing columns"):
        casreg.load(path)


def test_a_duplicate_cas_fails_the_load(tmp_path):
    path = registry(
        tmp_path,
        {"cas": "64-17-5", "molecular_formula": "C2H6O"},
        {"cas": "64-17-5", "molecular_formula": "C2H6O2"},
    )
    with pytest.raises(casreg.RegistryError, match="duplicate cas"):
        casreg.load(path)


def test_an_empty_cas_fails_the_load(tmp_path):
    with pytest.raises(casreg.RegistryError, match="empty cas"):
        casreg.load(registry(tmp_path, {"cas": "", "molecular_formula": "C2H6O"}))


def test_the_empty_registry_is_usable_as_a_no_op():
    """So a run with no CAS evidence takes the same code path, reporting unchecked."""
    assert len(casreg.EMPTY) == 0
    assert casreg.EMPTY.get("64-17-5") is None
    assert casreg.EMPTY.coverage({"64-17-5"})["without_registry_row"] == 1
