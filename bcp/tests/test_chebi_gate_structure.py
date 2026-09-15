"""Structural facts derived from a molfile, on which every other check depends.

The gate compares a record's *name* and its *asserted class* against what was
actually drawn. If this module is wrong, every one of those comparisons is wrong
in the same direction and the report still looks clean, so these tests pin each
derived field against a fixture whose chemistry is known by construction.
"""

from __future__ import annotations

import pytest

from chebi_gate import sdf, structure
from tests.chebi_gate_helpers import molblock, record_text, sdf_bytes


def analyse(mol: str, **fields: str) -> structure.Structure:
    """Parse a one-record SDF built from a molfile fixture and analyse it."""
    text = record_text(molblock(mol), {"NAME": mol, "CAS_NO": "64-17-5", **fields})
    record = sdf.parse_bytes(sdf_bytes(text)).records[0]
    return structure.analyse(record)


# ------------------------------------------------------------------- fragments


def test_single_fragment_neutral_compound():
    s = analyse("ethanol")
    assert s.parse is True
    assert s.formula == "C2H6O"
    assert s.parent_formula == "C2H6O"
    assert s.n_frag == 1
    assert s.counter == []
    assert s.frags == {"C2H6O": 1}
    assert s.net_charge == 0


def test_salt_splits_into_parent_and_counterion():
    s = analyse("amine_hcl")
    assert s.n_frag == 2
    assert s.counter == ["HCl"]
    assert s.parent_formula == "C2H7N"
    assert s.formula == "C2H8ClN"
    assert s.frags == {"C2H7N": 1, "HCl": 1}


def test_repeated_counterion_is_counted_twice():
    """A dihydrochloride has two HCl fragments, and both have to survive.

    The counterion list is built by object identity, not by formula, precisely so
    the second copy is not mistaken for the parent and dropped -- which would make
    a dihydrochloride's stoichiometry read as 1:1 and pass CON-02.
    """
    s = analyse("amine_2hcl")
    assert s.n_frag == 3
    assert s.counter == ["HCl", "HCl"]
    assert s.frags == {"C2H7N": 1, "HCl": 2}


def test_water_and_ethanol_counterions_are_identified():
    hydrate = analyse("amine_hydrate")
    assert hydrate.water == 1
    assert hydrate.ethanol == 0
    assert hydrate.counter == ["H2O"]

    ethanolate = analyse("amine_ethanolate")
    assert ethanolate.ethanol == 1
    assert ethanolate.water == 0
    assert ethanolate.counter == ["C2H6O"]


def test_anhydrous_record_reports_no_water():
    assert analyse("amine_hcl").water == 0


# ------------------------------------------------------- maleate versus fumarate


def test_maleate_double_bond_reads_as_z():
    """Maleate and fumarate share a formula and differ only in this geometry."""
    assert analyse("decylamine_maleate").geoms == ["Z"]


def test_fumarate_double_bond_reads_as_e():
    assert analyse("decylamine_fumarate").geoms == ["E"]


def test_geometry_is_only_reported_for_butenedioate_counterions():
    assert analyse("amine_hcl").geoms == []
    assert analyse("ethanol").geoms == []


# ----------------------------------------------------------- cationic nitrogen


def test_quaternary_ammonium_is_counted_as_a_cation_with_no_hydrogen():
    """A quaternary cation cannot be a hydrohalide: it has no H to have received."""
    assert analyse("quat_ammonium_br").cation_n_noh == 1


def test_protonated_ammonium_with_explicit_hydrogens_is_not_quaternary():
    """Regression: explicit H atoms must still count as hydrogens on the nitrogen.

    ``GetTotalNumHs()`` counts implicit and property hydrogens only, so a molfile
    that draws its hydrogens as atoms reports zero and a legitimate ionic
    hydrochloride is misread as a quaternary salt, failing CON-01. The fixture
    draws them explicitly on purpose.
    """
    s = analyse("ionic_amine_chloride")
    assert s.counter == ["Cl-"]
    assert s.net_charge == 0
    assert s.cation_n_noh == 0


# -------------------------------------------------------------- structure keys


def test_inchikey_and_skeleton_agree():
    s = analyse("ethanol")
    assert s.inchikey is not None
    assert s.skeleton == s.inchikey[:14]
    assert len(s.skeleton) == 14


def test_two_salt_forms_of_one_base_share_a_parent_skeleton():
    """This is how the gate finds a record's parent and spots duplicate salt forms."""
    mono = analyse("amine_hcl")
    di = analyse("amine_2hcl")
    assert mono.inchikey != di.inchikey
    assert mono.base_inchikey == di.base_inchikey
    assert mono.parent_skeleton == di.parent_skeleton


def test_neutralising_changes_the_key_of_an_ionic_drawing():
    s = analyse("ionic_amine_chloride")
    assert s.base_inchikey is not None
    assert s.base_inchikey != s.inchikey


# ------------------------------------------------------------- counts line bits


def test_chiral_flag_zero_by_default():
    assert analyse("ethanol").chiral_flag == 0


def test_chiral_flag_one_is_read():
    assert analyse("chiral_flag_absolute").chiral_flag == 1


def test_blank_chiral_flag_field_reads_as_zero_without_raising():
    """Regression: ``int(" ")`` raises, and ``" " or 0`` does not save you.

    The prototype wrote ``int(counts[12:15] or 0)``. A blank-padded chiral-flag
    field is a non-empty string of spaces, which is truthy, so the ``or 0`` never
    fires and the whole run dies on one malformed counts line.
    """
    assert analyse("chiral_flag_blank").chiral_flag == 0


def test_unparseable_chiral_flag_warns_and_reads_as_zero(caplog):
    # Mutate the chiral-flag field in place: columns 13-15 of the counts line,
    # which is a fixed-width format, so a substring replace would land elsewhere.
    block = molblock("ethanol")
    parts = block.split("\n")
    parts[3] = parts[3][:12] + " xx" + parts[3][15:]
    block = "\n".join(parts)
    text = record_text(block, {"NAME": "broken flag"})
    record = sdf.parse_bytes(sdf_bytes(text)).records[0]

    assert structure._chiral_flag(record.counts_line) == 0
    assert "unparseable chiral flag" in caplog.text


def test_zero_coordinates_are_detected():
    assert analyse("zero_coords").zero_coords is True
    assert analyse("ethanol").zero_coords is False


# ------------------------------------------------------------- failure handling


def test_unparseable_molfile_reports_parse_false_and_no_facts():
    """Every field must stay at its default so no check reads a fact into silence."""
    s = analyse("unparseable")
    assert s.parse is False
    assert s.formula is None
    assert s.inchikey is None
    assert s.base_inchikey is None
    assert s.n_frag == 0
    assert s.frags == {}
    assert s.counter == []
    assert s.net_charge is None
    assert s.skeleton is None
    assert s.parent_skeleton is None


def test_unparseable_molfile_still_reports_the_chiral_flag():
    """The counts line is text, so it is readable even when the atom block is not."""
    assert analyse("unparseable").chiral_flag == 0


@pytest.mark.parametrize(
    "mol",
    [
        "ethanol",
        "amine_hcl",
        "amine_2hcl",
        "decylamine_maleate",
        "decylamine_fumarate",
        "quat_ammonium_br",
        "ionic_amine_chloride",
        "amine_hydrate",
        "amine_ethanolate",
    ],
)
def test_every_parseable_fixture_yields_a_formula_and_a_key(mol):
    """EXT-04 compares the registry formula against this one, so it cannot be None.

    The prototype never set it: ``analyse_structure`` returned no ``formula`` key,
    the gate passed ``s.get("formula")`` to the comparator, and the comparator's
    ``if not a or not b: return ""`` silently dropped every comparison. The gate
    printed "checks run: ... EXT-04" and produced zero EXT-04 findings on all 195
    salt records with a full evidence cache.
    """
    s = analyse(mol)
    assert s.formula
    assert s.parent_formula
    assert s.inchikey
