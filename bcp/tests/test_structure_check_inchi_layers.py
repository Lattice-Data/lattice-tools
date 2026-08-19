"""Unit tests for structure_check.inchi layer comparison (no network)."""

from __future__ import annotations

import pytest

from structure_check.inchi import (
    COMPARED_LAYERS,
    CONSTITUTION_LAYERS,
    DIFFERENT_COMPOUND,
    IONISATION_LAYERS,
    ISOTOPE_DIFFERS,
    ISOTOPE_LAYER,
    FORM_DIFFERS,
    INCHI_VERDICTS,
    LAYERS_IDENTICAL,
    LAYERS_NOT_COMPARABLE,
    STEREO_DIFFERS,
    STEREO_LAYERS,
    STEREO_UNDEFINED_ON_ONE_SIDE,
    classify_pair,
    defined_side,
    defined_stereo,
    inchi_layers,
    is_multi_component,
    principal_component,
)

# Real structures from the curation run that exposed the bug this module fixes.
# Every pair below was verified as one compound by the previous implementation.

# CCMI / AVL-3288. /b present on both sides with opposite values: a genuine E/Z
# pair, which the old comparison could not see because it never read /b.
CCMI_Z = (
    "InChI=1S/C19H15Cl2N3O2/c1-12-10-18(26-24-12)17(11-22-15-6-2-13(20)3-7-15)"
    "19(25)23-16-8-4-14(21)5-9-16/h2-11,22H,1H3,(H,23,25)/b17-11-"
)
CCMI_E = CCMI_Z.replace("/b17-11-", "/b17-11+")

# Dihydrosphingosine. Identical /t, opposite /m -- the definition of an enantiomer
# pair. The old comparison reached its final `else` and called this a PubChem
# cataloguing duplicate.
DL_ERYTHRO_DIHYDROSPHINGOSINE = (
    "InChI=1S/C18H39NO2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-18(21)17(19)16-20/"
    "h17-18,20-21H,2-16,19H2,1H3/t17-,18+/m1/s1"
)
SPHINGANINE = DL_ERYTHRO_DIHYDROSPHINGOSINE.replace("/m1/", "/m0/")

# Formoterol fumarate against (S,S)-formoterol fumarate: the same /m flip inside a
# multi-component salt, so it exercises the formula split at the same time. The
# identifier this was given was arformoterol's -- a single enantiomer's entry
# assigned to a racemate.
FORMOTEROL_FUMARATE = (
    "InChI=1S/2C19H24N2O4.C4H4O4/c2*1-13(9-14-3-6-16(25-2)7-4-14)20-11-19(24)"
    "15-5-8-18(23)17(10-15)21-12-22;5-3(6)1-2-4(7)8/h2*3-8,10,12-13,19-20,23-24H,"
    "9,11H2,1-2H3,(H,21,22);1-2H,(H,5,6)(H,7,8)/b;;2-1+/t2*13-,19+;/m11./s1"
)
SS_FORMOTEROL_FUMARATE = FORMOTEROL_FUMARATE.replace("/m11./", "/m00./")

# SU 4312. The name side specifies /b, the CAS side specifies no stereo at all:
# PubChem holding two records for one substance, which is agreement.
SU4312_UNDEFINED = (
    "InChI=1S/C17H16N2O/c1-19(2)13-9-7-12(8-10-13)11-15-14-5-3-4-6-16(14)18-"
    "17(15)20/h3-11H,1-2H3,(H,18,20)"
)
SU4312_Z = SU4312_UNDEFINED + "/b15-11-"

# Doxorubicin free base against its hydrochloride: identical principal component,
# different formula layer.
DOXORUBICIN = (
    "InChI=1S/C27H29NO11/c1-10-22(31)13(28)6-17(38-10)39-15-8-27(36,16(30)9-29)"
    "7-12-19(15)26(35)21-20(24(12)33)23(32)11-4-3-5-14(37-2)18(11)25(21)34/"
    "h3-5,10,13,15,17,22,29,31,33,35-36H,6-9,28H2,1-2H3/t10-,13-,15-,17-,22+,27-/"
    "m0/s1"
)
DOXORUBICIN_HCL = (
    "InChI=1S/C27H29NO11.ClH/c1-10-22(31)13(28)6-17(38-10)39-15-8-27(36,16(30)"
    "9-29)7-12-19(15)26(35)21-20(24(12)33)23(32)11-4-3-5-14(37-2)18(11)25(21)34;/"
    "h3-5,10,13,15,17,22,29,31,33,35-36H,6-9,28H2,1-2H3;1H/"
    "t10-,13-,15-,17-,22+,27-;/m0./s1"
)

ETHANOL = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

# ABT-702 as its mono- and di-hydrochloride: the same base, differing only in how
# many equivalents of HCl. Real PubChem records, and exactly the salt distinction
# the strict policy turns on.
ABT702_HCL = (
    "InChI=1S/C22H19BrN6O.ClH/c23-16-3-1-2-14(10-16)17-11-18(28-22-20(17)21(24)"
    "26-13-27-22)15-4-5-19(25-12-15)29-6-8-30-9-7-29;/h1-5,10-13H,6-9H2,"
    "(H2,24,26,27,28);1H"
)
ABT702_2HCL = (
    "InChI=1S/C22H19BrN6O.2ClH/c23-16-3-1-2-14(10-16)17-11-18(28-22-20(17)21(24)"
    "26-13-27-22)15-4-5-19(25-12-15)29-6-8-30-9-7-29;;/h1-5,10-13H,6-9H2,"
    "(H2,24,26,27,28);2*1H"
)

# Apomorphine hydrochloride against its hemihydrate. Real PubChem records, and the
# only fixture here whose *principal* component carries a stoichiometric multiplier
# on one side and not the other -- 2C17H17NO2 against C17H17NO2. The ABT-702 pair
# does not pin that, because there the multiplier sits on the counterion.
APOMORPHINE_HCL = (
    "InChI=1S/2C17H17NO2.2ClH.H2O/c2*1-18-8-7-10-3-2-4-12-15(10)13(18)9-11-5-6-"
    "14(19)17(20)16(11)12;;;/h2*2-6,13,19-20H,7-9H2,1H3;2*1H;1H2/t2*13-;;;/m11.../s1"
)
APOMORPHINE_HCL_HEMIHYDRATE = (
    "InChI=1S/C17H17NO2.ClH.H2O/c1-18-8-7-10-3-2-4-12-15(10)13(18)9-11-5-6-14(19)"
    "17(20)16(11)12;;/h2-6,13,19-20H,7-9H2,1H3;1H;1H2/t13-;;/m1../s1"
)

# Pyrvinium pamoate: the drug cation and the pamoate counterion have formulas of
# equal character length, so ranking components by string length picks the
# counterion. This is the same wrong answer PubChem's desalting gives.
PYRVINIUM_PAMOATE_FORMULA = "2C26H28N3.C23H16O6"


# --------------------------------------------------------------------------
# Layer parsing
# --------------------------------------------------------------------------


def test_inchi_layers_splits_formula_from_layers() -> None:
    layers, formula = inchi_layers(ETHANOL)
    assert formula == "C2H6O"
    assert layers["c"] == "1-2-3"
    assert layers["h"] == "3H,2H2,1H3"


@pytest.mark.parametrize("raw", ["", "not an inchi", "InChI=1S", "InChI=1S/"])
def test_inchi_layers_returns_empty_rather_than_raising(raw: str) -> None:
    """Callers are on paths that promise not to raise; unparseable must be inert.

    Both halves are asserted separately: an `or` here would pass while one of them
    silently returned junk.
    """
    layers, formula = inchi_layers(raw)
    assert layers == {}
    assert formula == ""


def test_defined_stereo_counts_all_four_layers() -> None:
    assert defined_stereo(inchi_layers(SU4312_UNDEFINED)[0]) == 0
    assert defined_stereo(inchi_layers(SU4312_Z)[0]) == 1
    assert defined_stereo(inchi_layers(SPHINGANINE)[0]) == 3
    # Spelled 4, not len(STEREO_LAYERS): comparing against the constant under test
    # would keep passing if a layer were dropped from it.
    assert defined_stereo(inchi_layers(FORMOTEROL_FUMARATE)[0]) == 4
    assert len(STEREO_LAYERS) == 4


# --------------------------------------------------------------------------
# Principal component
# --------------------------------------------------------------------------


def test_principal_component_strips_counterion() -> None:
    assert principal_component("C27H29NO11.ClH") == "C27H29NO11"


def test_principal_component_picks_the_drug_not_the_pamoate() -> None:
    """Ranking by atom count, not formula-string length.

    Pyrvinium is C26H28N3 and pamoate is C23H16O6: both eight characters once the
    stoichiometric multiplier is dropped, so a string tiebreak returns the
    counterion. Atom counts are 57 against 45.
    """
    assert principal_component(PYRVINIUM_PAMOATE_FORMULA) == "2C26H28N3"


def test_principal_component_of_single_component_is_itself() -> None:
    assert principal_component("C2H6O") == "C2H6O"
    assert principal_component("") == ""


def test_is_multi_component_detects_salts() -> None:
    assert is_multi_component(DOXORUBICIN_HCL) is True
    assert is_multi_component(DOXORUBICIN) is False


@pytest.mark.parametrize(
    "raw", ["", "   ", "C2H6O", "not an inchi", "a/b", "InChI=1S/??/c1", "InChI=1S"]
)
def test_an_unreadable_structure_has_no_component_count(raw: str) -> None:
    """`None`, not `False`: exactly the inputs `classify_pair()` refuses.

    Each of these has no parseable formula layer, and `False` would report that as
    "one component" -- indistinguishable from a genuine free base. The answer
    decides whether a compound absent from ChEBI needs a new entry or a salt on an
    existing parent, so the wrong direction here is a wrong submission.
    """
    assert is_multi_component(raw) is None
    assert classify_pair(raw, DOXORUBICIN) == LAYERS_NOT_COMPARABLE


# --------------------------------------------------------------------------
# The comparison, on the structures that shipped wrong
# --------------------------------------------------------------------------


def test_opposite_double_bond_geometry_is_a_stereoisomer_difference() -> None:
    """CCMI: /b on both sides with different values, so two different molecules.

    Verified as one compound before /b was read at all.
    """
    assert classify_pair(CCMI_Z, CCMI_E) == STEREO_DIFFERS


def test_identical_t_with_opposite_m_is_a_stereoisomer_difference() -> None:
    """Dihydrosphingosine: /t equal, /m opposite -- an enantiomer pair.

    This fell through to the wrong branch and was reported as a cataloguing
    duplicate, which promoted one enantiomer's ChEBI entry onto the other.
    """
    assert classify_pair(DL_ERYTHRO_DIHYDROSPHINGOSINE, SPHINGANINE) == STEREO_DIFFERS


def test_m_flip_inside_a_multi_component_salt_is_still_a_stereo_difference() -> None:
    """Formoterol fumarate: the salt must not mask the /m11 against /m00 flip."""
    assert classify_pair(FORMOTEROL_FUMARATE, SS_FORMOTEROL_FUMARATE) == STEREO_DIFFERS


def test_layer_absent_on_one_side_is_a_pubchem_duplicate_not_a_difference() -> None:
    """SU 4312: one record says (Z), the other says nothing. Same substance."""
    verdict = classify_pair(SU4312_UNDEFINED, SU4312_Z)
    assert verdict == STEREO_UNDEFINED_ON_ONE_SIDE


def test_defined_side_prefers_the_more_specific_structure() -> None:
    """The anchor must be the side that specifies the stereochemistry.

    Anchoring on the vaguer side matches a stereo-unspecified database entry to a
    compound whose geometry is specified, which is the substitution the strict
    policy exists to reject.
    """
    assert defined_side(SU4312_UNDEFINED, SU4312_Z) == SU4312_Z
    assert defined_side(SU4312_Z, SU4312_UNDEFINED) == SU4312_Z


def test_salt_and_free_base_differ_by_form_not_by_compound() -> None:
    assert classify_pair(DOXORUBICIN, DOXORUBICIN_HCL) == FORM_DIFFERS


def test_different_stoichiometry_of_one_base_is_a_form_difference() -> None:
    """ABT-702 mono-HCl against di-HCl: two forms of one compound.

    Here the multiplier sits on the *counterion*, so this pins the formula-layer
    comparison but NOT the multiplier stripping -- see the apomorphine test below.
    """
    assert classify_pair(ABT702_HCL, ABT702_2HCL) == FORM_DIFFERS


def test_a_multiplier_on_the_principal_component_is_ignored() -> None:
    """Apomorphine hydrochloride against its hemihydrate: 2C17H17NO2 vs C17H17NO2.

    The stoichiometric multiplier has to be stripped before the principal
    components are compared, or one apomorphine and two compare as *different
    molecules* and a hydrate question is reported as a wrong-compound question.

    This is the only fixture that pins that stripping. Deleting both
    `_LEADING_COUNT.sub()` calls from `classify_pair` left every other test in this
    module passing.
    """
    assert classify_pair(APOMORPHINE_HCL, APOMORPHINE_HCL_HEMIHYDRATE) == FORM_DIFFERS


def test_unrelated_compounds_differ_by_compound() -> None:
    assert classify_pair(ETHANOL, DOXORUBICIN) == DIFFERENT_COMPOUND


def test_a_structure_is_identical_to_itself() -> None:
    for inchi in (ETHANOL, CCMI_Z, FORMOTEROL_FUMARATE, DOXORUBICIN_HCL):
        assert classify_pair(inchi, inchi) == LAYERS_IDENTICAL


@pytest.mark.parametrize(
    ("a", "b"),
    [
        ("", ""),
        (ETHANOL, ""),
        ("", ETHANOL),
        (ETHANOL, "not an inchi"),
        (ETHANOL, "a/b"),
        ("InChI=1S//c1-2-3", "InChI=1S//c1-2-3"),
    ],
)
def test_missing_structure_is_not_comparable_never_identical(a: str, b: str) -> None:
    """Silence from a database must never render as agreement between structures.

    The predecessor returned "identical" for two empty strings and escaped the
    consequences only because its one caller guarded the call site.
    """
    assert classify_pair(a, b) == LAYERS_NOT_COMPARABLE


def test_comparison_is_symmetric_under_argument_order() -> None:
    pairs = [
        (CCMI_Z, CCMI_E),
        (DL_ERYTHRO_DIHYDROSPHINGOSINE, SPHINGANINE),
        (SU4312_UNDEFINED, SU4312_Z),
        (DOXORUBICIN, DOXORUBICIN_HCL),
        (ETHANOL, DOXORUBICIN),
        (ETHANOL, ""),
    ]
    for a, b in pairs:
        assert classify_pair(a, b) == classify_pair(b, a)


def test_every_declared_verdict_is_reachable() -> None:
    """Exhaustiveness, not membership.

    Asserting each result is *in* INCHI_VERDICTS is nearly a tautology -- the
    function only ever returns those constants. What is worth pinning is the other
    direction: every verdict the module declares is produced by some real pair, so a
    constant cannot be declared and then never returned, and a consumer's mapping
    cannot silently go stale.
    """
    observed = {
        classify_pair(a, b)
        for a, b in [
            (ETHANOL, ETHANOL),
            (DOXORUBICIN, DOXORUBICIN_HCL),
            (CCMI_Z, CCMI_E),
            (SU4312_UNDEFINED, SU4312_Z),
            (ETHANOL, ETHANOL_D6),
            (ETHANOL, DIMETHYL_ETHER),
            (ETHANOL, ""),
        ]
    }
    assert observed == set(INCHI_VERDICTS)


# --------------------------------------------------------------------------
# Same formula, different molecule: what the stereo-only comparison could not see
# --------------------------------------------------------------------------

# Constitutional isomers. Each pair shares a formula layer and differs only in /c
# or /h, which the first version of this module never read -- so all of these
# compared LAYERS_IDENTICAL, the verdict meaning "the same compound".
DIMETHYL_ETHER = "InChI=1S/C2H6O/c1-3-2/h1-2H3"
N_BUTANE = "InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3"
ISOBUTANE = "InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3"
CATECHOL = "InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H"
HYDROQUINONE = "InChI=1S/C6H6O2/c7-5-1-2-6(8)4-3-5/h1-4,7-8H"

ACETIC_ACID = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
ACETATE = ACETIC_ACID + "/p-1"
ETHANOL_D6 = ETHANOL + "/i1D3,2D2,3D"

# Diatrizoate meglumine: a hydrogen-rich counterion against an iodinated parent.
# Meglumine wins on total atoms 30 to 29; diatrizoate wins on heavy atoms 20 to 13.
DIATRIZOATE_MEGLUMINE_FORMULA = "C11H9I3N2O4.C7H17NO5"


@pytest.mark.parametrize(
    ("a", "b", "label"),
    [
        (ETHANOL, DIMETHYL_ETHER, "ethanol / dimethyl ether"),
        (N_BUTANE, ISOBUTANE, "n-butane / isobutane"),
        (CATECHOL, HYDROQUINONE, "catechol / hydroquinone"),
    ],
)
def test_constitutional_isomers_are_different_compounds(
    a: str, b: str, label: str
) -> None:
    """One formula, two molecules -- separated only by /c and /h.

    Comparing the formula and the stereo layers while never reading the
    connectivity meant every one of these pairs reported as the same compound. This
    is the module's worst possible failure direction, and it mattered most in its
    intended slot: the InChIKey path only asks for a finer comparison once it has
    already established that connectivity differs.
    """
    assert classify_pair(a, b) == DIFFERENT_COMPOUND, label


def test_a_truncated_inchi_does_not_compare_equal_to_the_intact_string() -> None:
    """A spreadsheet cell limit clips the tail of an InChI.

    The clipped string keeps its formula and loses layers, so a comparison that
    reads only the formula and the stereo layers called it the same compound.
    Refusing is the safe direction.
    """
    intact = "InChI=1S/C6H12O6/c1-2-3(7)4(8)5(9)6(10)11/h2-6,8-11H,7H2"
    clipped = "InChI=1S/C6H12O6/c1-2-3(7)4(8)5(9)6(10)11/h2-6"
    assert classify_pair(intact, clipped) == DIFFERENT_COMPOUND


def test_an_isotopologue_is_not_the_unlabelled_compound() -> None:
    """/i leaves the formula layer alone, so it slipped past a formula gate.

    Ethanol-d6 has its own InChIKey, its own CAS number and its own ChEBI entry.
    It is neither a salt form nor a stereoisomer of ethanol.
    """
    assert classify_pair(ETHANOL, ETHANOL_D6) == ISOTOPE_DIFFERS


def test_a_conjugate_base_is_a_different_form_not_the_same_structure() -> None:
    """/p also leaves the formula layer as the neutral parent's.

    ChEBI holds acetic acid and acetate as separate entries, so matching one to
    the other is the same class of error as matching a salt to its free base.
    """
    assert classify_pair(ACETIC_ACID, ACETATE) == FORM_DIFFERS
    assert classify_pair(ACETIC_ACID, ACETIC_ACID + "/q+1") == FORM_DIFFERS


def test_isotopic_sublayers_do_not_overwrite_the_real_stereo_layers() -> None:
    """InChI re-emits /t /m /s after /i to give the isotopic stereochemistry.

    Keying layers by first character and assigning meant the isotopic /m replaced
    the real one, so a deuterium label could manufacture a stereo difference.
    """
    labelled = "InChI=1S/C3H8O/c1-2-3/h3H,2H2,1H3/t3-/m0/s1/i1D3/t3-/m1/s1"
    layers, _ = inchi_layers(labelled)
    assert layers["m"] == "0"
    assert layers["i"] == "1D3"


def test_principal_component_ranks_heavy_atoms_before_hydrogens() -> None:
    """Diatrizoate meglumine: the counterion has more atoms but fewer heavy ones.

    Ranking on total atoms alone let hydrogen outvote iodine, so the meglumine
    counterion was returned as the parent and the salt compared as a different
    compound from its own free acid.
    """
    assert principal_component(DIATRIZOATE_MEGLUMINE_FORMULA) == "C11H9I3N2O4"
    # and the pyrvinium tie still breaks the right way on total atoms
    assert principal_component(PYRVINIUM_PAMOATE_FORMULA) == "2C26H28N3"


@pytest.mark.parametrize(
    ("a", "b"), [("a/b", "x/b"), ("InChI=1S/??/c1", "InChI=1S/??/c1")]
)
def test_a_string_that_merely_contains_a_slash_is_not_comparable(
    a: str, b: str
) -> None:
    """The gate matches a formula pattern rather than checking for non-emptiness."""
    assert classify_pair(a, b) == LAYERS_NOT_COMPARABLE


# --------------------------------------------------------------------------
# Layers this comparison does not know how to read
# --------------------------------------------------------------------------

# Non-standard InChI. Standard InChI is `InChI=1S`; the non-standard flavour is
# `InChI=1` and can carry two layers that do not exist in the standard one -- `/f`,
# the fixed-H layer, which names a specific tautomer, and `/r`, the reconnected
# layer, which restores the metal bonds standard InChI breaks. Both change which
# molecule the string denotes, and neither is in COMPARED_LAYERS.
NONSTANDARD_ACETIC_ACID = "InChI=1/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
NONSTANDARD_ACETIC_ACID_FIXED_H = NONSTANDARD_ACETIC_ACID + "/f/h3H"


def test_a_non_standard_inchi_gets_no_verdict() -> None:
    """The tiers read a fixed layer list, complete for standard InChI only.

    Nothing in the parser rejects `InChI=1`, so before this gate a non-standard
    string was compared on the layers that happened to be recognised and the rest
    treated as absent on both sides -- a confident answer from an admittedly partial
    reading, which is the defect this module exists to remove.
    """
    assert classify_pair(NONSTANDARD_ACETIC_ACID, ACETIC_ACID) == LAYERS_NOT_COMPARABLE
    assert (
        classify_pair(NONSTANDARD_ACETIC_ACID, NONSTANDARD_ACETIC_ACID)
        == LAYERS_NOT_COMPARABLE
    )
    assert (
        classify_pair(NONSTANDARD_ACETIC_ACID_FIXED_H, NONSTANDARD_ACETIC_ACID)
        == LAYERS_NOT_COMPARABLE
    )


@pytest.mark.parametrize("layer", ["/f/h3H", "/rC2H6O", "/x1"])
def test_an_unrecognised_layer_gets_no_verdict(layer: str) -> None:
    """Refused on the layer itself, not only on the version string.

    Version and residue are checked separately so that a layer this module has
    never heard of is refused even when it arrives on a string claiming to be
    standard -- which is what makes COMPARED_LAYERS enforced rather than merely
    declared.
    """
    assert classify_pair(ETHANOL + layer, ETHANOL) == LAYERS_NOT_COMPARABLE
    assert classify_pair(ETHANOL, ETHANOL + layer) == LAYERS_NOT_COMPARABLE


def test_the_declared_layer_set_covers_every_layer_the_tiers_read() -> None:
    """COMPARED_LAYERS is the union of the tiers, not a hand-maintained copy.

    A tier that grew a layer without it being declared would refuse every real
    structure carrying that layer; a declared layer no tier reads would be silently
    ignored. Both are caught by building the set from the tiers themselves, which is
    what this pins.
    """
    assert COMPARED_LAYERS == set(
        CONSTITUTION_LAYERS + IONISATION_LAYERS + (ISOTOPE_LAYER,) + STEREO_LAYERS
    )
