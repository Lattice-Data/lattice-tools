"""ChEBI relationship classes, and whether a drawn structure supports the one asserted.

A salt record submitted to ChEBI carries a ``RELATIONSHIP`` field naming the class
it belongs to -- ``ISA36807`` for a hydrochloride, ``ISA50221`` for a maleate. The
code is what ChEBI's curators load; the check here is that the structure actually
drawn contains the counterion the code claims.

The table is deliberately a closed set, and an unrecognised code is a *finding*
rather than a pass. That is the behaviour handoff case 17 is about: adding
``ISA48369`` to the submission broke two records until the gate learned it, and the
right failure mode is "this gate does not know that code, someone confirm it"
rather than silent acceptance of a code nobody has checked.
"""

from __future__ import annotations

import re

from .formula import parse as parse_formula
from . import structure as structure_mod
from .structure import Structure

# Classes the two-file validator was written against.
CLASSES = {
    "ISA36807": "hydrochloride",
    "ISA48367": "hydrobromide",
    "ISA50221": "maleate",
    "ISA50921": "fumarate",
    "ISA64148": "oxalate",
    "ISA38700": "sodium",
    "ISA50356": "iodide",
    "ISA51337": "sulfate",
    "ISA64382": "sulfonate",
    "ISA38037": "methanesulfonate",
    "ISA50394": "potassium",
    "ISA50562": "tartrate",
}

# Classes added after the first submission round, kept separate so it stays visible
# which codes the original validator was never tested against.
#
# ISA48369 is the parent class for quaternary and pyridinium bromides. Those cannot
# be hydrobromides: their cation has no hydrogen to have received, which is the
# whole reason for reclassifying them. So the hydrohalide condition must NOT apply
# to this class -- a point worth stating, because applying the stricter test would
# look like the safer choice and would reject every correctly reclassified record.
EXTRA_CLASSES = {
    "ISA48369": "organic-bromide",
    "ISA35273": "quaternary-ammonium",
}

ALL_CLASSES = {**CLASSES, **EXTRA_CLASSES}

# Counterion formulae per class, as RDKit spells them ("Cl-", not "[Cl-]").
_HALIDES = {
    "hydrochloride": ("Cl", "HCl", "Cl-"),
    "hydrobromide": ("Br", "HBr", "Br-"),
}

# Solvent of crystallisation. Ignored when testing which counterions a class
# requires, because a hydrated hydrochloride is still a hydrochloride: the
# original rule demanded that *every* counterion be a halide, so any hydrate or
# ethanolate of a hydrohalide could never satisfy CON-01 and was reported as a
# class the structure does not support. Efonidipine hydrochloride monoethanolate
# is exactly that shape.
SOLVATES = (structure_mod.WATER_FORMULA, structure_mod.ETHANOL_FORMULA)
_IODIDE = ("I", "HI", "I-")
# The same ion list as the hydrobromide class; the classes differ in what they
# require of the *cation*, not in which counterion counts.
_ORGANIC_BROMIDE = _HALIDES["hydrobromide"]
_SODIUM = ("Na", "Na+")
_POTASSIUM = ("K", "K+")
_SULFATE = ("H2O4S", "HO4S-", "O4S-2")
_METHANESULFONATE = ("CH4O3S", "CH3O3S-")
# Protonation states of one counterion, as whole formulae with an optional charge
# suffix. They were prefix tests, and a prefix is not a formula: "C4H6O6" is a
# prefix of "C4H6O6S", so a fragment that merely begins like a tartrate counted as
# one. The charge suffix is matched explicitly instead, which is the only thing
# the prefix test was really there to allow for.
_OXALATE = ("C2H2O4", "C2HO4", "C2O4")
_TARTRATE = ("C4H6O6", "C4H5O6", "C4H4O6")

# ISA64382 is the parent class of the sulfonate counterions -- tosylate, mesylate,
# besylate, napsylate, napadisylate -- and cannot be a closed formula list the way
# the others are. The rule it replaces tested the substrings "S" and "O3" against
# every fragment including the parent, which is two separate mistakes. A drug
# parent can carry a sulfur and three oxygens of its own with no counterion in
# sight: six records of the reference batch have exactly that shape, and FR 122047
# hydrochloride -- whose only counterion is HCl -- reports the sulfonate class as
# supported if its RELATIONSHIP is mistyped as ISA64382 instead of ISA36807. And
# "S" is a substring of Si, Se, Sn and Sb, so the test was never about sulfur.
#
# Testing the counterions alone is not the fix either, and 3-Methyl-GABA
# napadisylate is why: the naphthalenedisulfonate is half again the size of the
# base, so it *is* the parent fragment by heavy-atom count and `counter` holds two
# copies of the base. What identifies a sulfonate counterion is its composition --
# carbon, hydrogen, oxygen and sulfur only, with three oxygens per sulfur -- and
# that reads the same whichever fragment happens to be largest.
_SULFONATE_ELEMENTS = frozenset({"C", "H", "O", "S"})


def _is_sulfonate(formula: str) -> bool:
    """Whether one fragment formula is a sulfonic acid or sulfonate anion."""
    counts = parse_formula(formula)
    if not counts or set(counts) - _SULFONATE_ELEMENTS:
        return False
    sulfur = counts.get("S", 0)
    return bool(sulfur) and counts.get("O", 0) == 3 * sulfur


# Charge suffix as RDKit spells a fragment: "Cl-", "Na+", "O4S-2".
_CHARGE_SUFFIX = re.compile(r"[+-]\d*$")


def _uncharged(formula: str) -> str:
    """A fragment formula without its charge suffix: ``"O4S-2"`` -> ``"O4S"``."""
    return _CHARGE_SUFFIX.sub("", formula)


# Butenedioate has no class table of its own: maleate and fumarate share a formula
# and are told apart by geometry, not composition. The stoichiometry counter still
# has to know it is a counterion, so it is listed here -- acid and both anions.
_BUTENEDIOATE = ("C4H4O4", "C4H3O4", "C4H2O4")

# Every formula the gate counts as a counterion, charge-stripped, assembled from
# the per-class tables above instead of restated beside them.
#
# CON-02 used to keep its own list and it held almost no anionic spellings, so an
# ionically drawn salt had its counterion counted as *base*: "X mono tosylate"
# drawn as a cation plus C7H7O3S- gave salt_n 0 and base_n 2, and CON-02 reported
# "implies 1 counterion(s) per base, structure has 0" at high on a record whose
# stoichiometry was exactly right. Hydrogen maleate, oxalate, tartrate and mesylate
# had the same shape, and the halide and alkali-metal rules were unaffected only
# because "Cl-" and "Na+" happened to be in both lists.
#
# Reachable but *latent* on the 290-record batch, and that is worth recording so
# nobody re-derives it as urgent: every ionic drawing in there is a halide or an
# alkali metal (Na+ 7, I- 5, Br- 2, K+ 2, all in both lists already) or a charged
# parent, which is correctly not a counterion. Reclassifying every fragment of all
# 290 records moves zero of them, so the baseline anchor is unchanged by this fix
# -- it is a latent false `high` on a drawing convention this batch happens not to
# use, not a finding that was firing.
#
# Having a second table *was* the defect, so this is deliberately the only one: a
# counterion added for a class rule is one the counter recognises for free.
COUNTERIONS = frozenset(
    _uncharged(x)
    for x in (
        *(ion for ions in _HALIDES.values() for ion in ions),
        *_IODIDE,
        *_ORGANIC_BROMIDE,
        *_SODIUM,
        *_POTASSIUM,
        *_SULFATE,
        *_METHANESULFONATE,
        *_OXALATE,
        *_TARTRATE,
        *_BUTENEDIOATE,
    )
)


def is_counterion(formula: str) -> bool:
    """Whether a fragment formula is a counterion rather than base or solvent.

    Indifferent to how the charge is spelled, because the same salt is drawn
    neutral by some depositors and ionic by others and the stoichiometry is
    identical either way. The sulfonates are recognised by composition rather than
    by list, exactly as :func:`class_supported` does it, so the whole family --
    tosylate, mesylate, besylate, napsylate, napadisylate -- is covered in both
    spellings without enumerating them.

    Solvates are not counterions and are not counted here; a caller that needs to
    exclude them from a base count checks :data:`SOLVATES` as well, because "is
    this solvent" and "is this a counterion" are two questions and CON-03 owns the
    first one.
    """
    return _uncharged(formula) in COUNTERIONS or _is_sulfonate(formula)


def class_name(code: str) -> str | None:
    """The class a ``RELATIONSHIP`` code names, or None if this gate does not know it."""
    return ALL_CLASSES.get(code.strip())


def class_supported(name: str, structure: Structure) -> bool:
    """Whether the drawn structure contains the counterion the class claims.

    Returns False for a class name this function has no rule for, so a code added
    to :data:`ALL_CLASSES` without a matching rule fails loudly instead of passing
    by default.
    """
    counter = structure.counter
    frags = structure.frags
    # Every rule below that demands *all* counterions match tests this rather than
    # `counter`, so that solvent of crystallisation cannot make the class fail. It
    # is computed once, for all of them: the exemption first existed only inside
    # the hydrohalide branch, and the two halide rules that kept testing `counter`
    # held a hydrated iodide and a hydrated quaternary bromide on a reason no
    # chemist could act on. The `any` rules further down need no equivalent -- an
    # extra solvate fragment cannot falsify "one of these is present".
    salt_ions = [x for x in counter if x not in SOLVATES]

    if name in _HALIDES:
        # A hydrohalide protonates a basic nitrogen, so the cation keeps its
        # hydrogen. A quaternary or pyridinium cation has none, and is a different
        # class however well the formula adds up.
        allowed = _HALIDES[name]
        return (
            bool(salt_ions)
            and all(x in allowed for x in salt_ions)
            and structure.cation_n_noh == 0
        )
    if name == "iodide":
        return bool(salt_ions) and all(x in _IODIDE for x in salt_ions)
    if name == "organic-bromide":
        return bool(salt_ions) and all(x in _ORGANIC_BROMIDE for x in salt_ions)
    if name == "quaternary-ammonium":
        return structure.cation_n_noh > 0
    if name == "maleate":
        # Maleate and fumarate are the same formula. Geometry is the only evidence,
        # and "undefined" is not treated as either: a maleate drawn with an
        # unspecified double bond is a defect for a chemist, not a maleate.
        #
        # `counter` is the "there is a counterion at all" test, as it is for the
        # sulfonate rule: the geometry is read from every fragment, so a lone
        # maleic acid record would otherwise satisfy the class on its own C=C.
        return (
            bool(salt_ions)
            and bool(structure.geoms)
            and all(g == "Z" for g in structure.geoms)
        )
    if name == "fumarate":
        return (
            bool(salt_ions)
            and bool(structure.geoms)
            and all(g == "E" for g in structure.geoms)
        )
    # `bool(salt_ions)` on all four, for the reason maleate, fumarate and sulfonate
    # already carry it: these scan every fragment, so a single-fragment record --
    # oxalic acid on its own, asserting ISA64148 -- satisfied its own class. And
    # nothing else caught it: `is_salt` is True because a RELATIONSHIP is set,
    # `has_counterion` is False so INT-03 asks for nothing, and CON-07 only fires
    # on a *non*-salt carrying a class. The record cleared.
    if name == "oxalate":
        return bool(salt_ions) and any(_uncharged(x) in _OXALATE for x in frags)
    if name == "tartrate":
        return bool(salt_ions) and any(_uncharged(x) in _TARTRATE for x in frags)
    if name == "sodium":
        return any(x in _SODIUM for x in counter)
    if name == "potassium":
        return any(x in _POTASSIUM for x in counter)
    if name == "sulfate":
        return bool(salt_ions) and any(x in _SULFATE for x in frags)
    if name == "methanesulfonate":
        return bool(salt_ions) and any(x in _METHANESULFONATE for x in frags)
    if name == "sulfonate":
        # `counter` being non-empty is the "there is a counterion at all" test; the
        # sulfonate itself is looked for across every fragment, because it may be
        # the largest one.
        return bool(salt_ions) and any(_is_sulfonate(x) for x in frags)
    return False


def unsupported_reason(name: str, structure: Structure) -> str:
    """Plain-English evidence for why a class is not supported, for the finding text."""
    if name in ("hydrochloride", "hydrobromide") and structure.cation_n_noh:
        return "quaternary/pyridinium cation, not a hydrohalide"
    # `counter` is every fragment except the largest, which on the shape this
    # module goes out of its way to handle -- a counterion bigger than its base --
    # is the *base*. Saying "counterions ['C5H11NO2', 'C5H11NO2']" to a chemist is
    # worse than saying nothing, so the fragments are named as what they are.
    return (
        f"fragments besides the largest {structure.counter}, geometry {structure.geoms}"
    )
