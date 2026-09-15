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
_IODIDE = ("I", "HI", "I-")
_ORGANIC_BROMIDE = ("Br", "HBr", "Br-")
_SODIUM = ("Na", "Na+")
_POTASSIUM = ("K", "K+")
_SULFATE = ("H2O4S", "HO4S-", "O4S-2")
_METHANESULFONATE = ("CH4O3S", "CH3O3S-")
_OXALATE_PREFIXES = ("C2H2O4", "C2HO4", "C2O4")
_TARTRATE_PREFIXES = ("C4H6O6", "C4H5O6", "C4H4O6")


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

    if name in _HALIDES:
        # A hydrohalide protonates a basic nitrogen, so the cation keeps its
        # hydrogen. A quaternary or pyridinium cation has none, and is a different
        # class however well the formula adds up.
        allowed = _HALIDES[name]
        return (
            bool(counter)
            and all(x in allowed for x in counter)
            and structure.cation_n_noh == 0
        )
    if name == "iodide":
        return bool(counter) and all(x in _IODIDE for x in counter)
    if name == "organic-bromide":
        return bool(counter) and all(x in _ORGANIC_BROMIDE for x in counter)
    if name == "quaternary-ammonium":
        return structure.cation_n_noh > 0
    if name == "maleate":
        # Maleate and fumarate are the same formula. Geometry is the only evidence,
        # and "undefined" is not treated as either: a maleate drawn with an
        # unspecified double bond is a defect for a chemist, not a maleate.
        return bool(structure.geoms) and all(g == "Z" for g in structure.geoms)
    if name == "fumarate":
        return bool(structure.geoms) and all(g == "E" for g in structure.geoms)
    if name == "oxalate":
        return any(x.startswith(_OXALATE_PREFIXES) for x in frags)
    if name == "tartrate":
        return any(x.startswith(_TARTRATE_PREFIXES) for x in frags)
    if name == "sodium":
        return any(x in _SODIUM for x in counter)
    if name == "potassium":
        return any(x in _POTASSIUM for x in counter)
    if name == "sulfate":
        return any(x in _SULFATE for x in frags)
    if name == "methanesulfonate":
        return any(x in _METHANESULFONATE for x in frags)
    if name == "sulfonate":
        return any(("S" in x and "O3" in x) or ("S2" in x and "O6" in x) for x in frags)
    return False


def unsupported_reason(name: str, structure: Structure) -> str:
    """Plain-English evidence for why a class is not supported, for the finding text."""
    if name in ("hydrochloride", "hydrobromide") and structure.cation_n_noh:
        return "quaternary/pyridinium cation, not a hydrohalide"
    return f"counterions {structure.counter}, geometry {structure.geoms}"
