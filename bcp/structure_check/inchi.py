"""
Compare two structures by reading their InChI layers, with no network calls.

The rest of this package compares **InChIKeys**: `skeleton()` takes the first 14
characters and `compare_structures()` tests full-key equality. That is cheap and
right for "same or not same", but an InChIKey is a hash, so when two keys differ it
cannot say *how*. `refine_skeleton_difference()` answers that by asking PubChem for
a desalted parent -- three requests per structure, and documented as unreliable for
multi-component salts.

An InChI **string** carries the same information directly and for free. The formula
layer distinguishes a salt from its free base; the `/b`, `/t`, `/m` and `/s` layers
distinguish stereoisomers. This module is that reading.

Why the distinction earns its own module: matching a compound to a ChEBI entry that
differs only in salt form is a curation error of a completely different kind from
matching it to the wrong stereoisomer, and a curator triages them differently. One
is "ChEBI has this compound but not this form"; the other is "we do not know which
molecule this row means".

`refine_skeleton_difference()` reaches this route because the resolvers carry the
InChI string alongside the key, at no extra request -- PubChem and ChEBI both return
the two in the same payload. A formula cannot be recovered from a 14-character
connectivity hash, so a caller holding only keys has nothing to read here, and the
parent-lookup route decides for a record with no InChI.

`defined_stereo()`, `defined_side()` and `is_multi_component()` are library surface
for callers outside this repo, so they have no production caller here.

**A comparison must declare the layers it is willing to ignore, not the ones it
remembers to look at.** That is the rule the whole module is built to, and every
refusal below follows from it: a layer outside `COMPARED_LAYERS`, a sublayer outside
`ISOTOPE_SUBLAYERS`, a non-standard version string and an unparseable formula all
return `LAYERS_NOT_COMPARABLE` rather than a verdict drawn from the layers that
happen to be recognised. `tests/test_structure_check_inchi_layers.py` pins each
tier against real compounds -- see `CHEBI_IDENTIFICATION.md` for how the method this
module came from arrived at them.
"""

from __future__ import annotations

import re

# Every InChI layer that carries stereochemistry.
#
#   /b  double-bond (E/Z) geometry
#   /t  tetrahedral parity, relative
#   /m  mirror flag -- which enantiomer of the /t arrangement
#   /s  stereo type (absolute, relative, racemic)
#
# Enumerated as a constant rather than checked ad hoc: all four are compared, and a
# comparison that reads only some of them cannot say so.
STEREO_LAYERS = ("b", "t", "m", "s")

# The layers that say which molecule this is, as opposed to which form or which
# isomer of it. `/c` is the connectivity table and `/h` the hydrogen positions;
# together they are the constitution, and a comparison that skips them calls
# constitutional isomers identical: ethanol and dimethyl ether are both C2H6O and
# differ only in `/c`, as do n-butane and isobutane, and catechol and hydroquinone.
# So does an InChI whose `/h` layer a cell limit truncated, against the intact
# string -- the shape a spreadsheet pipeline actually produces.
CONSTITUTION_LAYERS = ("c", "h")

# Charge and added/removed protons. These leave the formula layer as the neutral
# parent's, so they slip past a formula-equality gate: acetic acid and acetate
# differ only by `/p-1` and ChEBI holds them as separate entries.
IONISATION_LAYERS = ("q", "p")

# Isotopic substitution. A deuterated compound is its own ChEBI entry with its own
# CAS number, so it is neither a form nor a stereoisomer of the unlabelled one.
ISOTOPE_LAYER = "i"

# The sublayers a standard InChI may re-emit *inside* the isotopic block. They need
# their own set because `inchi_layers()` keeps everything from `/i` onwards under
# one key, so a layer arriving after `/i` never becomes a key and the
# COMPARED_LAYERS residue check below cannot see it. Without this set, `…/i1D3/x1`
# gets a confident verdict while the same `/x1` before `/i` is refused -- position
# deciding whether an unknown layer counts.
ISOTOPE_SUBLAYERS = frozenset("bhmst")

# Every layer this comparison reads. A layer outside this set is not ignored --
# `classify_pair()` refuses the pair, which is what makes the set enforced rather
# than declared. Non-standard InChI adds `/f` (fixed-H) and `/r` (reconnected
# metal), either of which can change which molecule the string denotes.
COMPARED_LAYERS = frozenset(
    CONSTITUTION_LAYERS + IONISATION_LAYERS + (ISOTOPE_LAYER,) + STEREO_LAYERS
)

# Verdicts. Finer-grained than this package's COMPARISON_VERDICTS, which cannot
# express "one side simply does not say" and cannot tell an isotopologue from a
# stereoisomer. `structure_check.client.comparison_verdict_from_inchi()` is the
# lossy projection onto the coarse vocabulary, and asserts at import time that it
# maps every verdict named here.
LAYERS_IDENTICAL = "layers_identical"
FORM_DIFFERS = "form_differs"
STEREO_DIFFERS = "stereo_differs"
STEREO_UNDEFINED_ON_ONE_SIDE = "stereo_undefined_on_one_side"
ISOTOPE_DIFFERS = "isotope_differs"
DIFFERENT_COMPOUND = "different_compound"
LAYERS_NOT_COMPARABLE = "not_comparable"

INCHI_VERDICTS = (
    LAYERS_IDENTICAL,
    FORM_DIFFERS,
    STEREO_DIFFERS,
    STEREO_UNDEFINED_ON_ONE_SIDE,
    ISOTOPE_DIFFERS,
    DIFFERENT_COMPOUND,
    LAYERS_NOT_COMPARABLE,
)

# A formula layer: optional stoichiometric multiplier, then element symbols and
# counts, with "." between components. A string that merely contains a "/" is not a
# formula layer: without this, "a/b" against "x/b" compares identical.
#
# Every digit run is bounded, because `int()` raises past 4300 digits and this
# module must return LAYERS_NOT_COMPARABLE rather than raise. Six digits is four
# orders of magnitude past any real stoichiometry. Modelled as "optional
# multiplier, then one or more (element symbol, optional count)" rather than as a
# character class, because a bound inside a repeated alternation is no bound at
# all: `(?:[A-Za-z]|[0-9]{0,6})*` matches 4400 digits as several shorter runs.
_COUNT = "[0-9]{0,6}"
_COMPONENT = rf"{_COUNT}(?:[A-Z][a-z]?{_COUNT})+"
_FORMULA_RE = re.compile(rf"^{_COMPONENT}(?:\.{_COMPONENT})*\Z")

# The version field of a **standard** InChI. The layer set above is complete for
# this version only, so `classify_pair()` requires it rather than assuming it: a
# non-standard string parses cleanly through `inchi_layers()` and would otherwise
# be handed a confident verdict drawn from an incomplete reading.
_STANDARD_VERSION = "InChI=1S"

# Bounded here too, so `principal_component()` is safe called directly with a
# formula that never went through _FORMULA_RE.
_ELEMENT = re.compile(rf"([A-Z][a-z]?)({_COUNT})")
_LEADING_COUNT = re.compile(r"^[0-9]+")


def inchi_layers(inchi: str) -> tuple[dict[str, str], str]:
    """
    Split an InChI into its layers.

    Returns `(layers, formula)`, where `layers` maps a single-letter prefix to the
    rest of that layer and `formula` is the formula layer. An InChI is
    `InChI=1S/<formula>/<layer>/<layer>/...`, so the formula is element 1 after
    splitting on "/" and the layers follow.

    Returns `({}, "")` for empty or malformed input rather than raising: every
    caller in this package is on a path that promises not to raise, and a
    structure we cannot parse must read as *unknown*, never as *equal*.

    **Everything from `/i` onwards is one value under `"i"`, not separate layers.**
    An InChI re-emits `/h`, `/t`, `/m` and `/s` after `/i` to give the *isotopic*
    hydrogen positions and stereochemistry, in the same single-letter namespace as
    the ordinary ones. The two namespaces must not merge, in either direction:

    - Heavy water and tritiated water differ only in that block -- `/i/hD2` against
      `/i/hT2`, both carrying an ordinary `/h1H2` first -- so letting the ordinary
      `/h` win drops the only difference there is and reports a false `match`.
    - A structure whose only `/t` is the isotopic one must not have that value read
      as its ordinary stereochemistry, which invents a stereo difference against a
      structure that really specifies one.

    Keeping the block whole in `layers["i"]` does both, and keeps every key a single
    letter. It is `ISOTOPE_SUBLAYERS` rather than `COMPARED_LAYERS` that then covers
    the residue check inside it.
    """
    # Stripped: a trailing newline from a text file or a CSV cell would otherwise
    # land in the final layer's value and make a structure differ from itself.
    parts = inchi.strip().split("/") if inchi else []
    layers: dict[str, str] = {}
    for index, part in enumerate(parts[2:]):
        if not part:
            continue
        if part[0] == ISOTOPE_LAYER:
            # The rest of the string, rejoined as it arrived, so the isotopic block
            # is compared as a unit by the isotope tier.
            layers[ISOTOPE_LAYER] = "/".join(parts[2 + index :])[1:]
            break
        # setdefault rather than assignment for the ordinary layers too: a
        # malformed string repeating a prefix must not have the later value win.
        layers.setdefault(part[0], part[1:])
    return layers, (parts[1] if len(parts) > 1 else "")


def _isotope_residue(layers: dict[str, str]) -> set[str]:
    """Prefixes inside the isotopic block that this module does not recognise.

    The block is `<substitution spec>[/<sublayer>…]`, and the spec carries no prefix
    letter of its own -- a leading `/` means there is no spec at all, as in
    `…/i/hD2`. So the sublayers are everything after the first `/`.
    """
    block = layers.get(ISOTOPE_LAYER)
    if not block:
        return set()
    return {part[0] for part in block.split("/")[1:] if part} - ISOTOPE_SUBLAYERS


def inchi_version(inchi: str) -> str:
    """
    The version field of an InChI string -- `"InChI=1S"` for standard -- or `""`.

    Re-splits rather than widening `inchi_layers()`'s documented two-tuple: one
    gate in one function needs this, and every other caller would have had to grow
    a third name it never reads.
    """
    parts = inchi.strip().split("/") if inchi else []
    return parts[0] if parts else ""


def defined_stereo(layers: dict[str, str]) -> int:
    """
    How many stereo layers this structure actually specifies.

    A measure of how specific a structure is, for reporting. **Not** the way to pick
    an anchor: use `defined_side()`, which tests containment of the defined-layer
    sets. A count makes a side defining only `/b` and a side defining only `/t` look
    equally specific, and whichever was passed first then won -- the tie that
    containment exists to refuse.

    All four layers are counted, `/b` included: a count over `/t` and `/m` alone
    reads a pair whose only difference is `/b` on one side as equally specific --
    both absent on both sides -- and an anchor chosen from that tie matches a
    stereo-unspecified database entry to a compound whose geometry is specified.
    """
    return sum(1 for k in STEREO_LAYERS if layers.get(k) is not None)


def _atom_counts(component: str) -> tuple[int, int]:
    """`(total atoms, heavy atoms)` for one component of a formula layer.

    The stoichiometric multiplier needs no stripping here: `_ELEMENT` only matches
    from `[A-Z]`, so `findall` skips a leading digit run either way. It does have to
    be stripped where components are *compared*, which is what the `_LEADING_COUNT`
    calls in `_rank()` and `classify_pair()` are for.
    """
    total = heavy = 0
    for element, count in _ELEMENT.findall(component):
        n = int(count) if count else 1
        total += n
        if element != "H":
            heavy += n
    return total, heavy


def _rank(component: str) -> tuple[int, int, str]:
    """`principal_component`'s ordering: heavy atoms, then total atoms, then text.

    The multiplier is dropped from the text tiebreak only, so `"2C26H28N3"` ranks on
    the atoms of one unit and ties break on the formula rather than on the count.
    """
    total, heavy = _atom_counts(component)
    return heavy, total, _LEADING_COUNT.sub("", component)


def principal_component(formula: str) -> str:
    """
    The largest component of a multi-component formula layer.

    `"C27H29NO11.ClH"` -> `"C27H29NO11"`. Two structures sharing this differ only
    by counterion or water. The stoichiometric multiplier is kept in the returned
    string but ignored when choosing, so ranking is on the atoms of one unit. That
    makes pyrvinium a heavy-atom tie -- 29 against pamoate's 29 -- decided by the
    hydrogen-inclusive total, 57 against 45. Counting the multiplier would decide
    it 58 to 29 and pick the same component; it is ignored because a component's
    identity is its formula, not how many equivalents of it the salt carries.

    Components are ranked by **heavy-atom count, then total atoms**, not by the
    length of the formula string and not by total atoms alone. Length ties on
    exactly the compound this is most often needed for:
    pyrvinium pamoate is `2C26H28N3.C23H16O6`, whose two components are both eight
    characters once the multiplier is dropped, so a string tiebreak returned the
    *pamoate counterion* -- the same wrong answer PubChem's desalting gives, which
    is the reason for preferring this route in the first place. Atom counts are
    57 against 45 and pick the drug. Across 407 multi-component structures from a
    real curation run the two rankings disagree on three.

    Total atoms alone is not enough either, because hydrogen then outvotes
    everything: diatrizoate meglumine is `C11H9I3N2O4.C7H17NO5`, where the
    hydrogen-rich meglumine counterion wins on total atoms 30 to 29 and the
    iodinated parent wins on heavy atoms 20 to 13. Heavy atoms decide first and
    total atoms break the tie, which is what keeps pyrvinium correct as well.

    **Limit, stated precisely:** this is only meaningful for an organic salt or
    hydrate, where one component is the compound and the rest are counterions or
    water. For a coordination complex or an organometallic that InChI has split
    into fragments -- `6ClH.14H3N.2O.3Ru` for ruthenium red, `C6H5.C2H4O2.Hg` for
    phenylmercuric acetate -- no component is the parent and the answer is
    meaningless whichever key is used. `classify_pair()` is unaffected, because it
    applies the same ranking to both sides and only ever compares the results.
    """
    components = [c for c in (formula or "").split(".") if c]
    if not components:
        return ""
    return max(components, key=_rank)


def classify_pair(inchi_a: str, inchi_b: str) -> str:
    """
    How two structures stand to one another, from their layers alone.

    Returns one of `INCHI_VERDICTS`:

    - `LAYERS_NOT_COMPARABLE` -- one or both structures are missing or unparseable,
      or either carries a layer outside `COMPARED_LAYERS` (or a sublayer outside
      `ISOTOPE_SUBLAYERS` inside the isotopic block), or **either** is not standard
      InChI. Either, not neither: a standard string against an `InChI=1` one is
      refused, and that mixed pair is the likely one here -- ChEBI's
      `standard_inchi` against a vendor or legacy string. The last two conditions
      are the same refusal as the first: a verdict reached by quietly skipping a
      layer whose meaning this module does not know would be the defect it exists
      to remove.
    - `LAYERS_IDENTICAL` -- same formula, same value in every stereo layer.
    - `FORM_DIFFERS` -- two ways in, and they are different questions with one
      answer. Either the principal components match and the **formula layers
      differ**, which is a salt, a hydrate, or a different stoichiometry of the same
      cation; or the formulas are *equal*, the constitution matches, and a `/q` or
      `/p` layer differs, which is a charge or protonation state -- acetic acid
      against acetate. Both are "the same compound in another form", which is why
      they share a verdict, but a caller must not read `FORM_DIFFERS` as implying the
      formula layers differ.

      The first route returns before the constitution, isotope and stereo
      comparisons run, so a pair that differs in **both** form and stereochemistry
      reports only the form difference. The tier order is forced -- a salt's extra
      component also changes `/c`, so the formula relationship has to be settled
      first -- but it means one call names one kind of disagreement, not every kind.
      A caller needing both must compare the stereo layers itself.

      **Two pairs get this wrong**, because the verdict is reached as soon as the
      stripped principal components match, so it is exact only where that shared
      component is the same molecule on both sides *and* is the compound rather than
      the counterion:

      1. *Two unrelated drugs sharing a heavy counterion.* `principal_component()`
         picks the counterion when it is the largest fragment, so hydroxyzine pamoate
         and pyrantel pamoate share one. Not fixable from the formula: `{drug, ClH}`
         against `{drug, BrH}` and `{hydroxyzine, pamoate}` against
         `{pyrantel, pamoate}` are one shape -- one shared component, one differing
         component per side -- leaving size as the only discriminator, and size is
         what failed. Counterions run from chloride's 1 heavy atom to pamoate's 29,
         overlapping the drugs (pyrantel 14, hydroxyzine 26). A subset-relation rule
         looks like the fix and is not: it reports a hydrochloride against a
         hydrobromide of the *same* drug as `DIFFERENT_COMPOUND`, the false positive
         this module exists to remove.
      2. *Two constitutional isomers carrying different counterions.*
         `C21H27ClN2O2.ClH` against an isomer's `C21H27ClN2O2.BrH` matches on the
         principal component's **formula** while the molecules differ in `/c`, which
         the formula tier returns before reading. With the *same* counterion the pair
         reaches the constitution tier and is correctly `DIFFERENT_COMPOUND`, so only
         the counterion difference hides it.

      Case 2 is the one the per-component `/c` sublayers could decide, and aligning
      them is harder than it looks: `/c` enumerates component **instances** and
      compresses runs with `N*`, so apomorphine hemihydrate's three formula
      components `2C17H17NO2.2ClH.H2O` give four `;`-separated sublayers standing for
      five instances. A mistake there corrupts `FORM_DIFFERS` for every ordinary
      salt. Case 1 needs to know which fragment is the drug, and is not decidable
      here at all.

      The direction is safe for both -- `salt_differs` surfaces the row as `check`
      rather than `ok`, so nothing is waved through, though it is a systematic move
      out of `investigate`. `test_two_salts_of_the_same_counterion` and
      `test_two_isomers_with_different_counterions` pin them, so neither reads as an
      oversight.
    - `STEREO_DIFFERS` -- same formula, and some stereo layer is present on **both**
      sides with **different** values. A genuine stereoisomer difference.
    - `STEREO_UNDEFINED_ON_ONE_SIDE` -- same formula, and every differing layer is
      absent on one side. PubChem holds both a stereo-defined and a
      stereo-undefined record for many substances, so this is a cataloguing
      artifact rather than a disagreement; anchor on `defined_side()`, which is
      containment of the defined-layer sets rather than a count of them.
    - `DIFFERENT_COMPOUND` -- the principal components differ.

    The `STEREO_DIFFERS` / `STEREO_UNDEFINED_ON_ONE_SIDE` split is the whole point.
    "Absent on one side" is safe to treat as agreement; "present on both and
    disagreeing" is two different molecules and must never be waved through.

    Missing input returns `LAYERS_NOT_COMPARABLE`, and in particular two empty
    strings are **not** identical. A caller may guard its own call site, but a
    library cannot rely on that: reporting silence from a database as agreement
    between two structures is the most damaging thing a comparison here can do.
    """
    layers_a, formula_a = inchi_layers(inchi_a)
    layers_b, formula_b = inchi_layers(inchi_b)
    if not _FORMULA_RE.match(formula_a or "") or not _FORMULA_RE.match(formula_b or ""):
        return LAYERS_NOT_COMPARABLE
    # Both gates below refuse rather than ignore. The tiers that follow read a fixed
    # list of layers, which is complete for standard InChI and for nothing else, so
    # a non-standard string or an unrecognised layer must not be handed a verdict
    # drawn from the layers that happen to be recognised.
    versions = (inchi_version(inchi_a), inchi_version(inchi_b))
    if any(version != _STANDARD_VERSION for version in versions):
        return LAYERS_NOT_COMPARABLE
    if (set(layers_a) | set(layers_b)) - COMPARED_LAYERS:
        return LAYERS_NOT_COMPARABLE
    # The isotopic block is one value, so the check above cannot see inside it.
    if _isotope_residue(layers_a) | _isotope_residue(layers_b):
        return LAYERS_NOT_COMPARABLE

    # Tiered, coarsest first, so a difference is never misfiled as the wrong *kind*
    # of difference. The order is forced: a salt's extra component also changes /c,
    # so the formula relationship must be settled before the constitution is
    # compared, or every salt pair would report as two different molecules.
    principal_a = _LEADING_COUNT.sub("", principal_component(formula_a))
    principal_b = _LEADING_COUNT.sub("", principal_component(formula_b))
    if not principal_a or principal_a != principal_b:
        return DIFFERENT_COMPOUND
    if formula_a != formula_b:
        return FORM_DIFFERS

    # Formulas now equal, so a /c or /h difference is a constitutional isomer, a
    # tautomer, or a truncated string -- all three "not the same molecule".
    if any(layers_a.get(k) != layers_b.get(k) for k in CONSTITUTION_LAYERS):
        return DIFFERENT_COMPOUND

    if any(layers_a.get(k) != layers_b.get(k) for k in IONISATION_LAYERS):
        return FORM_DIFFERS
    if layers_a.get(ISOTOPE_LAYER) != layers_b.get(ISOTOPE_LAYER):
        return ISOTOPE_DIFFERS

    for layer in STEREO_LAYERS:
        a, b = layers_a.get(layer), layers_b.get(layer)
        if a != b and a is not None and b is not None:
            return STEREO_DIFFERS
    if any(layers_a.get(k) != layers_b.get(k) for k in STEREO_LAYERS):
        return STEREO_UNDEFINED_ON_ONE_SIDE
    return LAYERS_IDENTICAL


def defined_side(inchi_a: str, inchi_b: str) -> str:
    """
    Of two structures, the one specifying more stereochemistry.

    "More specific" means its set of defined stereo layers is a strict **superset**,
    not merely larger. A side defining only `/b` and a side defining only `/t` both
    define one layer and neither dominates; counting would have picked whichever was
    passed first. Such ties, and exact ties, return `inchi_a`. A caller using this
    as an identification anchor should treat a tie as no preference: an
    exact-structure match against the returned side then fails to match any ChEBI
    entry, which is the safe direction.
    """
    layers_a, _ = inchi_layers(inchi_a)
    layers_b, _ = inchi_layers(inchi_b)
    set_a = {k for k in STEREO_LAYERS if layers_a.get(k) is not None}
    set_b = {k for k in STEREO_LAYERS if layers_b.get(k) is not None}
    # Containment, not a count: a side defining {b} and a side defining {t} both
    # count 1, and neither is more specific than the other. Only a strict superset
    # is grounds for preferring one.
    if set_b > set_a:
        return inchi_b
    return inchi_a


def is_multi_component(inchi: str) -> bool | None:
    """
    Whether the formula layer has more than one component -- a salt or a hydrate.

    ChEBI registers a salt as its own entry, separate from the free base, so this
    is what decides whether a compound absent from ChEBI needs a new entry or only
    a salt added to an existing parent.

    Returns `None` on exactly one input: a missing string, or a formula layer that
    does not parse. `False` there would report "I could not read this" as "one
    component" -- the collapse of unknown onto an answer that the rest of this
    module refuses -- and it would send a curator the wrong way round, to a new
    ChEBI entry rather than a salt on an existing parent.

    **A narrower refusal than `classify_pair()`'s, deliberately.** That function also
    refuses non-standard InChI and any layer outside `COMPARED_LAYERS`, because it
    reads a fixed list of layers and one it does not recognise could change which
    molecule the string denotes. This function reads the formula layer only, and the
    formula layer means the same thing in `InChI=1` as in `InChI=1S` and is not
    changed by an `/f` or `/r` layer sitting after it -- so `InChI=1/C27H29NO11.ClH/…`
    answers `True` here and `LAYERS_NOT_COMPARABLE` there, and both are right.
    """
    _, formula = inchi_layers(inchi)
    if not _FORMULA_RE.match(formula or ""):
        return None
    return "." in formula
