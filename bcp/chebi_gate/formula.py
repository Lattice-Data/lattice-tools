"""Comparing a registry molecular formula against the stoichiometry actually drawn.

This is the arithmetic behind EXT-04, the only stoichiometry check whose evidence
does not come from the step that generated the SDFs, and therefore the only one
that can see a missing counterion or a missing solvate.

Three things make it more than string equality, and each is a bug that happened:

**Ratios reduce.** A sesquifumarate registered as ``C12H18N2O.3/2C4H4O4`` and drawn
as two base units with three fumarates are the same substance. Formulae are
compared as reduced whole-number element ratios, so 1:1.5 equals 2:3.

**An indefinite multiplier is not 1.** The registry writes
``C19H23N.xC4H4O4`` with ``(1:?)`` in the name when it declines to fix the
component ratio. Reading the ``x`` as 1 made three records report "formula agrees"
when nothing had been confirmed. The gate for that runs *before* any parsing,
because by the time a multiplier has been parsed the information that it was
indefinite is gone.

**A non-element is not an element.** The element pattern ``[A-Z][a-z]?`` matches
any capital-plus-optional-lowercase run, so a vendor code that has lost its space
-- ``WAG 994`` becoming ``WAG994`` -- parses as W, A, G994 and yields a confident,
meaningless formula. Symbols are checked against the periodic table, so garbage is
reported as unparseable rather than compared.

The hydrogen-only case is deliberately not a defect. A registry that lists the
neutral acid against a drawing that shows the ion differs by exactly one hydrogen
per counterion; 10 of 290 records on the reference batch were this, and calling
them defects would bury the three real disagreements.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd

from rdkit import Chem

log = logging.getLogger(__name__)

# Every element symbol RDKit knows, so a formula token that is not an element is a
# parse failure instead of a silently invented element.
_PERIODIC_TABLE = Chem.GetPeriodicTable()
ELEMENTS = frozenset(
    _PERIODIC_TABLE.GetElementSymbol(number) for number in range(1, 119)
)

# CAS exports carry HTML markup in formulae: C<sub>12</sub>H<sub>18</sub>, and
# <sup>3</sup>/<sub>2</sub> for a sesqui ratio.
_MARKUP = re.compile(r"</?su[bp]>")

# One trailing charge, as RDKit writes it on a whole-molecule formula: "C25H37O4-",
# "C6H5O7-3". Only one, and only at the end.
_TRAILING_CHARGE = re.compile(r"[+-]\d*$")

# A component's optional leading multiplier, whole or fractional.
_MULTIPLIER = re.compile(r"^(\d+/\d+|\d+)?(.*)$")

# An element symbol and its optional count.
_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")

# The registry declining to fix a component ratio. Three spellings, all seen:
# a letter multiplier in a dotted formula ("C19H23N.xC4H4O4"), and two ways of
# writing an unknown ratio in the registry name ("(1:?)", "(?:1)").
INDEFINITE = re.compile(r"\.[xn][A-Z]|\(\d+:\?\)|\(\?:")

# The registry's own marker, when the parsed table already carries the verdict.
RATIO_UNKNOWN_MARKER = "unknown-to-CAS"

AGREE = "agree"
CONVENTION = "convention"
DISAGREE = "DISAGREE"
RATIO_UNKNOWN = "ratio-unknown-to-CAS"
NO_RECORD = "no-registry-record"
UNPARSEABLE = "unparseable"

# Report order: worst first, so a summary reads top-down.
STATUS_ORDER = (DISAGREE, UNPARSEABLE, RATIO_UNKNOWN, NO_RECORD, CONVENTION, AGREE)


@dataclass(frozen=True)
class Comparison:
    """The verdict of one registry-formula against one drawn formula."""

    status: str
    note: str = ""
    registry_formula: str = ""
    drawn_formula: str = ""
    diff: dict[str, int] = field(default_factory=dict)

    @property
    def confirms(self) -> bool:
        """Whether the registry positively confirmed the drawn stoichiometry."""
        return self.status in (AGREE, CONVENTION)

    @property
    def checked(self) -> bool:
        """Whether a comparison happened at all.

        False means the record cannot be judged on this axis, which the report has
        to say rather than leave to look like a pass.
        """
        return self.status not in (NO_RECORD, RATIO_UNKNOWN, UNPARSEABLE)


def strip_markup(text: str) -> str:
    """Remove the HTML sub/sup markup a SciFinder export carries."""
    return _MARKUP.sub("", text or "").strip()


def is_indefinite(*texts: str) -> bool:
    """Whether any text says the registry does not fix the component ratio.

    Called on the formula *and* the registry name together, because the evidence
    appears in either: the formula carries ``.xC4H4O4`` and the name carries
    ``(1:?)``, and three records on the reference batch had both.
    """
    joined = " ".join(t or "" for t in texts)
    return bool(INDEFINITE.search(joined)) or RATIO_UNKNOWN_MARKER in joined


def parse(text: str) -> Counter | None:
    """A formula as element counts, honouring dotted components and multipliers.

    Returns None when there is nothing to parse or when a token is not an element
    symbol. Counts are :class:`~fractions.Fraction` so a ``3/2`` multiplier stays
    exact until :func:`reduce_ratio` turns the whole thing into whole numbers.
    """
    if not text:
        return None
    cleaned = _TRAILING_CHARGE.sub("", strip_markup(text))
    if not cleaned:
        return None

    total: Counter = Counter()
    for part in cleaned.split("."):
        if not part:
            continue
        match = _MULTIPLIER.match(part)
        if match is None:  # pragma: no cover - the pattern always matches
            return None
        multiplier = Fraction(match.group(1)) if match.group(1) else Fraction(1)
        body = match.group(2)
        if not body:
            continue
        consumed = 0
        for symbol, count in _ELEMENT.findall(body):
            if symbol not in ELEMENTS:
                log.debug("formula %r has non-element token %r", text, symbol)
                return None
            total[symbol] += multiplier * (int(count) if count else 1)
            consumed += len(symbol) + len(count)
        if consumed != len(body):
            # Characters the element pattern could not account for, e.g. a stray
            # bracket or a lowercase run. Better unparseable than partly read.
            log.debug("formula %r has unreadable text in component %r", text, part)
            return None
    return total or None


def reduce_ratio(counts: Counter | None) -> dict[str, int] | None:
    """Element counts as the smallest whole-number ratio, or None.

    This is what makes a sesquifumarate written 1:1.5 compare equal to one drawn
    2:3. Absolute size is deliberately discarded: the registry states a
    composition, not how many units someone chose to draw.
    """
    if not counts:
        return None
    denominator = 1
    for value in counts.values():
        denominator = (
            denominator * value.denominator // gcd(denominator, value.denominator)
        )
    integers = {symbol: int(value * denominator) for symbol, value in counts.items()}
    divisor = 0
    for value in integers.values():
        divisor = gcd(divisor, value)
    if divisor > 1:
        integers = {symbol: value // divisor for symbol, value in integers.items()}
    return integers


def compare(
    registry_formula: str | None,
    drawn_formula: str | None,
    *,
    registry_name: str = "",
    registry_ratio: str = "",
) -> Comparison:
    """Compare a registry formula against a drawn one, as reduced ratios.

    The order of the guards is the substance of this function. The indefinite
    check runs before any parsing, because ``x`` parses as nothing and the fact
    that a multiplier was indefinite cannot be recovered afterwards. That ordering
    is what three records on the reference batch needed: without it they came back
    "formula agrees" while the registry had confirmed no ratio at all.
    """
    if not registry_formula:
        return Comparison(
            status=NO_RECORD,
            note="no CAS registry row for this number, so stoichiometry is unchecked",
            drawn_formula=drawn_formula or "",
        )

    if is_indefinite(registry_formula, registry_name, registry_ratio):
        return Comparison(
            status=RATIO_UNKNOWN,
            note="the registry itself does not fix the component ratio",
            registry_formula=strip_markup(registry_formula),
            drawn_formula=drawn_formula or "",
        )

    left = reduce_ratio(parse(registry_formula))
    right = reduce_ratio(parse(drawn_formula or ""))
    clean = strip_markup(registry_formula)
    if not left or not right:
        which = "registry" if not left else "drawn"
        return Comparison(
            status=UNPARSEABLE,
            note=f"the {which} formula could not be parsed",
            registry_formula=clean,
            drawn_formula=drawn_formula or "",
        )

    if left == right:
        return Comparison(
            status=AGREE,
            note=f"registry formula {clean} agrees with the drawn structure "
            "(same ratio)",
            registry_formula=clean,
            drawn_formula=drawn_formula or "",
        )

    # Positive means the registry has more of that element than the drawing.
    diff = {
        symbol: left.get(symbol, 0) - right.get(symbol, 0)
        for symbol in set(left) | set(right)
        if left.get(symbol, 0) != right.get(symbol, 0)
    }
    if set(diff) == {"H"}:
        return Comparison(
            status=CONVENTION,
            note=f"registry formula {clean} differs from drawn {drawn_formula} by "
            f"{diff['H']:+d} H only (neutral-acid versus ionic drawing convention)",
            registry_formula=clean,
            drawn_formula=drawn_formula or "",
            diff=diff,
        )

    detail = ", ".join(f"{symbol}{value:+d}" for symbol, value in sorted(diff.items()))
    return Comparison(
        status=DISAGREE,
        note=f"registry formula {clean} disagrees with drawn {drawn_formula}: {detail}",
        registry_formula=clean,
        drawn_formula=drawn_formula or "",
        diff=diff,
    )
