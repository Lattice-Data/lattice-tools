"""
Validate and repair CAS Registry Numbers before spending a request on them.

Shared by `chebi_lookup` and `structure_check`. It lives at the `bcp/` top level
rather than inside either package so neither has to import the other for a check
that belongs to both. Without it, `chebi_lookup.cas_to_cid()` and
`structure_check.client.cas_structure()` both send whatever the spreadsheet contained
straight to PubChem. A malformed number costs a request and comes back as "PubChem
does not know this", which reads identically to "this compound is not indexed" --
so the defect is attributed to the database instead of the cell.

**The load-bearing fact about the check digit: it proves the number is well-formed,
not that it is the right compound's number.** A valid checksum is necessary and
nowhere near sufficient, and two whole classes of corruption are invisible to it:

- **A spurious leading zero.** `0362-07-2` for `362-07-2`. Digits are weighted by
  position from the right, so a leading zero contributes `0 x weight` and shifts
  nothing -- the two strings produce the *identical* check digit. No arithmetic can
  separate them; only asking PubChem can.
- **A zero-padded check digit.** `42971-09-05` for `42971-09-5`. This one is worse
  than cosmetic, because it makes the string look like `YYYY-MM-DD` and a
  spreadsheet will then silently convert it to a date: `2113-05-05 00:00:00` is
  m-chlorophenylbiguanide hydrochloride's `2113-05-5`, two corruptions deep. The
  two are one causal chain, not two coincidences.

Two further shapes seen in real submissions, both repairable:

- **Rotated segments.** `3-4-7689` for `7689-03-4`. Offered only when the rotated
  reading passes its check digit, so a number that merely resembles the pattern is
  never mangled further.
- **Mojibake separators.** `864461?31?4`, where the hyphens did not survive an
  encoding round trip.

**One ambiguity is irreducible and is not guarded against.** Because a padded check
digit makes a CAS look like a date, the reverse is also true: a genuine date with a
single-digit day satisfies the CAS checksum about one time in ten, so
`2020-01-01 00:00:00` is repaired to `2020-01-1` and reported valid. No string rule
separates that from `2113-05-05 00:00:00`, which really is m-chlorophenylbiguanide
hydrochloride. A plausible-year test does not help -- the real compound's number *is*
year-shaped. The repair is applied because in a CAS column it is right far more often
than not, it is always reported in `repairs` so a curator can see it, and PubChem
resolution is what settles it. A date in a CAS column is a data-entry error either
way; this only means the error surfaces as an unresolvable number rather than as a
malformed one.

And one shape that is deliberately **not** repaired: a bare digit run such as
`6857789`. Reading it as `6857-78-9` is tempting and wrong -- such values are
PubChem **CIDs**, and `6857-78-9` is a real registry number belonging to an
unrelated compound. Guessing does not recover a lost number, it invents a wrong
one. `classify_cas()` returns `CAS_INVALID_FORMAT` and the row goes to a human.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Any

# A CAS Registry Number: two to seven digits, two digits, one check digit.
#
# `[0-9]` throughout, never `\d`: `\d` matches every Unicode decimal, and int()
# parses those happily, so an Arabic-Indic "\u0665\u0660-\u0660\u0660-\u0660"
# matched, satisfied the checksum and was reported CAS_VALID -- a string PubChem
# can never resolve. `\Z` rather than `$`, because `$` also matches before a
# trailing newline.
CAS_RE = re.compile(r"^([0-9]{2,7})-([0-9]{2})-([0-9])\Z")

# The Unicode dashes NFKC does not already fold, plus the soft hyphen, which is a
# format character but in a CAS column is a hyphen that survived a round trip as a
# discretionary break -- the same corruption as 864461?31?4. U+FF0D is in the
# class even though NFKC already folds it to ASCII hyphen: redundant, but the
# class is the documented set rather than the remainder after folding.
_SEPARATORS = re.compile(r"[‐-―−－?\u00ad\u2043\u02d7\u2e3a-\u2e3b]")
_EXCEL_TIME = re.compile(r"^(.*?)\s+00:00:00\Z")
# First group bounded like a real CAS body, so an 8-digit first segment is not
# "repaired" into something that can never be a registry number anyway.
_ZERO_PADDED_CHECK = re.compile(r"^([0-9]{2,7})-([0-9]{2})-0([0-9])\Z")
# The remainder must not itself begin with 0, or "007-00-1" becomes "07-00-1" --
# still leading-zero corrupted, while reporting the defect as repaired.
_LEADING_ZERO = re.compile(r"^0+([1-9][0-9]{1,6})-([0-9]{2})-([0-9])\Z")
_ROTATED = re.compile(r"^([0-9]{1,2})-([0-9]{1,2})-([0-9]{3,7})\Z")


def _has_leading_zero(cas: str) -> bool:
    """Whether the first segment still begins with 0 after every repair.

    A real CAS Registry Number never does. This is the residue the check digit
    cannot see: `00-00-0` sums to 0 and its check digit *is* 0, so an all-zero
    placeholder cell satisfies the checksum. `007-00-1` is the other shape --
    `_LEADING_ZERO` requires the remainder to start `[1-9]`, precisely so the repair
    does not turn it into `07-00-1` and call it fixed, so it arrives here untouched
    with both zeros still on the front. Without this check both are well-formed by
    the checksum and worth a request.
    """
    match = CAS_RE.match(cas or "")
    return bool(match) and match.group(1).startswith("0")


CAS_VALID = "valid"
CAS_INVALID_FORMAT = "invalid_format"
CAS_INVALID_CHECKSUM = "invalid_checksum"
CAS_MISSING = "missing"

CAS_CLASSES = (CAS_VALID, CAS_INVALID_FORMAT, CAS_INVALID_CHECKSUM, CAS_MISSING)

REPAIR_UNICODE_FOLD = "unicode_fold"
REPAIR_EXCEL_DATE = "excel_date"
REPAIR_SEPARATOR_MOJIBAKE = "separator_mojibake"
REPAIR_ZERO_PADDED_CHECK_DIGIT = "zero_padded_check_digit"
REPAIR_LEADING_ZERO = "leading_zero"
REPAIR_SEGMENT_ROTATION = "segment_rotation"


def cas_check_digit_ok(cas: str) -> bool | None:
    """
    Whether a CAS number's check digit is consistent with its body.

    Digits are taken right to left starting from the one before the check digit,
    weighted 1, 2, 3, ..., summed, and the total mod 10 must equal the check digit.

    Returns `None`, not `False`, when the string is not shaped like a CAS number at
    all -- "the question does not apply" and "the answer is no" are different
    findings, and collapsing them would let a blank cell be reported as a bad
    number.
    """
    match = CAS_RE.match(cas or "")
    if not match:
        return None
    body = match.group(1) + match.group(2)
    total = sum(int(d) * i for i, d in enumerate(reversed(body), start=1))
    return total % 10 == int(match.group(3))


def repair_rotation(raw: str) -> str:
    """
    Undo rotated segments: `"3-4-7689"` was `"7689-03-4"`.

    Returns `""` unless the rotated reading passes its check digit, so a number
    that merely resembles the pattern is never mangled further. Twelve of twelve
    real rows validated after this repair; twelve correct check digits by chance is
    about one in 10^12, which is what makes the rule safe to apply mechanically.
    """
    match = _ROTATED.match(raw or "")
    if not match:
        return ""
    first, second, third = match.groups()
    candidate = f"{third}-{int(first):02d}-{second}"
    return candidate if cas_check_digit_ok(candidate) else ""


def strip_leading_zero(raw: str) -> str:
    """
    Drop a spurious leading zero: `"0362-07-2"` -> `"362-07-2"`.

    Returns `""` when there is no leading zero to drop. The check digit cannot
    detect this defect -- both forms produce the same one -- so a caller that only
    validates will pass the broken value straight to PubChem, where it 404s. A real
    registry number has no leading zero in its first segment.

    `_LEADING_ZERO` is the single statement of what this accepts, and this is its
    single implementation -- including for `normalize_cas`, which calls it rather
    than matching the pattern itself. Two implementations diverge: a `CAS_RE` gate
    caps the first segment at seven digits and so declines `"01234567-07-2"`, a
    leading zero on a legitimate seven-digit body, which the pattern repairs.
    """
    match = _LEADING_ZERO.match(raw or "")
    if not match:
        return ""
    return f"{match.group(1)}-{match.group(2)}-{match.group(3)}"


def normalize_cas(raw: Any) -> tuple[str, str]:
    """
    Repair the mechanical corruptions of a CAS number.

    Returns `(cas, repairs)`, where `repairs` is a `+`-joined record of what was
    applied, empty when the value was already clean. Every repair is reported
    rather than applied silently: a curator has to be able to see that the value
    in the cell is not the value that was looked up.

    Repairs are recorded in the order they are applied, because they compose --
    a spreadsheet can only read a CAS number as a date *after* the check digit has
    been padded to two digits, and `0362-07-02` carries both a leading zero and a
    padded check digit.

    Segment rotation is deliberately excluded: it is a guess that has to be
    validated against the check digit, so it belongs in `classify_cas()`'s
    invalid-format branch via `repair_rotation()`, not in unconditional cleanup.
    """
    if raw is None or isinstance(raw, float) and raw != raw:  # NaN
        return "", ""
    text = (raw if isinstance(raw, str) else str(raw)).strip()
    if not text:
        return "", ""
    repairs: list[str] = []

    # Compatibility folding first, so a mojibake time suffix ("00\uff1a00\uff1a00")
    # is stripped by the step below rather than surviving the only pass that could
    # have removed it. Separator substitution is kept apart from the fold so the
    # two are reported honestly: a full-width digit is not a mangled hyphen.
    folded = unicodedata.normalize("NFKC", text)
    if folded != text:
        text = folded
        repairs.append(REPAIR_UNICODE_FOLD)

    match = _EXCEL_TIME.match(text)
    if match:
        text = match.group(1)
        repairs.append(REPAIR_EXCEL_DATE)

    separated = _SEPARATORS.sub("-", text)
    if separated != text:
        text = separated
        repairs.append(REPAIR_SEPARATOR_MOJIBAKE)

    match = _ZERO_PADDED_CHECK.match(text)
    if match:
        # Unconditional, including on a string that is a date rather than a CAS.
        # The date/CAS ambiguity is irreducible and is argued for in the module
        # docstring (the one-in-ten collision, `2020-01-01` vs `2113-05-05`); the
        # repair is reported in `repairs`.
        text = f"{match.group(1)}-{match.group(2)}-{match.group(3)}"
        repairs.append(REPAIR_ZERO_PADDED_CHECK_DIGIT)

    stripped = strip_leading_zero(text)
    if stripped:
        text = stripped
        repairs.append(REPAIR_LEADING_ZERO)

    return text, "+".join(repairs)


def classify_cas(raw: Any) -> tuple[str, str, str]:
    """
    Normalise a CAS number and say whether it is usable.

    Returns `(cas, cas_class, repairs)` with `cas_class` drawn from `CAS_CLASSES`:

    - `CAS_MISSING` -- the cell was empty. The name is the only handle.
    - `CAS_VALID` -- well-formed and the check digit agrees. Worth a request.
    - `CAS_INVALID_CHECKSUM` -- well-formed but the check digit disagrees. Usually a
      wholesale wrong number rather than a typo, so it is not worth guessing at.
    - `CAS_INVALID_FORMAT` -- no repair produced a registry number. Two shapes, and
      they return different things on purpose. Where nothing CAS-shaped came out at
      all, the value comes back **as the cell held it** with `repairs` empty: a
      repair that produced nothing usable is not information a curator needs, and
      nothing on that path is looked up, so there is no looked-up value that could
      differ from the cell. Where a repair *did* produce a CAS-shaped string that is
      still not a registry number -- a first segment that keeps a leading zero, such
      as `"0-1-007"` rotating to `"007-00-1"` -- the **repaired** value is returned
      and every repair is reported, because otherwise a curator gets a string that
      appears in no spreadsheet with no account of where it came from.

    `CAS_VALID` means *worth asking about*, never *correct*. The number can be
    perfectly formed and belong to an entirely different compound -- which is why
    the method that uses this resolves the compound's **name** independently and
    compares structures, rather than trusting a number that validates.

    `raw` is annotated `Any` deliberately: these values come from spreadsheet cells,
    so `None`, an `int` and `float("nan")` all arrive here and are all handled.
    """
    cas, repairs = normalize_cas(raw)
    if not cas:
        return "", CAS_MISSING, repairs
    # What the cell held, for the refusal path at the end. Recomputed rather than
    # returned from normalize_cas, which has two callers and one documented tuple.
    text = (raw if isinstance(raw, str) else str(raw)).strip()

    if CAS_RE.match(cas):
        if _has_leading_zero(cas):
            return cas, CAS_INVALID_FORMAT, repairs
        if cas_check_digit_ok(cas):
            return cas, CAS_VALID, repairs
        return cas, CAS_INVALID_CHECKSUM, repairs

    rotated = repair_rotation(cas)
    if rotated:
        # Rotation can expose a leading zero that normalize_cas never saw, because
        # in the rotated form it sat in the last segment: "3-4-07689" rotates to
        # "07689-03-4", which passes the check digit (leading zeros are invisible
        # to it) and would otherwise be reported CAS_VALID. The rotation label
        # comes first because the rotation is what exposed the leading zero --
        # repairs are recorded in the order they are applied.
        repairs = "+".join([r for r in (repairs, REPAIR_SEGMENT_ROTATION) if r])
        trimmed = strip_leading_zero(rotated)
        if trimmed:
            rotated = trimmed
            repairs = "+".join([r for r in (repairs, REPAIR_LEADING_ZERO) if r])
        # Reported on both exits, because the value returned is the rotated one
        # either way. Omitting the label on the refusal path handed a curator a
        # string that appears in no spreadsheet with no account of where it came
        # from -- exactly the silent repair normalize_cas() promises never to make.
        if _has_leading_zero(rotated):
            return rotated, CAS_INVALID_FORMAT, repairs
        return rotated, CAS_VALID, repairs

    # No repair produced anything CAS-shaped, so the cell is reported as it was and
    # no repair is claimed. Substituting separators before classifying means a cell
    # that is plainly a name comes through here -- "what?" folded to "what-",
    # "N-methylamine" from an em dash -- and reporting `separator_mojibake` on those
    # is noise in the one column whose job is auditability. Nothing is hidden by
    # dropping it: the promise is that a value which *was* looked up is never
    # silently different from the cell, and nothing on this path is looked up at all
    # (`_cas_for_lookup` refuses CAS_INVALID_FORMAT without a request).
    return text, CAS_INVALID_FORMAT, ""
