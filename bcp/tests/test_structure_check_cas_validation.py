"""Unit tests for structure_check.cas validation and repair (no network)."""

from __future__ import annotations

import pytest

from structure_check.cas import (
    CAS_INVALID_CHECKSUM,
    CAS_INVALID_FORMAT,
    CAS_MISSING,
    CAS_VALID,
    REPAIR_EXCEL_DATE,
    REPAIR_LEADING_ZERO,
    REPAIR_SEGMENT_ROTATION,
    REPAIR_SEPARATOR_MOJIBAKE,
    REPAIR_ZERO_PADDED_CHECK_DIGIT,
    cas_check_digit_ok,
    classify_cas,
    normalize_cas,
    repair_rotation,
    strip_leading_zero,
)

# Real registry numbers, each load-bearing for one rule.
ETHANOL = "64-17-5"
CAMPTOTHECIN = "7689-03-4"
M_CPBG_HCL = "2113-05-5"  # m-chlorophenylbiguanide hydrochloride
TWO_METHOXYESTRADIOL = "362-07-2"
VINPOCETINE = "42971-09-5"
BENZOXIQUINE = "86-75-9"
DINACICLIB_WRONG = "779353-01-3"  # off by one in the check digit; real is ...-01-4


# --------------------------------------------------------------------------
# The check digit
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cas", [ETHANOL, CAMPTOTHECIN, M_CPBG_HCL, VINPOCETINE, BENZOXIQUINE]
)
def test_real_registry_numbers_pass_their_check_digit(cas: str) -> None:
    assert cas_check_digit_ok(cas) is True


def test_a_wrong_check_digit_is_detected() -> None:
    assert cas_check_digit_ok(DINACICLIB_WRONG) is False


@pytest.mark.parametrize("raw", ["", "not a cas", "6857789", "1-2-3-4"])
def test_check_digit_is_none_when_the_question_does_not_apply(raw: str) -> None:
    """Not-shaped-like-a-CAS is a different finding from a failing check digit.

    Collapsing the two would report an empty cell as a bad registry number.
    """
    assert cas_check_digit_ok(raw) is None


def test_a_leading_zero_does_not_change_the_check_digit() -> None:
    """Why PubChem is the only thing that can catch this defect.

    Digits are weighted by position from the right, so a leading zero contributes
    zero and shifts nothing. The two strings are arithmetically indistinguishable.
    """
    assert cas_check_digit_ok("0362-07-2") == cas_check_digit_ok("362-07-2") is True


# --------------------------------------------------------------------------
# Mechanical repairs
# --------------------------------------------------------------------------


def test_spreadsheet_read_the_cas_as_a_date() -> None:
    """`2113-05-05 00:00:00` is `2113-05-5`, two corruptions deep.

    The date reading only becomes possible after the check digit has been padded to
    two digits, which is what makes the string look like `YYYY-MM-DD`.
    """
    cas, repairs = normalize_cas("2113-05-05 00:00:00")
    assert cas == M_CPBG_HCL
    assert REPAIR_EXCEL_DATE in repairs
    assert REPAIR_ZERO_PADDED_CHECK_DIGIT in repairs


def test_zero_padded_check_digit_is_unpadded() -> None:
    cas, repairs = normalize_cas("42971-09-05")
    assert cas == VINPOCETINE
    assert repairs == REPAIR_ZERO_PADDED_CHECK_DIGIT


def test_leading_zero_and_padded_check_digit_compose() -> None:
    cas, repairs = normalize_cas("0362-07-02")
    assert cas == TWO_METHOXYESTRADIOL
    assert REPAIR_ZERO_PADDED_CHECK_DIGIT in repairs
    assert REPAIR_LEADING_ZERO in repairs


def test_mojibake_separators_are_restored() -> None:
    cas, repairs = normalize_cas("86?75?9")
    assert cas == BENZOXIQUINE
    assert repairs == REPAIR_SEPARATOR_MOJIBAKE


def test_a_clean_number_is_left_alone_and_reports_no_repair() -> None:
    assert normalize_cas(ETHANOL) == (ETHANOL, "")
    assert normalize_cas(f"  {ETHANOL}  ") == (ETHANOL, "")


def test_rotated_segments_are_repaired_only_when_the_result_validates() -> None:
    assert repair_rotation("3-4-7689") == CAMPTOTHECIN
    assert repair_rotation("1-1-1111") == ""
    assert repair_rotation(ETHANOL) == ""


def test_strip_leading_zero_repairs_or_declines() -> None:
    """All three branches, including the one that declines a real leading zero.

    `05-00-0` does begin with a zero, but stripping it leaves a single-digit first
    segment, which is not CAS-shaped -- so the repair is refused rather than
    producing a malformed number. That branch had no coverage.
    """
    assert strip_leading_zero("0362-07-2") == TWO_METHOXYESTRADIOL
    assert strip_leading_zero("09-17-5") == ""  # len(trimmed) < 2
    assert strip_leading_zero("00-17-5") == ""  # not trimmed
    assert strip_leading_zero("007-00-1") == ""  # a zero survives the strip
    assert strip_leading_zero(ETHANOL) == ""
    assert strip_leading_zero("") == ""


def test_a_rotation_that_can_only_yield_a_leading_zero_is_refused() -> None:
    """`0-1-007` rotates to a first segment that is all zeros but one.

    The rotation gate passes it, because the check digit cannot see leading zeros;
    the final guard is what refuses it.

    All three fields, not just the class: the value returned is the **rotated**
    string, so omitting the label handed a curator `007-00-1` -- a number that
    appears in no spreadsheet -- with no record of how it got there.
    """
    assert classify_cas("0-1-007") == (
        "007-00-1",
        CAS_INVALID_FORMAT,
        REPAIR_SEGMENT_ROTATION,
    )


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------


def test_a_clean_number_classifies_valid() -> None:
    assert classify_cas(ETHANOL) == (ETHANOL, CAS_VALID, "")


def test_an_empty_cell_is_missing_not_invalid() -> None:
    assert classify_cas("") == ("", CAS_MISSING, "")
    assert classify_cas("   ") == ("", CAS_MISSING, "")


def test_a_failing_check_digit_classifies_invalid_checksum() -> None:
    cas, cas_class, _ = classify_cas(DINACICLIB_WRONG)
    assert (cas, cas_class) == (DINACICLIB_WRONG, CAS_INVALID_CHECKSUM)


def test_a_repaired_rotation_classifies_valid_and_says_so() -> None:
    cas, cas_class, repairs = classify_cas("3-4-7689")
    # Exact, not `in`: a substring assertion keeps passing if classify_cas drops
    # the repairs it inherited from normalize_cas and reports only its own.
    assert (cas, cas_class, repairs) == (
        CAMPTOTHECIN,
        CAS_VALID,
        REPAIR_SEGMENT_ROTATION,
    )


def test_a_bare_digit_run_is_never_guessed_into_a_cas() -> None:
    """`6857789` is a PubChem CID, and `6857-78-9` is a different compound.

    Something upstream once made exactly this guess, which is why the value must
    reach a human as unusable rather than as a plausible registry number.
    """
    cas, cas_class, _ = classify_cas("6857789")
    assert cas_class == CAS_INVALID_FORMAT
    assert cas == "6857789"


def test_repairs_are_reported_through_classification() -> None:
    """A curator has to see that the value looked up is not the value in the cell."""
    cas, cas_class, repairs = classify_cas("0362-07-02")
    assert (cas, cas_class) == (TWO_METHOXYESTRADIOL, CAS_VALID)
    assert repairs == (f"{REPAIR_ZERO_PADDED_CHECK_DIGIT}+{REPAIR_LEADING_ZERO}")


# --------------------------------------------------------------------------
# Defects the check digit cannot see, found by fuzzing the repairs
# --------------------------------------------------------------------------


def test_rotation_can_expose_a_leading_zero_and_it_is_repaired() -> None:
    """`3-4-07689` rotates to `07689-03-4`, which passes its check digit.

    Leading zeros are invisible to the checksum, so returning the rotated string
    unexamined reported a number PubChem cannot resolve as CAS_VALID. The rotated
    form has to go back through the leading-zero repair.
    """
    cas, cas_class, repairs = classify_cas("3-4-07689")
    assert (cas, cas_class) == (CAMPTOTHECIN, CAS_VALID)
    # Both labels, in order, exactly: asserting membership of each would pass if
    # classify_cas truncated the joined record to its first element.
    assert repairs == f"{REPAIR_LEADING_ZERO}+{REPAIR_SEGMENT_ROTATION}"


@pytest.mark.parametrize("raw", ["00-00-0", "0000-00-0", "007-00-1", "0-00-0"])
def test_a_first_segment_still_starting_with_zero_is_not_a_registry_number(
    raw: str,
) -> None:
    """`00-00-0` sums to zero and its check digit *is* zero, so it validates.

    A zero-filled or placeholder cell was therefore classed as well-formed and
    worth a request. A real CAS Registry Number has no leading zero.
    """
    _, cas_class, _ = classify_cas(raw)
    assert cas_class == CAS_INVALID_FORMAT


@pytest.mark.parametrize("raw", ["٥٠-٠٠-٠", "۵۰-۰۰-۰"])
def test_non_ascii_digits_are_not_a_valid_cas(raw: str) -> None:
    """`\\d` matches every Unicode decimal and int() parses them.

    Arabic-Indic and Extended Arabic-Indic digits therefore matched the pattern and
    satisfied the checksum, yielding CAS_VALID for a string PubChem can never
    resolve. The pattern uses `[0-9]`.
    """
    _, cas_class, _ = classify_cas(raw)
    assert cas_class == CAS_INVALID_FORMAT


def test_full_width_digits_are_folded_rather_than_rejected() -> None:
    """NFKC folds these to ASCII, so the number is recoverable -- unlike the above."""
    cas, cas_class, repairs = classify_cas("６４-１７-５")
    assert (cas, cas_class) == (ETHANOL, CAS_VALID)
    assert repairs


def test_a_trailing_newline_does_not_pass_as_a_cas() -> None:
    """`$` matches before a final newline; the pattern anchors with `\\Z`."""
    assert cas_check_digit_ok("64-17-5\n") is None
    assert classify_cas("64-17-5\n")[1] == CAS_VALID  # classify strips first


def test_an_eight_digit_first_segment_is_not_partially_repaired() -> None:
    """No CAS body is that long, so reporting a repair on it is noise."""
    cas, cas_class, repairs = classify_cas("12345678-12-01")
    assert cas_class == CAS_INVALID_FORMAT
    assert repairs == ""
    assert cas == "12345678-12-01"


def test_a_mojibake_time_suffix_is_still_stripped() -> None:
    """Compatibility folding has to run before the time suffix is matched.

    Otherwise full-width colons survive the only pass that could remove them and
    the value is reported malformed with a separator repair claimed.
    """
    cas, cas_class, _ = classify_cas("2113-05-05 00：00：00")
    assert (cas, cas_class) == (M_CPBG_HCL, CAS_VALID)


@pytest.mark.parametrize("raw", [6857789, 64175, None, float("nan")])
def test_non_string_input_classifies_rather_than_crashing(raw: object) -> None:
    """A spreadsheet reader hands back int64 and float('nan'), not str."""
    cas, cas_class, _ = classify_cas(raw)  # type: ignore[arg-type]
    assert cas_class in (CAS_INVALID_FORMAT, CAS_MISSING)


def test_a_genuine_date_can_still_pass_and_that_is_documented() -> None:
    """The one ambiguity this module cannot resolve, pinned so it is not a surprise.

    A padded check digit makes a CAS look like a date, so a real date with a
    single-digit day satisfies the checksum roughly one time in ten. No string rule
    separates `2020-01-01` from `2113-05-05`, which genuinely is
    m-chlorophenylbiguanide hydrochloride's number. The repair is always reported,
    and PubChem is what settles it.
    """
    cas, cas_class, repairs = classify_cas("2020-01-01 00:00:00")
    assert (cas, cas_class) == ("2020-01-1", CAS_VALID)
    assert REPAIR_EXCEL_DATE in repairs
