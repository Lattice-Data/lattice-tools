"""Byte fidelity of SDF parsing and re-emission.

These tests exist because the prototype this was ported from failed the invariant
they pin. Its cleared file was written as ``chunk.rstrip("\\n") + "\\n"``, and every
record in the real submission files ends with a blank line, so every cleared record
came back one byte short -- measured 0 of 13 byte-identical on the 95-record novel
file. Nothing caught it, because the only assertion made about the cleared file was
that it contained no ``GATE_*`` field.
"""

from __future__ import annotations

import pytest

from chebi_gate import sdf
from tests.chebi_gate_helpers import (
    TEMPLATE_TAGS,
    molblock,
    neutral_record,
    record_text,
    salt_record,
    sdf_bytes,
)

GATE_FIELDS = {
    "GATE_STATUS": "HELD",
    "GATE_CHECKS_FAILED": "CON-02",
    "GATE_RUN": "chebi_gate test",
}


@pytest.fixture
def two_record_sdf() -> bytes:
    return sdf_bytes(salt_record(), neutral_record())


# ------------------------------------------------------------------- parsing


def test_parses_every_record(two_record_sdf):
    parsed = sdf.parse_bytes(two_record_sdf)
    assert len(parsed.records) == 2
    assert [r.index for r in parsed.records] == [1, 2]
    assert [r.name for r in parsed.records] == [
        "ethylamine hydrochloride",
        "ethanol",
    ]


def test_data_fields_keep_template_order(two_record_sdf):
    salt = sdf.parse_bytes(two_record_sdf).records[0]
    assert salt.tag_order == list(TEMPLATE_TAGS)
    assert salt.data["CAS_NO"] == "557-66-4"
    assert salt.data["RELATIONSHIP"] == "ISA36807"


def test_absent_tag_is_absent_not_empty():
    """An omitted tag and an empty tag are different findings, so parse differs."""
    parsed = sdf.parse_bytes(sdf_bytes(salt_record(synonym=None)))
    record = parsed.records[0]
    assert "SYNONYM" not in record.data
    assert record.data["IUPAC_NAME"] == "ethanamine;hydrochloride"


def test_repeated_tag_keeps_the_first_value(caplog):
    text = record_text(molblock("ethanol"), {"NAME": "ethanol", "CAS_NO": "64-17-5"})
    text += "> <NAME>\nsomething else\n\n"
    parsed = sdf.parse_bytes(sdf_bytes(text))
    assert parsed.records[0].data["NAME"] == "ethanol"
    assert "repeats tag" in caplog.text


def test_name_falls_back_to_the_title_line():
    """A record with no NAME still has to be nameable in a finding message."""
    text = record_text(molblock("ethanol"), {"CAS_NO": "64-17-5"})
    record = sdf.parse_bytes(sdf_bytes(text)).records[0]
    assert record.data.get("NAME") is None
    assert record.name == "ethanol"


def test_cas_is_stripped_and_empty_when_absent():
    text = record_text(molblock("ethanol"), {"NAME": "ethanol"})
    assert sdf.parse_bytes(sdf_bytes(text)).records[0].cas == ""
    padded = record_text(
        molblock("ethanol"), {"NAME": "ethanol", "CAS_NO": "  64-17-5  "}
    )
    assert sdf.parse_bytes(sdf_bytes(padded)).records[0].cas == "64-17-5"


def test_mol_block_ends_at_m_end(two_record_sdf):
    record = sdf.parse_bytes(two_record_sdf).records[0]
    assert record.mol_block.endswith("M  END\n")
    assert "> <NAME>" not in record.mol_block


def test_counts_line_is_the_fourth_line(two_record_sdf):
    record = sdf.parse_bytes(two_record_sdf).records[0]
    assert "V2000" in record.counts_line
    assert record.counts_line.startswith(" 12 10")


def test_blank_chunk_between_terminators_is_skipped(caplog):
    raw = sdf_bytes(salt_record()) + b"\n" + sdf.RECORD_TERMINATOR
    parsed = sdf.parse_bytes(raw)
    assert len(parsed.records) == 1
    assert "empty record" in caplog.text


@pytest.mark.parametrize(
    "where",
    ["before the first record", "between two records", "after the last record"],
)
def test_a_blank_chunk_round_trips_where_it_was_written(where):
    """`dumps()` says it reproduces the input bytes exactly, and did not.

    Skipped chunks were a flat tuple concatenated after every record, so a blank
    chunk between two records came back at the end of the file. The one test of
    this put the chunk last, where the reordering cannot be seen -- and the
    290-record round-trip anchor leans on `dumps()`.
    """
    first = salt_record().encode()
    second = salt_record(name="second", cas="64-17-5").encode()
    blank = sdf.RECORD_TERMINATOR
    body = {
        "before the first record": blank + first + blank + second + blank,
        "between two records": first + blank + blank + second + blank,
        "after the last record": first + blank + second + blank + blank,
    }[where]

    parsed = sdf.parse_bytes(body)
    assert len(parsed.records) == 2
    assert parsed.dumps() == body


def test_two_blank_chunks_in_a_row_keep_their_place():
    first = salt_record().encode()
    blank = sdf.RECORD_TERMINATOR
    body = (
        first
        + blank
        + blank
        + blank
        + salt_record(name="s", cas="64-17-5").encode()
        + blank
    )
    parsed = sdf.parse_bytes(body)
    assert len(parsed.records) == 2
    assert parsed.dumps() == body


# --------------------------------------------------------------- byte fidelity


def test_whole_file_round_trip_is_byte_identical(two_record_sdf):
    assert sdf.parse_bytes(two_record_sdf).dumps() == two_record_sdf


def test_each_record_reemits_its_own_bytes(two_record_sdf):
    parsed = sdf.parse_bytes(two_record_sdf)
    joined = b"".join(r.raw + sdf.RECORD_TERMINATOR for r in parsed.records)
    assert joined == two_record_sdf


def test_cleared_file_is_byte_identical_to_input(tmp_path, two_record_sdf):
    """Handoff invariant 19: a cleared record carries no annotation and no drift."""
    source = tmp_path / "in.sdf"
    source.write_bytes(two_record_sdf)
    out = tmp_path / "cleared.sdf"

    parsed = sdf.parse_file(source)
    written = sdf.write_records(parsed.records, out)

    assert written == 2
    assert out.read_bytes() == two_record_sdf
    assert b"GATE_" not in out.read_bytes()


def test_cleared_file_is_written_even_when_empty(tmp_path):
    out = tmp_path / "cleared.sdf"
    assert sdf.write_records([], out) == 0
    assert out.exists()
    assert out.read_bytes() == b""


def test_non_ascii_byte_survives_the_round_trip():
    """The record that most needs handing back for repair must not be corrupted.

    The prototype decoded with ``errors="replace"``, so a non-ASCII byte became
    U+FFFD and could never be re-emitted. The gate refuses such a file (INT-06);
    it must still be able to give the bytes back unchanged.
    """
    # Built by byte substitution rather than by encoding a str, so the fixture
    # builder stays ASCII-only and the non-ASCII byte is unambiguously one byte.
    raw = sdf_bytes(neutral_record(synonym="cafe_ alcohol")).replace(
        b"cafe_", b"caf\xe9"
    )
    assert b"\xe9" in raw

    parsed = sdf.parse_bytes(raw)
    assert parsed.dumps() == raw
    assert parsed.records[0].raw + sdf.RECORD_TERMINATOR == raw


def test_crlf_line_endings_are_preserved():
    raw = sdf_bytes(neutral_record()).replace(b"\n", b"\r\n")
    parsed = sdf.parse_bytes(raw)
    assert parsed.dumps() == raw


def test_unterminated_final_record_is_flagged_and_preserved():
    raw = sdf_bytes(salt_record(), neutral_record(), terminate_last=False)
    parsed = sdf.parse_bytes(raw)

    assert len(parsed.records) == 2
    assert parsed.records[0].terminated is True
    assert parsed.records[1].terminated is False
    assert parsed.dumps() == raw


def test_trailing_whitespace_after_the_last_terminator_is_kept():
    raw = sdf_bytes(neutral_record()) + b"\n\n"
    parsed = sdf.parse_bytes(raw)
    assert len(parsed.records) == 1
    assert parsed.trailer == b"\n\n"
    assert parsed.dumps() == raw


# ------------------------------------------------------------- annotated output


def test_strip_fields_is_an_identity_when_the_fields_are_absent(two_record_sdf):
    for record in sdf.parse_bytes(two_record_sdf).records:
        assert sdf.strip_fields(record.raw, GATE_FIELDS) == record.raw


def test_annotate_then_strip_returns_the_original_bytes(two_record_sdf):
    """The held-back file is an SDF so a fixed record can be fed straight back."""
    for record in sdf.parse_bytes(two_record_sdf).records:
        annotated = record.with_fields(GATE_FIELDS)
        assert annotated != record.raw
        assert sdf.strip_fields(annotated, GATE_FIELDS) == record.raw


def test_annotating_is_idempotent(two_record_sdf):
    """Re-running the gate on its own output replaces annotations, not stacks them."""
    record = sdf.parse_bytes(two_record_sdf).records[0]
    once = record.with_fields(GATE_FIELDS)

    reparsed = sdf.parse_bytes(once + sdf.RECORD_TERMINATOR).records[0]
    twice = reparsed.with_fields(GATE_FIELDS)

    assert twice == once
    assert once.count(b"> <GATE_STATUS>") == 1
    assert twice.count(b"> <GATE_STATUS>") == 1


def test_annotated_record_keeps_its_original_fields(two_record_sdf):
    record = sdf.parse_bytes(two_record_sdf).records[0]
    annotated = record.with_fields(GATE_FIELDS)
    reparsed = sdf.parse_bytes(annotated + sdf.RECORD_TERMINATOR).records[0]

    for tag in TEMPLATE_TAGS:
        assert reparsed.data[tag] == record.data[tag]
    assert reparsed.data["GATE_STATUS"] == "HELD"
    assert reparsed.data["GATE_CHECKS_FAILED"] == "CON-02"


def test_write_annotated_writes_one_record_per_entry(tmp_path, two_record_sdf):
    records = sdf.parse_bytes(two_record_sdf).records
    out = tmp_path / "held.sdf"

    count = sdf.write_annotated([(r, GATE_FIELDS) for r in records], out)

    assert count == 2
    written = sdf.parse_file(out)
    assert len(written.records) == 2
    assert all(r.data["GATE_STATUS"] == "HELD" for r in written.records)


# ------------------------------------------- multi-line values and line endings


MULTILINE_GATE = {
    "GATE_STATUS": "HELD",
    "GATE_REASONS": "CON-02 [high] HOLDS the name says di\n\nEXT-04 [high] registry disagrees",
}


def test_annotating_is_idempotent_when_a_value_contains_a_blank_line(two_record_sdf):
    """Regression. A finding's text is prose and prose has paragraph breaks.

    Removal used a second regex that stopped at the first blank line, while the
    reader treats a blank line as a field terminator only when "> <" follows. So
    the value was half-removed and its orphaned tail grafted onto the preceding
    real field: feeding the held file back in turned RELATIONSHIP into
    "ISA36807\n\nEXT-04 [high] registry disagrees", which the next run then
    reported as an unknown relationship. The gate corrupted a record and faulted
    it for being corrupt.
    """
    record = sdf.parse_bytes(two_record_sdf).records[0]
    once = record.with_fields(MULTILINE_GATE)

    reparsed = sdf.parse_bytes(once + sdf.RECORD_TERMINATOR).records[0]
    twice = reparsed.with_fields(MULTILINE_GATE)

    assert twice == once
    assert reparsed.data["RELATIONSHIP"] == record.data["RELATIONSHIP"]
    assert reparsed.data["GATE_REASONS"] == MULTILINE_GATE["GATE_REASONS"]


def test_append_then_strip_is_an_identity_for_a_multi_line_value(two_record_sdf):
    for record in sdf.parse_bytes(two_record_sdf).records:
        annotated = record.with_fields(MULTILINE_GATE)
        assert sdf.strip_fields(annotated, MULTILINE_GATE) == record.raw


def test_a_field_removed_from_the_middle_leaves_no_gap(two_record_sdf):
    record = sdf.parse_bytes(two_record_sdf).records[0]
    without = record.without_fields(["IUPAC_NAME"])
    reparsed = sdf.parse_bytes(without + sdf.RECORD_TERMINATOR).records[0]

    assert "IUPAC_NAME" not in reparsed.data
    assert reparsed.data["SYNONYM"] == record.data["SYNONYM"]
    assert reparsed.data["CAS_NO"] == record.data["CAS_NO"]
    assert b"\n\n\n" not in without


def test_a_non_ascii_value_can_still_be_annotated():
    """The record the held file most needs to carry must not crash the write.

    with_fields encoded with strict UTF-8 while the rest of the module uses
    ascii+surrogateescape, so a detail echoing a non-ASCII byte raised
    UnicodeEncodeError -- on exactly the record the module exists to hand back.
    """
    raw = sdf_bytes(salt_record(iupac="x")).replace(b"ISA36807", b"\xe9SA36807")
    record = sdf.parse_bytes(raw).records[0]

    annotated = record.with_fields(
        {"GATE_STATUS": "HELD", "GATE_REASONS": f"class {record.data['RELATIONSHIP']}"}
    )
    assert b"\xe9SA36807" in annotated
    reparsed = sdf.parse_bytes(annotated + sdf.RECORD_TERMINATOR).records[0]
    assert reparsed.data["GATE_STATUS"] == "HELD"


def test_crlf_terminators_still_separate_records():
    """Regression, and the reason the old CRLF test was vacuous.

    With a bare-LF terminator, a CRLF file parsed as ONE unterminated record, so
    every per-record check ran against a concatenation of the whole file -- and the
    round-trip assertion still passed, because re-emitting one giant record
    reproduces the input exactly.
    """
    raw = sdf_bytes(salt_record(iupac="x"), neutral_record(cas="64-17-5", iupac="y"))
    crlf = raw.replace(b"\n", b"\r\n")
    parsed = sdf.parse_bytes(crlf)

    assert len(parsed.records) == 2
    assert [r.terminator for r in parsed.records] == [b"$$$$\r\n"] * 2
    assert all(r.terminated for r in parsed.records)
    assert parsed.dumps() == crlf
    assert parsed.malformed == ()


def test_crlf_records_parse_their_data_fields():
    """The half the terminator regression left behind, and it hid the same way.

    Splitting records on either line ending was not enough: the molfile-end marker
    was the plain string "\\nM  END\\n", which a CRLF file spells
    "\\r\\nM  END\\r\\n", so the separator came back empty and the data-field loop
    never ran. Every record parsed with ``data == {}``, and the checks then held a
    sound file for defects it did not have -- INT-03 "NAME missing or empty" on a
    record whose NAME was present, INT-02 a title mismatch against a NAME read as
    absent. INT-06 reports CRLF itself at ``low``, which holds nothing, so there
    was no finding naming the real problem and no route out of it.

    Asserting a *field value* is the point. The test above asserts boundaries and
    round-trip bytes, and both were already correct while every field was missing.
    """
    lf = sdf_bytes(salt_record(iupac="x"), neutral_record(cas="64-17-5", iupac="y"))
    crlf = lf.replace(b"\n", b"\r\n")

    from_lf = [r.data for r in sdf.parse_bytes(lf).records]
    from_crlf = [r.data for r in sdf.parse_bytes(crlf).records]

    assert from_crlf == from_lf
    assert from_crlf[0]["NAME"] == "ethylamine hydrochloride"
    assert from_crlf[0]["CAS_NO"] == "557-66-4"
    assert from_crlf[0]["RELATIONSHIP"] == "ISA36807"
    # No stray carriage return survives into a value, or every CAS and every
    # RELATIONSHIP code would miss its join key by one character.
    assert not any("\r" in v for d in from_crlf for v in d.values())
    # The tag order the annotator writes against, and the spans it removes by, both
    # have to survive too.
    assert [r.tag_order for r in sdf.parse_bytes(crlf).records] == [
        r.tag_order for r in sdf.parse_bytes(lf).records
    ]


def test_the_exact_terminator_bytes_are_preserved_per_record():
    mixed = sdf_bytes(salt_record(iupac="x")).replace(
        b"$$$$\n", b"$$$$\r\n"
    ) + sdf_bytes(neutral_record(cas="64-17-5", iupac="y"))
    parsed = sdf.parse_bytes(mixed)
    assert [r.terminator for r in parsed.records] == [b"$$$$\r\n", b"$$$$\n"]
    assert parsed.dumps() == mixed


# --------------------------------------------------------- malformed detection


@pytest.mark.parametrize(
    "label,build",
    [
        (
            "no final terminator",
            lambda: sdf_bytes(salt_record(iupac="x"), terminate_last=False),
        ),
        ("dollars with no newline", lambda: sdf_bytes(salt_record(iupac="x"))[:-1]),
        (
            "content after the last terminator",
            lambda: sdf_bytes(salt_record(iupac="x")) + b"leftover\n",
        ),
    ],
)
def test_a_malformed_file_says_so_and_still_round_trips(label, build):
    """A truncated input used to clear silently, and its cleared file was five
    bytes longer than the input because the gate supplied the missing terminator.
    """
    data = build()
    parsed = sdf.parse_bytes(data)
    assert parsed.malformed, label
    assert parsed.dumps() == data, label


def test_a_well_formed_file_reports_nothing_malformed(two_record_sdf):
    assert sdf.parse_bytes(two_record_sdf).malformed == ()


def test_write_records_always_terminates_even_an_unterminated_record(tmp_path):
    """So the output is a valid SDF. Callers refuse malformed input separately."""
    data = sdf_bytes(salt_record(iupac="x"), terminate_last=False)
    records = sdf.parse_bytes(data).records
    out = tmp_path / "o.sdf"
    sdf.write_records(records, out)

    written = out.read_bytes()
    assert written.endswith(sdf.RECORD_TERMINATOR)
    assert written.count(b"$$$$") == 1, "no doubled terminator"


def test_a_record_with_no_trailing_blank_line_can_still_be_annotated():
    """Regression, and it was live on the prototype's own held file.

    A data field is separated from the next by a blank line, so appending to a
    record that ends with a single newline runs the new field onto the previous
    value. without_fields returned raw untouched when no requested tag was
    present -- the case on every record's *first* annotation -- so the first
    appended field was swallowed: RELATIONSHIP came back as
    "ISA36807\n> <GATE_STATUS>\nHELD" and GATE_STATUS vanished as a field.

    All seven records of the prototype's HOLD_questionable.sdf are shaped this
    way, because it wrote records as rstrip("\n") + "\n". Only the *first*
    appended field is affected, which is why reading GATE_CHECKS_FAILED looked
    fine.
    """
    single = salt_record(iupac="x").rstrip("\n") + "\n"
    record = sdf.parse_bytes((single + "$$$$\n").encode()).records[0]
    assert not record.raw.endswith(b"\n\n"), "fixture must lack the blank line"

    annotated = record.with_fields(GATE_FIELDS)
    reparsed = sdf.parse_bytes(annotated + sdf.RECORD_TERMINATOR).records[0]

    assert reparsed.data["RELATIONSHIP"] == "ISA36807"
    for tag, value in GATE_FIELDS.items():
        assert reparsed.data[tag] == value, tag
    assert reparsed.with_fields(GATE_FIELDS) == annotated, "still idempotent"


def test_without_fields_normalises_even_when_it_removes_nothing():
    """The normalisation is the mechanism; asserting it directly pins it."""
    single = salt_record(iupac="x").rstrip("\n") + "\n"
    record = sdf.parse_bytes((single + "$$$$\n").encode()).records[0]
    assert record.without_fields(["GATE_STATUS"]).endswith(b"\n\n")
    assert record.without_fields([]).endswith(b"\n\n")


def test_a_crlf_record_gives_a_title_that_matches_its_name():
    """Matching the molfile marker with \r?\n is only half of reading CRLF.

    `text.split("\n")` left the title as "ethylamine hydrochloride\r" while the
    data-field pattern stopped before the carriage return, so INT-02 reported
    "mol title differs from NAME" on every record of a file whose only defect is
    its line endings -- and INT-06 rates those `low`, which holds nothing.
    """
    raw = sdf_bytes(salt_record()).replace(b"\n", b"\r\n")
    record = sdf.parse_bytes(raw).records[0]

    assert record.newline == "\r\n"
    assert "\r" not in record.title
    assert record.title == record.data["NAME"]
    assert "\r" not in record.counts_line


def test_a_crlf_record_survives_the_annotate_and_strip_round_trip():
    """The re-run workflow, on CRLF. The span loop stopped at the "\r".

    So a removed field left its blank line behind, the re-emitted record glued the
    separator onto the previous value, and CAS_NO came back as "557-66-4\r\n" --
    the gate corrupting a record and then faulting it for being corrupt.
    """
    raw = sdf_bytes(salt_record()).replace(b"\n", b"\r\n")
    record = sdf.parse_bytes(raw).records[0]

    annotated = record.with_fields({"GATE_STATUS": "HELD"})
    reparsed = sdf.parse_bytes(annotated + record.terminator).records[0]
    assert reparsed.data["CAS_NO"] == "557-66-4"
    assert reparsed.data["GATE_STATUS"] == "HELD"
    assert b"\r\n" in annotated and b"\n\n" not in annotated.replace(b"\r\n", b"")

    stripped = reparsed.without_fields(["GATE_STATUS"])
    assert stripped == record.raw


def test_an_lf_record_is_unchanged_by_the_line_ending_handling():
    raw = sdf_bytes(salt_record())
    record = sdf.parse_bytes(raw).records[0]
    assert record.newline == "\n"
    annotated = record.with_fields({"GATE_STATUS": "HELD"})
    reparsed = sdf.parse_bytes(annotated + record.terminator).records[0]
    assert reparsed.without_fields(["GATE_STATUS"]) == record.raw


def test_a_dollar_run_at_the_end_of_a_value_does_not_split_the_record():
    """The terminator was unanchored, so "$$$$" ending any line split a record.

    A synonym is the obvious way in, and a value echoed back into a re-fed held
    file is the one that matters -- this module's re-run workflow.
    """
    raw = sdf_bytes(salt_record(synonym="weird name $$$$"))
    parsed = sdf.parse_bytes(raw)

    assert len(parsed.records) == 1
    assert parsed.records[0].data["SYNONYM"] == "weird name $$$$"
    assert parsed.dumps() == raw
