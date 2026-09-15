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
