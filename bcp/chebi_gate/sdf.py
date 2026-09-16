"""Byte-faithful SDF parsing and re-emission for the ChEBI submission gate.

The gate splits one input SDF into records cleared for submission and records held
back, and the cleared file has to be byte-identical to its input: ChEBI loads it
with the field template they sent us, and an added field or a dropped blank line
breaks that load. So this module keeps every record's *exact bytes* and decodes
only a throwaway copy for analysis.

That is not how the prototype this was ported from worked, and the difference is
not cosmetic. It decoded the whole file with ``errors="replace"`` before splitting,
then re-emitted records as ``chunk.rstrip("\\n") + "\\n"``. Every record in both
submission files ends with a blank line, so ``rstrip`` removed it: measured on
chebi_bulk_group_B_novel.sdf, 0 of the 13 cleared records came back
byte-identical, each one byte short. Nothing reported this, because the only check
on the cleared file was that it contained no ``GATE_*`` field.

**One notion of where a field ends.** Fields are removed by the byte span the
reader recorded for them, never by a second regex. Two notions is how this broke
once already: the removal pattern was non-greedy up to the first blank line, while
the reader treats a blank line as a field terminator only when ``> <`` follows it.
A ``GATE_REASONS`` value containing an ordinary paragraph break was therefore
half-removed and its orphaned tail grafted onto the preceding real field. Feeding
the held file back in -- the workflow this module exists to support -- turned
``RELATIONSHIP`` into ``ISA36807\\n\\nline two after a blank``, which the next run
reported as an unknown relationship. The gate corrupted a record and then faulted
it for being corrupt.

**Both line endings terminate a record.** A CRLF file whose terminator is
``$$$$\\r\\n`` used to parse as one unterminated record, so every per-record check
ran against a concatenation of the whole file -- and a round-trip test still
passed, because re-emitting one giant record reproduces the input exactly. The
exact terminator bytes are kept per record so re-emission stays faithful.

Non-ASCII input survives a round trip instead of being replaced by U+FFFD, so the
gate can *refuse* a non-ASCII file (INT-06) without having corrupted the one
record it most needs to hand back for repair.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

RECORD_TERMINATOR = b"$$$$\n"

# A record terminator with either line ending.
_TERMINATOR = re.compile(rb"\$\$\$\$\r?\n")

# The molfile ends at this marker; everything after it is the data-field block.
#
# Matched with either line ending, because a CRLF file ends its molfile
# "\r\nM  END\r\n" and a bare-LF marker found none of it: `partition` returned an
# empty separator, the data-field loop never ran, and every record of a CRLF file
# parsed with no fields at all. INT-06 reports CRLF at `low`, which holds nothing,
# so such a file was instead held by INT-02 and INT-03 findings naming defects it
# did not have -- a missing NAME on a record whose NAME was right there. This is
# the record-terminator bug one layer up, and it hid for the same reason: the old
# test asserted boundaries and round-trip bytes, never a parsed field value.
_MOL_END = re.compile(r"\r?\nM  END\r?\n")


def _partition_mol_end(text: str) -> tuple[str, str, str]:
    """``str.partition`` over :data:`_MOL_END`, so both callers split identically.

    Returning the matched text as the separator keeps the span arithmetic in
    :func:`_build` correct: ascii+surrogateescape maps one byte to one character,
    so ``len(head) + len(sep)`` still indexes the original bytes whether the marker
    arrived with a carriage return or without one.
    """
    match = _MOL_END.search(text)
    if match is None:
        return text, "", ""
    return text[: match.start()], match.group(0), text[match.end() :]


# One data field: "> <TAG>", the value on following lines, terminated by a blank
# line that precedes the next "> <" or by end of record.
#
# Kept character-for-character from the prototype validator on purpose. It decides
# which tags exist and what their values are, so every finding about a tag or a
# value is downstream of it; a "tidier" regex silently rewrites the findings table
# and the baseline stops reproducing. In particular the value is non-greedy and
# the terminator is a lookahead, so a value containing a blank line runs on until
# a blank line followed by "> <".
_DATA_FIELD = re.compile(r"> <([^>]+)>\r?\n(.*?)(?=\r?\n\r?\n> <|(?:\r?\n)*\Z)", re.S)

# Text is decoded from bytes with ascii+surrogateescape, which maps each byte to
# exactly one character. That one-to-one mapping is what lets a span measured in
# the decoded text index the original bytes.
_CODEC = "ascii"
_ERRORS = "surrogateescape"


class SdfError(Exception):
    """Raised when an SDF cannot be parsed or a re-emission would not be faithful."""


@dataclass(frozen=True)
class DataField:
    """One data field and where it sits in the record's bytes."""

    tag: str
    value: str
    start: int
    end: int


@dataclass(frozen=True)
class SdfRecord:
    """One SDF record, with its original bytes and a decoded view for analysis.

    ``raw`` is the exact bytes between record terminators, including the trailing
    blank line, and excluding the terminator itself. ``terminator`` holds the exact
    terminator bytes that followed it, or empty when the record was unterminated.
    """

    raw: bytes
    index: int
    text: str
    title: str
    counts_line: str
    data: dict[str, str]
    fields: tuple[DataField, ...] = ()
    terminator: bytes = RECORD_TERMINATOR

    @property
    def tag_order(self) -> list[str]:
        return [f.tag for f in self.fields]

    @property
    def terminated(self) -> bool:
        return bool(self.terminator)

    @property
    def mol_block(self) -> str:
        """The molfile, up to and including ``M  END``."""
        head, sep, _ = _partition_mol_end(self.text)
        return head + sep if sep else self.text

    @property
    def name(self) -> str:
        """The NAME field, falling back to the molfile title line.

        The fallback matters: INT-03 refuses a record whose NAME is missing, but
        the other checks still have to say *which* record they are talking about,
        and a record with no NAME still has a title.
        """
        return self.data.get("NAME", self.title).strip()

    @property
    def cas(self) -> str:
        """The CAS_NO field, stripped. Empty string when absent.

        The gate's join key to waivers, quarantine, the CAS registry table and the
        PubChem cache. NAME is display text and changes when a record is renamed.
        """
        return self.data.get("CAS_NO", "").strip()

    def without_fields(self, tags) -> bytes:
        """This record's bytes with the named data fields removed.

        Removal is by recorded span, so it agrees exactly with what the reader
        parsed however many blank lines a value contains.

        The result always ends with exactly one blank line, whether or not any
        field was actually removed. That normalisation is not cosmetic: a data
        field is separated from the next by a blank line, so appending to a record
        that ends with a single newline runs the new field straight onto the
        previous value. An earlier version returned ``self.raw`` untouched when no
        requested tag was present -- which is the case on every record's *first*
        annotation -- and the first appended field was then swallowed. Measured on
        a record ending ``> <RELATIONSHIP>\nISA36807\n``: RELATIONSHIP came back as
        ``"ISA36807\n> <GATE_STATUS>\nHELD"`` and GATE_STATUS vanished as a field.
        All seven records of the prototype's own held file are shaped that way,
        because it wrote records as ``rstrip("\n") + "\n"``.
        """
        wanted = {str(t) for t in tags}
        spans = [(f.start, f.end) for f in self.fields if f.tag in wanted]
        keep = bytearray()
        cursor = 0
        for start, end in sorted(spans):
            keep += self.raw[cursor:start]
            cursor = max(cursor, end)
        keep += self.raw[cursor:]
        return bytes(keep).rstrip(b"\n") + b"\n\n"

    def with_fields(self, fields: dict[str, str]) -> bytes:
        """This record's bytes with ``fields`` appended, replacing any already there.

        Used for the held file, which is an SDF rather than a report so a fixed
        record can be fed straight back into the gate. Appending is idempotent:
        re-running the gate on its own output replaces the annotations instead of
        stacking a second copy.
        """
        body = self.without_fields(fields)
        for tag, value in fields.items():
            # surrogateescape, not the default strict UTF-8: a value echoed from a
            # non-ASCII record carries lone surrogates, and encoding them strictly
            # raised UnicodeEncodeError on exactly the record this file exists to
            # hand back.
            body += f"> <{tag}>\n{value}\n\n".encode(_CODEC, errors=_ERRORS)
        return body


@dataclass(frozen=True)
class SdfFile:
    """A parsed SDF: its records, its original bytes, and any unterminated tail."""

    records: list[SdfRecord]
    raw: bytes
    trailer: bytes = b""
    path: Path | None = None
    # (index of the record this chunk followed, bytes). Position matters: the
    # chunks used to be a flat tuple concatenated after every record, so a blank
    # chunk between two records came back at the end of the file and `dumps()` was
    # not byte-exact despite saying it was. The only test of it put the blank chunk
    # last, where the reordering is invisible, and the 290-record round-trip anchor
    # leans on `dumps()`.
    skipped_chunks: tuple[tuple[int, bytes], ...] = field(default_factory=tuple)

    def dumps(self) -> bytes:
        """Reproduce the input bytes exactly."""
        after: dict[int, list[bytes]] = {}
        for index, chunk in self.skipped_chunks:
            after.setdefault(index, []).append(chunk)
        out = b"".join(after.get(0, ()))
        for record in self.records:
            out += record.raw + record.terminator
            out += b"".join(after.get(record.index, ()))
        return out + self.trailer

    @property
    def malformed(self) -> tuple[str, ...]:
        """Why this file is not a well-formed SDF, if it is not.

        A caller must refuse to judge records from a malformed file rather than
        work around it: if a record boundary is in doubt, so is every finding
        attributed to a record, and re-emitting invents bytes. A truncated input
        used to clear silently, with a cleared file five bytes longer than its
        input because the missing terminator was supplied for it.
        """
        problems = []
        unterminated = [r.index for r in self.records if not r.terminated]
        if unterminated:
            problems.append(
                f"record {unterminated[0]} is not terminated by $$$$; the file is "
                "truncated or was hand-edited"
            )
        if self.trailer.strip():
            problems.append(
                f"{len(self.trailer)} bytes follow the last terminator: "
                f"{self.trailer[:40]!r}"
            )
        return tuple(problems)


def parse_bytes(data: bytes, *, path: Path | None = None) -> SdfFile:
    """Parse SDF bytes into records, keeping each record's bytes verbatim."""
    records: list[SdfRecord] = []
    skipped: list[bytes] = []
    cursor = 0
    for match in _TERMINATOR.finditer(data):
        chunk = data[cursor : match.start()]
        terminator = data[match.start() : match.end()]
        cursor = match.end()
        if not chunk.strip():
            log.warning("skipping empty record at byte %d", match.start())
            # Tagged with the record it follows, so dumps() can put it back where
            # it was rather than at the end of the file. 0 means "before the first".
            skipped.append((len(records), chunk + terminator))
            continue
        records.append(_build(chunk, len(records) + 1, terminator))

    tail = data[cursor:]
    trailer = b""
    if tail.strip():
        # A record with no terminator. Kept as a record so it can be reported,
        # rather than dropped, which would make cleared+held disagree with the
        # input count.
        records.append(_build(tail, len(records) + 1, b""))
    else:
        trailer = tail

    return SdfFile(
        records=records,
        raw=data,
        trailer=trailer,
        path=path,
        skipped_chunks=tuple(skipped),
    )


def parse_file(path: str | Path) -> SdfFile:
    """Parse an SDF from disk."""
    path = Path(path)
    return parse_bytes(path.read_bytes(), path=path)


def _build(chunk: bytes, index: int, terminator: bytes) -> SdfRecord:
    """Decode one record chunk into an :class:`SdfRecord`.

    Decoding uses ``surrogateescape`` rather than ``replace`` so that a non-ASCII
    byte survives into ``text`` as a lone surrogate and can be encoded back to the
    original byte. ``replace`` is lossy, and a lossy decode in the parser means the
    gate can never re-emit a non-ASCII record unchanged.
    """
    text = chunk.decode(_CODEC, errors=_ERRORS)
    head, sep, rest = _partition_mol_end(text)

    data: dict[str, str] = {}
    fields: list[DataField] = []
    if sep:
        offset = len(head) + len(sep)
        for match in _DATA_FIELD.finditer(rest):
            tag, value = match.group(1), match.group(2)
            if tag in data:
                log.warning("record %d repeats tag %r; keeping the first", index, tag)
                continue
            data[tag] = value
            # Span covers the tag line, the value and the blank line that follows,
            # so removing a field leaves neither a gap nor a doubled blank line.
            start = offset + match.start()
            end = offset + match.end()
            while end < len(text) and text[end] == "\n":
                end += 1
            fields.append(DataField(tag=tag, value=value, start=start, end=end))

    lines = text.split("\n")
    return SdfRecord(
        raw=chunk,
        index=index,
        text=text,
        title=lines[0] if lines else "",
        counts_line=lines[3] if len(lines) > 3 else "",
        data=data,
        fields=tuple(fields),
        terminator=terminator,
    )


def strip_fields(raw: bytes, tags) -> bytes:
    """Remove the named data fields from a record's bytes.

    Parses the record so removal uses the same field boundaries the reader does.
    Prefer :meth:`SdfRecord.without_fields` when a record is already to hand.
    """
    parsed = parse_bytes(raw + RECORD_TERMINATOR)
    if not parsed.records:
        return raw
    return parsed.records[0].without_fields(tags)


def write_records(records: list[SdfRecord], path: str | Path) -> int:
    """Write records to an SDF, each as its original bytes. Returns the count.

    Every record is terminated, including one that arrived unterminated -- an SDF
    requires it. Callers must refuse a malformed input rather than rely on this to
    repair one: see :attr:`SdfFile.malformed`.
    """
    path = Path(path)
    payload = b"".join(r.raw + (r.terminator or RECORD_TERMINATOR) for r in records)
    path.write_bytes(payload)
    return len(records)


def write_annotated(
    annotated: list[tuple[SdfRecord, dict[str, str]]], path: str | Path
) -> int:
    """Write records with per-record fields appended. Returns the count."""
    path = Path(path)
    payload = b"".join(
        record.with_fields(fields) + (record.terminator or RECORD_TERMINATOR)
        for record, fields in annotated
    )
    path.write_bytes(payload)
    return len(annotated)
