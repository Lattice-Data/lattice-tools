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
byte-identical, each one byte short. Nothing reported this, because the only
check on the cleared file was that it contained no ``GATE_*`` field.

Two consequences of working in bytes, both deliberate:

- non-ASCII input survives a round trip instead of being replaced by U+FFFD, so
  the gate can *refuse* a non-ASCII file (INT-06) without having corrupted it
- the trailing bytes after the final ``$$$$`` are kept verbatim in
  :attr:`SdfFile.trailer`, so :meth:`SdfFile.dumps` reproduces the input exactly
  even when the file is malformed
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

RECORD_TERMINATOR = b"$$$$\n"

# The molfile ends at this marker; everything after it is the data-field block.
MOL_END = "\nM  END\n"

# One data field: "> <TAG>", the value on following lines, terminated by a blank
# line that precedes the next "> <" or by end of record.
#
# Kept character-for-character from the prototype validator on purpose. It decides
# which tags exist and what their values are, so every finding about a tag or a
# value is downstream of it; a "tidier" regex silently rewrites the findings table
# and the baseline stops reproducing. In particular the value is non-greedy and
# the terminator is a lookahead, so a value containing a blank line is truncated
# at that blank line rather than swallowing the next field.
_DATA_FIELD = re.compile(r"> <([^>]+)>\n(.*?)(?=\n\n> <|\n*\Z)", re.S)

# A data field block to remove, by tag name. Matches the "> <TAG>\nvalue" run and
# the blank line that separates it from whatever follows, so removing a field from
# the middle of a record does not leave a double blank line behind.
_FIELD_BLOCK = r"> <{tag}>\n.*?(?:\n\n|\n*\Z)"


class SdfError(Exception):
    """Raised when an SDF cannot be parsed or a re-emission would not be faithful."""


@dataclass(frozen=True)
class SdfRecord:
    """One SDF record, with its original bytes and a decoded view for analysis.

    ``raw`` is the exact bytes between record terminators, including the trailing
    blank line, and excluding the ``$$$$\\n`` terminator itself. Re-emitting
    ``raw + RECORD_TERMINATOR`` reproduces the input for any record that was
    terminated in the input.
    """

    raw: bytes
    index: int
    text: str
    title: str
    counts_line: str
    data: dict[str, str]
    tag_order: list[str] = field(default_factory=list)
    terminated: bool = True

    @property
    def mol_block(self) -> str:
        """The molfile, up to and including ``M  END``."""
        head, sep, _ = self.text.partition(MOL_END)
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

        This is the gate's join key -- to waivers, to quarantine, to the CAS
        registry table and to the PubChem cache. NAME is display text and changes
        when a record is renamed; CAS does not.
        """
        return self.data.get("CAS_NO", "").strip()

    def with_fields(self, fields: dict[str, str]) -> bytes:
        """This record's bytes with ``fields`` appended, replacing any already there.

        Used for the held-back file, which is an SDF rather than a report so that a
        fixed record can be fed straight back into the gate. Appending is
        idempotent: re-running the gate on its own output replaces the annotations
        instead of stacking a second copy.
        """
        body = strip_fields(self.raw, fields)
        for tag, value in fields.items():
            block = f"> <{tag}>\n{value}\n\n".encode()
            body += block
        return body


@dataclass(frozen=True)
class SdfFile:
    """A parsed SDF: its records, its original bytes, and any unterminated tail."""

    records: list[SdfRecord]
    raw: bytes
    trailer: bytes = b""
    path: Path | None = None

    def dumps(self) -> bytes:
        """Reproduce the input bytes exactly."""
        out = b"".join(
            r.raw + (RECORD_TERMINATOR if r.terminated else b"") for r in self.records
        )
        return out + self.trailer


def strip_fields(raw: bytes, tags: object) -> bytes:
    """Remove the named data fields from a record's bytes.

    ``tags`` is any iterable of tag names; a mapping is accepted so a caller can
    pass the same dict it is about to append.

    The result always ends with exactly one blank line, which is the shape every
    record in the submission files has. Enforcing it here rather than at each call
    site is what keeps :meth:`SdfRecord.with_fields` idempotent.
    """
    text = raw.decode("ascii", errors="surrogateescape")
    for tag in tags:
        text = re.sub(
            _FIELD_BLOCK.format(tag=re.escape(str(tag))), "", text, flags=re.S
        )
    return text.rstrip("\n").encode("ascii", errors="surrogateescape") + b"\n\n"


def parse_bytes(data: bytes, *, path: Path | None = None) -> SdfFile:
    """Parse SDF bytes into records, keeping each record's bytes verbatim."""
    records: list[SdfRecord] = []
    chunks = data.split(RECORD_TERMINATOR)

    # split() on a terminated file leaves one empty tail chunk; on an unterminated
    # file the tail holds a real record. Distinguishing the two is what lets the
    # gate refuse a truncated file instead of quietly inventing a terminator.
    tail = chunks.pop()
    trailer = b""
    if tail.strip():
        records_tail: bytes | None = tail
    else:
        records_tail = None
        trailer = tail

    for position, chunk in enumerate(chunks, start=1):
        if not chunk.strip():
            log.warning("skipping empty record at position %d", position)
            continue
        records.append(_build(chunk, len(records) + 1, terminated=True))

    if records_tail is not None:
        records.append(_build(records_tail, len(records) + 1, terminated=False))

    return SdfFile(records=records, raw=data, trailer=trailer, path=path)


def parse_file(path: str | Path) -> SdfFile:
    """Parse an SDF from disk."""
    path = Path(path)
    return parse_bytes(path.read_bytes(), path=path)


def _build(chunk: bytes, index: int, *, terminated: bool) -> SdfRecord:
    """Decode one record chunk into an :class:`SdfRecord`.

    Decoding uses ``surrogateescape`` rather than ``replace`` so that a non-ASCII
    byte survives into ``text`` as a lone surrogate and can be encoded back to the
    original byte. ``replace`` is lossy, and a lossy decode in the parser means the
    gate can never re-emit a non-ASCII record unchanged -- which is precisely the
    record it most needs to hand back for repair.
    """
    text = chunk.decode("ascii", errors="surrogateescape")
    _, sep, rest = text.partition(MOL_END)

    data: dict[str, str] = {}
    tag_order: list[str] = []
    if sep:
        for match in _DATA_FIELD.finditer(rest):
            tag, value = match.group(1), match.group(2)
            if tag in data:
                log.warning("record %d repeats tag %r; keeping the first", index, tag)
                continue
            data[tag] = value
            tag_order.append(tag)

    lines = text.split("\n")
    return SdfRecord(
        raw=chunk,
        index=index,
        text=text,
        title=lines[0] if lines else "",
        counts_line=lines[3] if len(lines) > 3 else "",
        data=data,
        tag_order=tag_order,
        terminated=terminated,
    )


def write_records(records: list[SdfRecord], path: str | Path) -> int:
    """Write records to an SDF, each as its original bytes. Returns the count.

    The file is always created, even with no records, so a downstream step can rely
    on it existing rather than branching on whether the gate found anything.
    """
    path = Path(path)
    payload = b"".join(r.raw + RECORD_TERMINATOR for r in records)
    path.write_bytes(payload)
    return len(records)


def write_annotated(
    annotated: list[tuple[SdfRecord, dict[str, str]]], path: str | Path
) -> int:
    """Write records with per-record fields appended. Returns the count."""
    path = Path(path)
    payload = b"".join(
        record.with_fields(fields) + RECORD_TERMINATOR for record, fields in annotated
    )
    path.write_bytes(payload)
    return len(annotated)
