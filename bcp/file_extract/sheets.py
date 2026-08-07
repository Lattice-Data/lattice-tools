"""Lattice submission-sheet shaping for extracted file metadata.

Turns enriched S3 records into two TSVs whose columns mirror the Lattice
``SequenceFile`` and ``SequenceFileSet`` tabs, so a run's output pastes into the
submission workbook at A2 without reordering. Alias derivation lives here too,
since the sheets reference each other by alias rather than by uuid.

Pure by design: nothing in this module talks to S3.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Iterable, Sequence

from .constants import (
    CRAM_SLOT_COLUMNS,
    FILE_FORMAT_CRAM,
    FILE_FORMAT_FASTQ,
    NO_FILE_AVAILABLE_FALSE,
    READ_SLOT_COLUMNS,
    READ_SLOT_ORDER,
    RUN_CARDINALITY_PAIRED,
    RUN_CARDINALITY_SINGLE,
    SEQUENCE_FILE_COLUMNS,
    SEQUENCE_FILE_SET_COLUMNS,
    STATUS_CURRENT,
)

# Trailing read designator on a deliverable FASTQ. The chunk group is optional so
# both `_R1_001.fastq.gz` and `_R1.fastq.gz` resolve; the set stem is whatever
# precedes the match, which is the alias convention Lattice submissions use.
_FASTQ_SLOT_RE = re.compile(r"_(R[123]|I[12])(?:_(\d+))?\.fastq\.gz$")

_LAB_PATH_PREFIX = "/labs/"
_LAB_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_FORBIDDEN_CELL_CHARS = ("\t", "\n", "\r")

MAX_REPORTED_COLLISIONS = 5


class SheetBuildError(ValueError):
    """Raised when records cannot produce a valid submission sheet."""


@dataclass(frozen=True)
class LabIdentity:
    """Lab as the sheets need it: a linkTo path plus an alias namespace."""

    name: str
    namespace: str

    @property
    def path(self) -> str:
        return f"{_LAB_PATH_PREFIX}{self.name}/"

    @classmethod
    def parse(cls, lab: str, *, namespace: str | None = None) -> "LabIdentity":
        """Accept `heather-marlow` or `/labs/heather-marlow/` interchangeably."""
        name = (lab or "").strip()
        if name.startswith(_LAB_PATH_PREFIX):
            name = name[len(_LAB_PATH_PREFIX) :]
        name = name.strip("/")
        if not _LAB_NAME_RE.match(name):
            raise SheetBuildError(
                f"Invalid --lab {lab!r}: expected a lab name such as "
                "'heather-marlow' or '/labs/heather-marlow/'"
            )
        alias_ns = (namespace or name).strip()
        if not _LAB_NAME_RE.match(alias_ns):
            raise SheetBuildError(
                f"Invalid --alias-namespace {namespace!r}: expected a bare "
                "namespace such as 'heather-marlow'"
            )
        return cls(name=name, namespace=alias_ns)


@dataclass(frozen=True)
class SequenceFileRecord:
    """One deliverable file, ready for a SequenceFile row.

    ``set_stem`` empty means the file carries no parseable read designator, so it
    gets a SequenceFile row but joins no set.
    """

    key: str
    s3_uri: str
    file_format: str
    file_size: int
    set_stem: str = ""
    slot: str = ""
    chunk: str = ""
    crc64nvme_base64: str = ""
    read_count: int | None = None

    @property
    def filename(self) -> str:
        return self.key.rsplit("/", 1)[-1]

    @property
    def directory(self) -> str:
        return self.key.rsplit("/", 1)[0] if "/" in self.key else ""

    @property
    def group_key(self) -> tuple[str, str]:
        """Files sharing one SequenceFileSet row.

        Chunk is deliberately absent: `_R1_001` and `_R1_002` are two pieces of
        one lane's read 1, so they belong to the same set -- where they collide
        on the read1 slot and are rejected, rather than being split into two sets
        that would claim two sequencing runs happened.
        """
        return (self.directory, self.set_stem)

    @property
    def sort_key(self) -> tuple[str, str, str, int, str]:
        slot_index = (
            READ_SLOT_ORDER.index(self.slot) if self.slot in READ_SLOT_ORDER else 0
        )
        return (self.directory, self.set_stem, self.chunk, slot_index, self.filename)

    def file_alias(self, namespace: str) -> str:
        return f"{namespace}:{self.filename}"

    def set_alias(self, namespace: str) -> str:
        return f"{namespace}:{self.set_stem}"


@dataclass(frozen=True)
class FileSetGroup:
    """Files that were sequenced together and share one SequenceFileSet row."""

    set_stem: str
    directory: str
    members: tuple[SequenceFileRecord, ...] = field(default_factory=tuple)

    @property
    def sort_key(self) -> tuple[str, str]:
        return (self.directory, self.set_stem)

    def set_alias(self, namespace: str) -> str:
        return f"{namespace}:{self.set_stem}"


@dataclass(frozen=True)
class SheetOptions:
    """Submitter-supplied values the sheets cannot derive from S3."""

    lab: LabIdentity
    cro_order: str
    is_pilot_order: bool
    sequencing_platform: str
    sequence_file_path: Path
    sequence_file_set_path: Path
    cram_slot: str | None = None


def parse_fastq_slot(filename: str) -> tuple[str, str, str]:
    """Split a FASTQ basename into (set_stem, slot, chunk).

    Returns empty strings when no read designator is present.
    """
    match = _FASTQ_SLOT_RE.search(filename)
    if not match:
        return ("", "", "")
    return (filename[: match.start()], match.group(1), match.group(2) or "")


def cram_set_stem(filename: str) -> str:
    """Set stem for a CRAM: the basename without its extension."""
    if filename.endswith(".cram"):
        return filename[: -len(".cram")]
    return filename


def sample_dir_for(key: str) -> str:
    """Directory naming the CRO group -- the segment ahead of the last /raw/.

    Falls back to the file's own parent when the key has no /raw/ (--no-require-raw).
    """
    head = key.rsplit("/raw/", 1)[0] if "/raw/" in key else key.rsplit("/", 1)[0]
    return head.rsplit("/", 1)[-1]


def build_fastq_record(*, key: str, s3_uri: str, file_size: int) -> SequenceFileRecord:
    set_stem, slot, chunk = parse_fastq_slot(key.rsplit("/", 1)[-1])
    return SequenceFileRecord(
        key=key,
        s3_uri=s3_uri,
        file_format=FILE_FORMAT_FASTQ,
        file_size=file_size,
        set_stem=set_stem,
        slot=slot,
        chunk=chunk,
    )


def build_cram_record(*, key: str, s3_uri: str, file_size: int) -> SequenceFileRecord:
    return SequenceFileRecord(
        key=key,
        s3_uri=s3_uri,
        file_format=FILE_FORMAT_CRAM,
        file_size=file_size,
        set_stem=cram_set_stem(key.rsplit("/", 1)[-1]),
    )


def enrich_record(
    record: SequenceFileRecord,
    *,
    crc: str | None,
    read_count: int | None,
) -> SequenceFileRecord:
    """Attach fetched checksum and read count, leaving failures blank."""
    return replace(
        record,
        crc64nvme_base64=crc or "",
        read_count=read_count,
    )


def validate_plan(
    records: Sequence[SequenceFileRecord],
    *,
    namespace: str,
) -> None:
    """Reject unsubmittable listings before any per-file S3 work is spent.

    Both checks read filenames only, so an order that cannot be submitted fails
    in the first second rather than after a request per file.
    """
    validate_read_slots(records)
    validate_aliases(records, namespace=namespace)


def validate_read_slots(records: Iterable[SequenceFileRecord]) -> None:
    """Reject a set that would need two files in one read slot.

    In practice this means chunked FASTQs -- `_R1_001` beside `_R1_002` for one
    lane. A SequenceFileSet holds a single file per slot, so such a delivery
    cannot be represented, and splitting it into two sets would claim two
    sequencing runs where there was one. Neither Novogene nor Psomagen is
    expected to ship split FASTQs, so this means something upstream changed and
    wants a curator's eyes rather than a guess.
    """
    slots_by_group: dict[tuple[str, str], dict[str, list[SequenceFileRecord]]] = {}
    for record in records:
        if not record.set_stem or not record.slot:
            continue
        by_slot = slots_by_group.setdefault(record.group_key, {})
        by_slot.setdefault(record.slot, []).append(record)

    problems: list[str] = []
    for (_directory, set_stem), by_slot in sorted(slots_by_group.items()):
        for slot, members in sorted(by_slot.items()):
            if len(members) > 1:
                chunks = ", ".join(
                    sorted(record.chunk or "(none)" for record in members)
                )
                problems.append(f"{set_stem} slot {slot} <- chunks {chunks}")

    if problems:
        shown = problems[:MAX_REPORTED_COLLISIONS]
        more = len(problems) - len(shown)
        detail = "\n  ".join(shown)
        suffix = f"\n  ... and {more} more" if more else ""
        raise SheetBuildError(
            "Chunked reads are not supported: a SequenceFileSet holds one file "
            "per read slot, so these must be concatenated into a single file per "
            f"read before submission:\n  {detail}{suffix}"
        )


def validate_aliases(
    records: Iterable[SequenceFileRecord],
    *,
    namespace: str,
) -> None:
    """Reject alias collisions before any per-file S3 work is spent.

    Lattice aliases are unique per lab, so two same-named files under different
    sample directories cannot both be submitted -- and would silently overwrite
    one another's set membership if they were.
    """
    file_alias_keys: dict[str, list[str]] = {}
    set_alias_groups: dict[str, set[tuple[str, str]]] = {}
    for record in records:
        file_alias_keys.setdefault(record.file_alias(namespace), []).append(record.key)
        if record.set_stem:
            set_alias_groups.setdefault(record.set_alias(namespace), set()).add(
                record.group_key
            )

    problems = [
        f"{alias} <- {', '.join(sorted(keys))}"
        for alias, keys in sorted(file_alias_keys.items())
        if len(keys) > 1
    ]
    problems += [
        f"{alias} <- {', '.join(sorted(directory for directory, _ in groups))}"
        for alias, groups in sorted(set_alias_groups.items())
        if len(groups) > 1
    ]
    if problems:
        shown = problems[:MAX_REPORTED_COLLISIONS]
        more = len(problems) - len(shown)
        detail = "\n  ".join(shown)
        suffix = f"\n  ... and {more} more" if more else ""
        raise SheetBuildError(
            "Alias collisions would produce duplicate Lattice aliases; "
            f"resolve these before submitting:\n  {detail}{suffix}"
        )


def group_records(
    records: Sequence[SequenceFileRecord],
) -> tuple[list[FileSetGroup], list[str]]:
    """Bucket records into file sets, returning (groups, warnings).

    Assumes validate_read_slots has already run: a set holding two files for one
    slot raises here too, as a backstop for callers that skip validation.
    """
    warnings: list[str] = []
    buckets: dict[tuple[str, str], list[SequenceFileRecord]] = {}
    for record in sorted(records, key=lambda r: r.sort_key):
        if not record.set_stem:
            warnings.append(
                f"{record.filename}: no read designator parsed; the file gets a "
                "SequenceFile row but joins no SequenceFileSet"
            )
            continue
        buckets.setdefault(record.group_key, []).append(record)

    groups: list[FileSetGroup] = []
    for (directory, set_stem), members in buckets.items():
        validate_read_slots(members)
        groups.append(
            FileSetGroup(
                set_stem=set_stem,
                directory=directory,
                members=tuple(members),
            )
        )

    groups.sort(key=lambda g: g.sort_key)
    warnings.extend(_set_warnings(groups))
    return groups, warnings


def _set_warnings(groups: Sequence[FileSetGroup]) -> list[str]:
    warnings: list[str] = []
    for group in groups:
        by_slot = {member.slot: member for member in group.members}
        if "R2" in by_slot and "R1" not in by_slot:
            warnings.append(
                f"{group.set_stem}: R2 present without R1; read1 will be empty"
            )
        r1, r2 = by_slot.get("R1"), by_slot.get("R2")
        if (
            r1 is not None
            and r2 is not None
            and r1.read_count is not None
            and r2.read_count is not None
            and r1.read_count != r2.read_count
        ):
            warnings.append(
                f"{group.set_stem}: R1 read_count ({r1.read_count}) and R2 "
                f"read_count ({r2.read_count}) differ -- possible truncated file"
            )
    return warnings


def run_cardinality_for(slots: Iterable[str]) -> str:
    """Paired when a second read is present, single otherwise.

    Index and R3 files occupy their own slots and do not change the term; a CRAM
    set has no read designator at all and is single-ended.
    """
    return RUN_CARDINALITY_PAIRED if "R2" in set(slots) else RUN_CARDINALITY_SINGLE


def alias_array(alias: str) -> str:
    """Render a one-element alias array the way the sheet stores it."""
    return json.dumps([alias])


def sequence_file_row(
    record: SequenceFileRecord,
    *,
    lab: LabIdentity,
) -> list[object]:
    """Project a record onto SEQUENCE_FILE_COLUMNS."""
    values = {
        "aliases": alias_array(record.file_alias(lab.namespace)),
        "lab": lab.path,
        "file_format": record.file_format,
        "s3_uri": record.s3_uri,
        "crc64nvme_base64": record.crc64nvme_base64,
        "no_file_available": NO_FILE_AVAILABLE_FALSE,
        "status": STATUS_CURRENT,
        "read_count": record.read_count,
        "file_size": record.file_size,
    }
    return _project(values, SEQUENCE_FILE_COLUMNS)


def sequence_file_set_row(
    group: FileSetGroup,
    *,
    options: SheetOptions,
) -> list[object]:
    """Project a group onto SEQUENCE_FILE_SET_COLUMNS.

    ``library`` is left blank on purpose: it links to a Library object whose alias
    is not derivable from the S3 layout. The diagnostic TSV carries sample_dir and
    set_alias so the column can be filled by lookup.
    """
    lab = options.lab
    values: dict[str, object] = {
        "aliases": alias_array(group.set_alias(lab.namespace)),
        "lab": lab.path,
        "run_cardinality": run_cardinality_for(m.slot for m in group.members),
        "status": STATUS_CURRENT,
        "CRO_order": options.cro_order,
        "is_pilot_order": "TRUE" if options.is_pilot_order else "FALSE",
        "sequencing_platform": options.sequencing_platform,
    }
    for member in group.members:
        values[_slot_column(member, options)] = member.file_alias(lab.namespace)
    return _project(values, SEQUENCE_FILE_SET_COLUMNS)


def _slot_column(record: SequenceFileRecord, options: SheetOptions) -> str:
    if record.slot in READ_SLOT_COLUMNS:
        return READ_SLOT_COLUMNS[record.slot]
    if record.file_format == FILE_FORMAT_CRAM:
        if options.cram_slot not in CRAM_SLOT_COLUMNS:
            raise SheetBuildError(
                f"CRAM {record.filename} needs a slot; pass --cram-slot "
                f"{{{','.join(sorted(CRAM_SLOT_COLUMNS))}}}"
            )
        return CRAM_SLOT_COLUMNS[options.cram_slot]
    raise SheetBuildError(
        f"{record.filename}: read designator {record.slot!r} maps to no set slot"
    )


def _project(values: dict[str, object], columns: Sequence[str]) -> list[object]:
    unknown = set(values) - set(columns)
    if unknown:
        raise SheetBuildError(f"Unknown sheet column(s): {', '.join(sorted(unknown))}")
    return [values.get(column, "") for column in columns]


def _cell(value: object, column: str) -> str:
    """Render one cell, refusing anything that would break the TSV grid."""
    if value is None:
        text = ""
    elif isinstance(value, bool):
        text = "TRUE" if value else "FALSE"
    else:
        text = str(value)
    if any(char in text for char in _FORBIDDEN_CELL_CHARS):
        raise SheetBuildError(
            f"Value for column {column!r} contains a tab or newline and cannot be "
            f"written to a TSV cell: {text!r}"
        )
    return text


def write_sheet(
    path: str | Path,
    columns: Sequence[str],
    rows: Iterable[Sequence[object]],
) -> Path:
    """Write a submission sheet TSV.

    Quoting is disabled: alias arrays hold double quotes, and an unquoted
    ``["lab:file"]`` survives both a spreadsheet import and a raw-text paste,
    whereas the quoted form only survives the former.
    """
    destination = Path(path)
    columns = list(columns)
    with destination.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(
            fh,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            quotechar=None,
            escapechar=None,
            lineterminator="\n",
        )
        writer.writerow(columns)
        for row in rows:
            if len(row) != len(columns):
                raise SheetBuildError(
                    f"Row has {len(row)} values but the sheet has {len(columns)} "
                    "columns"
                )
            writer.writerow(
                [_cell(value, column) for value, column in zip(row, columns)]
            )
    return destination


def write_sheets(
    records: Sequence[SequenceFileRecord],
    groups: Sequence[FileSetGroup],
    *,
    options: SheetOptions,
) -> None:
    """Write both submission sheets for one extraction run."""
    write_sheet(
        options.sequence_file_path,
        SEQUENCE_FILE_COLUMNS,
        (
            sequence_file_row(record, lab=options.lab)
            for record in sorted(records, key=lambda r: r.sort_key)
        ),
    )
    write_sheet(
        options.sequence_file_set_path,
        SEQUENCE_FILE_SET_COLUMNS,
        (sequence_file_set_row(group, options=options) for group in groups),
    )


def validate_cro_order(value: str) -> str:
    """Reject a CRO order that cannot be trusted in an output filename."""
    order = (value or "").strip()
    if not order:
        raise SheetBuildError("--cro-order must not be empty")
    if any(char in order for char in ("/", "\\")) or any(
        char in order for char in _FORBIDDEN_CELL_CHARS
    ):
        raise SheetBuildError(
            f"Invalid --cro-order {value!r}: expected an order identifier such as "
            "'NVUS2024101701-15'"
        )
    return order


def cro_order_mismatch_warning(cro_order: str, prefix: str) -> str | None:
    """Warn when the supplied order does not appear in the S3 prefix."""
    if cro_order.lower() in prefix.lower():
        return None
    return (
        f"--cro-order {cro_order!r} does not appear in the S3 prefix {prefix!r}; "
        "check for a typo or a path from a different order"
    )


def default_sheet_output_names(
    cro_order: str,
    *,
    out_dir: str | Path = ".",
) -> tuple[Path, Path]:
    """Default sheet paths: <order>_SequenceFile.tsv and <order>_SequenceFileSet.tsv."""
    directory = Path(out_dir)
    return (
        directory / f"{cro_order}_SequenceFile.tsv",
        directory / f"{cro_order}_SequenceFileSet.tsv",
    )
