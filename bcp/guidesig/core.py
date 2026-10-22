"""Core library for protospacer set signatures."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path

_ACGT = frozenset("ACGT")
DEFAULT_COLUMN = "guide_protospacer"
ID_COLUMN = "guide_id"

# Role/target consistency: a control guide has no target gene, so its
# guide_target_gene_id must be one of NA_VALUES. The reverse is not enforced:
# a non-control guide may carry NA (e.g. "Targeting a Region" guides). The role
# is compared case-insensitively; NA_VALUES are exact spellings.
ROLE_COLUMN = "guide_role"
TARGET_COLUMN = "guide_target_gene_id"
CONTROL_ROLE = "control"
NA_VALUES = ("NA", "na", "N/A", "n/a")
_NA_EXPECTED = "one of " + ", ".join(repr(v) for v in NA_VALUES)
# Held to the same control rule as guide_target_gene_id, but only when present.
TARGET_NAME_COLUMN = "guide_target_gene_name"

# The all-trailing-G / all-leading-C check needs this many unique protospacers:
# in a smaller set every guide can share an end by chance.
END_CHECK_MIN_PROTOSPACERS = 10

# Columns whose absence is warned about but never fails a run.
RECOMMENDED_COLUMNS = (TARGET_NAME_COLUMN,)

# Single source of truth for the accepted input formats: the parser reads the
# delimiter from here and the CLI derives its --format choices from the keys.
FILE_FORMATS: dict[str, str] = {"tsv": "\t", "csv": ","}


class GuideSigError(ValueError):
    """Raised when a guide-template file cannot yield a valid signature."""


def _delimiter_for(file_format: str) -> str:
    """Return the field delimiter for ``file_format``."""
    # Neither except arm is reachable from the CLI, since argparse restricts
    # --format to the keys of FILE_FORMATS. They serve library callers: KeyError
    # for a name that is merely wrong (including None), TypeError for a value
    # that is not hashable at all. Both stay GuideSigError, like every other
    # failure this module reports.
    try:
        return FILE_FORMATS[file_format]
    except (KeyError, TypeError) as exc:
        accepted = ", ".join(sorted(FILE_FORMATS))
        raise GuideSigError(
            f"unknown file format {file_format!r} (accepted: {accepted})"
        ) from exc


def _wrong_format_hint(header: list[str], file_format: str) -> str:
    """Return a hint naming a better-fitting format, or an empty string.

    A header that collapsed to a single field while still holding another
    format's delimiter is the fingerprint of a ``file_format`` that does not
    match the file. Saying so turns the resulting "column absent" report into
    something the caller can act on, which matters most for ``--compare``: it
    reads both files under one ``--format``, so two files in different formats
    always fail here.
    """
    if len(header) != 1:
        return ""
    fits = sorted(
        name
        for name, delimiter in FILE_FORMATS.items()
        if name != file_format and delimiter in header[0]
    )
    if not fits:
        return ""
    return (
        f"; the whole header parsed as one {file_format} field, so this file "
        f"looks like {' or '.join(fits)} rather than {file_format} "
        "(check the format argument)"
    )


def header_warnings(path: str | Path, *, file_format: str) -> list[str]:
    """Return warnings for recommended columns absent from the header.

    These never fail a run. A file that cannot be read yields no warnings here;
    the signature read reports that failure itself.
    """
    file_path = Path(path)
    try:
        with file_path.open(encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter=_delimiter_for(file_format)))
    except (GuideSigError, OSError, UnicodeDecodeError, csv.Error, StopIteration):
        return []
    return [
        f"{file_path}: column {name!r} absent from header"
        for name in RECOMMENDED_COLUMNS
        if name not in header
    ]


def protospacer_set(
    path: str | Path,
    *,
    file_format: str,
    column: str = DEFAULT_COLUMN,
) -> set[str]:
    """Return the set of uppercase ACGT protospacer strings from a guide file.

    Args:
        path: Guide-template file to read.
        file_format: Input format, one of the keys of ``FILE_FORMATS``. Required;
            the format is never inferred from the file name or contents.
        column: Header name of the protospacer column.

    Raises:
        GuideSigError: If the format is unknown, the file cannot yield a set,
            a control row's ``guide_target_gene_id`` (or, when that column is
            present, ``guide_target_gene_name``) is not ``NA``, or
            every protospacer ends with G or every one starts with C (checked
            only for sets of at least 10 unique protospacers).
    """
    sequences, _ids = _read_guide_file(
        path, file_format=file_format, column=column, id_column=None
    )
    return sequences


def protospacers_and_guide_ids(
    path: str | Path,
    *,
    file_format: str,
    column: str = DEFAULT_COLUMN,
    id_column: str = ID_COLUMN,
) -> tuple[set[str], set[str]]:
    """Return the protospacer set and the guide_id set from one read of a file.

    Every check ``protospacer_set`` applies runs here too. Guide IDs are taken
    from the same rows the protospacers are (comment and blank rows skipped);
    an empty ``guide_id`` cell is not collected.

    Raises:
        GuideSigError: For any ``protospacer_set`` failure, if ``id_column`` is
            absent from the header, or if no row carries a guide ID.
    """
    sequences, ids = _read_guide_file(
        path, file_format=file_format, column=column, id_column=id_column
    )
    if not ids:
        raise GuideSigError(f"{Path(path)}: zero {id_column} values found")
    return sequences, ids


def _read_guide_file(
    path: str | Path,
    *,
    file_format: str,
    column: str,
    id_column: str | None,
) -> tuple[set[str], set[str]]:
    """Read and check a guide file; return (protospacers, guide IDs).

    Guide IDs are collected only when ``id_column`` is given, so a plain
    signature run does not require that column.
    """
    delimiter = _delimiter_for(file_format)
    label = file_format.upper()
    file_path = Path(path)
    try:
        with file_path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise GuideSigError(
                    f"{file_path}: file is empty or not {label}-parseable"
                ) from exc

            try:
                col_idx = header.index(column)
            except ValueError as exc:
                raise GuideSigError(
                    f"{file_path}: protospacer column {column!r} absent from "
                    f"header{_wrong_format_hint(header, file_format)}"
                ) from exc

            for required in (ROLE_COLUMN, TARGET_COLUMN):
                if required not in header:
                    raise GuideSigError(
                        f"{file_path}: required column {required!r} absent from header"
                    )
            role_idx = header.index(ROLE_COLUMN)
            target_idx = header.index(TARGET_COLUMN)
            # Optional column: checked only when present (its absence is a
            # header warning, see header_warnings).
            name_idx = (
                header.index(TARGET_NAME_COLUMN)
                if TARGET_NAME_COLUMN in header
                else None
            )

            id_idx: int | None = None
            if id_column is not None:
                if id_column not in header:
                    raise GuideSigError(
                        f"{file_path}: required column {id_column!r} absent from header"
                    )
                id_idx = header.index(id_column)

            sequences: set[str] = set()
            ids: set[str] = set()
            # Reported row numbers are physical file lines, derived from
            # reader.line_num rather than a record counter: a quoted cell holding
            # an embedded newline spans several lines, and a record counter would
            # then under-report every later row against what the user sees when
            # they open the file. line_num is the last line of the record just
            # read, so the next record begins on the line after it.
            record_line = reader.line_num + 1
            for row in reader:
                physical_row = record_line
                record_line = reader.line_num + 1

                # Fully blank lines (common trailing Excel/Sheets artifact) are
                # skipped; a row with cells but an empty first column is retained.
                if not row or all(cell == "" for cell in row):
                    continue

                # The comment marker is read from the first column only, whatever
                # `column` is set to, because these templates carry their
                # instruction rows in the leading `property` column and a `#`
                # elsewhere is data. A template that moved `property` off the
                # front would stop having its comment rows skipped, and they
                # would then be dropped as empty or rejected as non-ACGT at a
                # named line: a loud failure, never a false match.
                first_col = row[0]
                if first_col.lstrip().startswith("#"):
                    continue

                if len(row) < len(header):
                    raise GuideSigError(
                        f"{file_path}: row {physical_row}: ragged row shorter "
                        "than header"
                    )

                # A trailing delimiter leaves one surplus empty field and is a
                # routine export artifact, so surplus fields are tolerated while
                # they hold nothing. Values in columns the header never declared
                # mean a truncated header or a shifted export, which makes the
                # protospacer read from this row untrustworthy.
                surplus = row[len(header) :]
                if any(cell.strip() for cell in surplus):
                    raise GuideSigError(
                        f"{file_path}: row {physical_row}: ragged row with data "
                        "in fields beyond the header"
                    )

                role = row[role_idx].strip()
                target = row[target_idx].strip()
                is_control = role.lower() == CONTROL_ROLE
                if is_control and target not in NA_VALUES:
                    raise GuideSigError(
                        f"{file_path}: row {physical_row}: {ROLE_COLUMN} is "
                        f"{role!r} but {TARGET_COLUMN} is {target!r} "
                        f"(expected {_NA_EXPECTED})"
                    )
                if is_control and name_idx is not None:
                    target_name = row[name_idx].strip()
                    if target_name not in NA_VALUES:
                        raise GuideSigError(
                            f"{file_path}: row {physical_row}: {ROLE_COLUMN} is "
                            f"{role!r} but {TARGET_NAME_COLUMN} is "
                            f"{target_name!r} (expected {_NA_EXPECTED})"
                        )

                if id_idx is not None:
                    guide_id = row[id_idx].strip()
                    if guide_id:
                        ids.add(guide_id)

                value = row[col_idx].strip().upper()
                if not value:
                    continue

                if not set(value) <= _ACGT:
                    raise GuideSigError(
                        f"{file_path}: row {physical_row}: non-ACGT character "
                        f"in protospacer value {value!r}"
                    )

                sequences.add(value)
    except UnicodeDecodeError as exc:
        raise GuideSigError(
            f"{file_path}: not valid UTF-8 (expected utf-8/utf-8-sig): {exc}"
        ) from exc
    except OSError as exc:
        raise GuideSigError(f"{file_path}: file unreadable: {exc}") from exc
    except csv.Error as exc:
        raise GuideSigError(f"{file_path}: not {label}-parseable: {exc}") from exc

    if not sequences:
        raise GuideSigError(f"{file_path}: zero usable sequences after filtering")

    # Set-wide end checks: a single guide ending in G or starting with C is
    # normal, but every guide doing so suggests a systematic artifact in how the
    # sequences were exported. Each condition is checked on its own, and only
    # for sets of at least END_CHECK_MIN_PROTOSPACERS unique protospacers.
    count = len(sequences)
    if count < END_CHECK_MIN_PROTOSPACERS:
        return sequences, ids
    if all(seq.endswith("G") for seq in sequences):
        raise GuideSigError(
            f"{file_path}: Probable trailing G issue, all {count} protospacers "
            "end with 'G'"
        )
    if all(seq.startswith("C") for seq in sequences):
        raise GuideSigError(
            f"{file_path}: Probable leading C issue, all {count} protospacers "
            "start with 'C'"
        )

    return sequences, ids


def signature(
    path: str | Path,
    *,
    file_format: str,
    column: str = DEFAULT_COLUMN,
) -> str:
    """Return the gsig1 set signature for the protospacers in ``path``.

    The signature covers only the protospacer set, so the same sequences stored
    as TSV and as CSV produce an identical signature.
    """
    sequences = protospacer_set(path, file_format=file_format, column=column)
    return signature_of(sequences)


def signature_of(sequences: set[str]) -> str:
    """Return the gsig1 set signature for an already-checked protospacer set."""
    sorted_sequences = sorted(sequences)
    payload = "\n".join(sorted_sequences)
    digest = sha256(payload.encode("utf-8")).hexdigest()[:32]
    return f"gsig1:set:n={len(sequences)}:{digest}"
