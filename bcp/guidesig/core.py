"""Core library for protospacer set signatures."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path

_ACGT = frozenset("ACGT")
DEFAULT_COLUMN = "guide_protospacer"

# Single source of truth for the accepted input formats: the parser reads the
# delimiter from here and the CLI derives its --format choices from the keys.
FILE_FORMATS: dict[str, str] = {"tsv": "\t", "csv": ","}


class GuideSigError(ValueError):
    """Raised when a guide-template file cannot yield a valid signature."""


def _delimiter_for(file_format: str) -> str:
    """Return the field delimiter for ``file_format``."""
    try:
        return FILE_FORMATS[file_format]
    except (KeyError, TypeError) as exc:
        accepted = ", ".join(sorted(FILE_FORMATS))
        raise GuideSigError(
            f"unknown file format {file_format!r} (accepted: {accepted})"
        ) from exc


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
        GuideSigError: If the format is unknown, or the file cannot yield a set.
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
                    f"{file_path}: protospacer column {column!r} absent from header"
                ) from exc

            sequences: set[str] = set()
            for physical_row, row in enumerate(reader, start=2):
                # Fully blank lines (common trailing Excel/Sheets artifact) are
                # skipped; a row with cells but an empty first column is retained.
                if not row or all(cell == "" for cell in row):
                    continue

                first_col = row[0]
                if first_col.lstrip().startswith("#"):
                    continue

                if len(row) < len(header):
                    raise GuideSigError(
                        f"{file_path}: row {physical_row}: ragged row shorter "
                        "than header"
                    )

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

    return sequences


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
    sorted_sequences = sorted(sequences)
    payload = "\n".join(sorted_sequences)
    digest = sha256(payload.encode("utf-8")).hexdigest()[:32]
    return f"gsig1:set:n={len(sequences)}:{digest}"
