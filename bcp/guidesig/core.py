"""Core library for protospacer set signatures."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path

_ACGT = frozenset("ACGT")
DEFAULT_COLUMN = "guide_protospacer"


class GuideSigError(ValueError):
    """Raised when a guide-template TSV cannot yield a valid signature."""


def protospacer_set(path: str | Path, column: str = DEFAULT_COLUMN) -> set[str]:
    """Return the set of uppercase ACGT protospacer strings from a TSV."""
    file_path = Path(path)
    try:
        with file_path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise GuideSigError(
                    f"{file_path}: file is empty or not TSV-parseable"
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

                if len(row) <= col_idx:
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
    except OSError as exc:
        raise GuideSigError(f"{file_path}: file unreadable: {exc}") from exc
    except csv.Error as exc:
        raise GuideSigError(f"{file_path}: not TSV-parseable: {exc}") from exc

    if not sequences:
        raise GuideSigError(f"{file_path}: zero usable sequences after filtering")

    return sequences


def signature(path: str | Path, column: str = DEFAULT_COLUMN) -> str:
    """Return the gsig1 set signature for the protospacers in ``path``."""
    sequences = protospacer_set(path, column=column)
    sorted_sequences = sorted(sequences)
    payload = "\n".join(sorted_sequences)
    digest = sha256(payload.encode("utf-8")).hexdigest()[:32]
    return f"gsig1:set:n={len(sequences)}:{digest}"
