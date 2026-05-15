from __future__ import annotations

from pathlib import Path
from typing import Iterator


def load_identifiers(
    path: Path,
    *,
    id_column: str | None = None,
) -> Iterator[str]:
    """Load identifiers from a text, CSV, or TSV file.

    Args:
        path: Input file path.
        id_column: Optional column name for tabular input.

    Yields:
        Stripped non-empty identifier strings.
    """
    raise NotImplementedError
