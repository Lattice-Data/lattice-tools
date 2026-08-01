"""Shared helpers for guidesig tests."""

from __future__ import annotations

import csv
import io
from pathlib import Path

from guidesig import FILE_FORMATS

FIXTURES = Path(__file__).parent / "fixtures"


def write_delimited(
    path: Path,
    rows: list[list[str]],
    *,
    file_format: str,
    bom: bool = True,
) -> Path:
    """Write ``rows`` as a CRLF-terminated delimited file.

    Fields are written through ``csv.writer`` so that values containing the
    delimiter are quoted, which is what makes CSV fixtures meaningful.

    Args:
        path: Destination file.
        rows: Header row followed by data rows.
        file_format: Key of ``FILE_FORMATS`` selecting the delimiter.
        bom: Prefix the file with a UTF-8 BOM, as Excel exports do.

    Returns:
        The path that was written.
    """
    buffer = io.StringIO()
    writer = csv.writer(
        buffer,
        delimiter=FILE_FORMATS[file_format],
        lineterminator="\r\n",
    )
    writer.writerows(rows)
    body = buffer.getvalue()
    path.write_bytes((("\ufeff" + body) if bom else body).encode("utf-8"))
    return path
