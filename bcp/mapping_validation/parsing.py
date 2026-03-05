from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass(frozen=True)
class MappingRow:
    """Single row from a mapping file."""

    s3_path: str
    local_path: str
    line_num: int


def _split_mapping_line(line: str) -> Tuple[str, str] | None:
    """
    Split a raw line from a mapping file into (s3_path, local_path).

    Supports the common formats in existing mapping files:
    - Comma-separated:  s3://...,/ORPROJ1/...
    - Tab-separated:    s3://...\t/ORPROJ1/...
    - Fallback: first “s3://...” token and the rest treated as local path.
    """
    stripped = line.strip()
    if not stripped:
        return None

    # Preferred: CSV-style with a single comma separating S3 and local
    if "," in stripped:
        s3, local = stripped.split(",", 1)
        return s3.strip(), local.strip()

    # TSV-style
    if "\t" in stripped:
        parts = stripped.split("\t")
        if len(parts) >= 2:
            return parts[0].strip(), parts[1].strip()

    # Generic “s3:// ... <whitespace> /path” pattern
    if stripped.startswith("s3://"):
        # Find the first space that precedes an absolute local path
        # e.g. "s3://...   /ORPROJ1/..." or "s3://... /mnt/..."
        for idx in range(len(stripped)):
            if stripped[idx] == " " and idx + 1 < len(stripped) and stripped[idx + 1] == "/":
                s3 = stripped[:idx]
                local = stripped[idx + 1 :]
                return s3.strip(), local.strip()

    return None


def parse_mapping_file(path: str | Path) -> List[MappingRow]:
    """
    Parse a mapping file into a list of MappingRow objects.

    Each non-empty line is expected to contain exactly two columns:
    an S3 path and a local filesystem path. Lines that cannot be
    parsed into two columns are skipped, but their line numbers are
    preserved for diagnostics in downstream validators.
    """
    rows: List[MappingRow] = []
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        for line_num, raw in enumerate(fh, 1):
            parts = _split_mapping_line(raw)
            if parts is None:
                continue
            s3, local = parts
            if not s3 or not local:
                continue
            rows.append(MappingRow(s3_path=s3, local_path=local, line_num=line_num))
    return rows


__all__ = [
    "MappingRow",
    "parse_mapping_file",
]

