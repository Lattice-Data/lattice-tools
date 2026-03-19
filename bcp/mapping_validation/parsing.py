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


def _looks_like_s3(s: str) -> bool:
    return s.startswith("s3://")


def _should_skip_line(stripped: str) -> bool:
    """Heuristics for mapping-file header/meta lines.

    In practice, exported mapping files sometimes include:
    - header rows: `Local Path,S3 Path` or `S3 Path,Local Path`
    - metadata rows starting with `@` (e.g. `@NVUS...mapping_processed.csv`)
    - comments/notes starting with `#`
    """

    if stripped.startswith("#") or stripped.startswith("@"):
        return True

    lower = stripped.lower()
    # Common header variants
    if lower.startswith("local path") and "s3 path" in lower:
        return True
    if lower.startswith("s3 path") and "local path" in lower:
        return True

    return False


def _split_mapping_line(line: str) -> Tuple[str, str] | None:
    """
    Split a raw mapping line into (token1, token2).

    The caller will decide which token is S3 vs local (depending on
    provider or heuristics).

    Supports the common formats in existing mapping files:
    - Comma-separated:  s3://...,/ORPROJ1/...
    - Tab-separated:    s3://...\t/ORPROJ1/...
    - Fallback: first “s3://...” token and the rest treated as local path.
    """
    stripped = line.strip()
    if not stripped:
        return None
    if _should_skip_line(stripped):
        return None

    # Preferred: CSV-style with a single comma separating S3 and local
    if "," in stripped:
        col1, col2 = stripped.split(",", 1)
        return col1.strip(), col2.strip()

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
            if (
                stripped[idx] == " "
                and idx + 1 < len(stripped)
                and stripped[idx + 1] == "/"
            ):
                s3 = stripped[:idx]
                local = stripped[idx + 1 :]
                return s3.strip(), local.strip()

    return None


def parse_mapping_file(
    path: str | Path, provider: str | None = None
) -> List[MappingRow]:
    """
    Parse a mapping file into a list of MappingRow objects.

    Each non-empty line is expected to contain exactly two columns:
    two path-like fields (one S3, one local). Lines that cannot be
    parsed into two non-empty fields are skipped, but their line numbers
    are preserved for diagnostics in downstream validators.

    Provider-dependent exports:
    - `psomagen` exports are sometimes `Local Path,S3 Path`
    - `novogene` exports are sometimes `S3 Path,Local Path`

    We use `provider` as a hint, but fall back to heuristics based on
    which token looks like an `s3://...` URL.
    """
    provider = provider.lower() if provider else None
    rows: List[MappingRow] = []
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        for line_num, raw in enumerate(fh, 1):
            parts = _split_mapping_line(raw)
            if parts is None:
                continue

            tok1, tok2 = parts
            if not tok1 or not tok2:
                continue

            # Choose which token is S3. Provider is a hint, not a guarantee.
            if provider == "novogene":
                expected_s3_idx = 0
            elif provider == "psomagen":
                expected_s3_idx = 1
            else:
                expected_s3_idx = 0

            expected_s3 = tok1 if expected_s3_idx == 0 else tok2
            other = tok2 if expected_s3_idx == 0 else tok1

            if _looks_like_s3(expected_s3):
                s3, local = expected_s3, other
            elif _looks_like_s3(other):
                s3, local = other, expected_s3
            else:
                # Neither token looks like an S3 URL. Skip this line.
                continue

            rows.append(MappingRow(s3_path=s3, local_path=local, line_num=line_num))
    return rows


__all__ = [
    "MappingRow",
    "parse_mapping_file",
]
