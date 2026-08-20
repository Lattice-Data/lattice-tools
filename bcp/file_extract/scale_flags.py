"""Scale CLI list-flag validation for scale_h5ad."""

from __future__ import annotations

from typing import Sequence

from .s3_utils import parse_s3_uri
from .scale_wells import ScaleExtractError

_FORBIDDEN = frozenset("/\\\t\n\r")


def validate_raw_subdirs(values: Sequence[str]) -> list[str]:
    """Strip raw folder names or ``s3://`` prefixes.

    Tokens may be space-separated (``nargs="+"``) and/or comma-separated
    in a single argument.
    """
    cleaned: list[str] = []
    tokens: list[str] = []
    for value in values:
        tokens.extend(part.strip() for part in (value or "").split(","))
    for item in tokens:
        item = item.rstrip("/")
        if not item:
            raise ScaleExtractError("--raw-subdirs must not contain empty values")
        if item.startswith("s3://"):
            try:
                parse_s3_uri(item + "/")
            except ValueError as exc:
                raise ScaleExtractError(
                    f"Invalid --raw-subdirs {item!r}: {exc}"
                ) from exc
            cleaned.append(item)
            continue
        if any(char in item for char in _FORBIDDEN):
            raise ScaleExtractError(
                f"Invalid --raw-subdirs {item!r}: expected a raw folder "
                "name or s3:// URI"
            )
        cleaned.append(item)
    if not cleaned:
        raise ScaleExtractError("--raw-subdirs must not be empty")
    return cleaned
