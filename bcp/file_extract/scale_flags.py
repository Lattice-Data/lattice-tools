"""Shared Scale CLI list-flag validation for scale_h5ad and later scale_cram."""

from __future__ import annotations

from typing import Sequence

from .scale_wells import ScaleExtractError

_FORBIDDEN = frozenset("/\\\t\n\r")


def validate_id_list(values: Sequence[str], *, flag: str) -> list[str]:
    """Strip and reject empty or path-like identifiers."""
    cleaned: list[str] = []
    for value in values:
        item = (value or "").strip()
        if not item:
            raise ScaleExtractError(f"{flag} must not contain empty values")
        if any(char in item for char in _FORBIDDEN):
            raise ScaleExtractError(
                f"Invalid {flag} {value!r}: expected an identifier without path characters"
            )
        cleaned.append(item)
    if not cleaned:
        raise ScaleExtractError(f"{flag} must not be empty")
    return cleaned
