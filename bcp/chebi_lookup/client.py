from __future__ import annotations

from typing import Any


def lookup_identifier(identifier: str) -> dict[str, Any]:
    """Query ChEBI (or related services) for a single identifier.

    Args:
        identifier: ChEBI ID or other supported lookup key.

    Returns:
        Parsed result fields for the identifier.
    """
    raise NotImplementedError
