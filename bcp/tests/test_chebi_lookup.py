"""Unit tests for chebi_lookup (io, client)."""

from __future__ import annotations

import pytest


@pytest.mark.skip(reason="chebi_lookup implementation pending")
def test_lookup_identifier_placeholder() -> None:
    from chebi_lookup.client import lookup_identifier

    lookup_identifier("CHEBI:15365")
