"""
Live PubChem + ChEBI tests for structure_check.

Opt-in with `pytest -m "pubchem and chebi"`, or either marker alone. These pin the
real findings the tool was built for, so a change in either upstream API that
would silently flip a verdict shows up here rather than in a curation sheet.
"""

from __future__ import annotations

import pytest

from structure_check.client import (
    MATCH,
    REVIEW_INVESTIGATE,
    SKELETON_DIFFERS,
    STEREO_DIFFERS,
    check_row,
)

pytestmark = [pytest.mark.pubchem, pytest.mark.chebi]


def _skip_if_unresolved(result: dict) -> None:
    """Skip rather than fail when an upstream API is unreachable."""
    unresolved = {"cas_unresolved", "chebi_unresolved", "name_unresolved"}
    if (
        result["id_cas_verdict"] in unresolved
        or result["name_cas_verdict"] in unresolved
    ):
        pytest.skip(f"upstream unresolved: {result}")


def test_live_healthy_row_agrees_three_ways() -> None:
    result = check_row(cas="64-17-5", chebi_id="CHEBI:16236", name="Ethanol")
    _skip_if_unresolved(result)
    assert result["id_cas_verdict"] == MATCH
    assert result["name_cas_verdict"] == MATCH


def test_live_wrong_cas_for_the_named_compound() -> None:
    """
    Alexidine: the row's CAS belongs to an unrelated carboxylic acid.

    The ChEBI ID is correct *for that CAS*, which is why comparing names alone
    could never explain this row.
    """
    result = check_row(
        cas="22573-88-2", chebi_id="CHEBI:27391", name="Alexidine_dihydrochloride"
    )
    _skip_if_unresolved(result)
    assert result["id_cas_verdict"] == MATCH
    assert result["name_cas_verdict"] == SKELETON_DIFFERS
    assert result["review"] == REVIEW_INVESTIGATE


def test_live_stereoisomer_difference_is_not_a_different_molecule() -> None:
    """UCN-01 vs the epimer its CAS resolves to: same skeleton, different stereo."""
    result = check_row(cas="112953-11-4", chebi_id="CHEBI:221840", name="UCN-01")
    _skip_if_unresolved(result)
    assert result["name_cas_verdict"] == STEREO_DIFFERS


def test_live_chebi_id_agrees_with_its_cas_across_a_naming_mismatch() -> None:
    """
    Scriptaid: ChEBI's name is the long systematic one, which reads as a mismatch
    to any string comparison. The structures agree, which is the point.
    """
    result = check_row(cas="287383-59-9", chebi_id="CHEBI:92401", name="Scriptaid")
    _skip_if_unresolved(result)
    assert result["id_cas_verdict"] == MATCH
    assert result["name_cas_verdict"] == MATCH


def test_live_token_fallback_resolves_a_dual_alias_cell() -> None:
    result = check_row(cas="149647-78-9", name="Vorinostat_SAHA")
    _skip_if_unresolved(result)
    assert result["name_cas_verdict"] == MATCH
    assert "tokens:" in result["name_query"]
