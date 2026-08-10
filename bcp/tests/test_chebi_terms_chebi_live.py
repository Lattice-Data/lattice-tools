"""Live ChEBI API tests for chebi_terms (opt-in with pytest -m chebi)."""

from __future__ import annotations

import pytest
import requests

from chebi_terms.client import (
    CAS_CONFIRMED,
    STATUS_NOT_FOUND,
    STATUS_OK,
    STATUS_SECONDARY,
    describe_chebi_id,
    verify_chebi_id,
)

pytestmark = pytest.mark.chebi

REFERENCE_ID = "CHEBI:16236"  # ethanol
REFERENCE_CAS = "64-17-5"
SECONDARY_ID = "CHEBI:18484"  # resolves to CHEBI:90, (−)-epicatechin


def _skip_on_network_error(exc: BaseException) -> None:
    if isinstance(exc, requests.exceptions.RequestException):
        pytest.skip(f"ChEBI unreachable: {exc}")


def test_describe_live_returns_name_and_synonyms() -> None:
    try:
        result = describe_chebi_id(REFERENCE_ID)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["id_status"] == STATUS_OK
    assert result["chebi_accession"] == REFERENCE_ID
    assert result["chebi_name"] == "ethanol"
    assert result["chebi_synonyms"]
    # Markup must never reach callers.
    assert "<" not in result["chebi_synonyms"]


def test_verify_live_confirms_cas_from_chebi_xrefs() -> None:
    try:
        result = verify_chebi_id(REFERENCE_ID, expected_cas=REFERENCE_CAS)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["cas_verdict"] == CAS_CONFIRMED


def test_verify_live_flags_secondary_id_and_reports_primary() -> None:
    try:
        result = verify_chebi_id(SECONDARY_ID)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["id_status"] == STATUS_SECONDARY
    assert result["chebi_accession"] != SECONDARY_ID
    assert result["chebi_accession"].startswith("CHEBI:")


def test_verify_live_unknown_id() -> None:
    try:
        result = verify_chebi_id("CHEBI:99999999")
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["id_status"] == STATUS_NOT_FOUND
