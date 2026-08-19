"""Live PubChem API tests for chebi_lookup (opt-in with pytest -m pubchem)."""

from __future__ import annotations

import pytest
import requests

from chebi_lookup.client import cas_to_cid, lookup_cas

pytestmark = pytest.mark.pubchem

REFERENCE_CAS = "64-17-5"


def _skip_on_network_error(exc: BaseException) -> None:
    if isinstance(exc, requests.exceptions.RequestException):
        pytest.skip(f"PubChem unreachable: {exc}")


def test_cas_to_cid_live() -> None:
    try:
        cid = cas_to_cid(REFERENCE_CAS)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert isinstance(cid, int)
    assert cid > 0


def test_lookup_cas_live_has_cid_and_chebi() -> None:
    try:
        result = lookup_cas(REFERENCE_CAS)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["pubchem_cid"]
    assert str(result["chebi_id"]).upper().startswith("CHEBI:")


def test_lookup_cas_live_returns_both_smiles_fields() -> None:
    """The property rename that emptied these returned no error to notice.

    Asking by `IsomericSMILES`/`CanonicalSMILES` while PubChem answered by
    `SMILES`/`ConnectivitySMILES` gave "" on every lookup with no exception, no
    non-200 and no log line. The mocked suite cannot catch that recurring, because
    its fixtures are recordings and answer with whatever names were current when they
    were made. Only a live call proves PubChem still accepts the names in PROPERTIES
    in the request URL *and* still answers with them.
    """
    try:
        result = lookup_cas(REFERENCE_CAS)
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["isomeric_smiles"]
    assert result["canonical_smiles"]


def test_lookup_cas_live_unknown_cas() -> None:
    try:
        result = lookup_cas("00-00-0")
    except requests.exceptions.RequestException as exc:
        _skip_on_network_error(exc)
        raise
    assert result["pubchem_cid"] == ""
    assert result["chebi_id"] == ""
