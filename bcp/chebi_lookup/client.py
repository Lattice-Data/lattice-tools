"""PubChem PUG REST client: CAS Registry Number → CID → properties + ChEBI xref."""

from __future__ import annotations

import logging
import time
from typing import Any

import requests

# PubChem rate limits: 5 req/s, 400 req/min.
# 3 calls per compound at 0.25s each ≈ 1.3 req/s — well under the limit.
REQUEST_DELAY = 0.25
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

PROPERTIES = ",".join(
    [
        "Title",
        "IUPACName",
        "MolecularFormula",
        "MolecularWeight",
        "IsomericSMILES",
        "CanonicalSMILES",
        "InChI",
        "InChIKey",
        "XLogP",
        "TPSA",
    ]
)

MAX_SYNONYMS = 20

log = logging.getLogger(__name__)

OUTPUT_FIELDS_APPENDED = [
    "pubchem_cid",
    "chebi_id",
    "preferred_name",
    "iupac_name",
    "molecular_formula",
    "molecular_weight",
    "isomeric_smiles",
    "canonical_smiles",
    "inchi",
    "inchikey",
    "xlogp",
    "tpsa",
    "synonyms",
]


def empty_result() -> dict[str, Any]:
    """Return a result dict with all output fields set to empty strings."""
    return {field: "" for field in OUTPUT_FIELDS_APPENDED}


def get_with_retry(url: str, params: dict | None = None) -> requests.Response | None:
    """GET with exponential backoff on 429/503. Returns None on 404 or exhausted retries."""
    delay = RETRY_BACKOFF
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=15)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 404:
                return None
            if resp.status_code in (429, 503):
                log.warning(
                    "Rate limited (%s), waiting %ss [%s/%s]",
                    resp.status_code,
                    delay,
                    attempt,
                    MAX_RETRIES,
                )
                time.sleep(delay)
                delay *= 2
            else:
                log.warning(
                    "HTTP %s [%s/%s]: %s",
                    resp.status_code,
                    attempt,
                    MAX_RETRIES,
                    url,
                )
                time.sleep(delay)
                delay *= 2
        except requests.exceptions.RequestException as exc:
            log.warning(
                "Request error: %s [%s/%s]",
                exc,
                attempt,
                MAX_RETRIES,
            )
            time.sleep(delay)
            delay *= 2
    log.error("Retries exhausted: %s", url)
    return None


def cas_to_cid(cas: str) -> int | None:
    """
    Resolve CAS → PubChem CID via the dedicated name/CID endpoint.

    When multiple CIDs are returned, the first is PubChem's canonical choice.
    """
    resp = get_with_retry(f"{BASE}/compound/name/{cas}/cids/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return None
    try:
        cids = resp.json().get("IdentifierList", {}).get("CID", [])
        if len(cids) > 1:
            log.debug(
                "CAS %s → %s CIDs, using first (canonical): %s",
                cas,
                len(cids),
                cids[0],
            )
        return cids[0] if cids else None
    except (ValueError, KeyError):
        return None


def cid_to_properties(cid: int) -> dict[str, Any]:
    """Fetch scalar properties and ChEBI xref for a PubChem CID."""
    result = empty_result()
    result["pubchem_cid"] = cid

    resp = get_with_retry(f"{BASE}/compound/cid/{cid}/property/{PROPERTIES}/JSON")
    time.sleep(REQUEST_DELAY)
    if resp:
        try:
            props = resp.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            result["preferred_name"] = props.get("Title", "")
            result["iupac_name"] = props.get("IUPACName", "")
            result["molecular_formula"] = props.get("MolecularFormula", "")
            result["molecular_weight"] = props.get("MolecularWeight", "")
            result["isomeric_smiles"] = props.get("IsomericSMILES", "")
            result["canonical_smiles"] = props.get("CanonicalSMILES", "")
            result["inchi"] = props.get("InChI", "")
            result["inchikey"] = props.get("InChIKey", "")
            result["xlogp"] = props.get("XLogP", "")
            result["tpsa"] = props.get("TPSA", "")
        except (ValueError, KeyError, IndexError):
            pass

    resp2 = get_with_retry(f"{BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON")
    time.sleep(REQUEST_DELAY)
    if resp2:
        try:
            reg_ids = (
                resp2.json()
                .get("InformationList", {})
                .get("Information", [{}])[0]
                .get("RegistryID", [])
            )
            for rid in reg_ids:
                if rid.upper().startswith("CHEBI:"):
                    result["chebi_id"] = rid
                    break
        except (ValueError, KeyError, IndexError):
            pass

    return result


def cid_to_synonyms(cid: int) -> str:
    """Fetch synonyms as a pipe-separated string (capped at MAX_SYNONYMS)."""
    resp = get_with_retry(f"{BASE}/compound/cid/{cid}/synonyms/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return ""
    try:
        syns = (
            resp.json()
            .get("InformationList", {})
            .get("Information", [{}])[0]
            .get("Synonym", [])
        )
        return " | ".join(syns[:MAX_SYNONYMS])
    except (ValueError, KeyError, IndexError):
        return ""


def lookup_cas(cas: str) -> dict[str, Any]:
    """Resolve one CAS to PubChem properties, ChEBI xref, and synonyms."""
    cas = cas.strip()
    if not cas:
        return empty_result()

    cid = cas_to_cid(cas)
    if cid is None:
        return empty_result()

    result = cid_to_properties(cid)
    result["synonyms"] = cid_to_synonyms(cid)
    return result
