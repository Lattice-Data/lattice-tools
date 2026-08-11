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


# What a request turned out to be, for callers that need to tell an answer from a
# silence. get_with_retry collapses the last two into None, which is fine when the
# caller only wants the payload and wrong when it wants to know whether PubChem was
# asked at all.
OUTCOME_OK = "ok"
OUTCOME_NOT_FOUND = "not_found"
OUTCOME_UNREACHABLE = "unreachable"


def get_with_retry_status(
    url: str, params: dict | None = None
) -> tuple[requests.Response | None, str]:
    """
    GET with exponential backoff on 429/503, reporting *why* it returned nothing.

    Same request behaviour as get_with_retry — which is now a wrapper over this —
    but a 404 comes back as OUTCOME_NOT_FOUND and exhausted retries as
    OUTCOME_UNREACHABLE. PubChem answering "no such compound" is evidence about the
    compound; PubChem never answering is not, and a caller that cannot tell them
    apart has to treat every outage as a finding.
    """
    delay = RETRY_BACKOFF
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=15)
            if resp.status_code == 200:
                return resp, OUTCOME_OK
            if resp.status_code == 404:
                return None, OUTCOME_NOT_FOUND
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
    return None, OUTCOME_UNREACHABLE


def get_with_retry(url: str, params: dict | None = None) -> requests.Response | None:
    """GET with exponential backoff on 429/503. Returns None on 404 or exhausted retries."""
    resp, _outcome = get_with_retry_status(url, params)
    return resp


def cas_to_cid_status(cas: str) -> tuple[int | None, str]:
    """
    Resolve CAS → PubChem CID, reporting whether PubChem was reachable.

    A 200 whose body is not the JSON this endpoint is documented to return means
    something is answering that is not this API — a maintenance page, a proxy, a
    moved endpoint. That is an outage wearing a success code, so it reports
    UNREACHABLE rather than being mistaken for "no such CAS". An empty CID list is
    a real answer and stays NOT_FOUND.
    """
    resp, outcome = get_with_retry_status(f"{BASE}/compound/name/{cas}/cids/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return None, outcome
    try:
        cids = resp.json().get("IdentifierList", {}).get("CID", [])
        if len(cids) > 1:
            log.debug(
                "CAS %s → %s CIDs, using first (canonical): %s",
                cas,
                len(cids),
                cids[0],
            )
        # Inside the guard: a CID field that is an int or a dict rather than a
        # list is the same kind of malformed payload, and must not escape as a
        # traceback when the docstring above promises UNREACHABLE.
        return (cids[0], OUTCOME_OK) if cids else (None, OUTCOME_NOT_FOUND)
    except (ValueError, KeyError, AttributeError, TypeError, IndexError):
        log.error("Unparseable CID payload for %s — treating as unreachable", cas)
        return None, OUTCOME_UNREACHABLE


def cas_to_cid(cas: str) -> int | None:
    """
    Resolve CAS → PubChem CID via the dedicated name/CID endpoint.

    When multiple CIDs are returned, the first is PubChem's canonical choice.
    """
    cid, _outcome = cas_to_cid_status(cas)
    return cid


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
