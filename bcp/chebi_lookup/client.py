"""PubChem PUG REST client: CAS Registry Number → CID → properties + ChEBI xref."""

from __future__ import annotations

import logging
import time
import urllib.parse
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

# Statuses a bad *cell value* can genuinely provoke, so retrying is pointless and
# the answer is about the value: PUG REST returns 400 for input it cannot parse,
# 414 for a URL made too long by one, and 422 for one it parses but cannot use.
#
# Everything else 4xx is deliberately excluded, because no cell can cause it.
# 401/403/407 is a blocked client; 405 is the endpoint refusing the method; 410
# is the endpoint retired; 413 on a bodyless GET is a proxy, since URL length is
# 414. Classing any of those as "no such compound" would mark every value in the
# column missing, leave `outages` empty because nothing was unreachable, and
# trip `suspect_columns` — sending the operator to fix a column that was never
# wrong. They fall through to the retry-then-unreachable path instead.
MALFORMED_REQUEST_CODES = frozenset({400, 414, 422})


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
                # Logged inside the guard so no message promises a delay that is
                # not about to happen — the final attempt does not wait.
                if attempt < MAX_RETRIES:
                    log.warning(
                        "Rate limited (%s), waiting %ss [%s/%s]",
                        resp.status_code,
                        delay,
                        attempt,
                        MAX_RETRIES,
                    )
                    time.sleep(delay)
                delay *= 2
            elif resp.status_code in MALFORMED_REQUEST_CODES:
                # These say the *request* was wrong, which is a statement about
                # the cell: PUG REST returns 400 PUGREST.BadRequest for input it
                # cannot parse and 414 for an over-long URL, which is exactly what
                # a column of free-text notes produces. Retrying gets the same
                # reply twice more, and calling it unreachable would blame the
                # network for a bad cell — and suppress the wrong-column
                # diagnosis, since an outage disqualifies suspect_columns.
                #
                # Deliberately *not* every 4xx. A 401/403 is a blocked or
                # unauthenticated client, which says nothing about the compound;
                # mapping it here would report a whole run of good CAS numbers as
                # "no such compound" and send the operator to check their column
                # while the real cause was access. Those fall through to the
                # retry-then-unreachable path below, where chebi_terms also puts
                # them (it raises ChebiUnavailableError rather than reporting a
                # miss).
                #
                # chebi_terms treats 400 as unavailable rather than as a miss,
                # and that is right *there*: its URL is built from an ID already
                # reduced to digits, so a 400 cannot have been caused by a sheet
                # cell and must be the server. Here the path carries free text
                # straight from the column, so a 400 usually is the cell.
                log.warning("HTTP %s (not retryable): %s", resp.status_code, url)
                return None, OUTCOME_NOT_FOUND
            else:
                log.warning(
                    "HTTP %s [%s/%s]: %s",
                    resp.status_code,
                    attempt,
                    MAX_RETRIES,
                    url,
                )
                if attempt < MAX_RETRIES:
                    time.sleep(delay)
                delay *= 2
        except requests.exceptions.RequestException as exc:
            log.warning(
                "Request error: %s [%s/%s]",
                exc,
                attempt,
                MAX_RETRIES,
            )
            # No sleep after the final attempt: the backoff exists to space out
            # the *next* try, and there isn't one. chebi_terms guards this too;
            # it is 8 of the ~14s an exhausted request otherwise costs, paid on
            # every row of a run against a dead upstream.
            if attempt < MAX_RETRIES:
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

    The CAS is quoted here rather than by the caller, since it reaches the URL
    path and these values come from untrusted sheets: a stray "/" or "?" would
    rewrite or truncate the request.
    """
    resp, outcome = get_with_retry_status(
        f"{BASE}/compound/name/{urllib.parse.quote(str(cas), safe='')}/cids/JSON"
    )
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return None, outcome
    try:
        # Subscripted, not .get()-with-a-default: a body that parses as JSON but
        # carries no IdentifierList at all is a different server answering
        # (`{"Fault": ...}`), not PubChem saying "no such CAS". Defaulting it to
        # [] would report NOT_FOUND, and five such values would then blame the
        # CAS column for what the network did — the very thing this function's
        # outcome exists to prevent. Only an *empty* list is a real answer.
        cids = resp.json()["IdentifierList"]["CID"]
        if not isinstance(cids, list):
            raise TypeError("CID is not a list")
        if len(cids) > 1:
            log.debug(
                "CAS %s → %s CIDs, using first (canonical): %s",
                cas,
                len(cids),
                cids[0],
            )
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
