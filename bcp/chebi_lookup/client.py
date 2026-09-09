"""PubChem PUG REST client: CAS Registry Number → CID → properties + ChEBI xref."""

from __future__ import annotations

import logging
import time
import urllib.parse
from typing import Any

import requests

from cas_registry import CAS_INVALID_FORMAT, CAS_MISSING, classify_cas

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
        # PubChem renamed these: IsomericSMILES -> SMILES (stereo-bearing) and
        # CanonicalSMILES -> ConnectivitySMILES (connectivity only). The old names
        # are still ACCEPTED in the request URL, but the response comes back keyed
        # by the new ones, so asking with the old names and reading with them
        # returned "" on every lookup without any error, non-200 or log line.
        "SMILES",
        "ConnectivitySMILES",
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
    # What was actually asked, and what was wrong with what the cell held.
    # cas_registry repairs mechanical corruption before the request, and a
    # repair applied without being reported is the one thing that module refuses to
    # do: a curator has to be able to see that the value looked up is not the value
    # in the spreadsheet. cas_queried is empty when no request was made.
    "cas_queried",
    "cas_class",
    "cas_repairs",
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
    url: str,
    params: dict | None = None,
    *,
    malformed_is_answer: bool = False,
) -> tuple[requests.Response | None, str]:
    """
    GET with exponential backoff on 429/503, reporting *why* it returned nothing.

    Same request behaviour as get_with_retry — which is now a wrapper over this —
    but a 404 comes back as OUTCOME_NOT_FOUND and exhausted retries as
    OUTCOME_UNREACHABLE. PubChem answering "no such compound" is evidence about the
    compound; PubChem never answering is not, and a caller that cannot tell them
    apart has to treat every outage as a finding.

    `malformed_is_answer` says this URL carries a value taken from a spreadsheet,
    so a status meaning "your request was wrong" is a statement about that value.
    Only the callers that interpolate sheet text set it. On a URL built from an
    id this code produced itself — a CID, say — a 400 cannot be the cell's fault
    and must be the server, so it retries and reports UNREACHABLE like any other
    server fault. Getting that backwards would report a good CAS column as full
    of unknown compounds and send the operator to fix the column.
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
            elif malformed_is_answer and resp.status_code in MALFORMED_REQUEST_CODES:
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


def _cas_for_lookup(raw: str) -> tuple[str, str, str]:
    """
    `(value to send, class, repairs)` for a spreadsheet's CAS cell.

    The value is `""` when the cell is not worth a request -- blank, or nothing a
    mechanical repair can turn into a registry number. Returned that way rather than
    as a class the caller has to interpret, so "do not ask" is one truth test and the
    `cas_queried` column is literally what was asked.

    CAS validation lives in `cas_registry`, a top-level module that belongs to
    neither package, so both can import it at module scope without a cycle.
    """
    queried, cas_class, repairs = classify_cas(raw)
    if cas_class in (CAS_MISSING, CAS_INVALID_FORMAT):
        return "", cas_class, repairs
    return queried, cas_class, repairs


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

    **The value is validated and mechanically repaired before it is sent**, by
    cas_registry, and this is the one place both call paths pass through --
    chebi_lookup.lookup_cas() and structure_check.client.cas_structure(). It matters here
    more than anywhere else because the endpoint is PubChem's **name** endpoint: it
    resolves anything, so an unchecked cell holding a compound name resolves *as a
    name*, and a row's CAS side then agrees with its name side because both asked
    the same question. A cell no repair can turn into a registry number reports
    NOT_FOUND without a request -- not UNREACHABLE, since nothing was unreachable.

    A well-formed number whose check digit disagrees is still sent. It is not a
    registry number, but PubChem indexes vendor synonyms verbatim, so what it
    answers is evidence about the row; refusing to ask would be a claim about
    PubChem's index rather than about the number. None of the 11 such numbers in the
    second batch resolved, so this costs requests and loses nothing.
    """
    queried, cas_class, repairs = _cas_for_lookup(cas)
    if not queried:
        log.warning(
            "Not a CAS Registry Number, not asking PubChem: %r (%s)", cas, cas_class
        )
        return None, OUTCOME_NOT_FOUND
    if repairs:
        log.warning("CAS %r repaired to %r before lookup (%s)", cas, queried, repairs)

    resp, outcome = get_with_retry_status(
        f"{BASE}/compound/name/{urllib.parse.quote(queried, safe='')}/cids/JSON",
        malformed_is_answer=True,
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
        if cids and not isinstance(cids[0], int):
            # The return is annotated int | None and is interpolated unquoted
            # into a URL path, so anything else is a malformed payload.
            raise TypeError("CID is not an integer")
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
            result["isomeric_smiles"] = props.get("SMILES", "")
            result["canonical_smiles"] = props.get("ConnectivitySMILES", "")
            result["inchi"] = props.get("InChI", "")
            result["inchikey"] = props.get("InChIKey", "")
            result["xlogp"] = props.get("XLogP", "")
            result["tpsa"] = props.get("TPSA", "")
        except (ValueError, KeyError, AttributeError, TypeError, IndexError) as exc:
            # A 200 whose body is not this endpoint's JSON is an outage wearing a
            # success code, exactly as cas_to_cid_status documents. This function
            # returns no outcome, so the least it can do is refuse to be silent:
            # swallowing it bare would leave the row indistinguishable from
            # "PubChem has no properties for this CID".
            log.warning("Unusable property payload for CID %s: %s", cid, exc)

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
        except (ValueError, KeyError, AttributeError, TypeError, IndexError) as exc:
            # Same reasoning: without this line an unusable xref body reads as
            # "this compound has no ChEBI cross-reference", which is a finding
            # about chemistry rather than about the network.
            log.warning("Unusable registry-ID payload for CID %s: %s", cid, exc)

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
        # str() each synonym: PubChem has shipped numeric-looking synonyms as bare
        # numbers, and join() raises TypeError on the first one it meets. Widen the
        # except tuple too -- a valid-JSON body of the wrong shape (a proxy or
        # maintenance page) raises AttributeError here, not ValueError.
        return " | ".join(str(s) for s in syns[:MAX_SYNONYMS])
    except (ValueError, KeyError, AttributeError, TypeError, IndexError) as exc:
        log.warning("Unusable synonym payload for CID %s: %s", cid, exc)
        return ""


def lookup_cas(cas: str) -> dict[str, Any]:
    """
    Resolve one CAS to PubChem properties, ChEBI xref, and synonyms.

    Every row carries `cas_class` and `cas_repairs` whatever the outcome, so an
    empty result says which of "the cell was blank", "no repair made this a registry
    number" and "PubChem does not index it" happened. `cas_queried` is the value
    actually sent, empty when nothing was.
    """
    cas = (cas or "").strip()
    queried, cas_class, repairs = _cas_for_lookup(cas)
    report = {"cas_queried": queried, "cas_class": cas_class, "cas_repairs": repairs}
    if not queried:
        # cas_to_cid_status refuses these too. Returning here keeps the refusal out of
        # the log twice and puts the reason in the row rather than only in the log.
        return {**empty_result(), **report}

    cid = cas_to_cid(queried)
    if cid is None:
        return {**empty_result(), **report}

    result = cid_to_properties(cid)
    result["synonyms"] = cid_to_synonyms(cid)
    return {**result, **report}
