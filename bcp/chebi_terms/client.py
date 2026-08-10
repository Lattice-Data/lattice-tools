"""ChEBI REST client: ChEBI ID → authoritative name, synonyms, and validation verdicts."""

from __future__ import annotations

import html
import logging
import re
import time
from typing import Any, Iterator

import requests

# The ChEBI backend publishes no documented rate limit. Mirror the polite
# spacing chebi_lookup uses for PubChem rather than hammering it. Kept local
# rather than shared with chebi_lookup so tuning one tool's retry behaviour
# cannot silently change the other's.
REQUEST_DELAY = 0.25
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0

BASE = "https://www.ebi.ac.uk/chebi/backend/api/public"

# Name types emitted in chebi_synonyms. BRAND NAME and INN are deliberately
# excluded from output (acetylsalicylic acid alone carries 64 brand names) but
# are still consulted when matching a supplied name.
SYNONYM_TYPES = ("SYNONYM", "IUPAC NAME")
SYNONYM_LANGUAGE = "en"

log = logging.getLogger(__name__)

OUTPUT_FIELDS_APPENDED = [
    "chebi_accession",
    "chebi_name",
    "chebi_synonyms",
    "chebi_definition",
    "chebi_stars",
    "id_status",
    "name_verdict",
    "cas_verdict",
]


class ChebiUnavailableError(Exception):
    """
    ChEBI could not be reached, or gave no usable answer after MAX_RETRIES.

    Deliberately distinct from a 404: it says nothing about whether the compound
    exists, so callers must never report it as not_found.
    """


# id_status — what we learned about the compound. Highest precedence first when
# more than one could apply.
STATUS_NOT_FOUND = "not_found"
STATUS_NOT_RELEASED = "not_released"
STATUS_SECONDARY = "secondary"
STATUS_OK = "ok"

# id_status — what we learned about the input instead of the compound.
STATUS_MISSING = "missing"
STATUS_INVALID = "invalid"

# id_status — we never got to ask. Not a statement about the compound.
STATUS_LOOKUP_FAILED = "lookup_failed"

NOT_CHECKED = "not_checked"
NAME_MATCH = "match"
NAME_SYNONYM_MATCH = "synonym_match"
NAME_MISMATCH = "mismatch"
CAS_CONFIRMED = "confirmed"
CAS_NOT_IN_CHEBI = "not_in_chebi"

# ChEBI embeds inline markup in names: "(3<i>R</i>,4<i>S</i>)-...".
_TAG_RE = re.compile(r"</?[A-Za-z][^>]*>")
# Curators type an ASCII hyphen where ChEBI stores U+2212 et al: "(−)-epicatechin".
_DASH_MAP = dict.fromkeys(
    map(ord, "‐‑‒–—―−"),
    "-",
)
_ACCESSION_RE = re.compile(r"^(?:CHEBI:)?(\d+)$", re.IGNORECASE)


def empty_result() -> dict[str, Any]:
    """Return a result dict with facts blank and both name/CAS verdicts unchecked."""
    result: dict[str, Any] = {field: "" for field in OUTPUT_FIELDS_APPENDED}
    result["name_verdict"] = NOT_CHECKED
    result["cas_verdict"] = NOT_CHECKED
    return result


def clean_name(value: Any) -> str:
    """Strip ChEBI's inline HTML markup and normalize whitespace."""
    if value is None:
        return ""
    text = html.unescape(str(value))
    text = _TAG_RE.sub("", text)
    return " ".join(text.split())


def match_key(value: Any) -> str:
    """Normalize a name for comparison: markup, dash variants, case, whitespace."""
    return clean_name(value).translate(_DASH_MAP).casefold()


def normalize_chebi_id(chebi_id: Any) -> tuple[str, str] | None:
    """
    Parse a ChEBI identifier into (numeric_id, accession).

    Accepts "16236", "CHEBI:16236", and "chebi:16236". Returns None when the
    value is not a ChEBI identifier at all, so malformed input costs no request.
    """
    if chebi_id is None:
        return None
    match = _ACCESSION_RE.match(str(chebi_id).strip())
    if match is None:
        return None
    numeric = match.group(1).lstrip("0")
    if not numeric:
        return None
    return numeric, f"CHEBI:{numeric}"


def get_with_retry(url: str) -> requests.Response | None:
    """
    GET with exponential backoff on 429/503.

    Returns None for a 404 — an answer, meaning no such compound. Raises
    ChebiUnavailableError when ChEBI cannot be reached or keeps failing, so that
    "we could not ask" is never mistaken for "ChEBI says no".
    """
    delay = RETRY_BACKOFF
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, timeout=15)
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
            log.warning("Request error: %s [%s/%s]", exc, attempt, MAX_RETRIES)
            time.sleep(delay)
            delay *= 2
    raise ChebiUnavailableError(
        f"ChEBI unreachable after {MAX_RETRIES} attempts: {url}"
    )


def fetch_compound(numeric_id: str) -> dict[str, Any] | None:
    """
    Fetch the ChEBI compound payload for a numeric ID.

    Returns None only when ChEBI answers that no such compound exists. A garbled
    or unreachable ChEBI raises ChebiUnavailableError instead — a 200 carrying a
    non-JSON body is a server fault, not evidence about the compound.

    ChEBI resolves secondary IDs to their primary compound server-side, so the
    returned chebi_accession may differ from the ID requested.
    """
    resp = get_with_retry(f"{BASE}/compound/{numeric_id}/")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return None
    try:
        payload = resp.json()
    except ValueError as exc:
        raise ChebiUnavailableError(
            f"Non-JSON response for ChEBI ID {numeric_id}"
        ) from exc
    if not isinstance(payload, dict):
        raise ChebiUnavailableError(
            f"Unexpected response shape for ChEBI ID {numeric_id}: "
            f"{type(payload).__name__}"
        )
    return payload


def _iter_name_entries(payload: dict[str, Any]) -> Iterator[tuple[str, dict[str, Any]]]:
    """Yield (name_type, entry) for every name in the payload's names dict."""
    names = payload.get("names")
    if not isinstance(names, dict):
        return
    for name_type, entries in names.items():
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if isinstance(entry, dict):
                yield str(name_type).upper(), entry


def extract_synonyms(
    payload: dict[str, Any],
    *,
    max_synonyms: int | None = None,
) -> list[str]:
    """
    Collect English SYNONYM and IUPAC NAME values, deduplicated, order preserved.

    The compound's own name is omitted — it is reported separately as chebi_name.
    """
    official = match_key(payload.get("name"))
    seen = {official} if official else set()
    synonyms: list[str] = []

    for name_type, entry in _iter_name_entries(payload):
        if name_type not in SYNONYM_TYPES:
            continue
        language = entry.get("language_code")
        if language and str(language).lower() != SYNONYM_LANGUAGE:
            continue
        value = clean_name(entry.get("name"))
        if not value:
            continue
        key = match_key(value)
        if key in seen:
            continue
        seen.add(key)
        synonyms.append(value)
        if max_synonyms is not None and len(synonyms) >= max_synonyms:
            break

    return synonyms


def _name_candidates(payload: dict[str, Any]) -> set[str]:
    """
    Every name this compound is known by, as comparison keys.

    Deliberately permissive — brand names, INNs, and ChEBI's ASCII spellings all
    count as a match, so an accented or trade name in the input is not flagged.
    """
    candidates = set()
    for value in (payload.get("name"), payload.get("ascii_name")):
        key = match_key(value)
        if key:
            candidates.add(key)
    for _name_type, entry in _iter_name_entries(payload):
        for value in (entry.get("name"), entry.get("ascii_name")):
            key = match_key(value)
            if key:
                candidates.add(key)
    return candidates


def _cas_key(value: Any) -> str:
    """
    Normalize a CAS Registry Number for comparison.

    Absorbs the same Unicode dash variants match_key() handles, plus stray
    whitespace, so a curator's '64‑17‑5' is not reported as a disagreement with
    ChEBI. Internal dashes are kept: collapsing '64-17-5' and '64175' onto one
    key would risk a false 'confirmed', which is the more dangerous direction.
    """
    return "".join(clean_name(value).translate(_DASH_MAP).split())


def extract_cas_numbers(payload: dict[str, Any]) -> list[str]:
    """Collect CAS Registry Numbers that ChEBI itself records for this compound."""
    accessions = payload.get("database_accessions")
    if not isinstance(accessions, dict):
        return []
    entries = accessions.get("CAS")
    if not isinstance(entries, list):
        return []
    cas_numbers = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        value = str(entry.get("accession_number") or "").strip()
        if value and value not in cas_numbers:
            cas_numbers.append(value)
    return cas_numbers


def _id_status(requested_accession: str, payload: dict[str, Any]) -> str:
    """Classify the requested ID against the payload ChEBI returned for it."""
    if payload.get("is_released") is False:
        return STATUS_NOT_RELEASED
    returned = normalize_chebi_id(payload.get("chebi_accession"))
    requested = normalize_chebi_id(requested_accession)
    if returned is not None and requested is not None and returned[0] != requested[0]:
        return STATUS_SECONDARY
    return STATUS_OK


def verify_payload(
    requested_accession: str,
    payload: dict[str, Any] | None,
    *,
    expected_name: str | None = None,
    expected_cas: str | None = None,
    max_synonyms: int | None = None,
) -> dict[str, Any]:
    """
    Build a result row from an already-fetched payload. Pure — makes no requests.

    Splitting this out lets batch mode cache one payload per unique ID and still
    evaluate per-row expectations against it.
    """
    result = empty_result()

    if payload is None:
        result["id_status"] = STATUS_NOT_FOUND
        return result

    result["chebi_accession"] = clean_name(payload.get("chebi_accession"))
    result["chebi_name"] = clean_name(payload.get("name"))
    result["chebi_definition"] = clean_name(payload.get("definition"))
    stars = payload.get("stars")
    result["chebi_stars"] = "" if stars is None else stars
    result["chebi_synonyms"] = " | ".join(
        extract_synonyms(payload, max_synonyms=max_synonyms)
    )
    result["id_status"] = _id_status(requested_accession, payload)

    if expected_name and expected_name.strip():
        expected_key = match_key(expected_name)
        if expected_key and expected_key == match_key(payload.get("name")):
            result["name_verdict"] = NAME_MATCH
        elif expected_key and expected_key in _name_candidates(payload):
            result["name_verdict"] = NAME_SYNONYM_MATCH
        else:
            result["name_verdict"] = NAME_MISMATCH

    if expected_cas and expected_cas.strip():
        expected_key = _cas_key(expected_cas)
        known = {_cas_key(cas) for cas in extract_cas_numbers(payload)}
        if expected_key and expected_key in known:
            result["cas_verdict"] = CAS_CONFIRMED
        else:
            result["cas_verdict"] = CAS_NOT_IN_CHEBI

    return result


def verify_chebi_id(
    chebi_id: str,
    *,
    expected_name: str | None = None,
    expected_cas: str | None = None,
    max_synonyms: int | None = None,
) -> dict[str, Any]:
    """
    Resolve one ChEBI ID and check it against an expected name and/or CAS.

    Never raises: an unreachable ChEBI is reported as id_status lookup_failed
    rather than as a claim about the compound.
    """
    parsed = normalize_chebi_id(chebi_id)
    if parsed is None:
        result = empty_result()
        blank = not str(chebi_id or "").strip()
        result["id_status"] = STATUS_MISSING if blank else STATUS_INVALID
        return result

    numeric, accession = parsed
    try:
        payload = fetch_compound(numeric)
    except ChebiUnavailableError as exc:
        log.error("%s", exc)
        result = empty_result()
        result["id_status"] = STATUS_LOOKUP_FAILED
        return result

    return verify_payload(
        accession,
        payload,
        expected_name=expected_name,
        expected_cas=expected_cas,
        max_synonyms=max_synonyms,
    )


def describe_chebi_id(
    chebi_id: str,
    *,
    max_synonyms: int | None = None,
) -> dict[str, Any]:
    """Resolve one ChEBI ID to its name, synonyms, and definition — no expectations."""
    return verify_chebi_id(chebi_id, max_synonyms=max_synonyms)
