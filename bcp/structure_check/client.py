"""
Cross-check a curation row's identifiers by structure rather than by name.

Comparing chemical names as strings does not work: ChEBI calls CHEBI:92401
"6-(1,3-dioxo-2-benzo[de]isoquinolinyl)-N-hydroxyhexanamide" where a sheet says
"Scriptaid", and both are correct. InChIKeys do work, because they are derived
from the structure itself.

So each identifier in a row is resolved to a structure independently — the CAS
number and the compound name through PubChem, the ChEBI ID through ChEBI — and
those structures are compared. Two questions get separate answers:

  id_cas_verdict    does the ChEBI ID hold the same molecule as the CAS number?
                    ("has ChEBI got something else under this ID?")
  name_cas_verdict  is the CAS number right for the compound the row claims to be?

PubChem access is reused from chebi_lookup and ChEBI access from chebi_terms, so
this module adds no new HTTP client and inherits their retry and rate-limit
behaviour.
"""

from __future__ import annotations

import logging
import re
import time
import urllib.parse
from typing import Any

from chebi_lookup.client import BASE, REQUEST_DELAY, cas_to_cid, get_with_retry
from chebi_terms.client import (
    ChebiUnavailableError,
    fetch_compound,
    normalize_chebi_id,
)

log = logging.getLogger(__name__)

# An InChIKey is SKELETON(14) "-" STEREO+FLAGS(10) "-" PROTONATION(1).
# The first block hashes molecular connectivity, so a difference there is a
# genuinely different molecule. A difference only after it is the same skeleton
# with different stereochemistry, isotopes, or charge state — worth a look, but
# not the same class of error.
SKELETON_LEN = 14

# Underscore-separated tokens shorter than this, or with no letters, are not
# plausible compound names ("263", "28", "D").
MIN_TOKEN_LEN = 3

# A loosely named cell can resolve to many CIDs, and each one costs three more
# requests in refine_skeleton_difference. Unbounded per-row cost shows up as a
# run that never finishes, so cap it and say so when truncating.
MAX_NAME_CANDIDATES = 10

# Counterions, salt forms, and hydration states. Their presence means the cell is
# one compound plus its form ("PIK-75_HCl"), not two aliases ("Vorinostat_SAHA"),
# so the tokens must never be queried separately: PIK-75 and PIK-75 HCl have
# different InChIKeys, and comparing the free base would manufacture a difference
# that is purely our own doing.
SALT_TOKENS = frozenset(
    {
        "hcl",
        "hbr",
        "hydrochloride",
        "dihydrochloride",
        "chloride",
        "dichloride",
        "bromide",
        "dibromide",
        "iodide",
        "sodium",
        "disodium",
        "potassium",
        "calcium",
        "magnesium",
        "lithium",
        "ammonium",
        "sulfate",
        "sulphate",
        "bisulfate",
        "phosphate",
        "diphosphate",
        "nitrate",
        "acetate",
        "citrate",
        "tartrate",
        "bitartrate",
        "maleate",
        "fumarate",
        "succinate",
        "malate",
        "oxalate",
        "mesylate",
        "besylate",
        "tosylate",
        "triflate",
        "pamoate",
        "palmitate",
        "stearate",
        "gluconate",
        "lactate",
        "salicylate",
        "benzoate",
        "hydrate",
        "dihydrate",
        "trihydrate",
        "semihydrate",
        "hemihydrate",
        "monohydrate",
        "anhydrous",
        "free",
        "base",
        "salt",
        "hydrobromide",
        "hydroiodide",
        "trifluoroacetate",
        "embonate",
        "tfa",
    }
)

# Multiplicity prefixes a counterion can carry: "2HCl", "3HCl", "dihydrate".
_COUNT_PREFIXES = (
    "sesqui",
    "hemi",
    "mono",
    "di",
    "tri",
    "tetra",
    "penta",
    "hexa",
    "bis",
    "tris",
)

# Counterion words common enough as a *suffix* that matching the whole fragment
# would miss them ("trifluoroacetate", "hydroiodide", and any "…HCl" spelling).
_SALT_SUFFIXES = (
    "hcl",
    "hbr",
    "tfa",
    "hydrochloride",
    "hydrobromide",
    "hydroiodide",
    "trifluoroacetate",
    "acetate",
    "sulfate",
    "sulphate",
    "phosphate",
    "maleate",
    "mesylate",
    "besylate",
    "tosylate",
    "citrate",
    "tartrate",
    "fumarate",
    "succinate",
    "oxalate",
    "pamoate",
    "embonate",
    "hydrate",
    "chloride",
    "bromide",
    "iodide",
)


def _is_salt_word(fragment: str) -> bool:
    """
    True when a name fragment denotes a counterion, salt form, or hydration state.

    Deliberately generous. Firing wrongly only disables the token fallback, whose
    worst outcome is name_unresolved; failing to fire lets the fallback query a
    free base against a salt's CAS and invent a difference.
    """
    word = "".join(fragment.split()).casefold().lstrip("0123456789")
    if not word:
        return False
    if word in SALT_TOKENS:
        return True
    if any(
        word.startswith(prefix) and word[len(prefix) :] in SALT_TOKENS
        for prefix in _COUNT_PREFIXES
    ):
        return True
    return word.endswith(_SALT_SUFFIXES)


# Structure comparison outcomes.
MATCH = "match"
STEREO_DIFFERS = "stereo_differs"
SALT_DIFFERS = "salt_differs"
SKELETON_DIFFERS = "skeleton_differs"
# Neither a match nor a difference: there was nothing on one side to compare.
NOT_COMPARABLE = "not_comparable"

# Reasons a comparison could not be made. None of these is a finding about the
# compound — the distinction chebi_terms draws between an answer and a silence.
NAME_UNRESOLVED = "name_unresolved"
# The name resolved to more structures than were compared, and none matched. Not a
# difference: the evidence is incomplete, so saying "different molecule" would be a
# finding manufactured by our own truncation.
NAME_AMBIGUOUS = "name_ambiguous"
CAS_UNRESOLVED = "cas_unresolved"
# Kept apart so a run-wide outage is distinguishable from a sheet of bad IDs. The
# ChEBI side can draw this line because fetch_compound raises rather than
# returning None; the PubChem side cannot.
CHEBI_UNREACHABLE = "chebi_unreachable"
CHEBI_UNRESOLVED = "chebi_unresolved"
CHEBI_NO_STRUCTURE = "chebi_no_structure"
NOT_CHECKED = "not_checked"

UNRESOLVED_VERDICTS = (
    NAME_UNRESOLVED,
    NAME_AMBIGUOUS,
    CAS_UNRESOLVED,
    CHEBI_UNREACHABLE,
    CHEBI_UNRESOLVED,
    CHEBI_NO_STRUCTURE,
)

# The single sortable column: what a curator should do with this row.
REVIEW_INVESTIGATE = "investigate"
REVIEW_CHECK = "check"
REVIEW_OK = "ok"
REVIEW_UNVERIFIED = "unverified"

# Alphabetically "check" < "investigate", so sorting the review column alone puts
# the softer flag first. review_rank exists so "sort by it" is literally true.
REVIEW_RANK = {
    REVIEW_INVESTIGATE: 1,
    REVIEW_CHECK: 2,
    REVIEW_UNVERIFIED: 3,
    REVIEW_OK: 4,
}

OUTPUT_FIELDS_APPENDED = [
    "review_rank",
    "review",
    "id_cas_verdict",
    "name_cas_verdict",
    "cas_inchikey",
    "cas_pubchem_name",
    "chebi_inchikey",
    "name_query",
    "name_inchikey",
]


def empty_result() -> dict[str, Any]:
    """A result row with both verdicts unchecked and nothing resolved."""
    result: dict[str, Any] = {field: "" for field in OUTPUT_FIELDS_APPENDED}
    result["id_cas_verdict"] = NOT_CHECKED
    result["name_cas_verdict"] = NOT_CHECKED
    result["review"] = REVIEW_UNVERIFIED
    result["review_rank"] = REVIEW_RANK[REVIEW_UNVERIFIED]
    return result


def skeleton(inchikey: str) -> str:
    """The connectivity block of an InChIKey."""
    return (inchikey or "")[:SKELETON_LEN]


def compare_structures(reference: str, candidates: list[str]) -> str:
    """
    Compare one InChIKey against every key an identifier resolved to.

    A match against any candidate is a match: a name that resolves to several
    entries (a free base and its salts, say) is not evidence of an error. A stray
    candidate can therefore only mask a finding, never invent one — the safe
    direction for a tool whose output drives curation edits.

    Returns NOT_COMPARABLE when either side is missing. That is not a difference,
    and callers must not render it as one — this function is public, so it says so
    rather than leaving the ambiguity to whoever reads the return value.
    """
    present = [key for key in candidates if key]
    if not present or not reference:
        return NOT_COMPARABLE
    if reference in present:
        return MATCH
    if any(skeleton(key) == skeleton(reference) for key in present):
        return STEREO_DIFFERS
    return SKELETON_DIFFERS


_parent_cache: dict[str, str] = {}


def parent_inchikey(inchikey: str) -> str:
    """
    PubChem's desalted parent structure for an InChIKey, or "" if unavailable.

    A counterion is part of an InChIKey's connectivity block, so doxorubicin and
    doxorubicin hydrochloride hash to different skeletons even though a sheet
    naming one and citing the other is ordinary looseness rather than an error.
    Comparing parents collapses that difference.

    Used only to *downgrade* a skeleton difference, never to claim a match:
    PubChem's parent assignment is unreliable for multi-component salts — for
    pyrvinium pamoate it returns the pamoic acid counterion — so a wrong answer
    here can only leave a row flagged, never wave one through.
    """
    if not inchikey:
        return ""
    if inchikey in _parent_cache:
        return _parent_cache[inchikey]

    result = ""
    key = urllib.parse.quote(inchikey, safe="")
    resp = get_with_retry(f"{BASE}/compound/inchikey/{key}/cids/JSON")
    time.sleep(REQUEST_DELAY)
    cid = _first_cid(resp)
    if cid is not None:
        resp = get_with_retry(f"{BASE}/compound/cid/{cid}/cids/JSON?cids_type=parent")
        time.sleep(REQUEST_DELAY)
        parent_cid = _first_cid(resp)
        if parent_cid is not None:
            resp = get_with_retry(
                f"{BASE}/compound/cid/{parent_cid}/property/InChIKey/JSON"
            )
            time.sleep(REQUEST_DELAY)
            if resp is not None:
                try:
                    result = resp.json()["PropertyTable"]["Properties"][0]["InChIKey"]
                except (ValueError, KeyError, IndexError, TypeError):
                    result = ""

    # Only successes are cached: one transient failure must not disable the salt
    # demotion for this key for the rest of the run.
    if result:
        _parent_cache[inchikey] = result
    return result


def _first_cid(resp: Any) -> int | None:
    if resp is None:
        return None
    try:
        cids = resp.json()["IdentifierList"]["CID"]
    except (ValueError, KeyError, TypeError):
        return None
    return cids[0] if cids else None


def refine_skeleton_difference(reference: str, candidates: list[str]) -> str:
    """
    Decide whether a skeleton difference is only a salt-form difference.

    Returns SALT_DIFFERS when the two sides share a desalted parent, otherwise
    leaves SKELETON_DIFFERS. Costs three requests per structure and so is only
    worth calling once a difference has already been found.
    """
    reference_parent = parent_inchikey(reference)
    if not reference_parent:
        return SKELETON_DIFFERS
    for candidate in candidates:
        if candidate and parent_inchikey(candidate) == reference_parent:
            return SALT_DIFFERS
    return SKELETON_DIFFERS


def name_candidates(raw: Any) -> tuple[list[str], list[str]]:
    """
    Query strings to try for a name, as (whole_string_forms, fallback_tokens).

    Salt and suffix information is never stripped. "PIK-75 HCl" and "PIK-75"
    resolve to different InChIKeys, so dropping a "_HCl" would quietly compare the
    wrong chemistry and report a difference that is entirely our own doing.

    Tokens are a fallback only, for the "Vorinostat_SAHA" convention where a sheet
    packs two aliases for one compound into one cell. A cell naming a salt form is
    not that, so any SALT_TOKENS match disables the fallback entirely: reporting
    name_unresolved is better than comparing a molecule the row never claimed.
    """
    text = " ".join(str(raw or "").split())
    if not text:
        return [], []

    whole = [text]
    spaced = text.replace("_", " ")
    if spaced != text:
        whole.append(spaced)

    # Guard on both separators: "Doxorubicin HCl_Adriamycin" hides its counterion
    # behind a space, and splitting only on "_" would miss it.
    if any(_is_salt_word(piece) for piece in re.split(r"[_\s]+", text)):
        return whole, []

    tokens = [
        token
        for token in text.split("_")
        if len(token) >= MIN_TOKEN_LEN and any(char.isalpha() for char in token)
    ]
    # A lone token is already covered by the whole-string attempts above.
    return whole, tokens if len(tokens) > 1 else []


def inchikeys_for_name(name: str) -> list[str]:
    """Every InChIKey PubChem resolves a name to. Empty when it resolves none."""
    query = urllib.parse.quote(name, safe="")
    resp = get_with_retry(f"{BASE}/compound/name/{query}/property/InChIKey/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return []
    try:
        properties = resp.json()["PropertyTable"]["Properties"]
    except (ValueError, KeyError, TypeError):
        return []
    keys = []
    for entry in properties:
        key = entry.get("InChIKey") if isinstance(entry, dict) else None
        if key and key not in keys:
            keys.append(key)
    return keys


def name_structure(raw: Any) -> tuple[str, list[str], int]:
    """
    Resolve a name to (query_actually_used, inchikeys, total_found).

    Whole-string forms are tried first and the first hit wins, so a salt form stays
    intact. Only when none resolve are tokens tried, and then their results are
    unioned: "Vorinostat_SAHA" means one compound under two aliases, and a token
    that resolves to something unrelated can only add a key.

    Returns (query_actually_used, capped_keys, total_found). query_actually_used is
    reported so the comparison is auditable — a difference should never be
    traceable to a query the caller cannot see — and total_found lets the caller
    tell a complete comparison from a truncated one.
    """
    whole, tokens = name_candidates(raw)

    for candidate in whole:
        keys = inchikeys_for_name(candidate)
        if keys:
            return candidate, _cap(candidate, keys), len(keys)

    union: list[str] = []
    used: list[str] = []
    for token in tokens:
        keys = inchikeys_for_name(token)
        if keys:
            used.append(token)
            union.extend(key for key in keys if key not in union)
    if union:
        query = "tokens: " + "|".join(used)
        return query, _cap(query, union), len(union)

    return "", [], 0


def _cap(query: str, keys: list[str]) -> list[str]:
    """Limit how many structures one cell can cost, loudly."""
    if len(keys) <= MAX_NAME_CANDIDATES:
        return keys
    log.warning(
        "Name %r resolved to %s structures; comparing the first %s.",
        query,
        len(keys),
        MAX_NAME_CANDIDATES,
    )
    return keys[:MAX_NAME_CANDIDATES]


def chebi_structure(chebi_id: Any) -> tuple[str, str]:
    """
    Resolve a ChEBI ID to (inchikey, problem), where problem is "" when fine.

    CHEBI_UNREACHABLE means ChEBI could not be asked; CHEBI_UNRESOLVED means it was
    asked and has no such record (or the ID was malformed). Both mean "no structure
    to compare" for this row, but only the first one being universal implies an
    outage, so the run-level check can tell them apart. CHEBI_NO_STRUCTURE means
    ChEBI has the entry and records no structure, as class terms and R-group
    entries do not.
    """
    parsed = normalize_chebi_id(chebi_id)
    if parsed is None:
        return "", CHEBI_UNRESOLVED
    try:
        payload = fetch_compound(parsed[0])
    except ChebiUnavailableError as exc:
        log.error("%s", exc)
        return "", CHEBI_UNREACHABLE
    if payload is None:
        return "", CHEBI_UNRESOLVED
    key = (payload.get("default_structure") or {}).get("standard_inchi_key") or ""
    return (key, "") if key else ("", CHEBI_NO_STRUCTURE)


def cas_structure(cas: Any) -> tuple[str, str]:
    """
    Resolve a CAS number to (inchikey, pubchem_preferred_name).

    Two requests, not chebi_lookup.lookup_cas()'s four: only the InChIKey and Title
    are wanted, and fetching xrefs and synonyms for every row spent about a minute
    of REQUEST_DELAY on a 117-row sheet producing data that was discarded.

    Caveat inherited from chebi_lookup.get_with_retry: it returns None for a 404
    *and* for exhausted retries, so an empty result here cannot distinguish
    "PubChem has no such CAS" from "PubChem was unreachable". The ChEBI side does
    draw that distinction; this side does not, which is why a run where nothing
    resolved is reported as degraded rather than as a clean sheet.
    """
    text = str(cas or "").strip()
    if not text:
        return "", ""
    # chebi_lookup.cas_to_cid interpolates straight into the URL path, and these
    # cells come from the same untrusted sheet as the names: a stray "/" or "?"
    # would rewrite or truncate the request. Every other URL here is quoted too.
    cid = cas_to_cid(urllib.parse.quote(text, safe=""))
    if cid is None:
        return "", ""
    resp = get_with_retry(f"{BASE}/compound/cid/{cid}/property/InChIKey,Title/JSON")
    time.sleep(REQUEST_DELAY)
    if resp is None:
        return "", ""
    try:
        properties = resp.json()["PropertyTable"]["Properties"][0]
    except (ValueError, KeyError, IndexError, TypeError):
        return "", ""
    return properties.get("InChIKey") or "", properties.get("Title") or ""


def review_level(id_cas_verdict: str, name_cas_verdict: str) -> str:
    """
    Collapse both verdicts into one sortable answer.

    Ranked by what it would cost to be wrong: a different skeleton means the row
    points at another molecule, so it outranks a stereochemical difference, which
    outranks a check that simply could not be made.
    """
    verdicts = (id_cas_verdict, name_cas_verdict)
    if SKELETON_DIFFERS in verdicts:
        return REVIEW_INVESTIGATE
    if STEREO_DIFFERS in verdicts or SALT_DIFFERS in verdicts:
        return REVIEW_CHECK
    if MATCH in verdicts:
        return REVIEW_OK
    return REVIEW_UNVERIFIED


def check_row(
    *,
    cas: Any = None,
    chebi_id: Any = None,
    name: Any = None,
) -> dict[str, Any]:
    """
    Cross-check one row's identifiers by structure.

    The CAS number is the pivot for both comparisons, so without it neither can be
    made. Never raises: anything unresolvable is reported as such rather than as a
    finding about the compound.
    """
    result = empty_result()

    wants_id_check = bool(str(chebi_id or "").strip())
    wants_name_check = bool(str(name or "").strip())
    if not wants_id_check and not wants_name_check:
        # Nothing to compare the CAS against, so resolving it would cost requests
        # for a result that is discarded either way.
        return result

    cas_key, cas_name = cas_structure(cas)
    result["cas_inchikey"] = cas_key
    result["cas_pubchem_name"] = cas_name
    # A blank cell was never asked about; only a non-blank CAS that PubChem
    # could not resolve is a failed check. cas_structure() cannot draw this line
    # itself (it returns ("", "") for both), so the caller does it here.
    cas_verdict_if_missing = (
        NOT_CHECKED if not str(cas or "").strip() else CAS_UNRESOLVED
    )

    if wants_id_check:
        chebi_key, problem = chebi_structure(chebi_id)
        result["chebi_inchikey"] = chebi_key
        if problem:
            result["id_cas_verdict"] = problem
        elif not cas_key:
            result["id_cas_verdict"] = cas_verdict_if_missing
        else:
            verdict = compare_structures(cas_key, [chebi_key])
            if verdict == SKELETON_DIFFERS:
                verdict = refine_skeleton_difference(cas_key, [chebi_key])
            result["id_cas_verdict"] = verdict

    if wants_name_check:
        query, name_keys, total_found = name_structure(name)
        truncated = total_found > len(name_keys)
        if truncated:
            query = f"{query} (truncated: {total_found} structures)"
        result["name_query"] = query
        result["name_inchikey"] = " | ".join(name_keys)
        if not name_keys:
            result["name_cas_verdict"] = NAME_UNRESOLVED
        elif not cas_key:
            result["name_cas_verdict"] = cas_verdict_if_missing
        else:
            verdict = compare_structures(cas_key, name_keys)
            if verdict == SKELETON_DIFFERS:
                verdict = refine_skeleton_difference(cas_key, name_keys)
            # PubChem returns structures in CID order, not relevance order, so the
            # true match may be among the ones we did not compare. Reporting a
            # difference here would be a finding of our own making.
            if truncated and verdict != MATCH:
                verdict = NAME_AMBIGUOUS
            result["name_cas_verdict"] = verdict

    result["review"] = review_level(
        result["id_cas_verdict"], result["name_cas_verdict"]
    )
    result["review_rank"] = REVIEW_RANK[result["review"]]
    return result
