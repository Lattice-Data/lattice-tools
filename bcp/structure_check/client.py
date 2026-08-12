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
from collections.abc import Collection
from typing import Any

from chebi_lookup.client import (
    BASE,
    OUTCOME_UNREACHABLE,
    REQUEST_DELAY,
    cas_to_cid_status,
    get_with_retry_status,
)
from chebi_terms.client import (
    ChebiUnavailableError,
    fetch_compound,
    normalize_chebi_id,
)

log = logging.getLogger(__name__)


class PubChemUnavailableError(RuntimeError):
    """
    PubChem could not be asked — not "PubChem has no such compound".

    The mirror of chebi_terms.ChebiUnavailableError, and raised for the same
    reason: a resolver that returns "nothing" for both cases forces every caller
    to treat an outage as evidence about the compound. Raised rather than
    returned so the signal can be added without widening the return tuple of
    every resolver between here and check_row.
    """


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
        # Counterions written as formulae rather than words: "Diclofenac_2Na",
        # "Foo_Na2", "Bar_NaCl". Spelled-out forms were covered from the start,
        # and these read as ordinary tokens without them — so the fallback would
        # query the free base against a salt's CAS, the difference this guard
        # exists to stop the tool inventing.
        "na",
        "na2",
        "k",
        "k2",
        "ca",
        "mg",
        "li",
        "cl",
        "br",
        "nacl",
        "kcl",
        "hi",
        "tsoh",
        "msoh",
        "oh",
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

    Deliberately generous, and wider than "counterion" strictly implies: after
    digits are stripped, bare element symbols match too, so "Vitamin_K3" and
    "CA_074_Me" also lose the fallback. That is the safe direction — firing
    wrongly only disables the token fallback, whose worst outcome is
    name_unresolved, and the whole-string forms still run — while failing to fire
    lets the fallback query a free base against a salt's CAS and invent a
    difference.
    """
    # Stripped from both ends: a multiplicity prefix ("2Na") and a stoichiometry
    # suffix ("Na2") are the same counterion written two ways.
    word = "".join(fragment.split()).casefold().strip("0123456789")
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
# Kept apart so a run-wide outage is distinguishable from a sheet of bad IDs. Both
# sides can draw this line now: fetch_compound raises, and the PubChem side raises
# PubChemUnavailableError off get_with_retry_status.
CHEBI_UNREACHABLE = "chebi_unreachable"
CAS_UNREACHABLE = "cas_unreachable"
NAME_UNREACHABLE = "name_unreachable"
CHEBI_UNRESOLVED = "chebi_unresolved"
CHEBI_NO_STRUCTURE = "chebi_no_structure"
NOT_CHECKED = "not_checked"

# An upstream was asked and never answered. The only verdicts that say something
# about the *network* rather than about the compound, which is why they are the
# only ones a run-level outage verdict is allowed to count.
OUTAGE_VERDICTS = (
    CHEBI_UNREACHABLE,
    CAS_UNREACHABLE,
    NAME_UNREACHABLE,
)

UNRESOLVED_VERDICTS = (
    NAME_UNRESOLVED,
    NAME_AMBIGUOUS,
    CAS_UNRESOLVED,
    CAS_UNREACHABLE,
    NAME_UNREACHABLE,
    CHEBI_UNREACHABLE,
    CHEBI_UNRESOLVED,
    CHEBI_NO_STRUCTURE,
)

# What happened when one identifier was resolved, independent of any comparison.
# check_file's run-level accounting reads these rather than re-deriving intent from
# the sheet cell: the cell says what was asked for, only the resolver knows what
# actually happened, and every previous version of the run verdict was wrong
# because it measured one against the other.
IDENTIFIER_NOT_CHECKED = "not_checked"
IDENTIFIER_RESOLVED = "resolved"
IDENTIFIER_MISSING = "missing"
IDENTIFIER_UNREACHABLE = "unreachable"

# Which requested check went unasked because an upstream could not be reached.
# Kept out of `review`, whose four values are a verdict on the chemistry; this is
# a verdict on the run, and it belongs in its own column so a re-run after the
# outage clears can be diffed against this one.
UNASKED_NONE = ""
UNASKED_ID = "id"
UNASKED_NAME = "name"
UNASKED_BOTH = "both"

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
    "unasked",
    "id_cas_verdict",
    "name_cas_verdict",
    "cas_inchikey",
    "cas_pubchem_name",
    "chebi_inchikey",
    "name_query",
    "name_inchikey",
]

# Per-identifier outcomes. Returned by check_row for the run-level accounting but
# deliberately not written to the CSV: the verdict columns already say what
# happened to each row, and three more columns would be restating them.
STATUS_FIELDS = ["cas_status", "chebi_status", "name_status"]


def empty_result() -> dict[str, Any]:
    """A result row with both verdicts unchecked and nothing resolved."""
    result: dict[str, Any] = {field: "" for field in OUTPUT_FIELDS_APPENDED}
    result["id_cas_verdict"] = NOT_CHECKED
    result["name_cas_verdict"] = NOT_CHECKED
    result["review"] = REVIEW_UNVERIFIED
    result["review_rank"] = REVIEW_RANK[REVIEW_UNVERIFIED]
    # `unasked` is already "" from the comprehension above, which is UNASKED_NONE.
    for field in STATUS_FIELDS:
        result[field] = IDENTIFIER_NOT_CHECKED
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

    Deliberately stays on the lenient get_with_retry rather than raising, and its
    failures never reach the run-level outage verdict. This is the one path where
    an outage cannot produce a false clean sheet — it can only leave a row flagged
    for a human — and counting it would exit 1 on a run whose only content is a
    genuine finding, since it is by far the highest-volume caller here. Failures
    are logged so a flagged row is still traceable to the network.
    """
    if not inchikey:
        return ""
    if inchikey in _parent_cache:
        return _parent_cache[inchikey]

    unreachable = False

    def _step(url: str) -> Any:
        """One hop of the chain, noting an outage without letting it escape."""
        nonlocal unreachable
        resp, outcome = get_with_retry_status(url)
        if outcome == OUTCOME_UNREACHABLE:
            unreachable = True
            log.warning(
                "Parent lookup unreachable for %s; leaving the difference as found.",
                inchikey,
            )
        time.sleep(REQUEST_DELAY)
        return resp

    result = ""
    key = urllib.parse.quote(inchikey, safe="")
    cid = _first_cid(_step(f"{BASE}/compound/inchikey/{key}/cids/JSON"))
    if cid is not None:
        parent_cid = _first_cid(
            _step(f"{BASE}/compound/cid/{cid}/cids/JSON?cids_type=parent")
        )
        if parent_cid is not None:
            resp = _step(f"{BASE}/compound/cid/{parent_cid}/property/InChIKey/JSON")
            if resp is not None:
                try:
                    result = resp.json()["PropertyTable"]["Properties"][0]["InChIKey"]
                except (ValueError, KeyError, IndexError, TypeError):
                    result = ""

    # Answers are cached, silences are not. "PubChem has no parent for this
    # structure" is an answer and will not change during the run, so caching the
    # empty string saves re-walking three requests for every repeat. A transient
    # failure must still not disable the salt demotion for the rest of the run.
    if result or not unreachable:
        _parent_cache[inchikey] = result
    return result


def _first_cid(resp: Any) -> int | None:
    """
    First CID in a PUG list response, or None for anything unusable.

    The subscript is inside the guard: a 200 whose CID field is an int or a dict
    rather than a list would otherwise raise out through parent_inchikey and
    refine_skeleton_difference into check_row, which promises never to raise and
    would take the whole sheet down with it.
    """
    if resp is None:
        return None
    try:
        cids = resp.json()["IdentifierList"]["CID"]
        if not isinstance(cids, list):
            # A string payload would otherwise subscript to its first character
            # — "5291" quietly becoming "5" — the one malformed shape the int
            # and dict guards let through.
            raise TypeError("CID is not a list")
        return cids[0] if cids else None
    except (ValueError, KeyError, TypeError, IndexError):
        return None


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
    """
    Every InChIKey PubChem resolves a name to. Empty when it resolves none.

    Raises PubChemUnavailableError when PubChem could not be asked, so "this name
    is not a compound" stays distinguishable from "we never found out".
    """
    query = urllib.parse.quote(name, safe="")
    resp, outcome = get_with_retry_status(
        f"{BASE}/compound/name/{query}/property/InChIKey/JSON"
    )
    time.sleep(REQUEST_DELAY)
    if outcome == OUTCOME_UNREACHABLE:
        raise PubChemUnavailableError(f"PubChem unreachable resolving name {name!r}")
    if resp is None:
        return []
    keys: list[str] = []
    try:
        # The iteration belongs inside the guard too: `Properties: null` or a bare
        # number is as malformed as a missing key, and a dict would iterate its
        # own keys and quietly yield nothing — reported as name_unresolved, which
        # would be a claim about the compound rather than about the payload.
        properties = resp.json()["PropertyTable"]["Properties"]
        if not isinstance(properties, list):
            raise TypeError("Properties is not a list")
        for entry in properties:
            key = entry.get("InChIKey") if isinstance(entry, dict) else None
            if key and key not in keys:
                keys.append(key)
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        # A 200 that is not this endpoint's JSON is something other than this API
        # answering — a maintenance page or a moved endpoint. Not evidence that
        # the name is unknown.
        raise PubChemUnavailableError(
            f"PubChem returned an unparseable payload for name {name!r}"
        ) from exc
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

    PubChemUnavailableError from any candidate propagates rather than being caught,
    which deliberately suppresses the token fallback. The fallback's premise is
    that no whole-string form resolves; an unreachable whole-string query never
    established that, and falling through would compare the row against a free
    base the sheet never named while the audit trail read like a normal fallback.
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
    # `or {}` only rescues a falsy default_structure; a list or a string is truthy
    # and has no .get, and check_payload_shape does not validate this field. An
    # AttributeError here would escape check_row's never-raises contract.
    structure = payload.get("default_structure")
    key = (
        structure.get("standard_inchi_key") or "" if isinstance(structure, dict) else ""
    )
    return (key, "") if key else ("", CHEBI_NO_STRUCTURE)


def cas_structure(cas: Any) -> tuple[str, str]:
    """
    Resolve a CAS number to (inchikey, pubchem_preferred_name).

    Two requests, not chebi_lookup.lookup_cas()'s four: only the InChIKey and Title
    are wanted, and fetching xrefs and synonyms for every row spent about a minute
    of REQUEST_DELAY on a 117-row sheet producing data that was discarded.

    Raises PubChemUnavailableError when *either* request could not be answered, so
    "PubChem has no such CAS" stays distinguishable from "PubChem was unreachable".
    Both requests matter: a property call that fails after the CID resolved is
    still a CAS this run never saw the structure of.
    """
    text = str(cas or "").strip()
    if not text:
        return "", ""
    # Not quoted here: cas_to_cid_status quotes the value itself, so that its
    # other caller (chebi_lookup.lookup_cas) is covered by the same guard rather
    # than by each call site remembering. Quoting twice would encode the '%'.
    cid, outcome = cas_to_cid_status(text)
    if outcome == OUTCOME_UNREACHABLE:
        raise PubChemUnavailableError(f"PubChem unreachable resolving CAS {text!r}")
    if cid is None:
        return "", ""
    resp, outcome = get_with_retry_status(
        f"{BASE}/compound/cid/{cid}/property/InChIKey,Title/JSON"
    )
    time.sleep(REQUEST_DELAY)
    if outcome == OUTCOME_UNREACHABLE:
        raise PubChemUnavailableError(
            f"PubChem unreachable fetching properties for CAS {text!r}"
        )
    if resp is None:
        return "", ""
    try:
        # The .get() calls stay inside the guard: a Properties list holding
        # strings rather than objects would otherwise raise AttributeError past
        # it, through check_row — which promises never to raise — and abort the
        # sheet at whatever row the bad payload landed on.
        properties = resp.json()["PropertyTable"]["Properties"][0]
        return properties.get("InChIKey") or "", properties.get("Title") or ""
    except (ValueError, KeyError, IndexError, TypeError, AttributeError) as exc:
        raise PubChemUnavailableError(
            f"PubChem returned an unparseable payload for CAS {text!r}"
        ) from exc


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
    skip: Collection[str] = (),
) -> dict[str, Any]:
    """
    Cross-check one row's identifiers by structure.

    The CAS number is the pivot for both comparisons, so without it neither can be
    made. Never raises: anything unresolvable is reported as such rather than as a
    finding about the compound.

    `skip` names sides ("cas", "chebi", "name") whose upstream a caller already
    knows to be down, so the request is not made at all and the side is reported
    unreachable directly. That is what it is — the row's question goes unanswered
    either way — and it lets a batch run finish the columns that still work
    instead of paying full retry backoff per row for a settled answer.
    """
    result = empty_result()

    wants_id_check = bool(str(chebi_id or "").strip())
    wants_name_check = bool(str(name or "").strip())
    if not wants_id_check and not wants_name_check:
        # Nothing to compare the CAS against, so resolving it would cost requests
        # for a result that is discarded either way. cas_status stays not_checked:
        # a CAS this run never asked about is not a CAS that failed.
        return result

    has_cas = bool(str(cas or "").strip())
    cas_unreachable = False
    cas_key, cas_name = "", ""
    if has_cas:
        if "cas" in skip:
            cas_unreachable = True
        else:
            try:
                cas_key, cas_name = cas_structure(cas)
            except PubChemUnavailableError as exc:
                log.error("%s", exc)
                cas_unreachable = True
        result["cas_status"] = (
            IDENTIFIER_UNREACHABLE
            if cas_unreachable
            else IDENTIFIER_RESOLVED
            if cas_key
            else IDENTIFIER_MISSING
        )
    result["cas_inchikey"] = cas_key
    result["cas_pubchem_name"] = cas_name
    # A blank cell was never asked about, and an unreachable one was asked without
    # an answer. Only a non-blank CAS that PubChem answered "no" about is a failed
    # check. cas_structure() cannot draw these lines itself, so the caller does.
    if cas_unreachable:
        cas_verdict_if_missing = CAS_UNREACHABLE
    elif has_cas:
        cas_verdict_if_missing = CAS_UNRESOLVED
    else:
        cas_verdict_if_missing = NOT_CHECKED

    if not cas_key:
        # The pivot is gone, so both verdicts are already decided and resolving
        # the other identifiers would spend 2-3 requests per row on structures
        # nothing can be compared against. This is the mirror of the early return
        # above, and it matters most on the sheet where --cas-column is pointed at
        # the wrong header: every row would otherwise pay full price to produce
        # `unverified`. Their statuses stay not_checked, which is accurate — this
        # run learned nothing about those columns.
        if wants_id_check:
            result["id_cas_verdict"] = cas_verdict_if_missing
        if wants_name_check:
            result["name_cas_verdict"] = cas_verdict_if_missing
        return _finish(result)

    if wants_id_check:
        if "chebi" in skip:
            chebi_key, problem = "", CHEBI_UNREACHABLE
        else:
            chebi_key, problem = chebi_structure(chebi_id)
        result["chebi_inchikey"] = chebi_key
        result["chebi_status"] = (
            IDENTIFIER_UNREACHABLE
            if problem == CHEBI_UNREACHABLE
            else IDENTIFIER_MISSING
            if problem == CHEBI_UNRESOLVED
            # chebi_no_structure means ChEBI held the record and it legitimately
            # carries no structure, as class terms and R-group entries do. That is
            # a real ChEBI ID, so it must not count as evidence of a wrong column.
            else IDENTIFIER_RESOLVED
        )
        if problem:
            result["id_cas_verdict"] = problem
        else:
            verdict = compare_structures(cas_key, [chebi_key])
            if verdict == SKELETON_DIFFERS:
                verdict = refine_skeleton_difference(cas_key, [chebi_key])
            result["id_cas_verdict"] = verdict

    if wants_name_check:
        if "name" in skip:
            result["name_status"] = IDENTIFIER_UNREACHABLE
            result["name_cas_verdict"] = NAME_UNREACHABLE
            return _finish(result)
        try:
            query, name_keys, total_found = name_structure(name)
        except PubChemUnavailableError as exc:
            log.error("%s", exc)
            result["name_status"] = IDENTIFIER_UNREACHABLE
            result["name_cas_verdict"] = NAME_UNREACHABLE
            return _finish(result)
        result["name_status"] = IDENTIFIER_RESOLVED if name_keys else IDENTIFIER_MISSING
        truncated = total_found > len(name_keys)
        if truncated:
            query = f"{query} (truncated: {total_found} structures)"
        result["name_query"] = query
        result["name_inchikey"] = " | ".join(name_keys)
        if not name_keys:
            result["name_cas_verdict"] = NAME_UNRESOLVED
        else:
            verdict = compare_structures(cas_key, name_keys)
            # A truncated row that lands here is heading for name_ambiguous
            # below regardless of what refinement finds, so skip a refinement
            # that costs up to 3 requests per candidate for a result that is
            # about to be discarded.
            if verdict == SKELETON_DIFFERS and not truncated:
                verdict = refine_skeleton_difference(cas_key, name_keys)
            # PubChem returns structures in CID order, not relevance order, so the
            # true match may be among the ones we did not compare. Reporting a
            # difference here would be a finding of our own making.
            if truncated and verdict != MATCH:
                verdict = NAME_AMBIGUOUS
            result["name_cas_verdict"] = verdict

    return _finish(result)


def _finish(result: dict[str, Any]) -> dict[str, Any]:
    """Fill in the derived columns: review, review_rank, unasked."""
    result["review"] = review_level(
        result["id_cas_verdict"], result["name_cas_verdict"]
    )
    result["review_rank"] = REVIEW_RANK[result["review"]]

    # `review` stays a verdict on the chemistry — `ok` still means every comparison
    # that was made agreed, which is true even when the other side was unreachable.
    # Which question went unasked is a different fact about a different subject, so
    # it gets its own column rather than a fifth review level: overloading the sort
    # key would bury a real skeleton_differs finding behind a transient 503.
    id_unasked = result["id_cas_verdict"] in OUTAGE_VERDICTS
    name_unasked = result["name_cas_verdict"] in OUTAGE_VERDICTS
    if id_unasked and name_unasked:
        result["unasked"] = UNASKED_BOTH
    elif id_unasked:
        result["unasked"] = UNASKED_ID
    elif name_unasked:
        result["unasked"] = UNASKED_NAME
    return result
