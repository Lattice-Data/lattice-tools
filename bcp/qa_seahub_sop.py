"""
SOP validation for SeaHub lab raw uploads (``raw_assay='seahub_sci'``).

Validates the full S3 location and filename structure against the SOP::

    s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/{sublibrary}/{wafer}/
        {wafer}-{sublibrary}[_{well}]_{sublibrary type}-{UG}-{barcode}.trim.*

Known-good examples (all three must produce zero violations)::

    czi-hamazaki  hamazaki-seahub-bcp/CHEM3-R100/raw/R100E/441389/
        441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram
    czi-trapnell  trapnell-seahub-bcp/REF3/raw/REF3_P05_2/436830/
        436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram
    czi-trapnell  trapnell-seahub-bcp/CHEM16/raw/P03/432640/
        432640-CHEM16_P03_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram

Design notes
------------
* Every rule is evaluated independently, and the token-level rules run even
  when the stem does not match :data:`SEAHUB_STEM_RE`.  A wholly non-compliant
  name therefore reports each of its defects by name instead of collapsing into
  a single ``unparseable_stem``.
* There is deliberately no "ExperimentID must not appear in the filename" rule.
  The ExperimentID legitimately appears inside sublibrary names such as
  ``REF3_P05_2``, and is absent entirely from the hamazaki example, so
  sublibrary agreement is the only rule governing that token.
* The sublibrary folder may be spelled either ``{sublibrary}`` or with the
  redundant ``{ExperimentID}_`` elided, and neither is a defect: the third
  example above is the elided form.  The ExperimentID is already an ancestor
  path segment, so the prefix carries nothing the path does not, and every real
  trimmed upload measured -- REF3, GENE7, CHEM16 -- elides it on every
  sublibrary; not one mixes the two spellings.  Nothing downstream depends on
  which is used: the filename is the authoritative sublibrary name, and the
  cross-bucket identity is ``(wafer, UG)``.  Demanding the full form reported
  every GENE7 sublibrary as broken and proposed a move for all 5184 of its
  objects, on an upload that is fine.  ``sublibrary_mismatch`` is reserved for a
  filename neither spelling explains.
* A repeated leading wafer token is normalized away *before* any token rule
  runs, and the rules then see only the normalized stem.  On the raw stem the
  type-less pattern swallows the second wafer into the group
  (``437120-REF3_P04_1_A1``), which no folder can be reconciled against;
  validating both forms would instead double every row.  ``repeated_token``
  stays basename-scoped and is now the fallback for repeats that are not the
  wafer -- comparing filename tokens against folder segments would flag the
  valid ``REF3_P05_2_A10`` under ExperimentID folder ``REF3`` as a duplicate.
* ``expected_name`` repairs only the defect its own rule describes, computed on
  the normalized stem.  Composing several into a single corrected key is
  :mod:`qa_seahub_rename`'s job, because that needs the vendor index to resolve
  a missing sublibrary type.
* ``scope`` controls how far a fact reaches, and therefore how it dedupes:
  ``object`` stays per object, ``stem`` collapses per well,
  ``suffix`` collapses per distinct unrecognised extension, and ``upload``
  collapses to one row for the whole listing, since the bucket and project are
  the same fact however many wells sit under them.  The authoritative list is
  :data:`qa_constants.SEAHUB_VIOLATION_SCOPES`; no count is written down here,
  because the last one went stale the moment ``suffix`` was added.
* There is no ``bad_ug`` or ``bad_barcode`` rule.  Both stem patterns pin
  ``(?P<ug>Z\\d{4})`` and ``(?P<barcode>[ACGT]+)$``, so reaching such a rule
  would already guarantee it passes; a malformed token fails both patterns and
  is reported as ``unparseable_stem`` instead.
* Violations are non-fatal: QA continues and reports them as a table.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping

from qa_constants import (
    SEAHUB_BARE_TO_TRIM_SUFFIX,
    SEAHUB_DOUBLED_WAFER_RE,
    SEAHUB_DOWNLOAD_DUP_SUFFIX_RE,
    SEAHUB_NON_SEQ_BASENAMES,
    SEAHUB_NON_SEQ_EXTENSIONS,
    SEAHUB_NON_SEQ_NAME_RES,
    SEAHUB_STEM_NO_TYPE_RE,
    SEAHUB_STEM_RE,
    SEAHUB_SUBLIBRARY_TYPES,
    SEAHUB_WELL_RE,
    raw_expected,
)
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family

__all__ = [
    "SeahubStemGroup",
    "SopViolation",
    "expected_sop_name",
    "group_seahub_keys",
    "is_non_sequencing_artifact",
    "normalize_doubled_wafer",
    "seahub_group_parts",
    "sop_violation_summary",
    "unrecognized_suffix",
    "validate_seahub_group",
    "validate_seahub_key",
    "validate_seahub_stems",
]


@dataclass(frozen=True)
class SopViolation:
    """One SOP rule broken by one object (or stem).

    ``expected_name`` is the corrected *basename* repairing only this rule's
    defect.  It is deliberately not the fully-corrected name: a well with four
    defects yields four rows, each naming its own repair, and composing them into
    one target key is :mod:`qa_seahub_rename`'s job.  So a doubled-wafer bare
    file gets ``437120-437120-....trim.cram`` from ``missing_trim_infix``, which
    repairs the infix and nothing else -- the doubled token is
    ``duplicated_wafer_token``'s row to repair, and it carries its own
    ``expected_name``.  Rules that run after normalization see the normalized
    stem; ``missing_trim_infix`` is emitted before it and works on the basename
    as delivered.

    ``scope`` says how widely the fact applies and drives dedup when reporting:
    ``object`` for a rule about one S3 object, ``stem`` for one about a well,
    ``suffix`` for one about a distinct unrecognised extension, and ``upload``
    for one about the bucket or project, which is the same fact however many
    wells sit under it.  A suffix or upload defect reported per object would bury
    every other finding beneath it.

    There is no ``expected_folder``, and no ``folder`` scope.  Both existed for
    ``sublibrary_folder_truncated`` alone, and nothing else ever set or read
    them; once the elided-prefix folder became an accepted spelling they had no
    producer left.
    """

    type: str
    s3_path: str
    detail: str
    expected_name: str = ""
    scope: str = "stem"

    def as_dict(self) -> dict[str, str]:
        return asdict(self)


def expected_sop_name(basename: str, suffix: str, family: str) -> str:
    """Return the SOP-suffixed filename for a bare-family artifact."""
    if family != "bare":
        return basename
    trim_suffix = SEAHUB_BARE_TO_TRIM_SUFFIX.get(suffix)
    if trim_suffix is None:
        return basename
    return f"{basename[: -len(suffix)]}{trim_suffix}"


def _split_tokens(stem: str) -> list[str]:
    """Split a stem on hyphens and underscores, dropping empty tokens."""
    tokens: list[str] = []
    for hyphen_part in stem.split("-"):
        tokens.extend(part for part in hyphen_part.split("_") if part)
    return tokens


def _check_repeated_tokens(group: str, stem: str, s3_path: str) -> list[SopViolation]:
    """Catch a token repeated inside the sublibrary/well portion of a stem.

    Scoped to the ``group`` slot, and run only once the stem has matched, so the
    bag holds nothing but tokens the SOP puts in one place. Splitting the whole
    stem on both ``-`` and ``_`` before the match, as this did, put the type
    token in the same bag as the sublibrary name: a sublibrary legitimately named
    ``REF3_GEX_P01`` under type ``GEX_hash_oligo`` contributed ``GEX`` twice and
    reported ``repeated_token`` as its only violation on an otherwise clean name.
    Via the rename gate that withheld the proposal and flipped the well to
    ``UNKNOWN`` -- for a whole sublibrary at once, 288 to 576 objects, with no
    filename the tool would have accepted.

    The rule has never fired truthfully: zero across 81,936 real objects. It is
    kept because a repeat *within* the group (``REF3_P01_P01``) is a real naming
    error and nothing else reports it; the wafer case belongs to
    ``duplicated_wafer_token``.
    """
    seen: dict[str, int] = {}
    for token in _split_tokens(group):
        seen[token] = seen.get(token, 0) + 1
    repeated = sorted(token for token, count in seen.items() if count > 1)
    if not repeated:
        return []
    return [
        SopViolation(
            type="repeated_token",
            s3_path=s3_path,
            detail=(
                f"token(s) {', '.join(repeated)} appear more than once in the "
                f"sublibrary portion {group!r} of {stem!r}"
            ),
        )
    ]


def normalize_doubled_wafer(stem: str) -> tuple[str, bool]:
    """Drop a repeated leading wafer token, returning ``(stem, was_doubled)``.

    Only an *identical* repeat is collapsed.  ``437120-438514-...`` names two
    different wafers, which is a defect QA must not guess a repair for.
    """
    match = SEAHUB_DOUBLED_WAFER_RE.match(stem)
    if match is None or match.group("first") != match.group("second"):
        return stem, False
    return stem[len(match.group("first")) + 1 :], True


def is_non_sequencing_artifact(basename: str) -> bool:
    """True for browser/bulk-download leftovers rather than pipeline output.

    Only consulted once :func:`qa_mods.seahub_stem_and_family` has already
    declined the name, so a real artifact can never reach here.  A single
    download counter is stripped as an *additional* candidate, so both
    ``ug-icon.png`` and ``ug-icon.png.1`` classify while ``x.trim.cram.1``
    stays an unexpected suffix.
    """
    lowered = basename.lower()
    candidates = [lowered]
    stripped = SEAHUB_DOWNLOAD_DUP_SUFFIX_RE.sub("", lowered)
    if stripped != lowered:
        candidates.append(stripped)
    for candidate in candidates:
        if candidate in SEAHUB_NON_SEQ_BASENAMES:
            return True
        if any(pattern.match(candidate) for pattern in SEAHUB_NON_SEQ_NAME_RES):
            return True
        if any(candidate.endswith(ext) for ext in SEAHUB_NON_SEQ_EXTENSIONS):
            return True
    return False


def seahub_group_parts(
    group: str, sublibrary: str, experiment_id: str
) -> tuple[str, str, bool]:
    """Reconcile a filename's sublibrary token against its folder.

    Returns ``(sublibrary, trailing_token, matched)``, where ``sublibrary`` is
    the name *the filename* gives -- the authoritative one -- and ``matched`` is
    False only when the folder cannot explain the filename at all.

    Two folder spellings are accepted, and neither is a defect: the full
    ``{ExperimentID}_{sublibrary}`` (``REF3_P05_2/…REF3_P05_2_A10…``) and the
    same name with the redundant ``{ExperimentID}_`` elided
    (``P03/…CHEM16_P03_A1…``).  The ExperimentID is already an ancestor path
    segment, so the prefix carries nothing the path does not, and every real
    trimmed upload measured -- REF3, GENE7, CHEM16 -- elides it on every
    sublibrary.  A filename carrying no ExperimentID at all (``R100E/…R100E…``)
    matches the first form directly.

    Candidate order decides the answer only when *both* candidates can explain
    the group, which needs the folder name to be a leading token-prefix of its
    own ExperimentID -- folder ``P10`` under ExperimentID ``P10_2`` makes the full
    form ``P10_2_P10``, which itself starts with ``P10_``.  Trying the folder as
    delivered first prefers the literal reading of the path.  No real upload has
    that shape, and every other input reaches the same result either way, so this
    order is a tie-break rather than a rule; it is pinned by
    ``test_candidate_order_only_matters_when_both_can_explain_the_group``.
    """
    candidates = [sublibrary]
    if experiment_id:
        candidates.append(f"{experiment_id}_{sublibrary}")

    for candidate in candidates:
        if group == candidate:
            return candidate, "", True
        prefix = f"{candidate}_"
        if group.startswith(prefix):
            return candidate, group[len(prefix) :], True
    return sublibrary, "", False


def _check_group_token(
    group: str,
    sublibrary: str,
    experiment_id: str,
    s3_path: str,
    trailing_claimed: bool = False,
) -> list[SopViolation]:
    """Report a folder the filename's sublibrary token cannot be squared with.

    Only the genuine disagreement is a defect.  A folder that elides the
    redundant ``{ExperimentID}_`` prefix is one of two accepted spellings --
    see :func:`seahub_group_parts` -- so ``sublibrary_mismatch`` means what it
    says: neither spelling explains the filename.

    ``trailing_claimed`` says the caller has already reported the trailing token
    under another rule.  The relaxed-stem branch does exactly that: a trailing
    token that is not a well is reported as ``invalid_sublibrary_type``, and
    ``bad_well`` would then describe the same token again with contradictory
    advice -- one row saying it should be a type, the next saying it should be a
    well.  An upload that misspells its type token on every well doubled its SOP
    table that way, against a module built on one row per distinct fact.
    """
    name_sublibrary, trailing, matched = seahub_group_parts(
        group, sublibrary, experiment_id
    )

    if not matched:
        accepted = f"{sublibrary!r}"
        if experiment_id:
            accepted += f" nor {f'{experiment_id}_{sublibrary}'!r}"
        return [
            SopViolation(
                type="sublibrary_mismatch",
                s3_path=s3_path,
                detail=(
                    f"filename sublibrary {group!r} is neither {accepted}, with or "
                    "without a trailing _<well>; rename either the folder or the "
                    "files so they agree"
                ),
            )
        ]

    if trailing and not trailing_claimed and not SEAHUB_WELL_RE.match(trailing):
        return [
            SopViolation(
                type="bad_well",
                s3_path=s3_path,
                detail=(
                    f"token {trailing!r} after sublibrary {name_sublibrary!r} is "
                    "not a well of the form [A-H]<1-2 digits>"
                ),
            )
        ]
    return []


def _check_path(bucket: str, s3_key: str) -> tuple[list[SopViolation], dict | None]:
    """Validate bucket/project/depth, returning violations and folder info."""
    s3_path = f"s3://{bucket}/{s3_key}"
    violations: list[SopViolation] = []

    lab = ""
    if not bucket.startswith("czi-") or len(bucket) <= len("czi-"):
        violations.append(
            SopViolation(
                type="bad_bucket",
                s3_path=s3_path,
                detail=f"bucket {bucket!r} is not of the form czi-<lab>",
                scope="upload",
            )
        )
    else:
        lab = bucket[len("czi-") :]

    project = s3_key.split("/")[0] if "/" in s3_key else ""
    # An exact prefix, not the first hyphen-separated token: a lab whose name
    # contains a hyphen ("van-der-berg") splits to "van" and fails against its
    # own project. The rule is upload-scope, so the consequence is one row for
    # the whole listing plus the rename cell's "these destinations sit somewhere
    # the SOP rejects" banner -- on an upload that is correct. The trailing "-"
    # is what stops lab "lab" from accepting project "labalpha-seahub-bcp".
    if lab and project and not (project == lab or project.startswith(f"{lab}-")):
        violations.append(
            SopViolation(
                type="lab_project_mismatch",
                s3_path=s3_path,
                detail=(
                    f"project {project!r} must start with the lab name "
                    f"{lab!r} from bucket {bucket!r}"
                ),
                scope="upload",
            )
        )

    path_info = parse_seahub_raw_path(s3_key)
    if path_info is None:
        violations.append(
            SopViolation(
                type="bad_path_depth",
                s3_path=s3_path,
                detail=(
                    "key is not {project}/{ExperimentID}/raw/{sublibrary}/"
                    "{wafer}/{filename}"
                ),
            )
        )
    return violations, path_info


def validate_seahub_key(
    bucket: str,
    s3_key: str,
    *,
    assay_by_identity: Mapping[tuple[str, str], str] | None = None,
) -> list[SopViolation]:
    """Validate one SeaHub object against the SOP.

    Returns every rule broken by this object; an empty list means compliant.

    ``assay_by_identity`` maps ``(wafer, UG)`` to the sublibrary type read from
    the untrimmed vendor delivery.  Supplied, it lets ``invalid_sublibrary_type``
    name the corrected filename; omitted, that rule still fires but cannot say
    what the missing token should have been.  A plain mapping keeps this module a
    leaf that imports nothing but constants and parsers.
    """
    s3_path = f"s3://{bucket}/{s3_key}"
    violations, path_info = _check_path(bucket, s3_key)
    if path_info is None:
        return violations

    basename = s3_key.split("/")[-1]
    parsed = seahub_stem_and_family(basename)
    if parsed is None:
        if is_non_sequencing_artifact(basename):
            violations.append(
                SopViolation(
                    type="non_sequencing_artifact",
                    s3_path=s3_path,
                    detail=(
                        f"{basename!r} is likely a file unrelated to the CRO "
                        "sequencing — it matches bulk-download tool output "
                        "(browser pages, icons, manifest listings) rather than "
                        "any SOP artifact, and carries no sequencing data QA "
                        "can validate"
                    ),
                    scope="object",
                )
            )
            return violations
        violations.append(
            SopViolation(
                type="unexpected_suffix",
                s3_path=s3_path,
                detail=(
                    f"{basename!r} does not end in a SeaHub artifact suffix; "
                    "expected one of the six .trim.* artifacts"
                ),
                scope="object",
            )
        )
        return violations

    raw_stem, suffix, family = parsed
    if family == "bare":
        violations.append(
            SopViolation(
                type="missing_trim_infix",
                s3_path=s3_path,
                detail=(
                    f"artifact {suffix!r} is missing the SOP '.trim' infix, so "
                    "the upload is indistinguishable from untrimmed data"
                ),
                expected_name=expected_sop_name(basename, suffix, family),
            )
        )

    # Normalize the doubled wafer BEFORE any token rule runs, and then run them
    # on the normalized stem only.  On the raw stem the type-less pattern
    # swallows the second wafer into the group ('437120-REF3_P04_1_A1'), which no
    # folder can be reconciled against; validating both stems would instead
    # double every row for the ~1439 affected objects.
    stem, was_doubled = normalize_doubled_wafer(raw_stem)
    if was_doubled:
        violations.append(
            SopViolation(
                type="duplicated_wafer_token",
                s3_path=s3_path,
                detail=(
                    f"wafer {stem.split('-')[0]!r} appears twice at the head of "
                    f"{raw_stem!r}; the SOP name carries it once"
                ),
                expected_name=f"{stem}{suffix}",
            )
        )

    # False on the strict path: a stem that matched SEAHUB_STEM_RE has its type
    # token accounted for, so the trailing token is genuinely a well or nothing.
    type_token_claimed = False
    match = SEAHUB_STEM_RE.match(stem)
    if match is None:
        relaxed = SEAHUB_STEM_NO_TYPE_RE.match(stem)
        if relaxed is None:
            violations.append(
                SopViolation(
                    type="unparseable_stem",
                    s3_path=s3_path,
                    detail=(
                        f"{stem!r} does not decompose into {{wafer}}-"
                        "{sublibrary}[_{well}]_{type}-{UG}-{barcode}"
                    ),
                )
            )
            return violations
        match = relaxed
        vendor_assay = (assay_by_identity or {}).get(
            (match.group("wafer"), match.group("ug"))
        )
        # The relaxed pattern matches whether the type token is absent or merely
        # unrecognised, and those are different facts. A trailing token that is
        # not a well is a type the SOP does not know, so saying the stem "carries
        # no sublibrary type" describes the wrong defect -- and appending the
        # vendor type to it would propose a name with the unrecognised token
        # still buried in the sublibrary.
        _sublibrary, trailing, _matched = seahub_group_parts(
            match.group("group"), path_info["sublibrary"], path_info["experiment_id"]
        )
        unrecognized = bool(trailing) and not SEAHUB_WELL_RE.match(trailing)
        type_token_claimed = unrecognized
        if unrecognized:
            corrected = ""
            detail = (
                f"{stem!r} carries {trailing!r} where the sublibrary type "
                f"belongs; expected one of {', '.join(SEAHUB_SUBLIBRARY_TYPES)}"
            )
            if vendor_assay:
                detail += (
                    f" (the untrimmed vendor delivery gives {vendor_assay!r}, but "
                    "the name cannot be corrected without knowing whether "
                    f"{trailing!r} belongs to the sublibrary)"
                )
        elif vendor_assay:
            corrected = (
                f"{match.group('wafer')}-{match.group('group')}_{vendor_assay}"
                f"-{match.group('ug')}-{match.group('barcode')}{suffix}"
            )
            detail = (
                f"{stem!r} carries no sublibrary type; the untrimmed vendor "
                f"delivery gives {vendor_assay!r} for this well"
            )
        else:
            corrected = ""
            detail = (
                f"{stem!r} carries no sublibrary type; expected one of "
                f"{', '.join(SEAHUB_SUBLIBRARY_TYPES)}"
            )
        violations.append(
            SopViolation(
                type="invalid_sublibrary_type",
                s3_path=s3_path,
                detail=detail,
                expected_name=corrected,
            )
        )

    violations.extend(_check_repeated_tokens(match.group("group"), stem, s3_path))

    if match.group("wafer") != path_info["wafer"]:
        violations.append(
            SopViolation(
                type="wafer_mismatch",
                s3_path=s3_path,
                detail=(
                    f"filename wafer {match.group('wafer')!r} differs from "
                    f"folder wafer {path_info['wafer']!r}"
                ),
            )
        )

    violations.extend(
        _check_group_token(
            match.group("group"),
            path_info["sublibrary"],
            path_info["experiment_id"],
            s3_path,
            trailing_claimed=type_token_claimed,
        )
    )

    # No bad_ug / bad_barcode rule: both stem patterns pin (?P<ug>Z\d{4}) and
    # (?P<barcode>[ACGT]+)$, so a match guarantees both tokens are well formed
    # and the rules could never fire. A malformed one fails both patterns and is
    # reported as unparseable_stem, whose detail names the expected shape.
    return violations


@dataclass(frozen=True)
class SeahubStemGroup:
    """The artifacts of one well, as delivered.

    ``identity`` is ``(wafer, UG)`` when the normalized stem parses -- the same
    key :mod:`qa_seahub_source` indexes both sides of the cross-bucket comparison
    on, so a group joins the vendor index directly.  Otherwise it falls back to
    ``(raw_dir, stem)``, which at least keeps the artifacts of one unparseable
    well together.

    Carries what a caller or its tests read. ``experiment_id``, ``sublibrary``,
    ``wafer_folder``, ``stems``, ``wafer`` and ``ug`` were populated here and
    read nowhere -- not by a caller and not by a test -- and each is already
    recoverable from ``keys`` via :func:`parse_seahub_raw_path` or from
    ``normalized_stem``.
    """

    bucket: str
    raw_dir: str
    identity: tuple[str, str]
    keys: tuple[str, ...]
    suffixes: tuple[str, ...]
    families: frozenset[str]
    normalized_stem: str
    well_id: str

    @property
    def has_cram(self) -> bool:
        return any(s in (".cram", ".trim.cram") for s in self.suffixes)


def group_seahub_keys(
    bucket: str, s3_keys: list[str]
) -> tuple[list[SeahubStemGroup], list[str]]:
    """Bucket a listing into per-well groups, returning ``(groups, unparsed)``.

    ``unparsed`` holds keys that carry no recognizable artifact suffix -- junk and
    genuinely unexpected names -- which are reported per object rather than per
    well.
    """
    unparsed: list[str] = []
    buckets: dict[tuple[str, tuple[str, str]], list[tuple[str, str, str, str]]] = {}

    for s3_key in sorted(s3_keys):
        parsed = seahub_stem_and_family(s3_key.split("/")[-1])
        if parsed is None:
            unparsed.append(s3_key)
            continue
        raw_stem, suffix, family = parsed
        raw_dir = "/".join(s3_key.split("/")[:-1])
        normalized, _doubled = normalize_doubled_wafer(raw_stem)
        match = SEAHUB_STEM_RE.match(normalized) or SEAHUB_STEM_NO_TYPE_RE.match(
            normalized
        )
        identity = (
            (match.group("wafer"), match.group("ug"))
            if match is not None
            else (raw_dir, raw_stem)
        )
        buckets.setdefault((raw_dir, identity), []).append(
            (s3_key, raw_stem, suffix, family)
        )

    groups: list[SeahubStemGroup] = []
    for (raw_dir, identity), members in sorted(buckets.items()):
        path_info = parse_seahub_raw_path(f"{raw_dir}/x") or {}
        first_stem = members[0][1]
        normalized, _doubled = normalize_doubled_wafer(first_stem)
        match = SEAHUB_STEM_RE.match(normalized) or SEAHUB_STEM_NO_TYPE_RE.match(
            normalized
        )
        well_id = ""
        if match is not None:
            _sublibrary, trailing, _matched = seahub_group_parts(
                match.group("group"),
                path_info.get("sublibrary", ""),
                path_info.get("experiment_id", ""),
            )
            well_id = trailing
        groups.append(
            SeahubStemGroup(
                bucket=bucket,
                raw_dir=raw_dir,
                identity=identity,
                keys=tuple(m[0] for m in members),
                suffixes=tuple(m[2] for m in members),
                families=frozenset(m[3] for m in members),
                normalized_stem=normalized,
                well_id=well_id,
            )
        )
    return groups, unparsed


def validate_seahub_group(
    group: SeahubStemGroup,
    *,
    assay_by_identity: Mapping[tuple[str, str], str] | None = None,
) -> list[SopViolation]:
    """Validate every artifact of one well, then collapse duplicate facts.

    Every key is validated, not just the first: a well delivered as a compliant
    ``.trim.cram`` plus a bare ``_fail.csv`` has a real defect that validating
    only the alphabetically-first artifact silently drops.

    Dedup is by ``(type, stem, family)``.  Deduping on ``expected_name`` instead
    would not collapse anything -- three of the rules carry a suffix-dependent
    name -- turning a six-artifact well into nineteen rows.
    """
    violations: list[SopViolation] = []
    seen: set[tuple[str, str, str]] = set()
    for s3_key in group.keys:
        parsed = seahub_stem_and_family(s3_key.split("/")[-1])
        stem, family = (parsed[0], parsed[2]) if parsed is not None else ("", "")
        for violation in validate_seahub_key(
            group.bucket, s3_key, assay_by_identity=assay_by_identity
        ):
            fingerprint = (violation.type, stem, family)
            if fingerprint in seen:
                continue
            seen.add(fingerprint)
            violations.append(violation)
    return violations


def validate_seahub_stems(
    bucket: str,
    s3_keys: list[str],
    *,
    assay_by_identity: Mapping[tuple[str, str], str] | None = None,
) -> list[SopViolation]:
    """Validate a listing, reporting one row per distinct fact.

    The scopes collapse differently: ``object`` rules stay per object, ``stem``
    rules collapse to one row per well, ``suffix`` rules collapse to one row per
    distinct unrecognised extension, and ``upload`` rules collapse to one row for
    the whole listing.  Without that a single wrong bucket would report once per
    well, which on a 288-well upload buries every other finding.
    """
    groups, unparsed = group_seahub_keys(bucket, s3_keys)

    violations: list[SopViolation] = []
    # Deduped across both loops below: the bucket is fixed for the whole listing,
    # so a bucket or project defect is one fact about the upload, not per well.
    # Keyed on the detail so two projects under one bucket still report each.
    seen_upload_facts: set[tuple[str, str]] = set()

    def _upload_fact_is_new(violation: SopViolation) -> bool:
        if violation.scope != "upload":
            return True
        fingerprint = (violation.type, violation.detail)
        if fingerprint in seen_upload_facts:
            return False
        seen_upload_facts.add(fingerprint)
        return True

    for s3_key in unparsed:
        for violation in validate_seahub_key(
            bucket, s3_key, assay_by_identity=assay_by_identity
        ):
            if _upload_fact_is_new(violation):
                violations.append(violation)

    for group in groups:
        for violation in validate_seahub_group(
            group, assay_by_identity=assay_by_identity
        ):
            if _upload_fact_is_new(violation):
                violations.append(violation)
    return _collapse_unknown_suffixes(violations, s3_keys)


def unrecognized_suffix(basename: str) -> str:
    """The trailing extension of a name that matched no artifact suffix.

    Reported so a submitter is told *which* spelling was used rather than being
    handed one row per object. A SeaHub stem ends in ``-{barcode}``, so the
    extension is everything from the first dot after the final hyphen; names
    with no hyphen fall back to the first dot, and names with no dot at all
    report as empty.
    """
    tail = basename.rsplit("-", 1)[-1]
    dot = tail.find(".")
    if dot >= 0:
        return tail[dot:]
    dot = basename.find(".")
    return basename[dot:] if dot >= 0 else ""


def _collapse_unknown_suffixes(
    violations: list[SopViolation], s3_keys: list[str]
) -> list[SopViolation]:
    """One row per distinct unknown suffix, plus a headline if none are known.

    An upload that misspells its artifacts misspells every one of them, so
    ``unexpected_suffix`` per object is one row per object -- 480 identical rows
    on a real upload, which buries the finding it is trying to report. Collapse
    to one row naming the spelling, its count and an example, which is what a
    submitter needs in order to fix it.
    """
    collapsed: list[SopViolation] = []
    by_suffix: dict[str, list[SopViolation]] = {}
    for violation in violations:
        if violation.type != "unexpected_suffix":
            collapsed.append(violation)
            continue
        suffix = unrecognized_suffix(violation.s3_path.rsplit("/", 1)[-1])
        by_suffix.setdefault(suffix, []).append(violation)

    expected = ", ".join(raw_expected["seahub_sci"])
    for suffix, group in sorted(by_suffix.items()):
        shown = suffix or "(no extension)"
        collapsed.append(
            SopViolation(
                type="unexpected_suffix",
                s3_path=group[0].s3_path,
                detail=(
                    f"{len(group)} object(s) end in {shown!r}, which is not a "
                    f"SeaHub artifact suffix; the SOP artifacts are {expected}. "
                    f"Example: {group[0].s3_path}"
                ),
                scope="suffix",
            )
        )

    if s3_keys and not any(
        seahub_stem_and_family(key.rsplit("/", 1)[-1]) for key in s3_keys
    ):
        collapsed.append(
            SopViolation(
                type="no_recognized_artifacts",
                s3_path=f"{len(s3_keys)} object(s)",
                detail=(
                    f"none of the {len(s3_keys)} object(s) in this upload end in a "
                    f"SeaHub artifact suffix ({expected}), so no well could be "
                    "identified and no completeness or trimming check could run"
                ),
                scope="upload",
            )
        )
    return collapsed


def sop_violation_summary(violations: list[SopViolation]) -> dict[str, int]:
    """Count violations by type, for a compact notebook summary."""
    counts: dict[str, int] = {}
    for violation in violations:
        counts[violation.type] = counts.get(violation.type, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
