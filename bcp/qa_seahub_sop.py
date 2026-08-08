"""
SOP validation for SeaHub lab raw uploads (``raw_assay='seahub_sci'``).

Validates the full S3 location and filename structure against the SOP::

    s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/{sublibrary}/{wafer}/
        {wafer}-{sublibrary}[_{well}]_{sublibrary type}-{UG}-{barcode}.trim.*

Known-good examples (both must produce zero violations)::

    czi-hamazaki  hamazaki-seahub-bcp/CHEM3-R100/raw/R100E/441389/
        441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram
    czi-trapnell  trapnell-seahub-bcp/REF3/raw/REF3_P05_2/436830/
        436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram

Design notes
------------
* Every rule is evaluated independently, and the token-level rules run even
  when the stem does not match :data:`SEAHUB_STEM_RE`.  A wholly non-compliant
  name therefore reports each of its defects by name instead of collapsing into
  a single ``unparseable_stem``.
* There is deliberately no "ExperimentID must not appear in the filename" rule.
  The ExperimentID legitimately appears inside sublibrary names such as
  ``REF3_P05_2``, and is absent entirely from the hamazaki example, so
  sublibrary agreement is the only rule governing that token.  That asymmetry is
  also why :func:`seahub_group_parts` tries the folder as-is *before* trying
  ``{ExperimentID}_{folder}``.
* A folder/filename disagreement is attributed to the **folder**, and only the
  folder, when the filename token is exactly ``{ExperimentID}_{folder}``.  The
  untrimmed vendor delivery is the source of truth for the sublibrary name, and
  it carries the full form -- so a correctly-named file under a plate-only
  folder is a folder defect (``sublibrary_folder_truncated``), not a filename
  one.  ``sublibrary_mismatch`` is reserved for names the folder cannot explain
  either way.
* A repeated leading wafer token is normalized away *before* any token rule
  runs, and the rules then see only the normalized stem.  On the raw stem the
  type-less pattern swallows the second wafer into the group
  (``437120-REF3_P04_1_A1``), which no folder can be reconciled against;
  validating both forms would instead double every row.  ``repeated_token``
  stays basename-scoped and is now the fallback for repeats that are not the
  wafer -- comparing filename tokens against folder segments would flag the
  valid ``REF3_P05_2_A10`` under ExperimentID folder ``REF3`` as a duplicate.
* ``expected_name`` repairs only the defect its own rule describes, computed on
  the normalized stem; ``expected_folder`` is set only by the folder rule.
  Composing them into a single corrected key is :mod:`qa_seahub_rename`'s job,
  because that needs the vendor index to resolve a missing sublibrary type.
* ``scope`` controls how far a fact reaches, and therefore how it dedupes:
  ``object`` stays per object, ``stem`` collapses per well, ``folder``
  collapses per sublibrary directory across every wafer and well beneath it,
  and ``upload`` collapses to one row for the whole listing, since the bucket
  and project are the same fact however many wells sit under them.
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
    defect, computed on the wafer-normalized stem.  It is deliberately not the
    fully-corrected name: a well with four defects yields four rows, each naming
    its own repair, and composing them into one target key is
    :mod:`qa_seahub_rename`'s job.

    ``expected_folder`` is the corrected ``{sublibrary}`` path segment, set only
    by ``sublibrary_folder_truncated``.

    ``scope`` says how widely the fact applies and drives dedup when reporting:
    ``object`` for a rule about one S3 object, ``stem`` for one about a well,
    ``folder`` for one about a whole sublibrary directory, and ``upload`` for one
    about the bucket or project, which is the same fact however many wells sit
    under it.  A folder or upload defect reported per object would bury every
    other finding beneath it.
    """

    type: str
    s3_path: str
    detail: str
    expected_name: str = ""
    expected_folder: str = ""
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


def _check_repeated_tokens(stem: str, s3_path: str) -> list[SopViolation]:
    tokens = _split_tokens(stem)
    seen: dict[str, int] = {}
    for token in tokens:
        seen[token] = seen.get(token, 0) + 1
    repeated = sorted(token for token, count in seen.items() if count > 1)
    if not repeated:
        return []
    return [
        SopViolation(
            type="repeated_token",
            s3_path=s3_path,
            detail=(
                f"token(s) {', '.join(repeated)} appear more than once in {stem!r}"
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
) -> tuple[str, str, str]:
    """Reconcile a filename's sublibrary token against its folder.

    Returns ``(implied_folder, trailing_token, state)`` with ``state`` one of:

    * ``ok`` — the folder already matches the filename
    * ``truncated`` — the filename says ``{ExperimentID}_{folder}``, so the
      *folder* is missing its ExperimentID prefix
    * ``mismatch`` — the two cannot be reconciled either way

    Candidate order is load-bearing: the compliant folder is tried first, so a
    name whose folder is already correct never reports as truncated.  That is
    what separates ``REF3_P05_2/…REF3_P05_2_A10…`` (ok) and ``R100E/…R100E…``
    (ok, no ExperimentID in the name at all) from ``P05_1/…REF3_P05_1_A10…``
    (truncated).
    """
    candidates = [(sublibrary, "ok")]
    if experiment_id:
        candidates.append((f"{experiment_id}_{sublibrary}", "truncated"))

    for candidate, state in candidates:
        if group == candidate:
            return candidate, "", state
        prefix = f"{candidate}_"
        if group.startswith(prefix):
            return candidate, group[len(prefix) :], state
    return sublibrary, "", "mismatch"


def _check_group_token(
    group: str, sublibrary: str, experiment_id: str, s3_path: str
) -> list[SopViolation]:
    """Attribute a folder/filename disagreement to whichever side is wrong.

    The vendor delivery is the source of truth for the sublibrary name, so when
    the filename carries the full ``{ExperimentID}_{sublibrary}`` and the folder
    carries only part of it, the folder is at fault.  ``sublibrary_mismatch``
    then means what it says: the two genuinely disagree.
    """
    implied_folder, trailing, state = seahub_group_parts(
        group, sublibrary, experiment_id
    )
    violations: list[SopViolation] = []

    if state == "mismatch":
        return [
            SopViolation(
                type="sublibrary_mismatch",
                s3_path=s3_path,
                detail=(
                    f"filename sublibrary {group!r} is neither folder sublibrary "
                    f"{sublibrary!r} nor {sublibrary!r}_<well>; rename either the "
                    "folder or the files so they agree"
                ),
            )
        ]

    # Orthogonal to the folder question, so it can co-occur with truncation.
    if trailing and not SEAHUB_WELL_RE.match(trailing):
        violations.append(
            SopViolation(
                type="bad_well",
                s3_path=s3_path,
                detail=(
                    f"token {trailing!r} after sublibrary {implied_folder!r} is "
                    "not a well of the form [A-H]<1-2 digits>"
                ),
            )
        )

    if state == "truncated":
        violations.append(
            SopViolation(
                type="sublibrary_folder_truncated",
                s3_path=s3_path,
                detail=(
                    f"folder sublibrary {sublibrary!r} is missing its "
                    f"ExperimentID prefix; the filename says {implied_folder!r}, "
                    f"so rename the folder to {implied_folder!r}"
                ),
                expected_folder=implied_folder,
                scope="folder",
            )
        )
    return violations


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
    if lab and project and project.split("-")[0] != lab:
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

    violations.extend(_check_repeated_tokens(stem, s3_path))

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
        if vendor_assay:
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

    ``stems`` can hold more than one value: a well whose artifacts disagree about
    the name is itself the finding.
    """

    bucket: str
    raw_dir: str
    experiment_id: str
    sublibrary: str
    wafer_folder: str
    identity: tuple[str, str]
    keys: tuple[str, ...]
    suffixes: tuple[str, ...]
    families: frozenset[str]
    stems: tuple[str, ...]
    normalized_stem: str
    wafer: str
    ug: str
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
            _folder, trailing, _state = seahub_group_parts(
                match.group("group"),
                path_info.get("sublibrary", ""),
                path_info.get("experiment_id", ""),
            )
            well_id = trailing
        groups.append(
            SeahubStemGroup(
                bucket=bucket,
                raw_dir=raw_dir,
                experiment_id=path_info.get("experiment_id", ""),
                sublibrary=path_info.get("sublibrary", ""),
                wafer_folder=path_info.get("wafer", ""),
                identity=identity,
                keys=tuple(m[0] for m in members),
                suffixes=tuple(m[2] for m in members),
                families=frozenset(m[3] for m in members),
                stems=tuple(sorted({m[1] for m in members})),
                normalized_stem=normalized,
                wafer=match.group("wafer") if match else "",
                ug=match.group("ug") if match else "",
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

    Dedup is by ``(type, stem, family, expected_folder)``.  Deduping on
    ``expected_name`` instead would not collapse anything -- three of the rules
    carry a suffix-dependent name -- turning a six-artifact well into nineteen
    rows.
    """
    violations: list[SopViolation] = []
    seen: set[tuple[str, str, str, str]] = set()
    for s3_key in group.keys:
        parsed = seahub_stem_and_family(s3_key.split("/")[-1])
        stem, family = (parsed[0], parsed[2]) if parsed is not None else ("", "")
        for violation in validate_seahub_key(
            group.bucket, s3_key, assay_by_identity=assay_by_identity
        ):
            fingerprint = (violation.type, stem, family, violation.expected_folder)
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

    Four scopes collapse differently: ``object`` rules stay per object, ``stem``
    rules collapse to one row per well, ``folder`` rules collapse to one row per
    sublibrary directory -- deliberately ignoring the wafer, since a truncated
    folder name is one fact about a sublibrary however many wafers and wells sit
    beneath it -- and ``upload`` rules collapse to one row for the whole listing.
    Without that, REF3's seven truncated folders would report as several hundred
    rows, and a single wrong bucket would report once per well, either of which
    buries everything else.
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

    seen_folder_facts: set[tuple[str, str, str]] = set()
    for group in groups:
        for violation in validate_seahub_group(
            group, assay_by_identity=assay_by_identity
        ):
            if not _upload_fact_is_new(violation):
                continue
            if violation.scope == "folder":
                fingerprint = (
                    violation.type,
                    f"{group.raw_dir.rsplit('/', 1)[0]}",
                    violation.expected_folder,
                )
                if fingerprint in seen_folder_facts:
                    continue
                seen_folder_facts.add(fingerprint)
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
