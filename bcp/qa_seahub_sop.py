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
  sublibrary agreement is the only rule governing that token.
* ``repeated_token`` is scoped to the basename.  Comparing filename tokens
  against folder segments would flag the valid ``REF3_P05_2_A10`` under
  ExperimentID folder ``REF3`` as a duplicate.
* Violations are non-fatal: QA continues and reports them as a table.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from qa_constants import (
    SEAHUB_BARCODE_RE,
    SEAHUB_BARE_TO_TRIM_SUFFIX,
    SEAHUB_STEM_NO_TYPE_RE,
    SEAHUB_STEM_RE,
    SEAHUB_SUBLIBRARY_TYPES,
    SEAHUB_UG_RE,
    SEAHUB_WELL_RE,
)
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family

__all__ = [
    "SopViolation",
    "expected_sop_name",
    "sop_violation_summary",
    "validate_seahub_key",
    "validate_seahub_stems",
]


@dataclass(frozen=True)
class SopViolation:
    """One SOP rule broken by one object (or stem)."""

    type: str
    s3_path: str
    detail: str
    expected_name: str = ""

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


def _check_group_token(group: str, sublibrary: str, s3_path: str) -> list[SopViolation]:
    """Group must be the folder sublibrary, optionally plus a well token."""
    if group == sublibrary:
        return []
    prefix = f"{sublibrary}_"
    if group.startswith(prefix):
        well = group[len(prefix) :]
        if SEAHUB_WELL_RE.match(well):
            return []
        return [
            SopViolation(
                type="bad_well",
                s3_path=s3_path,
                detail=(
                    f"token {well!r} after sublibrary {sublibrary!r} is not a "
                    "well of the form [A-H]<1-2 digits>"
                ),
            )
        ]
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


def validate_seahub_key(bucket: str, s3_key: str) -> list[SopViolation]:
    """Validate one SeaHub object against the SOP.

    Returns every rule broken by this object; an empty list means compliant.
    """
    s3_path = f"s3://{bucket}/{s3_key}"
    violations, path_info = _check_path(bucket, s3_key)
    if path_info is None:
        return violations

    basename = s3_key.split("/")[-1]
    parsed = seahub_stem_and_family(basename)
    if parsed is None:
        violations.append(
            SopViolation(
                type="unexpected_suffix",
                s3_path=s3_path,
                detail=(
                    f"{basename!r} does not end in a SeaHub artifact suffix; "
                    "expected one of the six .trim.* artifacts"
                ),
            )
        )
        return violations

    stem, suffix, family = parsed
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
        violations.append(
            SopViolation(
                type="invalid_sublibrary_type",
                s3_path=s3_path,
                detail=(
                    f"{stem!r} carries no sublibrary type; expected one of "
                    f"{', '.join(SEAHUB_SUBLIBRARY_TYPES)}"
                ),
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
        _check_group_token(match.group("group"), path_info["sublibrary"], s3_path)
    )

    if not SEAHUB_UG_RE.match(match.group("ug")):
        violations.append(
            SopViolation(
                type="bad_ug",
                s3_path=s3_path,
                detail=f"UG {match.group('ug')!r} is not of the form Z####",
            )
        )
    if not SEAHUB_BARCODE_RE.match(match.group("barcode")):
        violations.append(
            SopViolation(
                type="bad_barcode",
                s3_path=s3_path,
                detail=(
                    f"barcode {match.group('barcode')!r} contains characters "
                    "outside A/C/G/T"
                ),
            )
        )
    return violations


def validate_seahub_stems(bucket: str, s3_keys: list[str]) -> list[SopViolation]:
    """Validate a listing, reporting each stem once rather than each artifact.

    A wholly misnamed well has six artifacts sharing one stem; without this
    collapsing, one defect would be reported six times.  Artifact-specific
    rules (``missing_trim_infix``, ``unexpected_suffix``) still need per-object
    evaluation, so the first object of each stem is validated in full and its
    ``missing_trim_infix`` is reported once for the whole stem.
    """
    violations: list[SopViolation] = []
    seen_stems: set[tuple[str, str]] = set()
    for s3_key in sorted(s3_keys):
        parsed = seahub_stem_and_family(s3_key.split("/")[-1])
        if parsed is None:
            violations.extend(validate_seahub_key(bucket, s3_key))
            continue
        stem_id = ("/".join(s3_key.split("/")[:-1]), parsed[0])
        if stem_id in seen_stems:
            continue
        seen_stems.add(stem_id)
        violations.extend(validate_seahub_key(bucket, s3_key))
    return violations


def sop_violation_summary(violations: list[SopViolation]) -> dict[str, int]:
    """Count violations by type, for a compact notebook summary."""
    counts: dict[str, int] = {}
    for violation in violations:
        counts[violation.type] = counts.get(violation.type, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
