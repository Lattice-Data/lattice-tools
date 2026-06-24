"""
SOP-driven per-prefix file-completeness validators for raw mappings.

Each validator walks a list of :class:`MappingRow` objects, classifies each
row's S3 suffix against the SOP-required artifact tables in
:mod:`mapping_validation.sop_artifacts`, and reports:

- ``unexpected_sample_file`` for S3 suffixes not in the SOP table (this is
  the signal that catches truncated names like ``_codes.csv`` instead of
  ``_trimmer-failure_codes.csv``)
- ``s3_local_artifact_mismatch`` when the local basename's SOP suffix
  contradicts the S3 row's suffix
- ``missing_sample_artifacts`` when a per-prefix bundle is missing any
  SOP-required artifact

Per-RunID metadata files (e.g. ``merged_trimmer-failure_codes.csv``) are
recognised by :data:`mapping_validation.constants._RUN_METADATA_RE` and
skipped here; their presence is enforced by the family-specific
S3-SOP regex check in :mod:`mapping_validation.validators`.
"""

from __future__ import annotations

import os
import re
from collections import defaultdict
from typing import Iterable, List, Tuple

from .parsing import MappingRow
from .sop_artifacts import (
    SCALE_PER_UG_SUFFIX_TO_ARTIFACT,
    SCALE_PER_WELL_EXT_TO_ARTIFACT,
    SCALE_REQUIRED_PER_UG_ARTIFACTS,
    SCALE_REQUIRED_PER_WELL_ARTIFACTS,
    SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT,
    SCI_REQUIRED_PER_PREFIX_ARTIFACTS,
    TENX_FASTQ_OPTIONAL_SUFFIX_TO_ARTIFACT,
    TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT,
    TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT,
    TENX_FASTQ_REQUIRED_PER_PREFIX_ARTIFACTS,
    TENX_FASTQ_REQUIRED_PER_TAIL_ARTIFACTS,
    classify_basename,
    classify_suffix,
    expected_suffixes_message,
    is_scale_per_well_optional_suffix,
    scale_classify_suffix,
)
from .validators import (
    _build_seahub_s3_patterns,
    _check_s3_local_suffix_consistency,
    _get_10x_s3_regex,
    _is_run_metadata,
)


# ---------------------------------------------------------------------------
# Shared helpers for reconstructing the concrete names of missing files
# ---------------------------------------------------------------------------


def _reverse_artifact_suffixes(table: dict[str, str]) -> dict[str, str]:
    """Map each artifact name back to its canonical SOP suffix.

    The forward tables map ``suffix -> artifact``; several suffix variants can
    share an artifact only for non-required metadata files, so keeping the
    first occurrence yields the canonical suffix for every required artifact.
    """
    reverse: dict[str, str] = {}
    for suffix, artifact in table.items():
        reverse.setdefault(artifact, suffix)
    return reverse


def _expected_missing_files(
    stem: str | None, missing: Iterable[str], artifact_to_suffix: dict[str, str]
) -> list[str]:
    """Reconstruct the SOP-expected S3 names for a bundle's missing artifacts.

    ``stem`` is the shared S3 path of an existing sibling row with its own
    suffix stripped, so appending the missing artifact's canonical suffix
    yields the exact path that should exist but does not.
    """
    if stem is None:
        return []
    return [
        stem + artifact_to_suffix[artifact]
        for artifact in missing
        if artifact in artifact_to_suffix
    ]


# ---------------------------------------------------------------------------
# Scale completeness
# ---------------------------------------------------------------------------


def _scale_expected_suffixes_msg() -> str:
    """Diagnostic message listing every SOP-recognised Scale suffix shape."""
    aggregate = expected_suffixes_message(SCALE_PER_UG_SUFFIX_TO_ARTIFACT)
    per_well = ", ".join(
        f"_<wellcode>{ext}" for ext in sorted(SCALE_PER_WELL_EXT_TO_ARTIFACT)
    )
    return f"per-well: {per_well}; per-UG: {aggregate}"


def _classify_scale_local_basename(
    basename: str,
) -> tuple[str, str, str | None] | None:
    """Classify a local Scale basename as per-well or per-UG aggregate.

    Returns the same shape as :func:`scale_classify_suffix` so the two sides
    can be compared.
    """
    well_match = re.search(
        r"[_-](?P<wellcode>\d{1,2}[A-H])(?P<ext>\.(?:cram|csv|json))$", basename
    )
    if well_match:
        artifact = SCALE_PER_WELL_EXT_TO_ARTIFACT.get(well_match.group("ext"))
        if artifact is not None:
            return ("per_well", artifact, well_match.group("wellcode"))

    aggregate = classify_basename(basename, SCALE_PER_UG_SUFFIX_TO_ARTIFACT)
    if aggregate is not None:
        return ("per_ug", aggregate, None)
    return None


def validate_scale_raw_completeness(mappings: Iterable[MappingRow]) -> dict:
    """Enforce SOP-listed file completeness for Scale raw S3 mappings.

    For each Scale S3 path matched by the SOP regex, classify the suffix
    after the captured ``ug_rt`` token against
    :data:`mapping_validation.sop_artifacts.SCALE_PER_WELL_EXT_TO_ARTIFACT`
    and :data:`SCALE_PER_UG_SUFFIX_TO_ARTIFACT`.  Aggregate per
    ``(group_id, runid, assay, ug_rt)`` for per-UG aggregates and
    ``(group_id, runid, assay, ug_rt, well_code)`` for per-well artifacts.

    After the scan, every per-well bundle must contain ``cram``/``csv``/
    ``json`` and every per-UG bundle must contain the five aggregate
    artifacts.  Missing items become ``missing_sample_artifacts`` errors.

    Cross-checks each row with :func:`_check_s3_local_suffix_consistency`
    so e.g. local ``_trimmer-failure_codes.csv`` paired with S3
    ``_codes.csv`` is flagged.
    """
    patterns = _build_seahub_s3_patterns("scale")

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0

    well_artifacts: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(set)
    aggregate_artifacts: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    well_lines: dict[tuple[str, str, str, str, str], set[int]] = defaultdict(set)
    aggregate_lines: dict[tuple[str, str, str, str], set[int]] = defaultdict(set)
    well_stem: dict[tuple[str, str, str, str, str], str] = {}
    aggregate_stem: dict[tuple[str, str, str, str], str] = {}
    well_rev = _reverse_artifact_suffixes(SCALE_PER_WELL_EXT_TO_ARTIFACT)

    for row in mappings:
        total += 1
        if _is_run_metadata(row.s3_path):
            metadata_count += 1
            continue

        m: re.Match[str] | None = None
        form = ""
        for pat, pat_form in patterns:
            m = pat.match(row.s3_path)
            if m:
                form = pat_form
                break
        if not m:
            continue

        # Only the SOP form (``_QSR-N[-SCALEPLEX]<suffix>``) carries SOP
        # artifact suffixes.  The ``index`` form uses index-sequence-based
        # naming and is out of scope for SOP completeness.
        if form != "sop":
            continue

        gd = m.groupdict()
        # The suffix capture group includes everything after the ug_rt token.
        suffix = gd.get("suffix", "")
        if suffix == "":
            errors.append(
                {
                    "type": "unexpected_sample_file",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": (
                        "missing Scale sample file suffix; expected one of "
                        f"{_scale_expected_suffixes_msg()}"
                    ),
                }
            )
            continue

        classification = scale_classify_suffix(suffix)
        if classification is None:
            if is_scale_per_well_optional_suffix(suffix):
                # Recognised provider-specific sidecar (e.g. _10A.cram-metadata.json);
                # accept silently rather than flooding errors with one row per well.
                continue
            errors.append(
                {
                    "type": "unexpected_sample_file",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": (
                        f"unexpected Scale sample file suffix '{suffix}'; "
                        f"expected one of {_scale_expected_suffixes_msg()}"
                    ),
                }
            )
            mismatch = _check_s3_local_suffix_consistency(
                row,
                None,
                _scale_local_suffix_table(),
            )
            if mismatch is not None:
                errors.append(mismatch)
            continue

        matched += 1
        scope, artifact, well_code = classification
        group_id = gd["group_id"]
        runid = gd["runid"]
        assay = gd["assay"]
        ug_rt = gd["ug_rt"]
        ug_key = (group_id, runid, assay, ug_rt)

        if scope == "per_well":
            assert well_code is not None
            well_key = ug_key + (well_code,)
            well_artifacts[well_key].add(artifact)
            well_lines[well_key].add(row.line_num)
            ext = well_rev.get(artifact)
            if ext:
                well_stem.setdefault(
                    well_key, row.s3_path[: len(row.s3_path) - len(ext)]
                )
        else:
            aggregate_artifacts[ug_key].add(artifact)
            aggregate_lines[ug_key].add(row.line_num)
            if suffix:
                aggregate_stem.setdefault(
                    ug_key, row.s3_path[: len(row.s3_path) - len(suffix)]
                )

        local_check = _check_scale_local_consistency(row, classification)
        if local_check is not None:
            errors.append(local_check)

    missing_per_well: dict[str, list[str]] = {}
    for key, found in well_artifacts.items():
        missing = sorted(SCALE_REQUIRED_PER_WELL_ARTIFACTS - found)
        if missing:
            label = f"{key[1]}-{key[0]}_{key[2]}_{key[3]}_{key[4]} (per-well artifacts)"
            missing_per_well[label] = missing
            errors.append(
                {
                    "type": "missing_sample_artifacts",
                    "s3_path": label,
                    "detail": (
                        "missing required per-well Scale artifacts: "
                        + ", ".join(missing)
                    ),
                    "missing": missing,
                    "present": sorted(found),
                    "lines": sorted(well_lines.get(key, set())),
                    "missing_files": _expected_missing_files(
                        well_stem.get(key), missing, well_rev
                    ),
                }
            )

    ug_rev = _reverse_artifact_suffixes(SCALE_PER_UG_SUFFIX_TO_ARTIFACT)
    missing_per_ug: dict[str, list[str]] = {}
    for key, found in aggregate_artifacts.items():
        missing = sorted(SCALE_REQUIRED_PER_UG_ARTIFACTS - found)
        if missing:
            label = f"{key[1]}-{key[0]}_{key[2]}_{key[3]} (per-UG aggregate artifacts)"
            missing_per_ug[label] = missing
            errors.append(
                {
                    "type": "missing_sample_artifacts",
                    "s3_path": label,
                    "detail": (
                        "missing required per-UG Scale artifacts: " + ", ".join(missing)
                    ),
                    "missing": missing,
                    "present": sorted(found),
                    "lines": sorted(aggregate_lines.get(key, set())),
                    "missing_files": _expected_missing_files(
                        aggregate_stem.get(key), missing, ug_rev
                    ),
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "total": total,
        "matched": matched,
        "well_prefixes_checked": len(well_artifacts),
        "ug_aggregates_checked": len(aggregate_artifacts),
        "missing_per_well": missing_per_well,
        "missing_per_ug": missing_per_ug,
    }


def _scale_local_suffix_table() -> dict[str, str]:
    """Combined Scale suffix table (well-code shapes can't be tabulated literally).

    Used only for cross-side mismatch detection on the local basename; the
    well-code variants are handled by :func:`_classify_scale_local_basename`.
    """
    return dict(SCALE_PER_UG_SUFFIX_TO_ARTIFACT)


def _check_scale_local_consistency(
    row: MappingRow,
    s3_classification: tuple[str, str, str | None],
) -> dict | None:
    """Compare classified S3 vs classified local for a Scale row.

    Returns an ``s3_local_artifact_mismatch`` error when the local basename
    is recognised as a SOP artifact that disagrees with the S3 side, taking
    the per-well / per-UG scope into account.
    """
    s3_scope, s3_artifact, s3_well_code = s3_classification
    local_basename = os.path.basename(row.local_path)
    local_classification = _classify_scale_local_basename(local_basename)
    if local_classification is None:
        return None
    local_scope, local_artifact, local_well_code = local_classification

    if s3_scope != local_scope:
        return {
            "type": "s3_local_artifact_mismatch",
            "line": row.line_num,
            "s3_path": row.s3_path,
            "local_path": row.local_path,
            "detail": (
                f"S3 is {s3_scope} ('{s3_artifact}') but local is "
                f"{local_scope} ('{local_artifact}')"
            ),
        }
    if s3_artifact != local_artifact:
        return {
            "type": "s3_local_artifact_mismatch",
            "line": row.line_num,
            "s3_path": row.s3_path,
            "local_path": row.local_path,
            "detail": (
                f"S3 artifact '{s3_artifact}' maps to local artifact '{local_artifact}'"
            ),
        }
    if s3_scope == "per_well" and s3_well_code != local_well_code:
        return {
            "type": "s3_local_artifact_mismatch",
            "line": row.line_num,
            "s3_path": row.s3_path,
            "local_path": row.local_path,
            "detail": (
                f"S3 well code '{s3_well_code}' does not match local well "
                f"code '{local_well_code}'"
            ),
        }
    return None


# ---------------------------------------------------------------------------
# sci completeness
# ---------------------------------------------------------------------------


def _sci_expected_suffixes_msg() -> str:
    return expected_suffixes_message(SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT)


def validate_sci_raw_completeness(mappings: Iterable[MappingRow]) -> dict:
    """Enforce SOP-listed file completeness for sci raw S3 mappings.

    Each S3 path matched by the sci SOP regex is classified by its suffix.
    The bundle key is ``(group_id, runid, assay, ug, barcode)``; every key
    must include the seven SOP-required artifacts (cram + csv/json + four
    trimmer-related files).
    """
    patterns = _build_seahub_s3_patterns("sci")

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0
    bundle_artifacts: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(set)
    bundle_lines: dict[tuple[str, str, str, str, str], set[int]] = defaultdict(set)
    bundle_stem: dict[tuple[str, str, str, str, str], str] = {}

    for row in mappings:
        total += 1
        if _is_run_metadata(row.s3_path):
            metadata_count += 1
            continue

        m: re.Match[str] | None = None
        for pat, _ in patterns:
            m = pat.match(row.s3_path)
            if m:
                break
        if not m:
            continue

        gd = m.groupdict()
        suffix = gd.get("suffix", "")
        artifact = classify_suffix(suffix, SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT)
        if artifact is None:
            errors.append(
                {
                    "type": "unexpected_sample_file",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": (
                        f"unexpected sci sample file suffix '{suffix}'; "
                        f"expected one of {_sci_expected_suffixes_msg()}"
                    ),
                }
            )
            mismatch = _check_s3_local_suffix_consistency(
                row,
                None,
                SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT,
            )
            if mismatch is not None:
                errors.append(mismatch)
            continue

        matched += 1
        key = (
            gd["group_id"],
            gd["runid"],
            gd["assay"],
            gd["ug"],
            gd["barcode"],
        )
        bundle_artifacts[key].add(artifact)
        bundle_lines[key].add(row.line_num)
        if suffix:
            bundle_stem.setdefault(key, row.s3_path[: len(row.s3_path) - len(suffix)])

        mismatch = _check_s3_local_suffix_consistency(
            row,
            artifact,
            SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT,
        )
        if mismatch is not None:
            errors.append(mismatch)

    sci_rev = _reverse_artifact_suffixes(SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT)
    missing_by_prefix: dict[str, list[str]] = {}
    for key, found in bundle_artifacts.items():
        missing = sorted(SCI_REQUIRED_PER_PREFIX_ARTIFACTS - found)
        if missing:
            label = f"{key[1]}-{key[0]}_{key[2]}-{key[3]}-{key[4]}"
            missing_by_prefix[label] = missing
            errors.append(
                {
                    "type": "missing_sample_artifacts",
                    "s3_path": label,
                    "detail": (
                        "missing required sci sample artifacts: " + ", ".join(missing)
                    ),
                    "missing": missing,
                    "present": sorted(found),
                    "lines": sorted(bundle_lines.get(key, set())),
                    "missing_files": _expected_missing_files(
                        bundle_stem.get(key), missing, sci_rev
                    ),
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "total": total,
        "matched": matched,
        "sample_prefixes_checked": len(bundle_artifacts),
        "missing_sample_artifacts": missing_by_prefix,
    }


# ---------------------------------------------------------------------------
# 10x raw FASTQ completeness
# ---------------------------------------------------------------------------


def _tenx_fastq_expected_per_prefix_msg() -> str:
    return expected_suffixes_message(TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT)


def _tenx_fastq_expected_per_tail_msg() -> str:
    return expected_suffixes_message(TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT)


def _split_tenx_fastq_suffix(
    suffix: str,
) -> Tuple[str | None, str | None, str | None]:
    """Split a 10x FASTQ row suffix into (per_tail_suffix, fastq_tail, per_prefix_suffix).

    Returns ``(tail_only_suffix, fastq_tail_text, None)`` for rows whose
    suffix starts with an Illumina tail (e.g. ``_S1_L001_R1_001.fastq.gz``);
    returns ``(None, None, suffix)`` for per-prefix rows (e.g.
    ``_trimmer-failure_codes.csv``); returns ``(None, None, None)`` when
    the suffix matches neither classification.
    """
    tail_match = re.match(
        r"^_(?P<tail>S\d+_L\d+_(?:R[123]|I[12])_\d+)(?P<rest>.+)$", suffix
    )
    if tail_match:
        return tail_match.group("rest"), tail_match.group("tail"), None
    return None, None, suffix


def validate_10x_raw_fastq_completeness(
    provider: str, mappings: Iterable[MappingRow]
) -> dict:
    """Enforce SOP-listed file completeness for 10x raw FASTQ mappings.

    For each S3 row matched by :func:`_get_10x_s3_regex`, classify the
    suffix that follows the ``-Z####-BARCODE`` prefix:

    - If the suffix is an Illumina tail variant (``_S*_L*_R*_*.csv`` etc.),
      record it under per-prefix → per-tail bookkeeping and require the
      four per-tail SOP artifacts (.fastq.gz, .csv, .json, _sample.fastq.gz).
    - Otherwise treat it as a per-prefix metadata suffix (.csv, .json,
      _trimmer-*.csv, _unmatched.*) and require the seven per-prefix SOP
      artifacts.

    CRAM rows (10x_cram mode) are ignored here — they have a separate
    validator (:func:`validate_s3_10x_cram_raw`).  Optional STAR Solo and
    application-QC suffixes are recognised silently and not required.
    """
    s3_regex, _ = _get_10x_s3_regex(provider)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0
    skipped_non_fastq = 0
    per_prefix_artifacts: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(
        set
    )
    per_tail_artifacts: dict[tuple[str, str, str, str, str, str], set[str]] = (
        defaultdict(set)
    )
    per_prefix_tails: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(set)
    per_prefix_lines: dict[tuple[str, str, str, str, str], set[int]] = defaultdict(set)
    per_tail_lines: dict[tuple[str, str, str, str, str, str], set[int]] = defaultdict(
        set
    )
    per_prefix_stem: dict[tuple[str, str, str, str, str], str] = {}
    per_tail_stem: dict[tuple[str, str, str, str, str, str], str] = {}

    for row in mappings:
        total += 1
        s3 = row.s3_path
        if _is_run_metadata(s3):
            metadata_count += 1
            continue

        base = os.path.basename(s3)
        # Per-sample ``.cram`` and ``.cram-metadata.json`` rows belong to the
        # 10x_cram mode and are skipped by the FASTQ-mode completeness check.
        # ``_unmatched.cram[-metadata.json]`` is part of the SOP-required
        # FASTQ-mode bundle and must NOT be skipped.
        if (
            base.endswith(".cram") or base.endswith(".cram-metadata.json")
        ) and "_unmatched.cram" not in base:
            skipped_non_fastq += 1
            continue

        m = s3_regex.match(s3)
        if not m:
            continue

        gd = m.groupdict()
        suffix = gd.get("suffix", "")
        if suffix == "":
            errors.append(
                {
                    "type": "unexpected_sample_file",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": (
                        "missing 10x FASTQ sample file suffix; expected per-tail "
                        f"({_tenx_fastq_expected_per_tail_msg()}) or per-prefix "
                        f"({_tenx_fastq_expected_per_prefix_msg()})"
                    ),
                }
            )
            continue

        tail_only, fastq_tail, per_prefix_suffix = _split_tenx_fastq_suffix(suffix)
        prefix_key = (
            gd["groupid"],
            gd["runid"],
            gd["assay"],
            gd["ug"],
            gd["barcode"],
        )

        if fastq_tail is not None and tail_only is not None:
            artifact = classify_suffix(
                tail_only, TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT
            )
            if artifact is None:
                errors.append(
                    {
                        "type": "unexpected_sample_file",
                        "line": row.line_num,
                        "s3_path": s3,
                        "detail": (
                            f"unexpected 10x FASTQ per-tail suffix '{tail_only}' "
                            f"after Illumina tail '_{fastq_tail}'; expected one "
                            f"of {_tenx_fastq_expected_per_tail_msg()}"
                        ),
                    }
                )
                continue
            matched += 1
            tail_key = prefix_key + (fastq_tail,)
            per_tail_artifacts[tail_key].add(artifact)
            per_tail_lines[tail_key].add(row.line_num)
            per_tail_stem.setdefault(tail_key, s3[: len(s3) - len(tail_only)])
            per_prefix_tails[prefix_key].add(fastq_tail)
            continue

        # per-prefix bucket
        assert per_prefix_suffix is not None
        artifact = classify_suffix(
            per_prefix_suffix, TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT
        )
        if artifact is None:
            optional_artifact = classify_suffix(
                per_prefix_suffix, TENX_FASTQ_OPTIONAL_SUFFIX_TO_ARTIFACT
            )
            if optional_artifact is not None:
                continue
            errors.append(
                {
                    "type": "unexpected_sample_file",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": (
                        f"unexpected 10x FASTQ per-prefix suffix '{per_prefix_suffix}'; "
                        f"expected one of {_tenx_fastq_expected_per_prefix_msg()}"
                    ),
                }
            )
            mismatch = _check_s3_local_suffix_consistency(
                row,
                None,
                TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT,
            )
            if mismatch is not None:
                errors.append(mismatch)
            continue

        matched += 1
        per_prefix_artifacts[prefix_key].add(artifact)
        per_prefix_lines[prefix_key].add(row.line_num)
        per_prefix_stem.setdefault(prefix_key, s3[: len(s3) - len(per_prefix_suffix)])

        mismatch = _check_s3_local_suffix_consistency(
            row,
            artifact,
            TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT,
        )
        if mismatch is not None:
            errors.append(mismatch)

    per_prefix_rev = _reverse_artifact_suffixes(
        TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT
    )
    per_tail_rev = _reverse_artifact_suffixes(TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT)

    missing_per_prefix: dict[str, list[str]] = {}
    for key, found in per_prefix_artifacts.items():
        missing = sorted(TENX_FASTQ_REQUIRED_PER_PREFIX_ARTIFACTS - found)
        if missing:
            label = (
                f"{key[1]}-{key[0]}_{key[2]}-{key[3]}-{key[4]} (per-prefix artifacts)"
            )
            missing_per_prefix[label] = missing
            errors.append(
                {
                    "type": "missing_sample_artifacts",
                    "s3_path": label,
                    "detail": (
                        "missing required 10x FASTQ per-prefix artifacts: "
                        + ", ".join(missing)
                    ),
                    "missing": missing,
                    "present": sorted(found),
                    "lines": sorted(per_prefix_lines.get(key, set())),
                    "missing_files": _expected_missing_files(
                        per_prefix_stem.get(key), missing, per_prefix_rev
                    ),
                }
            )

    missing_per_tail: dict[str, list[str]] = {}
    for prefix_key, tails in per_prefix_tails.items():
        for tail in sorted(tails):
            tail_key = prefix_key + (tail,)
            found = per_tail_artifacts.get(tail_key, set())
            missing = sorted(TENX_FASTQ_REQUIRED_PER_TAIL_ARTIFACTS - found)
            if missing:
                label = (
                    f"{prefix_key[1]}-{prefix_key[0]}_{prefix_key[2]}-"
                    f"{prefix_key[3]}-{prefix_key[4]}_{tail}"
                )
                missing_per_tail[label] = missing
                errors.append(
                    {
                        "type": "missing_sample_artifacts",
                        "s3_path": label,
                        "detail": (
                            "missing required 10x FASTQ per-tail artifacts: "
                            + ", ".join(missing)
                        ),
                        "missing": missing,
                        "present": sorted(found),
                        "lines": sorted(per_tail_lines.get(tail_key, set())),
                        "missing_files": _expected_missing_files(
                            per_tail_stem.get(tail_key), missing, per_tail_rev
                        ),
                    }
                )

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "skipped_non_fastq": skipped_non_fastq,
        "total": total,
        "matched": matched,
        "sample_prefixes_checked": len(per_prefix_artifacts),
        "tails_checked": sum(len(t) for t in per_prefix_tails.values()),
        "missing_per_prefix": missing_per_prefix,
        "missing_per_tail": missing_per_tail,
    }


__all__ = [
    "validate_scale_raw_completeness",
    "validate_sci_raw_completeness",
    "validate_10x_raw_fastq_completeness",
]
