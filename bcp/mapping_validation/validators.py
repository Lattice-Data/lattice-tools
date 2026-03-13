from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Set, Tuple

from .constants import (
    ASSAYS_BY_FAMILY,
    CANONICAL_ASSAY,
    PROVIDER_BUCKET,
    _RUN_METADATA_RE,
    build_assay_regex,
    get_assays,
    get_order_pattern,
)

_SEAHUB_FAMILIES: set[str] = {"scale", "sci"}
from .parsing import MappingRow
from .sif_io import load_sif_group_assays


def _is_run_metadata(s3_path: str) -> bool:
    """Return True if the S3 path's filename is a known run-level metadata file."""
    return bool(_RUN_METADATA_RE.match(os.path.basename(s3_path)))


def _validate_provider(provider: str) -> str:
    """Normalise and validate a provider name."""
    provider = provider.lower()
    if provider not in PROVIDER_BUCKET:
        raise ValueError(
            f"Unknown provider '{provider}', expected one of {sorted(PROVIDER_BUCKET)}"
        )
    return provider


# ---------------------------------------------------------------------------
# 10x raw S3 validation
# ---------------------------------------------------------------------------


def _build_10x_raw_s3_regex(provider: str, order_pattern: str, assay_re: str) -> re.Pattern[str]:
    """Build the canonical 10x raw S3 regex for a provider."""
    return re.compile(
        rf"^s3://(?P<bucket>czi-{re.escape(provider)})/"
        r"(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/"
        r"(?P<groupid>[^/]+)/raw/"
        rf"(?P<runid>\d+)-(?P<file_stem>.+?)_(?P<assay>{assay_re})-(?P<ug>Z\d{{4}})-(?P<barcode>[A-Za-z]+)"
        r"(?P<suffix>.*)$"
    )


def _get_10x_s3_regex(provider: str) -> Tuple[re.Pattern[str], set[str]]:
    """Return (compiled_regex, valid_assays) for 10x raw S3 paths."""
    provider = _validate_provider(provider)
    valid_assays = get_assays("10x", provider)
    assay_re = build_assay_regex(valid_assays)
    order_pattern = get_order_pattern(provider)
    regex = _build_10x_raw_s3_regex(provider, order_pattern, assay_re)
    return regex, valid_assays


def validate_s3_10x_raw(provider: str, mappings: Iterable[MappingRow]) -> dict:
    """
    Validate 10x-style raw S3 paths against the SOP for a given provider.

    The expected S3 layout is:

        s3://czi-{provider}/{lastname}-{project}/{order}/{GroupID}/raw/
            {RunID}-{GroupID}_{Assay}-{UG-BC}{suffix}

    Where:
    - provider: 'novogene' or 'psomagen'
    - order: provider-specific pattern
    - Assay: family-level assay set, optionally extended by provider
    - UG-BC: Z####-BARCODE (BARCODE is A/C/G/T only)
    """
    provider = _validate_provider(provider)
    s3_regex, valid_assays = _get_10x_s3_regex(provider)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0
    group_assays: dict[str, set[str]] = defaultdict(set)
    run_ids: set[str] = set()

    for row in mappings:
        total += 1

        if _is_run_metadata(row.s3_path):
            metadata_count += 1
            continue

        m = s3_regex.match(row.s3_path)
        if not m:
            warnings.append(
                {
                    "type": "parse_miss",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": "S3 path does not match expected 10x raw pattern",
                }
            )
            continue

        matched += 1
        gd = m.groupdict()
        group_assays[gd["groupid"]].add(gd["assay"])
        run_ids.add(gd["runid"])

        expected_prefix = f"/{gd['groupid']}/raw/{gd['runid']}-{gd['groupid']}_"
        if expected_prefix not in row.s3_path:
            errors.append(
                {
                    "type": "group_mismatch",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"groupid mismatch between path '{gd['groupid']}' and filename prefix",
                }
            )

        assay = gd["assay"]
        if assay not in valid_assays:
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"assay '{assay}' not allowed for provider '{provider}' "
                    f"(expected one of {sorted(valid_assays)})",
                }
            )

        if not re.fullmatch(r"[ACGT]+", gd["barcode"]):
            errors.append(
                {
                    "type": "invalid_barcode",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"barcode '{gd['barcode']}' contains characters outside A/C/G/T",
                }
            )

        project = gd["project"]
        if project != project.lower() or "_" in project:
            warnings.append(
                {
                    "type": "project_naming",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": "project should be lower-case with hyphen delimiters (per SOP)",
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "total": total,
        "matched": matched,
        "group_assays": dict(group_assays),
        "run_ids": run_ids,
    }


def find_unmatched_sif_paths_10x(
    mappings: Iterable[MappingRow],
    sif_groupids: set[str],
    provider: str,
) -> dict:
    """Find S3 paths not covered by any SIF GroupID.

    For each mapping row, extracts the GroupID from the S3 path and
    checks whether it appears in ``sif_groupids``.  Paths whose GroupID
    is absent, and paths that cannot be parsed at all (excluding known
    run-level metadata files), are collected and returned.
    """
    provider = _validate_provider(provider)
    s3_regex, _ = _get_10x_s3_regex(provider)

    unmatched_by_group: dict[str, List[MappingRow]] = defaultdict(list)
    unparsed: List[MappingRow] = []
    total = 0
    matched_sif = 0
    metadata = 0

    for row in mappings:
        total += 1
        if _is_run_metadata(row.s3_path):
            metadata += 1
            continue
        m = s3_regex.match(row.s3_path)
        if not m:
            unparsed.append(row)
            continue
        groupid = m.group("groupid")
        if groupid in sif_groupids:
            matched_sif += 1
        else:
            unmatched_by_group[groupid].append(row)

    return {
        "unmatched_by_group": dict(unmatched_by_group),
        "unparsed": unparsed,
        "total": total,
        "matched_sif": matched_sif,
        "metadata": metadata,
    }


# ---------------------------------------------------------------------------
# Seahub raw S3 validation (Scale + sci unified)
# ---------------------------------------------------------------------------


def _validate_seahub_family(assay_family: str) -> str:
    """Normalise and validate a seahub assay family name."""
    family = assay_family.lower()
    if family not in _SEAHUB_FAMILIES:
        raise ValueError(
            f"Unknown seahub assay family '{assay_family}', "
            f"expected one of {sorted(_SEAHUB_FAMILIES)}"
        )
    return family


def _build_seahub_s3_patterns(
    assay_family: str,
) -> list[Tuple[re.Pattern[str], str]]:
    """Build S3 regex patterns for a seahub assay family.

    Returns a list of (compiled_pattern, form_name) tuples.
    Scale has two forms (sop + index); sci has one (ug-barcode).
    """
    family = _validate_seahub_family(assay_family)
    order_re = get_order_pattern("novogene")
    valid_assays = ASSAYS_BY_FAMILY[family]
    assay_re = build_assay_regex(valid_assays)

    s3_prefix = (
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_re})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
    )

    if family == "scale":
        # Scale GroupIDs may contain underscores/hyphens (e.g. R096G, R112A-B).
        group_id_re = r"[A-Za-z0-9_-]+"
        stem = s3_prefix + rf"(?P<runid2>\d+)-(?P<group_id>{group_id_re})_(?P<assay>{assay_re})"
        sop = re.compile(
            stem
            + r"_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)"
            + r"(?P<suffix>.*)$"
        )
        idx = re.compile(
            stem
            + r"_(?P<index_seq>[ACGT]+)"
            + r"(?P<suffix>.*)$"
        )
        return [(sop, "sop"), (idx, "index")]

    # sci: use non-greedy group_id (.+?) so underscored IDs like
    # CHEM16_P07_F3 are parsed correctly.  The -Z\d{4}-[ACGT]+ anchor
    # after the assay prevents ambiguity.
    sci_pat = re.compile(
        s3_prefix
        + rf"(?P<runid2>\d+)-(?P<group_id>.+?)_(?P<assay>{assay_re})"
        + r"-(?P<ug>Z\d{4})-(?P<barcode>[A-Za-z]+)"
        + r"(?P<suffix>.*)$"
    )
    return [(sci_pat, "sci")]


def _seahub_post_match_checks(
    assay_family: str,
    form: str,
    gd: dict[str, str],
    assay: str,
    row: MappingRow,
    errors: List[dict],
) -> None:
    """Family-specific post-match validation hooks for seahub S3 paths."""
    if assay_family == "scale" and form == "sop":
        ug_rt = gd["ug_rt"]
        has_scaleplex = "SCALEPLEX" in ug_rt
        if assay in {"hash_oligo", "GEX_hash_oligo"} and not has_scaleplex:
            errors.append(
                {
                    "type": "hash_oligo_scaleplex_mismatch",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"assay '{assay}' requires SCALEPLEX in UG_RT, got '{ug_rt}'",
                }
            )
        if assay == "GEX" and has_scaleplex:
            errors.append(
                {
                    "type": "gex_scaleplex_violation",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"assay 'GEX' must not include SCALEPLEX in UG_RT, got '{ug_rt}'",
                }
            )

    elif assay_family == "sci":
        barcode = gd.get("barcode", "")
        if barcode and not re.fullmatch(r"[ACGT]+", barcode):
            errors.append(
                {
                    "type": "invalid_barcode",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"barcode '{barcode}' contains characters outside A/C/G/T",
                }
            )


def validate_s3_seahub_raw(
    assay_family: str, mappings: Iterable[MappingRow]
) -> dict:
    """
    Validate seahub-style (Scale or sci) raw S3 paths against the SOP.

    Shared S3 prefix layout:

        s3://czi-novogene/{project}/{order}/{experiment_id}/raw/
            {RunID}/{RunID}-{GroupID}_{Assay}<family-specific-suffix>

    The suffix differs by family:
    - scale SOP: ``_{Assay}_QSR-{N}[-SCALEPLEX]{suffix}``
    - scale index: ``_{Assay}_{IndexSequence}{suffix}``
    - sci: ``_{Assay}-Z{4digits}-{Barcode}{suffix}``
    """
    family = _validate_seahub_family(assay_family)
    valid_assays = ASSAYS_BY_FAMILY[family]
    patterns = _build_seahub_s3_patterns(family)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0
    group_assays: dict[str, set[str]] = defaultdict(set)
    run_ids: set[str] = set()

    for row in mappings:
        total += 1
        s3 = row.s3_path

        if _is_run_metadata(s3):
            metadata_count += 1
            continue

        if family == "scale" and "hash_oliga" in s3.lower():
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "assay 'hash_oliga' is not valid (did you mean 'hash_oligo'?)",
                }
            )

        m: re.Match[str] | None = None
        form = ""
        for pat, pat_form in patterns:
            m = pat.match(s3)
            if m:
                form = pat_form
                break

        if not m:
            warnings.append(
                {
                    "type": "parse_miss",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"S3 path does not match expected {family} raw pattern",
                }
            )
            continue

        matched += 1
        gd = m.groupdict()
        group_assays[gd["group_id"]].add(gd["assay"])
        run_ids.add(gd["runid"])

        literal_assay = gd["assay"]
        canonical = CANONICAL_ASSAY.get(literal_assay.lower())
        if canonical and literal_assay != canonical:
            errors.append(
                {
                    "type": "assay_casing",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"assay '{literal_assay}' violates SOP spelling '{canonical}'",
                }
            )

        assay = canonical or literal_assay.lower()
        if assay.lower() not in {a.lower() for a in valid_assays}:
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"assay '{literal_assay}' is not valid for {family} "
                    f"(expected one of {sorted(valid_assays)})",
                }
            )

        if gd["runid"] != gd["runid2"]:
            errors.append(
                {
                    "type": "runid_mismatch",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"runid mismatch between directory '{gd['runid']}' "
                    f"and filename '{gd['runid2']}'",
                }
            )

        project = gd["project"]
        if project != project.lower() or "_" in project:
            warnings.append(
                {
                    "type": "project_naming",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "project should be lower-case with hyphen delimiters (per SOP)",
                }
            )

        _seahub_post_match_checks(family, form, gd, assay, row, errors)

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "total": total,
        "matched": matched,
        "group_assays": dict(group_assays),
        "run_ids": run_ids,
    }


def validate_s3_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """Backward-compatible wrapper: Scale raw S3 validation."""
    return validate_s3_seahub_raw("scale", mappings)


def validate_s3_sci_raw(mappings: Iterable[MappingRow]) -> dict:
    """Backward-compatible wrapper: sci raw S3 validation."""
    return validate_s3_seahub_raw("sci", mappings)


# ---------------------------------------------------------------------------
# SIF completeness (shared)
# ---------------------------------------------------------------------------


def compare_groupid_assays(
    expected_group_assays: dict[str, set[str]],
    actual_group_assays: dict[str, set[str]],
) -> dict:
    """Compare expected vs actual GroupID -> assay-type mappings.

    Core comparison logic shared across all assay families.
    Both sides use lower-cased assay names.

    Returns a dict with:
    - expected_groupids, actual_groupids
    - missing_groupids, extra_groupids
    - missing_assays, extra_assays  (per-GroupID)
    - expected_group_assays, actual_group_assays  (pass-through)
    """
    expected_groupids = set(expected_group_assays.keys())
    actual_groupids = set(actual_group_assays.keys())

    missing_groupids = expected_groupids - actual_groupids
    extra_groupids = actual_groupids - expected_groupids

    missing_assays: dict[str, set[str]] = {}
    extra_assays: dict[str, set[str]] = {}

    for gid in expected_groupids & actual_groupids:
        expected_a = expected_group_assays.get(gid, set())
        actual_a = actual_group_assays.get(gid, set())
        miss = expected_a - actual_a
        extra = actual_a - expected_a
        if miss:
            missing_assays[gid] = miss
        if extra:
            extra_assays[gid] = extra

    return {
        "expected_groupids": expected_groupids,
        "actual_groupids": actual_groupids,
        "missing_groupids": missing_groupids,
        "extra_groupids": extra_groupids,
        "expected_group_assays": dict(expected_group_assays),
        "actual_group_assays": dict(actual_group_assays),
        "missing_assays": missing_assays,
        "extra_assays": extra_assays,
    }


def validate_sif_completeness_seahub(
    assay_family: str,
    mappings: Iterable[MappingRow],
    sif_path: str | Path,
) -> dict:
    """
    Compare seahub SIF GroupIDs and their assay types against S3 paths.

    Works for both Scale and sci families. Each SIF row represents one
    assay for a GroupID. This validator checks that every GroupID and
    assay type in the SIF is represented in the S3 paths and vice versa.
    """
    expected_group_assays = load_sif_group_assays(sif_path)
    patterns = _build_seahub_s3_patterns(assay_family)

    actual_group_assays: dict[str, set[str]] = defaultdict(set)
    for row in mappings:
        s3 = row.s3_path
        m: re.Match[str] | None = None
        for pat, _ in patterns:
            m = pat.match(s3)
            if m:
                break
        if not m:
            continue
        gd = m.groupdict()
        actual_group_assays[gd["group_id"]].add(gd["assay"].lower())

    return compare_groupid_assays(expected_group_assays, dict(actual_group_assays))


def validate_sif_completeness_scale(
    mappings: Iterable[MappingRow], sif_path: str | Path
) -> dict:
    """Backward-compatible wrapper: Scale SIF completeness."""
    return validate_sif_completeness_seahub("scale", mappings, sif_path)


def validate_local_paths_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate local paths for Novogene Scale/Quantum raw data.

    This focuses purely on local filesystem paths and checks for
    self-consistency of the QSR number, SCALEPLEX flag, and index /
    well-code within each path.
    """
    v1_pattern = re.compile(
        r"(?P<base_path>.+)/"
        r"(?P<wafer>\d+)-(?P<date>\d+_\d+)/"
        r"(?P<wafer2>\d+)-QSR-(?P<qsr_n>\d+)(?P<scaleplex>-SCALEPLEX)?-(?P<index>[ACGT]+)/"
        r"(?P<wafer3>\d+)-QSR-(?P<qsr_n2>\d+)(?P<scaleplex2>-SCALEPLEX)?-(?P<index2>[ACGT]+)"
        r"(?P<suffix>.*)"
    )

    v2_pattern = re.compile(
        r"(?P<base_path>.+)/"
        r"(?P<wafer>\d+)-(?P<date>\d+_\d+)/"
        r"(?P<dir_stem>(?P<wafer2>\d+)-QSR(?P<qsr_compact>\d+)(?P<sp_compact>SCALEPLEX)?"
        r"_QSR-(?P<qsr_dash>\d+)(?P<sp_dash>-SCALEPLEX)?)/"
        r"(?P<file_stem>(?P<wafer3>\d+)-QSR(?P<qsr_compact2>\d+)(?P<sp_compact2>SCALEPLEX)?"
        r"_QSR-(?P<qsr_dash2>\d+)(?P<sp_dash2>-SCALEPLEX)?)"
        r"_(?P<wellcode>[A-Za-z0-9]+)"
        r"(?P<suffix>.*)"
    )

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1
        local = row.local_path

        m = v1_pattern.match(local)
        if m:
            matched += 1
            gd = m.groupdict()

            if not (gd["wafer"] == gd["wafer2"] == gd["wafer3"]):
                errors.append(
                    {
                        "type": "wafer_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"wafer IDs differ within path: "
                        f"{gd['wafer']}, {gd['wafer2']}, {gd['wafer3']}",
                    }
                )
            if gd["qsr_n"] != gd["qsr_n2"]:
                errors.append(
                    {
                        "type": "qsr_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"QSR number mismatch between directory and filename: "
                        f"{gd['qsr_n']} vs {gd['qsr_n2']}",
                    }
                )
            if (gd["scaleplex"] is None) != (gd["scaleplex2"] is None):
                errors.append(
                    {
                        "type": "scaleplex_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": "SCALEPLEX present in one component of the path but not the other",
                    }
                )
            if gd["index"] != gd["index2"]:
                errors.append(
                    {
                        "type": "index_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"index sequence mismatch between directory and filename: "
                        f"{gd['index']} vs {gd['index2']}",
                    }
                )
            continue

        m = v2_pattern.match(local)
        if m:
            matched += 1
            gd = m.groupdict()

            if not (gd["wafer"] == gd["wafer2"] == gd["wafer3"]):
                errors.append(
                    {
                        "type": "wafer_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"wafer IDs differ within path: "
                        f"{gd['wafer']}, {gd['wafer2']}, {gd['wafer3']}",
                    }
                )

            if gd["qsr_compact"] != gd["qsr_dash"]:
                errors.append(
                    {
                        "type": "qsr_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"QSR number mismatch within directory: "
                        f"QSR{gd['qsr_compact']} vs QSR-{gd['qsr_dash']}",
                    }
                )

            if gd["qsr_compact2"] != gd["qsr_dash2"]:
                errors.append(
                    {
                        "type": "qsr_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"QSR number mismatch within filename: "
                        f"QSR{gd['qsr_compact2']} vs QSR-{gd['qsr_dash2']}",
                    }
                )

            if gd["dir_stem"] != gd["file_stem"]:
                errors.append(
                    {
                        "type": "dir_file_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"directory name '{gd['dir_stem']}' does not match "
                        f"filename prefix '{gd['file_stem']}'",
                    }
                )

            if bool(gd["sp_compact"]) != bool(gd["sp_dash"]):
                errors.append(
                    {
                        "type": "scaleplex_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": "SCALEPLEX present in one part of the directory name but not the other",
                    }
                )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
    }


def validate_s3_local_consistency_scale(mappings: Iterable[MappingRow]) -> dict:
    """
    Fuzzy consistency check between Scale S3 paths and local paths.

    This does not assume a precise local layout; instead it compares:
    - RunID (S3) vs wafer-like number (local)
    - QSR-N numbers seen in S3 vs those seen in local
    - SCALEPLEX presence in S3 vs local
    """
    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    qsr_pattern = re.compile(r"QSR-([0-9]+)", re.IGNORECASE)
    runid_pattern = re.compile(r"/raw/(\d+)/")
    local_wafer_pattern = re.compile(r"/(\d+)-\d+_\d+/")

    for row in mappings:
        total += 1
        s3 = row.s3_path
        local = row.local_path

        s3_lower = s3.lower()
        local_lower = local.lower()

        m_run = runid_pattern.search(s3)
        s3_runid = m_run.group(1) if m_run else None

        m_wafer = local_wafer_pattern.search(local)
        local_wafer = m_wafer.group(1) if m_wafer else None

        s3_qsr_nums = {int(n) for n in qsr_pattern.findall(s3)}
        local_qsr_nums = {int(n) for n in qsr_pattern.findall(local)}

        s3_scaleplex = "scaleplex" in s3_lower
        local_scaleplex = "scaleplex" in local_lower

        if any([s3_runid, local_wafer, s3_qsr_nums, local_qsr_nums, s3_scaleplex, local_scaleplex]):
            matched += 1
        else:
            continue

        if s3_runid and local_wafer and s3_runid != local_wafer:
            errors.append(
                {
                    "type": "runid_mismatch_s3_local",
                    "line": row.line_num,
                    "s3_path": s3,
                    "local_path": local,
                    "detail": f"S3 runid '{s3_runid}' does not match local wafer '{local_wafer}'",
                }
            )

        if s3_qsr_nums and local_qsr_nums:
            inter = s3_qsr_nums & local_qsr_nums
            if not inter:
                errors.append(
                    {
                        "type": "qsr_mismatch_s3_local",
                        "line": row.line_num,
                        "s3_path": s3,
                        "local_path": local,
                        "detail": f"QSR numbers disagree: S3 has {sorted(s3_qsr_nums)}, "
                        f"local has {sorted(local_qsr_nums)}",
                    }
                )
            elif s3_qsr_nums != local_qsr_nums:
                warnings.append(
                    {
                        "type": "qsr_partial_mismatch_s3_local",
                        "line": row.line_num,
                        "s3_path": s3,
                        "local_path": local,
                        "detail": f"QSR numbers overlap but differ: "
                        f"S3={sorted(s3_qsr_nums)}, local={sorted(local_qsr_nums)}",
                    }
                )
        elif s3_qsr_nums or local_qsr_nums:
            side = "S3" if s3_qsr_nums else "local"
            warnings.append(
                {
                    "type": "qsr_only_one_side",
                    "line": row.line_num,
                    "s3_path": s3,
                    "local_path": local,
                    "detail": f"QSR numbers present only on {side} side: "
                    f"{sorted(s3_qsr_nums or local_qsr_nums)}",
                }
            )

        if s3_scaleplex and not local_scaleplex:
            warnings.append(
                {
                    "type": "scaleplex_missing_in_local",
                    "line": row.line_num,
                    "s3_path": s3,
                    "local_path": local,
                    "detail": "S3 path includes SCALEPLEX but local path does not mention SCALEPLEX",
                }
            )
        elif local_scaleplex and not s3_scaleplex:
            warnings.append(
                {
                    "type": "scaleplex_missing_in_s3",
                    "line": row.line_num,
                    "s3_path": s3,
                    "local_path": local,
                    "detail": "Local path includes SCALEPLEX but S3 path does not mention SCALEPLEX",
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
    }


# ---------------------------------------------------------------------------
# sci-specific local path and S3/local consistency
# ---------------------------------------------------------------------------


def validate_local_paths_sci_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate local paths for Novogene sci raw data.

    Supports two local layouts:

    1. Double-UG (original):
        {base}/{RunID}-{date}_{time}/{RunID}-{GroupID}_{UG}-{UG}-{Barcode}/
            {RunID}-{GroupID}_{UG}-{UG}-{Barcode}{suffix}

    2. Single-UG (sci-plex style):
        {base}/{RunID}-{date}_{time}/{RunID}-{GroupID}-{UG}-{Barcode}/
            {RunID}-{GroupID}-{UG}-{Barcode}{suffix}

    GroupID may contain letters, digits, and underscores (e.g. CHEM16_P07_F3).
    Checks internal consistency of RunID, GroupID, UG, and Barcode between
    directory components and filename.
    """
    # Double-UG: GroupID_UG-UG-Barcode (underscore before UG, two UG codes)
    pattern_double_ug = re.compile(
        r"(?P<base>.+)/"
        r"(?P<runid>\d+)-(?P<date>\d+_\d+)/"
        r"(?P<runid2>\d+)-(?P<groupid>[A-Za-z0-9_]+)"
        r"_(?P<ug>Z\d{4})-(?P<ug2>Z\d{4})-(?P<barcode>[ACGT]+)/"
        r"(?P<runid3>\d+)-(?P<groupid2>[A-Za-z0-9_]+)"
        r"_(?P<ug3>Z\d{4})-(?P<ug4>Z\d{4})-(?P<barcode2>[ACGT]+)"
        r"(?P<suffix>.*)"
    )
    # Single-UG: GroupID-UG-Barcode (hyphen after GroupID, one UG)
    pattern_single_ug = re.compile(
        r"(?P<base>.+)/"
        r"(?P<runid>\d+)-(?P<date>\d+_\d+)/"
        r"(?P<runid2>\d+)-(?P<groupid>[A-Za-z0-9_]+)-(?P<ug>Z\d{4})-(?P<barcode>[ACGT]+)/"
        r"(?P<runid3>\d+)-(?P<groupid2>[A-Za-z0-9_]+)-(?P<ug3>Z\d{4})-(?P<barcode2>[ACGT]+)"
        r"(?P<suffix>.*)"
    )

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1
        local = row.local_path

        if _is_run_metadata(row.s3_path):
            continue

        m = pattern_double_ug.match(local)
        single_ug = False
        if not m:
            m = pattern_single_ug.match(local)
            if m:
                single_ug = True
        if not m:
            continue

        matched += 1
        gd = m.groupdict()

        if not (gd["runid"] == gd["runid2"] == gd["runid3"]):
            errors.append(
                {
                    "type": "runid_mismatch",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": f"RunID differs within path: {gd['runid']}, "
                    f"{gd['runid2']}, {gd['runid3']}",
                }
            )

        if gd["groupid"] != gd["groupid2"]:
            errors.append(
                {
                    "type": "groupid_mismatch",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": f"GroupID in directory '{gd['groupid']}' does not "
                    f"match filename '{gd['groupid2']}'",
                }
            )

        if not single_ug:
            if gd["ug"] != gd["ug2"]:
                errors.append(
                    {
                        "type": "ug_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"UG codes differ in directory: {gd['ug']} vs {gd['ug2']}",
                    }
                )

            if gd["ug3"] != gd["ug4"]:
                errors.append(
                    {
                        "type": "ug_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"UG codes differ in filename: {gd['ug3']} vs {gd['ug4']}",
                    }
                )

            if gd["ug"] != gd["ug3"]:
                errors.append(
                    {
                        "type": "ug_mismatch",
                        "line": row.line_num,
                        "local_path": local,
                        "detail": f"UG code in directory '{gd['ug']}' does not match "
                        f"filename '{gd['ug3']}'",
                    }
                )

        if gd["barcode"] != gd["barcode2"]:
            errors.append(
                {
                    "type": "barcode_mismatch",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": f"Barcode in directory '{gd['barcode']}' does not "
                    f"match filename '{gd['barcode2']}'",
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
    }


def validate_s3_local_consistency_sci(mappings: Iterable[MappingRow]) -> dict:
    """
    Cross-check sci S3 paths against local paths for consistency.

    Compares RunID, GroupID, UG code, and Barcode between the S3 path and
    the corresponding local path on each mapping row.
    """
    s3_patterns = _build_seahub_s3_patterns("sci")
    # Double-UG local: RunID-GroupID_UG-UG-Barcode
    local_regex_double_ug = re.compile(
        r"/(?P<runid>\d+)-\d+_\d+/"
        r"(?P<runid2>\d+)-(?P<groupid>[A-Za-z0-9_]+)"
        r"_(?P<ug>Z\d{4})-Z\d{4}-(?P<barcode>[ACGT]+)/"
    )
    # Single-UG local (sci-plex style): RunID-GroupID-UG-Barcode
    local_regex_single_ug = re.compile(
        r"/(?P<runid>\d+)-\d+_\d+/"
        r"(?P<runid2>\d+)-(?P<groupid>[A-Za-z0-9_]+)-(?P<ug>Z\d{4})-(?P<barcode>[ACGT]+)/"
    )

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1

        if _is_run_metadata(row.s3_path):
            continue

        s3_m: re.Match[str] | None = None
        for pat, _ in s3_patterns:
            s3_m = pat.match(row.s3_path)
            if s3_m:
                break
        local_m = local_regex_double_ug.search(row.local_path)
        if not local_m:
            local_m = local_regex_single_ug.search(row.local_path)

        if not s3_m or not local_m:
            continue

        matched += 1
        s3_gd = s3_m.groupdict()
        local_gd = local_m.groupdict()

        if s3_gd["runid"] != local_gd["runid"]:
            errors.append(
                {
                    "type": "runid_mismatch_s3_local",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "local_path": row.local_path,
                    "detail": f"S3 runid '{s3_gd['runid']}' does not match "
                    f"local runid '{local_gd['runid']}'",
                }
            )

        if s3_gd["group_id"] != local_gd["groupid"]:
            errors.append(
                {
                    "type": "groupid_mismatch_s3_local",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "local_path": row.local_path,
                    "detail": f"S3 groupid '{s3_gd['group_id']}' does not match "
                    f"local groupid '{local_gd['groupid']}'",
                }
            )

        if s3_gd["ug"] != local_gd["ug"]:
            errors.append(
                {
                    "type": "ug_mismatch_s3_local",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "local_path": row.local_path,
                    "detail": f"S3 UG '{s3_gd['ug']}' does not match "
                    f"local UG '{local_gd['ug']}'",
                }
            )

        if s3_gd["barcode"] != local_gd["barcode"]:
            errors.append(
                {
                    "type": "barcode_mismatch_s3_local",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "local_path": row.local_path,
                    "detail": f"S3 barcode '{s3_gd['barcode']}' does not match "
                    f"local barcode '{local_gd['barcode']}'",
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
    }


# ---------------------------------------------------------------------------
# Library-assay consistency (10x)
# ---------------------------------------------------------------------------


def validate_library_assay_consistency(
    mappings: Iterable[MappingRow],
    lib_assays: dict[str, str],
    provider: str,
) -> dict:
    """Cross-check library names in local paths against S3 paths."""
    provider = _validate_provider(provider)
    order_pattern = get_order_pattern(provider)
    valid_assays = get_assays("10x", provider)
    assay_re = build_assay_regex(valid_assays)

    s3_regex = re.compile(
        rf"^s3://czi-{re.escape(provider)}/[a-z0-9-]+/"
        rf"{order_pattern}/(?P<groupid>[^/]+)/raw/"
        rf"\d+-.+?_(?P<assay>{assay_re})-Z\d{{4}}-[ACGT]+"
    )

    libs_by_length = sorted(lib_assays.keys(), key=len, reverse=True)
    lib_patterns = {lib: f"-{lib}-" for lib in libs_by_length}

    assay_mismatches: List[dict] = []
    groupid_mismatches: List[dict] = []
    checked = 0
    skipped = 0

    for row in mappings:
        if _is_run_metadata(row.s3_path):
            continue
        m = s3_regex.match(row.s3_path)
        if not m:
            continue

        s3_assay = m.group("assay").lower()
        s3_groupid = m.group("groupid")

        found_lib: str | None = None
        for lib in libs_by_length:
            if lib_patterns[lib] in row.local_path:
                found_lib = lib
                break

        if found_lib is None:
            skipped += 1
            continue

        checked += 1

        if found_lib not in s3_groupid:
            groupid_mismatches.append(
                {
                    "line": row.line_num,
                    "library": found_lib,
                    "s3_groupid": s3_groupid,
                    "local_path": row.local_path,
                    "s3_path": row.s3_path,
                }
            )

        expected_assay = lib_assays[found_lib]
        if s3_assay != expected_assay:
            assay_mismatches.append(
                {
                    "line": row.line_num,
                    "library": found_lib,
                    "s3_assay": m.group("assay"),
                    "expected_assay": CANONICAL_ASSAY.get(expected_assay, expected_assay),
                    "local_path": row.local_path,
                    "s3_path": row.s3_path,
                }
            )

    return {
        "checked": checked,
        "assay_mismatches": assay_mismatches,
        "groupid_mismatches": groupid_mismatches,
        "skipped": skipped,
    }


# ---------------------------------------------------------------------------
# 10x processed S3 validation
# ---------------------------------------------------------------------------

_VALID_10X_PIPELINES: set[str] = {"cellranger"}
_RUN_DATE_RE: re.Pattern[str] = re.compile(r"^Run_\d{4}-\d{2}-\d{2}$")


def _build_10x_processed_s3_regex(
    provider: str, order_pattern: str
) -> re.Pattern[str]:
    """Build the regex for 10x processed S3 paths per the SOP.

    Expected layout::

        s3://czi-{provider}/{project}/{order}/{group_id}/processed/
            {pipeline}/{run_date}/outs/{file_path}
    """
    return re.compile(
        rf"^s3://(?P<bucket>czi-{re.escape(provider)})/"
        r"(?P<project>[^/]+)/"
        rf"(?P<order>{order_pattern})/"
        r"(?P<group_id>[^/]+)/processed/"
        r"(?P<pipeline>[^/]+)/"
        r"(?P<run_date>[^/]+)/"
        r"(?P<outs_marker>outs)/"
        r"(?P<file_path>.+)$"
    )


def validate_s3_10x_processed(
    provider: str, mappings: Iterable[MappingRow]
) -> dict:
    """Validate 10x processed S3 paths against the SOP.

    Expected S3 layout::

        s3://czi-{provider}/{project}/{order}/{GroupID}/processed/
            cellranger/{Run_YYYY-MM-DD}/outs/{file_path}
    """
    provider = _validate_provider(provider)
    order_pattern = get_order_pattern(provider)
    s3_regex = _build_10x_processed_s3_regex(provider, order_pattern)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    group_ids: set[str] = set()
    pipelines: set[str] = set()
    run_dates: set[str] = set()

    for row in mappings:
        total += 1
        s3 = row.s3_path

        m = s3_regex.match(s3)
        if not m:
            warnings.append(
                {
                    "type": "parse_miss",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "S3 path does not match expected 10x processed pattern",
                }
            )
            continue

        matched += 1
        gd = m.groupdict()
        group_ids.add(gd["group_id"])
        pipelines.add(gd["pipeline"])
        run_dates.add(gd["run_date"])

        pipeline = gd["pipeline"]
        if pipeline not in _VALID_10X_PIPELINES:
            errors.append(
                {
                    "type": "invalid_pipeline",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": (
                        f"pipeline '{pipeline}' is not valid for 10x processed "
                        f"(expected one of {sorted(_VALID_10X_PIPELINES)})"
                    ),
                }
            )

        run_date = gd["run_date"]
        if not _RUN_DATE_RE.match(run_date):
            warnings.append(
                {
                    "type": "run_date_format",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": (
                        f"run date '{run_date}' does not match expected "
                        "format Run_YYYY-MM-DD"
                    ),
                }
            )

        project = gd["project"]
        if project != project.lower() or "_" in project:
            warnings.append(
                {
                    "type": "project_naming",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": (
                        "project should be lower-case with hyphen "
                        "delimiters (per SOP)"
                    ),
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
        "group_ids": group_ids,
        "pipelines": pipelines,
        "run_dates": run_dates,
    }


# ---------------------------------------------------------------------------
# 10x processed SIF completeness
# ---------------------------------------------------------------------------


def validate_sif_completeness_10x_processed(
    provider: str,
    mappings: Iterable[MappingRow],
    sif_path: str,
) -> dict:
    """Compare GroupIDs found in processed S3 paths against SIF identifiers.

    Tries to load Group Identifiers from the SIF first.  If the SIF lacks
    a "Group Identifier" column (e.g. Ultima intake forms), falls back to
    the "Library name" column via :func:`load_sif_library_names`.
    """
    from .sif_io import _normalize_sif_groupid, load_sif_library_names

    s3_res = validate_s3_10x_processed(provider, mappings)
    s3_group_ids: set[str] = s3_res["group_ids"]

    sif_ga = load_sif_group_assays(sif_path)
    if sif_ga:
        sif_ids: set[str] = {_normalize_sif_groupid(k) for k in sif_ga}
    else:
        sif_ids = load_sif_library_names(sif_path)

    missing = sorted(sif_ids - s3_group_ids)
    extra = sorted(s3_group_ids - sif_ids)

    return {
        "sif_count": len(sif_ids),
        "s3_count": len(s3_group_ids),
        "missing": missing,
        "extra": extra,
    }


# ---------------------------------------------------------------------------
# 10x processed S3/local consistency
# ---------------------------------------------------------------------------


def validate_s3_local_consistency_10x_processed(
    provider: str, mappings: Iterable[MappingRow]
) -> dict:
    """Cross-check processed S3 and local paths for consistency.

    Checks:
    1. The GroupID from the S3 path appears in the local path.
    2. The file path after ``outs/`` matches between S3 and local.
    """
    provider = _validate_provider(provider)
    order_pattern = get_order_pattern(provider)
    s3_regex = _build_10x_processed_s3_regex(provider, order_pattern)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1
        s3 = row.s3_path
        local = row.local_path

        m = s3_regex.match(s3)
        if not m:
            continue

        matched += 1
        gd = m.groupdict()
        group_id = gd["group_id"]
        s3_file_path = gd["file_path"]

        if group_id not in local:
            errors.append(
                {
                    "type": "group_id_missing_local",
                    "line": row.line_num,
                    "s3_path": s3,
                    "local_path": local,
                    "detail": (
                        f"GroupID '{group_id}' from S3 not found in "
                        f"local path"
                    ),
                }
            )

        outs_idx = local.find("/outs/")
        if outs_idx >= 0:
            local_file_path = local[outs_idx + len("/outs/"):]
            if s3_file_path != local_file_path:
                errors.append(
                    {
                        "type": "file_path_mismatch",
                        "line": row.line_num,
                        "s3_path": s3,
                        "local_path": local,
                        "detail": (
                            f"file path after outs/ differs: "
                            f"S3='{s3_file_path}' vs local='{local_file_path}'"
                        ),
                    }
                )
        else:
            warnings.append(
                {
                    "type": "no_outs_in_local",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": "local path does not contain '/outs/' marker",
                }
            )

    return {
        "errors": errors,
        "warnings": warnings,
        "total": total,
        "matched": matched,
    }


__all__ = [
    "_is_run_metadata",
    "compare_groupid_assays",
    # 10x raw
    "validate_s3_10x_raw",
    "find_unmatched_sif_paths_10x",
    "validate_library_assay_consistency",
    # 10x processed
    "validate_s3_10x_processed",
    "validate_sif_completeness_10x_processed",
    "validate_s3_local_consistency_10x_processed",
    # Seahub (Scale + sci unified)
    "validate_s3_seahub_raw",
    "validate_sif_completeness_seahub",
    # Scale-specific (local path, S3/local consistency)
    "validate_local_paths_scale_raw",
    "validate_s3_local_consistency_scale",
    # sci-specific (local path, S3/local consistency)
    "validate_local_paths_sci_raw",
    "validate_s3_local_consistency_sci",
    # Backward-compatible wrappers
    "validate_s3_scale_raw",
    "validate_s3_sci_raw",
    "validate_sif_completeness_scale",
]
