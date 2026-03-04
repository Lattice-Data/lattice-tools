from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Tuple

from .constants import (
    ASSAYS_10X_NOVOGENE,
    ASSAYS_10X_PSOMAGEN,
    ASSAYS_SCALE,
    CANONICAL_ASSAY,
    _10X_NOVOGENE_ASSAY_RE,
    _10X_PSOMAGEN_ASSAY_RE,
    _NOVOGENE_ORDER_RE,
    _PSOMAGEN_ORDER_RE,
    _RUN_METADATA_RE,
    _SCALE_ASSAY_RE,
)
from .parsing import MappingRow
from .sif_io import load_sif_group_assays


def _is_run_metadata(s3_path: str) -> bool:
    """Return True if the S3 path's filename is a known run-level metadata file."""
    return bool(_RUN_METADATA_RE.match(os.path.basename(s3_path)))


def _build_10x_raw_s3_regex(provider: str, order_pattern: str, assay_re: str) -> re.Pattern[str]:
    """Build the canonical 10x raw S3 regex for a provider."""
    return re.compile(
        rf"^s3://(?P<bucket>czi-{re.escape(provider)})/"
        r"(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/"
        r"(?P<groupid>[^/]+)/raw/"
        rf"(?P<runid>\d+)-(?P<file_stem>.+?)_(?P<assay>{assay_re})-(?P<ug>Z\d{{4}})-(?P<barcode>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )


def validate_s3_10x_raw(provider: str, mappings: Iterable[MappingRow]) -> dict:
    """
    Validate 10x-style raw S3 paths against the SOP for a given provider.

    The expected S3 layout is:

        s3://czi-{provider}/{lastname}-{project}/{order}/{GroupID}/raw/
            {RunID}-{GroupID}_{Assay}-{UG-BC}{suffix}

    Where:
    - provider: 'novogene' or 'psomagen'
    - order:
        - novogene: NVUS#########-NN
        - psomagen: ANNNNN...
    - Assay:
        - novogene: GEX, CRI, ATAC
        - psomagen: GEX, CRI, ATAC, viral_ORF
    - UG-BC: Z####-BARCODE (BARCODE is A/C/G/T only)
    """
    provider = provider.lower()
    if provider not in {"novogene", "psomagen"}:
        raise ValueError(f"Unknown provider '{provider}', expected 'novogene' or 'psomagen'")

    if provider == "novogene":
        order_pattern = _NOVOGENE_ORDER_RE
        valid_assays = ASSAYS_10X_NOVOGENE
        assay_re = _10X_NOVOGENE_ASSAY_RE
    else:
        order_pattern = _PSOMAGEN_ORDER_RE
        valid_assays = ASSAYS_10X_PSOMAGEN
        assay_re = _10X_PSOMAGEN_ASSAY_RE

    s3_regex = _build_10x_raw_s3_regex(provider, order_pattern, assay_re)

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0
    group_assays: dict[str, set[str]] = defaultdict(set)

    for row in mappings:
        total += 1

        # Recognise run-level metadata files (not sample data)
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

        # GroupID consistency between path segment and filename prefix
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

        # UG format is enforced by the regex; double-check barcode is A/C/G/T only (already enforced)
        if not re.fullmatch(r"[ACGT]+", gd["barcode"]):
            errors.append(
                {
                    "type": "invalid_barcode",
                    "line": row.line_num,
                    "s3_path": row.s3_path,
                    "detail": f"barcode '{gd['barcode']}' contains characters outside A/C/G/T",
                }
            )

        # Project naming convention: lower-case and hyphen-delimited per SOP
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
    provider = provider.lower()
    if provider == "novogene":
        order_pattern = _NOVOGENE_ORDER_RE
        assay_re = _10X_NOVOGENE_ASSAY_RE
    else:
        order_pattern = _PSOMAGEN_ORDER_RE
        assay_re = _10X_PSOMAGEN_ASSAY_RE

    s3_regex = _build_10x_raw_s3_regex(provider, order_pattern, assay_re)

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

            # QSR number consistency within directory level
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

            # QSR number consistency within filename level
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

            # Directory stem must match filename stem
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

            # SCALEPLEX consistency within each level
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


def _build_scale_s3_patterns() -> Tuple[re.Pattern[str], re.Pattern[str]]:
    """Build and return the (sop_pattern, index_pattern) for Scale raw S3 paths."""
    sop = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{_NOVOGENE_ORDER_RE})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        rf"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>{_SCALE_ASSAY_RE})_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)"
        r"(?P<suffix>.*)$"
    )
    idx = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{_NOVOGENE_ORDER_RE})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        rf"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>{_SCALE_ASSAY_RE})_(?P<index_seq>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )
    return sop, idx


def validate_s3_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate Novogene Scale/Quantum raw S3 paths against the SOP.

    Two naming conventions are supported for the filename component:

    1. SOP / UG_RT form:
       {RunID}-{GroupID}_{Assay}_QSR-{N}[-SCALEPLEX]{suffix}

    2. Index-direct form:
       {RunID}-{GroupID}_{Assay}_{IndexSequence}{suffix}
    """
    sop_pattern, index_pattern = _build_scale_s3_patterns()

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0
    metadata_count = 0

    for row in mappings:
        total += 1
        s3 = row.s3_path

        # Recognise run-level metadata files (not sample data)
        if _is_run_metadata(s3):
            metadata_count += 1
            continue

        # Check for known typos early (before regex, since the typo prevents matching)
        if "hash_oliga" in s3.lower():
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "assay 'hash_oliga' is not valid (did you mean 'hash_oligo'?)",
                }
            )

        m = sop_pattern.match(s3)
        form = "sop"
        if not m:
            m = index_pattern.match(s3)
            form = "index"

        if not m:
            warnings.append(
                {
                    "type": "parse_miss",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "S3 path does not match expected Scale raw SOP or index-direct pattern",
                }
            )
            continue

        matched += 1
        gd = m.groupdict()
        literal_assay = gd["assay"]

        # The regex anchors guarantee the assay is a known name.
        # Check SOP canonical casing — any deviation is an SOP violation.
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

        # Normalise for downstream SOP rule checks
        assay = canonical or literal_assay.lower()

        # Validate assay is in the Scale-allowed set
        if assay.lower() not in {a.lower() for a in ASSAYS_SCALE}:
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"assay '{assay}' is not valid for Scale "
                    f"(expected one of {sorted(ASSAYS_SCALE)})",
                }
            )

        # RunID consistency between path segment and filename prefix
        if gd["runid"] != gd["runid2"]:
            errors.append(
                {
                    "type": "runid_mismatch",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": f"runid mismatch between path '{gd['runid']}' and filename '{gd['runid2']}'",
                }
            )

        # UG_RT-specific SOP checks
        if form == "sop":
            ug_rt = gd["ug_rt"]
            has_scaleplex = "SCALEPLEX" in ug_rt

            if assay in {"hash_oligo", "GEX_hash_oligo"} and not has_scaleplex:
                errors.append(
                    {
                        "type": "hash_oligo_scaleplex_mismatch",
                        "line": row.line_num,
                        "s3_path": s3,
                        "detail": f"assay '{assay}' requires SCALEPLEX in UG_RT, got '{ug_rt}'",
                    }
                )

            if assay == "GEX" and has_scaleplex:
                errors.append(
                    {
                        "type": "gex_scaleplex_violation",
                        "line": row.line_num,
                        "s3_path": s3,
                        "detail": f"assay 'GEX' must not include SCALEPLEX in UG_RT, got '{ug_rt}'",
                    }
                )

    return {
        "errors": errors,
        "warnings": warnings,
        "metadata_files": metadata_count,
        "total": total,
        "matched": matched,
    }


def validate_sif_completeness_scale(
    mappings: Iterable[MappingRow], sif_path: str | Path
) -> dict:
    """
    Compare Scale SIF GroupIDs and their assay types against S3 paths.

    Each SIF row represents one assay for a GroupID (e.g. R112A/GEX and
    R112A/Hash_oligo are two rows sharing the same GroupID).  This
    validator checks that:
    - Every GroupID in the SIF has at least one S3 path.
    - Every assay type listed in the SIF for a GroupID is represented
      in the S3 paths for that GroupID.
    - No unexpected GroupIDs or assay types appear in S3.
    """
    expected_group_assays = load_sif_group_assays(sif_path)
    expected_groupids = set(expected_group_assays.keys())

    sop_pattern, index_pattern = _build_scale_s3_patterns()

    actual_group_assays: dict[str, set[str]] = defaultdict(set)
    for row in mappings:
        s3 = row.s3_path
        m = sop_pattern.match(s3) or index_pattern.match(s3)
        if not m:
            continue
        gd = m.groupdict()
        group_id = gd["group_id"]
        assay = gd["assay"].lower()
        actual_group_assays[group_id].add(assay)

    actual_groupids = set(actual_group_assays.keys())
    missing_groupids = expected_groupids - actual_groupids
    extra_groupids = actual_groupids - expected_groupids

    # Per-GroupID assay completeness
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

        # Extract run ID from S3 and wafer-like ID from local
        m_run = runid_pattern.search(s3)
        s3_runid = m_run.group(1) if m_run else None

        m_wafer = local_wafer_pattern.search(local)
        local_wafer = m_wafer.group(1) if m_wafer else None

        # Extract all QSR numbers
        s3_qsr_nums = {int(n) for n in qsr_pattern.findall(s3)}
        local_qsr_nums = {int(n) for n in qsr_pattern.findall(local)}

        s3_scaleplex = "scaleplex" in s3_lower
        local_scaleplex = "scaleplex" in local_lower

        # Decide whether this pair is "matched" enough to reason about
        if any([s3_runid, local_wafer, s3_qsr_nums, local_qsr_nums, s3_scaleplex, local_scaleplex]):
            matched += 1
        else:
            continue

        # RunID vs wafer
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

        # QSR numbers
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
            # Only one side has QSR numbers – fuzzy warning.
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

        # SCALEPLEX presence
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


def validate_library_assay_consistency(
    mappings: Iterable[MappingRow],
    lib_assays: dict[str, str],
    provider: str,
) -> dict:
    """Cross-check library names in local paths against S3 paths."""
    provider = provider.lower()
    if provider == "novogene":
        order_pattern = _NOVOGENE_ORDER_RE
        assay_re = _10X_NOVOGENE_ASSAY_RE
    else:
        order_pattern = _PSOMAGEN_ORDER_RE
        assay_re = _10X_PSOMAGEN_ASSAY_RE

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

        # Check 1: library name must appear in the S3 GroupID
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

        # Check 2: assay in S3 must match SIF expectation for this library
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


__all__ = [
    "_is_run_metadata",
    "validate_s3_10x_raw",
    "find_unmatched_sif_paths_10x",
    "validate_local_paths_scale_raw",
    "validate_s3_scale_raw",
    "validate_sif_completeness_scale",
    "validate_s3_local_consistency_scale",
    "validate_library_assay_consistency",
]

