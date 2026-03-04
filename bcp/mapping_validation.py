"""
Mapping validation utilities for BCP workflows.

This module provides reusable primitives and a small CLI for validating
mapping CSV/TSV files that describe S3 ↔ local file path mappings.

Core helpers:
- parse_mapping_file: load a 2-column mapping (S3 path, local path)
- validate_uniqueness: ensure 1:1 mapping between S3 and local paths
- validate_s3_10x_raw: SOP checks for 10x-style raw S3 paths
- validate_local_paths_scale_raw / validate_s3_scale_raw: SOP checks for
  Scale/Quantum raw data (local and S3 sides)
- validate_sif_completeness_scale: SIF-vs-S3 completeness for Scale GroupIDs

CLI entrypoint (examples):

    python -m mapping_validation \\
        --mapping /path/to/mapping.csv \\
        --provider novogene \\
        --data raw --assay 10x

    python -m mapping_validation \\
        --mapping /path/to/mapping.csv \\
        --sif /path/to/SIF.csv \\
        --provider novogene \\
        --data raw --assay scale
"""

from __future__ import annotations

import os
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple
import argparse
import re
import sys

warnings.filterwarnings(
    "ignore",
    message="Data Validation extension is not supported",
    category=UserWarning,
    module="openpyxl",
)
import pandas as pd

# ---------------------------------------------------------------------------
# SOP-defined assay types per assay family
# ---------------------------------------------------------------------------
ASSAYS_10X_NOVOGENE: set[str] = {"GEX", "CRI", "ATAC"}
ASSAYS_10X_PSOMAGEN: set[str] = {"GEX", "CRI", "ATAC", "viral_ORF"}
ASSAYS_SCALE: set[str] = {"GEX", "hash_oligo", "GEX_hash_oligo"}
ASSAYS_SCI: set[str] = {"GEX", "hash_oligo", "GEX_hash_oligo"}

# Canonical (SOP) spellings – lower-cased key → canonical form
CANONICAL_ASSAY: dict[str, str] = {
    "gex": "GEX",
    "cri": "CRI",
    "atac": "ATAC",
    "hash_oligo": "hash_oligo",
    "gex_hash_oligo": "GEX_hash_oligo",
    "viral_orf": "viral_ORF",
}

# Regex alternation fragments used to anchor assay names inside S3 path
# regexes. Longest alternatives come first to avoid partial matches.
# Case-insensitive matching via (?i:...) so that ANY casing variant is
# captured (and then checked against the canonical SOP spelling).
_SCALE_ASSAY_RE = r"(?i:GEX_hash_oligo|hash_oligo|GEX)"
_10X_NOVOGENE_ASSAY_RE = r"(?i:ATAC|CRI|GEX)"
_10X_PSOMAGEN_ASSAY_RE = r"(?i:viral_ORF|ATAC|CRI|GEX)"

# Order number patterns (Novogene allows sub-order suffixes like -28-4)
_NOVOGENE_ORDER_RE = r"NVUS\d+-\d+(?:-\d+)*"
_PSOMAGEN_ORDER_RE = r"AN\d+"

# Run-level metadata filenames expected alongside sample data.
# The optional \d+_ prefix handles files placed at the order level
# with a RunID prefix (e.g. 439844_UploadCompleted.json).
_RUN_METADATA_RE = re.compile(
    r"^(?:\d+_)?(?:"
    r"LibraryInfo\.xml"
    r"|SequencingInfo\.json"
    r"|UploadCompleted\.json"
    r"|merged_trimmer-(?:failure_codes|stats)\.csv"
    r"|run_(?:SecondaryAnalysis|VariantCalling)\.txt"
    r")$"
)


def _is_run_metadata(s3_path: str) -> bool:
    """Return True if the S3 path's filename is a known run-level metadata file."""
    return bool(_RUN_METADATA_RE.match(os.path.basename(s3_path)))


def _normalize_sif_groupid(gid: str) -> str:
    """Normalise a SIF Group Identifier to match the S3 directory convention.

    SIF files for 10x use ``A + AF`` style (space-plus-space) while S3
    paths join the same parts with underscores: ``A_AF``.  Scale SIFs
    already use plain identifiers, so this is a no-op for them.
    """
    return gid.replace(" + ", "_")


@dataclass(frozen=True)
class MappingRow:
    """Single row from a mapping file."""

    s3_path: str
    local_path: str
    line_num: int


def _split_mapping_line(line: str) -> Tuple[str, str] | None:
    """
    Split a raw line from a mapping file into (s3_path, local_path).

    Supports the common formats in existing mapping files:
    - Comma-separated:  s3://...,/ORPROJ1/...
    - Tab-separated:    s3://...\t/ORPROJ1/...
    - Fallback: first “s3://...” token and the rest treated as local path.
    """
    stripped = line.strip()
    if not stripped:
        return None

    # Preferred: CSV-style with a single comma separating S3 and local
    if "," in stripped:
        s3, local = stripped.split(",", 1)
        return s3.strip(), local.strip()

    # TSV-style
    if "\t" in stripped:
        parts = stripped.split("\t")
        if len(parts) >= 2:
            return parts[0].strip(), parts[1].strip()

    # Generic “s3:// ... <whitespace> /path” pattern
    if stripped.startswith("s3://"):
        # Find the first space that precedes an absolute local path
        # e.g. "s3://...   /ORPROJ1/..." or "s3://... /mnt/..."
        for idx in range(len(stripped)):
            if stripped[idx] == " " and idx + 1 < len(stripped) and stripped[idx + 1] == "/":
                s3 = stripped[:idx]
                local = stripped[idx + 1 :]
                return s3.strip(), local.strip()

    return None


def parse_mapping_file(path: str | Path) -> List[MappingRow]:
    """
    Parse a mapping file into a list of MappingRow objects.

    Each non-empty line is expected to contain exactly two columns:
    an S3 path and a local filesystem path. Lines that cannot be
    parsed into two columns are skipped, but their line numbers are
    preserved for diagnostics in downstream validators.
    """
    rows: List[MappingRow] = []
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        for line_num, raw in enumerate(fh, 1):
            parts = _split_mapping_line(raw)
            if parts is None:
                continue
            s3, local = parts
            if not s3 or not local:
                continue
            rows.append(MappingRow(s3_path=s3, local_path=local, line_num=line_num))
    return rows


def validate_uniqueness(mappings: Iterable[MappingRow]) -> dict:
    """
    Check that each local path maps to exactly one S3 path and vice versa.

    Returns a dictionary with:
    - duplicate_locals: {local_path: [MappingRow, ...]}
    - duplicate_s3: {s3_path: [MappingRow, ...]}
    - total: total number of mappings inspected
    """
    local_to_rows: dict[str, List[MappingRow]] = defaultdict(list)
    s3_to_rows: dict[str, List[MappingRow]] = defaultdict(list)

    total = 0
    for row in mappings:
        total += 1
        local_to_rows[row.local_path].append(row)
        s3_to_rows[row.s3_path].append(row)

    duplicate_locals = {k: v for k, v in local_to_rows.items() if len(v) > 1}
    duplicate_s3 = {k: v for k, v in s3_to_rows.items() if len(v) > 1}

    return {
        "duplicate_locals": duplicate_locals,
        "duplicate_s3": duplicate_s3,
        "total": total,
    }


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

    Returns a dict with:
    - errors: list of issue dicts (hard violations)
    - warnings: list of issue dicts (soft SOP deviations)
    - total: number of mappings inspected
    - matched: number of S3 paths that matched the 10x raw pattern
    """
    provider = provider.lower()
    if provider not in {"novogene", "psomagen"}:
        raise ValueError(f"Unknown provider '{provider}', expected 'novogene' or 'psomagen'")

    if provider == "novogene":
        order_pattern = _NOVOGENE_ORDER_RE
        valid_assays = ASSAYS_10X_NOVOGENE
    else:
        order_pattern = _PSOMAGEN_ORDER_RE
        valid_assays = ASSAYS_10X_PSOMAGEN

    s3_regex = re.compile(
        rf"^s3://(?P<bucket>czi-{re.escape(provider)})/"
        r"(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/"
        r"(?P<groupid>[^/]+)/raw/"
        # After raw/ we expect <RunID>-<something>_<Assay>-Z####-BARCODE...
        r"(?P<runid>\d+)-(?P<file_stem>.+?)_(?P<assay>[A-Za-z_]+)-(?P<ug>Z\d{4})-(?P<barcode>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )

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


def validate_local_paths_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate local paths for Novogene Scale/Quantum raw data.

    This focuses purely on local filesystem paths and checks for
    self-consistency of the QSR number, SCALEPLEX flag, and index /
    well-code within each path.

    Two local layouts are recognised:

    V1 (index-direct):
        {base}/{wafer}-{date}/
            {wafer}-QSR-{N}[-SCALEPLEX]-{INDEX}/
                {wafer}-QSR-{N}[-SCALEPLEX]-{INDEX}{suffix}

    V2 (compact QSR):
        {base}/{wafer}-{date}/
            {wafer}-QSR{N}[SCALEPLEX]_QSR-{N}[-SCALEPLEX]/
                {wafer}-QSR{N}[SCALEPLEX]_QSR-{N}[-SCALEPLEX]_{WellCode}{suffix}

    Returns a dict with:
    - errors: list of issue dicts
    - warnings: list of issue dicts
    - total: number of mappings inspected
    - matched: number of local paths that matched either pattern
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


def validate_s3_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate Novogene Scale/Quantum raw S3 paths against the SOP.

    Two naming conventions are supported for the filename component:

    1. SOP / UG_RT form:
       {RunID}-{GroupID}_{Assay}_QSR-{N}[-SCALEPLEX]{suffix}

    2. Index-direct form:
       {RunID}-{GroupID}_{Assay}_{IndexSequence}{suffix}

    Where:
    - Assay is one of: GEX, hash_oligo, GEX_hash_oligo (per SOP).
    - hash_oliga (typo) is explicitly flagged as an error.
    - For UG_RT form:
        * hash_oligo and GEX_hash_oligo must include SCALEPLEX
        * GEX must NOT include SCALEPLEX

    Returns a dict with:
    - errors: list of issue dicts
    - warnings: list of issue dicts
    - metadata_files: count of recognized run-level metadata files
    - total: number of mappings inspected
    - matched: number of S3 paths that matched either pattern
    """
    sop_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{_NOVOGENE_ORDER_RE})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        rf"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>{_SCALE_ASSAY_RE})_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)"
        r"(?P<suffix>.*)$"
    )

    index_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{_NOVOGENE_ORDER_RE})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        rf"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>{_SCALE_ASSAY_RE})_(?P<index_seq>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )

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


def load_sif_group_assays(sif_path: str | Path) -> dict[str, set[str]]:
    """
    Load expected GroupID → {assay types} mapping from a SIF file.

    Each row in the SIF represents a sublibrary for a group identifier.
    The same group identifier appears once per assay type (e.g. GEX and
    CRI for 10x, or GEX and hash_oligo for Scale).  Returns a dict
    mapping each GroupID (as it appears in the SIF) to its set of assay
    types (normalised to lower-case).

    Falls back through three parsing strategies:
    1. Excel (.xlsx/.xlsm/.xls) with pandas
    2. Well-formed CSV with DictReader
    3. Heuristic Ultima intake-style CSV (Scale-specific)
    """
    import csv

    p = Path(sif_path)
    result: dict[str, set[str]] = defaultdict(set)

    def _find_col(cols_lower: dict[str, str], keyword: str) -> str | None:
        for key, name in cols_lower.items():
            if keyword in key:
                return name
        return None

    # Branch 1: Excel SIF
    if p.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        try:
            df = pd.read_excel(p)
        except Exception:
            df = None

        if df is not None and not df.empty:
            cols = {str(c).strip().lower(): c for c in df.columns}
            group_col = _find_col(cols, "group identifier")
            assay_col = _find_col(cols, "assay type")

            if group_col is not None:
                for _, row_data in df.iterrows():
                    gid = str(row_data.get(group_col, "")).strip()
                    if not gid or gid == "nan":
                        continue
                    assay_val = str(row_data.get(assay_col, "")).strip() if assay_col else ""
                    if assay_val and assay_val != "nan":
                        result[gid].add(assay_val.lower())
                    else:
                        result.setdefault(gid, set())

        if result:
            return dict(result)

    # Branch 2: well-formed CSV
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames:
            field_map = {name.lower(): name for name in reader.fieldnames}
            group_col_name = _find_col(field_map, "group identifier")
            assay_col_name = _find_col(field_map, "assay type")

            if group_col_name is not None:
                for row in reader:
                    gid = (row.get(group_col_name) or "").strip()
                    if not gid:
                        continue
                    assay_val = (row.get(assay_col_name) or "").strip() if assay_col_name else ""
                    if assay_val:
                        result[gid].add(assay_val.lower())
                    else:
                        result.setdefault(gid, set())

    if result:
        return dict(result)

    # Branch 3: Ultima intake-style heuristic (col 5=GroupID, col 6=Assay)
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith('"Library name'):
                continue
            if not re.match(r"^[A-Za-z0-9]+,QSR", stripped):
                continue
            parts = [c.strip() for c in stripped.split(",")]
            if len(parts) >= 7:
                gid = parts[5]
                assay_val = parts[6]
                if gid and assay_val:
                    result[gid].add(assay_val.lower())
                elif gid:
                    result.setdefault(gid, set())
            elif len(parts) >= 6:
                gid = parts[5]
                if gid:
                    result.setdefault(gid, set())

    return dict(result)


load_sif_scale_group_assays = load_sif_group_assays


def load_sif_scale_groupids(sif_path: str | Path) -> set[str]:
    """
    Load expected GroupIDs for Scale libraries from a SIF file.

    Convenience wrapper around :func:`load_sif_group_assays` that
    returns just the set of GroupIDs.
    """
    return set(load_sif_group_assays(sif_path).keys())


def load_sif_library_assays(sif_path: str | Path) -> dict[str, str]:
    """Load a Library-Name → assay-type mapping from a SIF file.

    Returns a dict mapping each individual library name to its assay type
    (lower-cased).  For 10x data the SIF typically has one row per
    library, so ``FTF1732A`` → ``gex`` and ``FTF1732AF`` → ``cri``.
    """
    p = Path(sif_path)
    result: dict[str, str] = {}

    def _find_col(cols_lower: dict[str, str], keyword: str) -> str | None:
        for key, name in cols_lower.items():
            if keyword in key:
                return name
        return None

    if p.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        try:
            df = pd.read_excel(p)
        except Exception:
            df = None
        if df is not None and not df.empty:
            cols = {str(c).strip().lower(): c for c in df.columns}
            lib_col = _find_col(cols, "library name")
            assay_col = _find_col(cols, "assay type")
            if lib_col is not None and assay_col is not None:
                for _, row_data in df.iterrows():
                    lib = str(row_data.get(lib_col, "")).strip()
                    assay = str(row_data.get(assay_col, "")).strip()
                    if lib and lib != "nan" and assay and assay != "nan":
                        result[lib] = assay.lower()
        if result:
            return result

    import csv
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames:
            field_map = {name.lower(): name for name in reader.fieldnames}
            lib_col_name = _find_col(field_map, "library name")
            assay_col_name = _find_col(field_map, "assay type")
            if lib_col_name is not None and assay_col_name is not None:
                for row in reader:
                    lib = (row.get(lib_col_name) or "").strip()
                    assay = (row.get(assay_col_name) or "").strip()
                    if lib and assay:
                        result[lib] = assay.lower()

    return result


def validate_library_assay_consistency(
    mappings: Iterable[MappingRow],
    lib_assays: dict[str, str],
    provider: str,
) -> dict:
    """Cross-check library names in local paths against S3 paths.

    For each mapping row the function:
    1. Parses the S3 path to extract the GroupID directory and assay.
    2. Searches the local path for a SIF library name (longest-first
       to avoid partial matches like ``FTF1732J`` inside ``FTF1732JF``).
    3. Checks that the library name appears in the S3 GroupID (e.g.
       ``FTF1732G`` must be part of ``FTF1732G_FTF1732GF``).
    4. Compares the S3 assay with the expected assay from the SIF.

    Returns a dict with ``checked``, ``assay_mismatches``,
    ``groupid_mismatches``, and ``skipped``.
    """
    provider = provider.lower()
    if provider == "novogene":
        order_pattern = _NOVOGENE_ORDER_RE
    else:
        order_pattern = _PSOMAGEN_ORDER_RE

    s3_regex = re.compile(
        rf"^s3://czi-{re.escape(provider)}/[a-z0-9-]+/"
        rf"{order_pattern}/(?P<groupid>[^/]+)/raw/"
        r"\d+-.+?_(?P<assay>[A-Za-z_]+)-Z\d{4}-[ACGT]+"
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
            groupid_mismatches.append({
                "line": row.line_num,
                "library": found_lib,
                "s3_groupid": s3_groupid,
                "local_path": row.local_path,
                "s3_path": row.s3_path,
            })

        # Check 2: assay in S3 must match SIF expectation for this library
        expected_assay = lib_assays[found_lib]
        if s3_assay != expected_assay:
            assay_mismatches.append({
                "line": row.line_num,
                "library": found_lib,
                "s3_assay": m.group("assay"),
                "expected_assay": CANONICAL_ASSAY.get(expected_assay, expected_assay),
                "local_path": row.local_path,
                "s3_path": row.s3_path,
            })

    return {
        "checked": checked,
        "assay_mismatches": assay_mismatches,
        "groupid_mismatches": groupid_mismatches,
        "skipped": skipped,
    }


def _build_scale_s3_patterns() -> Tuple[re.Pattern, re.Pattern]:
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

    Rules:
    - If both S3 runid and local wafer are found and differ → error.
    - If both S3 and local have QSR numbers and their sets are disjoint → error.
    - If QSR sets overlap but are not identical → warning.
    - If SCALEPLEX appears on one side but not the other → warning.
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


def _print_issue_examples(label: str, issues: List[dict], kind: str, max_examples: int = 5) -> None:
    """
    Print a short, grouped summary of errors or warnings with a few examples.

    `kind` is either 'errors' or 'warnings' for labeling.
    """
    if not issues:
        return

    by_type: dict[str, int] = defaultdict(int)
    for item in issues:
        by_type[item.get("type", "unknown")] += 1

    print(f"{label} {kind}:")
    for t, count in sorted(by_type.items(), key=lambda kv: kv[0]):
        print(f"  - {t}: {count}")

    print(f"  Showing up to {max_examples} example {kind}:")
    for item in issues[:max_examples]:
        line = item.get("line", "?")
        s3_path = item.get("s3_path")
        local_path = item.get("local_path")
        detail = item.get("detail")
        print(f"    line {line}:")
        if s3_path:
            print(f"      S3:    {s3_path}")
        if local_path:
            print(f"      local: {local_path}")
        if detail:
            print(f"      detail: {detail}")


def main() -> None:
    """
    CLI for validating mapping CSV/TSV files against SOP rules.

    Required arguments:
        --mapping PATH       Mapping file with two columns: S3, local
        --provider {novogene,psomagen}
        --data {raw,processed}
        --assay {10x,sci,scale}

    For Scale raw validation, a SIF file is strongly recommended:
        --sif PATH
    """

    parser = argparse.ArgumentParser(description="Validate S3/local mapping files against SOP rules.")
    parser.add_argument("--mapping", required=True, help="Path to mapping CSV/TSV file")
    parser.add_argument("--sif", help="Path to SIF CSV/XLSX for completeness checks (Scale)")
    parser.add_argument(
        "--provider",
        required=True,
        choices=["novogene", "psomagen"],
        help="Sequencing provider (affects SOP rules)",
    )
    parser.add_argument(
        "--data",
        required=True,
        choices=["raw", "processed"],
        help="Which kind of data the mapping represents",
    )
    parser.add_argument(
        "--assay",
        required=True,
        choices=["10x", "sci", "scale"],
        help="High-level assay family for SOP selection",
    )

    args = parser.parse_args()

    mappings = parse_mapping_file(args.mapping)
    if not mappings:
        print("No mappings parsed from file (check delimiter and content).", file=sys.stderr)
        raise SystemExit(1)

    exit_code = 0
    fail_reasons: List[str] = []

    # 1. Uniqueness
    uniq = validate_uniqueness(mappings)
    dup_local_count = len(uniq["duplicate_locals"])
    dup_s3_count = len(uniq["duplicate_s3"])
    print(f"Uniqueness: {uniq['total']} mappings, "
          f"{dup_local_count} duplicate locals, {dup_s3_count} duplicate S3 paths")
    if dup_local_count or dup_s3_count:
        exit_code = 1
        fail_reasons.append(f"duplicate mappings ({dup_local_count} local, {dup_s3_count} S3)")

    # 2. Mode-specific validation
    provider = args.provider
    if args.data == "raw" and args.assay == "10x":
        res = validate_s3_10x_raw(provider, mappings)
        meta_count = res.get("metadata_files", 0)
        s3_ga: dict[str, set[str]] = res.get("group_assays", {})
        print(
            f"10x raw SOP: matched {res['matched']} S3 paths, "
            f"{len(res['errors'])} errors, {len(res['warnings'])} warnings"
            + (f", {meta_count} run-metadata files" if meta_count else "")
        )
        _print_issue_examples("10x raw SOP", res["errors"], "errors")
        _print_issue_examples("10x raw SOP", res["warnings"], "warnings")
        if res["matched"] == 0:
            print(
                "  WARNING: none of the S3 paths matched the 10x raw pattern. "
                "Check that the S3 paths follow the SOP naming convention."
            )
            fail_reasons.append("no S3 paths matched the 10x raw pattern (0 matched)")
            exit_code = 1
        if res["errors"]:
            exit_code = 1
            fail_reasons.append(f"{len(res['errors'])} 10x S3 SOP errors")

        # --- sanity printouts: GroupIDs and assays from S3 ---
        if s3_ga:
            print(f"S3 GroupIDs found: {len(s3_ga)}")
            for gid in sorted(s3_ga):
                assays_str = ", ".join(sorted(s3_ga[gid]))
                print(f"  {gid}: {assays_str}")

        # --- SIF comparison (when provided) ---
        if args.sif:
            sif_ga = load_sif_group_assays(args.sif)
            sif_norm: dict[str, set[str]] = {
                _normalize_sif_groupid(k): v for k, v in sif_ga.items()
            }
            sif_ids = set(sif_norm.keys())
            s3_ids = set(s3_ga.keys())
            missing_ids = sif_ids - s3_ids
            extra_ids = s3_ids - sif_ids

            print(
                f"SIF completeness: expected={len(sif_ids)} GroupIDs, "
                f"found={len(s3_ids)} in S3, "
                f"missing={len(missing_ids)}, extra={len(extra_ids)}"
            )

            if missing_ids:
                print("  Missing GroupIDs from S3 (in SIF but not S3): "
                      + ", ".join(sorted(missing_ids)))
                exit_code = 1
                fail_reasons.append(
                    f"{len(missing_ids)} SIF GroupIDs missing from S3: "
                    + ", ".join(sorted(missing_ids))
                )
            if extra_ids:
                print("  Extra GroupIDs in S3 (not in SIF): "
                      + ", ".join(sorted(extra_ids)))
                exit_code = 1
                fail_reasons.append(
                    f"{len(extra_ids)} extra GroupIDs in S3 not in SIF: "
                    + ", ".join(sorted(extra_ids))
                )

            # Per-GroupID assay comparison
            missing_assays: dict[str, set[str]] = {}
            extra_assays: dict[str, set[str]] = {}
            for gid in sif_ids & s3_ids:
                expected = sif_norm[gid]
                actual = {a.lower() for a in s3_ga.get(gid, set())}
                m_diff = expected - actual
                e_diff = actual - expected
                if m_diff:
                    missing_assays[gid] = m_diff
                if e_diff:
                    extra_assays[gid] = e_diff

            if missing_assays:
                print("  GroupIDs with missing assay types in S3:")
                for gid in sorted(missing_assays):
                    print(f"    {gid}: missing {sorted(missing_assays[gid])}")
                exit_code = 1
                fail_reasons.append(
                    f"{len(missing_assays)} GroupIDs have missing assay types in S3"
                )
            if extra_assays:
                print("  GroupIDs with unexpected assay types in S3 (not in SIF):")
                for gid in sorted(extra_assays):
                    print(f"    {gid}: unexpected {sorted(extra_assays[gid])}")

            # Summary of expected SIF structure
            unique_combos: set[frozenset[str]] = set()
            for assays in sif_norm.values():
                unique_combos.add(frozenset(assays))
            combos_desc = [
                " + ".join(CANONICAL_ASSAY.get(a, a) for a in sorted(c))
                for c in sorted(unique_combos, key=sorted)
            ]
            print(f"  SIF assay combinations: {', '.join(combos_desc)}")

            # Library-level cross-checks (assay + GroupID consistency)
            lib_assays = load_sif_library_assays(args.sif)
            if lib_assays:
                lib_res = validate_library_assay_consistency(
                    mappings, lib_assays, provider
                )
                n_assay = len(lib_res["assay_mismatches"])
                n_gid = len(lib_res["groupid_mismatches"])
                print(
                    f"Library consistency: checked {lib_res['checked']} paths, "
                    f"{n_assay} assay mismatches, {n_gid} GroupID mismatches"
                    + (f", {lib_res['skipped']} skipped (no library name found)"
                       if lib_res["skipped"] else "")
                )
                if n_gid:
                    print("  GroupID mismatches (library in local path not in S3 GroupID):")
                    for item in lib_res["groupid_mismatches"][:10]:
                        print(
                            f"    line {item['line']}: local has '{item['library']}' "
                            f"but S3 GroupID is '{item['s3_groupid']}'"
                        )
                    if n_gid > 10:
                        print(f"    ... and {n_gid - 10} more")
                    exit_code = 1
                    fail_reasons.append(
                        f"{n_gid} GroupID mismatches "
                        "(local library name not found in S3 GroupID)"
                    )
                if n_assay:
                    print("  Assay mismatches (library in local path vs assay in S3):")
                    for item in lib_res["assay_mismatches"][:10]:
                        print(
                            f"    line {item['line']}: local has '{item['library']}' "
                            f"(SIF expects {item['expected_assay']}) "
                            f"but S3 assay is {item['s3_assay']}"
                        )
                    if n_assay > 10:
                        print(f"    ... and {n_assay - 10} more")
                    exit_code = 1
                    fail_reasons.append(
                        f"{n_assay} library-assay mismatches "
                        "(local library name doesn't match S3 assay per SIF)"
                    )
        else:
            print("10x mode: no --sif provided, skipping SIF completeness checks.")

    elif args.data == "raw" and args.assay == "scale" and provider == "novogene":
        # Local path sanity
        local_res = validate_local_paths_scale_raw(mappings)
        print(f"Scale raw (local): matched {local_res['matched']} paths, "
              f"{len(local_res['errors'])} errors, {len(local_res['warnings'])} warnings")
        _print_issue_examples("Scale raw (local)", local_res["errors"], "errors")
        _print_issue_examples("Scale raw (local)", local_res["warnings"], "warnings")
        if local_res["matched"] == 0:
            print(
                "  WARNING: none of the local paths matched any recognised Scale pattern. "
                "Local-path consistency checks were NOT applied."
            )
            fail_reasons.append("no local paths matched Scale pattern (0 matched)")
            exit_code = 1
        if local_res["errors"]:
            exit_code = 1
            fail_reasons.append(f"{len(local_res['errors'])} local-path errors")

        # S3 SOP checks
        s3_res = validate_s3_scale_raw(mappings)
        meta_count = s3_res.get("metadata_files", 0)
        print(
            f"Scale raw (S3): matched {s3_res['matched']} paths, "
            f"{len(s3_res['errors'])} errors, {len(s3_res['warnings'])} warnings"
            + (f", {meta_count} run-metadata files" if meta_count else "")
        )
        _print_issue_examples("Scale raw (S3)", s3_res["errors"], "errors")
        _print_issue_examples("Scale raw (S3)", s3_res["warnings"], "warnings")
        if s3_res["errors"]:
            exit_code = 1
            fail_reasons.append(f"{len(s3_res['errors'])} S3 SOP errors")

        # SIF completeness if provided
        if args.sif:
            sif_res = validate_sif_completeness_scale(mappings, args.sif)
            missing = sif_res["missing_groupids"]
            extra = sif_res["extra_groupids"]
            missing_assays = sif_res.get("missing_assays", {})
            extra_assays = sif_res.get("extra_assays", {})
            print(
                f"Scale SIF completeness: expected={len(sif_res['expected_groupids'])} GroupIDs, "
                f"found={len(sif_res['actual_groupids'])}, "
                f"missing={len(missing)}, extra={len(extra)}"
            )
            if missing:
                print("  Missing GroupIDs from S3 (present in SIF only): "
                      + ", ".join(sorted(missing)))
                exit_code = 1
                fail_reasons.append(
                    f"{len(missing)} SIF GroupIDs missing from S3: {', '.join(sorted(missing))}"
                )
            if extra:
                print("  Extra GroupIDs in S3 (not present in SIF): "
                      + ", ".join(sorted(extra)))
                exit_code = 1
                fail_reasons.append(
                    f"{len(extra)} extra GroupIDs in S3 not in SIF: {', '.join(sorted(extra))}"
                )

            # Per-GroupID assay completeness
            if missing_assays:
                print("  GroupIDs with missing assay types in S3:")
                for gid in sorted(missing_assays):
                    print(f"    {gid}: missing {sorted(missing_assays[gid])}")
                exit_code = 1
                fail_reasons.append(
                    f"{len(missing_assays)} GroupIDs have missing assay types in S3"
                )
            if extra_assays:
                print("  GroupIDs with unexpected assay types in S3 (not in SIF):")
                for gid in sorted(extra_assays):
                    print(f"    {gid}: unexpected {sorted(extra_assays[gid])}")

            # Summary of expected SIF structure
            expected_ga = sif_res.get("expected_group_assays", {})
            if expected_ga:
                unique_combos: set[frozenset[str]] = set()
                for gid, assays in expected_ga.items():
                    unique_combos.add(frozenset(assays))
                combos_desc = [
                    " + ".join(CANONICAL_ASSAY.get(a, a) for a in sorted(c))
                    for c in sorted(unique_combos, key=sorted)
                ]
                print(f"  SIF assay combinations: {', '.join(combos_desc)}")
        else:
            print("Scale mode: no --sif provided, skipping SIF completeness checks.")

        # Fuzzy S3/local consistency
        cons_res = validate_s3_local_consistency_scale(mappings)
        print(
            f"Scale S3/local consistency: matched {cons_res['matched']} pairs, "
            f"{len(cons_res['errors'])} errors, {len(cons_res['warnings'])} warnings"
        )
        _print_issue_examples("Scale S3/local consistency", cons_res["errors"], "errors")
        _print_issue_examples("Scale S3/local consistency", cons_res["warnings"], "warnings")
        if cons_res["errors"]:
            exit_code = 1
            fail_reasons.append(f"{len(cons_res['errors'])} S3/local consistency errors")

    else:
        print(
            f"Mode (provider={provider}, data={args.data}, assay={args.assay}) "
            "is not implemented yet.",
            file=sys.stderr,
        )
        exit_code = 1
        fail_reasons.append("unsupported validation mode")

    # Final verdict
    print()
    if exit_code == 0:
        print("VERDICT: PASS — mapping validates successfully against SOP rules.")
    else:
        print("VERDICT: FAIL — mapping has issues that need attention:")
        for reason in fail_reasons:
            print(f"  - {reason}")

    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
