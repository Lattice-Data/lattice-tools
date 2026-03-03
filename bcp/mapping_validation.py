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

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple
import argparse
import re
import sys


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
        # NVUS order numbers are of the form NVUS{digits}-{batch},
        # but the exact digit count can vary, so keep this permissive.
        order_pattern = r"NVUS\d+-\d+"
        valid_assays = {"GEX", "CRI", "ATAC"}
    else:
        order_pattern = r"AN\d+"
        valid_assays = {"GEX", "CRI", "ATAC", "viral_ORF"}

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

    for row in mappings:
        total += 1
        m = s3_regex.match(row.s3_path)
        if not m:
            # Not every mapping in a file must be 10x raw; record as a parse warning.
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
        "total": total,
        "matched": matched,
    }


def validate_local_paths_scale_raw(mappings: Iterable[MappingRow]) -> dict:
    """
    Validate local paths for Novogene Scale/Quantum raw data.

    This focuses purely on local filesystem paths and checks for
    self-consistency of the QSR number, SCALEPLEX flag, and index
    sequence within each path.

    Expected (index-direct) local layout:

        {base}/{wafer}-{date}/
            {wafer}-QSR-{N}[-SCALEPLEX]-{INDEX}/
                {wafer}-QSR-{N}[-SCALEPLEX]-{INDEX}{suffix}

    Where:
    - wafer: numeric run ID (e.g. 441969)
    - N: QSR number
    - INDEX: A/C/G/T sequence

    Returns a dict with:
    - errors: list of issue dicts
    - warnings: list of issue dicts (e.g. paths that do not match the pattern)
    - total: number of mappings inspected
    - matched: number of local paths that matched the expected pattern
    """
    pattern = re.compile(
        r"(?P<base_path>.+)/"
        r"(?P<wafer>\d+)-(?P<date>\d+_\d+)/"
        r"(?P<wafer2>\d+)-QSR-(?P<qsr_n>\d+)(?P<scaleplex>-SCALEPLEX)?-(?P<index>[ACGT]+)/"
        r"(?P<wafer3>\d+)-QSR-(?P<qsr_n2>\d+)(?P<scaleplex2>-SCALEPLEX)?-(?P<index2>[ACGT]+)"
        r"(?P<suffix>.*)"
    )

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1
        local = row.local_path
        m = pattern.match(local)
        if not m:
            warnings.append(
                {
                    "type": "parse_miss",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": "local path does not match expected Scale index-direct pattern",
                }
            )
            continue

        matched += 1
        gd = m.groupdict()

        # Wafer (run ID) should be consistent throughout the path
        if not (gd["wafer"] == gd["wafer2"] == gd["wafer3"]):
            errors.append(
                {
                    "type": "wafer_mismatch",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": f"wafer IDs differ within path: {gd['wafer']}, {gd['wafer2']}, {gd['wafer3']}",
                }
            )

        # QSR numbers must match
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

        # SCALEPLEX flags must be consistent
        if (gd["scaleplex"] is None) != (gd["scaleplex2"] is None):
            errors.append(
                {
                    "type": "scaleplex_mismatch",
                    "line": row.line_num,
                    "local_path": local,
                    "detail": "SCALEPLEX present in one component of the path but not the other",
                }
            )

        # Index sequences must match
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
    - Assay is one of: GEX, Hash_oligo, GEX_hash_oligo (others are treated as invalid)
    - Hash_oliga (typo) is explicitly flagged as an error.
    - For UG_RT form:
        * Hash_oligo and GEX_hash_oligo must include SCALEPLEX
        * GEX must NOT include SCALEPLEX

    Returns a dict with:
    - errors: list of issue dicts
    - warnings: list of issue dicts
    - total: number of mappings inspected
    - matched: number of S3 paths that matched either pattern
    """
    # Project and order components: trapnell-seahub-bcp / NVUS...
    order_pattern = r"NVUS\d+-\d+"

    sop_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        r"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>[A-Za-z_]+)_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)"
        r"(?P<suffix>.*)$"
    )

    index_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        r"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>[A-Za-z_]+)_(?P<index_seq>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )

    errors: List[dict] = []
    warnings: List[dict] = []
    total = 0
    matched = 0

    for row in mappings:
        total += 1
        s3 = row.s3_path

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

        # Explicitly flag common typo anywhere in the path
        if "Hash_oliga" in s3:
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "assay 'Hash_oliga' is not valid (did you mean 'Hash_oligo'?)",
                }
            )
        # Determine assay type using substrings to be more robust to parsing quirks
        if "GEX_hash_oligo" in s3:
            assay = "GEX_hash_oligo"
        elif "Hash_oligo" in s3:
            assay = "Hash_oligo"
        elif "GEX" in s3:
            assay = "GEX"
        else:
            errors.append(
                {
                    "type": "invalid_assay",
                    "line": row.line_num,
                    "s3_path": s3,
                    "detail": "unable to determine assay type from S3 path "
                    "(expected GEX, Hash_oligo, or GEX_hash_oligo)",
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

            if assay in {"Hash_oligo", "GEX_hash_oligo"} and not has_scaleplex:
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
        "total": total,
        "matched": matched,
    }


def load_sif_scale_groupids(sif_path: str | Path) -> set[str]:
    """
    Load expected GroupIDs for Scale libraries from a SIF CSV.

    The SIF used for Scale/Quantum includes a \"Group Identifier\" column
    that contains library-level identifiers such as R096A, R096G, etc.
    This helper extracts the non-empty values from that column.
    """
    import csv

    p = Path(sif_path)
    expected: set[str] = set()

    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return expected

        # Find the group identifier column in a case-insensitive way
        field_map = {name.lower(): name for name in reader.fieldnames}
        group_col_name = None
        for key, name in field_map.items():
            if "group identifier" in key:
                group_col_name = name
                break

        if group_col_name is None:
            return expected

        for row in reader:
            val = (row.get(group_col_name) or "").strip()
            if val:
                expected.add(val)

    return expected


def validate_sif_completeness_scale(
    mappings: Iterable[MappingRow], sif_path: str | Path
) -> dict:
    """
    Compare Scale SIF GroupIDs against GroupIDs observed in S3 Scale raw paths.

    - expected_groupids: from SIF \"Group Identifier\" column
    - actual_groupids: from S3 paths that match the Scale raw patterns
    """
    expected_groupids = load_sif_scale_groupids(sif_path)

    # Reuse the same patterns as validate_s3_scale_raw
    order_pattern = r"NVUS\d+-\d+"
    sop_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        r"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>[A-Za-z_]+)_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)"
        r"(?P<suffix>.*)$"
    )
    index_pattern = re.compile(
        r"^s3://(?P<bucket>czi-novogene)/(?P<project>[a-z0-9-]+)/"
        rf"(?P<order>{order_pattern})/(?P<experiment_id>[^/]+)/raw/"
        r"(?P<runid>\d+)/"
        r"(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>[A-Za-z_]+)_(?P<index_seq>[ACGT]+)"
        r"(?P<suffix>.*)$"
    )

    actual_groupids: set[str] = set()
    for row in mappings:
        s3 = row.s3_path
        m = sop_pattern.match(s3) or index_pattern.match(s3)
        if not m:
            continue
        gd = m.groupdict()
        actual_groupids.add(gd["group_id"])

    missing = expected_groupids - actual_groupids
    extra = actual_groupids - expected_groupids

    return {
        "expected_groupids": expected_groupids,
        "actual_groupids": actual_groupids,
        "missing_groupids": missing,
        "extra_groupids": extra,
    }


def main() -> None:
    """
    Simple CLI for validating mapping CSV/TSV files.

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
    parser.add_argument("--sif", help="Path to SIF CSV for completeness checks (Scale)")
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

    # 1. Uniqueness
    uniq = validate_uniqueness(mappings)
    dup_local_count = len(uniq["duplicate_locals"])
    dup_s3_count = len(uniq["duplicate_s3"])
    print(f"Uniqueness: {uniq['total']} mappings, "
          f"{dup_local_count} duplicate locals, {dup_s3_count} duplicate S3 paths")
    if dup_local_count or dup_s3_count:
        exit_code = 1

    # 2. Mode-specific validation
    provider = args.provider
    if args.data == "raw" and args.assay == "10x":
        res = validate_s3_10x_raw(provider, mappings)
        print(f"10x raw SOP: matched {res['matched']} S3 paths, "
              f"{len(res['errors'])} errors, {len(res['warnings'])} warnings")
        if res["errors"]:
            exit_code = 1

    elif args.data == "raw" and args.assay == "scale" and provider == "novogene":
        # Local path sanity
        local_res = validate_local_paths_scale_raw(mappings)
        print(f"Scale raw (local): matched {local_res['matched']} paths, "
              f"{len(local_res['errors'])} errors, {len(local_res['warnings'])} warnings")
        if local_res["errors"]:
            exit_code = 1

        # S3 SOP checks
        s3_res = validate_s3_scale_raw(mappings)
        print(f"Scale raw (S3): matched {s3_res['matched']} paths, "
              f"{len(s3_res['errors'])} errors, {len(s3_res['warnings'])} warnings")
        if s3_res["errors"]:
            exit_code = 1

        # SIF completeness if provided
        if args.sif:
            sif_res = validate_sif_completeness_scale(mappings, args.sif)
            missing = sif_res["missing_groupids"]
            extra = sif_res["extra_groupids"]
            print(
                f"Scale SIF completeness: expected={len(sif_res['expected_groupids'])}, "
                f"found={len(sif_res['actual_groupids'])}, "
                f"missing={len(missing)}, extra={len(extra)}"
            )
            if missing:
                print("  Missing GroupIDs from S3 (present in SIF only): "
                      + ", ".join(sorted(missing)))
                exit_code = 1
        else:
            print("Scale mode: no --sif provided, skipping SIF completeness checks.")

    else:
        print(
            f"Mode (provider={provider}, data={args.data}, assay={args.assay}) "
            "is not implemented yet.",
            file=sys.stderr,
        )
        exit_code = 1

    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
