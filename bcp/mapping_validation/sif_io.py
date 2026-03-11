from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Set

import warnings

import pandas as pd

from .constants import CANONICAL_ASSAY


warnings.filterwarnings(
    "ignore",
    message="Data Validation extension is not supported",
    category=UserWarning,
    module="openpyxl",
)


def _normalize_sif_groupid(gid: str) -> str:
    """Normalise a SIF Group Identifier to match the S3 directory convention.

    SIF files for 10x use ``A + AF`` style (space-plus-space) while S3
    paths join the same parts with underscores: ``A_AF``.  Scale SIFs
    already use plain identifiers, so this is a no-op for them.
    """
    return gid.replace(" + ", "_")


def _find_col(cols_lower: Dict[str, str], keyword: str) -> str | None:
    """Find a column whose lower-cased name contains ``keyword``."""
    for key, name in cols_lower.items():
        if keyword in key:
            return name
    return None


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
    import re

    p = Path(sif_path)
    result: dict[str, set[str]] = defaultdict(set)

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
    import csv

    p = Path(sif_path)
    result: dict[str, str] = {}

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


def load_sif_library_names(sif_path: str | Path) -> set[str]:
    """Load the set of unique Library Names from a SIF file.

    This is a fallback for SIFs that lack a "Group Identifier" column
    (e.g. Ultima intake forms where each library is a standalone sample).
    Searches for a column whose name contains "library name" and returns
    the set of non-empty values.
    """
    import csv

    p = Path(sif_path)
    result: set[str] = set()

    if p.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        try:
            df = pd.read_excel(p)
        except Exception:
            df = None
        if df is not None and not df.empty:
            cols = {str(c).strip().lower(): c for c in df.columns}
            lib_col = _find_col(cols, "library name")
            if lib_col is not None:
                for _, row_data in df.iterrows():
                    lib = str(row_data.get(lib_col, "")).strip()
                    if lib and lib != "nan":
                        result.add(lib)
        if result:
            return result

    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames:
            field_map = {name.lower(): name for name in reader.fieldnames}
            lib_col_name = _find_col(field_map, "library name")
            if lib_col_name is not None:
                for row in reader:
                    lib = (row.get(lib_col_name) or "").strip()
                    if lib:
                        result.add(lib)

    return result


__all__ = [
    "_normalize_sif_groupid",
    "load_sif_group_assays",
    "load_sif_scale_group_assays",
    "load_sif_scale_groupids",
    "load_sif_library_assays",
    "load_sif_library_names",
]

