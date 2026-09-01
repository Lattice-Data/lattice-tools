from __future__ import annotations

import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable

import re

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


def _group_id_match_keywords(provider: str | None) -> tuple[str, ...]:
    """Return ordered substrings used to locate the Group ID column for a vendor."""
    p = (provider or "").strip().lower()
    if p == "psomagen":
        return ("group id", "group identifier", "groupid")
    return ("group identifier", "group id", "groupid")


def _header_match_keys(orig: str) -> set[str]:
    """Lower-cased header variants for substring matching (footnotes, Unicode marks)."""
    raw = str(orig).strip().lower()
    keys = {raw, re.sub(r"[\*\^]+$", "", raw).strip()}
    plain = "".join(
        ch for ch in unicodedata.normalize("NFKD", raw) if not unicodedata.combining(ch)
    )
    keys.add(plain)
    keys.add(re.sub(r"[\*\^]+$", "", plain).strip())
    return {k for k in keys if k}


def _resolve_group_column(columns: Iterable[str], provider: str | None) -> str | None:
    """Pick the spreadsheet column that holds Group IDs (Psomagen vs Novogene headers)."""
    keywords = _group_id_match_keywords(provider)
    cols_list = list(columns)
    for kw in keywords:
        for orig in cols_list:
            for key in _header_match_keys(orig):
                if kw in key:
                    return str(orig)
    return None


def _assay_header_keywords(provider: str | None) -> tuple[str, ...]:
    """Ordered substrings to match spreadsheet *headers* for the assay column."""
    p = (provider or "").strip().lower()
    if p == "psomagen":
        return (
            # Long composite header on Psomagen order forms, e.g.
            # "Analysis Target / Feature Barcode / Supplemental libraries*"
            "analysis target",
            "supplemental libraries",
            "feature barcode",
            "assay type",
            "library type",
            "application type",
            "application",
            "experiment type",
            "sequencing type",
            "10x product",
            "product type",
            "data type",
            "service type",
            "assay",
        )
    return ("assay type", "assay")


def _resolve_assay_column(
    columns: Iterable[str], provider: str | None = None
) -> str | None:
    """Pick the assay / application column (Psomagen forms often omit 'Assay' in the name)."""
    cols = {str(c).strip().lower(): c for c in columns}
    for kw in _assay_header_keywords(provider):
        hit = _find_col(cols, kw)
        if hit:
            return hit
    for orig in columns:
        if str(orig).strip().lower() == "assay":
            return str(orig)
    for orig in columns:
        low = str(orig).strip().lower()
        if re.search(r"\bassay\b", low) and "assembly" not in low:
            return str(orig)
    return None


def _normalize_sif_assay_token(raw: str, provider: str | None = None) -> str:
    """Map SIF assay / application cell text to lower-case SOP tokens (gex, cri, …)."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return ""
    s = str(raw).strip().lower()
    if not s or s == "nan":
        return ""
    p = (provider or "").strip().lower()
    if s in CANONICAL_ASSAY:
        return CANONICAL_ASSAY[s].lower()
    if "hash" in s and "oligo" in s:
        return "hash_oligo"
    if "gex" in s and "hash" in s:
        return "gex_hash_oligo"
    if re.search(r"\batac\b", s):
        return "atac"
    if "viral" in s and "orf" in s:
        return "viral_orf"
    if re.search(r"\b(gex|gene\s*expression|scrna|sc\s*rna)\b", s) or (
        "gene expression" in s
    ):
        return "gex"
    if p in {"scale", "sci"} and "feature" in s and "barcode" in s:
        return "hash_oligo"
    if re.search(r"\b(cri|crispr)\b", s) or (
        "feature" in s and "barcode" in s and p not in {"scale", "sci"}
    ):
        return "cri"
    if "cell surface" in s and "protein" in s:
        return "cri"
    if len(s) <= 24:
        return s
    return ""


# Psomagen order forms ship with template example rows (e.g. AC_Sample1) in the
# Group ID column; these are not real samples and must not affect SIF completeness.
_PSOMAGEN_SIF_EXAMPLE_GROUPID_RE: re.Pattern[str] = re.compile(
    r"^AC_Sample\d+$", re.IGNORECASE
)


def _is_psomagen_sif_example_groupid(gid: str, provider: str | None) -> bool:
    """Return True for Psomagen order-form template Group IDs (not real samples)."""
    if (provider or "").strip().lower() != "psomagen":
        return False
    return bool(_PSOMAGEN_SIF_EXAMPLE_GROUPID_RE.match(gid.strip()))


# Psomagen / CZI order forms: long cover pages; header is often row ~20–40.
SIF_EXCEL_SCAN_NROWS: int = 600
SIF_EXCEL_MAX_HEADER_TRIES: int = 24


def _nonblank_row_cells(row: pd.Series) -> list[str]:
    cells: list[str] = []
    for v in row:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            continue
        s = str(v).strip()
        if s and s.lower() != "nan":
            cells.append(s)
    return cells


def _row_looks_like_sif_header(
    row: pd.Series, provider: str | None, *, require_group: bool
) -> bool:
    """Heuristic: row looks like a SIF column-header row (not a data or title line)."""
    cells = _nonblank_row_cells(row)
    if len(cells) < 2:
        return False
    lowered = [c.lower() for c in cells]
    has_assay = any(
        re.search(r"\b(assay|application)\b", low)
        or "analysis target" in low
        or ("supplemental" in low and "librar" in low)
        for low in lowered
    )
    has_lib = any("library" in low for low in lowered)
    kws_g = _group_id_match_keywords(provider)
    has_group = False
    for s in cells:
        for key in _header_match_keys(s):
            if any(kw in key for kw in kws_g):
                has_group = True
                break
        if has_group:
            break
    if require_group:
        return has_group and (has_assay or has_lib)
    return has_lib and has_assay


def _excel_sheet_names(path: Path) -> list[str] | None:
    for eng in ("openpyxl", None):
        try:
            kwargs = {"engine": eng} if eng else {}
            return list(pd.ExcelFile(path, **kwargs).sheet_names)
        except Exception:
            continue
    return None


def _read_excel_raw_preview(path: Path, sheet: str, nrows: int) -> pd.DataFrame | None:
    for eng in ("openpyxl", None):
        try:
            kwargs: dict = {"sheet_name": sheet, "header": None, "nrows": nrows}
            if eng:
                kwargs["engine"] = eng
            return pd.read_excel(path, **kwargs)
        except Exception:
            continue
    return None


def _read_excel_sif_sheet(
    path: Path, sheet: str, header_row: int
) -> pd.DataFrame | None:
    for eng in ("openpyxl", None):
        try:
            kwargs: dict = {"sheet_name": sheet, "header": header_row}
            if eng:
                kwargs["engine"] = eng
            df = pd.read_excel(path, **kwargs)
        except Exception:
            continue
        if df is not None and not df.empty:
            return df
    return None


def _header_row_candidates_for_sheet(
    path: Path,
    sheet: str,
    provider: str | None,
    *,
    require_group: bool,
) -> list[int]:
    """Return a small ordered list of 0-based ``header=`` indices to try for this sheet."""
    raw = _read_excel_raw_preview(path, sheet, SIF_EXCEL_SCAN_NROWS)
    if raw is None:
        return [0]

    n_preview = len(raw)
    scored: list[tuple[int, int]] = []
    for i in range(n_preview):
        if _row_looks_like_sif_header(
            raw.iloc[i], provider, require_group=require_group
        ):
            scored.append((int(raw.iloc[i].notna().sum()), i))

    if not scored and require_group:
        for i in range(n_preview):
            if _row_looks_like_sif_header(raw.iloc[i], provider, require_group=False):
                scored.append((int(raw.iloc[i].notna().sum()), i))

    if scored:
        scored.sort(key=lambda x: (-x[0], -x[1]))
        ordered: list[int] = []
        seen: set[int] = set()
        for _, idx in scored:
            if idx not in seen:
                ordered.append(idx)
                seen.add(idx)
        merged = [0] + [i for i in ordered if i != 0]
        seen.clear()
        out: list[int] = []
        for idx in merged:
            if idx not in seen:
                seen.add(idx)
                out.append(idx)
        return out[:SIF_EXCEL_MAX_HEADER_TRIES]

    return list(range(0, min(80, n_preview)))[:SIF_EXCEL_MAX_HEADER_TRIES]


def _parse_group_assays_from_dataframe(
    df: pd.DataFrame, provider: str | None
) -> dict[str, set[str]]:
    """Extract GroupID → assay types from one dataframe, or {} if no group column."""
    result: dict[str, set[str]] = defaultdict(set)
    group_col = _resolve_group_column(df.columns, provider)
    if group_col is None:
        return {}
    assay_col = _resolve_assay_column(df.columns, provider)
    assay_series = (
        df[assay_col].ffill()
        if assay_col is not None and assay_col in df.columns
        else None
    )
    for pos in range(len(df)):
        row_data = df.iloc[pos]
        gid = str(row_data.get(group_col, "")).strip()
        if not gid or gid == "nan":
            continue
        if _is_psomagen_sif_example_groupid(gid, provider):
            continue
        raw_assay = assay_series.iloc[pos] if assay_series is not None else ""
        token = _normalize_sif_assay_token(str(raw_assay), provider)
        if token:
            result[gid].add(token)
        else:
            result.setdefault(gid, set())
    return dict(result) if result else {}


def _load_group_assays_from_excel(
    path: Path, provider: str | None
) -> dict[str, set[str]]:
    """Scan sheets with preview-based header detection until Group IDs parse."""
    names = _excel_sheet_names(path)
    if names is None:
        for eng in ("openpyxl", None):
            try:
                kwargs: dict = {}
                if eng:
                    kwargs["engine"] = eng
                df = pd.read_excel(path, **kwargs)
            except Exception:
                continue
            if df is None or df.empty:
                continue
            return _parse_group_assays_from_dataframe(df, provider)
        return {}
    for sheet in names:
        for hr in _header_row_candidates_for_sheet(
            path, sheet, provider, require_group=True
        ):
            df = _read_excel_sif_sheet(path, sheet, hr)
            if df is None or df.empty:
                continue
            parsed = _parse_group_assays_from_dataframe(df, provider)
            if parsed:
                return parsed
    return {}


def load_sif_group_assays(
    sif_path: str | Path, provider: str | None = None
) -> dict[str, set[str]]:
    """
    Load expected GroupID → {assay types} mapping from a SIF file.

    Each row in the SIF represents a sublibrary for a group identifier.
    The same group identifier appears once per assay type (e.g. GEX and
    CRI for 10x, or GEX and hash_oligo for Scale).  Returns a dict
    mapping each GroupID (as it appears in the SIF) to its set of assay
    types (normalised to lower-case).

    ``provider`` affects which column header is preferred for locating Group IDs
    (e.g. Psomagen ``Group ID^`` vs Novogene ``Group Identifier``).  When omitted,
    ``Group Identifier`` is tried before ``Group ID``.

    Falls back through three parsing strategies:
    1. Excel (.xlsx/.xlsm/.xls) with pandas — preview the first ~600 rows of each
       sheet to find a SIF-like header row (Group + Library/Assay labels), then
       read with that ``header=`` index (forms with banners through row ~26).
    2. Well-formed CSV with DictReader
    3. Heuristic Ultima intake-style CSV (Scale-specific)
    """
    import csv

    p = Path(sif_path)
    result: dict[str, set[str]] = defaultdict(set)

    # Branch 1: Excel SIF
    if p.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        from_excel = _load_group_assays_from_excel(p, provider)
        if from_excel:
            return from_excel

    # Branch 2: well-formed CSV
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames:
            group_col_name = _resolve_group_column(reader.fieldnames, provider)
            assay_col_name = _resolve_assay_column(reader.fieldnames, provider)

            if group_col_name is not None:
                for row in reader:
                    gid = (row.get(group_col_name) or "").strip()
                    if not gid:
                        continue
                    if _is_psomagen_sif_example_groupid(gid, provider):
                        continue
                    assay_val = (
                        (row.get(assay_col_name) or "").strip()
                        if assay_col_name
                        else ""
                    )
                    token = _normalize_sif_assay_token(assay_val, provider)
                    if token:
                        result[gid].add(token)
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


def load_sif_scale_groupids(
    sif_path: str | Path, provider: str | None = None
) -> set[str]:
    """
    Load expected GroupIDs for Scale libraries from a SIF file.

    Convenience wrapper around :func:`load_sif_group_assays` that
    returns just the set of GroupIDs.
    """
    return set(load_sif_group_assays(sif_path, provider=provider).keys())


def _coerce_cell_to_str(value: object) -> str:
    """Stringify a spreadsheet cell without pandas' float artefacts.

    A column holding whole numbers alongside any blank cell is read as
    ``float64``, so a Group Identifier of ``1234`` would otherwise come back as
    ``"1234.0"`` and never match the ``1234`` S3 directory.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    s = str(value).strip()
    return "" if s.lower() == "nan" else s


def _library_records_from_dataframe(
    df: pd.DataFrame, provider: str | None
) -> tuple[dict[str, str], dict[str, set[str]]] | None:
    """Extract library → assay and library → GroupIDs from one SIF sheet.

    Returns ``None`` when the sheet does not carry the Library name and assay
    columns, so callers can move on to the next header row / sheet.  The group
    mapping is empty when the sheet has no Group Identifier column.
    """
    cols = {str(c).strip().lower(): c for c in df.columns}
    lib_col = _find_col(cols, "library name")
    assay_col = _resolve_assay_column(df.columns, provider)
    if lib_col is None or assay_col is None:
        return None
    group_col = _resolve_group_column(df.columns, provider)

    assay_ff = df[assay_col].ffill()
    assays: dict[str, str] = {}
    groups: dict[str, set[str]] = {}
    for pos in range(len(df)):
        row_data = df.iloc[pos]
        lib = _coerce_cell_to_str(row_data.get(lib_col))
        if not lib:
            continue
        assay = _normalize_sif_assay_token(
            _coerce_cell_to_str(assay_ff.iloc[pos]), provider
        )
        if assay:
            assays[lib] = assay
        if group_col is not None:
            # Deliberately not forward-filled: a blank Group Identifier means
            # "unknown", and inventing one would let an ungrouped library be
            # reported as belonging to the group above it.
            gid = _coerce_cell_to_str(row_data.get(group_col))
            if gid and not _is_psomagen_sif_example_groupid(gid, provider):
                groups.setdefault(lib, set()).add(_normalize_sif_groupid(gid))

    if not assays:
        return None
    return assays, groups


def load_sif_libraries(
    sif_path: str | Path, provider: str | None = None
) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Load ``(library → assay, library → GroupIDs)`` from a SIF file.

    Both mappings are read from the *same* sheet and header row, so they can
    never end up describing different parts of a workbook (a stale tab, or an
    example block above the real header).

    Assay types are lower-cased SOP tokens.  Group Identifiers are normalised
    with :func:`_normalize_sif_groupid` so they can be compared directly
    against S3 GroupID directory names.  A library maps to a *set* of groups
    because nothing guarantees a library name appears only once in a SIF.

    ``provider`` selects Psomagen-style assay/group column headers when needed.
    """
    import csv

    p = Path(sif_path)

    if p.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        names = _excel_sheet_names(p)
        if names:
            for sheet in names:
                for hr in _header_row_candidates_for_sheet(
                    p, sheet, None, require_group=False
                ):
                    df = _read_excel_sif_sheet(p, sheet, hr)
                    if df is None or df.empty:
                        continue
                    records = _library_records_from_dataframe(df, provider)
                    if records is not None:
                        return records

    assays: dict[str, str] = {}
    groups: dict[str, set[str]] = {}
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames:
            field_map = {name.lower(): name for name in reader.fieldnames}
            lib_col_name = _find_col(field_map, "library name")
            assay_col_name = _resolve_assay_column(reader.fieldnames, provider)
            group_col_name = _resolve_group_column(reader.fieldnames, provider)
            if lib_col_name is not None and assay_col_name is not None:
                for row in reader:
                    lib = (row.get(lib_col_name) or "").strip()
                    if not lib:
                        continue
                    assay = _normalize_sif_assay_token(
                        (row.get(assay_col_name) or "").strip(), provider
                    )
                    if assay:
                        assays[lib] = assay
                    if group_col_name is not None:
                        gid = (row.get(group_col_name) or "").strip()
                        if gid and not _is_psomagen_sif_example_groupid(gid, provider):
                            groups.setdefault(lib, set()).add(
                                _normalize_sif_groupid(gid)
                            )

    return assays, groups


def load_sif_library_assays(
    sif_path: str | Path, provider: str | None = None
) -> dict[str, str]:
    """Load a Library-Name → assay-type mapping from a SIF file.

    Returns a dict mapping each individual library name to its assay type
    (lower-cased).  For 10x data the SIF typically has one row per
    library, so ``FTF1732A`` → ``gex`` and ``FTF1732AF`` → ``cri``.

    ``provider`` selects Psomagen-style assay column headers when needed.
    """
    return load_sif_libraries(sif_path, provider)[0]


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
        sheet_names = _excel_sheet_names(p)
        if sheet_names:
            for sheet in sheet_names:
                for hr in _header_row_candidates_for_sheet(
                    p, sheet, None, require_group=False
                ):
                    df = _read_excel_sif_sheet(p, sheet, hr)
                    if df is None or df.empty:
                        continue
                    cols = {str(c).strip().lower(): c for c in df.columns}
                    lib_col = _find_col(cols, "library name")
                    if lib_col is None:
                        continue
                    names: set[str] = set()
                    for _, row_data in df.iterrows():
                        lib = str(row_data.get(lib_col, "")).strip()
                        if lib and lib != "nan":
                            names.add(lib)
                    if names:
                        return names

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


def load_sif_library_groups(
    sif_path: str | Path, provider: str | None = None
) -> dict[str, set[str]]:
    """Load a Library-Name → Group-Identifier(s) mapping from a SIF file.

    The SIF states this relationship explicitly: each library row carries the
    Group Identifier it belongs to.  Reading it avoids having to *infer* the
    relationship from the text of the names, which differs per convention:

    * paired 10x -- Group Identifier is a concatenation of its member library
      names (``LIB1`` and ``LIB1F`` both live under ``LIB1_LIB1F``);
    * multiome 10x -- library names extend the Group Identifier with an assay
      suffix (``CH01GEX`` and ``CH01ATAC`` both live under ``CH01``);
    * Psomagen 10x -- library name and Group Identifier are identical.

    Returns an empty dict when the SIF has no Group Identifier column.  See
    :func:`load_sif_libraries` for the details of how the sheet is chosen.
    """
    return load_sif_libraries(sif_path, provider)[1]


__all__ = [
    "_normalize_sif_groupid",
    "load_sif_group_assays",
    "load_sif_scale_group_assays",
    "load_sif_scale_groupids",
    "load_sif_libraries",
    "load_sif_library_assays",
    "load_sif_library_groups",
    "load_sif_library_names",
]
