"""
Generic DataFrame transforms for DB2 builders.

- Split embedded controlled-term dicts into _term_id / _term_name columns
- Collapse rows by a group key (scalar or list per cell)
"""

import pandas as pd

TERM_ID_SUFFIX = "_term_id"
TERM_NAME_SUFFIX = "_term_name"


def is_empty(val) -> bool:
    """Return True for None, NaN, empty string, or empty list."""
    if val is None:
        return True
    if isinstance(val, float) and pd.isna(val):
        return True
    return val in ("", [])


def to_items(val) -> list:
    """Normalize one cell into zero or more non-empty items for aggregation."""
    if is_empty(val):
        return []
    if isinstance(val, list):
        return [x for x in val if not is_empty(x)]
    return [val]


def single_or_list(series: pd.Series):
    """
    groupby().agg helper: return one scalar, a deduplicated list, or pd.NA.

    Used when multiple rows share the same key and should collapse to one row
    with either a single value or a list of distinct values.
    """
    items = [item for value in series for item in to_items(value)]
    if not items:
        return pd.NA
    unique = list(dict.fromkeys(items))
    return unique[0] if len(unique) == 1 else unique


def collapse_dataframe(
    df: pd.DataFrame,
    group_col: str,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Group by group_col and collapse other columns with single_or_list."""
    other_cols = columns or [c for c in df.columns if c != group_col]
    agg = {col: single_or_list for col in other_cols}
    return df.groupby(group_col, as_index=False).agg(agg)


def extract_term_id_from_ref(ref: str):
    """
    Parse semantic term ID from a controlled-term @id path.

    Example: '/controlled_terms/EFO:0920086/' -> 'EFO:0920086'
    """
    if not ref:
        return pd.NA
    if "/controlled_terms/" in ref:
        term_id = ref.split("/controlled_terms/")[-1]
        return term_id.rstrip("/")
    return ref


def cell_has_term_name_structure(val) -> bool:
    """
    Return True if val is a dict (or list of dicts) containing 'term_name'.

    Columns whose dicts lack term_name (e.g. embedded tissue samples) are skipped.
    """
    if is_empty(val):
        return False
    if isinstance(val, dict):
        return "term_name" in val
    if isinstance(val, list):
        return any(
            isinstance(item, dict) and "term_name" in item for item in val
        )
    return False


def detect_term_name_columns(df: pd.DataFrame) -> list[str]:
    """Return columns where any non-empty cell has term_name dict structure."""
    term_cols = []
    for col in df.columns:
        non_empty = df[col].dropna()
        if non_empty.empty:
            continue
        if non_empty.map(cell_has_term_name_structure).any():
            term_cols.append(col)
    return term_cols


def split_term_cell(val):
    """
    Convert one cell to parallel (term_id, term_name) values.

    Single dict -> scalars; list of dicts -> parallel lists; empty -> pd.NA pair.
    """
    if is_empty(val):
        return pd.NA, pd.NA

    if isinstance(val, dict):
        dicts = [val]
    elif isinstance(val, list):
        dicts = val
    else:
        return pd.NA, pd.NA

    dicts = [d for d in dicts if isinstance(d, dict) and d.get("term_name")]
    if not dicts:
        return pd.NA, pd.NA

    term_ids = [extract_term_id_from_ref(d.get("@id", "")) for d in dicts]
    term_names = [d["term_name"] for d in dicts]

    if len(dicts) == 1:
        return term_ids[0], term_names[0]
    return term_ids, term_names


def split_controlled_term_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace each controlled-term dict column with two columns:
        {col}_term_id   — e.g. EFO:0920086
        {col}_term_name — e.g. 10x gene expression flex v1

    The original column is dropped. Columns without term_name dicts are unchanged.
    """
    term_cols = detect_term_name_columns(df)
    if not term_cols:
        return df

    out = df.copy()
    for col in term_cols:
        id_col = f"{col}{TERM_ID_SUFFIX}"
        name_col = f"{col}{TERM_NAME_SUFFIX}"

        pairs = out[col].map(split_term_cell)
        out[id_col] = pairs.map(lambda p: p[0])
        out[name_col] = pairs.map(lambda p: p[1])
        out = out.drop(columns=[col])

    print(f"Split controlled-term columns: {term_cols}")
    return out
