"""
Generic functions for gathering DB2 object information, and DataFrame transforms for DB2 builder

- Split embedded controlled-term dicts into _term_id / _term_name columns
- Collapse rows by a group key (scalar or list per cell)
- Get url prefix from an object id
- Get API type from an object id
- Get object type from object id
- Get uuid from object id
- Get reference object ids from a field
"""

import re

import pandas as pd
import numpy as np
from constants import Configs

TERM_ID_SUFFIX = "_term_id"
TERM_NAME_SUFFIX = "_term_name"


def is_empty(val) -> bool:
    """Return True for None, NaN, empty string, or empty list."""
    if val is None:
        return True
    if isinstance(val, float) and pd.isna(val):
        return True
    return val in ("", [])


def strip_author_metadata_column_prefix(df: pd.DataFrame) -> pd.DataFrame:
    """Rename cols like tissues_author_metadata_foo -> foo."""
    rename = {
        c: re.sub(r'^.*?_author_metadata_', '', c)
        for c in df.columns
        if '_author_metadata_' in c
    }
    return df.rename(columns=rename) if rename else df


def collapse_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse columns that share the same name into one column.
    For each duplicated label, per row:
    - ignore empty values (is_empty)
    - if one distinct value remains, use it
    - if several distinct values remain, keep them as a list
      (same idea as single_or_list)
    - if none, use pd.NA
    Columns with unique names are left unchanged.
    """
    out = pd.DataFrame(index=df.index)
    seen: set[str] = set()
    for col in df.columns:
        if col in seen:
            continue
        seen.add(col)
        # df[col] is a Series if unique, DataFrame if duplicated labels
        block = df.loc[:, df.columns == col]
        if isinstance(block, pd.Series) or block.shape[1] == 1:
            out[col] = block.iloc[:, 0] if isinstance(block, pd.DataFrame) else block
            continue
        def _collapse_row(row: pd.Series):
            items = []
            for val in row.tolist():
                items.extend(to_items(val))
            if not items:
                return pd.NA
            unique = list(dict.fromkeys(items))
            return unique[0] if len(unique) == 1 else unique
        out[col] = block.apply(_collapse_row, axis=1)
    return out


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


def combine_bound_columns(
    df: pd.DataFrame,
    *,
    lower_col: str,
    upper_col: str,
    units_col: str | None = None,
    out_col: str,
    drop_source: bool = True,
) -> pd.DataFrame:
    """
    Combine lower/upper bound columns into one value, optionally with units.

    Equal bounds -> "{value}"
    Unequal     -> "{lower}-{upper}"
    With units  -> append " {units}"
    """
    missing = [c for c in (lower_col, upper_col) if c not in df.columns]
    if missing:
        return df

    lower = df[lower_col]
    upper = df[upper_col]
    equal = lower.eq(upper) | (lower.isna() & upper.isna())

    lower_s = lower.map(lambda v: "" if is_empty(v) else str(v))
    upper_s = upper.map(lambda v: "" if is_empty(v) else str(v))

    combined = np.where(equal, lower_s, lower_s + "-" + upper_s)

    if units_col and units_col in df.columns:
        units_s = df[units_col].map(lambda v: "" if is_empty(v) else str(v))
        combined = [
            f"{val} {unit}".strip() if unit else val
            for val, unit in zip(combined, units_s, strict=True)
        ]

    df[out_col] = combined

    if drop_source:
        cols = [lower_col, upper_col] + ([units_col] if units_col else [])
        df.drop(columns=[c for c in cols if c in df.columns], inplace=True)
    
    return df


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

    term_ids = [extract_controlled_term_id(d.get("@id", "")) for d in dicts]
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


# Gathering DB2 Objects generic functions

def extract_controlled_term_id(ref: str):
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


def match_object_config_from_id(object_id: str, configs: Configs):
    """Return (url_prefix, config) for an @id path, or (None, None)."""
    if not object_id:
        return None, None
    for path_key, config in configs.OBJECT_CONFIG.items():
        if f'/{path_key}/' in object_id:
            return path_key, config
    return None, None


def get_url_prefix_from_id(object_id: str, configs: Configs):
    """'/human_donors/uuid/' -> 'human_donors'"""
    prefix, _ = match_object_config_from_id(object_id, configs)
    return prefix


def get_api_type_from_id(object_id: str, configs: Configs):
    """'/human_donors/uuid/' -> 'HumanDonor'"""
    _, config = match_object_config_from_id(object_id, configs)
    return config['api_type'] if config else None


def get_config_obj_type(obj: dict, configs: Configs) -> str:
    """Map a resolved object to OBJECT_CONFIG url_prefix via @id, else @type."""
    obj_id = obj.get('@id', '')
    if obj_id:
        prefix = get_url_prefix_from_id(obj_id, configs)
        if prefix:
            return prefix
    for config_obj_type, config_obj_info in configs.OBJECT_CONFIG.items():
        if config_obj_info.get('api_type', '') == obj.get('@type', [''])[0]:
            return config_obj_type
    return ''


def extract_uuid_from_id(object_id: str) -> str:
    """'/tissues/uuid/' -> 'uuid'"""
    if '/' in object_id:
        return object_id.split('/')[-2] if object_id.endswith('/') else object_id.split('/')[-1]
    return object_id


def extract_references_from_field(field_value, field_name, configs: Configs) -> list:
    """Extract reference @id paths from a field using FIELD_TYPES."""
    if not field_value:
        return []
    field_spec = configs.FIELD_TYPES.get(field_name, {'type': 'string'})
    refs = []
    if field_spec['type'] == 'array':
        if isinstance(field_value, list):
            for item in field_value:
                if isinstance(item, dict):
                    ref_id = item.get('@id')
                    if ref_id:
                        refs.append(ref_id)
                elif isinstance(item, str):
                    refs.append(item)
    else:
        if isinstance(field_value, dict):
            ref_id = field_value.get('@id')
            if ref_id:
                refs.append(ref_id)
        elif isinstance(field_value, str):
            refs.append(field_value)
    return refs
