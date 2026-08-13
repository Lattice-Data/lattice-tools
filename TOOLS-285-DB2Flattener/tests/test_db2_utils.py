import numpy as np
import pandas as pd
import pytest

from DB2_utils import (
    collapse_dataframe,
    collapse_duplicate_columns,
    combine_bound_columns,
    extract_controlled_term_id,
    sort_ontology_term_id_column,
    split_term_cell,
    strip_author_metadata_column_prefix,
)


def test_combine_bound_columns_equal_and_unequal():
    df = pd.DataFrame(
        {
            "lower": [5, 1],
            "upper": [5, 5],
        }
    )
    result = combine_bound_columns(
        df, lower_col="lower", upper_col="upper", out_col="combined"
    )
    assert list(result["combined"]) == ["5", "1-5"]
    assert "lower" not in result.columns
    assert "upper" not in result.columns


def test_combine_bound_columns_with_units():
    df = pd.DataFrame(
        {
            "lower": [5, 1],
            "upper": [5, 5],
            "units": ["years", "years"],
        }
    )
    result = combine_bound_columns(
        df,
        lower_col="lower",
        upper_col="upper",
        units_col="units",
        out_col="combined",
    )
    assert list(result["combined"]) == ["5 years", "1-5 years"]
    assert "units" not in result.columns


def test_combine_bound_columns_empty_units():
    df = pd.DataFrame(
        {
            "lower": [5, 1],
            "upper": [5, 5],
            "units": ["", None],
        }
    )
    result = combine_bound_columns(
        df,
        lower_col="lower",
        upper_col="upper",
        units_col="units",
        out_col="combined",
    )
    assert list(result["combined"]) == ["5", "1-5"]


def test_combine_bound_columns_both_bounds_na():
    df = pd.DataFrame(
        {
            "lower": [np.nan],
            "upper": [np.nan],
        }
    )
    result = combine_bound_columns(
        df, lower_col="lower", upper_col="upper", out_col="combined"
    )
    assert list(result["combined"]) == [""]


def test_combine_bound_columns_missing_bound_column_returns_unchanged():
    df = pd.DataFrame({"lower": [1], "other": ["x"]})
    original = df.copy()
    result = combine_bound_columns(
        df, lower_col="lower", upper_col="upper", out_col="combined"
    )
    pd.testing.assert_frame_equal(result, original)
    assert "combined" not in result.columns


def test_combine_bound_columns_units_col_absent_still_combines():
    df = pd.DataFrame(
        {
            "lower": [1],
            "upper": [5],
        }
    )
    result = combine_bound_columns(
        df,
        lower_col="lower",
        upper_col="upper",
        units_col="units",
        out_col="combined",
    )
    assert list(result["combined"]) == ["1-5"]


def test_combine_bound_columns_drop_source_false_keeps_source_columns():
    df = pd.DataFrame(
        {
            "lower": [5],
            "upper": [5],
            "units": ["years"],
        }
    )
    result = combine_bound_columns(
        df,
        lower_col="lower",
        upper_col="upper",
        units_col="units",
        out_col="combined",
        drop_source=False,
    )
    assert list(result["combined"]) == ["5 years"]
    assert "lower" in result.columns
    assert "upper" in result.columns
    assert "units" in result.columns


def test_collapse_duplicate_columns_unique_unchanged():
    df = pd.DataFrame({"a": ["x"], "b": ["y"]})
    result = collapse_duplicate_columns(df)
    pd.testing.assert_frame_equal(result, df)


def test_collapse_duplicate_columns_same_value():
    df = pd.DataFrame([["x", "x"]], columns=["a", "a"])
    result = collapse_duplicate_columns(df)
    assert list(result.columns) == ["a"]
    assert result["a"].iloc[0] == "x"


def test_collapse_duplicate_columns_one_empty():
    df = pd.DataFrame([["x", ""], ["y", None]], columns=["a", "a"])
    result = collapse_duplicate_columns(df)
    assert list(result["a"]) == ["x", "y"]


def test_collapse_duplicate_columns_both_empty():
    df = pd.DataFrame([[None, ""]], columns=["a", "a"])
    result = collapse_duplicate_columns(df)
    assert pd.isna(result["a"].iloc[0])


def test_collapse_duplicate_columns_different_values():
    df = pd.DataFrame([["x", "y"]], columns=["a", "a"])
    result = collapse_duplicate_columns(df)
    assert result["a"].iloc[0] == ["x", "y"]


def test_collapse_duplicate_columns_nested_lists_flattened():
    df = pd.DataFrame([[["x"], ["y"]]], columns=["a", "a"])
    result = collapse_duplicate_columns(df)
    assert result["a"].iloc[0] == ["x", "y"]


def test_collapse_duplicate_columns_dedup_preserves_order():
    df = pd.DataFrame([["x", "x", "y"]], columns=["a", "a", "a"])
    result = collapse_duplicate_columns(df)
    assert result["a"].iloc[0] == ["x", "y"]


def test_collapse_duplicate_columns_mixed_with_unique():
    df = pd.DataFrame([["x", "y", "keep"]], columns=["a", "a", "b"])
    result = collapse_duplicate_columns(df)
    assert list(result.columns) == ["a", "b"]
    assert result["a"].iloc[0] == ["x", "y"]
    assert result["b"].iloc[0] == "keep"


def test_strip_author_metadata_column_prefix_renames():
    df = pd.DataFrame(
        {
            "tissues_author_metadata_mouse_litter_batch": ["A1"],
            "sample_name": ["s1"],
        }
    )
    result = strip_author_metadata_column_prefix(df)
    assert "mouse_litter_batch" in result.columns
    assert "tissues_author_metadata_mouse_litter_batch" not in result.columns
    assert list(result["mouse_litter_batch"]) == ["A1"]
    assert list(result["sample_name"]) == ["s1"]


def test_strip_author_metadata_column_prefix_noop_without_match():
    df = pd.DataFrame({"sample_name": ["s1"], "disease": ["normal"]})
    result = strip_author_metadata_column_prefix(df)
    pd.testing.assert_frame_equal(result, df)


def test_sort_ontology_term_id_column_reorders_paired_labels():
    df = pd.DataFrame(
        {
            "experimental_condition_ontology_term_id": [["EFO:2", "EFO:1"]],
            "experimental_condition": [["stim", "rest"]],
        }
    )
    result = sort_ontology_term_id_column(
        df, "experimental_condition_ontology_term_id"
    )
    assert result["experimental_condition_ontology_term_id"].iloc[0] == ["EFO:1", "EFO:2"]
    assert result["experimental_condition"].iloc[0] == ["rest", "stim"]


def test_sort_ontology_term_id_column_set_becomes_sorted_list():
    df = pd.DataFrame(
        {"self_reported_ethnicity_ontology_term_id": [{"HANCESTRO:2", "HANCESTRO:1"}]}
    )
    result = sort_ontology_term_id_column(
        df, "self_reported_ethnicity_ontology_term_id"
    )
    assert result["self_reported_ethnicity_ontology_term_id"].iloc[0] == [
        "HANCESTRO:1",
        "HANCESTRO:2",
    ]


def test_sort_ontology_term_id_column_no_paired_column():
    df = pd.DataFrame(
        {"experimental_condition_ontology_term_id": [["EFO:2", "EFO:1"]]}
    )
    result = sort_ontology_term_id_column(
        df, "experimental_condition_ontology_term_id"
    )
    assert result["experimental_condition_ontology_term_id"].iloc[0] == ["EFO:1", "EFO:2"]
    assert list(result.columns) == ["experimental_condition_ontology_term_id"]


def test_sort_ontology_term_id_column_missing_id_col_unchanged():
    df = pd.DataFrame({"sample_name": ["s1"]})
    result = sort_ontology_term_id_column(
        df, "experimental_condition_ontology_term_id"
    )
    pd.testing.assert_frame_equal(result, df)


def test_sort_ontology_term_id_column_length_mismatch_leaves_pair():
    df = pd.DataFrame(
        {
            "experimental_condition_ontology_term_id": [["EFO:2", "EFO:1"]],
            "experimental_condition": [["stim"]],
        }
    )
    result = sort_ontology_term_id_column(
        df, "experimental_condition_ontology_term_id"
    )
    assert result["experimental_condition_ontology_term_id"].iloc[0] == ["EFO:1", "EFO:2"]
    assert result["experimental_condition"].iloc[0] == ["stim"]


# --- controlled term splitting: missing values are None, not pd.NA ---


@pytest.mark.parametrize(
    'value',
    [
        pytest.param(None, id='none'),
        pytest.param([], id='empty_list'),
        pytest.param('not a dict', id='plain_string'),
        pytest.param([{'@id': '/controlled_terms/X:1/'}], id='dict_without_term_name'),
    ],
)
def test_split_term_cell_missing_values_are_none(value):
    """
    None rather than pd.NA: these values flow back into is_empty() via to_items()
    and single_or_list(), and `pd.NA == ""` evaluates to pd.NA rather than False,
    which raises "boolean value of NA is ambiguous".
    """
    term_id, term_name = split_term_cell(value)
    assert term_id is None
    assert term_name is None


def test_extract_controlled_term_id_empty_ref_is_none():
    assert extract_controlled_term_id('') is None


def test_collapse_dataframe_tolerates_missing_term_values():
    """
    Regression guard for the GEO build: a controlled term column populated for
    some rows and empty for others must still collapse.
    """
    empty = split_term_cell(None)[1]
    df = pd.DataFrame({'*library name': ['libA', 'libA'], '**tissue': ['lung', empty]})

    result = collapse_dataframe(df, group_col='*library name')

    assert result.iloc[0]['**tissue'] == 'lung'
