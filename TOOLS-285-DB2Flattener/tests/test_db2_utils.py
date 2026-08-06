import numpy as np
import pandas as pd

from DB2_utils import (
    collapse_duplicate_columns,
    combine_bound_columns,
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
