"""Tests for file_extract.scale_h5ad."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from file_extract.metadata_sheet import parse_sample_template
from file_extract.scale_flags import validate_raw_subdirs
from file_extract.scale_h5ad import (
    ScaleCram,
    _introspect_counts,
    _process_one,
    default_scale_h5ad_output_name,
    derived_from_by_filename,
    derived_from_filenames,
    derived_from_label,
    extract_scale_h5ad,
    file_belongs_to_sample,
    format_derived_from,
    format_feature_counts,
    format_samples_column,
    is_scale_h5ad_key,
    is_scale_mtx_key,
    leftover_cram_uris,
    leftover_cram_warning,
    empty_derived_from_warning,
    empty_raw_subdir_warning,
    is_cram_in_raw_search,
    list_raw_crams,
    parse_scale_cram_name,
    parse_samples_csv,
    raw_cram_search_prefixes,
    resolve_raw_subdir,
    sample_from_filename,
    tsv_filename,
    well_to_sample_map,
)
from file_extract.scale_wells import ScaleExtractError
from tests.file_extract_helpers import FIXTURES, MockS3Client

BUCKET = "czi-cro"
RUNDATE = "proj/ORD01/processed/Run_2001-01-13/"
OTHER = "proj/ORD01/processed/Run_1999-01-01/samples/OTHER.QSR-1_anndata.h5ad"
SAMPLES_CSV = FIXTURES / "scale_samples.csv"
SHEET_CSV = FIXTURES / "scale_sample_template.csv"

SAMP01_MERGED = f"{RUNDATE}samples/SAMP-01_anndata.h5ad"
SAMP01_QSR = f"{RUNDATE}samples/SAMP-01.QSR-1_anndata.h5ad"
SAMP01_QSR2 = f"{RUNDATE}samples/SAMP-01.QSR-2_anndata.h5ad"
SAMP02_QSR = f"{RUNDATE}samples/SAMP-02.QSR-1_anndata.h5ad"
CTRL_QSR = f"{RUNDATE}samples/CTRL-01.QSR-1_anndata.h5ad"
ALLCELLS = f"{RUNDATE}samples/SAMP-01.merged.allCells.csv"
CR_H5 = f"{RUNDATE}samples/sample_filtered_feature_bc_matrix.h5"
SAMP01_MTX_DIR = "SAMP-01.QSR-1-SCALEPLEX.filtered.matrix"
SAMP02_MTX_DIR = "SAMP-02.QSR-1-SCALEPLEX.filtered.matrix"
CTRL_MTX_DIR = "CTRL-01.QSR-1-SCALEPLEX.filtered.matrix"
SAMP01_MTX = f"{RUNDATE}scaleplex/{SAMP01_MTX_DIR}/matrix.mtx.gz"
SAMP02_MTX = f"{RUNDATE}scaleplex/{SAMP02_MTX_DIR}/matrix.mtx.gz"
CTRL_MTX = f"{RUNDATE}scaleplex/{CTRL_MTX_DIR}/matrix.mtx.gz"
RAW = "proj/ORD01/raw/426971/"
CRAM_SAMP01_GEX = f"{RAW}426971-RNA3-098C_GEX_QSR-1_1A.cram"
CRAM_SAMP01_GEX_B = f"{RAW}426971-RNA3-098C_GEX_QSR-1_2B.cram"
CRAM_SAMP02_GEX = f"{RAW}426971-RNA3-098C_GEX_QSR-1_3A.cram"
CRAM_SAMP01_PLX = f"{RAW}426971-RNA3-098C_hash_oligo_QSR-1-SCALEPLEX_1H.cram"
CRAM_SAMP02_PLX = f"{RAW}426971-RNA3-098C_hash_oligo_QSR-1-SCALEPLEX_3C.cram"
CRAM_CTRL = f"{RAW}426971-RNA3-098C_GEX_QSR-1_12H.cram"
CRAM_UNMATCHED = f"{RAW}426971-RNA3-098C_GEX_QSR-1_unmatched.cram"
CRAM_OTHER_QSR = f"{RAW}426971-RNA3-098C_GEX_QSR-2_1A.cram"
SAMP01_BARCODES = f"{RUNDATE}scaleplex/{SAMP01_MTX_DIR}/barcodes.tsv.gz"
NESTED_MTX = f"{RUNDATE}scaleplex/{SAMP01_MTX_DIR}/extra/matrix.mtx.gz"
OTHER_MTX = f"{RUNDATE}scaleplex/SAMP-01.noplex.filtered.matrix/matrix.mtx.gz"


def _sheet_text() -> str:
    return SHEET_CSV.read_text(encoding="utf-8")


def _samples_text() -> str:
    return SAMPLES_CSV.read_text(encoding="utf-8")


def test_sample_from_filename() -> None:
    assert sample_from_filename("SAMP-01.QSR-1_anndata.h5ad") == "SAMP-01"
    assert sample_from_filename("SAMP-01_anndata.h5ad") == "SAMP-01_anndata"


def test_file_belongs_to_sample() -> None:
    assert file_belongs_to_sample("SAMP-01_anndata.h5ad", "SAMP-01")
    assert file_belongs_to_sample("SAMP-01.QSR-1_anndata.h5ad", "SAMP-01")
    assert not file_belongs_to_sample("SAMP-02.QSR-1_anndata.h5ad", "SAMP-01")


def test_is_scale_h5ad_key() -> None:
    assert is_scale_h5ad_key(SAMP01_QSR, RUNDATE)
    assert not is_scale_h5ad_key(SAMP01_MERGED, RUNDATE)
    assert not is_scale_h5ad_key(ALLCELLS, RUNDATE)
    assert not is_scale_h5ad_key(CR_H5, RUNDATE)
    assert not is_scale_h5ad_key(OTHER, RUNDATE)
    assert not is_scale_h5ad_key(f"{RUNDATE}samples.csv", RUNDATE)


def test_is_scale_mtx_key() -> None:
    assert is_scale_mtx_key(SAMP01_MTX, RUNDATE)
    assert is_scale_mtx_key(f"{RUNDATE}scaleplex/{SAMP01_MTX_DIR}/matrix.mtx", RUNDATE)
    assert not is_scale_mtx_key(SAMP01_BARCODES, RUNDATE)
    assert not is_scale_mtx_key(NESTED_MTX, RUNDATE)
    assert not is_scale_mtx_key(OTHER_MTX, RUNDATE)
    assert not is_scale_mtx_key(SAMP01_QSR, RUNDATE)
    assert tsv_filename(SAMP01_MTX, RUNDATE) == f"{SAMP01_MTX_DIR}/matrix.mtx.gz"
    assert sample_from_filename(tsv_filename(SAMP01_MTX, RUNDATE)) == "SAMP-01"


def test_default_scale_h5ad_output_name() -> None:
    assert (
        default_scale_h5ad_output_name(RUNDATE) == "Run_2001-01-13_scale_h5ad_info.tsv"
    )


def test_parse_samples_csv() -> None:
    rows = parse_samples_csv(_samples_text())
    assert rows[0] == ("SAMP-01", "1A-2C")
    assert rows[-1] == ("CTRL-01", "12H")


def test_parse_sample_template_normalizes_wells() -> None:
    rows = parse_sample_template(_sheet_text())
    wells = {row["well"] for row in rows}
    assert "1A" in wells
    assert "2G" in wells
    assert "11A" not in wells
    assert not any("Plate only" in (row.get("RT_index") or "") for row in rows)


def test_parse_sample_template_skips_hash_comment_rows() -> None:
    csv_text = (
        "sample_name,RT_index\n"
        "# comment,Plate only: if multiple indices, list using a comma as separator\n"
        "tissue-1A,SCALEQUANT-A1\n"
    )
    rows = parse_sample_template(csv_text)
    assert [row["well"] for row in rows] == ["1A"]


def test_parse_sample_template_splits_comma_separated_rt_index() -> None:
    csv_text = (
        "sample_name,RT_index\n"
        'tissue-pool,"SCALEQUANT-A1,SCALEQUANT-B1,SCALEQUANT-C12"\n'
    )
    rows = parse_sample_template(csv_text)
    assert [row["well"] for row in rows] == ["1A", "1B", "12C"]
    assert {row["sample_name"] for row in rows} == {"tissue-pool"}


def test_parse_sample_template_invalid_data_row_still_errors() -> None:
    with pytest.raises(ScaleExtractError, match="SCALEQUANT-Z1"):
        parse_sample_template("sample_name,RT_index\ntissue-bad,SCALEQUANT-Z1\n")


def test_parse_sample_template_filters_by_experiment_name() -> None:
    csv_text = (
        "sample_name,RT_index,experiment_name\n"
        "keep-1A,SCALEQUANT-A1,RNA3_098\n"
        "drop-1B,SCALEQUANT-B1,CHEM13-R096\n"
        "keep-1C,SCALEQUANT-C1, RNA3_098 \n"
    )
    rows = parse_sample_template(csv_text, experiment="RNA3_098")
    assert [row["sample_name"] for row in rows] == ["keep-1A", "keep-1C"]
    assert {row["well"] for row in rows} == {"1A", "1C"}


def test_parse_sample_template_requires_experiment_name_when_filtering() -> None:
    with pytest.raises(ScaleExtractError, match="experiment_name"):
        parse_sample_template(
            "sample_name,RT_index\ntissue-1A,SCALEQUANT-A1\n",
            experiment="RNA3_098",
        )


def test_parse_sample_template_no_matching_experiment() -> None:
    csv_text = (
        "sample_name,RT_index,experiment_name\nother-1A,SCALEQUANT-A1,CHEM13-R096\n"
    )
    with pytest.raises(ScaleExtractError, match="RNA3_098"):
        parse_sample_template(csv_text, experiment="RNA3_098")


def test_format_samples_column_prefixes_lab() -> None:
    names = ("tissue-1A", "tissue-1B")
    assert json.loads(format_samples_column(names, "example-lab")) == [
        "example-lab:tissue-1A",
        "example-lab:tissue-1B",
    ]
    assert format_samples_column(names, "/labs/example-lab/") == format_samples_column(
        names, "example-lab"
    )


def test_validate_raw_subdirs_allows_names_and_s3_uris() -> None:
    assert validate_raw_subdirs([" 426971 ", "441969/"]) == ["426971", "441969"]
    assert validate_raw_subdirs(["426971,441969"]) == ["426971", "441969"]
    assert validate_raw_subdirs(["426971, 441969", "555000"]) == [
        "426971",
        "441969",
        "555000",
    ]
    assert validate_raw_subdirs(["s3://czi-novogene/proj/raw/426971/"]) == [
        "s3://czi-novogene/proj/raw/426971"
    ]
    assert validate_raw_subdirs(
        [
            "s3://czi-novogene/lab/order/RNA3_098,"
            "s3://czi-novogene/lab/order/CHEM13-R096/"
        ]
    ) == [
        "s3://czi-novogene/lab/order/RNA3_098",
        "s3://czi-novogene/lab/order/CHEM13-R096",
    ]
    with pytest.raises(ScaleExtractError, match="empty"):
        validate_raw_subdirs([""])
    with pytest.raises(ScaleExtractError, match="empty"):
        validate_raw_subdirs(["426971,"])
    with pytest.raises(ScaleExtractError, match="s3://"):
        validate_raw_subdirs(["proj/raw/426971"])


def test_validate_raw_subdirs_rejects_a_bucket_with_no_directory() -> None:
    """A bucket alone would be searched whole and match foreign crams."""
    for value in ["s3://novogene-delivery", "s3://novogene-delivery/"]:
        with pytest.raises(ScaleExtractError, match="name a directory"):
            validate_raw_subdirs([value])


def test_parse_scale_cram_name() -> None:
    gex = parse_scale_cram_name("426971-RNA3-098C_GEX_QSR-7_10A.cram")
    assert gex == ScaleCram(qsr="7", well="10A", scaleplex=False)
    hyphen = parse_scale_cram_name("440115-R115H_GEX_QSR-8-5B.cram")
    assert hyphen == ScaleCram(qsr="8", well="5B", scaleplex=False)
    plex = parse_scale_cram_name("426971-RNA3-098C_hash_oligo_QSR-7-SCALEPLEX_1E.cram")
    assert plex == ScaleCram(qsr="7", well="1E", scaleplex=True)
    plex_hyphen = parse_scale_cram_name(
        "440115-R115H_hash_oligo_QSR-8-SCALEPLEX-1A.cram"
    )
    assert plex_hyphen == ScaleCram(qsr="8", well="1A", scaleplex=True)
    assert parse_scale_cram_name("426971-RNA3-098C_GEX_QSR-7_unmatched.cram") is None
    assert (
        parse_scale_cram_name("426971-RNA3-098C_GEX_QSR-7_10C_unmatched.cram") is None
    )
    assert parse_scale_cram_name("notes.txt") is None
    assert parse_scale_cram_name("440115-R115H_GEX_QSR-8-13A.cram") is None


def test_derived_from_filenames() -> None:
    gex = ScaleCram(qsr="1", well="1A", scaleplex=False)
    plex = ScaleCram(qsr="1", well="1A", scaleplex=True)
    assert derived_from_filenames("SAMP-01", gex) == ("SAMP-01.QSR-1_anndata.h5ad",)
    assert derived_from_filenames("SAMP-01", plex) == (
        "SAMP-01.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx.gz",
        "SAMP-01.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx",
    )


def test_resolve_raw_subdir_uses_processed_sibling() -> None:
    assert resolve_raw_subdir("czi-cro", RUNDATE, "426971") == (
        "czi-cro",
        RAW,
    )
    assert resolve_raw_subdir(
        "czi-cro", RUNDATE, "s3://czi-novogene/lab/raw/441969"
    ) == ("czi-novogene", "lab/raw/441969/")
    assert resolve_raw_subdir(
        "czi-cro", RUNDATE, "s3://czi-novogene/lab/NVUS-04/RNA3_098"
    ) == ("czi-novogene", "lab/NVUS-04/RNA3_098/")


def test_raw_cram_search_prefixes_tries_raw_child_then_prefix() -> None:
    """The raw/ child is probed first; the prefix itself is the fallback."""
    assert raw_cram_search_prefixes("lab/NVUS-04/RNA3_098/") == (
        "lab/NVUS-04/RNA3_098/raw/",
        "lab/NVUS-04/RNA3_098/",
    )
    assert raw_cram_search_prefixes(RAW) == (f"{RAW}raw/", RAW)


def test_raw_cram_search_prefixes_does_not_probe_raw_twice() -> None:
    """A prefix already at the raw/ level needs no raw/raw/ probe."""
    assert raw_cram_search_prefixes("lab/NVUS-04/RNA3_098/raw") == (
        "lab/NVUS-04/RNA3_098/raw/",
    )
    assert raw_cram_search_prefixes("raw/") == ("raw/",)


def test_raw_cram_search_prefixes_tests_the_segment_not_the_suffix() -> None:
    """A name merely ending in "raw" is not itself the raw/ level."""
    assert raw_cram_search_prefixes("proj/ORD01/raw/ORD01_raw/") == (
        "proj/ORD01/raw/ORD01_raw/raw/",
        "proj/ORD01/raw/ORD01_raw/",
    )


def test_raw_cram_search_prefixes_never_yields_a_leading_slash() -> None:
    """An s3:// bucket root resolves to "", which is not a "/" prefix."""
    assert raw_cram_search_prefixes("/proj/ORD01/") == (
        "proj/ORD01/raw/",
        "proj/ORD01/",
    )


def test_raw_cram_search_prefixes_refuses_a_bucket_root() -> None:
    """Even raw/ at a bucket root spans deliveries, so there is no candidate."""
    for prefix in ["", "/"]:
        with pytest.raises(ScaleExtractError, match="not a bucket root"):
            raw_cram_search_prefixes(prefix)


def test_raw_cram_search_prefixes_covers_both_non_numeric_layouts() -> None:
    """No segment name separates these two shapes, so both must be tried.

    A bare name resolves under raw/, and an s3:// group directory can
    itself sit under a top-level raw/. Each keeps its own layout.
    """
    bare = resolve_raw_subdir("czi-cro", RUNDATE, "RNA3_098")[1]
    assert bare == "proj/ORD01/raw/RNA3_098/"
    assert bare in raw_cram_search_prefixes(bare)

    group = resolve_raw_subdir(
        "czi-cro", RUNDATE, "s3://novogene-delivery/raw/ORD01/RNA3_098"
    )[1]
    assert group == "raw/ORD01/RNA3_098/"
    assert f"{group}raw/" in raw_cram_search_prefixes(group)


def test_is_cram_in_raw_search_requires_numeric_folder() -> None:
    group_raw = "lab/NVUS-04/RNA3_098/raw/"
    assert is_cram_in_raw_search(f"{group_raw}426971/file.cram", group_raw)
    assert not is_cram_in_raw_search(f"{group_raw}file.cram", group_raw)
    assert not is_cram_in_raw_search(f"{group_raw}notes/file.cram", group_raw)
    assert is_cram_in_raw_search(CRAM_SAMP01_GEX, RAW)
    assert not is_cram_in_raw_search(CRAM_UNMATCHED, RAW)


def test_empty_derived_from_warning_names_the_file() -> None:
    """Literal text, because the docs promise the filename and --raw-subdirs."""
    assert empty_derived_from_warning("SAMP-01.QSR-2_anndata.h5ad") == (
        "derived_from is empty for 'SAMP-01.QSR-2_anndata.h5ad'; "
        "are any S3 raw cram directories missing from --raw-subdirs?"
    )


def test_empty_raw_subdir_warning_names_the_entry_and_every_prefix() -> None:
    """Literal text, because the docs promise the entry and all prefixes."""
    assert empty_raw_subdir_warning("999999", ("s3://b/p/raw/", "s3://b/p/")) == (
        "--raw-subdirs '999999' matched no crams under "
        "s3://b/p/raw/ or s3://b/p/; "
        "*.cram must sit in a numeric run folder -- one of those prefixes "
        "itself, or a folder directly under one"
    )
    assert empty_raw_subdir_warning("RNA3_098", ("s3://b/g/raw/",)) == (
        "--raw-subdirs 'RNA3_098' matched no crams under s3://b/g/raw/; "
        "*.cram must sit in a numeric run folder -- one of those prefixes "
        "itself, or a folder directly under one"
    )


def test_list_raw_crams_finds_numeric_children_of_group_uri() -> None:
    group = "lab/NVUS-04/RNA3_098/"
    nested = f"{group}raw/426971/426971-RNA3-098C_GEX_QSR-1_1A.cram"
    stray = f"{group}file.cram"
    client = MockS3Client(keys=[nested, stray, CRAM_SAMP01_GEX])
    found = list_raw_crams(
        client,
        BUCKET,
        RUNDATE,
        ["s3://czi-novogene/lab/NVUS-04/RNA3_098"],
    )
    assert found.crams == [("czi-novogene", nested)]
    assert found.empty_subdirs == []
    # The raw/ candidate hit, so the fallback prefix must not be walked.
    assert client.paginate_calls == 1


def test_list_raw_crams_walks_the_fallback_only_after_a_miss(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The fallback costs a second listing, so it is only reached on a miss.

    The debug line must name the prefix that actually supplied the crams,
    which here is the fallback rather than the raw/ child that missed.
    """
    nested = "proj/ORD01/raw/RNA3_098/426971/426971-RNA3-098C_GEX_QSR-1_1A.cram"
    client = MockS3Client(keys=[nested])
    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger="file_extract.scale_h5ad"):
        found = list_raw_crams(client, BUCKET, RUNDATE, ["RNA3_098"])
    assert found.crams == [(BUCKET, nested)]
    assert client.paginate_calls == 2
    assert f"1 crams under s3://{BUCKET}/proj/ORD01/raw/RNA3_098/" in caplog.text
    assert "RNA3_098/raw/" not in caplog.text


def test_list_raw_crams_logs_nothing_when_a_subdir_is_empty(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """No prefix supplied crams, so there is no winner to name."""
    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger="file_extract.scale_h5ad"):
        found = list_raw_crams(MockS3Client(keys=[]), BUCKET, RUNDATE, ["999999"])
    assert found.crams == []
    assert "crams under" not in caplog.text


def test_list_raw_crams_finds_numeric_children_of_bare_name() -> None:
    """A bare non-numeric name resolves to raw/{name}/{numeric}/."""
    nested = "proj/ORD01/raw/RNA3_098/426971/426971-RNA3-098C_GEX_QSR-1_1A.cram"
    found = list_raw_crams(MockS3Client(keys=[nested]), BUCKET, RUNDATE, ["RNA3_098"])
    assert found.crams == [(BUCKET, nested)]
    assert found.empty_subdirs == []


def test_list_raw_crams_finds_crams_under_a_raw_suffixed_name() -> None:
    """A subdir named e.g. ORD01_raw still gets its raw/ child searched."""
    nested = "proj/ORD01/raw/ORD01_raw/raw/426971/426971-R_GEX_QSR-1_1A.cram"
    found = list_raw_crams(MockS3Client(keys=[nested]), BUCKET, RUNDATE, ["ORD01_raw"])
    assert found.crams == [(BUCKET, nested)]
    assert found.empty_subdirs == []


def test_list_raw_crams_finds_numeric_children_under_top_level_raw() -> None:
    """A group dir that itself sits under raw/ still gets its raw/ level."""
    group = "raw/ORD01/RNA3_098/"
    nested = f"{group}raw/426971/426971-RNA3-098C_GEX_QSR-1_1A.cram"
    found = list_raw_crams(
        MockS3Client(keys=[nested]),
        BUCKET,
        RUNDATE,
        ["s3://novogene-delivery/raw/ORD01/RNA3_098"],
    )
    assert found.crams == [("novogene-delivery", nested)]
    assert found.empty_subdirs == []


def test_list_raw_crams_overlapping_subdirs_are_not_reported_empty() -> None:
    """A subdir whose crams an earlier entry claimed still counts as matched.

    Both values resolve to the same prefix, so the second finds only
    duplicates. The match counter therefore has to run before the dedup
    check -- counting after it would call the second entry empty and, under
    --strict, fail the run.
    """
    client = MockS3Client(keys=[CRAM_SAMP01_GEX])
    found = list_raw_crams(
        client,
        BUCKET,
        RUNDATE,
        ["426971", f"s3://{BUCKET}/proj/ORD01/raw/426971"],
    )
    assert found.crams == [(BUCKET, CRAM_SAMP01_GEX)]
    assert found.empty_subdirs == []


def test_list_raw_crams_reports_subdirs_that_matched_nothing() -> None:
    client = MockS3Client(keys=[CRAM_SAMP01_GEX])
    found = list_raw_crams(client, BUCKET, RUNDATE, ["426971", "999999"])
    assert found.crams == [(BUCKET, CRAM_SAMP01_GEX)]
    assert found.empty_subdirs == [
        (
            "999999",
            (
                f"s3://{BUCKET}/proj/ORD01/raw/999999/raw/",
                f"s3://{BUCKET}/proj/ORD01/raw/999999/",
            ),
        )
    ]


def test_well_to_sample_map_expands_barcodes() -> None:
    owners = well_to_sample_map([("SAMP-01", "1A-2C"), ("SAMP-02", "3A-3G")])
    assert owners["1A"] == "SAMP-01"
    assert owners["2C"] == "SAMP-01"
    assert owners["3A"] == "SAMP-02"
    assert "2D" not in owners


def test_derived_from_label_prefixes_basename() -> None:
    assert derived_from_label("example-lab", CRAM_SAMP01_GEX) == (
        "example-lab:426971-RNA3-098C_GEX_QSR-1_1A.cram"
    )
    assert derived_from_label("/labs/example-lab/", "file.cram") == (
        "example-lab:file.cram"
    )


def test_derived_from_by_filename_groups_crams() -> None:
    owners = well_to_sample_map(_parse_fixture_samples())
    grouped = derived_from_by_filename(
        [
            (BUCKET, CRAM_SAMP01_GEX),
            (BUCKET, CRAM_SAMP01_GEX_B),
            (BUCKET, CRAM_SAMP02_GEX),
            (BUCKET, CRAM_SAMP01_PLX),
            (BUCKET, CRAM_CTRL),
            (BUCKET, CRAM_UNMATCHED),
            (BUCKET, CRAM_OTHER_QSR),
        ],
        owners,
        {"CTRL-01"},
        "example-lab",
    )
    assert grouped["SAMP-01.QSR-1_anndata.h5ad"] == [
        derived_from_label("example-lab", CRAM_SAMP01_GEX),
        derived_from_label("example-lab", CRAM_SAMP01_GEX_B),
    ]
    assert grouped["SAMP-02.QSR-1_anndata.h5ad"] == [
        derived_from_label("example-lab", CRAM_SAMP02_GEX)
    ]
    assert grouped[f"{SAMP01_MTX_DIR}/matrix.mtx.gz"] == [
        derived_from_label("example-lab", CRAM_SAMP01_PLX)
    ]
    assert "SAMP-01.QSR-2_anndata.h5ad" in grouped
    assert "CTRL-01.QSR-1_anndata.h5ad" not in grouped
    leftover = leftover_cram_uris(
        [
            (BUCKET, CRAM_SAMP01_GEX),
            (BUCKET, CRAM_CTRL),
            (BUCKET, CRAM_UNMATCHED),
            (BUCKET, CRAM_OTHER_QSR),
        ],
        grouped,
        {"SAMP-01.QSR-1_anndata.h5ad"},
        "example-lab",
    )
    assert leftover == [
        f"s3://{BUCKET}/{CRAM_CTRL}",
        f"s3://{BUCKET}/{CRAM_UNMATCHED}",
        f"s3://{BUCKET}/{CRAM_OTHER_QSR}",
    ]
    assert leftover_cram_warning(leftover[0]).startswith("cram s3://")


def test_format_derived_from() -> None:
    assert json.loads(format_derived_from([])) == []
    assert json.loads(
        format_derived_from(
            [
                "example-lab:441479-R096A_GEX_QSR-1-10A.cram",
                "example-lab:441479-R096A_GEX_QSR-1-10B.cram",
            ]
        )
    ) == [
        "example-lab:441479-R096A_GEX_QSR-1-10A.cram",
        "example-lab:441479-R096A_GEX_QSR-1-10B.cram",
    ]


def _parse_fixture_samples() -> list[tuple[str, str]]:
    return parse_samples_csv(_samples_text())


def _h5ad_dims(bucket: str, key: str) -> tuple[int, int]:
    counts = {
        SAMP01_QSR: (11, 18129),
        SAMP01_QSR2: (15, 18129),
        SAMP02_QSR: (22, 900),
        CTRL_QSR: (3, 4),
    }
    if key not in counts:
        raise RuntimeError(f"unexpected h5ad key {key}")
    return counts[key]


def _mtx_dims(s3_client: object, bucket: str, key: str) -> tuple[int, int]:
    counts = {SAMP01_MTX: (101, 8), SAMP02_MTX: (202, 12), CTRL_MTX: (3, 2)}
    if key not in counts:
        raise RuntimeError(f"unexpected mtx key {key}")
    return counts[key]


def test_extract_scale_h5ad_writes_non_control_rows(tmp_path: Path) -> None:
    keys = [
        f"{RUNDATE}samples.csv",
        SAMP01_MERGED,
        SAMP01_QSR,
        SAMP02_QSR,
        CTRL_QSR,
        ALLCELLS,
        CR_H5,
        OTHER,
        SAMP01_MTX,
        SAMP02_MTX,
        CTRL_MTX,
        SAMP01_BARCODES,
        NESTED_MTX,
        OTHER_MTX,
        CRAM_SAMP01_GEX,
        CRAM_SAMP01_GEX_B,
        CRAM_SAMP02_GEX,
        CRAM_SAMP01_PLX,
        CRAM_SAMP02_PLX,
        CRAM_CTRL,
        CRAM_UNMATCHED,
        CRAM_OTHER_QSR,
    ]
    client = MockS3Client(
        keys=keys,
        sizes={
            SAMP01_QSR: 111,
            SAMP02_QSR: 222,
            CTRL_QSR: 33,
            SAMP01_MTX: 1001,
            SAMP02_MTX: 2002,
            CTRL_MTX: 33,
        },
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={
            SAMP01_MERGED: "crc-merged",
            SAMP01_QSR: "crc-qsr1",
            SAMP02_QSR: "crc-qsr2",
            CTRL_QSR: "crc-ctrl",
            SAMP01_MTX: "crc-mtx1",
            SAMP02_MTX: "crc-mtx2",
            CTRL_MTX: "crc-mtx-ctrl",
        },
    )
    out = tmp_path / "out.tsv"
    with (
        patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims),
        patch("file_extract.scale_h5ad.count_mtx_dims", side_effect=_mtx_dims),
    ):
        summary = extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )

    assert summary.total == 4
    assert summary.crc_ok == 4
    assert summary.enrichment_ok == 4
    assert any("CTRL-01" in warning for warning in summary.warnings)
    leftover = [w for w in summary.warnings if w.startswith("cram s3://")]
    leftover_uris = " ".join(leftover)
    assert CRAM_CTRL in leftover_uris
    assert CRAM_OTHER_QSR in leftover_uris
    assert CRAM_UNMATCHED not in leftover_uris
    assert CRAM_SAMP01_GEX not in leftover_uris
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert [row["filename"] for row in rows] == [
        "SAMP-01.QSR-1_anndata.h5ad",
        "SAMP-02.QSR-1_anndata.h5ad",
        f"{SAMP01_MTX_DIR}/matrix.mtx.gz",
        f"{SAMP02_MTX_DIR}/matrix.mtx.gz",
    ]
    assert rows[0]["sample"] == "SAMP-01"
    assert rows[0]["s3_uri"] == f"s3://{BUCKET}/{SAMP01_QSR}"
    assert rows[0]["crc64nvme_base64"] == "crc-qsr1"
    assert rows[0]["file_size"] == "111"
    assert rows[0]["observation_count"] == "11"
    assert json.loads(rows[0]["feature_counts"]) == [
        {"feature_type": "gene", "feature_count": 18129}
    ]
    assert rows[1]["crc64nvme_base64"] == "crc-qsr2"
    assert rows[1]["file_size"] == "222"
    assert rows[1]["observation_count"] == "22"
    assert json.loads(rows[1]["feature_counts"]) == [
        {"feature_type": "gene", "feature_count": 900}
    ]
    assert rows[2]["sample"] == "SAMP-01"
    assert rows[2]["s3_uri"] == f"s3://{BUCKET}/{SAMP01_MTX}"
    assert rows[2]["crc64nvme_base64"] == "crc-mtx1"
    assert rows[2]["file_size"] == "1001"
    assert rows[2]["observation_count"] == "101"
    assert json.loads(rows[2]["feature_counts"]) == [
        {"feature_type": "hash oligo", "feature_count": 8}
    ]
    assert rows[3]["sample"] == "SAMP-02"
    assert rows[3]["crc64nvme_base64"] == "crc-mtx2"
    assert rows[3]["file_size"] == "2002"
    assert rows[3]["observation_count"] == "202"
    assert json.loads(rows[3]["feature_counts"]) == [
        {"feature_type": "hash oligo", "feature_count": 12}
    ]
    assert all(row["crc64nvme_base64"] for row in rows)
    assert all(row["observation_count"] for row in rows)
    assert all(row["feature_counts"] for row in rows)
    assert json.loads(rows[0]["samples"]) == [
        "example-lab:tissue-1A",
        "example-lab:tissue-1B",
        "example-lab:tissue-1C",
        "example-lab:tissue-1D",
        "example-lab:tissue-1E",
        "example-lab:tissue-1F",
        "example-lab:tissue-1G",
        "example-lab:tissue-1H",
        "example-lab:tissue-2A",
        "example-lab:tissue-2B",
        "example-lab:tissue-2C",
    ]
    assert json.loads(rows[1]["samples"]) == [
        "example-lab:tissue-3A",
        "example-lab:tissue-3B",
        "example-lab:tissue-3C",
        "example-lab:tissue-3D",
        "example-lab:tissue-3E",
        "example-lab:tissue-3F",
        "example-lab:tissue-3G",
    ]
    assert json.loads(rows[2]["samples"]) == json.loads(rows[0]["samples"])
    assert json.loads(rows[3]["samples"]) == json.loads(rows[1]["samples"])
    assert json.loads(rows[0]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP01_GEX),
        derived_from_label("example-lab", CRAM_SAMP01_GEX_B),
    ]
    assert json.loads(rows[1]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP02_GEX)
    ]
    assert json.loads(rows[2]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP01_PLX)
    ]
    assert json.loads(rows[3]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP02_PLX)
    ]
    assert all(row["filename"] != "SAMP-01_anndata.h5ad" for row in rows)
    assert all(row["filename"] != "CTRL-01.QSR-1_anndata.h5ad" for row in rows)
    assert all(CTRL_MTX_DIR not in row["filename"] for row in rows)


def test_extract_scale_h5ad_derived_from_merges_raw_subdirs(tmp_path: Path) -> None:
    other_raw = "proj/ORD01/raw/441969/"
    other_cram = f"{other_raw}441969-RNA3-098C_GEX_QSR-1_1C.cram"
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR, CRAM_SAMP01_GEX, other_cram],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1"},
    )
    out = tmp_path / "merged.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971", "441969"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert json.loads(rows[0]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP01_GEX),
        derived_from_label("example-lab", other_cram),
    ]


def test_extract_scale_h5ad_derived_from_group_uri_numeric_raw(
    tmp_path: Path,
) -> None:
    group = "lab/NVUS-04/RNA3_098/"
    nested = f"{group}raw/426971/426971-RNA3-098C_GEX_QSR-1_1A.cram"
    stray = f"{group}file.cram"
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR, nested, stray],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1"},
    )
    out = tmp_path / "group.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["s3://czi-novogene/lab/NVUS-04/RNA3_098"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert json.loads(rows[0]["derived_from"]) == [
        derived_from_label("example-lab", nested)
    ]


def test_extract_scale_h5ad_warns_when_raw_subdir_matched_nothing(
    tmp_path: Path,
) -> None:
    """An empty derived_from is reported, not written out silently."""
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1"},
    )
    out = tmp_path / "no_crams.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        summary = extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["RNA3_098"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )
    assert summary.empty_raw_subdirs == ["RNA3_098"]
    warning = next(w for w in summary.warnings if w.startswith("--raw-subdirs"))
    # Equality, not substrings: the fallback uri is a substring of the raw/
    # one, so "uri in warning" is satisfied by the raw/ uri alone and pins
    # neither the subdir name nor that both prefixes are named.
    assert warning == empty_raw_subdir_warning(
        "RNA3_098",
        (
            f"s3://{BUCKET}/proj/ORD01/raw/RNA3_098/raw/",
            f"s3://{BUCKET}/proj/ORD01/raw/RNA3_098/",
        ),
    )
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert json.loads(rows[0]["derived_from"]) == []
    assert empty_derived_from_warning("SAMP-01.QSR-1_anndata.h5ad") in summary.warnings


def test_extract_scale_h5ad_warns_when_one_derived_from_is_empty(
    tmp_path: Path,
) -> None:
    """A QSR-2 row with no matching CRAM is named; filled QSR-1 rows are not."""
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR, SAMP01_QSR2, CRAM_SAMP01_GEX],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1", SAMP01_QSR2: "crc-qsr2"},
    )
    out = tmp_path / "partial_derived.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        summary = extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )
    empty_warnings = [
        w for w in summary.warnings if w.startswith("derived_from is empty")
    ]
    assert empty_warnings == [empty_derived_from_warning("SAMP-01.QSR-2_anndata.h5ad")]
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    by_name = {row["filename"]: row for row in rows}
    assert json.loads(by_name["SAMP-01.QSR-1_anndata.h5ad"]["derived_from"]) == [
        derived_from_label("example-lab", CRAM_SAMP01_GEX)
    ]
    assert json.loads(by_name["SAMP-01.QSR-2_anndata.h5ad"]["derived_from"]) == []


def test_format_feature_counts() -> None:
    assert json.loads(format_feature_counts("gene", 18129)) == [
        {"feature_type": "gene", "feature_count": 18129}
    ]


def test_process_one_gathers_crc_and_obs_for_h5ad_and_mtx() -> None:
    client = MockS3Client(crc_by_key={SAMP01_QSR: "crc-h5ad", SAMP01_MTX: "crc-mtx"})
    with (
        patch("file_extract.scale_h5ad.count_h5ad_dims", return_value=(11, 18129)),
        patch("file_extract.scale_h5ad.count_mtx_dims", return_value=(101, 8)),
    ):
        h5ad_crc, h5ad_crc_err, h5ad_obs, h5ad_fc, h5ad_err = _process_one(
            client, BUCKET, SAMP01_QSR, retries=1
        )
        mtx_crc, mtx_crc_err, mtx_obs, mtx_fc, mtx_err = _process_one(
            client, BUCKET, SAMP01_MTX, retries=1
        )
    assert (h5ad_crc, h5ad_crc_err, h5ad_obs, h5ad_err) == ("crc-h5ad", "", 11, "")
    assert json.loads(h5ad_fc or "") == [
        {"feature_type": "gene", "feature_count": 18129}
    ]
    assert (mtx_crc, mtx_crc_err, mtx_obs, mtx_err) == ("crc-mtx", "", 101, "")
    assert json.loads(mtx_fc or "") == [
        {"feature_type": "hash oligo", "feature_count": 8}
    ]
    with pytest.raises(RuntimeError, match="observation counter"):
        _introspect_counts(client, BUCKET, f"{RUNDATE}samples/notes.txt")


def test_extract_scale_h5ad_zero_matches_does_not_write(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", CTRL_QSR],
        object_bodies={f"{RUNDATE}samples.csv": "sample,barcodes\nCTRL-01,12H\n"},
        crc_by_key={CTRL_QSR: "crc-ctrl"},
    )
    out = tmp_path / "empty.tsv"
    summary = extract_scale_h5ad(
        client,
        BUCKET,
        RUNDATE,
        str(out),
        metadata_gid="sheet-uuid",
        metadata_experiment="RNA3_098",
        lab="example-lab",
        raw_subdirs=["426971"],
        show_progress=False,
        sheet_csv=_sheet_text(),
    )
    assert summary.total == 0
    assert not out.exists()


def test_extract_scale_h5ad_uses_fetch_sheet(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1"},
    )
    fetched: list[str] = []

    def fetch_sheet(sheet_id: str) -> str:
        fetched.append(sheet_id)
        return _sheet_text()

    out = tmp_path / "fetched.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="/labs/example-lab/",
            raw_subdirs=["426971"],
            show_progress=False,
            fetch_sheet=fetch_sheet,
        )
    assert fetched == ["sheet-uuid"]
    assert out.exists()
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert all(
        value.startswith("example-lab:") for value in json.loads(rows[0]["samples"])
    )


def test_extract_scale_h5ad_ignores_other_experiment_rows(tmp_path: Path) -> None:
    sheet = (
        _sheet_text().rstrip()
        + "\n"
        + "other-1A,SCALEQUANT-A1,CHEM13-R096\n"
        + "other-ctrl,SCALEQUANT-H12,CHEM13-R096\n"
    )
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR, CTRL_QSR],
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1", CTRL_QSR: "crc-ctrl"},
    )
    out = tmp_path / "filtered.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=_h5ad_dims):
        summary = extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971"],
            show_progress=False,
            sheet_csv=sheet,
        )
    assert any("CTRL-01" in warning for warning in summary.warnings)
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert [row["filename"] for row in rows] == ["SAMP-01.QSR-1_anndata.h5ad"]
    samples = json.loads(rows[0]["samples"])
    assert "example-lab:tissue-1A" in samples
    assert "example-lab:other-1A" not in samples
    assert "example-lab:other-ctrl" not in samples


def test_extract_scale_h5ad_obs_failure_leaves_count_empty(tmp_path: Path) -> None:
    client = MockS3Client(
        keys=[f"{RUNDATE}samples.csv", SAMP01_QSR, SAMP02_QSR],
        sizes={SAMP01_QSR: 111, SAMP02_QSR: 222},
        object_bodies={f"{RUNDATE}samples.csv": _samples_text()},
        crc_by_key={SAMP01_QSR: "crc-qsr1", SAMP02_QSR: "crc-qsr2"},
    )

    def count_dims(bucket: str, key: str) -> tuple[int, int]:
        if key == SAMP01_QSR:
            raise RuntimeError("No 'obs' group or dataset; not an AnnData h5ad")
        return 22, 900

    out = tmp_path / "partial.tsv"
    with patch("file_extract.scale_h5ad.count_h5ad_dims", side_effect=count_dims):
        summary = extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            str(out),
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )

    assert summary.total == 2
    assert summary.crc_ok == 2
    assert summary.enrichment_ok == 1
    assert len(summary.failures) == 1
    assert summary.failures[0][0] == SAMP01_QSR
    assert "obs" in summary.failures[0][2]
    rows = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    by_name = {row["filename"]: row for row in rows}
    assert by_name["SAMP-01.QSR-1_anndata.h5ad"]["observation_count"] == ""
    assert by_name["SAMP-01.QSR-1_anndata.h5ad"]["feature_counts"] == ""
    assert by_name["SAMP-01.QSR-1_anndata.h5ad"]["file_size"] == "111"
    assert by_name["SAMP-02.QSR-1_anndata.h5ad"]["observation_count"] == "22"
    assert json.loads(by_name["SAMP-02.QSR-1_anndata.h5ad"]["feature_counts"]) == [
        {"feature_type": "gene", "feature_count": 900}
    ]


def test_extract_scale_h5ad_missing_samples_csv() -> None:
    client = MockS3Client(keys=[SAMP01_QSR], crc_by_key={SAMP01_QSR: "crc"})
    with pytest.raises(ScaleExtractError, match="samples.csv"):
        extract_scale_h5ad(
            client,
            BUCKET,
            RUNDATE,
            "unused.tsv",
            metadata_gid="sheet-uuid",
            metadata_experiment="RNA3_098",
            lab="example-lab",
            raw_subdirs=["426971"],
            show_progress=False,
            sheet_csv=_sheet_text(),
        )


def test_metadata_sheet_parse_requires_rt_index() -> None:
    with pytest.raises(ScaleExtractError, match="RT_index"):
        parse_sample_template("sample_name\nfoo\n")
