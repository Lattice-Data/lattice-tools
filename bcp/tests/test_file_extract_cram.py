"""Tests for file_extract.cram."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from file_extract.cram import (
    CramListingWarnings,
    default_cram_output_name,
    extract_cram,
    is_target_file,
    is_unmatched_cram,
    only_unmatched_cram_warning,
    scan_cram_listing,
    ucram_found_warning,
)
from tests.file_extract_helpers import FIXTURES, MockS3Client

PREFIX = "test-order/NVUS0000000000-01/"
BUCKET = "example-bucket"
MATCHED_KEY = (
    f"{PREFIX}BO_AO1_3p_041/raw/442383-BO_AO1_3p_041_GEX-Z0077-CAATGACTATCTGAT.cram"
)
UNMATCHED_RAW_KEY = (
    f"{PREFIX}BO_AO1_3p_041/raw/"
    "442383-BO_AO1_3p_041_GEX-Z0077-CAATGACTATCTGAT_unmatched.cram"
)
UNMATCHED_TRIMMED_KEY = (
    f"{PREFIX}BO_AO1_3p_041/trimmed/"
    "442383-BO_AO1_3p_041_GEX-Z0077-CAATGACTATCTGAT_unmatched.cram"
)
UCRAM_KEY = f"{PREFIX}BO_AO1_3p_041/raw/sample.ucram"


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        (UNMATCHED_RAW_KEY, True),
        (UNMATCHED_TRIMMED_KEY, True),
        (MATCHED_KEY, False),
        (UCRAM_KEY, False),
    ],
)
def test_is_unmatched_cram(key: str, expected: bool) -> None:
    assert is_unmatched_cram(key) is expected


@pytest.mark.parametrize(
    ("key", "require_raw", "expected"),
    [
        (MATCHED_KEY, True, True),
        (MATCHED_KEY.replace("/raw/", "/processed/"), True, False),
        (UNMATCHED_RAW_KEY, True, False),
        (UNMATCHED_TRIMMED_KEY, True, False),
        (MATCHED_KEY.replace("/raw/", "/other/"), False, True),
        (UCRAM_KEY, True, False),
    ],
)
def test_is_target_file(key: str, require_raw: bool, expected: bool) -> None:
    assert is_target_file(key, require_raw=require_raw) is expected


def test_default_cram_output_name() -> None:
    assert default_cram_output_name("order/prefix/") == "prefix_cram_info.tsv"


def test_warning_helpers() -> None:
    assert only_unmatched_cram_warning(0) is None
    msg = only_unmatched_cram_warning(3)
    assert msg is not None
    assert "3 unmatched" in msg

    assert ucram_found_warning(0) is None
    ucram_msg = ucram_found_warning(2)
    assert ucram_msg is not None
    assert ".ucram" in ucram_msg


def test_only_unmatched_warning_drops_raw_clause_when_not_required() -> None:
    """--no-require-raw never restricted the search, so it must not be claimed."""
    assert "under /raw/" in str(only_unmatched_cram_warning(1, require_raw=True))
    relaxed = only_unmatched_cram_warning(1, require_raw=False)
    assert relaxed is not None
    assert "/raw/" not in relaxed
    assert "1 unmatched" in relaxed


def test_scan_cram_listing_classifies_in_one_pass() -> None:
    keys = [MATCHED_KEY, UNMATCHED_TRIMMED_KEY, UCRAM_KEY]
    client = MockS3Client(keys=keys, sizes={k: 100 for k in keys})
    listing = scan_cram_listing(client, BUCKET, PREFIX)
    assert listing.warnings == CramListingWarnings(
        unmatched_cram_count=1, ucram_count=1
    )
    assert [obj.key for obj in listing.targets] == [MATCHED_KEY]
    assert listing.targets[0].size_bytes == 100


def test_scan_cram_listing_honours_require_raw() -> None:
    outside_raw = MATCHED_KEY.replace("/raw/", "/other/")
    client = MockS3Client(keys=[outside_raw], sizes={outside_raw: 5})
    assert scan_cram_listing(client, BUCKET, PREFIX).targets == ()
    relaxed = scan_cram_listing(client, BUCKET, PREFIX, require_raw=False)
    assert [obj.key for obj in relaxed.targets] == [outside_raw]


def test_scan_cram_listing_walks_prefix_once() -> None:
    """The double walk this replaced cost one full listing per concern."""
    keys = [MATCHED_KEY, UCRAM_KEY]
    client = MockS3Client(keys=keys, sizes={k: 1 for k in keys})
    scan_cram_listing(client, BUCKET, PREFIX)
    assert client.paginate_calls == 1


def test_extract_cram_zero_matches_writes_nothing(tmp_path: Path) -> None:
    out = tmp_path / "out.tsv"
    summary = extract_cram(
        MockS3Client(), BUCKET, PREFIX, str(out), show_progress=False, inline=True
    )
    assert summary.total == 0
    assert not out.exists()


def test_extract_cram_records_failed_enrichment(tmp_path: Path) -> None:
    """A CRAM with no companion metadata still yields a row, with the error noted."""
    client = MockS3Client(
        keys=[MATCHED_KEY],
        sizes={MATCHED_KEY: 100},
        crc_by_key={MATCHED_KEY: "crc"},
        object_bodies={},
    )
    summary = extract_cram(
        client,
        BUCKET,
        PREFIX,
        str(tmp_path / "out.tsv"),
        show_progress=False,
        inline=True,
    )
    assert summary.crc_ok == 1
    assert summary.enrichment_ok == 0
    assert len(summary.failures) == 1
    assert summary.failures[0][0] == MATCHED_KEY


def test_extract_cram_scans_once_when_given_no_listing(tmp_path: Path) -> None:
    """Standalone callers keep working; the listing argument is an optimisation."""
    meta_body = (FIXTURES / "cram_metadata.json").read_text(encoding="utf-8")
    client = MockS3Client(
        keys=[MATCHED_KEY],
        sizes={MATCHED_KEY: 100},
        crc_by_key={MATCHED_KEY: "crc"},
        object_bodies={MATCHED_KEY + "-metadata.json": meta_body},
    )
    summary = extract_cram(
        client,
        BUCKET,
        PREFIX,
        str(tmp_path / "out.tsv"),
        show_progress=False,
        inline=True,
    )
    assert summary.total == 1
    assert client.paginate_calls == 1


def test_extract_cram_integration(tmp_path: Path) -> None:
    meta_body = (FIXTURES / "cram_metadata.json").read_text(encoding="utf-8")
    keys = [MATCHED_KEY, UNMATCHED_TRIMMED_KEY, UCRAM_KEY]
    sizes = {MATCHED_KEY: 2_000_000_000_000}
    crc = {MATCHED_KEY: "crc-matched"}
    bodies = {MATCHED_KEY + "-metadata.json": meta_body}
    client = MockS3Client(
        keys=keys,
        sizes=sizes,
        crc_by_key=crc,
        object_bodies=bodies,
    )

    out = tmp_path / "out.tsv"
    summary = extract_cram(
        client,
        BUCKET,
        PREFIX,
        str(out),
        inline=True,
        show_progress=False,
    )

    assert summary.total == 1
    assert summary.crc_ok == 1
    assert summary.enrichment_ok == 1
    assert not summary.failures

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["filename"].endswith(".cram")
    assert "unmatched" not in row["filename"]
    assert row["s3_uri"] == f"s3://{BUCKET}/{MATCHED_KEY}"
    assert row["size_bytes"] == "2000000000000"
    assert row["crc64nvme_base64"] == "crc-matched"
    assert row["read_count"] == "9876543"
    assert row["crc_error"] == ""
    assert row["metadata_error"] == ""
    assert "read" not in row
    assert "lane" not in row
