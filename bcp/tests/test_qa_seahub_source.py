"""
Tests for cross-bucket SeaHub trimming reconciliation (qa_seahub_source).

The vendor layout mirrored here is the real one measured from
``czi-novogene/labalpha-seahub-bcp/NVUS0000000000-02``: the group folder sits
*before* ``raw``, giving a different path depth from the lab's trimmed upload.
"""

from __future__ import annotations

import json

import pytest

from qa_seahub_recon import reconcile_trimming
from qa_seahub_source import (
    index_trimmed_upload,
    index_untrimmed_source,
    load_source_read_counts,
    parse_source_uri,
)

from tests.test_qa_gather import MockS3Client

SOURCE_URI = "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-02"
SOURCE_DIR = "labalpha-seahub-bcp/NVUS0000000000-02/REF5_P01/raw/448642"
TRIMMED_DIR = "labalpha-seahub-bcp/REF5/raw/REF5_P01/448642"

VENDOR_STEM_A = "448642-REF5_P01_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
VENDOR_STEM_B = "448642-REF5_P01_GEX_hash_oligo-Z0002-AGCTCGAATGCGATC"
TRIM_STEM_A = "448642-REF5_P01_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
TRIM_STEM_B = "448642-REF5_P01_A2_GEX_hash_oligo-Z0002-AGCTCGAATGCGATC"


def _vendor_keys(*stems: str) -> list[str]:
    keys: list[str] = []
    for stem in stems:
        keys.extend(
            [
                f"{SOURCE_DIR}/{stem}.cram",
                f"{SOURCE_DIR}/{stem}.cram-metadata.json",
                f"{SOURCE_DIR}/{stem}.csv",
                f"{SOURCE_DIR}/{stem}.json",
                f"{SOURCE_DIR}/{stem}_FlowQ.metric",
                f"{SOURCE_DIR}/{stem}_trimmer-stats.csv",
            ]
        )
    keys.extend(
        [
            f"{SOURCE_DIR}/448642_LibraryInfo.xml",
            f"{SOURCE_DIR}/448642_UploadCompleted.json",
            f"{SOURCE_DIR}/448642_merged_trimmer-stats.csv",
            f"{SOURCE_DIR}/448642_run_SecondaryAnalysis.txt",
            f"{SOURCE_DIR}/448642-unmatched.cram",
            "labalpha-seahub-bcp/NVUS0000000000-02/file-manifest.json",
        ]
    )
    return keys


def _trimmed_keys(*stems: str, omit_cram_for: str = "") -> list[str]:
    keys: list[str] = []
    for stem in stems:
        suffixes = [".trim.csv", ".trim.stderr", ".trim.stdout", ".trim_fail.csv"]
        if stem != omit_cram_for:
            suffixes += [".trim.cram", ".trim.cram-metadata.json"]
        keys.extend(f"{TRIMMED_DIR}/{stem}{suffix}" for suffix in suffixes)
    return keys


def _source_index(keys: list[str]) -> dict:
    return index_untrimmed_source(MockS3Client(keys=keys), SOURCE_URI)


class TestParseSourceUri:
    def test_splits_bucket_and_prefix(self):
        assert parse_source_uri(SOURCE_URI) == (
            "czi-novogene",
            "labalpha-seahub-bcp/NVUS0000000000-02/",
        )

    def test_rejects_non_s3_uri(self):
        with pytest.raises(ValueError):
            parse_source_uri("/local/path")


class TestIndexUntrimmedSource:
    def test_indexes_one_entry_per_well(self):
        index = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        assert sorted(index) == [("448642", "Z0001"), ("448642", "Z0002")]

    def test_captures_identity_fields_and_sidecar(self):
        entry = _source_index(_vendor_keys(VENDOR_STEM_A))[("448642", "Z0001")]
        assert entry.barcode == "CAGCTCGAATGCGAT"
        assert entry.group == "REF5_P01"
        assert entry.assay == "GEX_hash_oligo"
        assert entry.metadata_key.endswith(".cram-metadata.json")

    def test_skips_wafer_level_and_unmatched_files(self):
        index = _source_index(_vendor_keys(VENDOR_STEM_A))
        assert len(index) == 1
        assert "unmatched" not in index[("448642", "Z0001")].cram_key

    def test_reads_counts_from_sidecars(self):
        keys = _vendor_keys(VENDOR_STEM_A)
        contents = {
            f"{SOURCE_DIR}/{VENDOR_STEM_A}.cram-metadata.json": json.dumps(
                {"read_count": 260527531}
            )
        }
        s3 = MockS3Client(keys=keys, file_contents=contents)
        index = index_untrimmed_source(s3, SOURCE_URI)
        load_source_read_counts(s3, "czi-novogene", index)
        assert index[("448642", "Z0001")].read_count == 260527531


class TestIndexTrimmedUpload:
    def test_indexes_by_wafer_and_ug(self):
        index = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A, TRIM_STEM_B))
        assert sorted(index) == [("448642", "Z0001"), ("448642", "Z0002")]
        assert index[("448642", "Z0001")].sublibrary == "REF5_P01"
        assert index[("448642", "Z0001")].has_cram is True

    def test_flags_a_well_whose_cram_is_absent(self):
        index = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A, omit_cram_for=TRIM_STEM_A)
        )
        assert index[("448642", "Z0001")].has_cram is False

    def test_indexes_the_bare_family_too(self):
        bare = "448642-448642-REF5_P01_A1-Z0001-CAGCTCGAATGCGAT"
        keys = [f"{TRIMMED_DIR}/{bare}{s}" for s in (".cram", ".csv", "_fail.csv")]
        index = index_trimmed_upload(keys)
        entry = index[("448642", "Z0001")]
        assert entry.family == "bare"
        assert entry.assay is None
        assert entry.has_cram is True


class TestReconcileTrimming:
    def _fail_counts(self, stem: str, total: int) -> dict:
        return {stem: {"JumboSciGEX": {"failed": 100, "total": total}}}

    def test_consistent_pair_is_clean(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(TRIM_STEM_A, 260527531)
        )
        assert report.rows == []
        assert report.matched == 1

    def test_vendor_cram_with_no_trim_output(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed)
        not_trimmed = [r for r in report.rows if r["category"] == "not_trimmed"]
        assert [r["ug"] for r in not_trimmed] == ["Z0002"]

    def test_trimmed_well_whose_cram_is_absent_counts_as_not_trimmed(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A, omit_cram_for=TRIM_STEM_A)
        )
        report = reconcile_trimming(source, trimmed)
        assert [r["category"] for r in report.rows] == ["not_trimmed"]

    def test_orphan_trimmed_well(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A, TRIM_STEM_B))
        report = reconcile_trimming(source, trimmed)
        orphans = [r for r in report.rows if r["category"] == "orphan_trimmed"]
        assert [r["ug"] for r in orphans] == ["Z0002"]

    def test_identity_mismatch_on_barcode(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 100
        wrong_barcode = "448642-REF5_P01_A1_GEX_hash_oligo-Z0001-TTTTTTTTTTTTTTT"
        trimmed = index_trimmed_upload(_trimmed_keys(wrong_barcode))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(wrong_barcode, 100)
        )
        assert [r["category"] for r in report.rows] == ["identity_mismatch"]
        assert "barcode differs" in report.rows[0]["detail"]

    def test_read_count_mismatch_when_trimmer_saw_fewer_reads(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(TRIM_STEM_A, 158233602)
        )
        assert [r["category"] for r in report.rows] == ["read_count_mismatch"]
        assert "260527531" in report.rows[0]["detail"]

    def test_any_matching_format_total_is_accepted(self):
        # Format blocks declare different totals for the same well, so matching
        # any one of them means the trimmer consumed the delivered file.
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        fail_counts = {
            TRIM_STEM_A: {
                "JumboSciGEX": {"failed": 1, "total": 158233602},
                "JumboSciHash": {"failed": 2, "total": 260527531},
            }
        }
        report = reconcile_trimming(source, trimmed, fail_counts)
        assert report.rows == []

    def test_missing_vendor_read_count_is_reported_as_unavailable(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed, self._fail_counts(TRIM_STEM_A, 1))
        assert [r["category"] for r in report.rows] == ["metadata_unavailable"]

    def test_verdict_summarizes_both_sides(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed)
        assert "2 vendor CRAMs" in report.verdict()
        assert "1 trimmed wells" in report.verdict()
        assert report.counts["not_trimmed"] == 1
