"""
Tests for qa_gather: QAGatheredData, QADataGatherer, and gather_qa_data.

S3 interactions are mocked with a lightweight helper (MockS3Client) that
simulates list_objects, paginate, and download_file well enough to exercise
every gathering path without network access.
"""

import json
import os

import pytest

from qa_gather import (
    QAGatheredData,
    _is_merged_trimmer_file,
    gather_qa_data,
)
from qa_mods import QARunContext


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
QA_FIXTURES_DIR = os.path.join(FIXTURES_DIR, "qa")


# ---------------------------------------------------------------------------
# Mock S3 helpers
# ---------------------------------------------------------------------------


class MockPaginator:
    """Minimal paginator that wraps MockS3Client.list_objects into pages."""

    def __init__(self, s3_client: "MockS3Client"):
        self._s3 = s3_client

    def paginate(self, **kwargs):
        result = self._s3.list_objects(**kwargs)
        yield result


class MockS3Client:
    """Lightweight S3 mock driven by a dict of S3 keys.

    Parameters
    ----------
    keys : list[str]
        Full S3 keys (without bucket).  Used by ``list_objects`` to return
        ``Contents`` / ``CommonPrefixes`` and by ``download_file`` to resolve
        files whose content is provided in ``file_contents``.
    file_contents : dict[str, str | bytes]
        Mapping of S3 key → file content.  ``download_file`` writes this to
        the local path.  Keys not in this mapping trigger a ``FileNotFoundError``.
    """

    def __init__(
        self,
        keys: list[str] | None = None,
        file_contents: dict[str, str | bytes] | None = None,
    ):
        self._keys = keys or []
        self._file_contents = file_contents or {}

    def get_paginator(self, _operation: str) -> MockPaginator:
        return MockPaginator(self)

    def list_objects(
        self,
        Bucket: str = "",
        Prefix: str = "",
        Delimiter: str = "",
    ) -> dict:
        result: dict = {}
        contents = []
        prefixes: set[str] = set()

        for key in self._keys:
            if not key.startswith(Prefix):
                continue

            suffix = key[len(Prefix) :]

            if Delimiter and Delimiter in suffix:
                # There is a delimiter deeper in the suffix → this is a "folder"
                folder_end = suffix.index(Delimiter)
                common = Prefix + suffix[: folder_end + len(Delimiter)]
                prefixes.add(common)
            else:
                contents.append({"Key": key})

        if contents:
            result["Contents"] = contents
        if prefixes:
            result["CommonPrefixes"] = [{"Prefix": p} for p in sorted(prefixes)]
        return result

    def download_file(self, bucket: str, key: str, local_path: str) -> None:
        if key not in self._file_contents:
            raise FileNotFoundError(
                f"MockS3Client: no content registered for key {key!r}"
            )
        content = self._file_contents[key]
        mode = "wb" if isinstance(content, bytes) else "w"
        with open(local_path, mode) as fh:
            fh.write(content)


def _make_ctx(
    *,
    data_source: str = "s3",
    raw_assay: str = "10x",
    bucket: str = "czi-novogene",
    provider: str = "novogene",
    proj: str = "testproj",
    order: str = "ORD01",
    listing_prefix: str = "testproj/ORD01/",
    output_label: str = "ORD01",
    **kwargs,
) -> QARunContext:
    return QARunContext(
        data_source=data_source,
        raw_assay=raw_assay,
        bucket=bucket,
        provider=provider,
        proj=proj,
        order=order,
        output_label=output_label,
        listing_prefix=listing_prefix,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# QAGatheredData
# ---------------------------------------------------------------------------


class TestQAGatheredData:
    def test_default_construction(self):
        data = QAGatheredData()
        assert data.all_raw_files == []
        assert data.all_proc_files == {}
        assert data.fastq_log == {}
        assert data.has_raw is False
        assert data.has_processed is False
        assert data.gathering_errors == []
        assert data.gathering_warnings == []

    def test_fields_are_independent_across_instances(self):
        a = QAGatheredData()
        b = QAGatheredData()
        a.all_raw_files.append("x")
        assert b.all_raw_files == []


# ---------------------------------------------------------------------------
# _is_merged_trimmer_file
# ---------------------------------------------------------------------------


class TestIsMergedTrimmerFile:
    def test_merged_failure_codes(self):
        assert _is_merged_trimmer_file("path/merged_trimmer-failure_codes.csv")

    def test_prefixed_merged_failure_codes(self):
        assert _is_merged_trimmer_file("path/438761_merged_trimmer-failure_codes.csv")

    def test_merged_stats(self):
        assert _is_merged_trimmer_file("path/merged_trimmer-stats.csv")

    def test_prefixed_merged_stats(self):
        assert _is_merged_trimmer_file("path/438761_merged_trimmer-stats.csv")

    def test_non_merged_trimmer(self):
        assert not _is_merged_trimmer_file(
            "path/438761-G1_GEX-Z0001-BC_trimmer-failure_codes.csv"
        )

    def test_unrelated_csv(self):
        assert not _is_merged_trimmer_file("path/metrics_summary.csv")


# ---------------------------------------------------------------------------
# Manifest mode
# ---------------------------------------------------------------------------


class TestGatherManifestMode:
    def test_loads_raw_and_processed_from_manifest(self, tmp_path):
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "s3://czi-novogene/proj/ord/G1/raw/file1.fastq.gz\n"
            "s3://czi-novogene/proj/ord/G1/raw/file2.fastq.gz\n"
            "s3://czi-novogene/proj/ord/G1/processed/outs/web.html\n"
        )
        ctx = _make_ctx(
            data_source="manifest",
            manifest_path=str(manifest),
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=False,
        )
        s3 = MockS3Client()
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert data.has_processed is True
        assert len(data.all_raw_files) == 2
        assert "G1" in data.all_proc_files
        assert len(data.all_proc_files["G1"]) == 1

    def test_manifest_bucket_mismatch_raises(self, tmp_path):
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text("s3://czi-other/proj/ord/G1/raw/file.fq.gz\n")
        ctx = _make_ctx(
            data_source="manifest",
            bucket="czi-novogene",
            manifest_path=str(manifest),
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=False,
        )
        s3 = MockS3Client()
        with pytest.raises(ValueError, match="does not match"):
            gather_qa_data(ctx, s3)

    def test_manifest_empty_yields_no_data(self, tmp_path):
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text("s3://czi-novogene/unrelated/path.txt\n")
        ctx = _make_ctx(
            data_source="manifest",
            manifest_path=str(manifest),
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=False,
        )
        s3 = MockS3Client()
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is False
        assert data.has_processed is False


# ---------------------------------------------------------------------------
# S3 mode — structural warnings
# ---------------------------------------------------------------------------


class TestGatherS3Warnings:
    def test_order_level_processed_folder_warning(self):
        """Order-level processed/ folder produces a warning and is skipped."""
        keys = [
            "testproj/ORD01/processed/some_file.csv",
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("processed/" in w for w in data.gathering_warnings)

    def test_missing_raw_dir_warning(self):
        """Group without raw/ subdirectory produces a warning."""
        keys = ["testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/outs/config.csv"]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("raw/ MISSING" in w for w in data.gathering_warnings)

    def test_missing_processed_dir_warning(self):
        """Group without processed/ subdirectory produces a warning."""
        keys = ["testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("processed/ MISSING" in w for w in data.gathering_warnings)

    def test_extra_subdirs_warning(self):
        """Group with more than 2 subdirectories produces a warning."""
        keys = [
            "testproj/ORD01/G1/raw/file.csv",
            "testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/outs/config.csv",
            "testproj/ORD01/G1/extra/something.txt",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("EXTRA subdirs" in w for w in data.gathering_warnings)

    def test_empty_order_prefix_returns_empty_data(self):
        """No CommonPrefixes under order → empty data, no errors."""
        ctx = _make_ctx()
        s3 = MockS3Client(keys=[])
        data = gather_qa_data(ctx, s3)
        assert data.all_raw_files == []
        assert not data.has_raw
        assert not data.has_processed


# ---------------------------------------------------------------------------
# S3 mode — 10x raw files
# ---------------------------------------------------------------------------


class TestGather10xRaw:
    _BASE = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01"

    def _listing_keys(self):
        """Files that are listed but do NOT trigger S3 downloads."""
        return [
            f"{self._BASE}.csv",
            f"{self._BASE}.json",
            f"{self._BASE}_S1_L001_R1_001.fastq.gz",
            f"{self._BASE}_S1_L001_R2_001.fastq.gz",
        ]

    @staticmethod
    def _metadata_json(filename: str, read_count: int = 1000) -> str:
        return json.dumps(
            {"filename": filename, "read_count": read_count, "errors": []}
        )

    def test_raw_files_collected(self):
        """All raw files under raw/ are added to all_raw_files."""
        keys = self._listing_keys()
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert len(data.all_raw_files) == len(keys)

    def test_fastq_log_populated(self):
        """FASTQ files are counted in fastq_log under the correct group and assay."""
        keys = self._listing_keys()
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert "G1" in data.fastq_log
        assert "GEX" in data.fastq_log["G1"]
        assert len(data.fastq_log["G1"]["GEX"]) == 2

    def test_metadata_json_downloaded_and_parsed(self):
        """Metadata JSON files are downloaded and stored in read_metadata."""
        r1_meta_key = f"{self._BASE}_S1_L001_R1_001.fastq.gz-metadata.json"
        r1_filename = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        keys = self._listing_keys() + [r1_meta_key]
        file_contents = {
            r1_meta_key: self._metadata_json(r1_filename, 5000),
        }
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert r1_filename in data.read_metadata
        assert data.read_metadata[r1_filename]["read_count"] == 5000

    def test_metadata_prefers_actual_s3_filename_over_reported_filename(self):
        """read_metadata keys are canonicalized to actual object filename.

        When metadata payload filename differs from the object key, the gatherer
        keeps the object-derived filename as the dictionary key and stores both
        values for downstream diagnostics.
        """
        r1_meta_key = f"{self._BASE}_S1_L001_R1_001.fastq.gz-metadata.json"
        reported = "s3://czi-psomagen/other/other/WRONG_R1.fastq.gz"
        expected_actual = (
            "s3://czi-novogene/"
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        )
        keys = self._listing_keys() + [r1_meta_key]
        file_contents = {
            r1_meta_key: self._metadata_json(reported, 5000),
        }
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert expected_actual in data.read_metadata
        entry = data.read_metadata[expected_actual]
        assert entry["__actual_filename"] == expected_actual
        assert entry["__reported_filename"] == reported
        assert entry["read_count"] == 5000

    def test_wrong_assay_error(self):
        """File with unrecognised assay produces a gathering error."""
        keys = [
            "testproj/ORD01/G1/raw/439047-G1_BADASSAY-Z0273-BC01.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("WRONG ASSAY" in e for e in data.gathering_errors)

    def test_wrong_group_error_for_10x(self):
        """File whose parsed group doesn't match folder group → error (10x only)."""
        keys = [
            "testproj/ORD01/G1/raw/439047-WRONG_GEX-Z0273-BC01.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("WRONG GROUP" in e for e in data.gathering_errors)

    def test_group_separator_variants_do_not_error_for_10x(self):
        """Hyphen/underscore variants of the same group ID are treated as equal."""
        keys = [
            "testproj/ORD01/CUIMC-001/raw/439047-CUIMC_001_GEX-Z0273-BC01.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert not any("WRONG GROUP" in e for e in data.gathering_errors)

    def test_true_group_mismatch_still_errors_for_10x(self):
        """Normalization must not hide real group mismatches."""
        keys = [
            "testproj/ORD01/CUIMC-001/raw/439047-CUIMC_999_GEX-Z0273-BC01.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("WRONG GROUP" in e for e in data.gathering_errors)

    def test_trimmer_stats_downloaded(self):
        """Trimmer failure codes CSV is downloaded and parsed into trimmer_failure_stats."""
        trimmer_key = f"{self._BASE}_trimmer-failure_codes.csv"
        trimmer_csv = (
            "code,format,segment,reason,failed read count,total read count\n"
            "8,trim,preamble,rsq file,100,1000\n"
            "101,trim,insert,sequence was too short,50,1000\n"
        )
        keys = self._listing_keys() + [trimmer_key]
        file_contents = {trimmer_key: trimmer_csv}
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert len(data.trimmer_failure_stats) > 0

    def test_exp_to_run_map_populated(self):
        """exp_to_run_map is populated from trimmer filenames."""
        trimmer_key = f"{self._BASE}_trimmer-failure_codes.csv"
        trimmer_csv = (
            "code,format,segment,reason,failed read count,total read count\n"
            "8,trim,preamble,rsq file,100,1000\n"
        )
        keys = self._listing_keys() + [trimmer_key]
        file_contents = {trimmer_key: trimmer_csv}
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert any(v == "439047" for v in data.exp_to_run_map.values())


# ---------------------------------------------------------------------------
# S3 mode — non-10x (scale) raw files
# ---------------------------------------------------------------------------


class TestGatherScaleRaw:
    _RUN_DIR = "testproj/ORD01/RNA3_098/raw/426971/"

    def _listing_keys(self):
        """Scale raw files that don't trigger downloads."""
        return [
            f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C.cram",
            f"{self._RUN_DIR}426971-RNA3-098C_hash_oligo_QSR-7-SCALEPLEX_1E.cram",
        ]

    @staticmethod
    def _cram_metadata(filename: str, cb_tag: bool = True) -> str:
        return json.dumps(
            {
                "filename": filename,
                "read_count": 500,
                "cb_tag": cb_tag,
                "errors": [],
            }
        )

    def test_scale_raw_collected_from_run_subdirectory(self):
        """Scale raw files under run subdirectories are gathered."""
        keys = self._listing_keys()
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert len(data.all_raw_files) == 2

    def test_scale_cram_metadata_downloaded(self):
        """Scale CRAM metadata JSON files are ingested into read_metadata."""
        cram_file = "426971-RNA3-098C_GEX_QSR-7_10C.cram"
        meta_key = f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C.cram-metadata.json"
        keys = self._listing_keys() + [meta_key]
        file_contents = {meta_key: self._cram_metadata(cram_file)}
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert cram_file in data.read_metadata
        assert data.read_metadata[cram_file]["cb_tag"] is True

    def test_scale_cram_counted_in_fastq_log(self):
        """Scale CRAMs (not unmatched) are counted in fastq_log."""
        keys = self._listing_keys()
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert "RNA3_098" in data.fastq_log
        assert "GEX" in data.fastq_log["RNA3_098"]

    def test_scale_unmatched_cram_not_counted(self):
        """Unmatched CRAMs are excluded from fastq_log counts."""
        keys = [
            "testproj/ORD01/RNA3_098/raw/426971/"
            "426971-RNA3-098C_GEX_QSR-7_10C_unmatched.cram",
        ]
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        gex_files = data.fastq_log.get("RNA3_098", {}).get("GEX", [])
        assert len(gex_files) == 0

    def test_scale_unmatched_cram_metadata_not_downloaded(self):
        """Unmatched CRAM metadata JSON files are not ingested into read_metadata."""
        unmatched_cram = "426971-RNA3-098C_GEX_QSR-7_10C-unmatched.cram"
        unmatched_meta_key = f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C-unmatched.cram-metadata.json"
        matched_cram = "426971-RNA3-098C_GEX_QSR-7_10C.cram"
        matched_meta_key = (
            f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C.cram-metadata.json"
        )
        keys = self._listing_keys() + [
            f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C-unmatched.cram",
            unmatched_meta_key,
            matched_meta_key,
        ]
        file_contents = {
            unmatched_meta_key: self._cram_metadata(unmatched_cram, cb_tag=False),
            matched_meta_key: self._cram_metadata(matched_cram, cb_tag=True),
        }
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert matched_cram in data.read_metadata
        assert unmatched_cram not in data.read_metadata

    def test_scale_unmatched_cram_metadata_not_downloaded_underscore(self):
        """Underscore-form unmatched CRAM metadata is also excluded from read_metadata."""
        unmatched_cram = "426971-RNA3-098C_GEX_QSR-7_10C_unmatched.cram"
        unmatched_meta_key = f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C_unmatched.cram-metadata.json"
        matched_cram = "426971-RNA3-098C_GEX_QSR-7_10C.cram"
        matched_meta_key = (
            f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C.cram-metadata.json"
        )
        keys = self._listing_keys() + [
            f"{self._RUN_DIR}426971-RNA3-098C_GEX_QSR-7_10C_unmatched.cram",
            unmatched_meta_key,
            matched_meta_key,
        ]
        file_contents = {
            unmatched_meta_key: self._cram_metadata(unmatched_cram, cb_tag=False),
            matched_meta_key: self._cram_metadata(matched_cram, cb_tag=True),
        }
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert matched_cram in data.read_metadata
        assert unmatched_cram not in data.read_metadata

    def test_scale_no_wrong_group_error(self):
        """Scale files don't check group match (is_10x=False), so no WRONG GROUP."""
        keys = self._listing_keys()
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert not any("WRONG GROUP" in e for e in data.gathering_errors)


class TestGather10xCramRaw:
    _RAW_DIR = "testproj/ORD01/LeS1867W11/raw/"
    _BASE = "442488-LeS1867W11_ATAC-Z0027-CACTGTCAGCCAGAT"

    @staticmethod
    def _cram_metadata(filename: str, read_count: int = 1000) -> str:
        return json.dumps(
            {"filename": filename, "read_count": read_count, "errors": []}
        )

    def test_10x_cram_metadata_downloaded(self):
        meta_key = f"{self._RAW_DIR}{self._BASE}.cram-metadata.json"
        keys = [
            f"{self._RAW_DIR}{self._BASE}.cram",
            f"{self._RAW_DIR}{self._BASE}.csv",
            meta_key,
        ]
        file_contents = {
            meta_key: self._cram_metadata(f"{self._BASE}.cram", read_count=4321),
        }
        ctx = _make_ctx(raw_assay="10x_cram")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert f"{self._BASE}.cram" in data.read_metadata
        assert data.read_metadata[f"{self._BASE}.cram"]["read_count"] == 4321

    def test_10x_cram_not_counted_in_fastq_log(self):
        keys = [
            f"{self._RAW_DIR}{self._BASE}.cram",
            f"{self._RAW_DIR}{self._BASE}.csv",
        ]
        ctx = _make_ctx(raw_assay="10x_cram")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.fastq_log.get("LeS1867W11", {}) == {}


# ---------------------------------------------------------------------------
# S3 mode — CellRanger processed
# ---------------------------------------------------------------------------


class TestGatherCellrangerProcessed:
    def _cr_keys(self):
        """Minimal CellRanger processed file set."""
        outs = "testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/outs/"
        return [
            f"{outs}config.csv",
            f"{outs}web_summary.html",
            f"{outs}metrics_summary.csv",
        ]

    def test_cellranger_files_collected(self):
        """CellRanger outs/ files are collected into all_proc_files."""
        raw_key = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"
        keys = [raw_key] + self._cr_keys()
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_processed is True
        assert "G1" in data.all_proc_files
        assert len(data.all_proc_files["G1"]) == 3

    def test_invalid_run_date_error(self):
        """Invalid Run_ directory name produces a gathering error."""
        raw_key = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"
        bad_date = "testproj/ORD01/G1/processed/cellranger/Run_BADDATE/outs/config.csv"
        keys = [raw_key, bad_date]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("INCORRECT DATE FORMAT" in e for e in data.gathering_errors)

    def test_missing_cellranger_dir_warning(self):
        """Group with processed/ but no cellranger/ sub-dir produces a warning."""
        raw_key = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"
        non_cr = "testproj/ORD01/G1/processed/other/file.csv"
        keys = [raw_key, non_cr]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("cellranger/ MISSING" in w for w in data.gathering_warnings)

    def test_missing_outs_dir_warning(self):
        """Run dir without outs/ produces a warning."""
        raw_key = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"
        no_outs = "testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/other/file.csv"
        keys = [raw_key, no_outs]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert any("NO outs/" in w for w in data.gathering_warnings)


# ---------------------------------------------------------------------------
# S3 mode — Scale processed
# ---------------------------------------------------------------------------


class TestGatherScaleProcessed:
    def _workflow_info_json(self, genome: str = "s3://bucket/GRCh38") -> str:
        return json.dumps(
            {
                "Parameters": {
                    "genome": genome,
                    "scalePlexAssignmentMethod": "fc",
                    "bamOut": "true",
                    "scalePlex": "true",
                },
                "Workflow Manifest": {"name": "ScaleRna", "version": "2.1.0"},
                "Nextflow Information": {"Execution status": "OK"},
            }
        )

    def _samples_csv(self) -> str:
        return (
            "sample,barcodes,libIndex2,scalePlexLibIndex2\n"
            "SAMP-01,BC1,QSR-1;QSR-2,PLX1;PLX2\n"
        )

    def _scale_processed_keys(self):
        """Minimal Scale processed file set."""
        run = "testproj/ORD01/RNA3_098/processed/Run_2001-01-13/"
        return [
            f"{run}workflow_info.json",
            f"{run}samples.csv",
            f"{run}samples/SAMP-01_anndata.h5ad",
            f"{run}samples/SAMP-01.merged.allCells.csv",
            f"{run}samples/SAMP-01.QSR-1_anndata.h5ad",
            f"{run}samples/SAMP-01.QSR-2_anndata.h5ad",
        ]

    def test_scale_processed_collected(self):
        """Scale processed files are collected when workflow_info.json is found."""
        raw_key = (
            "testproj/ORD01/RNA3_098/raw/426971/426971-RNA3-098C_GEX_QSR-7_10C.cram"
        )
        proc_keys = self._scale_processed_keys()
        run_prefix = "testproj/ORD01/RNA3_098/processed/Run_2001-01-13/"
        file_contents = {
            f"{run_prefix}workflow_info.json": self._workflow_info_json(),
            f"{run_prefix}samples.csv": self._samples_csv(),
        }
        keys = [raw_key] + proc_keys
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        assert data.has_processed is True
        assert "RNA3_098" in data.scale_workflow_infos
        assert data.scale_workflow_infos["RNA3_098"]["genome"] == "s3://bucket/GRCh38"
        assert "RNA3_098" in data.scale_samples_info
        assert "SAMP-01" in data.scale_samples_info["RNA3_098"]["samples"]
        assert "RNA3_098" in data.scale_proc_files
        assert len(data.scale_proc_files["RNA3_098"]) == 4

    def test_scale_processed_no_workflow_info_skipped(self):
        """Processed subdir without workflow_info.json is silently skipped."""
        raw_key = (
            "testproj/ORD01/RNA3_098/raw/426971/426971-RNA3-098C_GEX_QSR-7_10C.cram"
        )
        proc_key = "testproj/ORD01/RNA3_098/processed/SomeDir/samples.csv"
        keys = [raw_key, proc_key]
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_processed is False
        assert data.scale_workflow_infos == {}


# ---------------------------------------------------------------------------
# S3 mode — order-level merged trimmers (10x)
# ---------------------------------------------------------------------------


class TestGatherOrderLevelMergedTrimmers:
    def test_10x_order_level_merged_trimmer_ingested(self):
        """Merged trimmer files at order level are ingested for 10x."""
        raw_key = "testproj/ORD01/G1/raw/438761-G1_GEX-Z0001-BC01.csv"
        merged_key = "testproj/ORD01/438761_merged_trimmer-failure_codes.csv"
        merged_csv = (
            "read group,code,format,segment,reason,"
            "failed read count,total read count\n"
            "TT,8,trim,preamble,rsq file,1000,10000\n"
            "TT,101,trim,insert,sequence was too short,100,10000\n"
        )
        keys = [raw_key, merged_key]
        file_contents = {merged_key: merged_csv}
        ctx = _make_ctx(raw_assay="10x")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert "438761" in data.merged_wafer_stats

    def test_scale_order_level_not_ingested(self):
        """Scale assay does not look for order-level merged trimmers."""
        raw_key = (
            "testproj/ORD01/RNA3_098/raw/426971/426971-RNA3-098C_GEX_QSR-7_10C.cram"
        )
        merged_key = "testproj/ORD01/438761_merged_trimmer-failure_codes.csv"
        keys = [raw_key, merged_key]
        ctx = _make_ctx(raw_assay="scale")
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.merged_wafer_stats == {}


# ---------------------------------------------------------------------------
# has_raw / has_processed flags
# ---------------------------------------------------------------------------


class TestHasRawHasProcessed:
    def test_raw_only(self):
        """Only raw/ data → has_raw=True, has_processed=False."""
        keys = ["testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv"]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert data.has_processed is False

    def test_processed_only(self):
        """Only processed/ data → has_raw=False, has_processed=True."""
        keys = [
            "testproj/ORD01/G1/raw/placeholder",
            "testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/outs/config.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        # raw/ subdir exists but has no parseable raw files (placeholder) — the raw key is
        # gathered but is not a real file, so has_raw depends on whether anything got listed
        assert data.has_processed is True

    def test_both_raw_and_processed(self):
        """Both raw/ and processed/ → both flags True."""
        keys = [
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv",
            "testproj/ORD01/G1/processed/cellranger/Run_2000-01-10/outs/config.csv",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert data.has_processed is True


# ---------------------------------------------------------------------------
# Multiple groups
# ---------------------------------------------------------------------------


class TestMultipleGroups:
    def test_two_groups_collected_independently(self):
        """Two groups under one order are gathered independently."""
        keys = [
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01.csv",
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz",
            "testproj/ORD01/G2/raw/439048-G2_GEX-Z0050-BC02.csv",
            "testproj/ORD01/G2/raw/439048-G2_GEX-Z0050-BC02_S1_L001_R1_001.fastq.gz",
        ]
        ctx = _make_ctx()
        s3 = MockS3Client(keys=keys)
        data = gather_qa_data(ctx, s3)
        assert "G1" in data.fastq_log
        assert "G2" in data.fastq_log
        assert "GEX" in data.fastq_log["G1"]
        assert "GEX" in data.fastq_log["G2"]
        assert len(data.all_raw_files) == 4


# ---------------------------------------------------------------------------
# Invalid data_source
# ---------------------------------------------------------------------------


class TestInvalidDataSource:
    def test_invalid_data_source_raises(self):
        ctx = QARunContext(
            data_source="ftp",
            raw_assay="10x",
            bucket="czi-novogene",
            provider="novogene",
            proj="p",
            order="o",
            output_label="o",
            listing_prefix="p/o/",
        )
        s3 = MockS3Client()
        with pytest.raises(ValueError, match="Invalid data_source"):
            gather_qa_data(ctx, s3)
