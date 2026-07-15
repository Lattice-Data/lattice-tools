"""
Tests for SeaHub lab raw QA support (seahub_sci assay).
"""

from __future__ import annotations

import os

from qa_checks import (
    check_expected_raw_files,
    check_extra_raw_files,
    validate_fastq_counts,
    validate_read_metadata,
)
from qa_gather import gather_qa_data
from qa_mods import (
    QARunContext,
    grab_seahub_trim_fail_csv,
    normalize_raw_assay,
    parse_raw_filename,
    parse_seahub_raw_path,
    resolve_qa_run_context,
    seahub_file_stem,
    seahub_trimmer_failure_storage_key,
    seahub_trimmer_group_storage_key,
)

from tests.test_qa_gather import MockS3Client, _make_ctx

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
QA_FIXTURES_DIR = os.path.join(FIXTURES_DIR, "qa")

SEAHUB_BASE = (
    "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
    "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
)
SEAHUB_KEY_CRAM = f"{SEAHUB_BASE}.trim.cram"
SEAHUB_KEY_FAIL = f"{SEAHUB_BASE}.trim_fail.csv"
SEAHUB_KEY_TRIM = f"{SEAHUB_BASE}.trim.csv"
SEAHUB_TRIM_FAIL = os.path.join(QA_FIXTURES_DIR, "seahub_trim_fail_sample.csv")
SEAHUB_TRIM = os.path.join(QA_FIXTURES_DIR, "seahub_trim_sample.csv")


class TestSeahubConstantsAndContext:
    def test_normalize_seahub_sci(self):
        assert normalize_raw_assay("seahub_sci") == "seahub_sci"

    def test_resolve_context_labalpha(self):
        ctx = resolve_qa_run_context(
            data_source="s3",
            raw_assay="seahub_sci",
            s3_path="s3://czi-labalpha/labalpha-seahub-bcp/REF3",
        )
        assert ctx.bucket == "czi-labalpha"
        assert ctx.provider == "labalpha"
        assert ctx.proj == "labalpha-seahub-bcp"
        assert ctx.order == "REF3"
        assert ctx.listing_prefix == "labalpha-seahub-bcp/REF3/"


class TestSeahubParsing:
    def test_parse_seahub_raw_path(self):
        info = parse_seahub_raw_path(SEAHUB_KEY_CRAM)
        assert info == {
            "experiment_id": "REF3",
            "sublibrary": "P05_1",
            "wafer": "430479",
        }

    def test_seahub_file_stem(self):
        assert (
            seahub_file_stem(
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim.cram"
            )
            == "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        )
        assert (
            seahub_file_stem(
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv"
            )
            == "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        )

    def test_parse_raw_filename_seahub_cram(self):
        assert parse_raw_filename(SEAHUB_KEY_CRAM, "seahub_sci") == (
            "430479",
            "REF3",
            "GEX_hash_oligo",
            "Z0097",
            "CAGTCAGTTGCAGAT",
        )

    def test_parse_raw_filename_seahub_trim_fail(self):
        assert parse_raw_filename(SEAHUB_KEY_FAIL, "seahub_sci") == (
            "430479",
            "REF3",
            "GEX_hash_oligo",
            "Z0097",
            "CAGTCAGTTGCAGAT",
        )

    def test_storage_keys(self):
        assert seahub_trimmer_group_storage_key(SEAHUB_KEY_FAIL) == "REF3/P05_1"
        assert seahub_trimmer_failure_storage_key(SEAHUB_KEY_FAIL) == (
            "430479",
            "430479",
        )


class TestSeahubTrimAdapters:
    def test_grab_seahub_trim_fail_csv_is_per_format(self):
        """One trimmer-fail value per format block, each < 100% (regression:
        the shared single-denominator parser summed both modalities' failures
        over one modality's total and produced >100%)."""
        stats: dict = {}
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        assert "430479" in stats
        # SeaHub inputs start from RSQ-passing reads -> no RSQ rows.
        assert stats["430479"]["rsq"] == []
        # Two format blocks (JumboSciHash, JumboSciGEX) -> two values.
        fails = stats["430479"]["trimmer_fail"]
        assert len(fails) == 2
        assert all(0 <= pct < 100 for pct in fails)

    def test_grab_seahub_trim_fail_csv_denominator_per_group(self):
        """Each modality's failures divide by its own total read count."""
        stats: dict = {}
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        hash_pct, gex_pct = stats["430479"]["trimmer_fail"]
        assert hash_pct == 100 * 157722381 / 260527531
        assert gex_pct == 100 * 100735571 / 158233602

    def test_grab_seahub_trim_fail_csv_appends(self):
        """Repeated wells accumulate into the same distribution."""
        stats: dict = {}
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        assert len(stats["430479"]["trimmer_fail"]) == 4


class TestSeahubChecks:
    def _full_stem_files(self, stem: str, raw_dir: str) -> list[str]:
        suffixes = [
            ".trim.cram",
            ".trim.csv",
            ".trim.stderr",
            ".trim.stdout",
            ".trim_fail.csv",
        ]
        return [f"{raw_dir}/{stem}{suf}" for suf in suffixes]

    def test_expected_raw_files_complete(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir)
        good, lost, found = check_expected_raw_files(all_raw, "seahub_sci")
        assert good == 1
        assert lost == []
        assert len(found) == 5

    def test_expected_raw_files_missing_suffix(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir)[:-1]
        _good, lost, _found = check_expected_raw_files(all_raw, "seahub_sci")
        assert len(lost) == 1
        assert ".trim_fail.csv" in lost[0]

    def test_extra_raw_allows_metadata(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir) + [
            f"{raw_dir}/{stem}.trim.cram-metadata.json"
        ]
        _good, _lost, found = check_expected_raw_files(all_raw, "seahub_sci")
        extra = check_extra_raw_files(all_raw, found, "seahub_sci")
        assert extra == []

    def test_fastq_counts_gex_hash_only(self):
        fastq_log = {
            "REF3": {
                "GEX_hash_oligo": ["a.trim.cram", "b.trim.cram"],
            }
        }
        assert validate_fastq_counts(fastq_log, "seahub_sci") == []

    def test_read_metadata_skips_r1_r2_pairing(self):
        meta = {
            SEAHUB_KEY_CRAM: {
                "read_count": 2_000_000,
                "__reported_filename": SEAHUB_KEY_CRAM.split("/")[-1],
            }
        }
        counts, errors, pairing = validate_read_metadata(meta, "seahub_sci")
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []


class TestSeahubGather:
    def _seahub_paginated_pages(self) -> dict:
        o = "labalpha-seahub-bcp/REF3/"
        raw = f"{o}raw/"
        sublib = f"{raw}P05_1/"
        wafer = f"{sublib}430479/"
        return {
            (o, "/"): [
                {"CommonPrefixes": [{"Prefix": raw}]},
            ],
            (raw, "/"): [
                {"CommonPrefixes": [{"Prefix": sublib}]},
            ],
            (sublib, "/"): [
                {"CommonPrefixes": [{"Prefix": wafer}]},
            ],
            (wafer, ""): [
                {
                    "Contents": [
                        {"Key": SEAHUB_KEY_CRAM},
                        {"Key": SEAHUB_KEY_FAIL},
                        {"Key": SEAHUB_KEY_TRIM},
                        {"Key": f"{SEAHUB_BASE}.trim.stderr"},
                        {"Key": f"{SEAHUB_BASE}.trim.stdout"},
                    ]
                }
            ],
        }

    def test_gather_from_s3_seahub_layout(self):
        pages = self._seahub_paginated_pages()
        paginated = dict(pages)
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            SEAHUB_KEY_TRIM: open(SEAHUB_TRIM).read(),
        }
        s3 = MockS3Client(paginated_pages=paginated, file_contents=file_contents)
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert len(data.all_raw_files) == 5
        assert "REF3" in data.fastq_log
        assert "GEX_hash_oligo" in data.fastq_log["REF3"]
        assert "430479" in data.trimmer_failure_stats
        assert "REF3/P05_1" in data.group_failure_stats

    def test_gather_manifest_enriches_seahub(self):
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        ctx = QARunContext(
            data_source="manifest",
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="",
            order="REF3",
            output_label="REF3",
            listing_prefix="",
            manifest_path=manifest,
            manifest_delimiter="\t",
            manifest_s3_column=6,
            manifest_has_header=True,
        )
        trim_fail_key = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv"
        )
        trim_key = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim.csv"
        )
        trim_fail_a2 = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim_fail.csv"
        )
        trim_a2 = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim.csv"
        )
        file_contents = {
            trim_fail_key: open(SEAHUB_TRIM_FAIL).read(),
            trim_key: open(SEAHUB_TRIM).read(),
            trim_fail_a2: open(SEAHUB_TRIM_FAIL).read(),
            trim_a2: open(SEAHUB_TRIM).read(),
        }
        s3 = MockS3Client(file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert len(data.all_raw_files) == 10
        assert "430479" in data.trimmer_failure_stats

    def _manifest_ctx(self, manifest: str) -> QARunContext:
        return QARunContext(
            data_source="manifest",
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="",
            order="REF3",
            output_label="REF3",
            listing_prefix="",
            manifest_path=manifest,
            manifest_delimiter="\t",
            manifest_s3_column=6,
            manifest_has_header=True,
        )

    def test_missing_metadata_json_warns(self):
        """Each .trim.cram lacking a .trim.cram-metadata.json sidecar warns."""
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        ctx = self._manifest_ctx(manifest)
        file_contents = {
            k: open(SEAHUB_TRIM_FAIL).read()
            for k in (
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv",
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim_fail.csv",
            )
        }
        s3 = MockS3Client(file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        missing = [
            w for w in data.gathering_warnings if w.startswith("METADATA MISSING")
        ]
        assert len(missing) == 2
        assert all(".trim.cram-metadata.json" in w for w in missing)

    def test_present_metadata_json_no_warn(self):
        """A .trim.cram with its .trim.cram-metadata.json sidecar does not warn."""
        pages = dict(self._seahub_paginated_pages())
        wafer_key = ("labalpha-seahub-bcp/REF3/raw/P05_1/430479/", "")
        meta_key = f"{SEAHUB_BASE}.trim.cram-metadata.json"
        pages[wafer_key] = [
            {"Contents": pages[wafer_key][0]["Contents"] + [{"Key": meta_key}]}
        ]
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            meta_key: '{"read_count": 1000, "filename": "%s"}'
            % SEAHUB_KEY_CRAM.split("/")[-1],
        }
        s3 = MockS3Client(paginated_pages=pages, file_contents=file_contents)
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, s3)
        missing = [
            w for w in data.gathering_warnings if w.startswith("METADATA MISSING")
        ]
        assert missing == []
