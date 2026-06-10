"""Tests for per-sample trimmer-stats parsing and wafer aggregation helpers."""

from __future__ import annotations

import os

import pytest

from qa_mods import (
    finalize_merged_wafer_stats,
    grab_sample_trimmer_stats_metrics,
    grab_trimmer_failure_codes_wafer_metrics,
    merge_partial_wafer_stats,
)

QA_FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "qa")
FLEX_STATS = os.path.join(QA_FIXTURES_DIR, "psomagen_trimmer_stats_flex_v2_sample.csv")


class TestGrabSampleTrimmerStatsMetrics:
    def test_psomagen_flex_v2_fixture_tt_and_read2s_q30(self):
        result = grab_sample_trimmer_stats_metrics(FLEX_STATS)
        assert result is not None
        assert result["tt_total_reads"] == 10695752970
        assert result["tt_failed_reads"] == 3559668377
        assert result["tt_pass_reads"] == 7136084593
        assert result["_q30_matched_bases"] == 148512897203
        assert result["_q30_failures"] == 73120475

    def test_missing_columns_returns_none(self, tmp_path):
        path = tmp_path / "bad.csv"
        path.write_text("read group,segment label\nTT,low quality bases\n")
        assert grab_sample_trimmer_stats_metrics(str(path)) is None


class TestGrabTrimmerFailureCodesWaferMetrics:
    def test_per_sample_rsq_totals(self, tmp_path):
        path = tmp_path / "failure.csv"
        path.write_text(
            "code,format,segment,reason,failed read count,total read count\n"
            "8,trim,preamble,rsq file,100,1000\n"
            "101,trim,insert,sequence was too short,50,1000\n"
        )
        result = grab_trimmer_failure_codes_wafer_metrics(str(path))
        assert result is not None
        assert result["rsq_total_reads"] == 1000
        assert result["rsq_failed_reads"] == 100
        assert result["rsq_pass_reads"] == 900


class TestMergePartialWaferStats:
    def test_sums_across_libraries_for_same_run_id(self):
        merged: dict = {}
        merge_partial_wafer_stats(
            merged,
            "441049",
            {
                "tt_total_reads": 1000,
                "tt_failed_reads": 100,
                "rsq_total_reads": 500,
                "rsq_failed_reads": 50,
                "_q30_matched_bases": 900,
                "_q30_failures": 100,
            },
        )
        merge_partial_wafer_stats(
            merged,
            "441049",
            {
                "tt_total_reads": 2000,
                "tt_failed_reads": 200,
                "rsq_total_reads": 800,
                "rsq_failed_reads": 80,
                "_q30_matched_bases": 1800,
                "_q30_failures": 200,
            },
        )
        finalize_merged_wafer_stats(merged)
        stats = merged["441049"]
        assert stats["tt_total_reads"] == 3000
        assert stats["tt_failed_reads"] == 300
        assert stats["rsq_total_reads"] == 1300
        assert stats["rsq_failed_reads"] == 130
        assert stats["tt_q30_pct"] == pytest.approx(
            100.0 * (900 + 1800) / (900 + 100 + 1800 + 200), rel=1e-9
        )

    def test_skips_tt_rsq_when_merged_failure_codes_present(self):
        merged = {
            "439047": {
                "_from_merged_failure_codes": True,
                "tt_total_reads": 10000,
                "tt_failed_reads": 1100,
                "rsq_total_reads": 1000,
                "rsq_failed_reads": 100,
            }
        }
        merge_partial_wafer_stats(
            merged,
            "439047",
            {"tt_total_reads": 999, "tt_failed_reads": 99, "rsq_total_reads": 999},
        )
        assert merged["439047"]["tt_total_reads"] == 10000
        assert merged["439047"]["rsq_total_reads"] == 1000
