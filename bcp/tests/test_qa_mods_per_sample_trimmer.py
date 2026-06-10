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
NOVO_SAMPLE_STATS = os.path.join(
    QA_FIXTURES_DIR, "novogene_per_sample_trimmer_stats.csv"
)
PSOM_SAMPLE_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "psomagen_per_sample_trimmer_failure_codes.csv"
)
NOVO_SAMPLE_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "novogene_per_sample_trimmer_failure_codes.csv"
)
NOVO_MERGED_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "novogene_merged_trimmer_failure_codes_sample.csv"
)


class TestGrabSampleTrimmerStatsMetrics:
    def test_flex_v2_yields_read2s_q30_and_no_tt(self):
        """Per-sample stats provide only Q30 components, never tt_* counts."""
        result = grab_sample_trimmer_stats_metrics(FLEX_STATS)
        assert result is not None
        # No TT-derived counts: per-sample files have no template read group.
        assert "tt_total_reads" not in result
        assert "tt_failed_reads" not in result
        # Q30 comes from the Read 2S segment for flex V2.
        assert result["_sample_q30_matched_bases"] == 148512897203
        assert result["_sample_q30_failures"] == 73120475

    def test_multiome_yields_low_quality_bases_q30(self):
        """Standard 10x layouts use the low quality bases segment for Q30."""
        result = grab_sample_trimmer_stats_metrics(NOVO_SAMPLE_STATS)
        assert result is not None
        assert "tt_total_reads" not in result
        assert result["_sample_q30_matched_bases"] == 1285752325
        assert result["_sample_q30_failures"] == 1473

    def test_missing_columns_returns_none(self, tmp_path):
        path = tmp_path / "bad.csv"
        path.write_text("read group,segment label\nTT,low quality bases\n")
        assert grab_sample_trimmer_stats_metrics(str(path)) is None


class TestGrabTrimmerFailureCodesWaferMetrics:
    def test_per_sample_rsq_totals_no_read_group_column(self, tmp_path):
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

    def test_psomagen_per_sample_barcode_read_group_routes_to_rsq(self):
        """A barcode read group (not TT) must use the per-sample RSQ branch."""
        result = grab_trimmer_failure_codes_wafer_metrics(PSOM_SAMPLE_FAILURE)
        assert result is not None
        assert result["rsq_total_reads"] == 10695752970
        assert result["rsq_failed_reads"] == 3261649289
        # Whole-sample volume must NOT be mislabeled as TT.
        assert "tt_total_reads" not in result

    def test_novogene_per_sample_barcode_read_group_routes_to_rsq(self):
        result = grab_trimmer_failure_codes_wafer_metrics(NOVO_SAMPLE_FAILURE)
        assert result is not None
        assert result["rsq_total_reads"] == 499361927
        assert result["rsq_failed_reads"] == 215538084
        assert "tt_total_reads" not in result

    def test_merged_with_tt_row_routes_to_merged_parser(self):
        """A real TT row marks a merged file and yields TT counts."""
        result = grab_trimmer_failure_codes_wafer_metrics(NOVO_MERGED_FAILURE)
        assert result is not None
        assert result["tt_total_reads"] == 269168026
        assert result["tt_failed_reads"] == 31013293 + 7854


class TestMergePartialWaferStats:
    def test_sums_rsq_and_q30_across_libraries_for_same_run_id(self):
        merged: dict = {}
        merge_partial_wafer_stats(
            merged,
            "441049",
            {
                "rsq_total_reads": 500,
                "rsq_failed_reads": 50,
                "_sample_q30_matched_bases": 900,
                "_sample_q30_failures": 100,
            },
        )
        merge_partial_wafer_stats(
            merged,
            "441049",
            {
                "rsq_total_reads": 800,
                "rsq_failed_reads": 80,
                "_sample_q30_matched_bases": 1800,
                "_sample_q30_failures": 200,
            },
        )
        finalize_merged_wafer_stats(merged)
        stats = merged["441049"]
        # TT never comes from per-sample inputs.
        assert "tt_total_reads" not in stats
        assert stats["rsq_total_reads"] == 1300
        assert stats["rsq_failed_reads"] == 130
        assert stats["rsq_pass_reads"] == 1170
        assert stats["sample_q30_pct"] == pytest.approx(
            100.0 * (900 + 1800) / (900 + 100 + 1800 + 200), rel=1e-9
        )

    def test_skips_rsq_when_merged_failure_codes_present(self):
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
            {"rsq_total_reads": 999, "rsq_failed_reads": 99},
        )
        # Authoritative merged RSQ totals must not be double-counted.
        assert merged["439047"]["rsq_total_reads"] == 1000
        assert merged["439047"]["tt_total_reads"] == 10000
