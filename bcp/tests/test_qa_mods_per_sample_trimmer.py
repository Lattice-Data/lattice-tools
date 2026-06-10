"""Tests for per-sample trimmer-stats parsing and wafer aggregation helpers."""

from __future__ import annotations

import os

from qa_mods import (
    finalize_merged_wafer_stats,
    grab_trimmer_failure_codes_wafer_metrics,
    merge_partial_wafer_stats,
)

QA_FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "qa")
PSOM_SAMPLE_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "psomagen_per_sample_trimmer_failure_codes.csv"
)
NOVO_SAMPLE_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "novogene_per_sample_trimmer_failure_codes.csv"
)
NOVO_MERGED_FAILURE = os.path.join(
    QA_FIXTURES_DIR, "novogene_merged_trimmer_failure_codes_sample.csv"
)
NOVO_MERGED_MULTIREASON = os.path.join(
    QA_FIXTURES_DIR, "novogene_merged_trimmer_failure_codes_multireason.csv"
)


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

    def test_merged_multireason_tt_sums_all_reasons_uses_rsq_total(self):
        """TT block with rsq file + too short + too long (divergent total).

        tt_total must come from the read-group total (the rsq-file row), not the
        smaller 'too long' segment total; tt_failed sums failures across all
        three TT reasons.
        """
        result = grab_trimmer_failure_codes_wafer_metrics(NOVO_MERGED_MULTIREASON)
        assert result is not None
        assert result["tt_total_reads"] == 369936374
        assert result["tt_failed_reads"] == 90182967 + 4653 + 2
        # RSQ totals pool every rsq-file row (TT + per-sample groups).
        assert result["rsq_total_reads"] == 369936374 + 32297 + 3100000 + 2800000
        assert result["rsq_failed_reads"] == 90182967 + 25633 + 1200000 + 900000

    def test_merged_tt_total_robust_to_row_ordering(self, tmp_path):
        """tt_total uses the max TT total even if the smaller segment is first."""
        path = tmp_path / "merged.csv"
        path.write_text(
            "read group,code,format,segment,reason,failed read count,total read count\n"
            # 'too long' (smaller total) listed FIRST to defeat positional logic.
            "TT,102,trim,insert,sequence was too long,2,175771633\n"
            "TT,8,trim,preamble,rsq file,90182967,369936374\n"
            "TT,101,trim,insert,sequence was too short,4653,369936374\n"
        )
        result = grab_trimmer_failure_codes_wafer_metrics(str(path))
        assert result is not None
        assert result["tt_total_reads"] == 369936374
        assert result["tt_failed_reads"] == 2 + 90182967 + 4653


class TestMergePartialWaferStats:
    def test_sums_rsq_across_libraries_for_same_run_id(self):
        merged: dict = {}
        merge_partial_wafer_stats(
            merged,
            "441049",
            {"rsq_total_reads": 500, "rsq_failed_reads": 50},
        )
        merge_partial_wafer_stats(
            merged,
            "441049",
            {"rsq_total_reads": 800, "rsq_failed_reads": 80},
        )
        finalize_merged_wafer_stats(merged)
        stats = merged["441049"]
        # TT never comes from per-sample inputs.
        assert "tt_total_reads" not in stats
        assert stats["rsq_total_reads"] == 1300
        assert stats["rsq_failed_reads"] == 130
        assert stats["rsq_pass_reads"] == 1170
        # Q30 is no longer derived from per-sample trimmer files.
        assert "sample_q30_pct" not in stats

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
