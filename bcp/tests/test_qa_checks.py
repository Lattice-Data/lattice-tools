"""
Tests for qa_checks validation functions.
"""

import logging

from qa_checks import (
    MIN_METADATA_READ_COUNT,
    build_wafer_failure_stats,
    check_expected_raw_files,
    check_extra_raw_files,
    validate_fastq_counts,
    summarize_fastq_count_validation,
    validate_processed_group,
    validate_read_metadata,
    validate_scale_cb_tag,
    validate_scale_processed_files,
    validate_scale_samples_csv,
    validate_scale_workflow_info,
)


class TestValidateFastqCounts:
    """Tests for validate_fastq_counts."""

    def test_no_errors_when_gex_cri_match(self):
        """When GEX and CRI counts match, no errors."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "CRI": ["c", "d"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

    def test_error_when_gex_cri_mismatch(self):
        """When GEX and CRI counts differ, one error per sample."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "CRI": ["c"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert len(errors) == 1
        assert "MISMATCH FQ COUNTS" in errors[0]
        assert "2 GEX" in errors[0] and "1 CRI" in errors[0]
        assert "G1" in errors[0]

    def test_scale_gex_hash_oligo_mismatch(self):
        """Scale: error when GEX and hash_oligo counts differ."""
        fastq_log = {"E1": {"GEX": ["a"], "hash_oligo": ["b", "c"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert len(errors) == 1
        assert "GEX" in errors[0] and "hash_oligo" in errors[0]

    def test_scale_gex_hash_oligo_match(self):
        """Scale: no error when GEX and hash_oligo counts match."""
        fastq_log = {"E1": {"GEX": ["a", "b"], "hash_oligo": ["c", "d"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert errors == []

    def test_scale_gex_only_no_error(self):
        """Scale: only GEX present → no comparison, no error."""
        fastq_log = {"E1": {"GEX": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert errors == []

    def test_scale_hash_oligo_only_no_error(self):
        """Scale: only hash_oligo present → no comparison, no error."""
        fastq_log = {"E1": {"hash_oligo": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert errors == []

    def test_sci_plex_gex_hash_oligo_match(self):
        """sci_plex: no error when GEX and hash_oligo counts match."""
        fastq_log = {"S1": {"GEX": ["a", "b"], "hash_oligo": ["c", "d"]}}
        errors = validate_fastq_counts(fastq_log, "sci_plex")
        assert errors == []

    def test_sci_plex_gex_hash_oligo_mismatch(self):
        """sci_plex: error when GEX and hash_oligo both present and counts differ."""
        fastq_log = {"S1": {"GEX": ["a"], "hash_oligo": ["b", "c"]}}
        errors = validate_fastq_counts(fastq_log, "sci_plex")
        assert len(errors) == 1
        assert "GEX" in errors[0] and "hash_oligo" in errors[0]
        assert "S1" in errors[0]

    def test_sci_plex_gex_hash_oligo_only_no_error(self):
        """sci_plex: only GEX_hash_oligo present → no comparison, no error."""
        fastq_log = {"S1": {"GEX_hash_oligo": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "sci_plex")
        assert errors == []

    def test_sci_plex_gex_only_no_error(self):
        """sci_plex: only GEX present → no error."""
        fastq_log = {"S1": {"GEX": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "sci_plex")
        assert errors == []

    def test_sci_plex_hash_oligo_only_no_error(self):
        """sci_plex: only hash_oligo present → no error."""
        fastq_log = {"S1": {"hash_oligo": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "sci_plex")
        assert errors == []

    def test_sci_jumbo_no_validation_returns_no_errors(self, caplog):
        """sci_jumbo: no validation; returns no errors and logs message."""
        fastq_log = {"J1": {"GEX": ["a"], "hash_oligo": ["b", "c"]}}
        with caplog.at_level(logging.WARNING):
            errors = validate_fastq_counts(fastq_log, "sci_jumbo")
        assert errors == []
        assert "sci_jumbo" in caplog.text
        assert "not validating" in caplog.text.lower()

    def test_10x_gex_atac_match(self):
        """10x: GEX and ATAC counts match → no errors."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "ATAC": ["c", "d"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

    def test_10x_gex_atac_mismatch(self):
        """10x: GEX and ATAC counts differ → error."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "ATAC": ["c"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert len(errors) == 1
        assert "GEX" in errors[0] and "ATAC" in errors[0]

    def test_10x_gex_only_no_error(self):
        """10x: GEX only → no comparison, no error."""
        fastq_log = {"G1": {"GEX": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

    def test_10x_atac_only_no_error(self):
        """10x: ATAC only → no comparison, no error."""
        fastq_log = {"G1": {"ATAC": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

    def test_10x_viral_orf_present_no_validation(self):
        """10x: viral_ORF present is no longer validated; no error (legacy)."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "viral_ORF": ["c", "d", "e"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

    def test_10x_cri_atac_present_logs_warning(self, caplog):
        """10x: CRI+ATAC present logs unexpected / QA expansion needed."""
        fastq_log = {
            "G1": {"GEX": ["a", "b"], "CRI": ["c", "d"], "ATAC": ["e", "f"]},
        }
        with caplog.at_level(logging.WARNING):
            errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []
        assert "CRI" in caplog.text and "ATAC" in caplog.text
        assert "expansion" in caplog.text.lower() or "unexpected" in caplog.text.lower()

    def test_10x_viral_orf_assay_no_validation_returns_no_errors(self, caplog):
        """10x_viral_ORF: no validation; returns no errors and logs legacy message."""
        fastq_log = {"G1": {"GEX": ["a", "b"], "viral_ORF": ["c"]}}
        with caplog.at_level(logging.WARNING):
            errors = validate_fastq_counts(fastq_log, "10x_viral_ORF")
        assert errors == []
        assert "10x_viral_ORF" in caplog.text
        assert (
            "not validating" in caplog.text.lower() or "legacy" in caplog.text.lower()
        )


class TestSummarizeFastqCountValidation:
    """Tests for summarize_fastq_count_validation."""

    def test_scale_all_matched(self):
        fastq_log = {
            "E1": {"GEX": ["a", "b"], "hash_oligo": ["c", "d"]},
            "E2": {"GEX": ["x"], "hash_oligo": ["y"]},
        }
        errors = validate_fastq_counts(fastq_log, "scale")
        assert errors == []

        summary = summarize_fastq_count_validation(fastq_log, "scale", errors)
        assert "checked 2" in summary
        assert "mismatches: 0" in summary
        assert "All matched." in summary

    def test_scale_only_gex_present(self):
        fastq_log = {"E1": {"GEX": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert errors == []

        summary = summarize_fastq_count_validation(fastq_log, "scale", errors)
        assert "checked 0" in summary
        assert "mismatches: 0" in summary
        assert "Only GEX present" in summary

    def test_scale_mismatch_reports_mismatches(self):
        fastq_log = {"E1": {"GEX": ["a", "b"], "hash_oligo": ["c"]}}
        errors = validate_fastq_counts(fastq_log, "scale")
        assert len(errors) == 1

        summary = summarize_fastq_count_validation(fastq_log, "scale", errors)
        assert "mismatches: 1" in summary
        assert "Mismatches found" in summary

    def test_10x_all_matched(self):
        fastq_log = {
            "G1": {"GEX": ["a", "b"], "CRI": ["c", "d"], "ATAC": ["e", "f"]},
            "G2": {"GEX": ["x"], "CRI": ["y"], "ATAC": ["z"]},
        }
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

        summary = summarize_fastq_count_validation(fastq_log, "10x", errors)
        assert "Fastq count validation (10x)" in summary
        assert "checked 4" in summary
        assert "mismatches: 0" in summary
        assert "All matched." in summary

    def test_10x_only_gex_present(self):
        fastq_log = {"G1": {"GEX": ["a", "b"]}}
        errors = validate_fastq_counts(fastq_log, "10x")
        assert errors == []

        summary = summarize_fastq_count_validation(fastq_log, "10x", errors)
        assert "checked 0" in summary
        assert "mismatches: 0" in summary
        assert "Only GEX present" in summary

    def test_10x_cram_summary_not_applicable(self):
        fastq_log = {}
        errors = validate_fastq_counts(fastq_log, "10x_cram")
        assert errors == []
        summary = summarize_fastq_count_validation(fastq_log, "10x_cram", errors)
        assert "10x_cram" in summary
        assert "not applicable" in summary
        assert "CRAM-only" in summary


class TestValidateReadMetadata:
    """Tests for validate_read_metadata."""

    def test_builds_group_read_counts(self):
        """Builds group_read_counts from R1 metadata (no R2)."""
        read_metadata = {
            "run-G1_GEX-UG01-R1.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
            },
        }
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert errors == []
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []
        assert "G1" in counts
        assert counts["G1"]["GEX"] == MIN_METADATA_READ_COUNT

    def test_r1_r2_mismatch_error(self):
        """Error when R1 and R2 read counts differ."""
        # Use 10x-style filename (run-group_assay-ug-barcode pattern, 4 hyphen parts)
        read_metadata = {
            "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT + 100,
                "errors": [],
            },
            "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT + 200,
                "errors": [],
            },
        }
        counts, errors, _pairing = validate_read_metadata(read_metadata, "10x")
        assert any("READ COUNT ERROR" in e for e in errors)

    def test_r1_r2_success_print_summary(self, capsys):
        """When R1/R2 match, optionally print a success summary to stdout."""
        read_metadata = {
            "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
            },
            "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
            },
        }
        validate_read_metadata(
            read_metadata,
            "10x",
            print_success=True,
            success_print_limit=5,
        )
        out = capsys.readouterr().out
        assert "validate_read_metadata(10x):" in out
        assert "r1_r2_pairs_compared=1" in out
        assert "matched=1" in out
        assert "mismatched=0" in out
        assert "r2_missing_r1_metadata=0" in out
        assert "Checked metadata read counts (2):" in out
        assert (
            "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz: "
            f"read_count={MIN_METADATA_READ_COUNT}" in out
        )
        assert (
            "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz: "
            f"read_count={MIN_METADATA_READ_COUNT}" in out
        )
        assert "MATCH:" in out

    def test_read_count_equal_to_minimum_passes_without_low_read_error(self):
        read_metadata = {
            "run-G1_GEX-UG01-R1.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
            },
        }
        _counts, errors, _pairing = validate_read_metadata(read_metadata, "10x")
        assert not any("below minimum" in e for e in errors)

    def test_read_count_below_minimum_raises_error(self):
        read_metadata = {
            "run-G1_GEX-UG01-R1.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT - 1,
                "errors": [],
            },
        }
        _counts, errors, _pairing = validate_read_metadata(read_metadata, "10x")
        assert any(
            "READ COUNT ERROR:" in e
            and "below minimum" in e
            and str(MIN_METADATA_READ_COUNT) in e
            for e in errors
        )

    def test_r2_below_minimum_also_raises_error(self):
        read_metadata = {
            "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
            },
            "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT - 1,
                "errors": [],
            },
        }
        _counts, errors, _pairing = validate_read_metadata(read_metadata, "10x")
        assert any(
            "READ COUNT ERROR:" in e and "R2_001.fastq.gz" in e and "below minimum" in e
            for e in errors
        )

    def test_r1_without_r2_lists_path_and_errors(self):
        r1 = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        read_metadata = {r1: {"read_count": MIN_METADATA_READ_COUNT, "errors": []}}
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert pairing["r1_without_r2_metadata"] == [r1]
        assert pairing["r2_without_r1_metadata"] == []
        assert any(r1 in e and "no R2 metadata" in e for e in errors)

    def test_r2_without_r1_lists_path_and_errors(self):
        r2 = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        read_metadata = {r2: {"read_count": MIN_METADATA_READ_COUNT, "errors": []}}
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert pairing["r2_without_r1_metadata"] == [r2]
        assert pairing["r1_without_r2_metadata"] == []
        assert any(r2 in e and "no R1 metadata" in e for e in errors)

    def test_metadata_errors_appended(self):
        """Metadata with errors list adds to errors."""
        read_metadata = {
            "run-G1_GEX-UG01_R1.fastq.gz": {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": ["bad"],
            },
        }
        counts, errors, _p = validate_read_metadata(read_metadata, "10x")
        assert any("METADATA.JSON ERROR" in e for e in errors)
        assert "G1" not in counts or counts.get("G1", {}).get("GEX") != 100

    def test_r2_in_group_id_pairs_correctly(self):
        """R2 in the group name (q_pcf_R2) must not confuse R1/R2 pairing."""
        r1 = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R1_001.fastq.gz"
        r2 = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R2_001.fastq.gz"
        read_metadata = {
            r1: {"read_count": MIN_METADATA_READ_COUNT + 500, "errors": []},
            r2: {"read_count": MIN_METADATA_READ_COUNT + 500, "errors": []},
        }
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []
        assert not any("PAIRING" in e for e in errors)
        assert "q_pcf_R2" in counts
        assert counts["q_pcf_R2"]["GEX"] == MIN_METADATA_READ_COUNT + 500

    def test_r1_in_group_id_pairs_correctly(self):
        """R1 in the group name (q_hf_R1) must not confuse R1/R2 pairing."""
        r1 = "439925-q_hf_R1_GEX-Z0004-CTGTGTAGGCATGAT_S1_L001_R1_001.fastq.gz"
        r2 = "439925-q_hf_R1_GEX-Z0004-CTGTGTAGGCATGAT_S1_L001_R2_001.fastq.gz"
        read_metadata = {
            r1: {"read_count": MIN_METADATA_READ_COUNT + 300, "errors": []},
            r2: {"read_count": MIN_METADATA_READ_COUNT + 300, "errors": []},
        }
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []
        assert not any("PAIRING" in e for e in errors)
        assert "q_hf_R1" in counts
        assert counts["q_hf_R1"]["GEX"] == MIN_METADATA_READ_COUNT + 300

    def test_r2_in_group_id_mismatch_detected(self):
        """Read count mismatch is still detected with R2 in the group name."""
        r1 = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R1_001.fastq.gz"
        r2 = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R2_001.fastq.gz"
        read_metadata = {
            r1: {"read_count": MIN_METADATA_READ_COUNT + 500, "errors": []},
            r2: {"read_count": MIN_METADATA_READ_COUNT + 999, "errors": []},
        }
        _counts, errors, _pairing = validate_read_metadata(read_metadata, "10x")
        assert any("READ COUNT ERROR" in e for e in errors)

    def test_metadata_filename_mismatch_warns_but_pairing_uses_actual_filename(self):
        """Filename mismatch in metadata content is warning-only.

        Pairing/count comparison continues using the canonical metadata key.
        """
        actual_r1 = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        actual_r2 = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        read_metadata = {
            actual_r1: {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
                "__actual_filename": actual_r1,
                "__reported_filename": "s3://bucket/WRONG_R1.fastq.gz",
            },
            actual_r2: {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
                "__actual_filename": actual_r2,
                "__reported_filename": actual_r2,
            },
        }
        counts, errors, pairing = validate_read_metadata(read_metadata, "10x")
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []
        assert any("METADATA FILENAME WARNING:" in e for e in errors)
        assert not any("READ METADATA PAIRING" in e for e in errors)
        assert counts["G1"]["GEX"] == MIN_METADATA_READ_COUNT

    def test_mismatch_warning_printed_in_summary(self, capsys):
        """Summary line includes mismatch warning count for observability."""
        actual_r1 = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        actual_r2 = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        read_metadata = {
            actual_r1: {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
                "__actual_filename": actual_r1,
                "__reported_filename": "s3://bucket/WRONG_R1.fastq.gz",
            },
            actual_r2: {
                "read_count": MIN_METADATA_READ_COUNT,
                "errors": [],
                "__actual_filename": actual_r2,
                "__reported_filename": actual_r2,
            },
        }
        validate_read_metadata(read_metadata, "10x", print_success=True)
        out = capsys.readouterr().out
        assert "metadata_filename_mismatch_warnings=1" in out

    def test_10x_cram_skips_r1_r2_pairing_and_summary(self, capsys):
        cram = "442488-LeS1867W11_ATAC-Z0027-CACTGTCAGCCAGAT.cram"
        read_metadata = {
            cram: {"read_count": MIN_METADATA_READ_COUNT, "errors": []},
        }
        counts, errors, pairing = validate_read_metadata(
            read_metadata, "10x_cram", print_success=True
        )
        out = capsys.readouterr().out
        assert "validate_read_metadata(10x_cram):" in out
        assert "CRAM-only layout" in out
        assert "r1_r2_pairs_compared" not in out
        assert "Checked metadata read counts (1):" in out
        assert f"{cram}: read_count={MIN_METADATA_READ_COUNT}" in out
        assert errors == []
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []
        assert counts["LeS1867W11"]["ATAC"] == MIN_METADATA_READ_COUNT


class Test10xCramRawFiles:
    """Tests for expected/extra raw files under 10x_cram policy."""

    def test_10x_cram_expected_and_optional_files(self):
        base = "proj/order/LeS1867W11/raw/442488-LeS1867W11_ATAC-Z0027-CACTGTCAGCCAGAT"
        all_raw = [
            f"{base}.cram",
            f"{base}.cram-metadata.json",
            f"{base}.csv",
            f"{base}.json",
            f"{base}_FlowQ.metric",
            f"{base}_SNVQ.metric",
            f"{base}_trimmer-failure_codes.csv",
            f"{base}_trimmer-stats.csv",
            f"{base}_extract_stats.h5",
        ]
        all_good, raw_lost, raw_found = check_expected_raw_files(all_raw, "10x_cram")
        assert all_good == 1
        assert raw_lost == []
        assert len(raw_found) == 8
        extra = check_extra_raw_files(all_raw, raw_found, "10x_cram")
        assert extra == []


class TestCheckExpectedRawFiles:
    """Tests for check_expected_raw_files."""

    def test_10x_all_present(self):
        """When all expected 10x endings are present, all_good=1."""
        from qa_mods import raw_expected

        base = "proj/order/G1/raw/439047-G1_GEX-Z0273-BC01"
        endings = raw_expected["10x"]
        all_raw = [f"{base}{e}" for e in endings]
        all_good, raw_lost, raw_found = check_expected_raw_files(all_raw, "10x")
        assert all_good == 1
        assert raw_lost == []
        assert len(raw_found) == len(all_raw)

    def test_10x_one_missing(self):
        """When one expected file is missing, raw_lost has one entry."""
        base = "proj/order/G1/raw/439047-G1_GEX-Z0273-BC01"
        all_raw = [
            f"{base}.csv",
            f"{base}_trimmer-failure_codes.csv",
            # omit .json and others
        ]
        all_good, raw_lost, raw_found = check_expected_raw_files(all_raw, "10x")
        assert all_good == 0 or len(raw_lost) >= 1
        assert len(raw_found) <= len(all_raw)


class TestCheckExtraRawFiles:
    """Tests for check_extra_raw_files."""

    def test_no_extra_when_all_in_raw_found(self):
        """When all_raw_files equals raw_found, extra is empty."""
        files = ["proj/order/G1/raw/439047-G1_GEX-Z0273-BC01.csv"]
        extra = check_extra_raw_files(files, files, "10x")
        assert extra == []

    def test_optional_10x_not_extra(self):
        """Optional 10x endings (e.g. .scRNA.applicationQC.h5) are not reported as extra."""
        base = "proj/order/G1/raw/439047-G1_GEX-Z0273-BC01"
        expected_file = f"{base}.csv"
        optional_file = f"{base}.scRNA.applicationQC.h5"
        all_raw = [expected_file, optional_file]
        raw_found = [expected_file]
        extra = check_extra_raw_files(all_raw, raw_found, "10x")
        assert optional_file not in extra

    def test_unknown_file_is_extra(self):
        """File not in raw_found and not optional is extra."""
        base = "proj/order/G1/raw/439047-G1_GEX-Z0273-BC01"
        all_raw = [f"{base}.csv", f"{base}_unknown_suffix.xyz"]
        raw_found = [f"{base}.csv"]
        extra = check_extra_raw_files(all_raw, raw_found, "10x")
        assert f"{base}_unknown_suffix.xyz" in extra

    # --- Scale-specific extra-file tests ---

    _SCALE_BASE = "proj/NVUS123/GENE9-R115/raw/440115"

    def test_scale_per_rt_gex_cram_not_extra(self):
        """Per-RT GEX CRAM files are recognized as known Scale files."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.cram"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_per_rt_gex_csv_not_extra(self):
        """Per-RT GEX CSV files are recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_per_rt_gex_json_not_extra(self):
        """Per-RT GEX JSON files are recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.json"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_per_rt_hash_oligo_cram_not_extra(self):
        """Per-RT hash_oligo CRAM files are recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-1A.cram"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_per_rt_double_digit_well_not_extra(self):
        """Per-RT files with double-digit row (e.g. 12H) are recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-12H.cram"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_per_rt_underscore_well_not_extra(self):
        """Production-style per-RT names use underscore before well (e.g. _5B)."""
        raw = [
            f"{self._SCALE_BASE}/426971-RNA3-098C_GEX_QSR-7_5B.cram",
            f"{self._SCALE_BASE}/426971-RNA3-098C_hash_oligo_QSR-7-SCALEPLEX_5B.json",
        ]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_trimmer_failure_codes_not_extra(self):
        """Aggregate trimmer-failure-codes.csv is recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-trimmer-failure-codes.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_trimmer_stats_not_extra(self):
        """Aggregate trimmer-stats.csv is recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-trimmer-stats.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_trimmer_failure_codes_underscore_form_not_extra(self):
        """Production uses _trimmer-failure_codes.csv (underscore before codes)."""
        raw = [
            f"{self._SCALE_BASE}/426971-RNA3-098C_GEX_QSR-7_trimmer-failure_codes.csv",
            f"{self._SCALE_BASE}/426971-RNA3-098C_GEX_QSR-7_trimmer-stats.csv",
            f"{self._SCALE_BASE}/426971-RNA3-098C_GEX_QSR-7_unmatched.cram",
        ]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_unmatched_cram_not_extra(self):
        """Aggregate unmatched.cram is recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-unmatched.cram"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_unmatched_csv_not_extra(self):
        """Aggregate unmatched.csv is recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-unmatched.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_aggregate_unmatched_json_not_extra(self):
        """Aggregate unmatched.json is recognized."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-unmatched.json"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_hash_oligo_aggregate_trimmer_not_extra(self):
        """hash_oligo aggregate trimmer files are recognized."""
        raw = [
            f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-trimmer-failure-codes.csv",
            f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-trimmer-stats.csv",
        ]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_hash_oligo_aggregate_unmatched_not_extra(self):
        """hash_oligo aggregate unmatched files are recognized."""
        raw = [
            f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.cram",
            f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.csv",
            f"{self._SCALE_BASE}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.json",
        ]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_wafer_sequencing_info_not_extra(self):
        """Wafer-level SequencingInfo.json is recognized."""
        raw = [f"{self._SCALE_BASE}/440115_SequencingInfo.json"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_wafer_library_info_not_extra(self):
        """Wafer-level LibraryInfo.xml is recognized."""
        raw = [f"{self._SCALE_BASE}/440115_LibraryInfo.xml"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_merged_trimmer_failure_codes_not_extra(self):
        """Merged trimmer-failure_codes.csv is recognized."""
        raw = [f"{self._SCALE_BASE}/merged_trimmer-failure_codes.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_merged_trimmer_stats_not_extra(self):
        """Merged trimmer-stats.csv is recognized."""
        raw = [f"{self._SCALE_BASE}/merged_trimmer-stats.csv"]
        assert check_extra_raw_files(raw, [], "scale") == []

    def test_scale_cram_metadata_sidecar_not_extra(self):
        """CRAM metadata sidecar is handled by existing sidecar logic."""
        cram = f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.cram"
        sidecar = f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.cram-metadata.json"
        assert check_extra_raw_files([cram, sidecar], [], "scale") == []

    def test_scale_unknown_extension_is_extra(self):
        """Unrecognized file extension IS flagged as extra."""
        raw = [f"{self._SCALE_BASE}/440115-R115H_GEX_QSR-8-5B.bam"]
        extra = check_extra_raw_files(raw, [], "scale")
        assert len(extra) == 1

    def test_scale_unknown_filename_is_extra(self):
        """Completely unrecognized filename IS flagged as extra."""
        raw = [f"{self._SCALE_BASE}/random_unknown_file.txt"]
        extra = check_extra_raw_files(raw, [], "scale")
        assert len(extra) == 1

    def test_scale_comprehensive_wafer(self):
        """A realistic set of Scale raw files: nothing should be extra."""
        base = self._SCALE_BASE
        files = [
            f"{base}/440115-R115H_GEX_QSR-8-1A.cram",
            f"{base}/440115-R115H_GEX_QSR-8-1A.csv",
            f"{base}/440115-R115H_GEX_QSR-8-1A.json",
            f"{base}/440115-R115H_GEX_QSR-8-1A.cram-metadata.json",
            f"{base}/440115-R115H_GEX_QSR-8-12H.cram",
            f"{base}/440115-R115H_GEX_QSR-8-12H.csv",
            f"{base}/440115-R115H_GEX_QSR-8-12H.json",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-1A.cram",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-1A.csv",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-1A.json",
            f"{base}/440115-R115H_GEX_QSR-8-trimmer-failure-codes.csv",
            f"{base}/440115-R115H_GEX_QSR-8-trimmer-stats.csv",
            f"{base}/440115-R115H_GEX_QSR-8-unmatched.cram",
            f"{base}/440115-R115H_GEX_QSR-8-unmatched.csv",
            f"{base}/440115-R115H_GEX_QSR-8-unmatched.json",
            f"{base}/440115-R115H_GEX_QSR-8-unmatched.cram-metadata.json",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-trimmer-failure-codes.csv",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-trimmer-stats.csv",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.cram",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.csv",
            f"{base}/440115-R115H_hash_oligo_QSR-8-SCALEPLEX-unmatched.json",
            f"{base}/440115_SequencingInfo.json",
            f"{base}/440115_LibraryInfo.xml",
            f"{base}/merged_trimmer-failure_codes.csv",
            f"{base}/merged_trimmer-stats.csv",
        ]
        assert check_extra_raw_files(files, [], "scale") == []

    def test_scale_comprehensive_wafer_with_unknown_file(self):
        """A realistic set plus one unknown file: only the unknown is extra."""
        base = self._SCALE_BASE
        unknown = f"{base}/unexpected_report.pdf"
        files = [
            f"{base}/440115-R115H_GEX_QSR-8-1A.cram",
            f"{base}/440115-R115H_GEX_QSR-8-trimmer-failure-codes.csv",
            f"{base}/440115_SequencingInfo.json",
            f"{base}/merged_trimmer-failure_codes.csv",
            unknown,
        ]
        extra = check_extra_raw_files(files, [], "scale")
        assert extra == [unknown]


class TestValidateProcessedGroup:
    """Tests for validate_processed_group."""

    def test_unknown_software_error(self):
        """Unknown software version produces an error."""
        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-99.0.0",
            "gex_alerts": [],
        }
        proc_files = ["proj/g/processed/cellranger/run/outs/config.csv"]
        result = validate_processed_group("G1", proc_files, report, {})
        assert any("CR ERROR" in e and "9.0.1 or 10.0.0" in e for e in result["errors"])

    def test_min_crispr_umi_error(self):
        """min-crispr-umi != 3 produces error."""
        report = {
            "chem": "flex",
            "extra": ["CRISPR"],
            "software": "cellranger-10.0.0",
            "gex_alerts": [],
            "min-crispr-umi": "5",
        }
        proc_files = ["proj/g/processed/cellranger/run/outs/config.csv"]
        result = validate_processed_group("G1", proc_files, report, {})
        assert any("min-crispr-umi" in e for e in result["errors"])

    def test_incl_int_error(self):
        """include-introns not true produces error."""
        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-9.0.1",
            "sub": "count",
            "gex_alerts": [],
            "incl_int": "false",
        }
        proc_files = ["proj/g/processed/cellranger/run/outs/config.csv"]
        result = validate_processed_group("G1", proc_files, report, {})
        assert any("include-introns" in e for e in result["errors"])

    def test_read_count_reconciliation(self):
        """When report GEX_reads != group_read_counts, error is added."""
        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-10.0.0",
            "gex_alerts": [],
            "GEX_reads": 1000,
        }
        proc_files = ["proj/g/processed/cellranger/run/outs/config.csv"]
        group_read_counts = {"G1": {"GEX": 500}}
        result = validate_processed_group("G1", proc_files, report, group_read_counts)
        assert any("READ COUNT ERROR" in e for e in result["errors"])

    def test_alerts_populated(self):
        """gex_alerts and crispr_alerts are returned as alerts with group/modality."""
        report = {
            "chem": "flex",
            "extra": ["CRISPR"],
            "software": "cellranger-10.0.0",
            "gex_alerts": [{"id": "gex1", "message": "low"}],
            "crispr_alerts": [{"id": "cri1"}],
        }
        proc_files = ["proj/g/processed/cellranger/run/outs/config.csv"]
        result = validate_processed_group("G1", proc_files, report, {})
        assert len(result["alerts"]) == 2
        gex_alert = next(a for a in result["alerts"] if a.get("modality") == "GEX")
        cri_alert = next(a for a in result["alerts"] if a.get("modality") == "CRI")
        assert gex_alert["group"] == "G1"
        assert cri_alert["group"] == "G1"

    def test_missing_expected_file(self):
        """When an expected outs file is missing, proc_missing is populated."""
        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-10.0.0",
            "gex_alerts": [],
        }
        # Only config.csv present; cellranger 10.0.0 nonflex expects many more
        proc_files = [
            "proj/g/processed/cellranger/run/outs/config.csv",
        ]
        result = validate_processed_group("G1", proc_files, report, {})
        assert len(result["proc_missing"]) >= 1
        assert any("group" in m and m["group"] == "G1" for m in result["proc_missing"])

    def test_count_index_csi_only_is_accepted(self):
        """When BAI is missing but CSI is present, the run is accepted."""
        from qa_constants import cellranger_expected

        bai = "possorted_genome_bam.bam.bai"
        csi = "possorted_genome_bam.bam.csi"
        expected_outs = list(cellranger_expected["count"]["outs"])
        assert bai in expected_outs
        assert csi not in expected_outs

        base = "proj/g/processed/cellranger/run/outs/"
        proc_files = [f"{base}{f}" for f in expected_outs if f != bai]
        proc_files.append(f"{base}{csi}")

        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-9.0.1",
            "sub": "count",
            "gex_alerts": [],
        }
        result = validate_processed_group("G1", proc_files, report, {})
        assert result["errors"] == []
        assert result["proc_missing"] == []
        assert result["process_extra"] == []

    def test_count_index_neither_bai_nor_csi_is_reported_as_missing(self):
        """When neither BAI nor CSI exists, both are flagged in proc_missing."""
        from qa_constants import cellranger_expected

        bai = "possorted_genome_bam.bam.bai"
        csi = "possorted_genome_bam.bam.csi"
        expected_outs = list(cellranger_expected["count"]["outs"])

        base = "proj/g/processed/cellranger/run/outs/"
        proc_files = [f"{base}{f}" for f in expected_outs if f != bai]
        # Do not add CSI either

        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-9.0.1",
            "sub": "count",
            "gex_alerts": [],
        }
        result = validate_processed_group("G1", proc_files, report, {})

        assert any(
            "missing both" in e and bai in e and csi in e for e in result["errors"]
        )
        assert len(result["proc_missing"]) == 1
        row = result["proc_missing"][0]
        assert row["group"] == "G1"
        assert row[bai] == "Y"
        assert row[csi] == "Y"

    def test_count_index_both_bai_and_csi_is_an_error(self):
        """When both BAI and CSI exist, it is an error (not accepted)."""
        from qa_constants import cellranger_expected

        bai = "possorted_genome_bam.bam.bai"
        csi = "possorted_genome_bam.bam.csi"
        expected_outs = list(cellranger_expected["count"]["outs"])
        assert bai in expected_outs

        base = "proj/g/processed/cellranger/run/outs/"
        proc_files = [f"{base}{f}" for f in expected_outs]
        proc_files.append(f"{base}{csi}")

        report = {
            "chem": "5p",
            "extra": [],
            "software": "cellranger-9.0.1",
            "sub": "count",
            "gex_alerts": [],
        }
        result = validate_processed_group("G1", proc_files, report, {})

        assert any("has both" in e and bai in e and csi in e for e in result["errors"])
        assert len(result["proc_missing"]) == 1
        row = result["proc_missing"][0]
        assert row["group"] == "G1"
        assert row[bai] == "Y"
        assert row[csi] == "Y"
        assert result["process_extra"] == []


class TestBuildWaferFailureStats:
    """Tests for build_wafer_failure_stats."""

    def test_aggregates_multiple_experiments_into_single_wafer(self):
        trimmer_failure_stats = {
            "order1/groupA": {"rsq": [1.0, 2.0], "trimmer_fail": [5.0]},
            "order1/groupB": {"rsq": [3.0], "trimmer_fail": [10.0]},
        }
        exp_to_run_map = {
            "order1/groupA": "439047",
            "order1/groupB": "439047",
        }
        wafer_stats = build_wafer_failure_stats(trimmer_failure_stats, exp_to_run_map)
        assert set(wafer_stats.keys()) == {"439047"}
        assert wafer_stats["439047"]["rsq"] == [1.0, 2.0, 3.0]
        assert wafer_stats["439047"]["trimmer_fail"] == [5.0, 10.0]

    def test_skips_experiments_without_run_mapping(self):
        trimmer_failure_stats = {
            "order1/groupA": {"rsq": [1.0], "trimmer_fail": [5.0]},
            "order1/groupB": {"rsq": [2.0], "trimmer_fail": [6.0]},
        }
        exp_to_run_map = {"order1/groupA": "439047"}  # groupB unmapped
        wafer_stats = build_wafer_failure_stats(trimmer_failure_stats, exp_to_run_map)
        assert set(wafer_stats.keys()) == {"439047"}
        assert wafer_stats["439047"]["rsq"] == [1.0]
        assert wafer_stats["439047"]["trimmer_fail"] == [5.0]


class TestValidateScaleWorkflowInfo:
    """Tests for validate_scale_workflow_info."""

    def test_all_params_correct_no_errors(self):
        """No errors when all required params match expected values."""
        params = {
            "bamOut": "true",
            "scalePlex": "true",
            "scalePlexAssignmentMethod": "fc",
            "genome": "s3://bucket/ref/GRCh38.json",
        }
        errors, info = validate_scale_workflow_info(params)
        assert errors == []
        assert any("SCALE GENOME" in m for m in info)

    def test_bamout_false_is_error(self):
        """Error when bamOut != 'true'."""
        params = {
            "bamOut": "false",
            "scalePlex": "true",
            "scalePlexAssignmentMethod": "fc",
            "genome": "s3://bucket/ref/GRCh38.json",
        }
        errors, _ = validate_scale_workflow_info(params)
        assert len(errors) == 1
        assert "bamOut" in errors[0]
        assert "'false'" in errors[0]

    def test_scaleplex_false_is_error(self):
        """Error when scalePlex != 'true'."""
        params = {
            "bamOut": "true",
            "scalePlex": "false",
            "scalePlexAssignmentMethod": "fc",
            "genome": "s3://bucket/ref/GRCh38.json",
        }
        errors, _ = validate_scale_workflow_info(params)
        assert len(errors) == 1
        assert "scalePlex" in errors[0]

    def test_assignment_method_wrong_is_error(self):
        """Error when scalePlexAssignmentMethod != 'fc'."""
        params = {
            "bamOut": "true",
            "scalePlex": "true",
            "scalePlexAssignmentMethod": "other",
            "genome": "s3://bucket/ref/GRCh38.json",
        }
        errors, _ = validate_scale_workflow_info(params)
        assert len(errors) == 1
        assert "scalePlexAssignmentMethod" in errors[0]

    def test_genome_in_info_messages(self):
        """Genome value appears in info messages."""
        params = {
            "bamOut": "true",
            "scalePlex": "true",
            "scalePlexAssignmentMethod": "fc",
            "genome": "s3://bucket/ref/GRCh38.json",
        }
        _, info = validate_scale_workflow_info(params)
        assert any("GRCh38" in m for m in info)


class TestValidateScaleSamplesCsv:
    """Tests for validate_scale_samples_csv."""

    def test_no_forbidden_columns_passes(self):
        """No errors when forbidden columns are absent."""
        columns = ["sample", "barcodes", "libIndex2", "scalePlexLibIndex2"]
        errors = validate_scale_samples_csv(columns)
        assert errors == []

    def test_scaleplexbarcodes_column_is_error(self):
        """Error when scalePlexBarcodes column is present."""
        columns = ["sample", "barcodes", "libIndex2", "scalePlexBarcodes"]
        errors = validate_scale_samples_csv(columns)
        assert len(errors) == 1
        assert "scalePlexBarcodes" in errors[0]


class TestValidateScaleCbTag:
    """Tests for validate_scale_cb_tag."""

    def test_all_true_passes(self):
        """No errors when all CRAM files have cb_tag=True."""
        read_metadata = {
            "sample1.cram": {"cb_tag": True},
            "sample2.cram": {"cb_tag": True},
        }
        errors = validate_scale_cb_tag(read_metadata)
        assert errors == []

    def test_false_cb_tag_is_error(self):
        """Error reported when cb_tag is False."""
        read_metadata = {
            "sample1.cram": {"cb_tag": True},
            "sample2.cram": {"cb_tag": False},
        }
        errors = validate_scale_cb_tag(read_metadata)
        assert len(errors) == 1
        assert "sample2.cram" in errors[0]
        assert "cb_tag=False" in errors[0]

    def test_missing_cb_tag_key_is_error(self):
        """Error when cb_tag key is absent from metadata."""
        read_metadata = {
            "sample1.cram": {"read_count": 100},
        }
        errors = validate_scale_cb_tag(read_metadata)
        assert len(errors) == 1
        assert "sample1.cram" in errors[0]
        assert "cb_tag=None" in errors[0]

    def test_unmatched_cram_skipped(self):
        """Unmatched CRAMs (both hyphen and underscore forms) are not checked for cb_tag."""
        read_metadata = {
            "sample1_GEX_QSR-5-unmatched.cram": {"cb_tag": False},
            "sample1_hash_oligo_QSR-5-SCALEPLEX-unmatched.cram": {"cb_tag": False},
            "sample1_GEX_QSR-5_unmatched.cram": {"cb_tag": False},
            "sample1_hash_oligo_QSR-5-SCALEPLEX_unmatched.cram": {"cb_tag": False},
            "sample1_GEX_QSR-5.cram": {"cb_tag": True},
        }
        errors = validate_scale_cb_tag(read_metadata)
        assert errors == []

    def test_cram_metadata_json_skipped(self):
        """Files ending in .cram-metadata.json are not checked."""
        read_metadata = {
            "sample1.cram-metadata.json": {"cb_tag": False},
        }
        errors = validate_scale_cb_tag(read_metadata)
        assert errors == []


class TestValidateScaleProcessedFiles:
    """Tests for validate_scale_processed_files."""

    def test_all_files_present(self):
        """No missing files when all expected per-sample and sublibrary files exist."""
        samples_info = {
            "samples": ["S1"],
            "sublibraries": {"S1": ["QSR-1", "QSR-2"]},
        }
        proc_files = [
            "proj/group/processed/run/samples/S1_anndata.h5ad",
            "proj/group/processed/run/samples/S1.merged.allCells.csv",
            "proj/group/processed/run/samples/S1.QSR-1_anndata.h5ad",
            "proj/group/processed/run/samples/S1.QSR-2_anndata.h5ad",
        ]
        result = validate_scale_processed_files(proc_files, samples_info)
        assert result["errors"] == []
        assert result["missing_files"] == []

    def test_merged_anndata_missing(self):
        """Merged anndata missing is flagged."""
        samples_info = {
            "samples": ["S1"],
            "sublibraries": {"S1": ["QSR-1"]},
        }
        proc_files = [
            "proj/group/processed/run/samples/S1.merged.allCells.csv",
            "proj/group/processed/run/samples/S1.QSR-1_anndata.h5ad",
        ]
        result = validate_scale_processed_files(proc_files, samples_info)
        assert any("S1_anndata.h5ad" in e for e in result["errors"])
        assert any(m["file"] == "S1_anndata.h5ad" for m in result["missing_files"])

    def test_sublibrary_anndata_missing(self):
        """Per-sublibrary anndata missing is flagged."""
        samples_info = {
            "samples": ["S1"],
            "sublibraries": {"S1": ["QSR-1", "QSR-2"]},
        }
        proc_files = [
            "proj/group/processed/run/samples/S1_anndata.h5ad",
            "proj/group/processed/run/samples/S1.merged.allCells.csv",
            "proj/group/processed/run/samples/S1.QSR-1_anndata.h5ad",
        ]
        result = validate_scale_processed_files(proc_files, samples_info)
        assert any("S1.QSR-2_anndata.h5ad" in e for e in result["errors"])
        assert any(
            m["file"] == "S1.QSR-2_anndata.h5ad" for m in result["missing_files"]
        )

    def test_allcells_csv_missing(self):
        """Merged allCells.csv missing is flagged."""
        samples_info = {
            "samples": ["S1"],
            "sublibraries": {"S1": ["QSR-1"]},
        }
        proc_files = [
            "proj/group/processed/run/samples/S1_anndata.h5ad",
            "proj/group/processed/run/samples/S1.QSR-1_anndata.h5ad",
        ]
        result = validate_scale_processed_files(proc_files, samples_info)
        assert any("S1.merged.allCells.csv" in e for e in result["errors"])
        assert any(
            m["file"] == "S1.merged.allCells.csv" for m in result["missing_files"]
        )
