"""
Tests for qa_checks validation functions.
"""

import logging

from qa_checks import (
    build_wafer_failure_stats,
    check_expected_raw_files,
    check_extra_raw_files,
    validate_fastq_counts,
    summarize_fastq_count_validation,
    validate_processed_group,
    validate_read_metadata,
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


class TestValidateReadMetadata:
    """Tests for validate_read_metadata."""

    def test_builds_group_read_counts(self):
        """Builds group_read_counts from R1 metadata (no R2)."""
        read_metadata = {
            "run-G1_GEX-UG01-R1.fastq.gz": {"read_count": 100, "errors": []},
        }
        counts, errors = validate_read_metadata(read_metadata, "10x")
        assert errors == []
        assert "G1" in counts
        assert counts["G1"]["GEX"] == 100

    def test_r1_r2_mismatch_error(self):
        """Error when R1 and R2 read counts differ."""
        # Use 10x-style filename (run-group_assay-ug-barcode pattern, 4 hyphen parts)
        read_metadata = {
            "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz": {
                "read_count": 100,
                "errors": [],
            },
            "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz": {
                "read_count": 200,
                "errors": [],
            },
        }
        counts, errors = validate_read_metadata(read_metadata, "10x")
        assert any("READ COUNT ERROR" in e for e in errors)

    def test_metadata_errors_appended(self):
        """Metadata with errors list adds to errors."""
        read_metadata = {
            "run-G1_GEX-UG01_R1.fastq.gz": {"read_count": 100, "errors": ["bad"]},
        }
        counts, errors = validate_read_metadata(read_metadata, "10x")
        assert any("METADATA.JSON ERROR" in e for e in errors)
        assert "G1" not in counts or counts.get("G1", {}).get("GEX") != 100


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
