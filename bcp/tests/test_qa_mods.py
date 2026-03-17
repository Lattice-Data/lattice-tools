"""
Tests for qa_mods parsing functions: parse_raw_filename, parse_met_summ, parse_web_summ.
Fixture files in tests/fixtures/qa/ are real outputs from CZI sequencing runs
and 10x Genomics CellRanger pipelines.
"""

import os

import pytest

from qa_mods import (
    extract_run_id_from_trimmer_filename,
    grab_trimmer_stats,
    parse_met_summ,
    parse_raw_filename,
    parse_web_summ,
)


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
QA_FIXTURES_DIR = os.path.join(FIXTURES_DIR, "qa")


class TestParseRawFilename:
    """Tests for parse_raw_filename function. Uses real S3-key style paths (bucket stripped)."""

    def test_10x_gex_full_path(self):
        """Parse 10x GEX path from fixture-style path."""
        path = "proj/order/Br1_A5/raw/439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT.csv"
        result = parse_raw_filename(path, "10x")
        assert result == ("439047", "Br1_A5", "GEX", "Z0273", "CTGCATGTTGCTGAGAT")

    def test_10x_gex_r1_fastq(self):
        """GEX R1 FASTQ file from a 10x run."""
        path = (
            "lange-human-embryogenesis/AN00028026/Br1_A5/raw/"
            "439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT_S1_L001_R1_001.fastq.gz"
        )
        assert parse_raw_filename(path, "10x") == (
            "439047",
            "Br1_A5",
            "GEX",
            "Z0273",
            "CTGCATGTTGCTGAGAT",
        )

    def test_10x_cri_full_path(self):
        """Parse 10x CRI path."""
        path = "proj/order/Treg_L01/raw/438523-Treg_L01_CRI-Z0012-CTGCCATAGCACGAT.csv"
        result = parse_raw_filename(path, "10x")
        assert result == ("438523", "Treg_L01", "CRI", "Z0012", "CTGCCATAGCACGAT")

    def test_10x_cri_r2_fastq(self):
        """CRI R2 FASTQ file from a 10x run."""
        path = (
            "marson-macrophages-tregs-pilot/AN00027127/Treg_L01/raw/"
            "438523-Treg_L01_CRI-Z0012-CTGCCATAGCACGAT_S1_L001_R2_001.fastq.gz"
        )
        assert parse_raw_filename(path, "10x") == (
            "438523",
            "Treg_L01",
            "CRI",
            "Z0012",
            "CTGCCATAGCACGAT",
        )

    def test_10x_gex_r2_json(self):
        """GEX R2 JSON metadata file (group name contains underscore)."""
        path = (
            "lange-human-embryogenesis/AN00028026/Br1_A5/raw/"
            "439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT_S1_L001_R2_001.json"
        )
        assert parse_raw_filename(path, "10x") == (
            "439047",
            "Br1_A5",
            "GEX",
            "Z0273",
            "CTGCATGTTGCTGAGAT",
        )

    def test_10x_gex_unmatched_cram(self):
        """GEX unmatched CRAM — suffix after barcode is _unmatched.cram."""
        path = (
            "marson-macrophages-tregs-pilot/AN00027127/Treg_L01/raw/"
            "438523-Treg_L01_GEX-Z0011-CACGCACTGCCAGAT_unmatched.cram"
        )
        assert parse_raw_filename(path, "10x") == (
            "438523",
            "Treg_L01",
            "GEX",
            "Z0011",
            "CACGCACTGCCAGAT",
        )

    def test_10x_atac_path(self):
        """Parse 10x ATAC path from multiome fixture."""
        path = "proj/order/CH13/raw/439048-CH13_ATAC-Z0050-CACATGGCAGCACAGAT_S1_L001_R1_001.csv"
        result = parse_raw_filename(path, "10x")
        assert result == ("439048", "CH13", "ATAC", "Z0050", "CACATGGCAGCACAGAT")

    def test_10x_atac_i2_fastq(self):
        """ATAC I2 FASTQ — assay ATAC, group has no underscore."""
        path = (
            "ucsf-killifish-atlas/NVUS2024101701-20/CH13/raw/"
            "439048-CH13_ATAC-Z0050-CACATGGCAGCACAGAT_S1_L001_I2_001.fastq.gz"
        )
        assert parse_raw_filename(path, "10x") == (
            "439048",
            "CH13",
            "ATAC",
            "Z0050",
            "CACATGGCAGCACAGAT",
        )

    def test_10x_gex_novogene_style_path(self):
        """Parse 10x GEX path (Novogene style from multiome_raw)."""
        path = "czi-novogene/ucsf-killifish-atlas/NVUS2024101701-20/CH13/raw/438586-CH13_GEX-Z0005-CATGTATCCTCTGAT.csv"
        result = parse_raw_filename(path, "10x")
        assert result == ("438586", "CH13", "GEX", "Z0005", "CATGTATCCTCTGAT")

    def test_sci_plex_gex_hash_oligo_cram(self):
        """sci_plex: GEX_hash_oligo CRAM in a run subdirectory."""
        path = (
            "hamazaki-seahub-bcp/NVUS2024101701-09/R097/raw/436012/"
            "436012-R097C_GEX_hash_oligo-Z0002-CATGTGCAGCCATCGAT.cram"
        )
        assert parse_raw_filename(path, "sci_plex") == (
            "436012",
            "R097C",
            "GEX_hash_oligo",
            "Z0002",
            "CATGTGCAGCCATCGAT",
        )

    def test_sci_plex_trimmer_stats_dash_in_suffix(self):
        """trimmer-stats.csv has a dash in suffix; barcode extraction still correct."""
        path = (
            "hamazaki-seahub-bcp/NVUS2024101701-09/R097/raw/436012/"
            "436012-R097C_GEX_hash_oligo-Z0046-CTCTCGCATGCAATGAT_trimmer-stats.csv"
        )
        assert parse_raw_filename(path, "sci_plex") == (
            "436012",
            "R097C",
            "GEX_hash_oligo",
            "Z0046",
            "CTCTCGCATGCAATGAT",
        )

    def test_scale_gex_path(self):
        """Parse scale GEX path; group from path, barcode is None."""
        path = "proj/order/EXP1/raw/run-EXP1_GEX-UG01.cram"
        result = parse_raw_filename(path, "scale")
        assert result == ("run", "EXP1", "GEX", "run-EXP1", None)

    def test_scale_gex_cram_real_path(self):
        """Scale GEX CRAM — group from path, assay from regex, ug from filename."""
        path = (
            "trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/"
            "426971-RNA3-098C_GEX_QSR-7_10C.cram"
        )
        run, group, assay, ug, barcode = parse_raw_filename(path, "scale")
        assert run == "426971"
        assert group == "RNA3_098"
        assert assay == "GEX"
        assert ug == "QSR-7"
        assert barcode is None

    def test_scale_hash_oligo_cram_real_path(self):
        """Scale hash_oligo CRAM — assay detected by hash_oligo regex."""
        path = (
            "trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/"
            "426971-RNA3-098C_hash_oligo_QSR-7-SCALEPLEX_1E.cram"
        )
        run, group, assay, ug, barcode = parse_raw_filename(path, "scale")
        assert run == "426971"
        assert group == "RNA3_098"
        assert assay == "hash_oligo"
        assert ug == "QSR-7-SCALEPLEX"
        assert barcode is None

    def test_scale_group_comes_from_path_not_filename(self):
        """Group for scale comes from path[2], independent of filename."""
        path = (
            "trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/"
            "426971-RNA3-098C_GEX_QSR-7_10C.cram"
        )
        _, group, _, _, _ = parse_raw_filename(path, "scale")
        assert group == "RNA3_098"

    def test_scale_hash_oligo_path(self):
        """Parse scale hash_oligo path."""
        path = "proj/order/EXP1/raw/run-EXP1_hash_oligo-UG02.cram"
        result = parse_raw_filename(path, "scale")
        assert result[1] == "EXP1"
        assert result[2] == "hash_oligo"
        assert result[4] is None

    def test_scale_gex_hash_oligo_path(self):
        """Parse scale GEX_hash_oligo path."""
        path = "proj/order/EXP1/raw/run-EXP1_GEX_hash_oligo-UG03.cram"
        result = parse_raw_filename(path, "scale")
        assert result[1] == "EXP1"
        assert result[2] == "GEX_hash_oligo"
        assert result[4] is None

    def test_short_path_returns_none(self):
        """Paths with fewer than 3 hyphen-split parts return None."""
        assert parse_raw_filename("a-b", "10x") is None
        assert parse_raw_filename("onlyonepart", "10x") is None
        assert parse_raw_filename("one-two", "10x") is None

    def test_returns_none_for_filename_with_too_few_dashes(self):
        """Filename with fewer than 3 dash-separated parts returns None."""
        path = "proj/order/GROUP/raw/filename-only.fastq.gz"
        assert parse_raw_filename(path, "10x") is None

    def test_returns_none_for_flat_filename_no_dashes(self):
        """Bare filename with no dashes returns None."""
        path = "proj/order/GROUP/raw/plainfile.csv"
        assert parse_raw_filename(path, "10x") is None

    def test_unknown_assay_falls_back_to_last_underscore_token(self):
        """When assay doesn't match valid_assays, fallback to last _ token of group_assay."""
        path = "proj/order/GROUP1/raw/438523-GROUP1_UNKNOWN-Z0012-CTGCC.csv"
        run, group, assay, ug, barcode = parse_raw_filename(path, "10x")
        assert run == "438523"
        assert group == "GROUP1"
        assert assay == "UNKNOWN"
        assert ug == "Z0012"
        assert barcode == "CTGCC"

    def test_uses_filename_only(self):
        """Only the last path component is used for parsing."""
        path = "any/prefix/439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT.csv"
        result = parse_raw_filename(path, "10x")
        assert result == ("439047", "Br1_A5", "GEX", "Z0273", "CTGCATGTTGCTGAGAT")


class TestExtractRunIdFromTrimmerFilename:
    """Tests for extract_run_id_from_trimmer_filename."""

    def test_valid_6_digit_run_id_examples(self):
        # Real-style examples from SOP / docs
        assert (
            extract_run_id_from_trimmer_filename(
                "434902-pilot_preandpostinj_tech_rep2_CRI-Z0028-CAGACTTGCTGCGAT_trimmer-failure_codes.csv"
            )
            == "434902"
        )
        assert (
            extract_run_id_from_trimmer_filename(
                "436073-R100A_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT_trimmer-failure_codes.csv"
            )
            == "436073"
        )
        assert (
            extract_run_id_from_trimmer_filename(
                "438761-t_fb_GEX-Z0003-CATCACACATGAATGAT_trimmer-failure_codes.csv"
            )
            == "438761"
        )

    def test_valid_8_digit_run_id_example(self):
        # 8-digit wafer id
        assert (
            extract_run_id_from_trimmer_filename(
                "43434720-Mac_L01_CRI-Z0002-CATGTGCAGCCATCGAT_trimmer-failure_codes.csv"
            )
            == "43434720"
        )

    def test_non_numeric_prefix_returns_none(self):
        name = "RUNX1-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT_trimmer-failure_codes.csv"
        assert extract_run_id_from_trimmer_filename(name) is None

    def test_too_short_or_too_long_prefix_returns_none(self):
        # 4 digits: too short
        assert (
            extract_run_id_from_trimmer_filename(
                "1234-Br1_A5_GEX-Z0273-CTGC_trimmer-failure_codes.csv"
            )
            is None
        )
        # 9 digits: too long
        assert (
            extract_run_id_from_trimmer_filename(
                "123456789-Br1_A5_GEX-Z0273-CTGC_trimmer-failure_codes.csv"
            )
            is None
        )


class TestGrabTrimmerStats:
    """Tests for grab_trimmer_stats."""

    def test_aggregates_rsq_and_trimmer_fail_percentages(self, tmp_path):
        # Use a small synthetic CSV fixture
        src = os.path.join(QA_FIXTURES_DIR, "trimmer_failure_codes_small.csv")
        dst = tmp_path / "trimmer_failure_codes_small.csv"
        with open(src, "r") as f_in, open(dst, "w") as f_out:
            f_out.write(f_in.read())

        trimmer_failure_stats: dict = {}
        grab_trimmer_stats(trimmer_failure_stats, "exp1", str(dst))

        assert "exp1" in trimmer_failure_stats
        stats = trimmer_failure_stats["exp1"]
        # rsq row: 100 failed of 1000 total => 10%
        assert stats["rsq"] == [10.0]
        # other reasons: 50 + 50 = 100 failed of 1000 total => 10%
        assert stats["trimmer_fail"] == [10.0]


class TestParseMetSumm:
    """
    Fixture files:
      metrics_summary_count.csv — single-row format from cellranger count (GEX-only)
      metrics_summary.csv       — multi-row format from cellranger multi (GEX + CRISPR)
    """

    def test_single_row_legacy_format(self):
        """Legacy metrics_summary with single row (count pipeline)."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary_count.csv")
        report = parse_met_summ(path)
        assert report == {"GEX_reads": 269_666_389}

    def test_count_format_no_cri_reads(self):
        """Single-row count CSV: no CRI_reads key present."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary_count.csv")
        report = parse_met_summ(path)
        assert "CRI_reads" not in report

    def test_count_format_comma_formatted_reads_parsed_correctly(self):
        """'Number of Reads' with commas (e.g. '269,666,389') is parsed as int."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary_count.csv")
        report = parse_met_summ(path)
        assert isinstance(report["GEX_reads"], int)

    def test_multi_row_with_gex_and_cri(self):
        """Multi-row metrics_summary with Gene Expression and CRISPR Guide Capture."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary.csv")
        report = parse_met_summ(path)
        assert "GEX_reads" in report
        assert "CRI_reads" in report
        assert isinstance(report["GEX_reads"], int)
        assert isinstance(report["CRI_reads"], int)
        assert report["GEX_reads"] > 0
        assert report["CRI_reads"] > 0

    def test_multi_format_gex_reads_value(self):
        """GEX_reads is the sum across all Gene Expression Fastq ID rows."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary.csv")
        report = parse_met_summ(path)
        assert report["GEX_reads"] == 13_397_292_766

    def test_multi_format_cri_reads_value(self):
        """CRI_reads is the sum across all CRISPR Guide Capture Fastq ID rows."""
        path = os.path.join(QA_FIXTURES_DIR, "metrics_summary.csv")
        report = parse_met_summ(path)
        assert report["CRI_reads"] == 3_157_567_476


class TestParseWebSumm:
    """
    Fixture files:
      web_summary.html    — cellranger 9.0.1 count, 3' chemistry, no CRISPR
      web_summary-2.html  — cellranger 9.0.1 multi, 5' chemistry, CRISPR + CellAnnotate
      web_summary-10.html — cellranger 10.0.0 multi, Flex, CRISPR + multiplexing
      web_summary-1.html  — Multiome (ARC) format — NOT supported (raises KeyError)
    """

    # --- cellranger 9.0.1 count (web_summary.html) ---

    def test_cr9_count_software_version(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["software"] == "cellranger-9.0.1"

    def test_cr9_count_subcommand(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["sub"] == "count"

    def test_cr9_count_chemistry_mapped_to_shortcode(self):
        """'Single Cell 3' v4 (polyA)' maps to '3p'."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["chem"] == "3p"

    def test_cr9_count_transcriptome(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["Transcriptome"] == "Coturnix_japonica-2.1"

    def test_cr9_count_ref_equals_transcriptome(self):
        """ref is set from gex_tab transcriptome value."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        report = parse_web_summ(path)
        assert report["ref"] == report["Transcriptome"]

    def test_cr9_count_include_introns(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["incl_int"] == "true"

    def test_cr9_count_no_extra_modalities(self):
        """Count run without CRISPR/Antibody: extra list is empty."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert parse_web_summ(path)["extra"] == []

    def test_cr9_count_gex_alerts_is_list(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary.html")
        assert isinstance(parse_web_summ(path)["gex_alerts"], list)

    # --- cellranger 9.0.1 multi + CRISPR + CellAnnotate (web_summary-2.html) ---

    def test_cr9_multi_crispr_software_version(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["software"] == "cellranger-9.0.1"

    def test_cr9_multi_crispr_subcommand(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["sub"] == "multi"

    def test_cr9_multi_crispr_chemistry(self):
        """'Single Cell 5' R2-only v3' maps to '5p'."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["chem"] == "5p"

    def test_cr9_multi_crispr_extra_contains_crispr(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert "CRISPR" in parse_web_summ(path)["extra"]

    def test_cr9_multi_cellannotate_in_extra(self):
        """skip-cell-annotation=false in experimental_design adds CellAnnotate to extra."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert "CellAnnotate" in parse_web_summ(path)["extra"]

    def test_cr9_multi_crispr_min_crispr_umi(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["min-crispr-umi"] == "3"

    def test_cr9_multi_crispr_create_bam(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["create-bam"] == "true"

    def test_cr9_multi_crispr_include_introns(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert parse_web_summ(path)["incl_int"] == "true"

    def test_cr9_multi_crispr_alerts_are_lists(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        report = parse_web_summ(path)
        assert isinstance(report["gex_alerts"], list)
        assert isinstance(report["crispr_alerts"], list)

    def test_cr9_multi_no_multiplex_flag(self):
        """No [samples] in experimental_design — multiplex key absent."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-2.html")
        assert "multiplex" not in parse_web_summ(path)

    # --- cellranger 10.0.0 flex + CRISPR + multiplexing (web_summary-10.html) ---

    def test_cr10_software_version(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert parse_web_summ(path)["software"] == "cellranger-10.0.0"

    def test_cr10_flex_chemistry(self):
        """'Flex Gene Expression' maps to 'flex'."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert parse_web_summ(path)["chem"] == "flex"

    def test_cr10_flex_no_incl_int(self):
        """Flex chemistry skips include-introns: incl_int key absent."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert "incl_int" not in parse_web_summ(path)

    def test_cr10_flex_probe_set_name(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert (
            parse_web_summ(path)["Probe Set Name"]
            == "Chromium Human Transcriptome Probe Set v1.1.0"
        )

    def test_cr10_flex_crispr_in_extra(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert "CRISPR" in parse_web_summ(path)["extra"]

    def test_cr10_flex_multiplex_flag(self):
        """[samples] in experimental_design sets multiplex=True."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert parse_web_summ(path).get("multiplex") is True

    def test_cr10_flex_crispr_alerts_present(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        report = parse_web_summ(path)
        assert isinstance(report["crispr_alerts"], list)
        assert len(report["crispr_alerts"]) > 0

    def test_cr10_flex_create_bam_false(self):
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        assert parse_web_summ(path)["create-bam"] == "false"

    def test_web_summary_10_returns_expected_keys(self):
        """parse_web_summ on web_summary-10.html returns expected structure."""
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-10.html")
        report = parse_web_summ(path)
        assert "ref" in report
        assert isinstance(report["extra"], list)
        assert "gex_alerts" in report

    # --- unsupported format ---

    def test_multiome_arc_format_raises_error(self):
        """
        Multiome (ARC) web_summary uses a different JSON schema.
        parse_web_summ raises KeyError (expects 'summary' or 'library' keys).
        Documents known limitation for refactoring.
        """
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-1.html")
        with pytest.raises(KeyError):
            parse_web_summ(path)
