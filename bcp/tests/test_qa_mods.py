"""
Tests for qa_mods parsing functions: parse_raw_filename, parse_met_summ, parse_web_summ.
Fixture files in tests/fixtures/qa/ are real outputs from CZI sequencing runs
and 10x Genomics CellRanger pipelines.
"""

import os

import pytest

from qa_mods import (
    extract_read_indicator,
    extract_run_id_from_merged_trimmer_path,
    extract_run_id_from_trimmer_filename,
    grab_merged_trimmer_q30,
    grab_merged_trimmer_stats,
    grab_trimmer_stats,
    is_order_level_processed_folder,
    is_valid_cellranger_run_dir_name,
    make_read_partner,
    normalize_raw_assay,
    parse_met_summ,
    parse_raw_filename,
    parse_scale_samples_csv,
    parse_scale_workflow_info,
    parse_web_summ,
    resolve_qa_run_context,
)


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
QA_FIXTURES_DIR = os.path.join(FIXTURES_DIR, "qa")


class TestIsOrderLevelProcessedFolder:
    def test_processed_under_order(self):
        o = "ny-biohub-califano/NVUS2024101701-17/"
        g = "ny-biohub-califano/NVUS2024101701-17/processed/"
        assert is_order_level_processed_folder(o, g) is True

    def test_raw_under_order_not_special_cased(self):
        assert is_order_level_processed_folder("a/b/", "a/b/raw/") is False

    def test_sample_group_not_skipped(self):
        assert (
            is_order_level_processed_folder("ny/NVUS1/", "ny/NVUS1/MS116A_MS116AF/")
            is False
        )

    def test_case_insensitive_processed(self):
        assert is_order_level_processed_folder("o/", "o/PROCESSED/") is True


class TestIsValidCellrangerRunDirName:
    def test_iso_date_only(self):
        assert is_valid_cellranger_run_dir_name("Run_2025-01-10") is True

    def test_iso_date_2026_with_suffix(self):
        assert is_valid_cellranger_run_dir_name("Run_2026-02-28_biohub") is True

    def test_underscore_separated_date(self):
        assert is_valid_cellranger_run_dir_name("Run_2025_12_31") is True

    def test_rejects_wrong_month(self):
        assert is_valid_cellranger_run_dir_name("Run_2025-13-01") is False

    def test_rejects_no_run_prefix(self):
        assert is_valid_cellranger_run_dir_name("2025-01-10") is False

    def test_rejects_path_with_slash(self):
        assert is_valid_cellranger_run_dir_name("Run_2025-01-10/extra") is False


class TestNormalizeRawAssay:
    def test_10x_case_insensitive(self):
        assert normalize_raw_assay("10X") == "10x"

    def test_10x_viral_orf_case(self):
        assert normalize_raw_assay("10x_viral_orf") == "10x_viral_ORF"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="raw_assay"):
            normalize_raw_assay("")
        with pytest.raises(ValueError, match="raw_assay"):
            normalize_raw_assay("   ")

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="not recognized"):
            normalize_raw_assay("unknown_assay")


class TestResolveQaRunContext:
    def test_s3_from_full_uri(self):
        ctx = resolve_qa_run_context(
            data_source="s3",
            raw_assay="10x",
            s3_path="s3://czi-novogene/myproj/NVUS2024101701-01/",
        )
        assert ctx.bucket == "czi-novogene"
        assert ctx.provider == "novogene"
        assert ctx.proj == "myproj"
        assert ctx.order == "NVUS2024101701-01"
        assert ctx.output_label == "NVUS2024101701-01"
        assert ctx.listing_prefix == "myproj/NVUS2024101701-01/"

    def test_s3_from_components(self):
        ctx = resolve_qa_run_context(
            data_source="s3",
            raw_assay="scale",
            provider="novogene",
            proj="p",
            order="o1",
        )
        assert ctx.bucket == "czi-novogene"
        assert ctx.output_label == "o1"

    def test_s3_run_label_override_output(self):
        ctx = resolve_qa_run_context(
            data_source="s3",
            raw_assay="10x",
            s3_path="s3://czi-novogene/myproj/NVUS2024101701-01/",
            run_label="my_run",
        )
        assert ctx.output_label == "my_run"

    def test_manifest_requires_run_label(self):
        path = os.path.join(FIXTURES_DIR, "test_manifest.tsv")
        with pytest.raises(ValueError, match="run_label"):
            resolve_qa_run_context(
                data_source="manifest",
                raw_assay="10x",
                manifest_path=path,
                manifest_delimiter="\t",
                manifest_s3_column=0,
                manifest_has_header=False,
                run_label="",
            )

    def test_manifest_resolves_bucket(self):
        path = os.path.join(FIXTURES_DIR, "test_manifest.tsv")
        ctx = resolve_qa_run_context(
            data_source="manifest",
            raw_assay="sci_plex",
            manifest_path=path,
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=False,
            run_label="NVUS2024101701-test",
        )
        assert ctx.bucket == "czi-novogene"
        assert ctx.output_label == "NVUS2024101701-test"
        assert ctx.data_source == "manifest"


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

    def test_10x_hyphenated_group_gex_csv(self):
        """10x parser accepts hyphens inside group IDs."""
        path = (
            "wang-tetrapod-atlas/NVUS2024101701-43/fbf_1-1/raw/"
            "440261-fbf_1-1_GEX-Z0052-CTATGCCACAGCATGAT.csv"
        )
        assert parse_raw_filename(path, "10x") == (
            "440261",
            "fbf_1-1",
            "GEX",
            "Z0052",
            "CTATGCCACAGCATGAT",
        )

    def test_10x_hyphenated_group_gex_variants(self):
        """Hyphenated-group parsing stays stable across common 10x suffix variants."""
        base = "440261-fbf_1-1_GEX-Z0052-CTATGCCACAGCATGAT"
        suffixes = [
            ".json",
            "_S1_L001_R1_001.csv",
            "_S1_L001_R1_001.fastq.gz",
            "_S1_L001_R2_001.fastq.gz",
            ".scRNA.applicationQC.h5",
            ".scRNA.applicationQC.html",
            "_Log.final.out",
            "_Log.out",
            "_Log.progress.out",
            "_ReadsPerGene.out.tab",
            "_unmatched.cram",
        ]
        for suffix in suffixes:
            path = f"proj/order/fbf_1-1/raw/{base}{suffix}"
            assert parse_raw_filename(path, "10x") == (
                "440261",
                "fbf_1-1",
                "GEX",
                "Z0052",
                "CTATGCCACAGCATGAT",
            )

    def test_10x_hyphenated_group_other_assays(self):
        """Hyphenated-group parsing works for CRI and ATAC assays."""
        cri = "proj/order/fbf_1-1/raw/440261-fbf_1-1_CRI-Z0052-CTATGCCACAGCATGAT.csv"
        atac = (
            "proj/order/fbf_1-1/raw/"
            "440261-fbf_1-1_ATAC-Z0052-CTATGCCACAGCATGAT_S1_L001_R1_001.fastq.gz"
        )
        assert parse_raw_filename(cri, "10x") == (
            "440261",
            "fbf_1-1",
            "CRI",
            "Z0052",
            "CTATGCCACAGCATGAT",
        )
        assert parse_raw_filename(atac, "10x") == (
            "440261",
            "fbf_1-1",
            "ATAC",
            "Z0052",
            "CTATGCCACAGCATGAT",
        )

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


class TestExtractRunIdFromMergedTrimmerPath:
    """Tests for extract_run_id_from_merged_trimmer_path."""

    def test_scale_sci_path_parent_run_id(self):
        """Scale/sci: run_id from path parent (raw/438761/merged_*)."""
        path = "proj/order/RNA3_098/raw/438761/merged_trimmer-failure_codes.csv"
        assert extract_run_id_from_merged_trimmer_path(path) == "438761"
        path2 = "proj/order/RNA3_098/raw/438761/merged_trimmer-stats.csv"
        assert extract_run_id_from_merged_trimmer_path(path2) == "438761"

    def test_scale_sci_filename_in_run_dir(self):
        """Scale/sci: run_id from path parent when filename has run_id prefix."""
        path = "proj/order/R097/raw/436012/436012_merged_trimmer-failure_codes.csv"
        assert extract_run_id_from_merged_trimmer_path(path) == "436012"

    def test_10x_order_level_filename_prefix(self):
        """10x: merged files under order/; run_id from filename prefix."""
        path = "proj/order/438761_merged_trimmer-failure_codes.csv"
        assert extract_run_id_from_merged_trimmer_path(path) == "438761"
        path2 = (
            "czi-novogene/weissman/NVUS2024101701-29/438761_merged_trimmer-stats.csv"
        )
        assert extract_run_id_from_merged_trimmer_path(path2) == "438761"

    def test_invalid_or_non_merged_returns_none(self):
        """Non-merged paths or missing run_id return None."""
        assert (
            extract_run_id_from_merged_trimmer_path("proj/order/raw/file.csv") is None
        )
        assert (
            extract_run_id_from_merged_trimmer_path(
                "proj/order/group/raw/438761/sample_trimmer-failure_codes.csv"
            )
            is None
        )
        # 10x order-level but no prefix (single wafer)
        assert (
            extract_run_id_from_merged_trimmer_path(
                "proj/order/merged_trimmer-failure_codes.csv"
            )
            is None
        )


class TestGrabMergedTrimmerStats:
    """Tests for grab_merged_trimmer_stats."""

    def test_tt_and_rsq_from_fixture(self, tmp_path):
        """Fixture has only TT rows: TT = sum all TT failures; RSQ = all rsq file rows."""
        src = os.path.join(QA_FIXTURES_DIR, "merged_trimmer_failure_tt_small.csv")
        dst = tmp_path / "merged_trimmer_failure.csv"
        with open(src, "r") as f_in, open(dst, "w") as f_out:
            f_out.write(f_in.read())
        result = grab_merged_trimmer_stats(str(dst))
        assert result is not None
        tt_total = 158267504
        tt_failed = 24911820 + 355 + 6
        assert result["tt_total_reads"] == tt_total
        assert result["tt_failed_reads"] == tt_failed
        assert result["tt_pass_reads"] == tt_total - tt_failed
        assert result["tt_fail_pct"] == pytest.approx(
            100.0 * tt_failed / tt_total, rel=1e-5
        )
        assert result["tt_pass_pct"] == pytest.approx(
            100.0 * (tt_total - tt_failed) / tt_total, rel=1e-5
        )
        # RSQ: only one "rsq file" row (TT) in this fixture
        assert result["rsq_total_reads"] == tt_total
        assert result["rsq_failed_reads"] == 24911820
        assert result["rsq_pass_reads"] == tt_total - 24911820
        assert result["rsq_fail_pct"] == pytest.approx(
            100.0 * 24911820 / tt_total, rel=1e-5
        )

    def test_rsq_sums_across_all_read_groups(self, tmp_path):
        """RSQ totals and failures are summed across all read groups with reason 'rsq file'."""
        path = tmp_path / "merged.csv"
        path.write_text(
            "read group,code,format,segment,reason,failed read count,total read count\n"
            "TT,8,trim,preamble,rsq file,100,1000\n"
            "TT,101,trim,insert,sequence was too short,10,1000\n"
            "UGAv3-1000,8,no trimming,insert,rsq file,50,500\n"
        )
        result = grab_merged_trimmer_stats(str(path))
        assert result is not None
        # TT: total from first TT row, failed = all TT rows
        assert result["tt_total_reads"] == 1000
        assert result["tt_failed_reads"] == 110  # 100 + 10
        assert result["tt_pass_reads"] == 890
        assert result["tt_fail_pct"] == pytest.approx(11.0, rel=1e-5)
        # RSQ: all rows with reason "rsq file" across all read groups
        assert result["rsq_total_reads"] == 1500  # 1000 + 500
        assert result["rsq_failed_reads"] == 150  # 100 + 50
        assert result["rsq_pass_reads"] == 1350
        assert result["rsq_fail_pct"] == pytest.approx(100.0 * 150 / 1500, rel=1e-5)

    def test_no_tt_row_returns_none(self, tmp_path):
        """CSV without TT read_group returns None."""
        path = tmp_path / "no_tt.csv"
        path.write_text(
            "read group,reason,failed read count,total read count\n"
            "UGAv3-1000,rsq file,100,1000\n"
        )
        assert grab_merged_trimmer_stats(str(path)) is None


class TestGrabMergedTrimmerQ30:
    """Tests for grab_merged_trimmer_q30."""

    def test_tt_low_quality_bases_from_fixture(self, tmp_path):
        """TT row with low quality bases segment returns Q30 percentage."""
        src = os.path.join(QA_FIXTURES_DIR, "merged_trimmer_stats_tt_small.csv")
        dst = tmp_path / "merged_trimmer_stats.csv"
        with open(src, "r") as f_in, open(dst, "w") as f_out:
            f_out.write(f_in.read())
        result = grab_merged_trimmer_q30(str(dst))
        assert result is not None
        # 3807451195 / (3807451195 + 7211)
        expected = 100.0 * 3807451195 / (3807451195 + 7211)
        assert abs(result - expected) < 0.0001

    def test_no_tt_row_returns_none(self, tmp_path):
        """CSV without TT read_group returns None."""
        path = tmp_path / "no_tt.csv"
        path.write_text(
            "read group,segment label,num matched bases,num failures\n"
            "Z0003,low quality bases,1000,50\n"
        )
        assert grab_merged_trimmer_q30(str(path)) is None

    def test_no_low_quality_segment_returns_none(self, tmp_path):
        """TT row without low quality bases segment returns None."""
        path = tmp_path / "no_segment.csv"
        path.write_text(
            "read group,segment label,num matched bases,num failures\n"
            "TT,insert,1000,50\n"
        )
        assert grab_merged_trimmer_q30(str(path)) is None


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
        Either no ``const data =`` payload (ValueError) or wrong shape (KeyError).
        """
        path = os.path.join(QA_FIXTURES_DIR, "web_summary-1.html")
        with pytest.raises((ValueError, KeyError)):
            parse_web_summ(path)


class TestParseScaleWorkflowInfo:
    """Tests for parse_scale_workflow_info."""

    def test_parses_parameters_and_manifest(self):
        """Happy path: all expected keys are present and correctly extracted."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_workflow_info_good.json")
        result = parse_scale_workflow_info(path)
        assert result["bamOut"] == "true"
        assert result["scalePlex"] == "true"
        assert result["scalePlexAssignmentMethod"] == "fc"
        assert result["workflow_name"] == "ScaleRna"
        assert result["workflow_version"] == "2.1.0"
        assert result["execution_status"] == "OK"

    def test_missing_parameters_raises(self):
        """ValueError raised when Parameters section is absent."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_workflow_info_bad_no_params.json")
        # Create a fixture missing Parameters
        import json

        data = {"Workflow Manifest": {"name": "ScaleRna", "version": "2.1.0"}}
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(data, f)
        try:
            with pytest.raises(ValueError, match="Parameters"):
                parse_scale_workflow_info(path)
        finally:
            os.remove(path)

    def test_genome_value_extracted(self):
        """Genome S3 path is extracted correctly from Parameters."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_workflow_info_good.json")
        result = parse_scale_workflow_info(path)
        assert result["genome"] is not None
        assert "GRCh38" in result["genome"]


class TestParseScaleSamplesCsv:
    """Tests for parse_scale_samples_csv."""

    def test_extracts_columns_and_samples(self):
        """Happy path: columns and samples are correctly extracted."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_samples_good.csv")
        result = parse_scale_samples_csv(path)
        assert "sample" in result["columns"]
        assert "libIndex2" in result["columns"]
        assert len(result["samples"]) > 0
        assert "SAMP-01-0001" in result["samples"]

    def test_sublibraries_parsed_from_libindex2(self):
        """Semicolon-delimited libIndex2 is correctly split into sublibrary IDs."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_samples_good.csv")
        result = parse_scale_samples_csv(path)
        sublibs = result["sublibraries"]["SAMP-01-0001"]
        assert sublibs == [
            "QSR-1",
            "QSR-2",
            "QSR-3",
            "QSR-4",
            "QSR-5",
            "QSR-6",
            "QSR-7",
            "QSR-8",
        ]

    def test_sample_count_matches(self):
        """Number of samples matches the number of rows in the CSV."""
        path = os.path.join(QA_FIXTURES_DIR, "scale_samples_good.csv")
        result = parse_scale_samples_csv(path)
        assert len(result["samples"]) == 4


class TestExtractReadIndicator:
    """Tests for extract_read_indicator — Illumina read indicator from filename tail."""

    def test_standard_r1(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        assert extract_read_indicator(f) == "R1"

    def test_standard_r2(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        assert extract_read_indicator(f) == "R2"

    def test_index_read_i1(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_I1_001.fastq.gz"
        assert extract_read_indicator(f) == "I1"

    def test_index_read_i2(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_I2_001.fastq.gz"
        assert extract_read_indicator(f) == "I2"

    def test_r3_read(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_R3_001.fastq.gz"
        assert extract_read_indicator(f) == "R3"

    def test_r2_in_group_id_returns_tail_r1(self):
        """R2 in group ID is ignored; the tail R1 indicator is returned."""
        f = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R1_001.fastq.gz"
        assert extract_read_indicator(f) == "R1"

    def test_r1_in_group_id_returns_tail_r2(self):
        """R1 in group ID is ignored; the tail R2 indicator is returned."""
        f = "439925-q_hf_R1_GEX-Z0004-CTGTGTAGGCATGAT_S1_L001_R2_001.fastq.gz"
        assert extract_read_indicator(f) == "R2"

    def test_no_indicator_returns_none(self):
        assert extract_read_indicator("439047-G1_GEX-Z0273-BC01.csv") is None

    def test_cram_file_returns_none(self):
        assert extract_read_indicator("426971-RNA3-098C_GEX_QSR-7_10C.cram") is None


class TestMakeReadPartner:
    """Tests for make_read_partner — swap Illumina read indicator at tail only."""

    def test_simple_r1_to_r2(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        expected = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        assert make_read_partner(f, "R1", "R2") == expected

    def test_simple_r2_to_r1(self):
        f = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"
        expected = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        assert make_read_partner(f, "R2", "R1") == expected

    def test_preserves_r2_in_group_id(self):
        """Only the tail indicator is swapped; R2 in the group ID is untouched."""
        f = "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R1_001.fastq.gz"
        result = make_read_partner(f, "R1", "R2")
        assert (
            result
            == "439925-q_pcf_R2_GEX-Z0028-CAGACTTGCTGCGAT_S1_L001_R2_001.fastq.gz"
        )
        assert result.count("_R2_") == 2  # group ID R2 + tail R2

    def test_preserves_r1_in_group_id(self):
        """Only the tail indicator is swapped; R1 in the group ID is untouched."""
        f = "439925-q_hf_R1_GEX-Z0004-CTGTGTAGGCATGAT_S1_L001_R2_001.fastq.gz"
        result = make_read_partner(f, "R2", "R1")
        assert (
            result == "439925-q_hf_R1_GEX-Z0004-CTGTGTAGGCATGAT_S1_L001_R1_001.fastq.gz"
        )
        assert result.count("_R1_") == 2  # group ID R1 + tail R1

    def test_no_indicator_returns_unchanged(self):
        f = "439047-G1_GEX-Z0273-BC01.csv"
        assert make_read_partner(f, "R1", "R2") == f
