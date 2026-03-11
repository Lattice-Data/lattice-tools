"""
Tests for qa_mods helper functions used in QA notebook.

Covers:
  - parse_raw_filename: filename parsing for all assay types
  - parse_met_summ: metrics_summary.csv parsing (count single-row and multi multi-row formats)
  - parse_web_summ: web_summary.html parsing (cellranger 9.0.1 count, 9.0.1 multi, 10.0.0 flex)

Fixture files in tests/fixtures/qa/ are real outputs from CZI sequencing runs
and 10x Genomics CellRanger pipelines.
"""
import os
import pytest
from qa_mods import parse_raw_filename, parse_met_summ, parse_web_summ


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), 'fixtures', 'qa')


# ---------------------------------------------------------------------------
# parse_raw_filename
# ---------------------------------------------------------------------------

class TestParseRawFilename:
    """
    Tests use real S3 keys (bucket prefix stripped, as the function receives them).
    Keys follow the pattern: {proj}/{order}/{group}/raw/[{run_dir}/]{filename}
    """

    # --- 10x assay ---

    def test_10x_cri_r2_fastq(self):
        """CRI R2 FASTQ file from a 10x run."""
        key = (
            'marson-macrophages-tregs-pilot/AN00027127/Treg_L01/raw/'
            '438523-Treg_L01_CRI-Z0012-CTGCCATAGCACGAT_S1_L001_R2_001.fastq.gz'
        )
        assert parse_raw_filename(key, '10x') == (
            '438523', 'Treg_L01', 'CRI', 'Z0012', 'CTGCCATAGCACGAT'
        )

    def test_10x_gex_r1_fastq(self):
        """GEX R1 FASTQ file from a 10x run."""
        key = (
            'lange-human-embryogenesis/AN00028026/Br1_A5/raw/'
            '439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT_S1_L001_R1_001.fastq.gz'
        )
        assert parse_raw_filename(key, '10x') == (
            '439047', 'Br1_A5', 'GEX', 'Z0273', 'CTGCATGTTGCTGAGAT'
        )

    def test_10x_gex_r2_json(self):
        """GEX R2 JSON metadata file from a 10x run (group name contains underscore)."""
        key = (
            'lange-human-embryogenesis/AN00028026/Br1_A5/raw/'
            '439047-Br1_A5_GEX-Z0273-CTGCATGTTGCTGAGAT_S1_L001_R2_001.json'
        )
        assert parse_raw_filename(key, '10x') == (
            '439047', 'Br1_A5', 'GEX', 'Z0273', 'CTGCATGTTGCTGAGAT'
        )

    def test_10x_gex_unmatched_cram(self):
        """GEX unmatched CRAM file — suffix after barcode is '_unmatched.cram'."""
        key = (
            'marson-macrophages-tregs-pilot/AN00027127/Treg_L01/raw/'
            '438523-Treg_L01_GEX-Z0011-CACGCACTGCCAGAT_unmatched.cram'
        )
        assert parse_raw_filename(key, '10x') == (
            '438523', 'Treg_L01', 'GEX', 'Z0011', 'CACGCACTGCCAGAT'
        )

    def test_10x_atac_fastq(self):
        """ATAC I2 FASTQ file — assay is ATAC, group has no underscore."""
        key = (
            'ucsf-killifish-atlas/NVUS2024101701-20/CH13/raw/'
            '439048-CH13_ATAC-Z0050-CACATGGCAGCACAGAT_S1_L001_I2_001.fastq.gz'
        )
        assert parse_raw_filename(key, '10x') == (
            '439048', 'CH13', 'ATAC', 'Z0050', 'CACATGGCAGCACAGAT'
        )

    # --- sci_plex assay ---
    # sci_plex files use an extra run subdirectory under raw/ but filename
    # parsing is identical to 10x (non-scale branch).

    def test_sci_plex_gex_hash_oligo_cram(self):
        """GEX_hash_oligo CRAM in a run subdirectory."""
        key = (
            'hamazaki-seahub-bcp/NVUS2024101701-09/R097/raw/436012/'
            '436012-R097C_GEX_hash_oligo-Z0002-CATGTGCAGCCATCGAT.cram'
        )
        assert parse_raw_filename(key, 'sci_plex') == (
            '436012', 'R097C', 'GEX_hash_oligo', 'Z0002', 'CATGTGCAGCCATCGAT'
        )

    def test_sci_plex_gex_hash_oligo_trimmer_stats(self):
        """
        trimmer-stats.csv has a dash in the suffix, causing extra elements when
        splitting the filename on '-'. Barcode extraction must still be correct.
        """
        key = (
            'hamazaki-seahub-bcp/NVUS2024101701-09/R097/raw/436012/'
            '436012-R097C_GEX_hash_oligo-Z0046-CTCTCGCATGCAATGAT_trimmer-stats.csv'
        )
        assert parse_raw_filename(key, 'sci_plex') == (
            '436012', 'R097C', 'GEX_hash_oligo', 'Z0046', 'CTCTCGCATGCAATGAT'
        )

    # --- scale assay ---
    # Group is extracted from path position [2], not the filename.
    # ug is the second-to-last '_'-delimited token in the filename.
    # barcode is always None.

    def test_scale_gex_cram(self):
        """Scale GEX CRAM — group from path, assay from regex, ug from filename."""
        key = (
            'trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/'
            '426971-RNA3-098C_GEX_QSR-7_10C.cram'
        )
        run, group, assay, ug, barcode = parse_raw_filename(key, 'scale')
        assert run == '426971'
        assert group == 'RNA3_098'
        assert assay == 'GEX'
        assert ug == 'QSR-7'
        assert barcode is None

    def test_scale_hash_oligo_cram(self):
        """Scale hash_oligo CRAM — assay detected by 'hash_oligo' regex, not 'GEX'."""
        key = (
            'trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/'
            '426971-RNA3-098C_hash_oligo_QSR-7-SCALEPLEX_1E.cram'
        )
        run, group, assay, ug, barcode = parse_raw_filename(key, 'scale')
        assert run == '426971'
        assert group == 'RNA3_098'
        assert assay == 'hash_oligo'
        assert ug == 'QSR-7-SCALEPLEX'
        assert barcode is None

    def test_scale_group_comes_from_path_not_filename(self):
        """Group for scale comes from path[2], independent of filename content."""
        key = (
            'trapnell-seahub-bcp/NVUS2024101701-04/RNA3_098/raw/426971/'
            '426971-RNA3-098C_GEX_QSR-7_10C.cram'
        )
        _, group, _, _, _ = parse_raw_filename(key, 'scale')
        assert group == 'RNA3_098'

    # --- edge cases ---

    def test_returns_none_for_filename_with_too_few_dashes(self):
        """Filenames with fewer than 3 dash-separated parts return None."""
        key = 'proj/order/GROUP/raw/filename-only.fastq.gz'
        assert parse_raw_filename(key, '10x') is None

    def test_returns_none_for_flat_filename_no_dashes(self):
        """A bare filename with no dashes at all returns None."""
        key = 'proj/order/GROUP/raw/plainfile.csv'
        assert parse_raw_filename(key, '10x') is None

    def test_unknown_assay_falls_back_to_last_underscore_token(self):
        """
        When the assay portion doesn't match any valid_assays entry, the function
        falls back to the last '_'-delimited token of group_assay as the assay name.
        """
        key = (
            'proj/order/GROUP1/raw/'
            '438523-GROUP1_UNKNOWN-Z0012-CTGCC.csv'
        )
        run, group, assay, ug, barcode = parse_raw_filename(key, '10x')
        assert run == '438523'
        assert group == 'GROUP1'
        assert assay == 'UNKNOWN'
        assert ug == 'Z0012'
        assert barcode == 'CTGCC'


# ---------------------------------------------------------------------------
# parse_met_summ
# ---------------------------------------------------------------------------

class TestParseMetSumm:
    """
    Fixture files:
      metrics_summary_count.csv  — single-row format from cellranger count
                                   (GEX-only, comma-formatted "Number of Reads")
      metrics_summary.csv        — multi-row format from cellranger multi
                                   (GEX + CRISPR Guide Capture rows)
    """

    def test_count_format_returns_gex_reads(self):
        """Single-row count CSV: only GEX_reads is returned."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary_count.csv')
        report = parse_met_summ(f)
        assert report == {'GEX_reads': 269666389}

    def test_count_format_no_cri_reads(self):
        """Single-row count CSV: no CRI_reads key present."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary_count.csv')
        report = parse_met_summ(f)
        assert 'CRI_reads' not in report

    def test_multi_format_gex_and_cri_reads(self):
        """Multi-row multi CSV (GEX + CRISPR): both GEX_reads and CRI_reads returned."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary.csv')
        report = parse_met_summ(f)
        assert 'GEX_reads' in report
        assert 'CRI_reads' in report

    def test_multi_format_gex_reads_value(self):
        """GEX_reads is the sum across all Gene Expression Fastq ID rows."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary.csv')
        report = parse_met_summ(f)
        assert report['GEX_reads'] == 13397292766

    def test_multi_format_cri_reads_value(self):
        """CRI_reads is the sum across all CRISPR Guide Capture Fastq ID rows."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary.csv')
        report = parse_met_summ(f)
        assert report['CRI_reads'] == 3157567476

    def test_count_format_comma_formatted_reads_parsed_correctly(self):
        """'Number of Reads' value with commas (e.g. '269,666,389') is parsed as int."""
        f = os.path.join(FIXTURES_DIR, 'metrics_summary_count.csv')
        report = parse_met_summ(f)
        assert isinstance(report['GEX_reads'], int)


# ---------------------------------------------------------------------------
# parse_web_summ
# ---------------------------------------------------------------------------

class TestParseWebSumm:
    """
    Fixture files in tests/fixtures/qa/:
      web_summary.html    — cellranger 9.0.1, 'count' subcommand, 3' chemistry, no CRISPR
      web_summary-2.html  — cellranger 9.0.1, 'multi' subcommand, 5' chemistry,
                            with CRISPR and CellAnnotate
      web_summary-10.html — cellranger 10.0.0, 'multi' subcommand, Flex chemistry,
                            with CRISPR and [samples] multiplexing
      web_summary-1.html  — Multiome (ARC) format — NOT supported by parse_web_summ
    """

    # --- cellranger 9.0.1 count (web_summary.html) ---

    def test_cr9_count_software_version(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['software'] == 'cellranger-9.0.1'

    def test_cr9_count_subcommand(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['sub'] == 'count'

    def test_cr9_count_chemistry_mapped_to_shortcode(self):
        """'Single Cell 3' v4 (polyA)' maps to '3p'."""
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['chem'] == '3p'

    def test_cr9_count_transcriptome(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['Transcriptome'] == 'Coturnix_japonica-2.1'

    def test_cr9_count_ref_equals_transcriptome(self):
        """ref is overwritten by gex_tab transcriptome value after experimental_design loop."""
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        report = parse_web_summ(f)
        assert report['ref'] == report['Transcriptome']

    def test_cr9_count_include_introns(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['incl_int'] == 'true'

    def test_cr9_count_no_extra_modalities(self):
        """count run without CRISPR/Antibody: extra list is empty."""
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert parse_web_summ(f)['extra'] == []

    def test_cr9_count_gex_alerts_is_list(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary.html')
        assert isinstance(parse_web_summ(f)['gex_alerts'], list)

    # --- cellranger 9.0.1 multi + CRISPR + CellAnnotate (web_summary-2.html) ---

    def test_cr9_multi_crispr_software_version(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['software'] == 'cellranger-9.0.1'

    def test_cr9_multi_crispr_subcommand(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['sub'] == 'multi'

    def test_cr9_multi_crispr_chemistry(self):
        """'Single Cell 5' R2-only v3' maps to '5p'."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['chem'] == '5p'

    def test_cr9_multi_crispr_extra_contains_crispr(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert 'CRISPR' in parse_web_summ(f)['extra']

    def test_cr9_multi_cellannotate_in_extra(self):
        """skip-cell-annotation=false in experimental_design adds CellAnnotate to extra."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert 'CellAnnotate' in parse_web_summ(f)['extra']

    def test_cr9_multi_crispr_min_crispr_umi(self):
        """min-crispr-umi parsed from experimental_design csv."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['min-crispr-umi'] == '3'

    def test_cr9_multi_crispr_create_bam(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['create-bam'] == 'true'

    def test_cr9_multi_crispr_include_introns(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert parse_web_summ(f)['incl_int'] == 'true'

    def test_cr9_multi_crispr_alerts_are_lists(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        report = parse_web_summ(f)
        assert isinstance(report['gex_alerts'], list)
        assert isinstance(report['crispr_alerts'], list)

    def test_cr9_multi_no_multiplex_flag(self):
        """No [samples] section in experimental_design — multiplex key absent."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-2.html')
        assert 'multiplex' not in parse_web_summ(f)

    # --- cellranger 10.0.0 flex + CRISPR + multiplexing (web_summary-10.html) ---

    def test_cr10_software_version(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert parse_web_summ(f)['software'] == 'cellranger-10.0.0'

    def test_cr10_flex_chemistry(self):
        """'Flex Gene Expression' maps to 'flex'."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert parse_web_summ(f)['chem'] == 'flex'

    def test_cr10_flex_no_incl_int(self):
        """Flex chemistry skips include-introns: incl_int key absent."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert 'incl_int' not in parse_web_summ(f)

    def test_cr10_flex_probe_set_name(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert parse_web_summ(f)['Probe Set Name'] == 'Chromium Human Transcriptome Probe Set v1.1.0'

    def test_cr10_flex_crispr_in_extra(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert 'CRISPR' in parse_web_summ(f)['extra']

    def test_cr10_flex_multiplex_flag(self):
        """[samples] section in experimental_design sets multiplex=True."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert parse_web_summ(f).get('multiplex') is True

    def test_cr10_flex_crispr_alerts_present(self):
        """cr10 CRISPR run has non-empty crispr_alerts list."""
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        report = parse_web_summ(f)
        assert isinstance(report['crispr_alerts'], list)
        assert len(report['crispr_alerts']) > 0

    def test_cr10_flex_create_bam_false(self):
        f = os.path.join(FIXTURES_DIR, 'web_summary-10.html')
        assert parse_web_summ(f)['create-bam'] == 'false'

    # --- unsupported format ---

    def test_multiome_arc_format_raises_error(self):
        """
        Multiome (ARC) web_summary.html uses a different JSON schema.
        parse_web_summ raises KeyError because it expects 'summary' or 'library' keys.
        This documents a known limitation to guard against during refactoring.
        """
        f = os.path.join(FIXTURES_DIR, 'web_summary-1.html')
        with pytest.raises(KeyError):
            parse_web_summ(f)
