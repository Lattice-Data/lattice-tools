"""
CLI tests for mapping_validation.main.

These tests exercise the argparse-based entrypoint with small temporary
mapping files and assert on exit codes and high-level verdict output.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

import mapping_validation


def _write_temp_mapping(tmp_path: Path, contents: str) -> Path:
    path = tmp_path / "mapping_cli.csv"
    path.write_text(contents)
    return path


_TENX_CRAM_CLI_PREFIX = (
    "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
    "442356-LeS188W1_GEX-Z0083-CAGTGTATTGCTGAT"
)


def _tenx_cram_cli_mapping(
    *,
    omit_suffixes: set[str] | None = None,
    local_overrides: dict[str, str] | None = None,
    extra_rows: list[str] | None = None,
) -> str:
    """Build CLI mapping text for a complete 10x_cram bundle."""
    omit_suffixes = omit_suffixes or set()
    local_overrides = local_overrides or {}
    lines: list[str] = []
    for suffix in (
        ".cram",
        ".csv",
        ".json",
        "_extract_stats.h5",
        "_SNVQ.metric",
        "_FlowQ.metric",
        "_trimmer-stats.csv",
        "_trimmer-failure_codes.csv",
    ):
        if suffix in omit_suffixes:
            continue
        local = local_overrides.get(suffix, f"/provider/export/sample{suffix}")
        lines.append(f"{_TENX_CRAM_CLI_PREFIX}{suffix},{local}")

    metadata_prefix = "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
    for filename in (
        "442356_LibraryInfo.xml",
        "442356_UploadCompleted.json",
        "run_SecondaryAnalysis.txt",
        "run_VariantCalling.txt",
        "442356_merged_trimmer-stats.csv",
        "442356_merged_trimmer-failure_codes.csv",
    ):
        lines.append(f"{metadata_prefix}{filename},/provider/export/{filename}")

    lines.extend(extra_rows or [])
    return "\n".join(lines) + "\n"


def _run_10x_cram_cli(tmp_path: Path, monkeypatch, mapping_text: str) -> None:
    """Execute the CLI in 10x_cram mode for a temporary mapping."""
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)
    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x_cram",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    mapping_validation.main()


_TENX_FASTQ_CLI_PREFIX = (
    "s3://czi-novogene/project-scaling-alpha/NVUS0000000000-29/CD4i_R1L01/raw/"
    "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT"
)


def _tenx_fastq_cli_mapping(
    *,
    s3_prefix: str = _TENX_FASTQ_CLI_PREFIX,
    assay_token: str = "GEX",
    omit_per_prefix: set[str] | None = None,
    omit_per_tail: set[tuple[str, str]] | None = None,
    extra_rows: list[str] | None = None,
) -> str:
    """Build mapping text for a SOP-complete 10x raw FASTQ bundle (R1 + R2)."""
    omit_per_prefix = omit_per_prefix or set()
    omit_per_tail = omit_per_tail or set()
    lines: list[str] = []
    for sfx in (
        ".csv",
        ".json",
        "_trimmer-stats.csv",
        "_trimmer-failure_codes.csv",
        "_unmatched.cram",
        "_unmatched.csv",
        "_unmatched.json",
    ):
        if sfx in omit_per_prefix:
            continue
        lines.append(f"{s3_prefix}{sfx},/local/{assay_token}{sfx}")
    for tail in ("S1_L001_R1_001", "S1_L001_R2_001"):
        for sfx in (".fastq.gz", ".csv", ".json", "_sample.fastq.gz"):
            if (tail, sfx) in omit_per_tail:
                continue
            lines.append(f"{s3_prefix}_{tail}{sfx},/local/{assay_token}_{tail}{sfx}")
    lines.extend(extra_rows or [])
    return "\n".join(lines) + "\n"


def test_cli_passes_for_valid_10x_mapping(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should exit 0 for a SOP-complete 10x Novogene FASTQ mapping."""
    mapping_text = _tenx_fastq_cli_mapping()
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "10x raw FASTQ completeness" in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_fails_on_duplicate_mappings(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should fail (exit 1) when duplicate mappings are present."""
    # The S3 paths here are not 10x SOP-compliant on purpose – we only
    # care that uniqueness fails and the CLI reports a non-zero status.
    mapping_text = "s3://bucket/a,/local/a\ns3://bucket/a,/local/b\n"
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "VERDICT: FAIL" in captured.out


_SCALE_CLI_S3_PREFIX = (
    "s3://czi-novogene/lab-seahub-alpha/NVUS0000000000-26/CHEM13-R096/raw/441969/"
    "441969-R096A_GEX_QSR-1"
)
_SCALE_CLI_LOCAL_PREFIX = (
    "/local_root/data1/V129/441969-20260220_2053/441969-QSR1_QSR-1/441969-QSR1_QSR-1"
)


def _scale_cli_mapping(
    *,
    omit_per_well: set[str] | None = None,
    omit_per_ug: set[str] | None = None,
    extra_rows: list[str] | None = None,
) -> str:
    """Build mapping text for a SOP-complete Scale per-well + per-UG bundle."""
    omit_per_well = omit_per_well or set()
    omit_per_ug = omit_per_ug or set()
    lines: list[str] = []
    for ext in (".cram", ".csv", ".json"):
        if ext in omit_per_well:
            continue
        lines.append(
            f"{_SCALE_CLI_S3_PREFIX}_7A{ext},{_SCALE_CLI_LOCAL_PREFIX}_7A{ext}"
        )
    for sfx in (
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
        "_unmatched.cram",
        "_unmatched.csv",
        "_unmatched.json",
    ):
        if sfx in omit_per_ug:
            continue
        lines.append(f"{_SCALE_CLI_S3_PREFIX}{sfx},{_SCALE_CLI_LOCAL_PREFIX}{sfx}")
    lines.extend(extra_rows or [])
    return "\n".join(lines) + "\n"


def test_cli_scale_raw_with_sif_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a SOP-complete Scale raw mapping with SIF."""
    mapping_path = _write_temp_mapping(tmp_path, _scale_cli_mapping())

    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM13-R096,CHEM13-R096_GEX,CTATGCACA,lab-seahub-alpha,"
        "CHEM13-R096,R096A,GEX\n"
    )
    sif_path = tmp_path / "scale_sif.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "scale",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "scale raw completeness" in captured.out
    assert "VERDICT: PASS" in captured.out


_SCI_CLI_S3_PREFIX = (
    "s3://czi-novogene/lab-seahub-beta/NVUS0000000000-32/CHEM3-R100/raw/441389/"
    "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT"
)
_SCI_CLI_LOCAL_PREFIX = (
    "/local_root/newsftp/S3/ultima/CR0-789/441389-20260224_2053/"
    "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT/"
    "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT"
)


def _sci_cli_mapping(
    *, omit: set[str] | None = None, extra_rows: list[str] | None = None
) -> str:
    """Build mapping text for a SOP-complete sci raw bundle."""
    omit = omit or set()
    suffixes = (
        ".cram",
        ".csv",
        ".json",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
        "_trimmer-stats_FlowQ.metric",
        "_trimmer-stats_SNVQ.metric",
    )
    lines: list[str] = []
    for sfx in suffixes:
        if sfx in omit:
            continue
        lines.append(f"{_SCI_CLI_S3_PREFIX}{sfx},{_SCI_CLI_LOCAL_PREFIX}{sfx}")
    lines.extend(extra_rows or [])
    return "\n".join(lines) + "\n"


def test_cli_sci_raw_with_sif_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a SOP-complete sci raw mapping with SIF."""
    mapping_path = _write_temp_mapping(tmp_path, _sci_cli_mapping())

    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM3-R100,R100E,Z0028,lab-seahub-beta,CHEM3-R100,R100E,GEX_hash_oligo\n"
    )
    sif_path = tmp_path / "sci_sif.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "sci",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "sci raw completeness" in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_processed_with_sif_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a minimal valid 10x processed mapping with SIF."""
    mapping_text = (
        "s3://czi-novogene/project-embryo-alpha/NVUS0000000000-19/"
        "e10_rep1_t13/processed/cellranger/Run_2000-03-10/outs/"
        "filtered_feature_bc_matrix.h5,"
        "/local/user_001/Ultima/projects_202602/X000SC00000000-Z00-F000_"
        "GRCm39-vM37_10xcellranger_v9.0/Data_process/sampleMatrix/e10_rep1_t13/outs/"
        "filtered_feature_bc_matrix.h5\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    sif_text = "Library name,Group Identifier,Assay Type\nLIB1,e10_rep1_t13,GEX\n"
    sif_path = tmp_path / "tenx_proc_sif.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "processed",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_processed_reports_normalized_groupid_warning(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """CLI should explicitly report '-'/'_' normalized GroupID warning counts."""
    mapping_text = (
        "s3://czi-novogene/test-project/NVUS0000000000-43/"
        "fbm_1-2/processed/cellranger/Run_2001-03-31/outs/"
        "filtered_feature_bc_matrix.h5,"
        "/local/user_002/projects_3_2026/X000SC00000000-Z00-F000/"
        "Data_process/sampleMatrix/fbm_1_2/outs/filtered_feature_bc_matrix.h5\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)
    sif_text = "Library name,Group Identifier,Assay Type\nLIB1,fbm_1-2,GEX\n"
    sif_path = tmp_path / "tenx_proc_sif_norm.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "processed",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "group_id_normalized_match: 1" in captured.out
    assert "matched after '-'/'_' normalization" in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_psomagen_10x_raw_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a SOP-complete 10x Psomagen raw mapping (viral_ORF assay)."""
    mapping_text = _tenx_fastq_cli_mapping(
        s3_prefix=(
            "s3://czi-psomagen/project-scaling-alpha/AN00000001/CD4i_R1L01/raw/"
            "416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT"
        ),
        assay_token="viral_ORF",
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "psomagen",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_raw_without_fastq_fails(tmp_path: Path, monkeypatch, capsys) -> None:
    """10x raw mode must fail when no FASTQ rows are present."""
    mapping_text = (
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356-LeS188W1_GEX-Z0083-CAGTGTATTGCTGAT.cram,/local/sample.cram\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "requires FASTQ" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_raw_allows_unmatched_cram_artifacts(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """10x raw mode passes when unmatched CRAMs accompany a SOP-complete bundle."""
    mapping_text = _tenx_fastq_cli_mapping(
        extra_rows=[
            f"{_TENX_FASTQ_CLI_PREFIX}_unmatched.cram-metadata.json,"
            "/local/sample_unmatched.cram-metadata.json"
        ]
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_cram_valid_bundle_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """10x_cram raw mode should pass for complete sample and run-level artifacts."""
    prefix = (
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356-LeS188W1_GEX-Z0083-CAGTGTATTGCTGAT"
    )
    mapping_text = (
        f"{prefix}.cram,/local/sample.cram\n"
        f"{prefix}.csv,/local/sample.csv\n"
        f"{prefix}.json,/local/sample.json\n"
        f"{prefix}_extract_stats.h5,/local/sample_extract_stats.h5\n"
        f"{prefix}_SNVQ.metric,/local/sample_SNVQ.metric\n"
        f"{prefix}_FlowQ.metric,/local/sample_FlowQ.metric\n"
        f"{prefix}_trimmer-stats.csv,/local/sample_trimmer-stats.csv\n"
        f"{prefix}_trimmer-failure-codes.csv,/local/sample_trimmer-failure-codes.csv\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356_LibraryInfo.xml,/local/442356_LibraryInfo.xml\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356_UploadCompleted.json,/local/442356_UploadCompleted.json\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "run_SecondaryAnalysis.txt,/local/run_SecondaryAnalysis.txt\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "run_VariantCalling.txt,/local/run_VariantCalling.txt\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356_merged_trimmer-stats.csv,/local/442356_merged_trimmer-stats.csv\n"
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356_merged_trimmer-failure_codes.csv,/local/442356_merged_trimmer-failure_codes.csv\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x_cram",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "10x_cram raw SOP" in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_cram_forbids_unmatched_files(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """10x_cram mode should fail if unmatched CRAM artifacts are present."""
    mapping_text = (
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356-LeS188W1_GEX-Z0083-CAGTGTATTGCTGAT_unmatched.cram,/local/unmatched.cram\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x_cram",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "forbidden_unmatched_cram" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_cram_missing_cram_fails(tmp_path: Path, monkeypatch, capsys) -> None:
    """10x_cram mode should fail clearly when the CRAM artifact is absent."""
    mapping_text = _tenx_cram_cli_mapping(omit_suffixes={".cram"})

    with pytest.raises(SystemExit) as excinfo:
        _run_10x_cram_cli(tmp_path, monkeypatch, mapping_text)

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "missing required sample artifacts: cram" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_cram_swapped_local_artifact_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """A row mapping an S3 JSON artifact to a local CRAM artifact should fail."""
    mapping_text = _tenx_cram_cli_mapping(
        local_overrides={".json": "/provider/export/sample.cram"}
    )

    with pytest.raises(SystemExit) as excinfo:
        _run_10x_cram_cli(tmp_path, monkeypatch, mapping_text)

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "s3_local_artifact_mismatch" in captured.out
    assert "json" in captured.out
    assert "cram" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_cram_different_local_basename_same_artifact_passes(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Different local basenames are allowed when artifact suffixes still agree."""
    mapping_text = _tenx_cram_cli_mapping(
        local_overrides={".json": "/provider/export/provider_renamed.json"}
    )

    with pytest.raises(SystemExit) as excinfo:
        _run_10x_cram_cli(tmp_path, monkeypatch, mapping_text)

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "s3_local_artifact_mismatch" not in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_10x_cram_groupid_mismatch_is_clear(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Filename GroupID mismatches should produce a direct group_mismatch error."""
    mismatched_row = (
        "s3://czi-novogene/project-alpha/NVUS0000000000-29/LeS188W1/raw/"
        "442356-LeS188_GEX-Z0083-CAGTGTATTGCTGAT.cram,/provider/export/mismatched.cram"
    )
    mapping_text = _tenx_cram_cli_mapping(extra_rows=[mismatched_row])

    with pytest.raises(SystemExit) as excinfo:
        _run_10x_cram_cli(tmp_path, monkeypatch, mapping_text)

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "group_mismatch" in captured.out
    assert "LeS188W1" in captured.out
    assert "LeS188" in captured.out
    assert "missing required sample artifacts" not in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_cram_extensionless_sample_message(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Extensionless sample rows should explain that the suffix is missing."""
    mapping_text = _tenx_cram_cli_mapping(
        extra_rows=[f"{_TENX_CRAM_CLI_PREFIX},/provider/export/sample"]
    )

    with pytest.raises(SystemExit) as excinfo:
        _run_10x_cram_cli(tmp_path, monkeypatch, mapping_text)

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "missing 10x_cram sample file suffix" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_scale_raw_truncated_codes_csv_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Reproduce the production bug: ``_codes.csv`` (S3) vs ``_trimmer-failure_codes.csv`` (local)."""
    bug_row = (
        "s3://czi-novogene/lab-seahub-alpha/NVUS0000000000-58/CHEM5-R114/raw/443546/"
        "443546-R114ABCD_hash_oligo_QSR-9-SCALEPLEX_codes.csv,"
        "/local_root/DATA1/V129/443546-20260418_0814/"
        "443546_1-QSR9SCALEPLEX_QSR-9-SCALEPLEX_trimmer-failure_codes.csv"
    )
    mapping_path = _write_temp_mapping(tmp_path, bug_row + "\n")

    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM5-R114,CHEM5-R114_GEX,CTATGCACA,lab-seahub-alpha,"
        "CHEM5-R114,R114ABCD,hash_oligo\n"
    )
    sif_path = tmp_path / "scale_sif_bug.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "scale",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "scale raw completeness" in captured.out
    assert "unexpected_sample_file" in captured.out
    assert "_codes.csv" in captured.out
    assert "s3_local_artifact_mismatch" in captured.out
    assert "scale raw completeness errors" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_scale_raw_hyphen_prefixed_per_ug_passes(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Scale mapping with hyphen-prefixed per-UG suffixes should pass."""
    s3_pfx = (
        "s3://czi-novogene/lab-seahub-alpha/NVUS0000000000-66/"
        "GENE12-R117/raw/443602/443602-R117E_GEX_QSR-5"
    )
    local_pfx = "/local_root/DATA1/V129/443602-20260509_0826/443602-QSR5_QSR-5/443602-QSR5_QSR-5"
    lines = []
    for ext in (".cram", ".csv", ".json"):
        lines.append(f"{s3_pfx}-6B{ext},{local_pfx}_6B{ext}")
    for sfx, local_sfx in (
        ("-trimmer-failure-codes.csv", "_trimmer-failure_codes.csv"),
        ("-trimmer-stats.csv", "_trimmer-stats.csv"),
        ("-unmatched.cram", "_unmatched.cram"),
        ("-unmatched.csv", "_unmatched.csv"),
        ("-unmatched.json", "_unmatched.json"),
    ):
        lines.append(f"{s3_pfx}{sfx},{local_pfx}{local_sfx}")
    mapping_text = "\n".join(lines) + "\n"
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "GENE12-R117,GENE12-R117_GEX,CTATGCACA,lab-seahub-alpha,"
        "GENE12-R117,R117E,GEX\n"
    )
    sif_path = tmp_path / "scale_sif_hyphen.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "scale",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "scale raw completeness" in captured.out
    assert "VERDICT: PASS" in captured.out


def test_cli_scale_raw_missing_per_well_artifact_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Dropping a per-well .json should be reported as missing per-well artifact."""
    mapping_path = _write_temp_mapping(
        tmp_path, _scale_cli_mapping(omit_per_well={".json"})
    )
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM13-R096,CHEM13-R096_GEX,CTATGCACA,lab-seahub-alpha,"
        "CHEM13-R096,R096A,GEX\n"
    )
    sif_path = tmp_path / "scale_sif_missing.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "scale",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "missing required per-well Scale artifacts: json" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_sci_raw_missing_trimmer_failure_codes_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Dropping ``_trimmer-failure_codes.csv`` should fail the sci CLI run."""
    mapping_path = _write_temp_mapping(
        tmp_path, _sci_cli_mapping(omit={"_trimmer-failure_codes.csv"})
    )
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM3-R100,R100E,Z0028,lab-seahub-beta,CHEM3-R100,R100E,GEX_hash_oligo\n"
    )
    sif_path = tmp_path / "sci_sif_missing.csv"
    sif_path.write_text(sif_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--sif",
        str(sif_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "sci",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert (
        "missing required sci sample artifacts: trimmer_failure_codes" in captured.out
    )
    assert "VERDICT: FAIL" in captured.out


def test_cli_10x_raw_missing_per_tail_sample_fastq_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Dropping the R1 ``_sample.fastq.gz`` should fail the 10x raw CLI run."""
    mapping_path = _write_temp_mapping(
        tmp_path,
        _tenx_fastq_cli_mapping(omit_per_tail={("S1_L001_R1_001", "_sample.fastq.gz")}),
    )

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "missing required 10x FASTQ per-tail artifacts: sample_fastq" in captured.out
    assert "10x FASTQ completeness errors" in captured.out
    assert "VERDICT: FAIL" in captured.out


_PER_TAIL_ERROR_DETAIL = "missing required 10x FASTQ per-tail artifacts: sample_fastq"


def _multi_bundle_10x_mapping(n_bundles: int) -> str:
    """Build N independent 10x raw bundles, each missing R1 ``_sample.fastq.gz``.

    Every bundle is otherwise SOP-complete, so each contributes exactly one
    ``missing_sample_artifacts`` (per-tail) completeness error.  Distinct
    GroupIDs keep both the S3 and local paths unique so the uniqueness check
    stays green and only the targeted completeness error is produced.
    """
    lines: list[str] = []
    for i in range(n_bundles):
        gid = f"CD4i_R1L{i:02d}"
        s3_prefix = (
            f"s3://czi-novogene/project-scaling-alpha/NVUS0000000000-29/{gid}/raw/"
            f"416640-{gid}_GEX-Z0238-CTGCACATTGTAGAT"
        )
        local_prefix = f"/local/{gid}/GEX"
        for sfx in (
            ".csv",
            ".json",
            "_trimmer-stats.csv",
            "_trimmer-failure_codes.csv",
            "_unmatched.cram",
            "_unmatched.csv",
            "_unmatched.json",
        ):
            lines.append(f"{s3_prefix}{sfx},{local_prefix}{sfx}")
        for tail in ("S1_L001_R1_001", "S1_L001_R2_001"):
            for sfx in (".fastq.gz", ".csv", ".json", "_sample.fastq.gz"):
                if tail == "S1_L001_R1_001" and sfx == "_sample.fastq.gz":
                    continue
                lines.append(f"{s3_prefix}_{tail}{sfx},{local_prefix}_{tail}{sfx}")
    return "\n".join(lines) + "\n"


def _run_10x_raw_cli(
    tmp_path: Path, monkeypatch, mapping_text: str, extra_args: list[str] | None = None
) -> int:
    """Run the CLI in 10x raw mode and return the captured exit code."""
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)
    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "novogene",
        "--data",
        "raw",
        "--assay",
        "10x",
    ]
    argv.extend(extra_args or [])
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()
    return int(excinfo.value.code)


def test_cli_default_caps_completeness_examples(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Default output caps per-tail completeness examples at five with a hint."""
    code = _run_10x_raw_cli(tmp_path, monkeypatch, _multi_bundle_10x_mapping(8))

    assert code == 1
    captured = capsys.readouterr()
    assert "missing_sample_artifacts: 8" in captured.out
    assert "Showing up to 5 example errors" in captured.out
    assert captured.out.count(_PER_TAIL_ERROR_DETAIL) == 5
    assert "... and 3 more errors (use --verbose to see all)" in captured.out
    assert "VERDICT: FAIL" in captured.out


def test_cli_verbose_shows_all_completeness_errors(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """--verbose prints every completeness error and drops the truncation hint."""
    code = _run_10x_raw_cli(
        tmp_path, monkeypatch, _multi_bundle_10x_mapping(8), extra_args=["--verbose"]
    )

    assert code == 1
    captured = capsys.readouterr()
    assert "missing_sample_artifacts: 8" in captured.out
    assert "Showing all 8 errors" in captured.out
    assert captured.out.count(_PER_TAIL_ERROR_DETAIL) == 8
    assert "use --verbose to see all" not in captured.out


def test_cli_verbose_includes_present_and_source_lines(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Verbose aggregate errors expose present artifacts and source line numbers."""
    _run_10x_raw_cli(
        tmp_path, monkeypatch, _multi_bundle_10x_mapping(1), extra_args=["--verbose"]
    )

    captured = capsys.readouterr()
    assert "present: " in captured.out
    assert "source lines: " in captured.out


def test_cli_verbose_lists_concrete_missing_file_names(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Verbose output reconstructs the concrete SOP-expected missing S3 paths."""
    _run_10x_raw_cli(
        tmp_path, monkeypatch, _multi_bundle_10x_mapping(1), extra_args=["--verbose"]
    )

    captured = capsys.readouterr()
    assert "missing files (SOP-expected, not found in mapping):" in captured.out
    assert (
        "s3://czi-novogene/project-scaling-alpha/NVUS0000000000-29/CD4i_R1L00/raw/"
        "416640-CD4i_R1L00_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001_sample.fastq.gz"
        in captured.out
    )


def test_cli_max_examples_respected(tmp_path: Path, monkeypatch, capsys) -> None:
    """--max-examples N limits examples to exactly N when more errors exist."""
    code = _run_10x_raw_cli(
        tmp_path,
        monkeypatch,
        _multi_bundle_10x_mapping(8),
        extra_args=["--max-examples", "2"],
    )

    assert code == 1
    captured = capsys.readouterr()
    assert "Showing up to 2 example errors" in captured.out
    assert captured.out.count(_PER_TAIL_ERROR_DETAIL) == 2
    assert "... and 6 more errors (use --verbose to see all)" in captured.out


def test_cli_max_examples_zero_shows_all(tmp_path: Path, monkeypatch, capsys) -> None:
    """--max-examples 0 is treated as unlimited, like --verbose."""
    code = _run_10x_raw_cli(
        tmp_path,
        monkeypatch,
        _multi_bundle_10x_mapping(8),
        extra_args=["--max-examples", "0"],
    )

    assert code == 1
    captured = capsys.readouterr()
    assert "Showing all 8 errors" in captured.out
    assert captured.out.count(_PER_TAIL_ERROR_DETAIL) == 8
    assert "use --verbose to see all" not in captured.out


def test_cli_unsupported_mode_fails(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should fail with an unsupported provider/data/assay combination."""
    mapping_text = "s3://bucket/a,/local/a\n"
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        "psomagen",
        "--data",
        "raw",
        "--assay",
        "scale",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "VERDICT: FAIL" in captured.out
    assert "unsupported validation mode" in captured.out
