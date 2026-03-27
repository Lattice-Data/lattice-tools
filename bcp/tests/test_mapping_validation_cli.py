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


def test_cli_passes_for_valid_10x_mapping(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should exit 0 for a simple valid 10x Novogene mapping."""
    mapping_text = (
        "s3://czi-novogene/project-scaling-alpha/"
        "NVUS0000000000-29/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz,"
        "/local/416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz\n"
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R2_001.fastq.gz,"
        "/local/416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R2_001.fastq.gz\n"
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


def test_cli_scale_raw_with_sif_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a minimal valid Scale raw mapping with SIF."""
    mapping_text = (
        "s3://czi-novogene/lab-seahub-alpha/NVUS0000000000-26/CHEM13-R096/raw/441969/"
        "441969-R096A_GEX_QSR-1-7A.json,"
        "/local_root/data1/V129/441969-20260220_2053/"
        "441969-QSR1_QSR-1/441969-QSR1_QSR-1_7A.json\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

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
    assert "VERDICT: PASS" in captured.out


def test_cli_sci_raw_with_sif_passes(tmp_path: Path, monkeypatch, capsys) -> None:
    """CLI should pass for a minimal valid sci raw mapping with SIF."""
    mapping_text = (
        "s3://czi-novogene/lab-seahub-beta/NVUS0000000000-32/CHEM3-R100/raw/441389/"
        "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric,"
        "/local_root/newsftp/S3/ultima/CR0-789/441389-20260224_2053/"
        "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT/"
        "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT_SNVQ.metric\n"
    )
    mapping_path = _write_temp_mapping(tmp_path, mapping_text)

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
    """CLI should pass for a minimal valid 10x Psomagen raw mapping."""
    mapping_text = (
        "s3://czi-psomagen/project-scaling-alpha/"
        "AN00000001/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz,"
        "/local/416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz\n"
        "s3://czi-psomagen/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "AN00012345/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT_S1_L001_R2_001.fastq.gz,"
        "/local/416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT_S1_L001_R2_001.fastq.gz\n"
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
