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
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz,"
        "/local/416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz\n"
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
    mapping_text = (
        "s3://bucket/a,/local/a\n"
        "s3://bucket/a,/local/b\n"
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
    assert "VERDICT: FAIL" in captured.out

