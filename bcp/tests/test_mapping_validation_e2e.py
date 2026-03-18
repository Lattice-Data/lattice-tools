from __future__ import annotations

from pathlib import Path
import sys

import pytest

import mapping_validation
import pandas as pd


FIXTURES_DIR = Path(__file__).parent / "fixtures" / "mapping_validation"


@pytest.mark.parametrize(
    "mapping_name,sif_name,provider,data,assay,expected_code",
    [
        # Happy paths
        (
            "novogene_10x_raw_valid.csv",
            "novogene_10x_sif.csv",
            "novogene",
            "raw",
            "10x",
            0,
        ),
        (
            "psomagen_10x_raw_valid.csv",
            "psomagen_10x_sif.csv",
            "psomagen",
            "raw",
            "10x",
            0,
        ),
        (
            "novogene_sci_raw_valid.csv",
            "novogene_sci_sif.csv",
            "novogene",
            "raw",
            "sci",
            0,
        ),
        (
            "novogene_scale_raw_sop_valid.csv",
            "novogene_scale_sif.csv",
            "novogene",
            "raw",
            "scale",
            0,
        ),
        (
            "novogene_10x_processed_valid.csv",
            "novogene_10x_processed_sif.csv",
            "novogene",
            "processed",
            "10x",
            0,
        ),
        # Error paths
        ("duplicates.csv", None, "novogene", "raw", "10x", 1),
        (
            "novogene_10x_raw_errors.csv",
            "novogene_10x_sif.csv",
            "novogene",
            "raw",
            "10x",
            1,
        ),
        (
            "novogene_sci_raw_errors.csv",
            "novogene_sci_sif.csv",
            "novogene",
            "raw",
            "sci",
            1,
        ),
        (
            "novogene_scale_raw_errors.csv",
            "novogene_scale_sif.csv",
            "novogene",
            "raw",
            "scale",
            1,
        ),
        (
            "novogene_10x_processed_errors.csv",
            "novogene_10x_processed_sif.csv",
            "novogene",
            "processed",
            "10x",
            1,
        ),
        (
            "novogene_sci_raw_missing_groupid.csv",
            "novogene_sci_sif.csv",
            "novogene",
            "raw",
            "sci",
            1,
        ),
    ],
)
def test_mapping_validation_e2e(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    mapping_name: str,
    sif_name: str | None,
    provider: str,
    data: str,
    assay: str,
    expected_code: int,
) -> None:
    """End-to-end CLI tests against small mapping/SIF fixtures."""
    # Copy mapping fixture into a temp file to avoid coupling to repo layout
    src_mapping = FIXTURES_DIR / "mappings" / mapping_name
    mapping_path = tmp_path / mapping_name
    mapping_path.write_text(src_mapping.read_text(encoding="utf-8"), encoding="utf-8")

    argv = [
        "mapping_validation",
        "--mapping",
        str(mapping_path),
        "--provider",
        provider,
        "--data",
        data,
        "--assay",
        assay,
    ]

    if sif_name is not None:
        src_sif = FIXTURES_DIR / "sif" / sif_name
        if src_sif.suffix == ".csv":
            # Real SIFs are typically Excel workbooks; convert the small CSV
            # fixtures into .xlsx files so we exercise the Excel loading path.
            df = pd.read_csv(src_sif)
            sif_path = tmp_path / src_sif.with_suffix(".xlsx").name
            df.to_excel(sif_path, index=False)
        else:
            sif_path = tmp_path / sif_name
            sif_path.write_text(src_sif.read_text(encoding="utf-8"), encoding="utf-8")
        argv.extend(["--sif", str(sif_path)])

    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == expected_code
    captured = capsys.readouterr()
    if expected_code == 0:
        assert "VERDICT: PASS" in captured.out
    else:
        assert "VERDICT: FAIL" in captured.out
