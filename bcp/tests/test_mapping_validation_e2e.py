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
        (
            "novogene_10x_cram_raw_valid.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            0,
        ),
        (
            "novogene_10x_cram_raw_different_local_basename_valid.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            0,
        ),
        (
            "psomagen_10x_illumina_raw_valid.csv",
            "psomagen_10x_illumina_sif.csv",
            "psomagen",
            "raw",
            "10x_illumina",
            0,
        ),
        (
            "psomagen_10x_processed_valid.csv",
            "psomagen_10x_processed_sif.csv",
            "psomagen",
            "processed",
            "10x",
            0,
        ),
        (
            "novogene_10x_raw_multiome_valid.csv",
            "novogene_10x_multiome_sif.csv",
            "novogene",
            "raw",
            "10x",
            0,
        ),
        # Multiome order whose local folders are named after the SIF Library
        # names (CH01GEX/CH01ATAC) while S3 uses the SOP stem (CH01_GEX) under
        # GroupID CH01.  The library name is *longer* than the GroupID here,
        # which is the reverse of the LIB1/LIB1_LIB1F layout.
        (
            "novogene_10x_raw_multiome_local_lib_names.csv",
            "novogene_10x_multiome_sif.csv",
            "novogene",
            "raw",
            "10x",
            0,
        ),
        # Same idea for the paired layout: local folders named LIB1/LIB1F while
        # S3 files both under GroupID LIB1_LIB1F.
        (
            "novogene_10x_raw_local_lib_names.csv",
            "novogene_10x_sif.csv",
            "novogene",
            "raw",
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
            "novogene_10x_cram_raw_unmatched_forbidden.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            1,
        ),
        (
            "novogene_10x_cram_raw_missing_cram.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            1,
        ),
        (
            "novogene_10x_cram_raw_swapped_local_artifact.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            1,
        ),
        (
            "novogene_10x_cram_raw_groupid_mismatch.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
            1,
        ),
        (
            "novogene_10x_cram_raw_extensionless_sample.csv",
            None,
            "novogene",
            "raw",
            "10x_cram",
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
        (
            "psomagen_10x_illumina_raw_missing_run_meta.csv",
            None,
            "psomagen",
            "raw",
            "10x_illumina",
            1,
        ),
        (
            "psomagen_10x_illumina_raw_incomplete_reads.csv",
            None,
            "psomagen",
            "raw",
            "10x_illumina",
            1,
        ),
        (
            "psomagen_10x_illumina_raw_ultima_stem.csv",
            None,
            "psomagen",
            "raw",
            "10x_illumina",
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


@pytest.mark.parametrize(
    "mapping_name,sif_name,expected_checked",
    [
        # Multiome: SIF library name extends the GroupID (CH01GEX -> CH01).
        (
            "novogene_10x_raw_multiome_local_lib_names.csv",
            "novogene_10x_multiome_sif.csv",
            34,
        ),
        # Paired: GroupID concatenates its members (LIB1, LIB1F -> LIB1_LIB1F).
        ("novogene_10x_raw_local_lib_names.csv", "novogene_10x_sif.csv", 30),
    ],
)
def test_mapping_validation_e2e_library_consistency_actually_runs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    mapping_name: str,
    sif_name: str,
    expected_checked: int,
) -> None:
    """Library consistency must inspect every row, not silently check nothing.

    A plain exit-code assertion cannot tell "all rows agree" from "no rows were
    looked at": every other 10x fixture reports ``checked 0 paths`` because its
    local paths embed the SOP stem rather than the bare library name.  Pinning
    the count keeps a regression back to a silent no-op visible.
    """
    mapping_path = tmp_path / mapping_name
    mapping_path.write_text(
        (FIXTURES_DIR / "mappings" / mapping_name).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    sif_path = tmp_path / "sif.xlsx"
    pd.read_csv(FIXTURES_DIR / "sif" / sif_name).to_excel(sif_path, index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
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
            "10x",
        ],
    )

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert (
        f"Library consistency: checked {expected_checked} paths, "
        "0 assay mismatches, 0 GroupID mismatches" in out
    )
    assert "VERDICT: PASS" in out


def test_mapping_validation_e2e_reports_groupid_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A misfiled row must name both the SIF and the S3 GroupID, and fail.

    One CH01GEX file is filed under GroupID CH02.  The run reports other
    problems too (S3 is now missing CH02's other assay), but the library
    consistency line has to say which library was expected where.
    """
    mapping_path = tmp_path / "m.csv"
    mapping_path.write_text(
        (
            FIXTURES_DIR / "mappings" / "novogene_10x_raw_multiome_groupid_mismatch.csv"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    sif_path = tmp_path / "sif.xlsx"
    pd.read_csv(
        FIXTURES_DIR / "sif" / "novogene_10x_multiome_two_group_sif.csv"
    ).to_excel(sif_path, index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
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
            "10x",
        ],
    )

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    out = capsys.readouterr().out
    assert "1 GroupID mismatches" in out
    assert "GroupID mismatches (library's GroupID does not match S3 GroupID):" in out
    assert (
        "local library 'CH01GEX' belongs to SIF GroupID 'CH01' "
        "but S3 GroupID is 'CH02'" in out
    )
    assert "VERDICT: FAIL" in out


def test_mapping_validation_e2e_notes_missing_sif_groupids(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A SIF with no Group Identifier column must say the check was weakened."""
    mapping_path = tmp_path / "m.csv"
    mapping_path.write_text(
        (
            FIXTURES_DIR / "mappings" / "novogene_10x_raw_multiome_local_lib_names.csv"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    sif_path = tmp_path / "sif.xlsx"
    pd.DataFrame(
        {"Library name": ["CH01GEX", "CH01ATAC"], "Assay Type": ["GEX", "ATAC"]}
    ).to_excel(sif_path, index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
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
            "10x",
        ],
    )

    with pytest.raises(SystemExit):
        mapping_validation.main()

    out = capsys.readouterr().out
    assert "Library consistency: checked 34 paths" in out
    assert "NOTE: 34 of 34 paths had no SIF GroupID" in out


def test_mapping_validation_e2e_10x_processed_multiome_pass(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Processed Multiome mapping with --tenx-profile multiome should pass."""
    mapping_path = tmp_path / "m.csv"
    mapping_path.write_text(
        (
            FIXTURES_DIR / "mappings" / "novogene_10x_processed_multiome_valid.csv"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    df = pd.read_csv(FIXTURES_DIR / "sif" / "novogene_10x_multiome_sif.csv")
    sif_path = tmp_path / "multiome_sif.xlsx"
    df.to_excel(sif_path, index=False)

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
        "--tenx-profile",
        "multiome",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 0
    assert "multiome processed outs" in capsys.readouterr().out.lower()


def test_mapping_validation_e2e_10x_processed_multiome_incomplete_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Missing required ATAC outs with --tenx-profile multiome should fail."""
    mapping_path = tmp_path / "m.csv"
    mapping_path.write_text(
        (
            FIXTURES_DIR / "mappings" / "novogene_10x_processed_multiome_incomplete.csv"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    df = pd.read_csv(FIXTURES_DIR / "sif" / "novogene_10x_multiome_sif.csv")
    sif_path = tmp_path / "multiome_sif.xlsx"
    df.to_excel(sif_path, index=False)

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
        "--tenx-profile",
        "multiome",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        mapping_validation.main()

    assert excinfo.value.code == 1
    out = capsys.readouterr().out
    assert "VERDICT: FAIL" in out
    assert "multiome" in out.lower()
