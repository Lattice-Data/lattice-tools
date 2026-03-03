"""
Tests for mapping_validation.py core helpers.
"""

from pathlib import Path

import pytest

from mapping_validation import (
    MappingRow,
    parse_mapping_file,
    validate_local_paths_scale_raw,
    validate_sif_completeness_scale,
    validate_s3_scale_raw,
    validate_s3_10x_raw,
    validate_uniqueness,
)


def _write_temp_mapping(tmp_path: Path, contents: str) -> Path:
    path = tmp_path / "mapping.csv"
    path.write_text(contents)
    return path


def test_parse_mapping_file_supports_comma_and_tab(tmp_path: Path) -> None:
    """parse_mapping_file should understand both comma- and tab-separated lines."""
    mapping_text = (
        "s3://bucket/project/order/group/raw/file1.fastq.gz,/local/path/file1.fastq.gz\n"
        "s3://bucket/project/order/group/raw/file2.fastq.gz\t/local/path/file2.fastq.gz\n"
        "\n"
        "# comment-like line that should be ignored\n"
    )
    path = _write_temp_mapping(tmp_path, mapping_text)

    rows = parse_mapping_file(path)

    assert len(rows) == 2
    assert rows[0].s3_path.endswith("file1.fastq.gz")
    assert rows[0].local_path.endswith("file1.fastq.gz")
    assert rows[1].s3_path.endswith("file2.fastq.gz")
    assert rows[1].local_path.endswith("file2.fastq.gz")


def test_parse_mapping_file_ignores_malformed_lines(tmp_path: Path) -> None:
    """Lines that cannot be split into two columns are skipped."""
    mapping_text = (
        "this is not a valid mapping line\n"
        "s3://bucket/a,/local/a\n"
        "s3://bucket/b /local/b\n"  # space-separated fallback
    )
    path = _write_temp_mapping(tmp_path, mapping_text)

    rows = parse_mapping_file(path)

    # only the last two lines should be parsed
    assert [r.s3_path for r in rows] == [
        "s3://bucket/a",
        "s3://bucket/b",
    ]


def test_validate_uniqueness_detects_duplicate_local_and_s3() -> None:
    """validate_uniqueness should flag many-to-one and one-to-many mappings."""
    rows = [
        MappingRow("s3://bucket/a", "/local/a", 1),
        MappingRow("s3://bucket/b", "/local/b", 2),
        # duplicate local path (two S3s -> one local)
        MappingRow("s3://bucket/c", "/local/a", 3),
        # duplicate S3 path (two locals -> one S3)
        MappingRow("s3://bucket/b", "/local/b2", 4),
    ]

    result = validate_uniqueness(rows)

    dup_locals = result["duplicate_locals"]
    dup_s3 = result["duplicate_s3"]

    assert "/local/a" in dup_locals
    assert {r.line_num for r in dup_locals["/local/a"]} == {1, 3}

    assert "s3://bucket/b" in dup_s3
    assert {r.line_num for r in dup_s3["s3://bucket/b"]} == {2, 4}


def test_validate_uniqueness_ok_for_one_to_one() -> None:
    """No duplicates when each S3/local pair is unique."""
    rows = [
        MappingRow("s3://bucket/a", "/local/a", 1),
        MappingRow("s3://bucket/b", "/local/b", 2),
    ]

    result = validate_uniqueness(rows)

    assert result["duplicate_locals"] == {}
    assert result["duplicate_s3"] == {}
    assert result["total"] == 2


def test_validate_s3_10x_raw_novogene_happy_path() -> None:
    """A correctly formatted Novogene 10x raw S3 path should produce no errors."""
    rows = [
        MappingRow(
            "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
            "NVUS2024101701-29/CD4i_R1L01/raw/"
            "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
            "/local/416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
            1,
        )
    ]

    result = validate_s3_10x_raw("novogene", rows)

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_s3_10x_raw_flags_invalid_assay_and_group_mismatch() -> None:
    """Invalid assay and group ID mismatch should be reported as errors."""
    # groupid in path is CD4i_R1L01, but filename uses CD4i_R1L02 and assay typo 'Hash_oliga'
    row = MappingRow(
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/CD4i_R1L01/raw/"
        "416640-CD4i_R1L02_Hash_oliga-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
        "/local/path",
        5,
    )

    result = validate_s3_10x_raw("novogene", [row])

    error_types = {e["type"] for e in result["errors"]}
    assert "group_mismatch" in error_types
    assert "invalid_assay" in error_types


def test_validate_s3_10x_raw_psomagen_allows_viral_orf() -> None:
    """Psomagen 10x validator should accept viral_ORF assay."""
    row = MappingRow(
        "s3://czi-psomagen/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "AN00012345/CD4i_R1L01/raw/"
        "416640-CD4i_R1L01_viral_ORF-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
        "/local/path",
        10,
    )

    result = validate_s3_10x_raw("psomagen", [row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_scale_raw_happy_path() -> None:
    """Well-formed Scale local path should produce no errors."""
    row = MappingRow(
        s3_path="s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
        "441969-R096G_GEX_CTATGCACA.json",
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441969-20260220_0135/"
            "441969-QSR-7-CTATGCACA/441969-QSR-7-CTATGCACA.json"
        ),
        line_num=1,
    )

    result = validate_local_paths_scale_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_scale_raw_detects_qsr_mismatch() -> None:
    """Mismatched QSR number between directory and filename should be flagged."""
    row = MappingRow(
        s3_path="s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
        "441969-R096G_GEX_CTATGCACA.json",
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441969-20260220_0135/"
            "441969-QSR-7-CTATGCACA/441969-QSR-8-CTATGCACA.json"
        ),
        line_num=2,
    )

    result = validate_local_paths_scale_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "qsr_mismatch" in error_types


def test_validate_s3_scale_raw_flags_hash_oliga_typo() -> None:
    """S3 Scale raw paths with 'Hash_oliga' should be treated as invalid assay."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096D_Hash_oliga_CACTGCTCA.json"
        ),
        local_path="/local/path",
        line_num=3,
    )

    result = validate_s3_scale_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "invalid_assay" in error_types


def test_validate_s3_scale_raw_requires_scaleplex_for_hash_oligo() -> None:
    """Hash_oligo SOP form without SCALEPLEX in UG_RT should be flagged."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_Hash_oligo_QSR-7.json"
        ),
        local_path="/local/path",
        line_num=4,
    )

    result = validate_s3_scale_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "hash_oligo_scaleplex_mismatch" in error_types


def test_validate_sif_completeness_scale_flags_missing_groupid(tmp_path: Path) -> None:
    """SIF completeness should report GroupIDs present in SIF but absent in S3."""
    # Create a tiny SIF with two group identifiers
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,Experiement Identifier,"
        "Group Identifier,Assay Type\n"
        "CHEM13-R096,R096A,QSR1,trapnell-seahub-bcp,CHEM13-R096,R096A,GEX\n"
        "CHEM13-R096,R096B,QSR2,trapnell-seahub-bcp,CHEM13-R096,R096B,GEX\n"
    )
    sif_path = tmp_path / "SIF_scale.csv"
    sif_path.write_text(sif_text)

    # Only R096B appears in S3 mappings; R096A should be reported missing
    mappings = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
                "441969-R096B_GEX_CTATGCACA.json"
            ),
            local_path="/local/path",
            line_num=1,
        )
    ]

    result = validate_sif_completeness_scale(mappings, sif_path)

    assert "R096B" in result["actual_groupids"]
    assert "R096A" in result["missing_groupids"]
