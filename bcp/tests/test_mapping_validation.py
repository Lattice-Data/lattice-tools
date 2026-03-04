"""
Tests for mapping_validation.py core helpers.
"""

from pathlib import Path

import pytest

from mapping_validation import (
    ASSAYS_SCALE,
    CANONICAL_ASSAY,
    MappingRow,
    _is_run_metadata,
    load_sif_scale_group_assays,
    parse_mapping_file,
    validate_local_paths_scale_raw,
    validate_sif_completeness_scale,
    validate_s3_scale_raw,
    validate_s3_10x_raw,
    validate_s3_local_consistency_scale,
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
    """S3 Scale raw paths with 'Hash_oliga' should be treated as invalid assay.

    The typo prevents regex matching, but the error must still be reported.
    """
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096D_Hash_oliga_CACTGCTCA.json"
        ),
        local_path="/local/path",
        line_num=3,
    )

    result = validate_s3_scale_raw([row])

    # The typo means the regex won't match, but the error is still captured
    assert result["matched"] == 0
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


def test_validate_s3_local_consistency_scale_qsr_mismatch() -> None:
    """Clearly different QSR numbers between S3 and local should be an error."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_GEX_QSR-1-7A.json"
        ),
        local_path=(
            "/ORPROJ1/DATA1/V129/441969-20260220_2053/"
            "441969-QSR7_QSR-7/441969-QSR7_QSR-7_7A.json"
        ),
        line_num=10,
    )

    res = validate_s3_local_consistency_scale([row])

    error_types = {e["type"] for e in res["errors"]}
    assert "qsr_mismatch_s3_local" in error_types


def test_validate_s3_local_consistency_scale_scaleplex_warning() -> None:
    """SCALEPLEX on S3 but not on local should be a warning, not an error."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096A_hash_oligo_QSR-1-SCALEPLEX-8G.cram"
        ),
        local_path=(
            "/ORPROJ1/DATA1/V129/441969-20260220_2053/"
            "441969-QSR1_QSR-1/441969-QSR1_QSR-1_8G.cram"
        ),
        line_num=11,
    )

    res = validate_s3_local_consistency_scale([row])

    assert not res["errors"]
    warn_types = {w["type"] for w in res["warnings"]}
    assert "scaleplex_missing_in_local" in warn_types


# --- V2 local path pattern tests ---


def test_validate_local_paths_scale_raw_v2_gex_happy_path() -> None:
    """V2 compact-QSR GEX local path should match and produce no errors."""
    row = MappingRow(
        s3_path="s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
        "439774-R112A_GEX_QSR-1-7A.json",
        local_path=(
            "/ORPROJ1/DATA1/V129/439774-20260220_2053/"
            "439774-QSR1_QSR-1/439774-QSR1_QSR-1_7A.json"
        ),
        line_num=1,
    )

    result = validate_local_paths_scale_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_scale_raw_v2_scaleplex_happy_path() -> None:
    """V2 compact-QSR SCALEPLEX local path should match and produce no errors."""
    row = MappingRow(
        s3_path="s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
        "439774-R112A_Hash_oligo_QSR-1-SCALEPLEX-8G.cram",
        local_path=(
            "/ORPROJ1/DATA1/V129/439774-20260220_2053/"
            "439774-QSR1SCALEPLEX_QSR-1-SCALEPLEX/"
            "439774-QSR1SCALEPLEX_QSR-1-SCALEPLEX_8G.cram"
        ),
        line_num=6,
    )

    result = validate_local_paths_scale_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_scale_raw_v2_detects_dir_file_mismatch() -> None:
    """V2 path where directory stem differs from filename stem should be flagged."""
    row = MappingRow(
        s3_path="s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
        "439774-R112A_GEX_QSR-1-7A.json",
        local_path=(
            "/ORPROJ1/DATA1/V129/439774-20260220_2053/"
            "439774-QSR1_QSR-1/439774-QSR2_QSR-2_7A.json"
        ),
        line_num=2,
    )

    result = validate_local_paths_scale_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "dir_file_mismatch" in error_types


# --- Assay-anchored regex tests ---


def test_validate_s3_scale_raw_correctly_parses_hash_oligo_group_id() -> None:
    """Hash_oligo paths should parse group_id as R112A, not R112A_Hash.

    The path uses 'Hash_oligo' (capital H) which violates SOP spelling
    'hash_oligo', so an assay_casing error is expected — but the group_id
    must still be correctly extracted as 'R112A'.
    """
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
            "439774-R112A_Hash_oligo_QSR-1-SCALEPLEX-8G.cram"
        ),
        local_path="/local/path",
        line_num=6,
    )

    result = validate_s3_scale_raw([row])

    assert result["matched"] == 1
    # Only error should be the casing violation — no structural errors
    assert len(result["errors"]) == 1
    assert result["errors"][0]["type"] == "assay_casing"
    assert "'Hash_oligo' violates SOP spelling 'hash_oligo'" in result["errors"][0]["detail"]


def test_validate_sif_completeness_scale_correct_group_ids(tmp_path: Path) -> None:
    """SIF completeness should correctly parse GroupIDs from Hash_oligo paths."""
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "R112A,QSR1,CGATGCGCA,trapnell-seahub-bcp,GENE8-R112,R112A,GEX\n"
        "R112A,QSR1SCALEPLEX,CACATCACA,trapnell-seahub-bcp,GENE8-R112,R112A,Hash_oligo\n"
    )
    sif_path = tmp_path / "SIF.csv"
    sif_path.write_text(sif_text)

    mappings = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_GEX_QSR-1-7A.json"
            ),
            local_path="/local/path",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_Hash_oligo_QSR-1-SCALEPLEX-8G.cram"
            ),
            local_path="/local/path2",
            line_num=6,
        ),
    ]

    result = validate_sif_completeness_scale(mappings, sif_path)

    # Both GEX and Hash_oligo should parse group_id as R112A (not R112A_Hash)
    assert result["actual_groupids"] == {"R112A"}
    assert result["missing_groupids"] == set()
    assert result["extra_groupids"] == set()
    # Both assay types present → no missing assays
    assert result["missing_assays"] == {}
    assert result["extra_assays"] == {}


# --- SIF per-GroupID assay completeness tests ---


def test_load_sif_scale_group_assays_returns_assay_mapping(tmp_path: Path) -> None:
    """SIF loader should return a dict of GroupID → set of assay types."""
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "R112A,QSR1,CGATGCGCA,trapnell-seahub-bcp,GENE8-R112,R112A,GEX\n"
        "R112A,QSR1SCALEPLEX,CACATCACA,trapnell-seahub-bcp,GENE8-R112,R112A,Hash_oligo\n"
        "R112B,QSR2,CGCATATCA,trapnell-seahub-bcp,GENE8-R112,R112B,GEX\n"
        "R112B,QSR2SCALEPLEX,CTGCAGTGA,trapnell-seahub-bcp,GENE8-R112,R112B,Hash_oligo\n"
    )
    sif_path = tmp_path / "SIF.csv"
    sif_path.write_text(sif_text)

    result = load_sif_scale_group_assays(sif_path)

    assert result == {
        "R112A": {"gex", "hash_oligo"},
        "R112B": {"gex", "hash_oligo"},
    }


def test_validate_sif_completeness_flags_missing_assay(tmp_path: Path) -> None:
    """If GEX exists but Hash_oligo is missing for a GroupID, flag it."""
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "R112A,QSR1,CGATGCGCA,trapnell-seahub-bcp,GENE8-R112,R112A,GEX\n"
        "R112A,QSR1SCALEPLEX,CACATCACA,trapnell-seahub-bcp,GENE8-R112,R112A,Hash_oligo\n"
    )
    sif_path = tmp_path / "SIF.csv"
    sif_path.write_text(sif_text)

    # Only GEX is in S3 — Hash_oligo is missing
    mappings = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_GEX_QSR-1-7A.json"
            ),
            local_path="/local/path",
            line_num=1,
        ),
    ]

    result = validate_sif_completeness_scale(mappings, sif_path)

    assert result["missing_groupids"] == set()
    assert "R112A" in result["missing_assays"]
    assert "hash_oligo" in result["missing_assays"]["R112A"]


def test_validate_sif_completeness_no_missing_when_all_assays_present(tmp_path: Path) -> None:
    """No missing assays when all SIF assay types are in S3."""
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "R112A,QSR1,CGATGCGCA,trapnell-seahub-bcp,GENE8-R112,R112A,GEX\n"
        "R112A,QSR1SCALEPLEX,CACATCACA,trapnell-seahub-bcp,GENE8-R112,R112A,Hash_oligo\n"
    )
    sif_path = tmp_path / "SIF.csv"
    sif_path.write_text(sif_text)

    mappings = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_GEX_QSR-1-7A.json"
            ),
            local_path="/local/path",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_hash_oligo_QSR-1-SCALEPLEX-8G.cram"
            ),
            local_path="/local/path2",
            line_num=6,
        ),
    ]

    result = validate_sif_completeness_scale(mappings, sif_path)

    assert result["missing_assays"] == {}
    assert result["extra_assays"] == {}


# --- Metadata file recognition tests ---


def test_is_run_metadata_recognizes_known_files() -> None:
    """Known run-level metadata files should be recognized."""
    assert _is_run_metadata("s3://bucket/path/raw/439774/439774_LibraryInfo.xml")
    assert _is_run_metadata("s3://bucket/path/raw/439774/439774_SequencingInfo.json")
    assert _is_run_metadata("s3://bucket/path/raw/439774/UploadCompleted.json")
    assert _is_run_metadata("s3://bucket/path/raw/439774/merged_trimmer-stats.csv")
    assert _is_run_metadata("s3://bucket/path/raw/439774/run_SecondaryAnalysis.txt")


def test_is_run_metadata_rejects_sample_files() -> None:
    """Sample data files should not be recognized as metadata."""
    assert not _is_run_metadata(
        "s3://bucket/path/raw/439774/439774-R112A_GEX_QSR-1-7A.json"
    )
    assert not _is_run_metadata(
        "s3://bucket/path/raw/439774/439774-R112A_Hash_oligo_QSR-1-SCALEPLEX-8G.cram"
    )


def test_validate_s3_scale_raw_separates_metadata_files() -> None:
    """Metadata files should be counted separately, not as parse misses."""
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774-R112A_GEX_QSR-1-7A.json"
            ),
            local_path="/local/path1",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-35/GENE8-R112/raw/439774/"
                "439774_LibraryInfo.xml"
            ),
            local_path="/local/path2",
            line_num=2,
        ),
    ]

    result = validate_s3_scale_raw(rows)

    assert result["matched"] == 1
    assert result["metadata_files"] == 1
    assert len(result["warnings"]) == 0
