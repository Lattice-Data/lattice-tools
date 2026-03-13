"""
Tests for mapping_validation.py core helpers.
"""

from pathlib import Path

import pytest

from mapping_validation import (
    ASSAYS_10X,
    ASSAYS_SCALE,
    CANONICAL_ASSAY,
    MappingRow,
    _is_run_metadata,
    _normalize_sif_groupid,
    build_assay_regex,
    compare_groupid_assays,
    find_unmatched_sif_paths_10x,
    get_assays,
    get_order_pattern,
    load_sif_library_assays,
    load_sif_scale_group_assays,
    parse_mapping_file,
    validate_library_assay_consistency,
    validate_local_paths_scale_raw,
    validate_local_paths_sci_raw,
    validate_s3_10x_processed,
    validate_s3_local_consistency_10x_processed,
    validate_s3_local_consistency_sci,
    validate_s3_seahub_raw,
    validate_sif_completeness_scale,
    validate_sif_completeness_seahub,
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
    """Invalid assay and group ID mismatch currently surface as a parse miss.

    The path has a mismatch between the GroupID directory (CD4i_R1L01) and the
    filename (CD4i_R1L02), and uses the non-SOP assay name 'Hash_oliga'. The
    validator does not attempt a fuzzy parse in this case and instead reports a
    generic parse_miss warning.
    """
    # groupid in path is CD4i_R1L01, but filename uses CD4i_R1L02 and assay typo 'Hash_oliga'
    row = MappingRow(
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/CD4i_R1L01/raw/"
        "416640-CD4i_R1L02_Hash_oliga-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
        "/local/path",
        5,
    )

    result = validate_s3_10x_raw("novogene", [row])

    # No hard errors are raised; instead we get a parse_miss warning.
    assert not result["errors"]
    warn_types = {w["type"] for w in result["warnings"]}
    assert "parse_miss" in warn_types


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
    """SIF completeness should correctly parse GroupIDs from Hash_oligo paths.

    This case intentionally uses the non-SOP assay spelling 'Hash_oligo';
    the neighbouring test with 'hash_oligo' in the SIF ensures both the
    casing error path and the group-id extraction behaviour are exercised.
    """
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


# ---------------------------------------------------------------------------
# Library-assay consistency (10x)
# ---------------------------------------------------------------------------


def test_validate_library_assay_consistency_happy_path() -> None:
    """Matching library name → assay and GroupID yields zero mismatches."""
    lib_assays = {"LIB1": "gex", "LIB1F": "cri"}
    rows = [
        MappingRow(
            s3_path="s3://czi-novogene/proj/NVUS2024101701-28/LIB1_LIB1F/raw/100-LIB1_LIB1F_GEX-Z0001-ACGT_R1.fastq.gz",
            local_path="/data/100-20260101_0000/100-LIB1-Z0001-ACGT/100-LIB1-Z0001-ACGT_R1.fastq.gz",
            line_num=1,
        ),
        MappingRow(
            s3_path="s3://czi-novogene/proj/NVUS2024101701-28/LIB1_LIB1F/raw/101-LIB1_LIB1F_CRI-Z0002-TGCA_R1.fastq.gz",
            local_path="/data/101-20260101_0000/101-LIB1F-Z0002-TGCA/101-LIB1F-Z0002-TGCA_R1.fastq.gz",
            line_num=2,
        ),
    ]
    res = validate_library_assay_consistency(rows, lib_assays, "novogene")
    assert res["checked"] == 2
    assert len(res["assay_mismatches"]) == 0
    assert len(res["groupid_mismatches"]) == 0


def test_validate_library_assay_consistency_detects_assay_mismatch() -> None:
    """CRI library appearing under GEX assay is flagged."""
    lib_assays = {"LIB1": "gex", "LIB1F": "cri"}
    rows = [
        MappingRow(
            s3_path="s3://czi-novogene/proj/NVUS2024101701-28/LIB1_LIB1F/raw/100-LIB1_LIB1F_GEX-Z0001-ACGT_R1.fastq.gz",
            local_path="/data/100-20260101_0000/100-LIB1F-Z0001-ACGT/100-LIB1F-Z0001-ACGT_R1.fastq.gz",
            line_num=1,
        ),
    ]
    res = validate_library_assay_consistency(rows, lib_assays, "novogene")
    assert res["checked"] == 1
    assert len(res["assay_mismatches"]) == 1
    m = res["assay_mismatches"][0]
    assert m["library"] == "LIB1F"
    assert m["s3_assay"] == "GEX"
    assert m["expected_assay"] == "CRI"


def test_validate_library_assay_consistency_detects_groupid_mismatch() -> None:
    """Library in local path must appear in the S3 GroupID."""
    lib_assays = {"LIB1": "gex", "LIB1F": "cri", "LIB2": "gex", "LIB2F": "cri"}
    rows = [
        MappingRow(
            s3_path="s3://czi-novogene/proj/NVUS2024101701-28/LIB2_LIB2F/raw/100-LIB2_LIB2F_GEX-Z0001-ACGT_R1.fastq.gz",
            local_path="/data/100-20260101_0000/100-LIB1-Z0001-ACGT/100-LIB1-Z0001-ACGT_R1.fastq.gz",
            line_num=1,
        ),
    ]
    res = validate_library_assay_consistency(rows, lib_assays, "novogene")
    assert res["checked"] == 1
    assert len(res["groupid_mismatches"]) == 1
    m = res["groupid_mismatches"][0]
    assert m["library"] == "LIB1"
    assert m["s3_groupid"] == "LIB2_LIB2F"


def test_normalize_sif_groupid_rewrites_space_plus() -> None:
    """_normalize_sif_groupid should map 'A + AF' → 'A_AF'."""
    assert _normalize_sif_groupid("A + AF") == "A_AF"
    # Plain IDs without ' + ' should be unchanged
    assert _normalize_sif_groupid("R112A") == "R112A"


def test_validate_s3_10x_raw_counts_metadata_files() -> None:
    """10x validator should count run-level metadata separately from data files."""
    base = (
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/CD4i_R1L01/raw/"
    )
    rows = [
        MappingRow(
            s3_path=base
            + "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
            local_path="/local/file1.fastq.gz",
            line_num=1,
        ),
        MappingRow(
            s3_path=base + "416640_UploadCompleted.json",
            local_path="/local/UploadCompleted.json",
            line_num=2,
        ),
    ]

    result = validate_s3_10x_raw("novogene", rows)

    assert result["matched"] == 1
    assert result["metadata_files"] == 1


def test_validate_s3_10x_raw_invalid_barcode_and_project_naming() -> None:
    """10x raw validator should treat malformed paths as parse_miss.

    With the current SOP-aligned regex, this path does not match the expected
    10x raw layout at all, so it is reported as a parse_miss warning rather
    than reaching project_naming or barcode checks.
    """
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/Weissman_Scaling/NVUS2024101701-29/CD4i_R1L01/raw/"
                "416640-CD4i_R1L01_GEX-Z0238-CTGCAXTATTGTAGAT_S1_L001_R1_001.fastq.gz"
            ),
            local_path="/local/file.fastq.gz",
            line_num=1,
        )
    ]

    result = validate_s3_10x_raw("novogene", rows)

    error_types = {e["type"] for e in result["errors"]}
    warn_types = {w["type"] for w in result["warnings"]}
    assert not error_types
    assert "parse_miss" in warn_types


def test_validate_s3_10x_raw_group_mismatch_error() -> None:
    """10x raw validator should flag group_mismatch when GroupID disagrees."""
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
                "NVUS2024101701-29/CD4i_R1L01/raw/"
                "416640-CD4i_R1L02_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz"
            ),
            local_path="/local/file.fastq.gz",
            line_num=2,
        )
    ]

    result = validate_s3_10x_raw("novogene", rows)

    error_types = {e["type"] for e in result["errors"]}
    assert "group_mismatch" in error_types


def test_find_unmatched_sif_paths_10x_reports_unmatched_groupids() -> None:
    """find_unmatched_sif_paths_10x should separate matched and unmatched GroupIDs."""
    base = (
        "s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/"
        "NVUS2024101701-29/"
    )
    # One mapping with GroupID present in SIF, one with GroupID missing from SIF
    rows = [
        MappingRow(
            s3_path=base
            + "CD4i_R1L01/raw/"
            "416640-CD4i_R1L01_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
            local_path="/local/ok.fastq.gz",
            line_num=1,
        ),
        MappingRow(
            s3_path=base
            + "CD4i_R1L02/raw/"
            "416640-CD4i_R1L02_GEX-Z0238-CTGCACATTGTAGAT_S1_L001_R1_001.fastq.gz",
            local_path="/local/extra.fastq.gz",
            line_num=2,
        ),
    ]
    sif_groupids = {"CD4i_R1L01"}

    res = find_unmatched_sif_paths_10x(rows, sif_groupids, "novogene")

    assert res["matched_sif"] == 1
    assert res["metadata"] == 0
    assert len(res["unparsed"]) == 0
    assert set(res["unmatched_by_group"].keys()) == {"CD4i_R1L02"}


def test_validate_s3_scale_raw_index_pattern_happy_path() -> None:
    """Scale S3 index-direct form should match and have no errors.

    This exercises the backward-compatible wrapper validate_s3_scale_raw;
    a separate test below covers the unified validate_s3_seahub_raw("scale", ...)
    entrypoint so both call paths stay wired up.
    """
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_GEX_CTATGCACA.json"
        ),
        local_path="/local/path.json",
        line_num=1,
    )

    result = validate_s3_scale_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_s3_local_consistency_scale_qsr_only_one_side_warning() -> None:
    """When QSR numbers appear only on S3 side, warn but do not error."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_GEX_QSR-1-7A.json"
        ),
        local_path=(
            "/ORPROJ1/DATA1/V129/441969-20260220_2053/"
            "441969-run/441969-run_7A.json"
        ),
        line_num=12,
    )

    res = validate_s3_local_consistency_scale([row])

    assert not res["errors"]
    warn_types = {w["type"] for w in res["warnings"]}
    assert "qsr_only_one_side" in warn_types


def test_validate_s3_local_consistency_scale_qsr_partial_overlap_warning() -> None:
    """Overlapping but non-identical QSR sets should yield a partial-mismatch warning."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_GEX_QSR-1-7A.json"
        ),
        local_path=(
            "/ORPROJ1/DATA1/V129/441969-20260220_2053/"
            "441969-QSR1_QSR-1_QSR-2/441969-QSR1_QSR-1_QSR-2_7A.json"
        ),
        line_num=13,
    )

    res = validate_s3_local_consistency_scale([row])

    assert not res["errors"]
    warn_types = {w["type"] for w in res["warnings"]}
    assert "qsr_partial_mismatch_s3_local" in warn_types


# ---------------------------------------------------------------------------
# Refactored constants / helpers tests
# ---------------------------------------------------------------------------


def test_get_assays_returns_base_set_for_10x() -> None:
    """get_assays('10x') should return the family base without provider extras."""
    assays = get_assays("10x")
    assert assays == ASSAYS_10X
    assert assays == {"GEX", "CRI", "ATAC"}


def test_get_assays_extends_with_provider_extras() -> None:
    """get_assays('10x', 'psomagen') should include viral_ORF."""
    assays = get_assays("10x", "psomagen")
    assert "viral_ORF" in assays
    assert assays == {"GEX", "CRI", "ATAC", "viral_ORF"}


def test_get_assays_novogene_equals_base_10x() -> None:
    """Novogene has no extras for 10x, so result equals the base set."""
    assert get_assays("10x", "novogene") == ASSAYS_10X


def test_get_assays_unknown_family_raises() -> None:
    """Unknown assay family should raise ValueError."""
    with pytest.raises(ValueError, match="Unknown assay family"):
        get_assays("unknown_family")


def test_get_order_pattern_unknown_provider_raises() -> None:
    """Unknown provider for order pattern should raise ValueError."""
    with pytest.raises(ValueError, match="Unknown provider"):
        get_order_pattern("unknown_provider")


def test_build_assay_regex_longest_first() -> None:
    """build_assay_regex should sort longest-first to avoid partial matches."""
    import re

    regex = build_assay_regex(ASSAYS_SCALE)
    m = re.match(regex, "GEX_hash_oligo")
    assert m is not None
    assert m.group(0) == "GEX_hash_oligo"


def test_compare_groupid_assays_detects_missing_and_extra() -> None:
    """compare_groupid_assays should report missing/extra GroupIDs and assays."""
    expected = {
        "G1": {"gex", "cri"},
        "G2": {"gex"},
    }
    actual = {
        "G1": {"gex"},
        "G3": {"atac"},
    }

    result = compare_groupid_assays(expected, actual)

    assert result["missing_groupids"] == {"G2"}
    assert result["extra_groupids"] == {"G3"}
    assert "G1" in result["missing_assays"]
    assert result["missing_assays"]["G1"] == {"cri"}


# ---------------------------------------------------------------------------
# Seahub unified validator tests
# ---------------------------------------------------------------------------


def test_validate_s3_seahub_raw_scale_happy_path() -> None:
    """Seahub unified validator with scale family should match SOP form."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441969-R096G_GEX_CTATGCACA.json"
        ),
        local_path="/local/path.json",
        line_num=1,
    )

    result = validate_s3_seahub_raw("scale", [row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_s3_seahub_raw_scale_runid_and_gex_scaleplex_mismatch() -> None:
    """Scale S3 validator should flag runid and GEX/SCALEPLEX violations."""
    # runid in directory (441969) disagrees with runid2 in filename (441970)
    # and GEX assay is incorrectly combined with SCALEPLEX in UG_RT.
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-26/CHEM13-R096/raw/441969/"
            "441970-R096A_GEX_QSR-1-SCALEPLEX-8G.cram"
        ),
        local_path="/local/path.cram",
        line_num=2,
    )

    result = validate_s3_seahub_raw("scale", [row])

    error_types = {e["type"] for e in result["errors"]}
    assert "runid_mismatch" in error_types
    assert "gex_scaleplex_violation" in error_types


def test_validate_s3_seahub_raw_sci_happy_path() -> None:
    """Seahub unified validator with sci family should match Z-barcode form."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        local_path="/local/path",
        line_num=1,
    )

    result = validate_s3_seahub_raw("sci", [row])

    assert result["matched"] == 1
    assert result["errors"] == []
    assert result["group_assays"] == {"R100E": {"GEX_hash_oligo"}}


def test_validate_s3_seahub_raw_sci_invalid_barcode_and_project_naming() -> None:
    """sci S3 validator should treat malformed paths as parse_miss.

    The current sci regex focuses on SOP layout; this fixture does not match
    that layout, so it results in a parse_miss warning rather than a more
    specific project_naming or barcode error.
    """
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/Hamazaki-Seahub-BCP/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCXAT.json"
        ),
        local_path="/local/path",
        line_num=2,
    )

    result = validate_s3_seahub_raw("sci", [row])

    error_types = {e["type"] for e in result["errors"]}
    warn_types = {w["type"] for w in result["warnings"]}
    assert not error_types
    assert "parse_miss" in warn_types


def test_validate_s3_seahub_raw_sci_tracks_group_assays() -> None:
    """Seahub sci mode should collect group_assays from parsed paths."""
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
                "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT.json"
            ),
            local_path="/local/a",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441908/"
                "441908-R100D_GEX_hash_oligo-Z0155-CTTCATATCTGAGAT.csv"
            ),
            local_path="/local/b",
            line_num=2,
        ),
    ]

    result = validate_s3_seahub_raw("sci", rows)

    assert result["matched"] == 2
    assert "R100E" in result["group_assays"]
    assert "R100D" in result["group_assays"]


def test_validate_s3_seahub_raw_sci_metadata_separated() -> None:
    """Metadata files in sci mode should be counted separately."""
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
                "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT.json"
            ),
            local_path="/local/data",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
                "441389_LibraryInfo.xml"
            ),
            local_path="/local/meta",
            line_num=2,
        ),
    ]

    result = validate_s3_seahub_raw("sci", rows)

    assert result["matched"] == 1
    assert result["metadata_files"] == 1


def test_validate_sif_completeness_seahub_sci(tmp_path: Path) -> None:
    """Seahub SIF completeness for sci should detect missing GroupIDs."""
    sif_text = (
        "Library name,Sublibrary name,Ultima Index Sequence,Project Identifier,"
        "Experiement Identifier,Group Identifier,Assay Type\n"
        "CHEM3-R100,R100E,Z0028,hamazaki-seahub-bcp,CHEM3-R100,R100E,GEX_hash_oligo\n"
        "CHEM3-R100,R100D,Z0155,hamazaki-seahub-bcp,CHEM3-R100,R100D,GEX_hash_oligo\n"
    )
    sif_path = tmp_path / "SIF_sci.csv"
    sif_path.write_text(sif_text)

    mappings = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
                "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT.json"
            ),
            local_path="/local/path",
            line_num=1,
        ),
    ]

    result = validate_sif_completeness_seahub("sci", mappings, sif_path)

    assert "R100E" in result["actual_groupids"]
    assert "R100D" in result["missing_groupids"]


def test_validate_local_paths_sci_raw_happy_path() -> None:
    """Well-formed sci local path should produce no errors."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441389-20260224_2053/"
            "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT/"
            "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        line_num=1,
    )

    result = validate_local_paths_sci_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_sci_raw_single_ug_sci_plex_style() -> None:
    """Single-UG sci-plex local path (RunID-GroupID-UG-Barcode) should match and validate."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-30/CHEM16/raw/441588/"
            "441588-CHEM16_P07_F3_GEX_hash_oligo-Z0310-CATGACAGTAATGAT_trimmer-stats.csv"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-754-001/441588-20260218_1318/"
            "441588-CHEM16_P07_F3-Z0310-CATGACAGTAATGAT/"
            "441588-CHEM16_P07_F3-Z0310-CATGACAGTAATGAT_trimmer-stats.csv"
        ),
        line_num=1,
    )

    result = validate_local_paths_sci_raw([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_local_paths_sci_raw_detects_barcode_mismatch() -> None:
    """Mismatched barcode between directory and filename should be flagged."""
    row = MappingRow(
        s3_path="s3://czi-novogene/proj/NVUS2024101701-32/EXP/raw/441389/441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT.json",
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441389-20260224_2053/"
            "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT/"
            "441389-R100E_Z0028-Z0028-AAAAAAAAAAAAGAT_SNVQ.metric"
        ),
        line_num=2,
    )

    result = validate_local_paths_sci_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "barcode_mismatch" in error_types


def test_validate_local_paths_sci_raw_detects_runid_groupid_ug_mismatch() -> None:
    """sci local validator should flag runid, groupid and UG inconsistencies."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441388-20260224_2053/"
            "441387-R100F_Z0028-Z0029-CAGACTTGCTGCGAT/"
            "441386-R100G_Z0030-Z0031-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        line_num=3,
    )

    result = validate_local_paths_sci_raw([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "runid_mismatch" in error_types
    assert "groupid_mismatch" in error_types
    assert "ug_mismatch" in error_types


def test_validate_s3_local_consistency_sci_happy_path() -> None:
    """Consistent S3 and local paths should produce no errors."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441389-20260224_2053/"
            "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT/"
            "441389-R100E_Z0028-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        line_num=1,
    )

    result = validate_s3_local_consistency_sci([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_s3_local_consistency_sci_single_ug_sci_plex_style() -> None:
    """Single-UG sci-plex local path should match S3 for consistency check."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-30/CHEM16/raw/441588/"
            "441588-CHEM16_P07_F3_GEX_hash_oligo-Z0310-CATGACAGTAATGAT_trimmer-stats.csv"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-754-001/441588-20260218_1318/"
            "441588-CHEM16_P07_F3-Z0310-CATGACAGTAATGAT/"
            "441588-CHEM16_P07_F3-Z0310-CATGACAGTAATGAT_trimmer-stats.csv"
        ),
        line_num=1,
    )

    result = validate_s3_local_consistency_sci([row])

    assert result["matched"] == 1
    assert result["errors"] == []


def test_validate_s3_local_consistency_sci_detects_groupid_mismatch() -> None:
    """Mismatched GroupID between S3 and local should be an error."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT.json"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441389-20260224_2053/"
            "441389-R100F_Z0028-Z0028-CAGACTTGCTGCGAT/"
            "441389-R100F_Z0028-Z0028-CAGACTTGCTGCGAT.json"
        ),
        line_num=3,
    )

    result = validate_s3_local_consistency_sci([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "groupid_mismatch_s3_local" in error_types


def test_validate_s3_local_consistency_sci_detects_runid_ug_barcode_mismatch() -> None:
    """sci S3/local consistency should flag runid, UG and barcode mismatches."""
    row = MappingRow(
        s3_path=(
            "s3://czi-novogene/hamazaki-seahub-bcp/NVUS2024101701-32/CHEM3-R100/raw/441389/"
            "441389-R100E_GEX_hash_oligo-Z0028-CAGACTTGCTGCGAT_SNVQ.metric"
        ),
        local_path=(
            "/ORPROJ1/NEWSFTP/S3/ultima/CR0-789/441388-20260224_2053/"
            "441388-R100E_Z0029-Z0029-AAAAAAAAAAAAAAA/"
            "441388-R100E_Z0029-Z0029-AAAAAAAAAAAAAAA_SNVQ.metric"
        ),
        line_num=4,
    )

    result = validate_s3_local_consistency_sci([row])

    assert result["matched"] == 1
    error_types = {e["type"] for e in result["errors"]}
    assert "runid_mismatch_s3_local" in error_types
    assert "ug_mismatch_s3_local" in error_types
    assert "barcode_mismatch_s3_local" in error_types


# ---------------------------------------------------------------------------
# 10x processed S3 validation
# ---------------------------------------------------------------------------

_PROC_S3_NOVO = (
    "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
    "e10_rep1_t13/processed/cellranger/Run_2025-03-10/outs/"
)
_PROC_LOCAL_NOVO = (
    "/ORPROJ1/GB/USER/liguo/Ultima/projects_202602/"
    "X202SC25127893-Z01-F001_GRCm39-vM37_10xcellranger_v9.0/"
    "Data_process/sampleMatrix/e10_rep1_t13/outs/"
)

_PROC_S3_PSOM = (
    "s3://czi-psomagen/weissman-perturb/AN00012345/"
    "CD4i_R1L01/processed/cellranger/Run_2025-02-01/outs/"
)
_PROC_LOCAL_PSOM = (
    "/data/process/CD4i_R1L01/outs/"
)


def test_validate_s3_10x_processed_novogene_basic() -> None:
    """Valid novogene processed paths should match and produce no errors."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_NOVO + "filtered_feature_bc_matrix.h5",
            local_path=_PROC_LOCAL_NOVO + "filtered_feature_bc_matrix.h5",
            line_num=1,
        ),
        MappingRow(
            s3_path=_PROC_S3_NOVO + "metrics_summary.csv",
            local_path=_PROC_LOCAL_NOVO + "metrics_summary.csv",
            line_num=2,
        ),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["matched"] == 2
    assert res["errors"] == []
    assert res["warnings"] == []
    assert "e10_rep1_t13" in res["group_ids"]
    assert "cellranger" in res["pipelines"]
    assert "Run_2025-03-10" in res["run_dates"]


def test_validate_s3_10x_processed_psomagen_basic() -> None:
    """Valid psomagen processed paths should match."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_PSOM + "multi/count/raw_feature_bc_matrix.h5",
            local_path=_PROC_LOCAL_PSOM + "multi/count/raw_feature_bc_matrix.h5",
            line_num=1,
        ),
    ]
    res = validate_s3_10x_processed("psomagen", rows)
    assert res["matched"] == 1
    assert res["errors"] == []
    assert "CD4i_R1L01" in res["group_ids"]


def test_validate_s3_10x_processed_invalid_pipeline() -> None:
    """An unexpected pipeline should produce an error."""
    bad_s3 = (
        "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
        "e10_rep1_t13/processed/starsolo/Run_2025-03-10/outs/file.h5"
    )
    rows = [
        MappingRow(s3_path=bad_s3, local_path="/local/file.h5", line_num=1),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(e["type"] == "invalid_pipeline" for e in res["errors"])


def test_validate_s3_10x_processed_bad_run_date_format() -> None:
    """A non-standard run date should produce a warning."""
    bad_s3 = (
        "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
        "e10_rep1_t13/processed/cellranger/2025-03-10/outs/file.h5"
    )
    rows = [
        MappingRow(s3_path=bad_s3, local_path="/local/file.h5", line_num=1),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(w["type"] == "run_date_format" for w in res["warnings"])


def test_validate_s3_10x_processed_project_naming_warning() -> None:
    """Non-lowercase project names should produce a warning."""
    bad_s3 = (
        "s3://czi-novogene/Weissman_Embryo/NVUS2024101701-19/"
        "e10_rep1_t13/processed/cellranger/Run_2025-03-10/outs/file.h5"
    )
    rows = [
        MappingRow(s3_path=bad_s3, local_path="/local/file.h5", line_num=1),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(w["type"] == "project_naming" for w in res["warnings"])


def test_validate_s3_10x_processed_unmatched_path() -> None:
    """A path that doesn't follow the processed pattern produces a warning."""
    bad_s3 = (
        "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
        "e10_rep1_t13/raw/440021-e10_rep1_t13_GEX-Z0035-CTGAATGATCTCGAT.csv"
    )
    rows = [
        MappingRow(s3_path=bad_s3, local_path="/local/file.csv", line_num=1),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["matched"] == 0
    assert any(w["type"] == "parse_miss" for w in res["warnings"])


def test_validate_s3_10x_processed_multiple_group_ids() -> None:
    """Multiple GroupIDs should all be collected."""
    rows = [
        MappingRow(
            s3_path=(
                "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
                "e9_rep1_t1/processed/cellranger/Run_2025-03-10/outs/file.h5"
            ),
            local_path="/local/e9_rep1_t1/outs/file.h5",
            line_num=1,
        ),
        MappingRow(
            s3_path=(
                "s3://czi-novogene/weissman-embryo-tracing/NVUS2024101701-19/"
                "e10_rep1_t13/processed/cellranger/Run_2025-03-10/outs/file.h5"
            ),
            local_path="/local/e10_rep1_t13/outs/file.h5",
            line_num=2,
        ),
    ]
    res = validate_s3_10x_processed("novogene", rows)
    assert res["group_ids"] == {"e9_rep1_t1", "e10_rep1_t13"}


# ---------------------------------------------------------------------------
# 10x processed S3/local consistency
# ---------------------------------------------------------------------------


def test_s3_local_consistency_10x_processed_match() -> None:
    """Consistent S3 and local paths should produce no errors."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_NOVO + "filtered_feature_bc_matrix.h5",
            local_path=_PROC_LOCAL_NOVO + "filtered_feature_bc_matrix.h5",
            line_num=1,
        ),
    ]
    res = validate_s3_local_consistency_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert res["errors"] == []
    assert res["warnings"] == []


def test_s3_local_consistency_10x_processed_group_id_missing() -> None:
    """GroupID from S3 not in local should produce an error."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_NOVO + "filtered_feature_bc_matrix.h5",
            local_path="/data/process/wrong_sample/outs/filtered_feature_bc_matrix.h5",
            line_num=1,
        ),
    ]
    res = validate_s3_local_consistency_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(e["type"] == "group_id_missing_local" for e in res["errors"])


def test_s3_local_consistency_10x_processed_file_path_mismatch() -> None:
    """Different file paths after outs/ should produce an error."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_NOVO + "filtered_feature_bc_matrix.h5",
            local_path=_PROC_LOCAL_NOVO + "raw_feature_bc_matrix.h5",
            line_num=1,
        ),
    ]
    res = validate_s3_local_consistency_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(e["type"] == "file_path_mismatch" for e in res["errors"])


def test_s3_local_consistency_10x_processed_no_outs_in_local() -> None:
    """Local path without /outs/ should produce a warning."""
    rows = [
        MappingRow(
            s3_path=_PROC_S3_NOVO + "filtered_feature_bc_matrix.h5",
            local_path="/data/e10_rep1_t13/filtered_feature_bc_matrix.h5",
            line_num=1,
        ),
    ]
    res = validate_s3_local_consistency_10x_processed("novogene", rows)
    assert res["matched"] == 1
    assert any(w["type"] == "no_outs_in_local" for w in res["warnings"])

