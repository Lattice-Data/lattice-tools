"""Tests for file_extract.sheets and the submission-sheet outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from file_extract.constants import (
    SEQUENCE_FILE_COLUMNS,
    SEQUENCE_FILE_SET_COLUMNS,
)
from file_extract.cram import extract_cram
from file_extract.fastq import extract_fastq
from file_extract.sheets import (
    FileSetGroup,
    _project,
    LabIdentity,
    SequenceFileRecord,
    SheetBuildError,
    SheetOptions,
    build_cram_record,
    build_fastq_record,
    cram_set_stem,
    cro_order_mismatch_warning,
    default_sheet_output_names,
    enrich_record,
    group_records,
    parse_fastq_slot,
    run_cardinality_for,
    sample_dir_for,
    sequence_file_row,
    sequence_file_set_row,
    validate_aliases,
    validate_cro_order,
    validate_plan,
    validate_read_slots,
    write_sheet,
)
from tests.file_extract_helpers import FIXTURES, MockS3Client

BUCKET = "example-bucket"
PREFIX = "test-order/NVUS0000000000-01/"

# Header rows copied verbatim from the Lattice submission workbooks. Kept literal
# here on purpose: this is the contract the emitted TSVs have to satisfy, so the
# test must not read it from the same constant the code does.
REFERENCE_SEQUENCE_FILE_HEADER = [
    "#response",
    "#response_time",
    "uuid",
    "aliases",
    "lab",
    "file_format",
    "s3_uri",
    "crc64nvme_base64",
    "no_file_available",
    "derived_from",
    "status",
    "submitter_comment",
    "description",
    "read_count",
    "file_size",
]

REFERENCE_SEQUENCE_FILE_SET_HEADER = [
    "#response",
    "#response_time",
    "uuid",
    "aliases",
    "lab",
    "library",
    "run_cardinality",
    "status",
    "submitter_comment",
    "description",
    "CRO_order",
    "is_pilot_order",
    "read1",
    "read2",
    "read3",
    "index1",
    "index2",
    "untrimmed_cram",
    "trimmed_cram",
    "sequencing_platform",
]


def _record(
    key: str,
    *,
    size: int = 100,
    read_count: int | None = None,
) -> SequenceFileRecord:
    record = build_fastq_record(key=key, s3_uri=f"s3://{BUCKET}/{key}", file_size=size)
    return enrich_record(record, crc="crc", read_count=read_count)


def _options(tmp_path: Path, **overrides: object) -> SheetOptions:
    defaults: dict[str, object] = {
        "lab": LabIdentity.parse("example-lab"),
        "cro_order": "NVUS0000000000-01",
        "is_pilot_order": False,
        "sequencing_platform": "Ultima Genomics UG 100",
        "sequence_file_path": tmp_path / "sf.tsv",
        "sequence_file_set_path": tmp_path / "sfs.tsv",
    }
    defaults.update(overrides)
    return SheetOptions(**defaults)  # type: ignore[arg-type]


def _read_tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


def test_sheet_columns_match_reference_workbooks() -> None:
    assert SEQUENCE_FILE_COLUMNS == REFERENCE_SEQUENCE_FILE_HEADER
    assert SEQUENCE_FILE_SET_COLUMNS == REFERENCE_SEQUENCE_FILE_SET_HEADER


@pytest.mark.parametrize(
    ("value", "name", "path"),
    [
        ("example-lab", "example-lab", "/labs/example-lab/"),
        ("/labs/example-lab/", "example-lab", "/labs/example-lab/"),
        ("/labs/example-lab", "example-lab", "/labs/example-lab/"),
        ("  other-lab  ", "other-lab", "/labs/other-lab/"),
    ],
)
def test_lab_identity_parse(value: str, name: str, path: str) -> None:
    lab = LabIdentity.parse(value)
    assert (lab.name, lab.path, lab.namespace) == (name, path, name)


def test_lab_identity_custom_namespace() -> None:
    lab = LabIdentity.parse("/labs/other-lab/", namespace="sims-lab")
    assert lab.path == "/labs/other-lab/"
    assert lab.namespace == "sims-lab"


@pytest.mark.parametrize(
    "value", ["", "   ", "has space", "a/b", "peter:sims", "/labs/"]
)
def test_lab_identity_rejects_junk(value: str) -> None:
    with pytest.raises(SheetBuildError):
        LabIdentity.parse(value)


def test_lab_identity_rejects_junk_namespace() -> None:
    with pytest.raises(SheetBuildError):
        LabIdentity.parse("other-lab", namespace="bad namespace")


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        (
            "435612-L_eri_3primeLib_GEX-Z0050-CACATGGCAGCACAGAT_S1_L001_R1_001.fastq.gz",
            ("435612-L_eri_3primeLib_GEX-Z0050-CACATGGCAGCACAGAT_S1_L001", "R1", "001"),
        ),
        ("sample_S1_L002_R2_003.fastq.gz", ("sample_S1_L002", "R2", "003")),
        ("sample_L001_R3_001.fastq.gz", ("sample_L001", "R3", "001")),
        ("sample_L001_I2_001.fastq.gz", ("sample_L001", "I2", "001")),
        ("chunkless_R1.fastq.gz", ("chunkless", "R1", "")),
        ("no_designator.fastq.gz", ("", "", "")),
        ("sample_R4_001.fastq.gz", ("", "", "")),
    ],
)
def test_parse_fastq_slot(filename: str, expected: tuple[str, str, str]) -> None:
    assert parse_fastq_slot(filename) == expected


def test_cram_set_stem() -> None:
    assert cram_set_stem("436387-R097D_GEX-Z0097.cram") == "436387-R097D_GEX-Z0097"
    assert cram_set_stem("no-extension") == "no-extension"


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        ("proj/NVUS-15/L_eri_3primeLib/raw/x_R1_001.fastq.gz", "L_eri_3primeLib"),
        ("proj/NVUS-09/R097/raw/436387/436387-x.cram", "R097"),
        ("proj/NVUS-15/L_eri/other/x_R1_001.fastq.gz", "other"),
    ],
)
def test_sample_dir_for(key: str, expected: str) -> None:
    assert sample_dir_for(key) == expected


@pytest.mark.parametrize(
    ("slots", "expected"),
    [
        (["R1", "R2"], "paired-end"),
        (["R1"], "single-end"),
        ([""], "single-end"),
        (["R1", "R2", "R3"], "paired-end"),
        (["R1", "I1", "I2"], "single-end"),
    ],
)
def test_run_cardinality_for(slots: list[str], expected: str) -> None:
    assert run_cardinality_for(slots) == expected


def test_group_records_pairs_reads() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
    ]
    groups, warnings = group_records(records)
    assert warnings == []
    assert len(groups) == 1
    assert groups[0].set_stem == "s_L001"
    assert [member.slot for member in groups[0].members] == ["R1", "R2"]


def test_group_records_splits_lanes() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L002_R1_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L002_R2_001.fastq.gz"),
    ]
    groups, _ = group_records(records)
    assert [g.set_stem for g in groups] == ["s_L001", "s_L002"]


def test_group_records_splits_sample_directories() -> None:
    records = [
        _record(f"{PREFIX}S2/raw/b_L001_R1_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/a_L001_R1_001.fastq.gz"),
    ]
    groups, _ = group_records(records)
    assert [(g.directory.rsplit("/", 2)[-2], g.set_stem) for g in groups] == [
        ("S1", "a_L001"),
        ("S2", "b_L001"),
    ]


def test_group_records_warns_on_unpaired_r2() -> None:
    groups, warnings = group_records(
        [_record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz")]
    )
    assert len(groups) == 1
    assert any("R2 present without R1" in w for w in warnings)


def test_group_records_skips_file_without_designator() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/orphan.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
    ]
    groups, warnings = group_records(records)
    assert len(groups) == 1
    assert any("no read designator" in w for w in warnings)


def test_group_records_warns_on_read_count_mismatch() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz", read_count=100),
        _record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz", read_count=99),
    ]
    _, warnings = group_records(records)
    assert any("possible truncated file" in w for w in warnings)


def _multi_chunk_records() -> list[SequenceFileRecord]:
    """One lane delivered as two chunks per read -- what we refuse to submit."""
    return [
        _record(f"{PREFIX}S1/raw/s_L001_{read}_{chunk}.fastq.gz")
        for chunk in ("001", "002")
        for read in ("R1", "R2")
    ]


def test_validate_read_slots_rejects_chunked_reads() -> None:
    """Chunks are pieces of one read, not separate sequencing runs.

    Splitting them into two sets would claim two runs happened; a set cannot hold
    both in one slot, so the delivery is refused and named instead.
    """
    with pytest.raises(SheetBuildError) as exc:
        validate_read_slots(_multi_chunk_records())
    message = str(exc.value)
    assert "Chunked reads are not supported" in message
    assert "concatenated" in message
    assert "s_L001 slot R1 <- chunks 001, 002" in message
    assert "s_L001 slot R2 <- chunks 001, 002" in message


def test_validate_read_slots_accepts_one_file_per_slot() -> None:
    validate_read_slots(
        [
            _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
            _record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz"),
            _record(f"{PREFIX}S1/raw/s_L001_I1_001.fastq.gz"),
            # Same read, different lane: a different set, not a duplicate slot.
            _record(f"{PREFIX}S1/raw/s_L002_R1_001.fastq.gz"),
        ]
    )


def test_validate_read_slots_rejects_chunkless_beside_chunked() -> None:
    with pytest.raises(SheetBuildError, match=r"chunks \(none\), 001"):
        validate_read_slots(
            [
                _record(f"{PREFIX}S1/raw/s_L001_R1.fastq.gz"),
                _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
            ]
        )


def test_validate_plan_rejects_chunks_before_alias_check() -> None:
    """Chunked input is named as such rather than as a confusing alias clash."""
    with pytest.raises(SheetBuildError, match="Chunked reads are not supported"):
        validate_plan(_multi_chunk_records(), namespace="example-lab")


def test_group_records_rejects_duplicate_slot() -> None:
    """Backstop for callers that group without validating first."""
    with pytest.raises(SheetBuildError, match="Chunked reads are not supported"):
        group_records(_multi_chunk_records())


def test_validate_aliases_accepts_distinct_files() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/a_L001_R1_001.fastq.gz"),
        _record(f"{PREFIX}S2/raw/b_L001_R1_001.fastq.gz"),
    ]
    validate_aliases(records, namespace="example-lab")


def test_validate_aliases_rejects_same_basename_in_two_dirs() -> None:
    records = [
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
        _record(f"{PREFIX}S2/raw/s_L001_R1_001.fastq.gz"),
    ]
    with pytest.raises(SheetBuildError) as exc:
        validate_aliases(records, namespace="example-lab")
    message = str(exc.value)
    assert "example-lab:s_L001_R1_001.fastq.gz" in message
    # Both the file alias and the set alias collide, so both are reported.
    assert "example-lab:s_L001" in message


def test_sequence_file_row_values() -> None:
    record = _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz", size=42, read_count=7)
    row = sequence_file_row(record, lab=LabIdentity.parse("example-lab"))
    assert len(row) == len(SEQUENCE_FILE_COLUMNS)
    values = dict(zip(SEQUENCE_FILE_COLUMNS, row))
    assert values["aliases"] == '["example-lab:s_L001_R1_001.fastq.gz"]'
    assert values["lab"] == "/labs/example-lab/"
    assert values["file_format"] == "fastq"
    assert values["no_file_available"] == "FALSE"
    assert values["status"] == "current"
    assert values["read_count"] == 7
    assert values["file_size"] == 42
    assert values["uuid"] == ""
    assert values["derived_from"] == ""
    assert values["#response"] == ""


def test_sequence_file_row_leaves_failed_enrichment_blank() -> None:
    record = enrich_record(
        build_fastq_record(
            key=f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz", s3_uri="s3://b/k", file_size=1
        ),
        crc=None,
        read_count=None,
    )
    values = dict(
        zip(
            SEQUENCE_FILE_COLUMNS, sequence_file_row(record, lab=LabIdentity.parse("x"))
        )
    )
    assert values["crc64nvme_base64"] == ""
    assert values["read_count"] is None


def test_sequence_file_set_row_values(tmp_path: Path) -> None:
    records = [
        _record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz"),
        _record(f"{PREFIX}S1/raw/s_L001_I1_001.fastq.gz"),
    ]
    groups, _ = group_records(records)
    row = sequence_file_set_row(
        groups[0], options=_options(tmp_path, is_pilot_order=True)
    )
    assert len(row) == len(SEQUENCE_FILE_SET_COLUMNS)
    values = dict(zip(SEQUENCE_FILE_SET_COLUMNS, row))
    assert values["aliases"] == '["example-lab:s_L001"]'
    assert values["run_cardinality"] == "paired-end"
    assert values["read1"] == "example-lab:s_L001_R1_001.fastq.gz"
    assert values["read2"] == "example-lab:s_L001_R2_001.fastq.gz"
    assert values["index1"] == "example-lab:s_L001_I1_001.fastq.gz"
    assert values["read3"] == ""
    assert values["CRO_order"] == "NVUS0000000000-01"
    assert values["is_pilot_order"] == "TRUE"
    assert values["sequencing_platform"] == "Ultima Genomics UG 100"
    assert values["status"] == "current"
    # Not derivable from the S3 layout; filled by hand from the diagnostic TSV.
    assert values["library"] == ""
    assert values["untrimmed_cram"] == ""
    assert values["trimmed_cram"] == ""


def test_sequence_file_set_row_is_pilot_false(tmp_path: Path) -> None:
    groups, _ = group_records([_record(f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz")])
    row = sequence_file_set_row(groups[0], options=_options(tmp_path))
    values = dict(zip(SEQUENCE_FILE_SET_COLUMNS, row))
    assert values["is_pilot_order"] == "FALSE"
    assert values["run_cardinality"] == "single-end"


@pytest.mark.parametrize(
    ("slot", "column", "other"),
    [
        ("trimmed", "trimmed_cram", "untrimmed_cram"),
        ("untrimmed", "untrimmed_cram", "trimmed_cram"),
    ],
)
def test_sequence_file_set_row_places_cram(
    tmp_path: Path, slot: str, column: str, other: str
) -> None:
    record = build_cram_record(
        key=f"{PREFIX}R097/raw/436387/436387-x.cram",
        s3_uri="s3://b/436387-x.cram",
        file_size=10,
    )
    groups, _ = group_records([record])
    row = sequence_file_set_row(groups[0], options=_options(tmp_path, cram_slot=slot))
    values = dict(zip(SEQUENCE_FILE_SET_COLUMNS, row))
    assert values[column] == "example-lab:436387-x.cram"
    assert values[other] == ""
    assert values["run_cardinality"] == "single-end"
    assert values["read1"] == ""


def test_sequence_file_set_row_requires_cram_slot(tmp_path: Path) -> None:
    record = build_cram_record(
        key=f"{PREFIX}R097/raw/x.cram", s3_uri="s3://b/x.cram", file_size=1
    )
    groups, _ = group_records([record])
    with pytest.raises(SheetBuildError, match="--cram-slot"):
        sequence_file_set_row(groups[0], options=_options(tmp_path))


def test_slot_column_rejects_unknown_read_designator(tmp_path: Path) -> None:
    """Not reachable from the CLI: the parser only yields R1-R3 and I1-I2."""
    group = FileSetGroup(
        set_stem="s_L001",
        directory="d",
        members=(
            SequenceFileRecord(
                key="d/s_L001_R9_001.fastq.gz",
                s3_uri="s3://b/k",
                file_format="fastq",
                file_size=1,
                set_stem="s_L001",
                slot="R9",
                chunk="001",
            ),
        ),
    )
    with pytest.raises(SheetBuildError, match="maps to no set slot"):
        sequence_file_set_row(group, options=_options(tmp_path))


def test_project_rejects_unknown_column() -> None:
    """Guards a typo in a row builder rather than any user input."""
    with pytest.raises(SheetBuildError, match="Unknown sheet column"):
        _project({"not_a_column": "x"}, ["aliases"])


def test_write_sheet_leaves_alias_arrays_unquoted(tmp_path: Path) -> None:
    path = tmp_path / "out.tsv"
    write_sheet(path, ["aliases", "read_count"], [['["lab:file.fastq.gz"]', None]])
    text = path.read_text(encoding="utf-8")
    assert text == 'aliases\tread_count\n["lab:file.fastq.gz"]\t\n'


def test_write_sheet_rejects_wrong_width(tmp_path: Path) -> None:
    with pytest.raises(SheetBuildError, match="2 values but the sheet has 3"):
        write_sheet(tmp_path / "o.tsv", ["a", "b", "c"], [["1", "2"]])


def test_write_sheet_rejects_embedded_tab(tmp_path: Path) -> None:
    with pytest.raises(SheetBuildError, match="contains a tab or newline"):
        write_sheet(tmp_path / "o.tsv", ["a"], [["has\ttab"]])


def test_write_sheet_renders_booleans(tmp_path: Path) -> None:
    path = tmp_path / "out.tsv"
    write_sheet(path, ["is_pilot_order"], [[True], [False]])
    assert _read_tsv(path)[1:] == [["TRUE"], ["FALSE"]]


@pytest.mark.parametrize("value", ["", "  ", "a/b", "a\\b", "with\ttab"])
def test_validate_cro_order_rejects_junk(value: str) -> None:
    with pytest.raises(SheetBuildError):
        validate_cro_order(value)


def test_validate_cro_order_strips() -> None:
    assert validate_cro_order("  NVUS0000000000-15 ") == "NVUS0000000000-15"


def test_cro_order_mismatch_warning() -> None:
    assert cro_order_mismatch_warning("NVUS-15", "bucket/path/NVUS-15/") is None
    assert cro_order_mismatch_warning("nvus-15", "bucket/path/NVUS-15/") is None
    warning = cro_order_mismatch_warning("NVUS-99", "bucket/path/NVUS-15/")
    assert warning is not None
    assert "does not appear in the S3 prefix" in warning


def test_default_sheet_output_names() -> None:
    sequence_file, file_set = default_sheet_output_names("NVUS-15", out_dir="out")
    assert sequence_file == Path("out/NVUS-15_SequenceFile.tsv")
    assert file_set == Path("out/NVUS-15_SequenceFileSet.tsv")


def test_file_set_group_alias() -> None:
    group = FileSetGroup(set_stem="s_L001", directory="d")
    assert group.set_alias("lab") == "lab:s_L001"


def test_extract_fastq_reproduces_submitted_sheets(tmp_path: Path) -> None:
    """End to end over a real 20-file order, pinned to the submitter's own sheet.

    The goldens were derived from a hand-filled SequenceFile TSV that a curator
    already accepted, with bucket, project, order and lab identifiers substituted
    for neutral ones. A diff here means the tool would no longer reproduce that
    submission.
    """
    listing = json.loads(
        (FIXTURES / "order15_listing.json").read_text(encoding="utf-8")
    )
    keys = sorted(listing)
    client = MockS3Client(
        keys=keys,
        sizes={key: listing[key]["size"] for key in keys},
        crc_by_key={key: listing[key]["crc"] for key in keys},
        object_bodies={
            f"{key}-metadata.json": json.dumps(
                {"read_count": listing[key]["read_count"]}
            )
            for key in keys
        },
    )
    options = _options(
        tmp_path,
        cro_order="NVUS0000000000-15",
        sequence_file_path=tmp_path / "sf.tsv",
        sequence_file_set_path=tmp_path / "sfs.tsv",
    )

    summary = extract_fastq(
        client,
        "example-bucket",
        "example-project/NVUS0000000000-15/",
        str(tmp_path / "diagnostics.tsv"),
        show_progress=False,
        inline=True,
        sheets=options,
    )

    assert summary.total == 20
    assert summary.set_count == 10
    assert summary.warnings == []
    assert options.sequence_file_path.read_text(encoding="utf-8") == (
        FIXTURES / "order15_sequence_file_sheet.tsv"
    ).read_text(encoding="utf-8")
    assert options.sequence_file_set_path.read_text(encoding="utf-8") == (
        FIXTURES / "order15_sequence_file_set_sheet.tsv"
    ).read_text(encoding="utf-8")


def test_extract_cram_writes_sheets(tmp_path: Path) -> None:
    key = f"{PREFIX}R097/raw/436387/436387-R097D_GEX-Z0097.cram"
    client = MockS3Client(
        keys=[key],
        sizes={key: 512},
        crc_by_key={key: "crc-cram"},
        object_bodies={f"{key}-metadata.json": json.dumps({"read_count": 139216024})},
    )
    options = _options(tmp_path, cram_slot="trimmed")

    summary = extract_cram(
        client,
        BUCKET,
        PREFIX,
        str(tmp_path / "diagnostics.tsv"),
        show_progress=False,
        inline=True,
        sheets=options,
    )

    assert summary.total == 1
    assert summary.set_count == 1

    file_rows = _read_tsv(options.sequence_file_path)
    assert file_rows[0] == SEQUENCE_FILE_COLUMNS
    file_values = dict(zip(SEQUENCE_FILE_COLUMNS, file_rows[1]))
    assert file_values["file_format"] == "cram"
    assert file_values["aliases"] == ('["example-lab:436387-R097D_GEX-Z0097.cram"]')
    assert file_values["read_count"] == "139216024"
    assert file_values["file_size"] == "512"

    set_rows = _read_tsv(options.sequence_file_set_path)
    assert set_rows[0] == SEQUENCE_FILE_SET_COLUMNS
    set_values = dict(zip(SEQUENCE_FILE_SET_COLUMNS, set_rows[1]))
    assert set_values["aliases"] == '["example-lab:436387-R097D_GEX-Z0097"]'
    assert set_values["trimmed_cram"] == "example-lab:436387-R097D_GEX-Z0097.cram"
    assert set_values["untrimmed_cram"] == ""
    assert set_values["run_cardinality"] == "single-end"


def test_extract_fastq_rejects_chunked_order_before_fetching(tmp_path: Path) -> None:
    """A chunked delivery is refused up front, not after a request per file."""
    names = [
        f"s_L001_{read}_{chunk}.fastq.gz"
        for chunk in ("001", "002")
        for read in ("R1", "R2")
    ]
    keys = [f"{PREFIX}S1/raw/{name}" for name in names]
    client = MockS3Client(
        keys=keys,
        sizes={key: 10 for key in keys},
        crc_by_key={key: "crc" for key in keys},
        object_bodies={f"{key}-metadata.json": '{"read_count": 5}' for key in keys},
    )

    with pytest.raises(SheetBuildError) as exc:
        extract_fastq(
            client,
            BUCKET,
            PREFIX,
            str(tmp_path / "diagnostics.tsv"),
            show_progress=False,
            inline=True,
            sheets=_options(tmp_path),
        )

    message = str(exc.value)
    assert "Chunked reads are not supported" in message
    assert "concatenated" in message
    # Refused during planning, so nothing was written and nothing was fetched.
    assert list(tmp_path.iterdir()) == []


def test_extract_fastq_without_sheets_writes_only_diagnostics(tmp_path: Path) -> None:
    key = f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz"
    client = MockS3Client(
        keys=[key],
        sizes={key: 1},
        crc_by_key={key: "crc"},
        object_bodies={f"{key}-metadata.json": json.dumps({"read_count": 5})},
    )
    summary = extract_fastq(
        client,
        BUCKET,
        PREFIX,
        str(tmp_path / "diagnostics.tsv"),
        show_progress=False,
        inline=True,
    )
    assert summary.set_count == 0
    assert list(tmp_path.iterdir()) == [tmp_path / "diagnostics.tsv"]
