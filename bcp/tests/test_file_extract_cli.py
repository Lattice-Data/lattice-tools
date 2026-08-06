"""CLI tests for file_extract."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from file_extract.cli import build_parser, main
from tests.file_extract_helpers import MockS3Client

PREFIX = "test-order/NVUS0000000000-01/"
BUCKET = "example-bucket"
H5_PREFIX = "proj/lib/processed/outs/per_sample_outs/"

# Submission-sheet inputs the fastq and cram subcommands require.
SHEET_ARGS = [
    "--lab",
    "example-lab",
    "--cro-order",
    "NVUS0000000000-01",
    "--is-pilot-order",
    "false",
    "--sequencing-platform",
    "Ultima Genomics UG 100",
]
CRAM_SHEET_ARGS = [*SHEET_ARGS, "--cram-slot", "trimmed"]


def test_build_parser_requires_subcommand() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "fastq" in out
    assert "cram" in out
    assert "h5" in out


def test_cli_invalid_uri() -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["fastq", "not-a-uri"])
    assert exc_info.value.code == 2


@patch("file_extract.cli.extract_fastq")
@patch("file_extract.cli.boto3.client")
def test_cli_fastq_success(
    mock_boto: MagicMock, mock_extract: MagicMock, tmp_path: Path
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(
        total=2,
        crc_ok=2,
        enrichment_ok=2,
        read_tally={"R1": 1, "R2": 1},
    )
    out = tmp_path / "fastq.tsv"

    code = main(
        [
            "fastq",
            f"s3://{BUCKET}/{PREFIX}",
            "-o",
            str(out),
            "--quiet",
            *SHEET_ARGS,
        ]
    )
    assert code == 0


@patch("file_extract.cli.extract_fastq")
@patch("file_extract.cli.boto3.client")
def test_cli_fastq_strict_on_failure(
    mock_boto: MagicMock, mock_extract: MagicMock, tmp_path: Path
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(
        total=1,
        crc_ok=1,
        enrichment_ok=0,
        failures=[("key", "", "metadata error")],
    )
    out = tmp_path / "fastq_strict.tsv"

    code = main(
        [
            "fastq",
            f"s3://{BUCKET}/{PREFIX}",
            "-o",
            str(out),
            "--quiet",
            "--strict",
            *SHEET_ARGS,
        ]
    )
    assert code == 1


@patch("file_extract.cli.extract_h5")
@patch("file_extract.cli.boto3.client")
def test_cli_h5_checksum_only(
    mock_boto: MagicMock, mock_extract: MagicMock, tmp_path: Path
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(total=1, crc_ok=1)
    out = tmp_path / "h5.tsv"

    code = main(
        [
            "h5",
            f"s3://{BUCKET}/{H5_PREFIX}",
            "-o",
            str(out),
            "--no-introspect",
            "--quiet",
        ]
    )
    assert code == 0


@patch("file_extract.cli.extract_fastq")
@patch("file_extract.cli.boto3.client")
def test_cli_fastq_zero_matches(
    mock_boto: MagicMock,
    mock_extract: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(total=0)
    code = main(["fastq", f"s3://{BUCKET}/{PREFIX}", "--quiet", *SHEET_ARGS])
    assert code == 0
    assert "Nothing to do" in capsys.readouterr().out


@patch("file_extract.cli.extract_cram")
@patch("file_extract.cli.scan_cram_listing_warnings")
@patch("file_extract.cli.boto3.client")
def test_cli_cram_success(
    mock_boto: MagicMock,
    mock_scan: MagicMock,
    mock_extract: MagicMock,
    tmp_path: Path,
) -> None:
    from file_extract.cram import CramListingWarnings
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_scan.return_value = CramListingWarnings()
    mock_extract.return_value = RunSummary(total=1, crc_ok=1, enrichment_ok=1)
    out = tmp_path / "cram.tsv"

    code = main(
        [
            "cram",
            f"s3://{BUCKET}/{PREFIX}",
            "-o",
            str(out),
            "--quiet",
            *CRAM_SHEET_ARGS,
        ]
    )
    assert code == 0


@patch("file_extract.cli.extract_cram")
@patch("file_extract.cli.scan_cram_listing_warnings")
@patch("file_extract.cli.boto3.client")
def test_cli_cram_zero_matches_only_unmatched(
    mock_boto: MagicMock,
    mock_scan: MagicMock,
    mock_extract: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from file_extract.cram import CramListingWarnings
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_scan.return_value = CramListingWarnings(unmatched_cram_count=2)
    mock_extract.return_value = RunSummary(total=0)
    code = main(["cram", f"s3://{BUCKET}/{PREFIX}", "--quiet", *CRAM_SHEET_ARGS])
    assert code == 0
    out = capsys.readouterr().out
    assert "Nothing to do" in out
    assert "2 unmatched" in out


@patch("file_extract.cli.extract_cram")
@patch("file_extract.cli.scan_cram_listing_warnings")
@patch("file_extract.cli.boto3.client")
def test_cli_cram_ucram_warning(
    mock_boto: MagicMock,
    mock_scan: MagicMock,
    mock_extract: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from file_extract.cram import CramListingWarnings
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_scan.return_value = CramListingWarnings(ucram_count=1)
    mock_extract.return_value = RunSummary(total=0)
    code = main(["cram", f"s3://{BUCKET}/{PREFIX}", "--quiet", *CRAM_SHEET_ARGS])
    assert code == 0
    out = capsys.readouterr().out
    assert ".ucram" in out


@pytest.mark.parametrize(
    "argv",
    [
        # Every sheet input is mandatory; dropping any one is a usage error.
        ["fastq", f"s3://{BUCKET}/{PREFIX}"],
        ["fastq", f"s3://{BUCKET}/{PREFIX}", "--lab", "example-lab"],
        ["cram", f"s3://{BUCKET}/{PREFIX}", *SHEET_ARGS],
        ["fastq", f"s3://{BUCKET}/{PREFIX}", *SHEET_ARGS, "--is-pilot-order", "maybe"],
    ],
)
def test_cli_rejects_incomplete_sheet_arguments(argv: list[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(argv)
    assert exc_info.value.code == 2


@pytest.mark.parametrize(
    ("value", "expected"),
    [("true", "TRUE"), ("TRUE", "TRUE"), ("false", "FALSE"), ("False", "FALSE")],
)
def test_cli_parses_pilot_flag_case_insensitively(value: str, expected: str) -> None:
    args = build_parser().parse_args(
        [
            "fastq",
            f"s3://{BUCKET}/{PREFIX}",
            "--lab",
            "example-lab",
            "--cro-order",
            "NVUS0000000000-01",
            "--sequencing-platform",
            "Ultima Genomics UG 100",
            "--is-pilot-order",
            value,
        ]
    )
    assert args.is_pilot_order is (expected == "TRUE")


def _inline_fastq(*args: object, **kwargs: object) -> object:
    """Run the real extractor in-process; the mock S3 client cannot cross a pool."""
    from file_extract.fastq import extract_fastq

    kwargs["inline"] = True
    return extract_fastq(*args, **kwargs)  # type: ignore[arg-type]


def _fastq_client() -> MockS3Client:
    keys = [
        f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz",
        f"{PREFIX}S1/raw/s_L001_R2_001.fastq.gz",
    ]
    return MockS3Client(
        keys=keys,
        sizes={key: 10 for key in keys},
        crc_by_key={key: "crc" for key in keys},
        object_bodies={f"{key}-metadata.json": '{"read_count": 5}' for key in keys},
    )


@patch("file_extract.cli.boto3.client")
def test_cli_fastq_writes_both_sheets(mock_boto: MagicMock, tmp_path: Path) -> None:
    mock_boto.return_value = _fastq_client()
    with patch("file_extract.cli.extract_fastq", _inline_fastq):
        code = main(
            [
                "fastq",
                f"s3://{BUCKET}/{PREFIX}",
                "-o",
                str(tmp_path / "diagnostics.tsv"),
                "--quiet",
                *SHEET_ARGS,
            ]
        )
    assert code == 0
    sequence_file = tmp_path / "NVUS0000000000-01_SequenceFile.tsv"
    file_set = tmp_path / "NVUS0000000000-01_SequenceFileSet.tsv"
    assert sequence_file.exists()
    assert file_set.exists()
    assert len(file_set.read_text(encoding="utf-8").splitlines()) == 2
    assert "example-lab:s_L001_R1_001.fastq.gz" in sequence_file.read_text(
        encoding="utf-8"
    )


@patch("file_extract.cli.boto3.client")
def test_cli_out_dir_receives_sheets(mock_boto: MagicMock, tmp_path: Path) -> None:
    mock_boto.return_value = _fastq_client()
    sheets_dir = tmp_path / "sheets"
    with patch("file_extract.cli.extract_fastq", _inline_fastq):
        code = main(
            [
                "fastq",
                f"s3://{BUCKET}/{PREFIX}",
                "-o",
                str(tmp_path / "diagnostics.tsv"),
                "--out-dir",
                str(sheets_dir),
                "--quiet",
                *SHEET_ARGS,
            ]
        )
    assert code == 0
    assert (sheets_dir / "NVUS0000000000-01_SequenceFile.tsv").exists()
    assert (sheets_dir / "NVUS0000000000-01_SequenceFileSet.tsv").exists()


@patch("file_extract.cli.boto3.client")
def test_cli_fastq_alias_collision_exits_one(
    mock_boto: MagicMock,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    keys = [
        f"{PREFIX}S1/raw/s_L001_R1_001.fastq.gz",
        f"{PREFIX}S2/raw/s_L001_R1_001.fastq.gz",
    ]
    mock_boto.return_value = MockS3Client(
        keys=keys,
        sizes={key: 10 for key in keys},
        crc_by_key={key: "crc" for key in keys},
    )
    with patch("file_extract.cli.extract_fastq", _inline_fastq):
        code = main(
            [
                "fastq",
                f"s3://{BUCKET}/{PREFIX}",
                "-o",
                str(tmp_path / "diagnostics.tsv"),
                "--quiet",
                *SHEET_ARGS,
            ]
        )
    assert code == 1
    assert "Alias collisions" in capsys.readouterr().err


@patch("file_extract.cli.boto3.client")
def test_cli_bad_lab_exits_one(
    mock_boto: MagicMock, capsys: pytest.CaptureFixture[str]
) -> None:
    mock_boto.return_value = MockS3Client()
    code = main(
        [
            "fastq",
            f"s3://{BUCKET}/{PREFIX}",
            "--quiet",
            "--lab",
            "not a lab",
            "--cro-order",
            "NVUS0000000000-01",
            "--is-pilot-order",
            "false",
            "--sequencing-platform",
            "Ultima Genomics UG 100",
        ]
    )
    assert code == 1
    assert "Invalid --lab" in capsys.readouterr().err


@patch("file_extract.cli.extract_fastq")
@patch("file_extract.cli.boto3.client")
def test_cli_warns_when_cro_order_absent_from_prefix(
    mock_boto: MagicMock,
    mock_extract: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(total=0)
    code = main(
        [
            "fastq",
            f"s3://{BUCKET}/{PREFIX}",
            "--quiet",
            "--lab",
            "example-lab",
            "--cro-order",
            "NVUS0000000000-99",
            "--is-pilot-order",
            "false",
            "--sequencing-platform",
            "Ultima Genomics UG 100",
        ]
    )
    assert code == 0
    assert "does not appear in the S3 prefix" in capsys.readouterr().out


@patch("file_extract.cli.extract_fastq")
@patch("file_extract.cli.boto3.client")
def test_cli_reports_sheet_warnings(
    mock_boto: MagicMock,
    mock_extract: MagicMock,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_extract.return_value = RunSummary(
        total=1,
        crc_ok=1,
        enrichment_ok=1,
        set_count=1,
        warnings=["s_L001: R2 present without R1; read1 will be empty"],
    )
    code = main(["fastq", f"s3://{BUCKET}/{PREFIX}", "--quiet", *SHEET_ARGS])
    assert code == 0
    out = capsys.readouterr().out
    assert "R2 present without R1" in out
    assert "SequenceFileSets: 1" in out


@patch("file_extract.cli.extract_cram")
@patch("file_extract.cli.scan_cram_listing_warnings")
@patch("file_extract.cli.boto3.client")
def test_cli_cram_strict_on_failure(
    mock_boto: MagicMock,
    mock_scan: MagicMock,
    mock_extract: MagicMock,
    tmp_path: Path,
) -> None:
    from file_extract.cram import CramListingWarnings
    from file_extract.models import RunSummary

    mock_boto.return_value = MockS3Client()
    mock_scan.return_value = CramListingWarnings()
    mock_extract.return_value = RunSummary(
        total=1,
        crc_ok=1,
        enrichment_ok=0,
        failures=[("key", "", "metadata error")],
    )
    out = tmp_path / "cram_strict.tsv"

    code = main(
        [
            "cram",
            f"s3://{BUCKET}/{PREFIX}",
            "-o",
            str(out),
            "--quiet",
            "--strict",
            *CRAM_SHEET_ARGS,
        ]
    )
    assert code == 1
