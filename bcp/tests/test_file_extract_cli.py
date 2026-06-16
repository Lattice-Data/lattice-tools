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


def test_build_parser_requires_subcommand() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def test_cli_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "fastq" in out
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
    code = main(["fastq", f"s3://{BUCKET}/{PREFIX}", "--quiet"])
    assert code == 0
    assert "Nothing to do" in capsys.readouterr().out
