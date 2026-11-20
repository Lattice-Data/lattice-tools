"""CLI tests for --resequencing-orders and the outside-raw/ warning (fastq, cram)."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from file_extract import sequence_sets
from file_extract.cli import _outside_raw_lines, main
from file_extract.constants import FASTQ_COLUMNS, SHEET_HELPER_COLUMNS
from file_extract.cram import extract_cram, extract_cram_sets, scan_cram_listing
from file_extract.fastq import (
    FastqListing,
    extract_fastq,
    extract_fastq_sets,
    scan_fastq_listing,
)
from file_extract.models import ListedObject
from tests.file_extract_helpers import MockS3Client

BUCKET = "czi-novogene"
PROJ = "smith-proj/"
MAIN_URI = f"s3://{BUCKET}/{PROJ}ORDER_A"
CONTAINS_NEW = ["--resequencing-contains-new", "true"]
SHEET_ARGS = [
    "--lab",
    "example-lab",
    "--cro-order",
    "ORDER_A",
    "--is-pilot-order",
    "false",
    "--sequencing-platform",
    "Ultima Genomics UG 100",
]


def _fq(order: str, group: str, stem: str) -> list[str]:
    base = f"{PROJ}{order}/{group}/raw/{stem}"
    return [f"{base}_R1_001.fastq.gz", f"{base}_R2_001.fastq.gz"]


def _client(keys: list[str]) -> MockS3Client:
    return MockS3Client(
        keys=keys,
        sizes={k: 10 for k in keys},
        crc_by_key={k: "crc==" for k in keys},
        object_bodies={f"{k}-metadata.json": '{"read_count": 7}' for k in keys},
    )


def _inline(func: Any):
    def run(*args: Any, **kwargs: Any) -> Any:
        kwargs["inline"] = True
        return func(*args, **kwargs)

    return run


def _run(
    keys: list[str],
    argv: list[str],
    tmp_path: Path,
    *,
    mode: str = "fastq",
    client: MockS3Client | None = None,
) -> int:
    """Run main() against a mock bucket, extraction in-process."""
    extra = ["--cram-slot", "trimmed"] if mode == "cram" else []
    with (
        patch("file_extract.cli.boto3.client", return_value=client or _client(keys)),
        patch("file_extract.cli.extract_fastq", _inline(extract_fastq)),
        patch("file_extract.cli.extract_cram", _inline(extract_cram)),
        patch("file_extract.cli.extract_fastq_sets", _inline(extract_fastq_sets)),
        patch("file_extract.cli.extract_cram_sets", _inline(extract_cram_sets)),
    ):
        return main(
            [
                mode,
                MAIN_URI,
                "-o",
                str(tmp_path / f"ORDER_A_{mode}_info.tsv"),
                "--out-dir",
                str(tmp_path),
                "-q",
                *SHEET_ARGS,
                *extra,
                *argv,
            ]
        )


def _rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()]


def _column(path: Path, name: str) -> list[str]:
    rows = _rows(path)
    index = rows[0].index(name)
    return [row[index] for row in rows[1:]]


def _reseq(*orders: str, contains_new: str | None = "true") -> list[str]:
    argv = [
        "--resequencing-orders",
        ",".join(f"s3://{BUCKET}/{PROJ}{o}" for o in orders),
    ]
    if contains_new is not None:
        argv += ["--resequencing-contains-new", contains_new]
    return argv


# --- routing through the CLI ---------------------------------------------


def test_known_group_is_added_and_original_rows_are_unchanged(tmp_path: Path) -> None:
    main_keys = _fq("ORDER_A", "GRP_001", "RUN1-GRP_001")
    reseq_keys = _fq("ORDER_B", "GRP_001", "RUN2-GRP_001")

    single = tmp_path / "single"
    single.mkdir()
    assert _run(main_keys, [], single) == 0
    multi = tmp_path / "multi"
    multi.mkdir()
    assert _run(main_keys + reseq_keys, _reseq("ORDER_B"), multi) == 0

    for sheet in ("ORDER_A_SequenceFile.tsv", "ORDER_A_SequenceFileSet.tsv"):
        before = (single / sheet).read_text().splitlines()
        after = (multi / sheet).read_text().splitlines()
        # Every original row, header included, survives byte for byte, first.
        assert after[: len(before)] == before
    assert _column(multi / "ORDER_A_SequenceFileSet.tsv", "CRO_order") == [
        "ORDER_A",
        "ORDER_B",
    ]
    assert _column(multi / "ORDER_A_SequenceFileSet.tsv", "aliases") == [
        '["example-lab:RUN1-GRP_001"]',
        '["example-lab:RUN2-GRP_001"]',
    ]
    assert not (multi / "ORDER_B_SequenceFile.tsv").exists()


def test_new_group_gets_its_own_output_set(tmp_path: Path) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1-GRP_001") + _fq(
        "ORDER_B", "GRP_003", "RUN2-GRP_003"
    )
    assert _run(keys, _reseq("ORDER_B"), tmp_path) == 0
    assert _column(tmp_path / "ORDER_A_SequenceFileSet.tsv", "aliases") == [
        '["example-lab:RUN1-GRP_001"]'
    ]
    assert _column(tmp_path / "ORDER_B_SequenceFileSet.tsv", "aliases") == [
        '["example-lab:RUN2-GRP_003"]'
    ]
    assert _column(tmp_path / "ORDER_B_SequenceFileSet.tsv", "CRO_order") == ["ORDER_B"]
    assert _column(tmp_path / "ORDER_B_fastq_info.tsv", "source_order") == [
        "ORDER_B",
        "ORDER_B",
    ]


def test_mixed_orders_split_and_two_orders_make_two_extra_sets(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = (
        _fq("ORDER_A", "GRP_001", "RUN1-GRP_001")
        + _fq("ORDER_A", "GRP_002", "RUN1-GRP_002")
        + _fq("ORDER_B", "GRP_001", "RUN2-GRP_001")
        + _fq("ORDER_B", "GRP_003", "RUN2-GRP_003")
        + _fq("ORDER_C", "GRP_002", "RUN3-GRP_002")
        + _fq("ORDER_C", "GRP_004", "RUN3-GRP_004")
    )
    assert _run(keys, _reseq("ORDER_B", "ORDER_C"), tmp_path) == 0
    assert _column(tmp_path / "ORDER_A_fastq_info.tsv", "source_order") == [
        "ORDER_A",
        "ORDER_A",
        "ORDER_A",
        "ORDER_A",
        "ORDER_B",
        "ORDER_B",
        "ORDER_C",
        "ORDER_C",
    ]
    assert _column(tmp_path / "ORDER_B_SequenceFileSet.tsv", "aliases") == [
        '["example-lab:RUN2-GRP_003"]'
    ]
    assert _column(tmp_path / "ORDER_C_SequenceFileSet.tsv", "aliases") == [
        '["example-lab:RUN3-GRP_004"]'
    ]
    out = capsys.readouterr().out
    assert "Added to ORDER_A: 1 GroupID(s), 2 file(s)" in out
    assert "New GroupIDs, written to ORDER_C outputs: 1 GroupID(s), 2 file(s)" in out
    assert "Output set ORDER_A (main order, plus 4 re-sequenced file(s))" in out
    assert "4 file(s); 2 GroupID(s)" in out
    assert "GRP_00" not in out


def test_info_tsv_has_source_order_only_with_the_flag(tmp_path: Path) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1-GRP_001")
    single = tmp_path / "single"
    single.mkdir()
    assert _run(keys, [], single) == 0
    assert _rows(single / "ORDER_A_fastq_info.tsv")[0] == list(FASTQ_COLUMNS) + list(
        SHEET_HELPER_COLUMNS
    )


def test_all_files_are_fetched_in_one_pass(tmp_path: Path) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1-GRP_001") + _fq(
        "ORDER_B", "GRP_003", "RUN2-GRP_003"
    )
    calls: list[list[str]] = []
    real = sequence_sets.fetch_results

    def spy(s3_client: Any, bucket: str, wanted: list[str], **kwargs: Any) -> Any:
        calls.append(list(wanted))
        return real(s3_client, bucket, wanted, **kwargs)

    with patch.object(sequence_sets, "fetch_results", spy):
        assert _run(keys, _reseq("ORDER_B"), tmp_path) == 0
    assert len(calls) == 1
    assert sorted(calls[0]) == sorted(keys)


def test_cram_reseq_routes_and_reports_order(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = [
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1-GRP_001.cram",
        f"{PROJ}ORDER_B/GRP_001/raw/RUN2-GRP_001.cram",
        f"{PROJ}ORDER_B/GRP_003/raw/RUN2-GRP_003.cram",
        f"{PROJ}ORDER_B/GRP_003/raw/RUN2-GRP_003.ucram",
    ]
    assert _run(keys, _reseq("ORDER_B"), tmp_path, mode="cram") == 0
    assert _column(tmp_path / "ORDER_A_SequenceFileSet.tsv", "CRO_order") == [
        "ORDER_A",
        "ORDER_B",
    ]
    assert _column(tmp_path / "ORDER_B_SequenceFileSet.tsv", "trimmed_cram") == [
        "example-lab:RUN2-GRP_003.cram"
    ]
    assert _column(tmp_path / "ORDER_A_cram_info.tsv", "source_order") == [
        "ORDER_A",
        "ORDER_B",
    ]
    assert "WARNING: ORDER_B: Found 1 .ucram file(s)" in capsys.readouterr().out


# --- errors ----------------------------------------------------------------


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["--resequencing-orders", "s3://other/smith-proj/ORDER_B", *CONTAINS_NEW],
            "same bucket",
        ),
        (["--resequencing-orders", f"s3://{BUCKET}", *CONTAINS_NEW], "a bucket alone"),
        (["--resequencing-orders", f"{MAIN_URI}/GRP_001", *CONTAINS_NEW], "overlaps"),
        (["--resequencing-orders", MAIN_URI, *CONTAINS_NEW], "is the main order"),
        (
            _reseq("ORDER_B") + ["--no-require-raw"],
            "--no-require-raw cannot be combined",
        ),
    ],
)
def test_usage_errors_exit_2_before_any_s3_call(
    argv: list[str], message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with patch("file_extract.cli.boto3.client") as boto:
        with pytest.raises(SystemExit) as info:
            main(["fastq", MAIN_URI, *SHEET_ARGS, *argv])
    assert info.value.code == 2
    assert message in capsys.readouterr().err
    boto.assert_not_called()


@pytest.mark.parametrize(
    ("bad_key", "message"),
    [
        (f"{PROJ}ORDER_B/RNA3_098/raw/426971/x_R1_001.fastq.gz", "only the 10x layout"),
        (f"{PROJ}ORDER_B/raw/x_R1_001.fastq.gz", "must be the order folder"),
    ],
)
def test_layout_errors_exit_1_before_any_fetch(
    bad_key: str,
    message: str,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1") + [bad_key]
    client = _client(keys)
    assert _run(keys, _reseq("ORDER_B"), tmp_path, client=client) == 1
    err = capsys.readouterr().err
    assert message in err
    assert f"s3://{BUCKET}/{bad_key}" in err
    assert client.get_object_calls == []
    assert list(tmp_path.iterdir()) == []


def test_file_directly_in_order_folder_is_an_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Only reachable with --no-require-raw, which re-sequencing refuses, but the
    # layout check still names it if a key ever slips through.
    loose = f"{PROJ}ORDER_B/loose_R1_001.fastq.gz"
    keys = _fq("ORDER_A", "GRP_001", "RUN1")
    main_listing = scan_fastq_listing(_client(keys), BUCKET, f"{PROJ}ORDER_A/")

    def scan(_client: Any, _bucket: str, prefix: str, **_kw: Any) -> FastqListing:
        if prefix.endswith("ORDER_B/"):
            return FastqListing(targets=(ListedObject(loose, 1),))
        return main_listing

    with patch("file_extract.cli.scan_fastq_listing", scan):
        assert _run(keys, _reseq("ORDER_B"), tmp_path) == 1
    err = capsys.readouterr().err
    assert "no GroupID folder" in err
    assert loose in err


def test_alias_collision_across_orders_names_both_orders(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = _fq("ORDER_A", "GRP_001", "SAME") + _fq("ORDER_B", "GRP_001", "SAME")
    client = _client(keys)
    assert _run(keys, _reseq("ORDER_B"), tmp_path, client=client) == 1
    err = capsys.readouterr().err
    assert "(ORDER_A)" in err
    assert "(ORDER_B)" in err
    assert client.get_object_calls == []


def test_empty_reseq_order_warns_and_fails_only_under_strict(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1")
    assert _run(keys, _reseq("ORDER_B"), tmp_path) == 0
    assert "ORDER_B" in capsys.readouterr().out
    strict = tmp_path / "strict"
    strict.mkdir()
    assert _run(keys, _reseq("ORDER_B") + ["--strict"], strict) == 1


def test_scale_h5ad_rejects_resequencing_orders(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as info:
        main(
            [
                "scale_h5ad",
                "s3://b/p/o/processed/run/",
                "--lab",
                "example-lab",
                "--metadata-gid",
                "gid",
                "--metadata-experiment",
                "exp",
                "--raw-subdirs",
                "426971",
                "--resequencing-orders",
                "s3://b/p/o2",
            ]
        )
    assert info.value.code == 2
    assert "--resequencing-orders" in capsys.readouterr().err


# --- outside-raw/ warning --------------------------------------------------


def test_fastq_outside_raw_warning_single_order(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    stray = f"{PROJ}ORDER_A/GRP_001/processed/stray_R1_001.fastq.gz"
    keys = _fq("ORDER_A", "GRP_001", "RUN1") + [stray]
    assert _run(keys, ["--strict"], tmp_path) == 0
    out = capsys.readouterr().out
    assert (
        "WARNING: 1 FASTQ file(s) were SKIPPED because they are not in a raw/ folder."
        in out
    )
    assert "These files are NOT in any output:" in out
    assert f"  s3://{BUCKET}/{stray}" in out
    assert "To include them, re-run with --no-require-raw." in out
    assert "1 file(s) outside raw/ were skipped; see above" in out
    for path in tmp_path.iterdir():
        assert "stray" not in path.read_text()


def test_no_require_raw_includes_the_file_and_does_not_warn(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    stray = f"{PROJ}ORDER_A/GRP_001/processed/stray_R1_001.fastq.gz"
    keys = _fq("ORDER_A", "GRP_001", "RUN1") + [stray]
    assert _run(keys, ["--no-require-raw"], tmp_path) == 0
    assert "SKIPPED" not in capsys.readouterr().out
    assert "stray_R1_001.fastq.gz" in (tmp_path / "ORDER_A_fastq_info.tsv").read_text()


def test_outside_raw_warning_with_zero_matches(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    stray = f"{PROJ}ORDER_A/GRP_001/stray_R1_001.fastq.gz"
    assert _run([stray], [], tmp_path) == 0
    out = capsys.readouterr().out
    assert "1 FASTQ file(s) were SKIPPED" in out
    assert "Nothing to do" in out


def test_outside_raw_warning_groups_by_order_without_hint(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = (
        _fq("ORDER_A", "GRP_001", "RUN1")
        + _fq("ORDER_B", "GRP_001", "RUN2")
        + [
            f"{PROJ}ORDER_A/GRP_001/processed/a_R1_001.fastq.gz",
            f"{PROJ}ORDER_B/GRP_001/processed/b_R1_001.fastq.gz",
            f"{PROJ}ORDER_B/GRP_001/processed/b_R2_001.fastq.gz",
        ]
    )
    assert _run(keys, _reseq("ORDER_B"), tmp_path) == 0
    out = capsys.readouterr().out
    assert "WARNING: 3 FASTQ file(s) were SKIPPED" in out
    assert "  ORDER_A (1):\n    s3://" in out
    assert "  ORDER_B (2):\n    s3://" in out
    assert "--no-require-raw" not in out


def test_cram_outside_raw_warning_does_not_double_report_unmatched(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = [
        f"{PROJ}ORDER_A/GRP_001/raw/c1.cram",
        f"{PROJ}ORDER_A/GRP_001/processed/c2.cram",
        f"{PROJ}ORDER_A/GRP_001/processed/c3_unmatched.cram",
    ]
    assert _run(keys, [], tmp_path, mode="cram") == 0
    out = capsys.readouterr().out
    assert "1 CRAM file(s) were SKIPPED" in out
    assert "processed/c2.cram" in out
    assert "c3_unmatched.cram" not in out


def test_scan_listings_collect_outside_raw() -> None:
    keys = [
        f"{PROJ}ORDER_A/G/raw/a_R1_001.fastq.gz",
        f"{PROJ}ORDER_A/G/out/b_R1_001.fastq.gz",
        f"{PROJ}ORDER_A/G/out/b_sample.fastq.gz",
        f"{PROJ}ORDER_A/G/raw/c.cram",
        f"{PROJ}ORDER_A/G/out/d.cram",
        f"{PROJ}ORDER_A/G/out/e_unmatched.cram",
    ]
    client = _client(keys)
    fq = scan_fastq_listing(client, BUCKET, f"{PROJ}ORDER_A/")
    assert [o.key for o in fq.targets] == [keys[0]]
    assert fq.outside_raw == (keys[1],)
    relaxed = scan_fastq_listing(client, BUCKET, f"{PROJ}ORDER_A/", require_raw=False)
    assert relaxed.outside_raw == ()
    cram = scan_cram_listing(client, BUCKET, f"{PROJ}ORDER_A/")
    assert cram.outside_raw == (keys[4],)
    assert cram.warnings.unmatched_cram_count == 1


def test_outside_raw_lines_cap_at_the_print_limit() -> None:
    uris = [f"s3://b/k{i}" for i in range(60)]
    lines = _outside_raw_lines(
        "FASTQ", [("ORDER_A", uris[:30]), ("ORDER_B", uris[30:])], hint=False, limit=50
    )
    assert sum(1 for line in lines if line.startswith("    s3://")) == 50
    assert lines[-1] == "  ... and 10 more"
    assert _outside_raw_lines("FASTQ", [(None, [])], hint=True) == []


def test_single_order_cli_mocks_still_drive_extraction(tmp_path: Path) -> None:
    """The single-order path still calls extract_fastq even on an empty listing."""
    from file_extract.models import RunSummary

    with (
        patch("file_extract.cli.boto3.client", return_value=MockS3Client()),
        patch("file_extract.cli.extract_fastq") as extract,
    ):
        extract.return_value = RunSummary(total=0)
        assert main(["fastq", MAIN_URI, "-q", *SHEET_ARGS]) == 0
    extract.assert_called_once()
    assert extract.call_args.kwargs["listing"].targets == ()


# --- h5 ----------------------------------------------------------------------

H5 = "sample_filtered_feature_bc_matrix.h5"


def _h5(order: str, group: str, sample: str = "S1") -> str:
    return (
        f"{PROJ}{order}/{group}/processed/cellranger/Run_2024-01-01/outs/"
        f"per_sample_outs/{sample}/{H5}"
    )


def _run_h5(keys: list[str], argv: list[str], tmp_path: Path) -> int:
    with patch("file_extract.cli.boto3.client", return_value=_client(keys)):
        return main(
            [
                "h5",
                MAIN_URI,
                "--lab",
                "example-lab",
                "--no-introspect",
                "--no-md5",
                "-o",
                str(tmp_path / "ORDER_A_h5_info.tsv"),
                "-q",
                *argv,
            ]
        )


def _h5_rows(path: Path) -> dict[str, dict[str, str]]:
    # The h5 info TSV uses ordinary CSV quoting, unlike the submission sheets.
    with path.open(newline="", encoding="utf-8") as fh:
        return {row["s3_uri"]: row for row in csv.DictReader(fh, delimiter="\t")}


def test_h5_raw_only_main_group_takes_the_reseq_h5(tmp_path: Path) -> None:
    # ORDER_A holds GRP_001 only as raw/; its h5 was made in ORDER_B.
    keys = [
        _h5("ORDER_A", "GRP_002"),
        f"{PROJ}ORDER_A/GRP_002/raw/RUN1-G2_R1_001.fastq.gz",
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1_R1_001.fastq.gz",
        _h5("ORDER_B", "GRP_001"),
        f"{PROJ}ORDER_B/GRP_001/raw/RUN2_R1_001.fastq.gz",
        _h5("ORDER_B", "GRP_003"),
        f"{PROJ}ORDER_B/GRP_003/raw/RUN2-G3_R1_001.fastq.gz",
    ]
    assert _run_h5(keys, _reseq("ORDER_B"), tmp_path) == 0

    main_rows = _h5_rows(tmp_path / "ORDER_A_h5_info.tsv")
    assert set(main_rows) == {
        f"s3://{BUCKET}/{_h5('ORDER_A', 'GRP_002')}",
        f"s3://{BUCKET}/{_h5('ORDER_B', 'GRP_001')}",
    }
    original = main_rows[f"s3://{BUCKET}/{_h5('ORDER_A', 'GRP_002')}"]
    resequenced = main_rows[f"s3://{BUCKET}/{_h5('ORDER_B', 'GRP_001')}"]
    assert original["source_order"] == "ORDER_A"
    assert resequenced["source_order"] == "ORDER_B"
    # GRP_002 is only in ORDER_A, so its derived_from is ORDER_A's raw/ alone.
    assert original["derived_from"] == '["example-lab:RUN1-G2_R1_001.fastq.gz"]'
    # GRP_001 is in both orders: the re-processed h5 lists both runs' FASTQs.
    assert json.loads(resequenced["derived_from"]) == [
        "example-lab:RUN1_R1_001.fastq.gz",
        "example-lab:RUN2_R1_001.fastq.gz",
    ]

    new_rows = _h5_rows(tmp_path / "ORDER_B_h5_info.tsv")
    assert list(new_rows) == [f"s3://{BUCKET}/{_h5('ORDER_B', 'GRP_003')}"]
    assert new_rows[f"s3://{BUCKET}/{_h5('ORDER_B', 'GRP_003')}"]["library"] == (
        "GRP_003"
    )


@pytest.mark.parametrize(
    "reseq_processed",
    [_h5("ORDER_B", "GRP_001"), f"{PROJ}ORDER_B/GRP_001/processed/other.h5"],
)
def test_h5_shared_group_with_two_processed_folders_errors(
    reseq_processed: str,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    keys = [
        _h5("ORDER_A", "GRP_001"),
        reseq_processed,
        f"{PROJ}ORDER_B/GRP_001/raw/RUN2_R1_001.fastq.gz",
        _h5("ORDER_B", "GRP_003"),
    ]
    assert _run_h5(keys, _reseq("ORDER_B"), tmp_path) == 1
    err = capsys.readouterr().err
    assert "GroupID GRP_001:" in err
    assert f"s3://{BUCKET}/{PROJ}ORDER_A/GRP_001/processed/" in err
    assert f"s3://{BUCKET}/{PROJ}ORDER_B/GRP_001/processed/" in err
    assert "GRP_003" not in err
    # Stopped before any output was written.
    assert list(tmp_path.iterdir()) == []


def test_h5_main_order_h5_lists_reseq_raw_fastqs_too(tmp_path: Path) -> None:
    # The h5 sits in the main order; ORDER_B only added reads for GRP_001.
    keys = [
        _h5("ORDER_A", "GRP_001"),
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1_R1_001.fastq.gz",
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1_R2_001.fastq.gz",
        f"{PROJ}ORDER_B/GRP_001/raw/RUN2_R1_001.fastq.gz",
        f"{PROJ}ORDER_B/GRP_001/raw/RUN2_R2_001.fastq.gz",
    ]
    assert _run_h5(keys, _reseq("ORDER_B"), tmp_path) == 0
    rows = _h5_rows(tmp_path / "ORDER_A_h5_info.tsv")
    row = rows[f"s3://{BUCKET}/{_h5('ORDER_A', 'GRP_001')}"]
    assert row["source_order"] == "ORDER_A"
    assert json.loads(row["derived_from"]) == [
        "example-lab:RUN1_R1_001.fastq.gz",
        "example-lab:RUN1_R2_001.fastq.gz",
        "example-lab:RUN2_R1_001.fastq.gz",
        "example-lab:RUN2_R2_001.fastq.gz",
    ]


def test_h5_single_order_has_no_source_order_column(tmp_path: Path) -> None:
    assert _run_h5([_h5("ORDER_A", "GRP_001")], [], tmp_path) == 0
    header = _rows(tmp_path / "ORDER_A_h5_info.tsv")[0]
    assert "source_order" not in header


@pytest.mark.parametrize(
    ("main_uri_suffix", "bad_key", "message"),
    [
        ("", f"{PROJ}ORDER_B/GRP_003/other/{H5}", "only the 10x layout"),
        ("/GRP_001", f"{PROJ}ORDER_B/GRP_003/raw/x", "must be the order folder"),
    ],
)
def test_h5_layout_errors(
    main_uri_suffix: str,
    bad_key: str,
    message: str,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    keys = [_h5("ORDER_A", "GRP_001"), bad_key]
    with patch("file_extract.cli.boto3.client", return_value=_client(keys)):
        code = main(
            [
                "h5",
                MAIN_URI + main_uri_suffix,
                "--lab",
                "example-lab",
                "--no-introspect",
                "--no-md5",
                "-o",
                str(tmp_path / "out.tsv"),
                "-q",
                *_reseq("ORDER_B"),
            ]
        )
    assert code == 1
    assert message in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


# --- --resequencing-contains-new ---------------------------------------------


def test_contains_new_false_writes_no_output_for_new_groups(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = (
        _fq("ORDER_A", "GRP_001", "RUN1-GRP_001")
        + _fq("ORDER_B", "GRP_001", "RUN2-GRP_001")
        + _fq("ORDER_B", "GRP_003", "RUN2-GRP_003")
    )
    calls: list[list[str]] = []
    real = sequence_sets.fetch_results

    def spy(s3_client: Any, bucket: str, wanted: list[str], **kwargs: Any) -> Any:
        calls.append(list(wanted))
        return real(s3_client, bucket, wanted, **kwargs)

    with patch.object(sequence_sets, "fetch_results", spy):
        code = _run(keys, _reseq("ORDER_B", contains_new="false"), tmp_path)
    assert code == 0
    # The shared GroupID is still added to the main order.
    assert _column(tmp_path / "ORDER_A_fastq_info.tsv", "source_order") == [
        "ORDER_A",
        "ORDER_A",
        "ORDER_B",
        "ORDER_B",
    ]
    # No ORDER_B outputs, and the left-out files were never fetched.
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "ORDER_A_SequenceFile.tsv",
        "ORDER_A_SequenceFileSet.tsv",
        "ORDER_A_fastq_info.tsv",
    ]
    assert len(calls) == 1
    assert not any("GRP_003" in key for key in calls[0])
    out = capsys.readouterr().out
    assert "New GroupIDs, NOT written" in out
    assert "1 GroupID(s), 2 file(s)" in out
    assert "GRP_003" not in out
    assert "Left out 2 file(s)" in out


def test_contains_new_defaults_to_false(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    keys = _fq("ORDER_A", "GRP_001", "RUN1-GRP_001") + _fq(
        "ORDER_B", "GRP_003", "RUN2-GRP_003"
    )
    assert _run(keys, _reseq("ORDER_B", contains_new=None), tmp_path) == 0
    assert not (tmp_path / "ORDER_B_fastq_info.tsv").exists()
    assert "New GroupIDs, NOT written" in capsys.readouterr().out


def test_contains_new_false_h5_writes_only_the_main_tsv(tmp_path: Path) -> None:
    keys = [
        _h5("ORDER_A", "GRP_002"),
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1_R1_001.fastq.gz",
        _h5("ORDER_B", "GRP_001"),
        _h5("ORDER_B", "GRP_003"),
    ]
    assert _run_h5(keys, _reseq("ORDER_B", contains_new="false"), tmp_path) == 0
    assert [p.name for p in tmp_path.iterdir()] == ["ORDER_A_h5_info.tsv"]
    assert set(_h5_rows(tmp_path / "ORDER_A_h5_info.tsv")) == {
        f"s3://{BUCKET}/{_h5('ORDER_A', 'GRP_002')}",
        f"s3://{BUCKET}/{_h5('ORDER_B', 'GRP_001')}",
    }


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            CONTAINS_NEW,
            "--resequencing-contains-new needs --resequencing-orders",
        ),
        (_reseq("ORDER_B", contains_new="maybe"), "expected 'true' or 'false'"),
    ],
)
def test_contains_new_usage_errors_exit_2(
    argv: list[str], message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with patch("file_extract.cli.boto3.client") as boto:
        with pytest.raises(SystemExit) as info:
            main(["fastq", MAIN_URI, *SHEET_ARGS, *argv])
    assert info.value.code == 2
    assert message in capsys.readouterr().err
    boto.assert_not_called()


def test_derived_from_cells_combines_raw_prefixes_and_warns_once_if_empty() -> None:
    from file_extract.h5 import derived_from_cells
    from file_extract.models import RunSummary

    keys = [
        f"{PROJ}ORDER_A/GRP_001/raw/RUN1_R1_001.fastq.gz",
        f"{PROJ}ORDER_B/GRP_001/raw/RUN1_R1_001.fastq.gz",
    ]
    h5_shared = _h5("ORDER_B", "GRP_001")
    h5_empty = _h5("ORDER_B", "GRP_009")
    summary = RunSummary(total=2)
    cells = derived_from_cells(
        _client(keys),
        BUCKET,
        [h5_shared, h5_empty],
        "example-lab",
        summary,
        raw_prefixes_by_key={
            h5_shared: [f"{PROJ}ORDER_A/GRP_001/raw/", f"{PROJ}ORDER_B/GRP_001/raw/"],
            h5_empty: [f"{PROJ}ORDER_A/GRP_009/raw/", f"{PROJ}ORDER_B/GRP_009/raw/"],
        },
    )
    # The same filename in both orders is one alias.
    assert json.loads(cells[h5_shared]) == ["example-lab:RUN1_R1_001.fastq.gz"]
    assert json.loads(cells[h5_empty]) == []
    assert len(summary.warnings) == 1
    assert "GRP_009/raw/" in summary.warnings[0]
