"""Tests for the guidesig --h5-uri guide_id check.

S3 is never contacted: listings come from file_extract's MockS3Client, and each
"S3" h5 key is served from a small local h5 built in pytest's temp folder.
"""

from __future__ import annotations

from pathlib import Path

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("botocore")

from guidesig import h5check  # noqa: E402
from guidesig.cli import main  # noqa: E402
from guidesig.h5check import (  # noqa: E402
    H5CheckError,
    check_guide_ids_against_h5,
    find_h5_layout,
    guide_ids_from_h5,
    read_sample_names,
)
from tests.file_extract_helpers import MockS3Client  # noqa: E402
from tests.guidesig_helpers import write_delimited  # noqa: E402

BUCKET = "bucket"
URI = f"s3://{BUCKET}/proj/order"
H5_NAME = "sample_filtered_feature_bc_matrix.h5"
OUTS = "proj/order/LIB1/processed/run1/outs/per_sample_outs"
KEY_S1 = f"{OUTS}/S1/count/{H5_NAME}"
KEY_S2 = f"{OUTS}/S2/count/{H5_NAME}"
CRISPR = "CRISPR Guide Capture"
GEX = "Gene Expression"


def _write_h5(path: Path, features: list[tuple[str, str]]) -> Path:
    """Write a minimal Cell Ranger-shaped h5 holding ``(id, feature_type)``."""
    with h5py.File(path, "w") as h5:
        matrix = h5.create_group("matrix")
        matrix.create_dataset("barcodes", data=[b"AAAC-1"])
        feats = matrix.create_group("features")
        feats.create_dataset("id", data=[fid.encode() for fid, _ in features])
        feats.create_dataset(
            "feature_type", data=[ftype.encode() for _, ftype in features]
        )
    return path


def _guides(*ids: str) -> list[tuple[str, str]]:
    return [(gid, CRISPR) for gid in ids]


def _local_reader(files: dict[str, Path]):
    """Serve S3 keys from local h5 files, as read_s3_h5_guide_ids would."""

    def reader(bucket: str, key: str) -> set[str]:
        assert bucket == BUCKET
        with h5py.File(files[key], "r") as h5:
            return guide_ids_from_h5(h5)

    return reader


def _sheet(tmp_path: Path, ids: list[str], *, id_column: str = "guide_id") -> Path:
    header = [
        "property",
        id_column,
        "guide_protospacer",
        "guide_role",
        "guide_target_gene_id",
    ]
    seqs = ["AAAA", "TTTA", "ATAT", "TATA"]
    rows = [
        ["", gid, seqs[i], "Targeting a Gene", "ENSGFAKE0000001"]
        for i, gid in enumerate(ids)
    ]
    return write_delimited(tmp_path / "sheet.tsv", [header, *rows], file_format="tsv")


# --- reading one h5 -----------------------------------------------------------


def test_guide_ids_from_h5_keeps_only_crispr_features(tmp_path: Path) -> None:
    path = _write_h5(
        tmp_path / "m.h5",
        [("ENSG1", GEX), ("g1", CRISPR), ("ab1", "Antibody Capture"), ("g2", CRISPR)],
    )
    with h5py.File(path, "r") as h5:
        assert guide_ids_from_h5(h5) == {"g1", "g2"}


def test_guide_ids_from_h5_requires_matrix_group(tmp_path: Path) -> None:
    path = tmp_path / "empty.h5"
    with h5py.File(path, "w"):
        pass
    with h5py.File(path, "r") as h5, pytest.raises(RuntimeError, match="'matrix'"):
        guide_ids_from_h5(h5)


def test_guide_ids_from_h5_requires_feature_ids(tmp_path: Path) -> None:
    path = tmp_path / "no_ids.h5"
    with h5py.File(path, "w") as h5:
        h5.create_group("matrix/features").create_dataset(
            "feature_type", data=[CRISPR.encode()]
        )
    with h5py.File(path, "r") as h5, pytest.raises(RuntimeError, match="features/id"):
        guide_ids_from_h5(h5)


# --- finding h5 files and folders that lack one ------------------------------


def test_find_h5_layout_flags_sample_and_library_folders_without_h5() -> None:
    keys = [
        KEY_S1,
        f"{OUTS}/S1/metrics_summary.csv",
        f"{OUTS}/S2/metrics_summary.csv",  # S2 has no h5
        "proj/order/LIB1/raw/x_R1_001.fastq.gz",
        "proj/order/LIB2/processed/run1/outs/web_summary.html",  # LIB2 has no h5
        "proj/order/LIB2/raw/y_R1_001.fastq.gz",
    ]
    layout = find_h5_layout(keys, H5_NAME)
    assert layout.h5_keys == [KEY_S1]
    assert layout.sample_dirs_without_h5 == [f"{OUTS}/S2/"]
    assert layout.library_dirs_without_h5 == ["proj/order/LIB2/"]


def test_find_h5_layout_matches_basename_exactly() -> None:
    keys = [f"{OUTS}/S1/count/raw_{H5_NAME}", KEY_S2]
    layout = find_h5_layout(keys, H5_NAME)
    assert layout.h5_keys == [KEY_S2]
    assert layout.sample_dirs_without_h5 == [f"{OUTS}/S1/"]


# --- checking guide_ids against every h5 -------------------------------------


def test_all_samples_pass_and_extras_are_notes(tmp_path: Path) -> None:
    files = {
        KEY_S1: _write_h5(tmp_path / "s1.h5", _guides("g1", "g2", "g9")),
        KEY_S2: _write_h5(tmp_path / "s2.h5", _guides("g1", "g2")),
    }
    report = check_guide_ids_against_h5(
        {"g1", "g2"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.ok
    assert (report.h5_total, report.h5_failed) == (2, 0)
    assert len(report.notes) == 1
    assert "sample S1 (library LIB1)" in report.notes[0]
    assert "1 guide(s) in h5 not in metadata; examples: g9" in report.notes[0]


def test_every_failing_sample_is_reported(tmp_path: Path) -> None:
    """One bad sample must not stop the others from being checked."""
    files = {
        KEY_S1: _write_h5(tmp_path / "s1.h5", _guides("g1")),
        KEY_S2: _write_h5(tmp_path / "s2.h5", _guides("g2")),
    }
    report = check_guide_ids_against_h5(
        {"g1", "g2", "g3"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert not report.ok
    assert (report.h5_total, report.h5_failed) == (2, 2)
    assert report.failures == [
        f"sample S1 (library LIB1) (s3://{BUCKET}/{KEY_S1}): 2 of 3 guide_ids "
        "missing from h5; examples: g2, g3",
        f"sample S2 (library LIB1) (s3://{BUCKET}/{KEY_S2}): 2 of 3 guide_ids "
        "missing from h5; examples: g1, g3",
    ]


def test_unreadable_h5_fails_that_sample_only(tmp_path: Path) -> None:
    files = {KEY_S2: _write_h5(tmp_path / "s2.h5", _guides("g1"))}
    good = _local_reader(files)

    def reader(bucket: str, key: str) -> set[str]:
        if key == KEY_S1:
            raise RuntimeError("truncated file")
        return good(bucket, key)

    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=[KEY_S1, KEY_S2]),
        reader=reader,
    )
    assert (report.h5_total, report.h5_failed) == (2, 1)
    assert report.failures == [
        f"sample S1 (library LIB1) (s3://{BUCKET}/{KEY_S1}): could not read h5: "
        "truncated file"
    ]


def test_h5_without_crispr_features_fails(tmp_path: Path) -> None:
    files = {KEY_S1: _write_h5(tmp_path / "s1.h5", [("ENSG1", GEX)])}
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"sample S1 (library LIB1) (s3://{BUCKET}/{KEY_S1}): "
        "no 'CRISPR Guide Capture' features in h5"
    ]


def test_missing_folders_are_failures(tmp_path: Path) -> None:
    files = {KEY_S1: _write_h5(tmp_path / "s1.h5", _guides("g1"))}
    keys = [*files, f"{OUTS}/S2/metrics_summary.csv"]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert (report.h5_total, report.h5_failed) == (1, 0)
    assert report.failures == [
        f"sample folder s3://{BUCKET}/{OUTS}/S2/ has no {H5_NAME}"
    ]


def test_no_h5_found_is_a_failure() -> None:
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=["proj/order/LIB1/raw/x_R1_001.fastq.gz"]),
        reader=lambda bucket, key: set(),
    )
    assert report.failures == [f"no {H5_NAME} found under {URI}"]


def test_custom_filename_is_honored(tmp_path: Path) -> None:
    raw_key = "proj/order/LIB1/processed/run1/outs/multi/count/raw.h5"
    files = {raw_key: _write_h5(tmp_path / "raw.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        filename="raw.h5",
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.ok
    assert report.h5_total == 1


def test_invalid_uri_raises() -> None:
    with pytest.raises(H5CheckError, match="Invalid S3 URI"):
        check_guide_ids_against_h5(
            {"g1"},
            "bucket/proj/order",
            s3_client=MockS3Client(),
            reader=lambda b, k: set(),
        )


def test_listing_failure_raises() -> None:
    class BrokenClient:
        def get_paginator(self, operation: str):
            raise RuntimeError("access denied")

    with pytest.raises(H5CheckError, match=f"could not list {URI}: access denied"):
        check_guide_ids_against_h5(
            {"g1"}, URI, s3_client=BrokenClient(), reader=lambda b, k: set()
        )


# --- CLI ----------------------------------------------------------------------


@pytest.fixture
def fake_s3(tmp_path: Path, monkeypatch):
    """Route the CLI's default S3 client and h5 reader to local stand-ins."""
    files: dict[str, Path] = {}
    state = {"client_made": False}

    def make_client(pool_size: int) -> MockS3Client:
        state["client_made"] = True
        return MockS3Client(keys=list(files))

    reader = _local_reader(files)
    monkeypatch.setattr(h5check, "check_dependencies", lambda: None)
    monkeypatch.setattr(h5check, "_make_s3_client", make_client)
    monkeypatch.setattr(
        h5check,
        "read_s3_h5_guide_ids",
        lambda bucket, key, *, pool_size: reader(bucket, key),
    )
    return files, state


def test_cli_pass_prints_signature_and_notes(tmp_path: Path, fake_s3, capsys) -> None:
    files, _state = fake_s3
    files[KEY_S1] = _write_h5(tmp_path / "s1.h5", _guides("g1", "g2", "g9"))
    sheet = _sheet(tmp_path, ["g1", "g2"])
    main(["--format", "tsv", "--h5-uri", URI, str(sheet)])
    captured = capsys.readouterr()
    sig, path = captured.out.strip().split("\t")
    assert sig.startswith("gsig1:set:n=2:")
    assert path == str(sheet)
    assert "note: sample S1 (library LIB1)" in captured.err
    assert "error:" not in captured.err


def test_cli_failure_lists_every_problem(tmp_path: Path, fake_s3, capsys) -> None:
    files, _state = fake_s3
    files[KEY_S1] = _write_h5(tmp_path / "s1.h5", _guides("g1"))
    files[KEY_S2] = _write_h5(tmp_path / "s2.h5", _guides("g2"))
    sheet = _sheet(tmp_path, ["g1", "g2"])
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--h5-uri", URI, str(sheet)])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "error: sample S1 (library LIB1)" in captured.err
    assert "error: sample S2 (library LIB1)" in captured.err
    assert (
        "error: guide_id check failed: 2 problem(s); 2 of 2 h5 file(s) failed"
        in captured.err
    )


def test_cli_sheet_failure_skips_s3(tmp_path: Path, fake_s3, capsys) -> None:
    _files, state = fake_s3
    sheet = _sheet(tmp_path, ["g1", "g2"], id_column="not_guide_id")
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--h5-uri", URI, str(sheet)])
    assert exc.value.code == 1
    assert "required column 'guide_id' absent from header" in capsys.readouterr().err
    assert state["client_made"] is False


@pytest.mark.parametrize(
    "extra_args, message",
    [
        (["--compare", "a.tsv", "b.tsv"], "not allowed with --compare"),
        (["a.tsv", "b.tsv"], "exactly one FILE"),
    ],
)
def test_cli_h5_usage_errors(capsys, extra_args: list[str], message: str) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--h5-uri", URI, *extra_args])
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


@pytest.mark.parametrize(
    "flag, value",
    [
        ("--workers", "4"),
        ("--retries", "2"),
        ("--h5-filename", "x"),
        ("--h5-samples-file", "samples.tsv"),
    ],
)
def test_cli_h5_flags_require_h5_uri(capsys, flag: str, value: str) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", flag, value, "a.tsv"])
    assert exc.value.code == 2
    assert "require --h5-uri" in capsys.readouterr().err


def test_cli_workers_must_be_positive(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--h5-uri", URI, "--workers", "0", "a.tsv"])
    assert exc.value.code == 2
    assert "expected a positive integer" in capsys.readouterr().err


# --- --h5-samples-file ----------------------------------------------------------

ES_OUTS = "proj/order/ES_m6_1/processed/run1/outs/per_sample_outs"
NES_OUTS = "proj/order/NES_m6_1/processed/run1/outs/per_sample_outs"
KEY_ES = f"{ES_OUTS}/ES_m6_1/count/{H5_NAME}"
KEY_NES = f"{NES_OUTS}/NES_m6_1/count/{H5_NAME}"


def test_read_sample_names_skips_blanks_comments_and_duplicates(
    tmp_path: Path,
) -> None:
    path = tmp_path / "samples.tsv"
    path.write_text("# ES samples\nES_m6_1\n\n  ES_m8_1  \nES_m6_1\n")
    assert read_sample_names(path) == ["ES_m6_1", "ES_m8_1"]


@pytest.mark.parametrize(
    "content, message",
    [
        ("ES_m6_1\nES_m8_1\tlibrary\n", "line 2: holds a tab"),
        ("# only a comment\n\n", "no sample names found"),
        ("", "no sample names found"),
    ],
)
def test_read_sample_names_rejects(tmp_path: Path, content: str, message: str) -> None:
    path = tmp_path / "samples.tsv"
    path.write_text(content)
    with pytest.raises(H5CheckError, match=message):
        read_sample_names(path)


def test_read_sample_names_unreadable(tmp_path: Path) -> None:
    with pytest.raises(H5CheckError, match="cannot read sample list"):
        read_sample_names(tmp_path / "missing.tsv")


def test_sample_list_reads_only_listed_h5(tmp_path: Path) -> None:
    """ES_m6_1 is an exact name: it must not pull in NES_m6_1."""
    files = {
        KEY_ES: _write_h5(tmp_path / "es.h5", _guides("g1")),
        KEY_NES: _write_h5(tmp_path / "nes.h5", _guides("n1")),
    }
    read: list[str] = []
    local = _local_reader(files)

    def reader(bucket: str, key: str) -> set[str]:
        read.append(key)
        return local(bucket, key)

    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["ES_m6_1"],
        samples_label="samples.tsv",
        s3_client=MockS3Client(keys=list(files)),
        reader=reader,
    )
    assert report.ok
    assert read == [KEY_ES]
    assert report.notes == [
        f"checked 1 of 2 {H5_NAME} file(s); samples from samples.tsv"
    ]


def test_listed_sample_without_h5_fails(tmp_path: Path) -> None:
    files = {KEY_ES: _write_h5(tmp_path / "es.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["ES_m6_1", "ES_m6_3"],
        samples_label="samples.tsv",
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"sample 'ES_m6_3' listed in samples.tsv has no {H5_NAME} under {URI}"
    ]
    assert report.h5_total == 1


def test_sample_list_matching_nothing_reports_each_name() -> None:
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["ES_m6_1"],
        samples_label="samples.tsv",
        s3_client=MockS3Client(keys=[KEY_NES]),
        reader=lambda bucket, key: set(),
    )
    assert report.failures == [
        f"sample 'ES_m6_1' listed in samples.tsv has no {H5_NAME} under {URI}"
    ]
    assert report.notes == [
        f"checked 0 of 1 {H5_NAME} file(s); samples from samples.tsv"
    ]


def test_missing_folder_checks_respect_sample_list(tmp_path: Path) -> None:
    """Only listed samples' folders are checked for a missing h5."""
    files = {KEY_ES: _write_h5(tmp_path / "es.h5", _guides("g1"))}
    keys = [
        *files,
        f"{NES_OUTS}/NES_m6_1/metrics_summary.csv",  # not listed, no h5
        "proj/order/NES_m9_1/processed/run1/outs/web_summary.html",
    ]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["ES_m6_1"],
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert report.ok


def test_cli_h5_samples_file(tmp_path: Path, fake_s3, capsys) -> None:
    files, _state = fake_s3
    files[KEY_ES] = _write_h5(tmp_path / "es.h5", _guides("g1", "g2"))
    files[KEY_NES] = _write_h5(tmp_path / "nes.h5", _guides("n1"))
    sheet = _sheet(tmp_path, ["g1", "g2"])
    samples = tmp_path / "samples.tsv"
    samples.write_text("ES_m6_1\n")
    main(
        [
            "--format",
            "tsv",
            "--h5-uri",
            URI,
            "--h5-samples-file",
            str(samples),
            str(sheet),
        ]
    )
    captured = capsys.readouterr()
    assert captured.out.startswith("gsig1:set:n=2:")
    assert f"note: checked 1 of 2 {H5_NAME} file(s); samples from {samples}" in (
        captured.err
    )
    assert "error:" not in captured.err


def test_cli_bad_samples_file_fails_before_s3(tmp_path: Path, fake_s3, capsys) -> None:
    _files, state = fake_s3
    sheet = _sheet(tmp_path, ["g1"])
    samples = tmp_path / "samples.tsv"
    samples.write_text("ES_m6_1\tES_m8_1\n")
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--format",
                "tsv",
                "--h5-uri",
                URI,
                "--h5-samples-file",
                str(samples),
                str(sheet),
            ]
        )
    assert exc.value.code == 1
    assert "line 1: holds a tab" in capsys.readouterr().err
    assert state["client_made"] is False
