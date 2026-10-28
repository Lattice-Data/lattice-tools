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
    find_multi_layout,
    find_per_sample_layout,
    guide_ids_from_h5,
    is_run_level_h5,
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

# Run-level layout: one outs/multi/count/ per library, under a run folder.
MULTI_NAME = "raw_feature_bc_matrix.h5"
LIB1 = "proj/order/LIB1/"
LIB2 = "proj/order/LIB2/"
RUN1 = f"{LIB1}processed/cellranger/Run_1/outs"
RUN2 = f"{LIB1}processed/cellranger/Run_2/outs"
LIB2_RUN1 = f"{LIB2}processed/cellranger/Run_1/outs"
MULTI_KEY = f"{RUN1}/multi/count/{MULTI_NAME}"
MULTI_KEY_RUN2 = f"{RUN2}/multi/count/{MULTI_NAME}"
MULTI_KEY_LIB2 = f"{LIB2_RUN1}/multi/count/{MULTI_NAME}"
# Cell Ranger 10.x writes the same basename flat at the top of outs/.
FLAT_KEY = f"{RUN1}/{MULTI_NAME}"
FLAT_KEY_LIB2 = f"{LIB2_RUN1}/{MULTI_NAME}"


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


def test_find_per_sample_layout_flags_sample_and_library_folders_without_h5() -> None:
    """LIB2 has processed output but no per_sample_outs/ at all: reported.

    Whether a library holds a run-level h5 is find_multi_layout's question, so
    this layout reports only the missing per_sample_outs/.
    """
    keys = [
        KEY_S1,
        f"{OUTS}/S1/metrics_summary.csv",
        f"{OUTS}/S2/metrics_summary.csv",  # S2 has no h5
        "proj/order/LIB1/raw/x_R1_001.fastq.gz",
        "proj/order/LIB2/processed/run1/outs/web_summary.html",  # LIB2 has no h5
        "proj/order/LIB2/raw/y_R1_001.fastq.gz",
    ]
    layout = find_per_sample_layout(keys, H5_NAME)
    assert layout.h5_keys == [KEY_S1]
    assert layout.sample_dirs_without_h5 == [f"{OUTS}/S2/"]
    assert layout.library_dirs_without_per_sample_outs == ["proj/order/LIB2/"]


def test_find_per_sample_layout_matches_basename_exactly() -> None:
    keys = [f"{OUTS}/S1/count/raw_{H5_NAME}", KEY_S2]
    layout = find_per_sample_layout(keys, H5_NAME)
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
        scope="per-sample",
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
        scope="per-sample",
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
        scope="per-sample",
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
        scope="per-sample",
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
        scope="per-sample",
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
        scope="per-sample",
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
    main(
        [
            "--format",
            "tsv",
            "--h5-uri",
            URI,
            "--h5-scope",
            "per-sample",
            str(sheet),
        ]
    )
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
        main(
            [
                "--format",
                "tsv",
                "--h5-uri",
                URI,
                "--h5-scope",
                "per-sample",
                str(sheet),
            ]
        )
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
        ("--h5-scope", "per-sample"),
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

OUTS_A = "proj/order/tb7k_q2/processed/run1/outs/per_sample_outs"
OUTS_B = "proj/order/Wtb7k_q2/processed/run1/outs/per_sample_outs"
KEY_A = f"{OUTS_A}/tb7k_q2/count/{H5_NAME}"
KEY_B = f"{OUTS_B}/Wtb7k_q2/count/{H5_NAME}"


def test_read_sample_names_skips_blanks_comments_and_duplicates(
    tmp_path: Path,
) -> None:
    path = tmp_path / "samples.tsv"
    path.write_text("# listed samples\ntb7k_q2\n\n  pm3x_v8  \ntb7k_q2\n")
    assert read_sample_names(path) == ["tb7k_q2", "pm3x_v8"]


@pytest.mark.parametrize(
    "content, message",
    [
        ("tb7k_q2\npm3x_v8\tlibrary\n", "line 2: holds a tab"),
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
    """tb7k_q2 is an exact name: it must not pull in Wtb7k_q2, which contains it."""
    files = {
        KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1")),
        KEY_B: _write_h5(tmp_path / "b.h5", _guides("n1")),
    }
    read: list[str] = []
    local = _local_reader(files)

    def reader(bucket: str, key: str) -> set[str]:
        read.append(key)
        return local(bucket, key)

    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2"],
        samples_label="samples.tsv",
        scope="per-sample",
        s3_client=MockS3Client(keys=list(files)),
        reader=reader,
    )
    assert report.ok
    assert read == [KEY_A]
    assert report.notes == [
        f"checked 1 of 2 {H5_NAME} file(s); samples from samples.tsv"
    ]


def test_listed_sample_without_h5_fails(tmp_path: Path) -> None:
    files = {KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2", "tb7k_q9"],
        samples_label="samples.tsv",
        scope="per-sample",
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"sample 'tb7k_q9' listed in samples.tsv has no {H5_NAME} under {URI}"
    ]
    assert report.h5_total == 1


def test_sample_list_matching_nothing_reports_each_name() -> None:
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2"],
        samples_label="samples.tsv",
        scope="per-sample",
        s3_client=MockS3Client(keys=[KEY_B]),
        reader=lambda bucket, key: set(),
    )
    assert report.failures == [
        f"sample 'tb7k_q2' listed in samples.tsv has no {H5_NAME} under {URI}"
    ]
    assert report.notes == [
        f"checked 0 of 1 {H5_NAME} file(s); samples from samples.tsv"
    ]


def test_missing_folder_checks_respect_sample_list(tmp_path: Path) -> None:
    """Only listed samples' folders are checked for a missing h5."""
    files = {KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1"))}
    keys = [
        *files,
        f"{OUTS_B}/Wtb7k_q2/metrics_summary.csv",  # not listed, no h5
        "proj/order/rz9d_h4/processed/run1/outs/web_summary.html",
    ]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2"],
        scope="per-sample",
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert report.ok


def test_cli_h5_samples_file(tmp_path: Path, fake_s3, capsys) -> None:
    files, _state = fake_s3
    files[KEY_A] = _write_h5(tmp_path / "a.h5", _guides("g1", "g2"))
    files[KEY_B] = _write_h5(tmp_path / "b.h5", _guides("n1"))
    sheet = _sheet(tmp_path, ["g1", "g2"])
    samples = tmp_path / "samples.tsv"
    samples.write_text("tb7k_q2\n")
    main(
        [
            "--format",
            "tsv",
            "--h5-uri",
            URI,
            "--h5-scope",
            "per-sample",
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
    samples.write_text("tb7k_q2\tpm3x_v8\n")
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--format",
                "tsv",
                "--h5-uri",
                URI,
                "--h5-scope",
                "per-sample",
                "--h5-samples-file",
                str(samples),
                str(sheet),
            ]
        )
    assert exc.value.code == 1
    assert "line 1: holds a tab" in capsys.readouterr().err
    assert state["client_made"] is False


def test_listed_sample_with_empty_folder_reported_once(tmp_path: Path) -> None:
    """A folder failure already names the sample; no second listed-sample line."""
    files = {KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1"))}
    keys = [*files, f"{OUTS_A}/tb7k_q5/metrics_summary.csv"]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2", "tb7k_q5"],
        samples_label="samples.tsv",
        scope="per-sample",
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"sample folder s3://{BUCKET}/{OUTS_A}/tb7k_q5/ has no {H5_NAME}"
    ]


def test_missing_h5_in_same_named_library_reported_once(tmp_path: Path) -> None:
    """Library and sample share a name; a missing h5 is one sample-folder line."""
    files = {KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1"))}
    keys = [
        *files,
        f"{OUTS_B}/Wtb7k_q2/metrics_summary.csv",
        "proj/order/Wtb7k_q2/raw/x_R1_001.fastq.gz",
    ]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        scope="per-sample",
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"sample folder s3://{BUCKET}/{OUTS_B}/Wtb7k_q2/ has no {H5_NAME}"
    ]


def test_library_with_sample_folders_not_reported() -> None:
    """LIB1 has sample folders (one with, one without an h5): only S2 is reported."""
    keys = [KEY_S1, f"{OUTS}/S2/metrics_summary.csv"]
    layout = find_per_sample_layout(keys, H5_NAME)
    assert layout.library_dirs_without_per_sample_outs == []
    assert layout.sample_dirs_without_h5 == [f"{OUTS}/S2/"]


def test_library_without_per_sample_outs_fails_the_run(tmp_path: Path) -> None:
    files = {KEY_S1: _write_h5(tmp_path / "s1.h5", _guides("g1"))}
    keys = [*files, "proj/order/LIB2/processed/run1/outs/web_summary.html"]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        scope="per-sample",
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    expected = f"library folder s3://{BUCKET}/proj/order/LIB2/ has no per_sample_outs/"
    assert report.failures == [expected]


def test_listed_sample_matching_library_failure_reported_once(tmp_path: Path) -> None:
    files = {KEY_A: _write_h5(tmp_path / "a.h5", _guides("g1"))}
    keys = [*files, "proj/order/pm3x_v8/processed/run1/outs/web_summary.html"]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        samples=["tb7k_q2", "pm3x_v8"],
        samples_label="samples.tsv",
        scope="per-sample",
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    expected = (
        f"library folder s3://{BUCKET}/proj/order/pm3x_v8/ has no per_sample_outs/"
    )
    assert report.failures == [expected]


# --- multi scope: one run-level h5 per library --------------------------------


def test_find_multi_layout_groups_one_h5_per_library() -> None:
    keys = [
        MULTI_KEY,
        MULTI_KEY_LIB2,
        f"{RUN1}/multi/count/feature_reference.csv",
        f"{RUN1}/per_sample_outs/S1/count/{H5_NAME}",
    ]
    layout = find_multi_layout(keys, MULTI_NAME)
    assert layout.h5_keys == sorted([MULTI_KEY, MULTI_KEY_LIB2])
    assert layout.libraries_without_h5 == []
    assert layout.ambiguous_libraries == []


def test_find_multi_layout_ignores_per_sample_h5() -> None:
    """A per-sample h5 is not under outs/multi/count/, so multi scope skips it."""
    layout = find_multi_layout([KEY_S1, KEY_S2], MULTI_NAME)
    assert layout.h5_keys == []
    assert layout.libraries_without_h5 == [LIB1]


def test_find_multi_layout_reports_library_without_h5() -> None:
    keys = [MULTI_KEY, f"{LIB2_RUN1}/web_summary.html"]
    layout = find_multi_layout(keys, MULTI_NAME)
    assert layout.h5_keys == [MULTI_KEY]
    assert layout.libraries_without_h5 == [LIB2]


def test_find_multi_layout_reports_two_run_folders_as_ambiguous() -> None:
    """Reprocessing left two run folders: which describes the delivery is undecidable."""
    layout = find_multi_layout([MULTI_KEY, MULTI_KEY_RUN2], MULTI_NAME)
    assert layout.h5_keys == []
    assert layout.libraries_without_h5 == []
    assert layout.ambiguous_libraries == [(LIB1, sorted([MULTI_KEY, MULTI_KEY_RUN2]))]


def test_multi_scope_is_the_default_and_labels_the_run(tmp_path: Path) -> None:
    files = {MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1", "g2", "g9"))}
    report = check_guide_ids_against_h5(
        {"g1", "g2"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.ok
    assert (report.h5_total, report.h5_failed) == (1, 0)
    assert report.notes == [
        f"library LIB1 (run Run_1) (s3://{BUCKET}/{MULTI_KEY}): 1 guide(s) in h5 "
        "not in metadata; examples: g9"
    ]


def test_multi_scope_reads_one_h5_per_library(tmp_path: Path) -> None:
    """Two libraries, one read each; the per-sample h5s alongside are not read."""
    files = {
        MULTI_KEY: _write_h5(tmp_path / "m1.h5", _guides("g1")),
        MULTI_KEY_LIB2: _write_h5(tmp_path / "m2.h5", _guides("g1")),
    }
    read: list[str] = []
    local = _local_reader(files)

    def reader(bucket: str, key: str) -> set[str]:
        read.append(key)
        return local(bucket, key)

    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=[*files, KEY_S1, KEY_S2]),
        reader=reader,
    )
    assert report.ok
    assert sorted(read) == sorted([MULTI_KEY, MULTI_KEY_LIB2])
    assert report.h5_total == 2


def test_multi_scope_missing_guides_fail(tmp_path: Path) -> None:
    files = {MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1", "g2"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"library LIB1 (run Run_1) (s3://{BUCKET}/{MULTI_KEY}): 1 of 2 guide_ids "
        "missing from h5; examples: g2"
    ]


def test_multi_scope_library_without_h5_fails(tmp_path: Path) -> None:
    files = {MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1"))}
    keys = [*files, f"{LIB2_RUN1}/web_summary.html"]
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=keys),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"library folder s3://{BUCKET}/{LIB2} has no run-level {MULTI_NAME} "
        "under outs/ and outs/multi/count/"
    ]
    assert report.h5_total == 1


def test_multi_scope_ambiguous_library_is_not_read(tmp_path: Path) -> None:
    files = {
        MULTI_KEY: _write_h5(tmp_path / "m1.h5", _guides("g1")),
        MULTI_KEY_RUN2: _write_h5(tmp_path / "m2.h5", _guides("g1")),
    }
    read: list[str] = []

    def reader(bucket: str, key: str) -> set[str]:
        read.append(key)
        return {"g1"}

    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=reader,
    )
    assert read == []
    assert report.h5_total == 0
    assert len(report.failures) == 1
    failure = report.failures[0]
    assert (
        f"library folder s3://{BUCKET}/{LIB1} has 2 run-level {MULTI_NAME} files"
        in failure
    )
    assert "expected exactly one, so none were read" in failure
    assert f"s3://{BUCKET}/{MULTI_KEY}" in failure
    assert f"s3://{BUCKET}/{MULTI_KEY_RUN2}" in failure


def test_multi_scope_with_no_multi_output_suggests_per_sample() -> None:
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=["proj/order/LIB1/raw/x_R1_001.fastq.gz"]),
        reader=lambda bucket, key: set(),
    )
    assert report.failures == [
        f"no run-level {MULTI_NAME} found under {URI} (looked under outs/ and "
        "outs/multi/count/); if this order has no cellranger multi output, use "
        "--h5-scope per-sample"
    ]


def test_scope_both_runs_each_scope(tmp_path: Path) -> None:
    files = {
        MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1")),
        f"{RUN1}/per_sample_outs/S1/count/{H5_NAME}": _write_h5(
            tmp_path / "s1.h5", _guides("g1")
        ),
    }
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        scope="both",
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.ok
    assert report.h5_total == 2


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"scope": "sample"}, "unknown scope 'sample'"),
        ({"scope": "both", "filename": "x.h5"}, "cannot apply to scope 'both'"),
        ({"samples": ["S1"]}, "does not apply to scope 'multi'"),
    ],
)
def test_scope_argument_errors(kwargs: dict, message: str) -> None:
    with pytest.raises(H5CheckError, match=message):
        check_guide_ids_against_h5(
            {"g1"}, URI, s3_client=MockS3Client(), reader=lambda b, k: set(), **kwargs
        )


def test_cli_multi_scope_is_the_default(tmp_path: Path, fake_s3, capsys) -> None:
    files, _state = fake_s3
    files[MULTI_KEY] = _write_h5(tmp_path / "m.h5", _guides("g1", "g2"))
    sheet = _sheet(tmp_path, ["g1", "g2"])
    main(["--format", "tsv", "--h5-uri", URI, str(sheet)])
    captured = capsys.readouterr()
    assert captured.out.startswith("gsig1:set:n=2:")
    assert "error:" not in captured.err


@pytest.mark.parametrize(
    "extra_args, message",
    [
        (
            ["--h5-samples-file", "samples.tsv"],
            "--h5-samples-file narrows the per-sample sweep",
        ),
        (
            ["--h5-scope", "both", "--h5-filename", "x.h5"],
            "--h5-filename is not allowed with --h5-scope both",
        ),
    ],
)
def test_cli_scope_usage_errors(capsys, extra_args: list[str], message: str) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--h5-uri", URI, *extra_args, "a.tsv"])
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


# --- both Cell Ranger layouts -------------------------------------------------


@pytest.mark.parametrize(
    "key, expected",
    [
        (MULTI_KEY, True),  # 9.0.1 nested
        (FLAT_KEY, True),  # 10.x flat
        (KEY_S1, False),  # per-sample, different basename
        (f"{RUN1}/per_sample_outs/{MULTI_NAME}", False),  # per_sample_outs is not outs
        (f"{RUN1}/multi/{MULTI_NAME}", False),  # multi/ without count/
        (f"{RUN1}/analysis/{MULTI_NAME}", False),  # some other outs/ subfolder
    ],
)
def test_is_run_level_h5_matches_segments_not_suffixes(
    key: str, expected: bool
) -> None:
    assert is_run_level_h5(key, MULTI_NAME) is expected


def test_find_multi_layout_accepts_the_flat_layout() -> None:
    """Cell Ranger 10.x writes no multi/count/ subtree at all."""
    keys = [
        FLAT_KEY,
        f"{RUN1}/filtered_feature_bc_matrix.h5",
        f"{RUN1}/feature_reference.csv",
        f"{RUN1}/per_sample_outs/S1/count/{H5_NAME}",
    ]
    layout = find_multi_layout(keys, MULTI_NAME)
    assert layout.h5_keys == [FLAT_KEY]
    assert layout.libraries_without_h5 == []
    assert layout.ambiguous_libraries == []


def test_find_multi_layout_mixes_layouts_across_libraries() -> None:
    """One order may hold a 9.0.1 library and a 10.x library side by side."""
    layout = find_multi_layout([MULTI_KEY, FLAT_KEY_LIB2], MULTI_NAME)
    assert layout.h5_keys == sorted([MULTI_KEY, FLAT_KEY_LIB2])
    assert layout.libraries_without_h5 == []


def test_find_multi_layout_both_layouts_in_one_library_is_ambiguous() -> None:
    """Reprocessing across versions leaves two run-level h5s under one library."""
    layout = find_multi_layout([MULTI_KEY, FLAT_KEY], MULTI_NAME)
    assert layout.h5_keys == []
    assert layout.ambiguous_libraries == [(LIB1, sorted([MULTI_KEY, FLAT_KEY]))]


def test_flat_layout_is_read_and_labelled(tmp_path: Path) -> None:
    files = {FLAT_KEY: _write_h5(tmp_path / "flat.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.ok
    assert report.h5_total == 1


def test_flat_layout_missing_guides_name_the_run(tmp_path: Path) -> None:
    files = {FLAT_KEY: _write_h5(tmp_path / "flat.h5", _guides("g1"))}
    report = check_guide_ids_against_h5(
        {"g1", "g2"},
        URI,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert report.failures == [
        f"library LIB1 (run Run_1) (s3://{BUCKET}/{FLAT_KEY}): 1 of 2 guide_ids "
        "missing from h5; examples: g2"
    ]


# --- per-scope worker defaults ------------------------------------------------


@pytest.mark.parametrize(
    "scope, expected",
    [
        ("multi", h5check.MULTI_SCOPE_WORKERS),
        ("per-sample", h5check.PER_SAMPLE_SCOPE_WORKERS),
        ("both", h5check.PER_SAMPLE_SCOPE_WORKERS),
    ],
)
def test_default_workers_for_scope(scope: str, expected: int) -> None:
    assert h5check.default_workers_for(scope) == expected


def test_multi_scope_opens_more_threads_than_per_sample() -> None:
    """The point of the split: multi opens one h5 per library, so it can widen."""
    assert h5check.MULTI_SCOPE_WORKERS > h5check.PER_SAMPLE_SCOPE_WORKERS


def test_cli_worker_help_defaults_match_h5check() -> None:
    """cli spells the defaults out for help text; keep the two in step."""
    from guidesig import cli

    assert cli.H5_MULTI_WORKERS_HELP_DEFAULT == h5check.MULTI_SCOPE_WORKERS
    assert cli.H5_PER_SAMPLE_WORKERS_HELP_DEFAULT == h5check.PER_SAMPLE_SCOPE_WORKERS


@pytest.mark.parametrize("scope", ["multi", "per-sample", "both"])
def test_explicit_workers_overrides_the_scope_default(
    tmp_path: Path, monkeypatch, scope: str
) -> None:
    seen: list[int] = []
    monkeypatch.setattr(
        h5check,
        "h5_worker_ceiling",
        lambda *, do_introspect, workers: seen.append(workers) or workers,
    )
    files = {MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1"))}
    check_guide_ids_against_h5(
        {"g1"},
        URI,
        scope=scope,
        workers=3,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert seen == [3]


@pytest.mark.parametrize(
    "scope, expected",
    [
        ("multi", 24),
        ("per-sample", 8),
        ("both", 8),
    ],
)
def test_unset_workers_uses_the_scope_default(
    tmp_path: Path, monkeypatch, scope: str, expected: int
) -> None:
    seen: list[int] = []
    monkeypatch.setattr(
        h5check,
        "h5_worker_ceiling",
        lambda *, do_introspect, workers: seen.append(workers) or workers,
    )
    files = {MULTI_KEY: _write_h5(tmp_path / "m.h5", _guides("g1"))}
    check_guide_ids_against_h5(
        {"g1"},
        URI,
        scope=scope,
        s3_client=MockS3Client(keys=list(files)),
        reader=_local_reader(files),
    )
    assert seen == [expected]


# --- S3 block size for the features read --------------------------------------


def test_features_block_size_is_well_below_the_s3fs_default() -> None:
    """The whole point of setting it: s3fs would otherwise fetch 50 MiB a seek."""
    assert h5check.FEATURES_BLOCK_SIZE == 1 * 2**20
    assert h5check.FEATURES_BLOCK_SIZE < 50 * 2**20


def test_read_s3_h5_guide_ids_passes_block_size_to_fs_open(
    tmp_path: Path, monkeypatch
) -> None:
    """Pin the block size, and the call shape that actually delivers it.

    fsspec.open forwards keywords as storage options to the filesystem
    constructor, where S3FileSystem takes default_block_size and not
    block_size, so that form cannot deliver one. A refactor back to it, or
    dropping the argument, would restore 50 MiB blocks and 25-second reads
    while leaving every other test green, so this asserts the kwarg reaches
    fs.open with the value we measured.
    """
    fsspec = pytest.importorskip("fsspec")
    local_h5 = _write_h5(tmp_path / "m.h5", _guides("g1", "g2"))
    captured: dict[str, object] = {}

    class FakeFS:
        def open(self, path, mode, *, block_size):
            captured.update(path=path, mode=mode, block_size=block_size)
            return open(local_h5, "rb")

    def fake_filesystem(protocol, **storage_options):
        captured.update(protocol=protocol, storage_options=storage_options)
        return FakeFS()

    monkeypatch.setattr(fsspec, "filesystem", fake_filesystem)

    guides = h5check.read_s3_h5_guide_ids(BUCKET, "lib/m.h5", pool_size=24)

    assert guides == {"g1", "g2"}
    assert captured["block_size"] == h5check.FEATURES_BLOCK_SIZE
    assert captured["protocol"] == "s3"
    assert captured["path"] == f"s3://{BUCKET}/lib/m.h5"
    assert captured["mode"] == "rb"
    # pool_size still reaches the client config, as file_extract passes it.
    assert captured["storage_options"] == {
        "config_kwargs": {"max_pool_connections": 24}
    }
