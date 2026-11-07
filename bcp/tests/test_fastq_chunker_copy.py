"""Copying files that need no split (`run --copy-singletons`), over file://.

No pigz or split needed: a plan with nothing to split only copies.
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest

from fastq_chunker import portal, s3io
from fastq_chunker import run as run_mod
from fastq_chunker.cli import main
from fastq_chunker.plan import (
    GB,
    FastqFile,
    load_plan,
    make_plan,
    plan_group,
    write_plan,
)
from fastq_chunker.run import (
    CHECKED_CRC,
    CHECKED_SIZE,
    RunError,
    RunOptions,
    check_copy,
    run_plan,
)


def make_world(tmp_path: Path, n_sets: int = 3):
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    groups = []
    for i in range(1, n_sets + 1):
        files = []
        for slot, mate in (("read1", "R1"), ("read2", "R2")):
            name = f"run{i}-LIB_{mate}_001.fastq.gz"
            data = f"{name}\n".encode() * (10 + i)
            (src / name).write_bytes(data)
            files.append(
                FastqFile(
                    group=f"G{i}",
                    slot=slot,
                    filename=name,
                    s3_uri=(src / name).as_uri(),
                    size_bytes=len(data),
                    reads=1000,
                )
            )
        groups.append(plan_group(files))
    plan = make_plan(groups, dst.as_uri() + "/", 80 * GB, 100 * GB, 1_000_000)
    return plan, src, dst


def opts(tmp_path, **kw) -> RunOptions:
    base = dict(
        copy_singletons=True,
        copy_workers=3,
        log_dir=tmp_path / "logs",
        manifest=tmp_path / "run_manifest.tsv",
    )
    base.update(kw)
    return RunOptions(**base)


def manifest(path: Path) -> dict[str, dict]:
    with open(path, newline="") as f:
        return {r["chunk"]: r for r in csv.DictReader(f, delimiter="\t")}


def names(plan):
    return sorted(fp.file.filename for g in plan.groups for fp in g.files)


def test_copies_every_unsplit_file_and_records_each(tmp_path):
    plan, src, dst = make_world(tmp_path)
    out = io.StringIO()
    assert run_plan(plan, opts(tmp_path), out=out) == 0
    assert sorted(p.name for p in dst.iterdir()) == names(plan)
    for name in names(plan):
        assert (dst / name).read_bytes() == (src / name).read_bytes()
    rows = manifest(tmp_path / "run_manifest.tsv")
    assert sorted(rows) == names(plan)
    assert {r["status"] for r in rows.values()} == {"ok"}
    assert {r["checked"] for r in rows.values()} == {CHECKED_SIZE}
    assert rows[names(plan)[0]]["bytes"] == str((src / names(plan)[0]).stat().st_size)


def test_rerun_skips_good_copies_and_redoes_a_wrong_one(tmp_path):
    plan, src, dst = make_world(tmp_path)
    assert run_plan(plan, opts(tmp_path), out=io.StringIO()) == 0
    bad = names(plan)[1]
    (dst / bad).write_bytes(b"truncated")
    out = io.StringIO()
    assert run_plan(plan, opts(tmp_path), out=out) == 0
    rows = manifest(tmp_path / "run_manifest.tsv")
    assert rows[bad]["status"] == "ok"
    assert (dst / bad).read_bytes() == (src / bad).read_bytes()
    assert sum(r["status"] == "skipped" for r in rows.values()) == len(rows) - 1


def test_force_copies_again(tmp_path):
    plan, _, _ = make_world(tmp_path, n_sets=1)
    run_plan(plan, opts(tmp_path), out=io.StringIO())
    run_plan(plan, opts(tmp_path, force=True), out=io.StringIO())
    rows = manifest(tmp_path / "run_manifest.tsv")
    assert {r["status"] for r in rows.values()} == {"ok"}


def test_missing_source_stops_before_any_copy(tmp_path):
    plan, src, dst = make_world(tmp_path)
    (src / names(plan)[0]).unlink()
    with pytest.raises(RunError, match="missing:"):
        run_plan(plan, opts(tmp_path), out=io.StringIO())
    assert list(dst.iterdir()) == []


def test_source_size_drift_stops_before_any_copy(tmp_path):
    plan, src, dst = make_world(tmp_path)
    (src / names(plan)[0]).write_bytes(b"changed")
    with pytest.raises(RunError, match="size mismatch"):
        run_plan(plan, opts(tmp_path), out=io.StringIO())
    assert list(dst.iterdir()) == []


class FlakyFS:
    """A filesystem whose copy fails for one file name, or writes it short."""

    def __init__(self, fs, name: str, mode: str):
        self._fs, self._name, self._mode = fs, name, mode

    def __getattr__(self, attr):
        return getattr(self._fs, attr)

    def copy(self, src, dst):
        if dst.endswith("/" + self._name):
            if self._mode == "raise":
                raise OSError("simulated copy failure")
            self._fs.pipe(dst, b"short")
            return
        self._fs.copy(src, dst)


@pytest.mark.parametrize(
    "mode, detail", [("raise", "simulated copy failure"), ("short", "size 5, plan")]
)
def test_one_failed_copy_does_not_stop_the_others(tmp_path, monkeypatch, mode, detail):
    plan, src, dst = make_world(tmp_path)
    bad = names(plan)[2]
    real = s3io.make_fs
    monkeypatch.setattr(s3io, "make_fs", lambda url: FlakyFS(real(url), bad, mode))
    out = io.StringIO()
    assert run_plan(plan, opts(tmp_path), out=out) == 1
    text = out.getvalue()
    assert f"FAILED copy {bad}: " in text and detail in text
    assert f"rerun with: --only {bad}" in text
    assert not (dst / bad).exists()  # a failed copy leaves nothing behind
    rows = manifest(tmp_path / "run_manifest.tsv")
    assert rows[bad]["status"] == "failed" and rows[bad]["bytes"] == ""
    assert sum(r["status"] == "ok" for r in rows.values()) == len(rows) - 1


def test_only_limits_the_copies(tmp_path):
    plan, _, dst = make_world(tmp_path)
    one = names(plan)[0]
    assert run_plan(plan, opts(tmp_path, only=[one]), out=io.StringIO()) == 0
    assert [p.name for p in dst.iterdir()] == [one]


def test_dry_run_lists_copies_and_needs_no_tools(tmp_path, monkeypatch):
    plan, src, dst = make_world(tmp_path, n_sets=1)

    def no_tools(*a, **k):
        raise AssertionError("tools looked up for a copy-only plan")

    monkeypatch.setattr(run_mod, "find_tools", no_tools)
    out = io.StringIO()
    assert run_plan(plan, opts(tmp_path, dry_run=True), out=out) == 0
    lines = out.getvalue().splitlines()
    assert lines == [f"copy {(src / n).as_uri()} {plan.dst}{n}" for n in names(plan)]
    assert list(dst.iterdir()) == []


def test_run_command_on_a_copy_only_plan_needs_no_tools(tmp_path, monkeypatch):
    plan, _, dst = make_world(tmp_path, n_sets=1)
    write_plan(plan, tmp_path / "plan.json")
    monkeypatch.setattr("fastq_chunker.cli.find_tools", lambda *a, **k: 1 / 0)
    rc = main(
        [
            "run",
            "--plan",
            str(tmp_path / "plan.json"),
            "--copy-singletons",
            "--log-dir",
            str(tmp_path / "logs"),
            "--manifest",
            str(tmp_path / "m.tsv"),
        ]
    )
    assert rc == 0
    assert sorted(p.name for p in dst.iterdir()) == names(plan)


# --- the CRC64NVME check -----------------------------------------------------


class FakeS3:
    def __init__(self, crc, kind):
        self.crc, self.kind, self.calls = crc, kind, 0

    def head_object(self, Bucket, Key, ChecksumMode):
        self.calls += 1
        assert ChecksumMode == "ENABLED"
        out = {}
        if self.crc:
            out["ChecksumCRC64NVME"] = self.crc
        if self.kind:
            out["ChecksumType"] = self.kind
        return out


def copied_file(tmp_path, crc="AAAAAAAAAAA="):
    plan, src, dst = make_world(tmp_path, n_sets=1)
    fp = plan.groups[0].files[0]
    fp.file.__dict__["crc64nvme"] = crc  # frozen dataclass; test setup only
    (dst / fp.file.filename).write_bytes((src / fp.file.filename).read_bytes())
    target = f"s3://dst/{fp.file.filename}"
    fs = s3io.make_fs(plan.dst)
    local_target = plan.dst + fp.file.filename
    return fs, local_target, target, fp


@pytest.mark.parametrize(
    "crc, kind, problem, checked",
    [
        ("AAAAAAAAAAA=", "FULL_OBJECT", None, CHECKED_CRC),
        (
            "BBBBBBBBBBB=",
            "FULL_OBJECT",
            "CRC64NVME BBBBBBBBBBB=, DB2 AAAAAAAAAAA=",
            CHECKED_CRC,
        ),
        ("AAAAAAAAAAA=", "COMPOSITE", None, CHECKED_SIZE),
        (None, None, None, CHECKED_SIZE),
    ],
)
def test_check_copy_compares_full_object_crc_only(
    tmp_path, monkeypatch, crc, kind, problem, checked
):
    fs, local_target, _, fp = copied_file(tmp_path)
    client = FakeS3(crc, kind)
    # head_crc64nvme parses an s3 URL; the size check runs on the local copy
    monkeypatch.setattr(
        s3io, "split_s3_url", lambda url: ("dst", url.rsplit("/", 1)[-1])
    )
    assert check_copy(fs, local_target, fp, client) == (problem, checked)
    assert client.calls == 1


def test_check_copy_skips_head_without_a_db2_crc(tmp_path):
    fs, local_target, _, fp = copied_file(tmp_path, crc=None)
    client = FakeS3("AAAAAAAAAAA=", "FULL_OBJECT")
    assert check_copy(fs, local_target, fp, client) == (None, CHECKED_SIZE)
    assert client.calls == 0


def test_check_copy_size_failure_comes_first(tmp_path):
    fs, local_target, _, fp = copied_file(tmp_path)
    Path(local_target.removeprefix("file://")).write_bytes(b"x")
    client = FakeS3("AAAAAAAAAAA=", "FULL_OBJECT")
    problem, _ = check_copy(fs, local_target, fp, client)
    assert problem.startswith("size 1, plan")
    assert client.calls == 0


# --- the crc travels from DB2 through the plan -------------------------------


def test_portal_carries_crc64nvme_into_the_plan(tmp_path):
    sets = [
        {
            "@id": "/sequence_file_sets/s/",
            "read1": "/sequence_files/f1/",
        }
    ]
    files = [
        {
            "@id": "/sequence_files/f1/",
            "file_format": "fastq",
            "file_size": 10,
            "read_count": 1,
            "s3_uri": "s3://b/a_R1.fastq.gz",
            "crc64nvme_base64": "Q1JDPQ==",
        }
    ]
    [members] = portal.build_groups(sets, files)
    assert members[0].crc64nvme == "Q1JDPQ=="
    plan = make_plan([plan_group(members)], "s3://d/", 80 * GB, 100 * GB, 1_000_000)
    write_plan(plan, tmp_path / "p.json")
    assert (
        load_plan(tmp_path / "p.json").groups[0].files[0].file.crc64nvme == "Q1JDPQ=="
    )


def test_a_plan_written_before_crc64nvme_still_loads(tmp_path):
    plan, _, _ = make_world(tmp_path, n_sets=1)
    write_plan(plan, tmp_path / "p.json")
    d = json.loads((tmp_path / "p.json").read_text())
    for g in d["groups"]:
        for fp in g["files"]:
            del fp["file"]["crc64nvme"]
    (tmp_path / "old.json").write_text(json.dumps(d))
    assert load_plan(tmp_path / "old.json").groups[0].files[0].file.crc64nvme is None
