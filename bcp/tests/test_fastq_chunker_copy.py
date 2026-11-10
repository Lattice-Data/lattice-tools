"""Copying files that need no split (`run --copy-singletons`), over file://.

No pigz or split needed: a plan with nothing to split only copies.
"""

from __future__ import annotations

import csv
import io
import json
import os
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


# --- S3 copies go through boto3's managed copy, parts in parallel ------------

MIB = 2**20


def test_copy_part_size_grows_to_stay_under_the_part_cap(monkeypatch):
    monkeypatch.delenv(s3io.COPY_PART_ENV, raising=False)
    assert s3io.copy_part_size(1) == s3io.COPY_PART
    assert s3io.copy_part_size(s3io.COPY_PART * s3io.MAX_PARTS) == s3io.COPY_PART
    huge = s3io.COPY_PART * s3io.MAX_PARTS + 1
    part = s3io.copy_part_size(huge)
    assert part > s3io.COPY_PART and -(-huge // part) <= s3io.MAX_PARTS
    monkeypatch.setenv(s3io.COPY_PART_ENV, str(5 * MIB))
    assert s3io.copy_part_size(1) == 5 * MIB


class RecordingClient:
    def __init__(self):
        self.calls = []

    def copy(self, source, bucket, key, ExtraArgs, Config):
        self.calls.append((source, bucket, key, ExtraArgs, Config))


def test_managed_copy_asks_for_crc64nvme_and_parallel_parts(monkeypatch):
    monkeypatch.delenv(s3io.COPY_PART_ENV, raising=False)
    client = RecordingClient()
    s3io.managed_copy(client, "s3://src/a/b.fastq.gz", "s3://dst/p/b.fastq.gz", 10, 7)
    [(source, bucket, key, extra, config)] = client.calls
    assert source == {"Bucket": "src", "Key": "a/b.fastq.gz"}
    assert (bucket, key) == ("dst", "p/b.fastq.gz")
    assert extra == {"ChecksumAlgorithm": "CRC64NVME"}
    assert config.max_concurrency == 7
    assert config.multipart_chunksize == s3io.COPY_PART
    assert config.multipart_threshold == s3io.COPY_THRESHOLD


def test_copy_options_must_be_positive(tmp_path):
    plan, _, _ = make_world(tmp_path, n_sets=1)
    for kw in ({"copy_parts": 0}, {"copy_workers": 0}):
        with pytest.raises(RunError, match="at least 1"):
            run_plan(plan, opts(tmp_path, **kw), out=io.StringIO())


def test_cli_passes_copy_parts_through(tmp_path, monkeypatch):
    plan, _, _ = make_world(tmp_path, n_sets=1)
    write_plan(plan, tmp_path / "plan.json")
    seen = {}

    def fake_run_plan(plan, opts, **kw):
        seen["opts"] = opts
        return 0

    monkeypatch.setattr("fastq_chunker.cli.run_plan", fake_run_plan)
    args = ["run", "--plan", str(tmp_path / "plan.json"), "--copy-singletons"]
    assert main(args + ["--copy-parts", "5", "--copy-workers", "2"]) == 0
    assert (seen["opts"].copy_parts, seen["opts"].copy_workers) == (5, 2)
    assert main(args) == 0
    assert seen["opts"].copy_parts == s3io.DEFAULT_COPY_PARTS
    assert seen["opts"].copy_workers == 4


@pytest.fixture(scope="module")
def moto_s3():
    pytest.importorskip("moto")
    import socket

    import boto3
    from moto.server import ThreadedMotoServer

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    server = ThreadedMotoServer(ip_address="127.0.0.1", port=port)
    server.start()
    env = {
        s3io.ENDPOINT_ENV: f"http://127.0.0.1:{port}",
        s3io.COPY_PART_ENV: str(5 * MIB),
        "AWS_ACCESS_KEY_ID": "testing",
        "AWS_SECRET_ACCESS_KEY": "testing",
        "AWS_DEFAULT_REGION": "us-east-1",
        "AWS_EC2_METADATA_DISABLED": "true",
    }
    old = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    client = boto3.client("s3", endpoint_url=env[s3io.ENDPOINT_ENV])
    for b in ("src", "dst"):
        client.create_bucket(Bucket=b)
    try:
        yield client
    finally:
        server.stop()
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def s3_world(client, prefix: str):
    """Three files in s3://src: one below the part threshold, two in parts."""
    files = []
    for i, size in enumerate((3 * MIB, 12 * MIB, 11 * MIB)):
        name = f"{prefix}{i}_R1_001.fastq.gz"
        data = bytes([i]) * size
        client.put_object(Bucket="src", Key=f"raw/{name}", Body=data)
        files.append(
            FastqFile(f"G{i}", "read1", name, f"s3://src/raw/{name}", size, 1000)
        )
    plan = make_plan(
        [plan_group([f]) for f in files], f"s3://dst/{prefix}/", 80 * GB, 100 * GB, 1
    )
    return plan, files


def test_s3_copies_run_in_parts_and_match_the_source(moto_s3, tmp_path, monkeypatch):
    plan, files = s3_world(moto_s3, "ok")
    pools, calls = [], []
    real_client, real_copy = s3io.make_s3_client, s3io.managed_copy
    monkeypatch.setattr(
        s3io, "make_s3_client", lambda n: pools.append(n) or real_client(n)
    )
    monkeypatch.setattr(
        s3io,
        "managed_copy",
        lambda *a: calls.append(a[1:]) or real_copy(*a),
    )
    out = io.StringIO()
    o = opts(tmp_path, copy_workers=2, copy_parts=3)
    assert run_plan(plan, o, out=out) == 0, out.getvalue()
    assert pools == [6]  # a connection for every part of every file in flight
    assert sorted(c[0] for c in calls) == sorted(f.s3_uri for f in files)
    assert all(c[3] == 3 for c in calls)
    for f in files:
        got = moto_s3.get_object(Bucket="dst", Key=f"ok/{f.filename}")
        etag = got["ETag"].strip('"')
        assert got["Body"].read() == bytes([int(f.group[1:])]) * f.size_bytes
        # below 5 MiB: one request; above: a multipart copy of 5 MiB parts
        assert ("-" in etag) == (f.size_bytes > 5 * MIB), (f.filename, etag)
    rows = manifest(tmp_path / "run_manifest.tsv")
    assert {r["status"] for r in rows.values()} == {"ok"}
    # a rerun finds every copy present and correct
    out = io.StringIO()
    assert run_plan(plan, o, out=out) == 0
    assert out.getvalue().count("present, skipping") == 3


def test_a_failed_s3_copy_is_aborted_and_reported(moto_s3, tmp_path, monkeypatch):
    plan, files = s3_world(moto_s3, "bad")
    victim = files[1]
    real_client = s3io.make_s3_client

    def flaky_client(n):
        client = real_client(n)
        real_part = client.upload_part_copy

        def part(**kw):
            if kw["Key"].endswith(victim.filename) and kw["PartNumber"] == 2:
                raise OSError("simulated part failure")
            return real_part(**kw)

        client.upload_part_copy = part
        return client

    monkeypatch.setattr(s3io, "make_s3_client", flaky_client)
    out = io.StringIO()
    assert run_plan(plan, opts(tmp_path, copy_parts=2), out=out) == 1
    text = out.getvalue()
    assert f"FAILED copy {victim.filename}: " in text
    assert "simulated part failure" in text
    listed = moto_s3.list_objects_v2(Bucket="dst", Prefix="bad/").get("Contents", [])
    assert sorted(o["Key"] for o in listed) == sorted(
        f"bad/{f.filename}" for f in files if f is not victim
    )
    # the multipart upload was aborted, not left behind to be billed
    uploads = moto_s3.list_multipart_uploads(Bucket="dst", Prefix="bad/")
    assert not uploads.get("Uploads")
