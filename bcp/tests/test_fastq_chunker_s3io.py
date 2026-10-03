"""Round trips for the streaming endpoints over file:// and a local moto server."""

from __future__ import annotations

import gzip
import hashlib
import io
import random
import threading
import time
import os
import socket
import subprocess
import sys
from pathlib import Path

import pytest

from fastq_chunker import s3io

BCP_DIR = str(Path(__file__).resolve().parents[1])
MIB = 2**20


def run_endpoint(args, stdin=None, env=None):
    full_env = {**os.environ, **(env or {}), "PYTHONPATH": BCP_DIR}
    return subprocess.run(
        [sys.executable, "-m", "fastq_chunker.s3io", *args],
        input=stdin,
        capture_output=True,
        env=full_env,
        check=True,
    )


def test_put_then_get_round_trip_local(tmp_path):
    data = os.urandom(3 * MIB + 17)
    url = (tmp_path / "out" / "x.bin").as_uri()
    put = run_endpoint(["put", url], stdin=data, env={s3io.BLOCK_ENV: str(MIB)})
    name, n, md5hex = put.stdout.decode().strip().split("\t")
    assert (name, int(n), md5hex) == ("x.bin", len(data), hashlib.md5(data).hexdigest())
    assert (tmp_path / "out" / "x.bin").read_bytes() == data
    s3io.write_sidecar(s3io.make_fs(url), url, md5hex)
    assert (tmp_path / "out" / "x.bin.md5").read_text() == f"{md5hex}  x.bin\n"

    got = run_endpoint(["get", url, "--md5"], env={s3io.BLOCK_ENV: str(MIB)})
    assert got.stdout == data
    assert got.stderr.decode().strip().splitlines()[-1] == f"md5\t{md5hex}"


def test_put_never_writes_a_sidecar(tmp_path):
    url = (tmp_path / "y.bin").as_uri()
    run_endpoint(["put", url], stdin=b"abc")
    assert (tmp_path / "y.bin").read_bytes() == b"abc"
    assert not (tmp_path / "y.bin.md5").exists()
    with pytest.raises(subprocess.CalledProcessError):
        run_endpoint(["put", url, "--sidecar"], stdin=b"abc")


def test_put_with_compressor(tmp_path):
    url = (tmp_path / "c.gz").as_uri()
    out = run_endpoint(
        ["put", url, "--compress", "gzip", "-c", "-1"], stdin=b"hello " * 1000
    )
    name, n, md5hex = out.stdout.decode().strip().split("\t")
    data = (tmp_path / "c.gz").read_bytes()
    assert gzip.decompress(data) == b"hello " * 1000
    assert (int(n), md5hex) == (len(data), hashlib.md5(data).hexdigest())


def test_failing_compressor_leaves_no_object(tmp_path):
    url = (tmp_path / "bad.gz").as_uri()
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "fastq_chunker.s3io",
            "put",
            url,
            "--compress",
            "sh",
            "-c",
            "head -c 100; exit 3",
        ],
        input=b"x" * 1000,
        capture_output=True,
        env={**os.environ, "PYTHONPATH": BCP_DIR},
    )
    assert proc.returncode != 0
    assert "exited 3" in proc.stderr.decode()
    assert not (tmp_path / "bad.gz").exists()


def test_empty_input_makes_empty_object(tmp_path):
    url = (tmp_path / "empty.bin").as_uri()
    out = run_endpoint(["put", url], stdin=b"").stdout.decode()
    assert out.split("\t")[1] == "0"
    assert (tmp_path / "empty.bin").stat().st_size == 0


def test_in_process_helpers(tmp_path):
    url = (tmp_path / "z.bin").as_uri()
    n, md5hex = s3io.put(url, src=io.BytesIO(b"hello"))
    assert (n, md5hex) == (5, hashlib.md5(b"hello").hexdigest())
    fs = s3io.make_fs(url)
    s3io.write_sidecar(fs, url, md5hex)
    assert s3io.read_sidecar(fs, url) == (md5hex, "z.bin")
    sink = io.BytesIO()
    assert s3io.get(url, out=sink) == md5hex
    assert sink.getvalue() == b"hello"


# ---- S3 through a local moto server -------------------------------------------


@pytest.fixture(scope="module")
def moto_endpoint():
    pytest.importorskip("moto")
    from moto.server import ThreadedMotoServer

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    server = ThreadedMotoServer(ip_address="127.0.0.1", port=port)
    server.start()
    env = {
        s3io.ENDPOINT_ENV: f"http://127.0.0.1:{port}",
        "AWS_ACCESS_KEY_ID": "testing",
        "AWS_SECRET_ACCESS_KEY": "testing",
        "AWS_DEFAULT_REGION": "us-east-1",
        "AWS_EC2_METADATA_DISABLED": "true",
    }
    old = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    import boto3

    boto3.client("s3", endpoint_url=env[s3io.ENDPOINT_ENV]).create_bucket(Bucket="dst")
    try:
        yield env
    finally:
        server.stop()
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_multipart_round_trip_via_s3(moto_endpoint):
    # 12 MiB with a 5 MiB block: three parts, so multipart is exercised
    block = 5 * MIB
    data = os.urandom(12 * MIB + 5)
    env = {**moto_endpoint, s3io.BLOCK_ENV: str(block)}
    url = "s3://dst/_test/big.bin"
    put = run_endpoint(["put", url], stdin=data, env=env)
    name, n, md5hex = put.stdout.decode().strip().split("\t")
    assert (name, int(n), md5hex) == (
        "big.bin",
        len(data),
        hashlib.md5(data).hexdigest(),
    )

    fs = s3io.make_fs(url)
    fs.invalidate_cache()
    assert fs.size(url) == len(data)
    s3io.write_sidecar(fs, url, md5hex)
    assert s3io.read_sidecar(fs, url) == (md5hex, "big.bin")
    got = run_endpoint(["get", url, "--md5"], env=env)
    assert got.stdout == data
    assert got.stderr.decode().strip().splitlines()[-1] == f"md5\t{md5hex}"


def test_failing_compressor_aborts_multipart_via_s3(moto_endpoint):
    env = {**moto_endpoint, s3io.BLOCK_ENV: str(5 * MIB)}
    url = "s3://dst/_test/aborted.bin"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "fastq_chunker.s3io",
            "put",
            url,
            "--compress",
            "sh",
            "-c",
            "head -c 6000000; exit 3",
        ],
        input=os.urandom(7 * MIB),
        capture_output=True,
        env={**os.environ, **env, "PYTHONPATH": BCP_DIR},
    )
    assert proc.returncode != 0
    fs = s3io.make_fs(url)
    fs.invalidate_cache()
    assert not fs.exists(url)
    import boto3

    client = boto3.client("s3", endpoint_url=env[s3io.ENDPOINT_ENV])
    pending = client.list_multipart_uploads(Bucket="dst").get("Uploads", [])
    assert pending == []


class SlowRangeFS:
    """A fake fs whose cat_file answers out of order, to prove get reassembles in order."""

    protocol = ("s3", "s3a")

    def __init__(self, data: bytes):
        self.data = data
        self.calls = 0
        self.in_flight = 0
        self.peak = 0
        self.completed = 0
        self.lock = threading.Lock()

    def size(self, url):
        return len(self.data)

    def cat_file(self, url, start=None, end=None):
        with self.lock:
            self.calls += 1
            self.in_flight += 1
            self.peak = max(self.peak, self.in_flight)
        time.sleep(random.Random(start).random() / 50)
        with self.lock:
            self.in_flight -= 1
            self.completed += 1
        return self.data[start:end]


def test_ranged_blocks_reassemble_in_order_and_bound_concurrency():
    data = os.urandom(1000 * 7 + 13)
    fs = SlowRangeFS(data)
    out = b"".join(s3io.ranged_blocks(fs, "s3://x/y", block=1000, concurrency=4))
    assert out == data
    assert fs.calls == 8
    assert 1 < fs.peak <= 4


def test_ranged_blocks_empty_object():
    fs = SlowRangeFS(b"")
    assert list(s3io.ranged_blocks(fs, "s3://x/y", block=10, concurrency=3)) == []


def test_get_via_s3_uses_ranges(moto_endpoint):
    block = 1 * MIB
    data = os.urandom(11 * MIB + 3)
    env = {**moto_endpoint, s3io.BLOCK_ENV: str(5 * MIB)}
    url = "s3://dst/_test/ranged.bin"
    run_endpoint(["put", url], stdin=data, env=env)
    for concurrency in ("1", "4"):
        got = run_endpoint(
            ["get", url, "--md5", "--concurrency", concurrency],
            env={**moto_endpoint, s3io.BLOCK_ENV: str(block)},
        )
        assert got.stdout == data
        assert (
            got.stderr.decode().strip().splitlines()[-1]
            == f"md5\t{hashlib.md5(data).hexdigest()}"
        )


def test_get_takes_the_ranged_path_for_s3(monkeypatch):
    data = os.urandom(2500)
    fs = SlowRangeFS(data)  # has cat_file but no open(): a sequential read would fail
    monkeypatch.setattr(s3io, "make_fs", lambda url: fs)
    monkeypatch.setenv(s3io.BLOCK_ENV, "1000")
    sink = io.BytesIO()
    assert (
        s3io.get("s3://x/y", out=sink, concurrency=3) == hashlib.md5(data).hexdigest()
    )
    assert sink.getvalue() == data
    assert fs.calls == 3


def test_ranged_blocks_do_not_run_ahead_of_a_slow_consumer():
    # memory is bounded only if fetching stalls once `concurrency` blocks are
    # complete but not yet consumed
    fs = SlowRangeFS(os.urandom(1000 * 12))
    concurrency = 3
    for i, _ in enumerate(
        s3io.ranged_blocks(fs, "s3://x/y", block=1000, concurrency=concurrency)
    ):
        time.sleep(0.05)
        assert fs.completed - (i + 1) <= concurrency, (
            f"{fs.completed} done after {i + 1} consumed"
        )


class FakeS3Client:
    """Records multipart calls; upload_part sleeps a random time so parts
    complete out of order and in-flight counts are observable."""

    def __init__(self, fail_part: int | None = None):
        self.calls: list[str] = []
        self.parts: list[tuple[int, bytes]] = []
        self.completed_with = None
        self.aborted = False
        self.in_flight = 0
        self.peak = 0
        self.lock = threading.Lock()
        self.fail_part = fail_part

    def create_multipart_upload(self, Bucket, Key):
        self.calls.append("create")
        return {"UploadId": "u1"}

    def upload_part(self, Bucket, Key, UploadId, PartNumber, Body):
        assert UploadId == "u1"
        with self.lock:
            self.in_flight += 1
            self.peak = max(self.peak, self.in_flight)
        time.sleep(random.Random(PartNumber).random() / 50)
        with self.lock:
            self.in_flight -= 1
            self.parts.append((PartNumber, bytes(Body)))
        if PartNumber == self.fail_part:
            raise RuntimeError("part failed")
        return {"ETag": f"etag-{PartNumber}"}

    def complete_multipart_upload(self, Bucket, Key, UploadId, MultipartUpload):
        self.calls.append("complete")
        self.completed_with = MultipartUpload["Parts"]

    def abort_multipart_upload(self, Bucket, Key, UploadId):
        self.calls.append("abort")
        self.aborted = True

    def put_object(self, Bucket, Key, Body):
        self.calls.append(f"put_object:{len(Body)}")


def test_multipart_writer_orders_parts_and_bounds_in_flight():
    client = FakeS3Client()
    w = s3io.S3MultipartWriter("s3://b/k", concurrency=3, client=client)
    blocks = [os.urandom(10) for _ in range(9)]
    for b in blocks:
        w.write(b)
        # memory bound: never more than `concurrency` parts held for upload
        assert len(w.pending) < 3
    w.commit()
    assert client.calls[0] == "create" and client.calls[-1] == "complete"
    assert [p["PartNumber"] for p in client.completed_with] == list(range(1, 10))
    assert [p["ETag"] for p in client.completed_with] == [
        f"etag-{i}" for i in range(1, 10)
    ]
    assert sorted(client.parts) == [(i + 1, b) for i, b in enumerate(blocks)]
    assert 1 < client.peak <= 3
    assert not client.aborted


def test_multipart_writer_empty_object_uses_put_object():
    client = FakeS3Client()
    w = s3io.S3MultipartWriter("s3://b/k", concurrency=2, client=client)
    w.commit()
    assert client.calls == ["put_object:0"]


def test_multipart_writer_discard_aborts():
    client = FakeS3Client()
    w = s3io.S3MultipartWriter("s3://b/k", concurrency=2, client=client)
    w.write(b"x")
    w.discard()
    assert client.aborted and "complete" not in client.calls


def test_multipart_writer_failed_part_surfaces_and_put_aborts(monkeypatch):
    client = FakeS3Client(fail_part=2)
    real = s3io.S3MultipartWriter
    monkeypatch.setattr(s3io, "S3MultipartWriter", lambda url, c: real(url, c, client))

    # bypass make_fs for a fake s3: only is_s3 and invalidate_cache are touched
    class FakeFS:
        protocol = ("s3", "s3a")

        def invalidate_cache(self):
            pass

    monkeypatch.setattr(s3io, "make_fs", lambda url: FakeFS())
    monkeypatch.setenv(s3io.BLOCK_ENV, "4")
    with pytest.raises(RuntimeError, match="part failed"):
        s3io.put("s3://b/k", src=io.BytesIO(b"0123456789abcdef"), upload_concurrency=2)
    assert client.aborted and "complete" not in client.calls


def test_split_s3_url():
    assert s3io.split_s3_url("s3://bucket/a/b/c.gz") == ("bucket", "a/b/c.gz")
    assert s3io.split_s3_url("s3://bucket/k") == ("bucket", "k")


@pytest.mark.parametrize("fail_part", [4, "complete"])
def test_failure_during_commit_still_aborts(monkeypatch, fail_part):
    # with 4 parts at concurrency 2, parts 3 and 4 are collected inside commit()
    client = FakeS3Client(fail_part=4 if fail_part == 4 else None)
    if fail_part == "complete":

        def boom(**kwargs):
            raise RuntimeError("complete failed")

        client.complete_multipart_upload = boom
    real = s3io.S3MultipartWriter
    monkeypatch.setattr(s3io, "S3MultipartWriter", lambda url, c: real(url, c, client))

    class FakeFS:
        protocol = ("s3", "s3a")

        def invalidate_cache(self):
            pass

    monkeypatch.setattr(s3io, "make_fs", lambda url: FakeFS())
    monkeypatch.setenv(s3io.BLOCK_ENV, "4")
    with pytest.raises(RuntimeError):
        s3io.put("s3://b/k", src=io.BytesIO(b"0123456789abcdef"), upload_concurrency=2)
    assert client.aborted
    assert client.calls.count("abort") == 1
