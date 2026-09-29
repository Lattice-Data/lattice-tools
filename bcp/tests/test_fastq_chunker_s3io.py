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
    put = run_endpoint(
        ["put", url, "--sidecar"], stdin=data, env={s3io.BLOCK_ENV: str(MIB)}
    )
    name, n, md5hex = put.stdout.decode().strip().split("\t")
    assert (name, int(n), md5hex) == ("x.bin", len(data), hashlib.md5(data).hexdigest())
    assert (tmp_path / "out" / "x.bin").read_bytes() == data
    assert (tmp_path / "out" / "x.bin.md5").read_text() == f"{md5hex}  x.bin\n"

    got = run_endpoint(["get", url, "--md5"], env={s3io.BLOCK_ENV: str(MIB)})
    assert got.stdout == data
    assert got.stderr.decode().strip().splitlines()[-1] == f"md5\t{md5hex}"


def test_put_writes_no_sidecar_by_default(tmp_path):
    url = (tmp_path / "y.bin").as_uri()
    run_endpoint(["put", url], stdin=b"abc")
    assert (tmp_path / "y.bin").read_bytes() == b"abc"
    assert not (tmp_path / "y.bin.md5").exists()


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
            "--sidecar",
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
    assert not (tmp_path / "bad.gz.md5").exists()


def test_empty_input_makes_empty_object(tmp_path):
    url = (tmp_path / "empty.bin").as_uri()
    out = run_endpoint(["put", url], stdin=b"").stdout.decode()
    assert out.split("\t")[1] == "0"
    assert (tmp_path / "empty.bin").stat().st_size == 0


def test_in_process_helpers(tmp_path):
    url = (tmp_path / "z.bin").as_uri()
    n, md5hex = s3io.put(url, src=io.BytesIO(b"hello"), sidecar=True)
    assert (n, md5hex) == (5, hashlib.md5(b"hello").hexdigest())
    fs = s3io.make_fs(url)
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
    put = run_endpoint(["put", url, "--sidecar"], stdin=data, env=env)
    name, n, md5hex = put.stdout.decode().strip().split("\t")
    assert (name, int(n), md5hex) == (
        "big.bin",
        len(data),
        hashlib.md5(data).hexdigest(),
    )

    fs = s3io.make_fs(url)
    fs.invalidate_cache()
    assert fs.size(url) == len(data)
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
            "--sidecar",
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
    assert not fs.exists(url + ".md5")
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
