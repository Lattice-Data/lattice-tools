"""Round trips for the streaming endpoints over file:// and a local moto server."""

from __future__ import annotations

import hashlib
import io
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
    assert (tmp_path / "out" / "x.bin.md5").read_text() == f"{md5hex}  x.bin\n"

    got = run_endpoint(["get", url, "--md5"], env={s3io.BLOCK_ENV: str(MIB)})
    assert got.stdout == data
    assert got.stderr.decode().strip().splitlines()[-1] == f"md5\t{md5hex}"


def test_put_without_sidecar(tmp_path):
    url = (tmp_path / "y.bin").as_uri()
    run_endpoint(["put", url, "--no-sidecar"], stdin=b"abc")
    assert (tmp_path / "y.bin").read_bytes() == b"abc"
    assert not (tmp_path / "y.bin.md5").exists()


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
    assert s3io.read_sidecar(fs, url) == (md5hex, "big.bin")
    got = run_endpoint(["get", url, "--md5"], env=env)
    assert got.stdout == data
    assert got.stderr.decode().strip().splitlines()[-1] == f"md5\t{md5hex}"
