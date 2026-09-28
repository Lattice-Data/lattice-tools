"""Streaming endpoints: ``get URL`` to stdout and ``put URL`` from stdin.

Both take any fsspec URL. In production that is ``s3://`` through s3fs, which
picks up IRSA credentials from the pod environment; tests use ``file://`` and a
local moto server reached through :data:`ENDPOINT_ENV`.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys

import fsspec

# 64 MiB multipart parts allow objects up to 640 GiB under S3's 10,000-part cap.
# S3 rejects parts below 5 MiB, so the override exists for tests, not tuning.
DEFAULT_BLOCK = 64 * 2**20
BLOCK_ENV = "FASTQ_CHUNKER_BLOCK"
ENDPOINT_ENV = "FASTQ_CHUNKER_S3_ENDPOINT"


def block_size() -> int:
    return int(os.environ.get(BLOCK_ENV, DEFAULT_BLOCK))


def make_fs(url: str) -> fsspec.AbstractFileSystem:
    protocol = fsspec.utils.get_protocol(url)
    if protocol == "s3":
        kwargs: dict = {"default_block_size": block_size()}
        endpoint = os.environ.get(ENDPOINT_ENV)
        if endpoint:
            kwargs["client_kwargs"] = {"endpoint_url": endpoint}
        return fsspec.filesystem("s3", **kwargs)
    return fsspec.filesystem(protocol, auto_mkdir=True)


def basename(url: str) -> str:
    return url.rstrip("/").rsplit("/", 1)[-1]


def get(url: str, out=None, md5: bool = False) -> str:
    """Stream ``url`` to ``out`` (default stdout); return the MD5 of the bytes."""
    out = out or sys.stdout.buffer
    fs = make_fs(url)
    block = block_size()
    digest = hashlib.md5()
    open_kwargs = {"block_size": block}
    if fs.protocol in ("s3", ("s3", "s3a")):
        open_kwargs["cache_type"] = "readahead"
    with fs.open(url, "rb", **open_kwargs) as f:
        while True:
            b = f.read(block)
            if not b:
                break
            digest.update(b)
            out.write(b)
    out.flush()
    if md5:
        print(f"md5\t{digest.hexdigest()}", file=sys.stderr, flush=True)
    return digest.hexdigest()


def put(
    url: str,
    src=None,
    sidecar: bool = False,
    compress: list[str] | None = None,
) -> tuple[int, str]:
    """Stream ``src`` (default stdin) to ``url``; return (bytes, md5).

    With ``compress``, that command is run with ``src`` on its stdin and its
    stdout is what gets uploaded. Running the compressor here rather than in a
    shell pipe means a compressor that dies mid-stream fails the upload instead
    of leaving a truncated object behind with a success exit code.
    """
    src = src or sys.stdin.buffer
    fs = make_fs(url)
    block = block_size()
    digest = hashlib.md5()
    n = 0
    proc = None
    if compress:
        proc = subprocess.Popen(compress, stdin=src, stdout=subprocess.PIPE)
        src = proc.stdout
    # autocommit=False is fsspec's transaction API: close() uploads the last
    # part but completes nothing, so commit() finishes the object and discard()
    # aborts the multipart upload (or removes the temp file locally)
    f = fs.open(url, "wb", block_size=block, autocommit=False)
    try:
        while True:
            b = src.read(block)
            if not b:
                break
            digest.update(b)
            f.write(b)
            n += len(b)
        if proc is not None:
            proc.stdout.close()
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"compressor {compress[0]} exited {rc}")
    except BaseException:
        try:
            f.close()
        finally:
            f.discard()
        raise
    f.close()
    f.commit()
    if sidecar:
        write_sidecar(fs, url, digest.hexdigest())
    return n, digest.hexdigest()


def sidecar_url(url: str) -> str:
    return url + ".md5"


def write_sidecar(fs: fsspec.AbstractFileSystem, url: str, md5hex: str) -> None:
    with fs.open(sidecar_url(url), "w") as f:
        f.write(f"{md5hex}  {basename(url)}\n")


def read_sidecar(fs: fsspec.AbstractFileSystem, url: str) -> tuple[str, str]:
    """Return (md5hex, basename) from a ``<url>.md5`` sidecar."""
    text = fs.open(sidecar_url(url), "r").read().strip()
    md5hex, _, name = text.partition("  ")
    return md5hex, name


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fastq_chunker.s3io")
    sub = parser.add_subparsers(dest="command", required=True)
    g = sub.add_parser("get")
    g.add_argument("url")
    g.add_argument("--md5", action="store_true", help="print the MD5 to stderr")
    p = sub.add_parser("put")
    p.add_argument("url")
    p.add_argument("--sidecar", action="store_true", help="also write <url>.md5")
    p.add_argument(
        "--compress",
        nargs=argparse.REMAINDER,
        help="compressor command to run on stdin; must come last",
    )
    args = parser.parse_args(argv)
    if args.command == "get":
        get(args.url, md5=args.md5)
    else:
        n, md5hex = put(args.url, sidecar=args.sidecar, compress=args.compress or None)
        # one line the orchestrator parses into the run manifest
        print(f"{basename(args.url)}\t{n}\t{md5hex}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
