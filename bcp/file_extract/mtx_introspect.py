"""Read observation and feature counts from a 10x-style Market Matrix directory."""

from __future__ import annotations

import gzip
import io
import zlib
from pathlib import Path
from typing import Any

from .s3_utils import get_object_bytes

BARCODE_BASENAMES = ("barcodes.tsv.gz", "barcodes.tsv")
FEATURE_BASENAMES = ("features.tsv.gz", "features.tsv")
MTX_HEADER_BYTES = 65536


def _count_tsv_lines(data: bytes, *, gzipped: bool) -> int:
    text = gzip.decompress(data).decode("utf-8") if gzipped else data.decode("utf-8")
    return sum(1 for line in text.splitlines() if line.strip())


def _mtx_header_dims(data: bytes, *, gzipped: bool) -> tuple[int, int]:
    """Return ``(n_rows, n_cols)`` from a Matrix Market header (features × cells)."""
    if gzipped:
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        raw = decoder.decompress(data)
        if not raw:
            raw = gzip.decompress(data)
    else:
        raw = data
    for line in io.TextIOWrapper(io.BytesIO(raw), encoding="utf-8"):
        stripped = line.strip()
        if not stripped or stripped.startswith("%"):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            raise RuntimeError("MTX header missing dimensions")
        return int(parts[0]), int(parts[1])
    raise RuntimeError("MTX header missing dimensions")


def _get_object_prefix(s3_client: Any, bucket: str, key: str) -> bytes:
    """Read the start of an object so MTX headers do not require a full download."""
    try:
        resp = s3_client.get_object(
            Bucket=bucket, Key=key, Range=f"bytes=0-{MTX_HEADER_BYTES - 1}"
        )
        return resp["Body"].read()
    except Exception:
        return get_object_bytes(s3_client, bucket, key)


def _count_sibling_lines(
    s3_client: Any, bucket: str, parent: str, basenames: tuple[str, ...]
) -> int | None:
    for name in basenames:
        sibling = f"{parent}/{name}"
        try:
            data = get_object_bytes(s3_client, bucket, sibling)
        except Exception:
            continue
        return _count_tsv_lines(data, gzipped=name.endswith(".gz"))
    return None


def _count_local_sibling_lines(parent: Path, basenames: tuple[str, ...]) -> int | None:
    for name in basenames:
        sibling = parent / name
        if sibling.is_file():
            return _count_tsv_lines(sibling.read_bytes(), gzipped=name.endswith(".gz"))
    return None


def count_mtx_dims(s3_client: Any, bucket: str, mtx_key: str) -> tuple[int, int]:
    """Return ``(n_obs, n_features)`` from siblings, else the MTX header."""
    parent = mtx_key.rsplit("/", 1)[0]
    n_obs = _count_sibling_lines(s3_client, bucket, parent, BARCODE_BASENAMES)
    n_features = _count_sibling_lines(s3_client, bucket, parent, FEATURE_BASENAMES)
    if n_obs is None or n_features is None:
        data = _get_object_prefix(s3_client, bucket, mtx_key)
        n_rows, n_cols = _mtx_header_dims(data, gzipped=mtx_key.endswith(".gz"))
        if n_obs is None:
            n_obs = n_cols
        if n_features is None:
            n_features = n_rows
    return n_obs, n_features


def count_mtx_observations(s3_client: Any, bucket: str, mtx_key: str) -> int:
    """Return n_obs from sibling barcodes, else the MTX header column count."""
    return count_mtx_dims(s3_client, bucket, mtx_key)[0]


def count_mtx_dims_local(path: str) -> tuple[int, int]:
    """Open a local matrix.mtx(.gz) for testing (same layout as S3)."""
    mtx_path = Path(path)
    n_obs = _count_local_sibling_lines(mtx_path.parent, BARCODE_BASENAMES)
    n_features = _count_local_sibling_lines(mtx_path.parent, FEATURE_BASENAMES)
    if n_obs is None or n_features is None:
        n_rows, n_cols = _mtx_header_dims(
            mtx_path.read_bytes(), gzipped=mtx_path.name.endswith(".gz")
        )
        if n_obs is None:
            n_obs = n_cols
        if n_features is None:
            n_features = n_rows
    return n_obs, n_features


def count_mtx_observations_local(path: str) -> int:
    """Open a local matrix.mtx(.gz) for testing (same layout as S3)."""
    return count_mtx_dims_local(path)[0]
