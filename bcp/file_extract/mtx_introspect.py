"""Read observation count from a 10x-style Market Matrix directory."""

from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Any

from .s3_utils import get_object_bytes

BARCODE_BASENAMES = ("barcodes.tsv.gz", "barcodes.tsv")


def _count_barcode_lines(data: bytes, *, gzipped: bool) -> int:
    text = gzip.decompress(data).decode("utf-8") if gzipped else data.decode("utf-8")
    return sum(1 for line in text.splitlines() if line.strip())


def _n_obs_from_mtx_header(data: bytes, *, gzipped: bool) -> int:
    raw = gzip.decompress(data) if gzipped else data
    for line in io.TextIOWrapper(io.BytesIO(raw), encoding="utf-8"):
        stripped = line.strip()
        if not stripped or stripped.startswith("%"):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            raise RuntimeError("MTX header missing dimensions")
        return int(parts[1])
    raise RuntimeError("MTX header missing dimensions")


def count_mtx_observations(s3_client: Any, bucket: str, mtx_key: str) -> int:
    """Return n_obs from sibling barcodes, else the MTX header column count."""
    parent = mtx_key.rsplit("/", 1)[0]
    for name in BARCODE_BASENAMES:
        sibling = f"{parent}/{name}"
        try:
            data = get_object_bytes(s3_client, bucket, sibling)
        except Exception:
            continue
        return _count_barcode_lines(data, gzipped=name.endswith(".gz"))
    data = get_object_bytes(s3_client, bucket, mtx_key)
    return _n_obs_from_mtx_header(data, gzipped=mtx_key.endswith(".gz"))


def count_mtx_observations_local(path: str) -> int:
    """Open a local matrix.mtx(.gz) for testing (same layout as S3)."""
    mtx_path = Path(path)
    for name in BARCODE_BASENAMES:
        sibling = mtx_path.parent / name
        if sibling.is_file():
            return _count_barcode_lines(
                sibling.read_bytes(), gzipped=name.endswith(".gz")
            )
    return _n_obs_from_mtx_header(
        mtx_path.read_bytes(), gzipped=mtx_path.name.endswith(".gz")
    )
