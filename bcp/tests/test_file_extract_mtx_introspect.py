"""Tests for file_extract.mtx_introspect (local 10x-style matrix dirs)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from file_extract.mtx_introspect import count_mtx_observations_local


def _write_mtx(path: Path, *, n_rows: int, n_cols: int, gzipped: bool) -> None:
    body = (
        "%%MatrixMarket matrix coordinate integer general\n"
        f"{n_rows} {n_cols} 2\n"
        "1 1 1\n"
        "2 1 1\n"
    ).encode("utf-8")
    path.write_bytes(gzip.compress(body) if gzipped else body)


def _write_barcodes(path: Path, *, n: int, gzipped: bool) -> None:
    body = "".join(f"bc-{i}\n" for i in range(n)).encode("utf-8")
    path.write_bytes(gzip.compress(body) if gzipped else body)


def test_count_mtx_observations_local_uses_barcodes(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx.gz"
    _write_mtx(mtx, n_rows=10, n_cols=99, gzipped=True)
    _write_barcodes(tmp_path / "barcodes.tsv.gz", n=4, gzipped=True)
    assert count_mtx_observations_local(str(mtx)) == 4


def test_count_mtx_observations_local_uncompressed_barcodes(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx"
    _write_mtx(mtx, n_rows=10, n_cols=99, gzipped=False)
    _write_barcodes(tmp_path / "barcodes.tsv", n=6, gzipped=False)
    assert count_mtx_observations_local(str(mtx)) == 6


def test_count_mtx_observations_local_falls_back_to_header(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx.gz"
    _write_mtx(mtx, n_rows=10, n_cols=7, gzipped=True)
    assert count_mtx_observations_local(str(mtx)) == 7


def test_count_mtx_observations_local_missing_header(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx"
    mtx.write_text(
        "%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="dimensions"):
        count_mtx_observations_local(str(mtx))
