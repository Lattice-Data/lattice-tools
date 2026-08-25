"""Tests for file_extract.mtx_introspect (local 10x-style matrix dirs)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from file_extract.mtx_introspect import (
    count_mtx_dims,
    count_mtx_dims_local,
    count_mtx_observations,
    count_mtx_observations_local,
)
from tests.file_extract_helpers import MockS3Client


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


def _write_features(path: Path, *, n: int, gzipped: bool) -> None:
    body = "".join(f"feat-{i}\n" for i in range(n)).encode("utf-8")
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


def test_count_mtx_dims_local_uses_features(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx.gz"
    _write_mtx(mtx, n_rows=99, n_cols=7, gzipped=True)
    _write_barcodes(tmp_path / "barcodes.tsv.gz", n=4, gzipped=True)
    _write_features(tmp_path / "features.tsv.gz", n=5, gzipped=True)
    assert count_mtx_dims_local(str(mtx)) == (4, 5)


def test_count_mtx_dims_local_uncompressed_features(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx"
    _write_mtx(mtx, n_rows=99, n_cols=7, gzipped=False)
    _write_features(tmp_path / "features.tsv", n=6, gzipped=False)
    assert count_mtx_dims_local(str(mtx)) == (7, 6)


def test_count_mtx_dims_local_falls_back_to_header_first_dim(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx.gz"
    _write_mtx(mtx, n_rows=10, n_cols=7, gzipped=True)
    assert count_mtx_dims_local(str(mtx)) == (7, 10)


def test_count_mtx_observations_s3_uses_barcodes() -> None:
    mtx = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx.gz"
    barcodes = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/barcodes.tsv.gz"
    client = MockS3Client(
        object_bodies={
            mtx: gzip.compress(
                b"%%MatrixMarket matrix coordinate integer general\n10 99 1\n"
            ),
            barcodes: gzip.compress(b"a\nb\nc\n"),
        }
    )
    assert count_mtx_observations(client, "bucket", mtx) == 3


def test_count_mtx_dims_s3_uses_features() -> None:
    mtx = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx.gz"
    barcodes = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/barcodes.tsv.gz"
    features = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/features.tsv.gz"
    client = MockS3Client(
        object_bodies={
            mtx: gzip.compress(
                b"%%MatrixMarket matrix coordinate integer general\n99 8 1\n"
            ),
            barcodes: gzip.compress(b"a\nb\nc\n"),
            features: gzip.compress(b"h1\nh2\nh3\nh4\n"),
        }
    )
    assert count_mtx_dims(client, "bucket", mtx) == (3, 4)


def test_count_mtx_observations_s3_falls_back_to_header() -> None:
    mtx = "run/scaleplex/S.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx.gz"
    client = MockS3Client(
        object_bodies={
            mtx: gzip.compress(
                b"%%MatrixMarket matrix coordinate integer general\n5 8 1\n1 1 1\n"
            )
        }
    )
    assert count_mtx_observations(client, "bucket", mtx) == 8
    assert count_mtx_dims(client, "bucket", mtx) == (8, 5)


def test_count_mtx_observations_local_missing_header(tmp_path: Path) -> None:
    mtx = tmp_path / "matrix.mtx"
    mtx.write_text(
        "%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="dimensions"):
        count_mtx_observations_local(str(mtx))
