"""Tests for file_extract.h5ad_introspect (local h5ad-shaped HDF5 files)."""

from __future__ import annotations

import pytest

pytest.importorskip("h5py")
import h5py  # noqa: E402
import numpy as np  # noqa: E402

from file_extract.h5ad_introspect import (  # noqa: E402
    count_h5ad_dims_local,
    count_h5ad_observations_local,
)


def _build_axes(
    path: str,
    *,
    n_obs: int,
    n_vars: int,
    obs_index: str = "_index",
    var_index: str = "_index",
) -> None:
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.create_dataset(
            obs_index,
            data=np.array([f"cell-{i}".encode() for i in range(n_obs)], dtype="S"),
        )
        var = h5.create_group("var")
        var.create_dataset(
            var_index,
            data=np.array([f"gene-{i}".encode() for i in range(n_vars)], dtype="S"),
        )


def _build_legacy_obs_dataset(path: str, *, n_obs: int, n_vars: int) -> None:
    with h5py.File(path, "w") as h5:
        h5.create_dataset(
            "obs",
            data=np.zeros(n_obs, dtype=[("index", "S8")]),
        )
        var = h5.create_group("var")
        var.create_dataset(
            "_index",
            data=np.array([f"gene-{i}".encode() for i in range(n_vars)], dtype="S"),
        )


def test_count_h5ad_observations_local_obs_index(tmp_path) -> None:
    path = tmp_path / "current.h5ad"
    _build_axes(str(path), n_obs=3, n_vars=7)
    assert count_h5ad_observations_local(str(path)) == 3


def test_count_h5ad_observations_local_obs_index_alias(tmp_path) -> None:
    path = tmp_path / "index.h5ad"
    _build_axes(str(path), n_obs=5, n_vars=2, obs_index="index")
    assert count_h5ad_observations_local(str(path)) == 5


def test_count_h5ad_observations_local_legacy_dataset(tmp_path) -> None:
    path = tmp_path / "legacy.h5ad"
    _build_legacy_obs_dataset(str(path), n_obs=4, n_vars=9)
    assert count_h5ad_observations_local(str(path)) == 4


def test_count_h5ad_dims_local_var_index(tmp_path) -> None:
    path = tmp_path / "vars.h5ad"
    _build_axes(str(path), n_obs=3, n_vars=12)
    assert count_h5ad_dims_local(str(path)) == (3, 12)


def test_count_h5ad_dims_local_var_index_alias(tmp_path) -> None:
    path = tmp_path / "var_index.h5ad"
    _build_axes(str(path), n_obs=2, n_vars=8, var_index="index")
    assert count_h5ad_dims_local(str(path)) == (2, 8)


def test_count_h5ad_observations_local_missing_obs(tmp_path) -> None:
    path = tmp_path / "empty.h5ad"
    with h5py.File(path, "w"):
        pass
    with pytest.raises(RuntimeError, match="obs"):
        count_h5ad_observations_local(str(path))


def test_count_h5ad_dims_local_missing_var(tmp_path) -> None:
    path = tmp_path / "no_var.h5ad"
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.create_dataset(
            "_index",
            data=np.array([b"cell-0"], dtype="S"),
        )
    with pytest.raises(RuntimeError, match="var"):
        count_h5ad_dims_local(str(path))


def test_count_h5ad_observations_local_obs_group_without_index(tmp_path) -> None:
    path = tmp_path / "no_index.h5ad"
    with h5py.File(path, "w") as h5:
        h5.create_group("obs")
        var = h5.create_group("var")
        var.create_dataset("_index", data=np.array([b"g0"], dtype="S"))
    with pytest.raises(RuntimeError, match="_index"):
        count_h5ad_observations_local(str(path))


def test_count_h5ad_dims_local_uses_attrs_index(tmp_path) -> None:
    path = tmp_path / "attrs.h5ad"
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        obs.create_dataset(
            "cell_id",
            data=np.array([f"c{i}".encode() for i in range(4)], dtype="S"),
        )
        var = h5.create_group("var")
        var.attrs["_index"] = "gene_id"
        var.create_dataset(
            "gene_id",
            data=np.array([f"g{i}".encode() for i in range(9)], dtype="S"),
        )
    assert count_h5ad_dims_local(str(path)) == (4, 9)


def test_count_h5ad_dims_local_attrs_index_bytes(tmp_path) -> None:
    path = tmp_path / "attrs_bytes.h5ad"
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.attrs["_index"] = b"barcode"
        obs.create_dataset(
            "barcode",
            data=np.array([b"bc0", b"bc1"], dtype="S"),
        )
        var = h5.create_group("var")
        var.create_dataset(
            "_index",
            data=np.array([b"g0"], dtype="S"),
        )
    assert count_h5ad_dims_local(str(path)) == (2, 1)


def test_count_h5ad_dims_local_categorical_codes(tmp_path) -> None:
    path = tmp_path / "categorical.h5ad"
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        cell_id = obs.create_group("cell_id")
        cell_id.create_dataset("categories", data=np.array([b"a", b"b"], dtype="S"))
        cell_id.create_dataset("codes", data=np.array([0, 1, 0, 1, 0], dtype="i4"))
        var = h5.create_group("var")
        gene_id = var.create_group("_index")
        gene_id.create_dataset("categories", data=np.array([b"g"], dtype="S"))
        gene_id.create_dataset("codes", data=np.array([0, 0, 0], dtype="i4"))
    assert count_h5ad_dims_local(str(path)) == (5, 3)
