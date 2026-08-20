"""Tests for file_extract.h5ad_introspect (local h5ad-shaped HDF5 files)."""

from __future__ import annotations

import pytest

pytest.importorskip("h5py")
import h5py  # noqa: E402
import numpy as np  # noqa: E402

from file_extract.h5ad_introspect import count_h5ad_observations_local  # noqa: E402


def _build_obs_group(path: str, *, n: int, index_name: str) -> None:
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        obs.create_dataset(
            index_name,
            data=np.array([f"cell-{i}".encode() for i in range(n)], dtype="S"),
        )


def _build_legacy_obs_dataset(path: str, *, n: int) -> None:
    with h5py.File(path, "w") as h5:
        h5.create_dataset(
            "obs",
            data=np.zeros(n, dtype=[("index", "S8")]),
        )


def test_count_h5ad_observations_local_obs_index(tmp_path) -> None:
    path = tmp_path / "current.h5ad"
    _build_obs_group(str(path), n=3, index_name="_index")
    assert count_h5ad_observations_local(str(path)) == 3


def test_count_h5ad_observations_local_obs_index_alias(tmp_path) -> None:
    path = tmp_path / "index.h5ad"
    _build_obs_group(str(path), n=5, index_name="index")
    assert count_h5ad_observations_local(str(path)) == 5


def test_count_h5ad_observations_local_legacy_dataset(tmp_path) -> None:
    path = tmp_path / "legacy.h5ad"
    _build_legacy_obs_dataset(str(path), n=4)
    assert count_h5ad_observations_local(str(path)) == 4


def test_count_h5ad_observations_local_missing_obs(tmp_path) -> None:
    path = tmp_path / "empty.h5ad"
    with h5py.File(path, "w"):
        pass
    with pytest.raises(RuntimeError, match="obs"):
        count_h5ad_observations_local(str(path))


def test_count_h5ad_observations_local_obs_group_without_index(tmp_path) -> None:
    path = tmp_path / "no_index.h5ad"
    with h5py.File(path, "w") as h5:
        h5.create_group("obs")
    with pytest.raises(RuntimeError, match="_index"):
        count_h5ad_observations_local(str(path))
