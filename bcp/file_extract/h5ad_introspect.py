"""Read observation and feature counts from AnnData (.h5ad) without loading the matrix."""

from __future__ import annotations

from typing import Any


def _n_axis_from_h5(h5: Any, axis: str) -> int:
    """Return n_obs or n_vars from an open AnnData HDF5 file."""
    import h5py

    if axis not in h5:
        raise RuntimeError(f"No '{axis}' group or dataset; not an AnnData h5ad")
    obj = h5[axis]
    if isinstance(obj, h5py.Dataset):
        return int(obj.shape[0])
    if isinstance(obj, h5py.Group):
        for name in ("_index", "index"):
            if name in obj:
                return int(obj[name].shape[0])
        raise RuntimeError(f"AnnData {axis} group has no _index or index")
    raise RuntimeError(f"Unexpected {axis} object; not an AnnData h5ad")


def count_h5ad_dims(bucket: str, key: str) -> tuple[int, int]:
    """Open an h5ad over S3 once and return ``(n_obs, n_vars)``."""
    import fsspec
    import h5py

    uri = f"s3://{bucket}/{key}"
    with fsspec.open(uri, "rb") as fobj:
        with h5py.File(fobj, "r") as h5:
            return _n_axis_from_h5(h5, "obs"), _n_axis_from_h5(h5, "var")


def count_h5ad_observations(bucket: str, key: str) -> int:
    """Open an h5ad over S3 and return the observation count."""
    return count_h5ad_dims(bucket, key)[0]


def count_h5ad_dims_local(path: str) -> tuple[int, int]:
    """Open a local h5ad for testing (same layout as count_h5ad_dims)."""
    import h5py

    with h5py.File(path, "r") as h5:
        return _n_axis_from_h5(h5, "obs"), _n_axis_from_h5(h5, "var")


def count_h5ad_observations_local(path: str) -> int:
    """Open a local h5ad for testing (same layout as count_h5ad_observations)."""
    return count_h5ad_dims_local(path)[0]
