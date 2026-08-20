"""Read observation count from AnnData (.h5ad) without loading the matrix."""

from __future__ import annotations

from typing import Any


def _n_obs_from_h5(h5: Any) -> int:
    """Return n_obs from an open AnnData HDF5 file."""
    import h5py

    if "obs" not in h5:
        raise RuntimeError("No 'obs' group or dataset; not an AnnData h5ad")
    obs = h5["obs"]
    if isinstance(obs, h5py.Dataset):
        return int(obs.shape[0])
    if isinstance(obs, h5py.Group):
        for name in ("_index", "index"):
            if name in obs:
                return int(obs[name].shape[0])
        raise RuntimeError("AnnData obs group has no _index or index")
    raise RuntimeError("Unexpected obs object; not an AnnData h5ad")


def count_h5ad_observations(bucket: str, key: str) -> int:
    """Open an h5ad over S3 and return the observation count."""
    import fsspec
    import h5py

    uri = f"s3://{bucket}/{key}"
    with fsspec.open(uri, "rb") as fobj:
        with h5py.File(fobj, "r") as h5:
            return _n_obs_from_h5(h5)


def count_h5ad_observations_local(path: str) -> int:
    """Open a local h5ad for testing (same layout as count_h5ad_observations)."""
    import h5py

    with h5py.File(path, "r") as h5:
        return _n_obs_from_h5(h5)
