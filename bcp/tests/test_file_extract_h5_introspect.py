"""Tests for file_extract.h5_introspect (local h5 files)."""

from __future__ import annotations

import pytest

pytest.importorskip("h5py")
import h5py  # noqa: E402

from file_extract.h5_introspect import introspect_h5_local  # noqa: E402


def _build_matrix_h5(path: str, *, with_genome: bool = True) -> None:
    import numpy as np

    with h5py.File(path, "w") as h5:
        matrix = h5.create_group("matrix")
        matrix.create_dataset(
            "barcodes",
            data=np.array([b"AAACCTG", b"AAACCTH"], dtype="S"),
        )
        features = matrix.create_group("features")
        types = ["Gene Expression", "Antibody Capture"]
        genomes = ["GRCh38", "GRCh38"]
        features.create_dataset(
            "feature_type",
            data=np.array(types, dtype="S"),
        )
        if with_genome:
            features.create_dataset(
                "genome",
                data=np.array(genomes, dtype="S"),
            )


def test_introspect_h5_local(tmp_path) -> None:
    h5_path = tmp_path / "matrix.h5"
    _build_matrix_h5(str(h5_path))
    obs, type_counts, genome_counts = introspect_h5_local(str(h5_path))
    assert obs == 2
    assert type_counts["Gene Expression"] == 1
    assert type_counts["Antibody Capture"] == 1
    assert genome_counts == {"GRCh38": 1}


def test_introspect_h5_local_invalid(tmp_path) -> None:
    h5_path = tmp_path / "bad.h5"
    with h5py.File(h5_path, "w"):
        pass
    with pytest.raises(RuntimeError, match="matrix"):
        introspect_h5_local(str(h5_path))
