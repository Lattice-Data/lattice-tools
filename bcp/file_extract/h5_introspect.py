from __future__ import annotations

from collections import Counter


def check_introspection_deps() -> None:
    """Import optional deps once, with an actionable error if missing."""
    try:
        import fsspec  # noqa: F401
        import h5py  # noqa: F401
        import s3fs  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            f"Introspection needs h5py + fsspec + s3fs (missing: {e.name}). "
            "Install them, or pass --no-introspect for a checksum-only run."
        )


def introspect_h5(
    bucket: str, key: str
) -> tuple[int, dict[str, int], dict[str, int] | None]:
    """
    Open the matrix h5 over S3 and return (observation_count, feature_type_counts,
    genome_gene_counts). Reads only headers + the small features table.
    """
    import h5py
    import fsspec

    uri = f"s3://{bucket}/{key}"
    with fsspec.open(uri, "rb") as fobj:
        with h5py.File(fobj, "r") as h5:
            if "matrix" not in h5:
                raise RuntimeError(
                    "No 'matrix' group; not a Cell Ranger feature-barcode matrix"
                )
            matrix = h5["matrix"]

            if "barcodes" not in matrix:
                raise RuntimeError("No 'matrix/barcodes' dataset")
            obs_count = int(matrix["barcodes"].shape[0])

            type_counts: dict[str, int] = {}
            types: list[str] | None = None
            if "features/feature_type" in matrix:
                types = [
                    t.decode("utf-8") if isinstance(t, (bytes, bytearray)) else str(t)
                    for t in matrix["features/feature_type"][:]
                ]
                type_counts = dict(Counter(types))

            genome_counts: dict[str, int] | None = None
            if types is not None and "features/genome" in matrix:
                genomes = [
                    g.decode("utf-8") if isinstance(g, (bytes, bytearray)) else str(g)
                    for g in matrix["features/genome"][:]
                ]
                gx_genomes = [
                    gen for gen, typ in zip(genomes, types) if typ == "Gene Expression"
                ]
                genome_counts = dict(Counter(gx_genomes))

            return obs_count, type_counts, genome_counts


def introspect_h5_local(path: str) -> tuple[int, dict[str, int], dict[str, int] | None]:
    """Open a local h5 file for testing (same layout as introspect_h5)."""
    import h5py

    with h5py.File(path, "r") as h5:
        if "matrix" not in h5:
            raise RuntimeError(
                "No 'matrix' group; not a Cell Ranger feature-barcode matrix"
            )
        matrix = h5["matrix"]
        if "barcodes" not in matrix:
            raise RuntimeError("No 'matrix/barcodes' dataset")
        obs_count = int(matrix["barcodes"].shape[0])

        type_counts: dict[str, int] = {}
        types: list[str] | None = None
        if "features/feature_type" in matrix:
            types = [
                t.decode("utf-8") if isinstance(t, (bytes, bytearray)) else str(t)
                for t in matrix["features/feature_type"][:]
            ]
            type_counts = dict(Counter(types))

        genome_counts: dict[str, int] | None = None
        if types is not None and "features/genome" in matrix:
            genomes = [
                g.decode("utf-8") if isinstance(g, (bytes, bytearray)) else str(g)
                for g in matrix["features/genome"][:]
            ]
            gx_genomes = [
                gen for gen, typ in zip(genomes, types) if typ == "Gene Expression"
            ]
            genome_counts = dict(Counter(gx_genomes))

        return obs_count, type_counts, genome_counts
