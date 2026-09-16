"""Per-row hashes of a dense ndarray or CSR matrix, computed concurrently.

Returns one digest per row, aligned to row index: out[i] is the hash of row i.

Hashing modes
-------------
canonical=False (default, fast)
    Hashes the raw CSR buffers (or dense row bytes). Requires canonical CSR
    form first -- sorted indices, no duplicate entries, and NO EXPLICIT STORED
    ZEROS. scipy's fancy-index __setitem__ leaves explicit zeros behind, which
    makes two logically identical rows hash differently; row_hashes() calls
    eliminate_zeros() to prevent that (see `copy`, since that mutates).
    Digests are storage-specific: f32 vs f64, or CSR vs dense, will differ.

canonical=True (slower, portable)
    Hashes (nonzero indices as int64, values as float32) per row. Digests
    match across CSR/CSC/dense and ignore explicit stored zeros and index
    order. Use when comparing matrices from different sources.

How workers get the matrix
--------------------------
sharing="inherit"   fork() only. Zero copies: children inherit the buffers
                    copy-on-write and only (lo, hi) tuples cross the boundary.
                    Requires start_method="fork", which CPython deprecates in
                    multi-threaded parents (3.12+) -- and a Jupyter kernel is
                    always multi-threaded. Unsafe if any thread holds a lock at
                    fork time; HDF5 in particular is not fork-safe, so avoid
                    this when adata is backed or h5py handles are open.

sharing="shm"       Works with any start method, including the safe ones.
                    Buffers are written once to /dev/shm (tmpfs on Linux) and
                    workers np.memmap them read-only; only path strings are
                    pickled. Costs one extra copy of the matrix, and on Linux
                    /dev/shm defaults to ~50% of RAM -- check `df -h /dev/shm`.
                    On platforms without /dev/shm (macOS) it falls back to a
                    temp dir, which is disk-backed but page-cached.

Python 3.14 changed the default start method on Linux from fork to forkserver;
this module always passes an explicit context, so behaviour does not shift
under you on upgrade.

NOT thread-safe / not reentrant: state lives in a module-level dict so forked
children inherit buffers instead of pickling them. Do not call row_hashes()
concurrently from multiple threads.
"""

import anndata as ad
import hashlib
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import scipy.sparse as sparse

_G = {}  # module-level so forked children inherit it copy-on-write

_SHM_DIR = "/dev/shm"


# --------------------------------------------------------------------------
# canonicalization
# --------------------------------------------------------------------------


def canonicalize_csr(X, copy=False, warn=True):
    """Put a CSR matrix into canonical form for raw-buffer hashing.

    Returns (matrix, changed). With copy=True the input is never mutated, and
    the copy is only taken when a fix is actually needed -- so an
    already-clean matrix costs nothing.
    """
    needs_dedupe = not X.has_canonical_format
    needs_sort = not X.has_sorted_indices
    # explicit zeros are the failure mode that silently breaks duplicate
    # detection, so always pay the cheap vectorized scan
    needs_zeros = bool((X.data == 0).any())
    changed = needs_dedupe or needs_sort or needs_zeros

    if not changed:
        return X, False

    if copy:
        X = X.copy()
    elif warn:
        warnings.warn(
            "canonicalizing CSR in place (explicit zeros removed / indices "
            "sorted); .nnz and buffer identity will change. Pass copy=True "
            "to leave the input untouched.",
            stacklevel=3,
        )

    if needs_dedupe:
        X.sum_duplicates()
    if needs_sort:
        X.sort_indices()
    if needs_zeros:
        X.eliminate_zeros()
    return X, True


# --------------------------------------------------------------------------
# sharing
# --------------------------------------------------------------------------


def shm_dir_available():
    """True if /dev/shm exists and is writable (Linux tmpfs)."""
    return os.path.isdir(_SHM_DIR) and os.access(_SHM_DIR, os.W_OK)


def shm_free_bytes():
    """Free bytes on the shm filesystem, or None if unavailable."""
    if not shm_dir_available():
        return None
    st = os.statvfs(_SHM_DIR)
    return st.f_bavail * st.f_frsize


def _export(arrays, nbytes_total):
    """Write arrays to a shm-backed dir. Returns (dirpath, {name: (path, dtype, shape)})."""
    if shm_dir_available():
        free = shm_free_bytes()
        if free is not None and nbytes_total > free:
            warnings.warn(
                f"{nbytes_total / 1e9:.1f} GB of buffers exceeds {free / 1e9:.1f} GB "
                f"free on {_SHM_DIR}; falling back to a temp dir. Remount "
                f"{_SHM_DIR} larger, or use sharing='inherit'.",
                stacklevel=3,
            )
            parent = None
        else:
            parent = _SHM_DIR
    else:
        parent = None

    dirpath = tempfile.mkdtemp(prefix="row_hashes_", dir=parent)
    spec = {}
    for name, arr in arrays.items():
        path = os.path.join(dirpath, f"{name}.bin")
        arr.tofile(path)  # raw contiguous bytes
        spec[name] = (path, arr.dtype.str, arr.shape)
    return dirpath, spec


def _attach(spec):
    """Reopen exported arrays read-only in a worker."""
    return {
        name: np.memmap(path, dtype=np.dtype(dt), mode="r", shape=tuple(shape))
        for name, (path, dt, shape) in spec.items()
    }


def _init_worker(state):
    """ProcessPoolExecutor initializer for sharing='shm'."""
    _G.clear()
    _G.update(state)
    arrays = _attach(state["spec"])
    _install_arrays(arrays, state["kind"], state["canonical"])


def _install_arrays(arrays, kind, canonical):
    """Populate _G's buffer views from typed arrays (parent or worker side)."""
    if kind == "dense":
        _G["A"] = arrays["A"]
        return
    _G["indptr"] = arrays["indptr"].tolist()  # list indexing beats np scalars
    if canonical:
        _G["data_arr"] = arrays["data"]
        _G["indices_arr"] = arrays["indices"]
    else:
        _G["itemsizes"] = (
            arrays["data"].dtype.itemsize,
            arrays["indices"].dtype.itemsize,
        )
        _G["data"] = memoryview(np.asarray(arrays["data"])).cast("B")
        _G["indices"] = memoryview(np.asarray(arrays["indices"])).cast("B")


# --------------------------------------------------------------------------
# hashing
# --------------------------------------------------------------------------


def _chunk(span):
    """Hash rows [lo, hi) and return their digests in row order."""
    lo, hi = span
    meta = _G["meta"]
    canonical = _G["canonical"]
    out = []

    if _G["kind"] == "dense":
        A = _G["A"]
        if canonical:
            for i in range(lo, hi):
                row = A[i]
                nz = np.flatnonzero(row)  # != 0, so -0.0 drops, NaN stays
                h = hashlib.sha256(meta)
                h.update(f"|{len(nz)}".encode())
                h.update(np.ascontiguousarray(nz, dtype=np.int64))
                h.update(np.ascontiguousarray(row[nz], dtype=np.float32))
                out.append(h.hexdigest())
        else:
            for i in range(lo, hi):
                h = hashlib.sha256(meta)
                h.update(np.ascontiguousarray(A[i]))  # memmap rows need a copy
                out.append(h.hexdigest())
        return out

    ptr = _G["indptr"]
    if canonical:
        data, idx = _G["data_arr"], _G["indices_arr"]
        for i in range(lo, hi):
            s, e = ptr[i], ptr[i + 1]
            ix, vs = idx[s:e], data[s:e]
            keep = vs != 0  # no mutation of the input needed
            ix = np.ascontiguousarray(ix[keep], dtype=np.int64)
            vs = np.ascontiguousarray(vs[keep], dtype=np.float32)
            order = np.argsort(ix, kind="stable")
            h = hashlib.sha256(meta)
            h.update(f"|{len(ix)}".encode())
            h.update(np.ascontiguousarray(ix[order]))
            h.update(np.ascontiguousarray(vs[order]))
            out.append(h.hexdigest())
    else:
        d, ix, (dw, iw) = _G["data"], _G["indices"], _G["itemsizes"]
        for i in range(lo, hi):
            s, e = ptr[i], ptr[i + 1]
            h = hashlib.sha256(meta)
            h.update(ix[s * iw : e * iw])
            h.update(d[s * dw : e * dw])
            out.append(h.hexdigest())
    return out


def default_start_method():
    """fork where it is the safe/only useful option, else forkserver."""
    available = mp.get_all_start_methods()
    if sys.platform == "darwin":
        # no /dev/shm on macOS, so inherit-via-fork is the only zero-copy path
        return "fork" if "fork" in available else "spawn"
    if shm_dir_available() and "forkserver" in available:
        return "forkserver"
    return "fork" if "fork" in available else available[0]


def _prepare(X, canonical, copy):
    """Normalize input into (kind, meta, arrays, bytes_per_row)."""
    n_rows = X.shape[0]
    if sparse.issparse(X):
        X = X.tocsr()
        if not canonical:
            # raw-buffer mode compares storage, so storage must be canonical
            X, _ = canonicalize_csr(X, copy=copy)
        meta = (
            f"canon|{X.shape[1]}".encode()
            if canonical
            else (f"csr|{X.data.dtype.str}|{X.indices.dtype.str}|{X.shape[1]}").encode()
        )
        arrays = {"data": X.data, "indices": X.indices, "indptr": X.indptr}
        return "csr", meta, arrays, X.data.nbytes / max(n_rows, 1)

    A = np.ascontiguousarray(X)  # strided rows would be rejected
    meta = (
        f"canon|{A.shape[1]}".encode()
        if canonical
        else f"dense|{A.dtype.str}|{A.shape[1]}".encode()
    )
    return "dense", meta, {"A": A}, A.shape[1] * A.itemsize


def _resolve(start_method, sharing):
    if start_method == "auto":
        start_method = default_start_method()
    if sharing == "auto":
        sharing = "inherit" if start_method == "fork" else "shm"
    if sharing == "inherit" and start_method != "fork":
        raise ValueError(
            f"sharing='inherit' requires start_method='fork', got "
            f"{start_method!r}. Use sharing='shm' with forkserver/spawn."
        )
    return start_method, sharing


def _get_workers():
    """Might as well use max number of workers"""
    return os.cpu_count()


class RowHasher:
    """Persistent worker pool for repeated per-row hashing of one matrix.

    Pool spin-up (~0.15s for forkserver/spawn) and the shm export are paid
    once here instead of per call. Measured on 8-core Linux, 200k x 30k CSR:
    one-shot forkserver+shm ran at 0.84x of serial, while the same work
    through a warm RowHasher ran at 5.6x. Use this for anything repeated;
    use row_hashes() for a single pass.

        with RowHasher(adata.X, n_workers=8) as rhr:
            first = rhr.hashes()
            again = rhr.hashes()      # no setup cost

    With sharing="inherit" the children snapshot the buffers copy-on-write at
    pool creation; mutating X afterwards will NOT be reflected. Re-create the
    RowHasher if the matrix changes.
    """

    def __init__(
        self,
        X,
        n_workers=None,
        batch=5000,
        canonical=False,
        copy=False,
        start_method="auto",
        sharing="auto",
    ):
        self.n_rows = X.shape[0]
        self.batch = batch
        self.canonical = canonical
        self.n_workers = n_workers if n_workers else _get_workers()
        self.start_method, self.sharing = _resolve(start_method, sharing)
        self.kind, self.meta, self._arrays, self.bytes_per_row = _prepare(
            X, canonical, copy
        )
        self._dirpath = None
        self._ex = None

    def _spans(self, batch=None):
        b = batch or self.batch
        return [(i, min(i + b, self.n_rows)) for i in range(0, self.n_rows, b)]

    def __enter__(self):
        _G.clear()
        _G["canonical"] = self.canonical
        _G["kind"] = self.kind
        _G["meta"] = self.meta
        ctx = mp.get_context(self.start_method)

        if self.sharing == "inherit":
            # install before forking so children inherit the buffers
            _install_arrays(self._arrays, self.kind, self.canonical)
            self._ex = ProcessPoolExecutor(self.n_workers, mp_context=ctx)
        else:
            total = sum(a.nbytes for a in self._arrays.values())
            self._dirpath, spec = _export(self._arrays, total)
            state = {
                "kind": self.kind,
                "meta": self.meta,
                "canonical": self.canonical,
                "spec": spec,
            }
            _install_arrays(self._arrays, self.kind, self.canonical)
            self._ex = ProcessPoolExecutor(
                self.n_workers,
                mp_context=ctx,
                initializer=_init_worker,
                initargs=(state,),
            )
        return self

    def __exit__(self, *exc):
        if self._ex is not None:
            self._ex.shutdown()
            self._ex = None
        if self._dirpath is not None:
            shutil.rmtree(self._dirpath, ignore_errors=True)
            self._dirpath = None
        return False

    def hashes(self, batch=None):
        """One digest per row, row-aligned. Requires an active context."""
        if self._ex is None:
            raise RuntimeError("use RowHasher as a context manager")
        return [h for c in self._ex.map(_chunk, self._spans(batch)) for h in c]

    def duplicates(self, batch=None):
        return _group(self.hashes(batch))


def row_hashes(
    X,
    n_workers=None,
    batch=5000,
    backend="auto",
    canonical=False,
    copy=False,
    start_method="auto",
    sharing="auto",
):
    """One sha256 hex digest per row of X (dense ndarray or scipy sparse).

    Single-pass convenience wrapper; for repeated passes use RowHasher, which
    keeps the pool and the shm export warm between calls.

    backend: "process", "thread" (capped at 2; best for fat dense rows),
        "serial", or "auto".
    start_method: "fork" | "forkserver" | "spawn" | "auto".
    sharing: "inherit" (fork only, zero-copy) | "shm" (any start method, one
        extra copy) | "auto".
    canonical, copy: see module docstring.
    """
    n_rows = X.shape[0]
    kind, meta, arrays, bytes_per_row = _prepare(X, canonical, copy)
    n_workers = n_workers if n_workers else _get_workers()

    if backend == "auto":
        backend = "thread" if bytes_per_row >= 16_384 else "process"
        if n_rows < 20_000:
            backend = "serial"
        elif canonical:
            backend = "process"  # allocation-heavy per row; threads do not help

    if backend in ("serial", "thread"):
        _G.clear()
        _G.update(canonical=canonical, kind=kind, meta=meta)
        _install_arrays(arrays, kind, canonical)
        spans = [(i, min(i + batch, n_rows)) for i in range(0, n_rows, batch)]
        if backend == "serial":
            chunks = [_chunk(s) for s in spans]
        else:
            with ThreadPoolExecutor(n_workers) as ex:
                chunks = list(ex.map(_chunk, spans))  # map preserves order
        return [h for c in chunks for h in c]

    # process backend: reuse RowHasher so there is one code path
    hasher = RowHasher.__new__(RowHasher)
    hasher.n_rows, hasher.batch, hasher.canonical = n_rows, batch, canonical
    hasher.n_workers = n_workers
    hasher.start_method, hasher.sharing = _resolve(start_method, sharing)
    hasher.kind, hasher.meta = kind, meta
    hasher._arrays, hasher.bytes_per_row = arrays, bytes_per_row
    hasher._dirpath = hasher._ex = None
    with hasher:
        return hasher.hashes()  # out[i] == hash of row i


def _group(digests):
    groups = {}
    for i, h in enumerate(digests):
        groups.setdefault(h, []).append(i)
    return {h: rows for h, rows in groups.items() if len(rows) > 1}


def find_duplicate_rows(X, **kw):
    """Map digest -> sorted list of row indices, for groups of size > 1."""
    return _group(row_hashes(X, **kw))


def evaluate_dup_counts(adata: ad.AnnData, **kw):
    """
    Entry point for row hashing
    Returns dataframe with duplicates or None

    Filters in_tissue==0 cells from spatial datasets (can result in false positive
    duplicate results)
    """
    if 'in_tissue' in adata.obs.columns:
        obs_to_keep = adata.obs[adata.obs['in_tissue'] != 0].index
        adata = adata[obs_to_keep, : ]

    matrix = adata.raw.X if adata.raw else adata.X

    hashes = row_hashes(matrix, **kw)
    
    hash_df = adata.obs.copy()
    hash_df['row_hash'] = hashes
    hash_df = hash_df[hash_df.duplicated(subset='row_hash',keep=False) == True]
    hash_df.sort_values('row_hash', inplace=True)

    if not hash_df.empty:
        print('duplicated raw counts', 'ERROR')
        return hash_df
    print('no duplicated raw counts', 'GOOD')
