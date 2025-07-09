import anndata as ad
import dask
import hashlib
import h5py
import numpy as np
import os
import pandas as pd
import psutil
from abc import ABC, abstractmethod
from anndata.compat import DaskArray
from dataclasses import dataclass
from dask.array import map_blocks
from scipy import sparse
from typing import Union, Callable
from cellxgene_mods import report, Sizes
from cellxgene_schema.utils import SPARSE_MATRIX_TYPES, read_h5ad


@dataclass
class IndicesResult:
    index: int
    barcode: str
    matrix_slice: sparse.csr_matrix
    data_array = None
    indices_array = None


class BaseMatrixFunctions(ABC):
    @staticmethod
    def get_matrix_functions(adata: ad.AnnData):
        if isinstance(adata.X, DaskArray):
            return DaskMatrixFunctions()
        return NonDaskMatrixFunctions()

    @abstractmethod
    def evaluate_sparsity(self, adata: ad.AnnData):
        pass

    @abstractmethod
    def evaluate_data(self, adata: ad.AnnData):
        pass

    @abstractmethod
    def evaluate_dup_counts(self, adata: ad.AnnData) -> pd.DataFrame | None:
        pass

    @abstractmethod
    def evaluate_all_zero_indices(self, adata: ad.AnnData):
        pass


class DaskMatrixFunctions(BaseMatrixFunctions):
    def _get_matrix_format(self, matrix: DaskArray) -> str:
        """
        Given a matrix, returns the format as one of: csc, csr, coo, dense
        or unknown.

        This mimics the scipy.sparse `format` property, but extends it to
        support ndarray and other classes AnnData may proxy the matrix with.
        """

        # Note: the AnnData proxy classes DO support the `format_str` property, but
        # doing a slice seemed safer, if less performant.  Using `format_str`, which
        # currently works, uses private API:
        #
        # >>> return getattr(matrix, "format_str", "dense)
        #
        matrix_format = "unknown"
        try:
            matrix_slice = matrix[0:1, 0:1].compute()
        except AttributeError:
            # compute() may fail on an unknown matrix value. if so, return "unknown"
            return matrix_format
        if isinstance(matrix_slice, sparse.spmatrix):
            matrix_format = matrix_slice.format
        elif isinstance(matrix_slice, np.ndarray):
            matrix_format = "dense"
        assert matrix_format in SPARSE_MATRIX_TYPES.union({"unknown", "dense"})
        return matrix_format

    def _count_matrix_nonzero(self, matrix: DaskArray) -> int:
        def count_nonzeros(
            matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool
        ) -> np.array:
            nnz = (
                matrix_chunk.nnz if is_sparse_matrix else np.count_nonzero(matrix_chunk)
            )
            return np.array([nnz])

        is_sparse_matrix = self._get_matrix_format(matrix) in SPARSE_MATRIX_TYPES
        # if matrix too small to chunk, then just call function and don't dask distribute
        if len(matrix.chunks[0]) > 1:
            nonzeros = (
                map_blocks(
                    count_nonzeros, matrix, is_sparse_matrix, drop_axis=1, dtype=int
                )
                .compute(scheduler="processes")
                .sum()
            )
        else:
            nonzeros = count_nonzeros(matrix.compute(), is_sparse_matrix)[0]
        return nonzeros

    def _determine_sparsity(self, matrix: DaskArray):
        nnz = self._count_matrix_nonzero(matrix)
        sparsity = 1 - nnz / np.prod(matrix.shape)
        return round(sparsity, 3)

    def _get_matrices_to_evaluate(
        self,
        adata: ad.AnnData,
    ) -> list[tuple[Union[np.ndarray, sparse.spmatrix], str]]:
        """
        Helper function to return list of tuples of matrix location and matrix name

        Other functions use this to iterate and report matrix name as error/warning message
        """
        # will always be adata.X, then add other matrices if they exist
        matrices_to_evaluate = [(adata.X, "X")]

        if adata.raw:
            matrices_to_evaluate.append((adata.raw.X, "raw.X"))
        for layer in adata.layers:
            matrices_to_evaluate.append((adata.layers[layer], f"layer '{layer}'"))

        return matrices_to_evaluate

    def evaluate_sparsity(self, adata: ad.AnnData):
        """
        Evaluate sparsity in qa notebook

        Expects matrices as dask array, so make sure to load Anndata with
        read_h5ad()
        """
        max_sparsity = 0.5
        valid = True

        matrices_to_evaluate = self._get_matrices_to_evaluate(adata)

        for matrix, matrix_name in matrices_to_evaluate:
            format = self._get_matrix_format(matrix)
            sparsity = self._determine_sparsity(matrix)
            report(f"{matrix_name} sparsity: {sparsity}")
            if sparsity > max_sparsity and format not in SPARSE_MATRIX_TYPES:
                report(
                    f"{matrix_name} should be converted to csr sparse, found to be {format}",
                    "ERROR",
                )
                valid = False

        if valid:
            report("All matrices have passed checks", "GOOD")

    def evaluate_data(self, adata: ad.AnnData):
        """
        3 other data checks on matrix: min, max, and if raw matrix, all integer check

        For better efficency, all-integer check only done on raw matrix
        Probably room for improvement to delay chunks, better call compute(), or
        other dask optimizations to not chunk through each matrix 2-3 times
        """
        min_maxs = {}
        matrices_to_evaluate = self._get_matrices_to_evaluate(adata)
        matrix_names = [name for _, name in matrices_to_evaluate]
        raw_matrix = "raw.X" if "raw.X" in matrix_names else "X"

        def get_max_chunk(matrix_chunk: sparse.spmatrix) -> np.array:
            return np.array([matrix_chunk.max()]).reshape(-1, 1)

        def get_min_chunk(matrix_chunk: sparse.spmatrix) -> np.array:
            return np.array([matrix_chunk.min()]).reshape(-1, 1)

        def all_integers_chunk(matrix_chunk: sparse.spmatrix) -> np.array:
            data_array = matrix_chunk.data
            is_all_integers = np.all(np.round(data_array) == data_array)
            return np.array([is_all_integers]).reshape(-1, 1)

        for matrix, matrix_name in matrices_to_evaluate:
            with dask.config.set(scheduler="processes"):
                max = (
                    map_blocks(get_max_chunk, matrix, dtype=int).compute().ravel().max()
                )
                min = (
                    map_blocks(get_min_chunk, matrix, dtype=int).compute().ravel().min()
                )
                if matrix_name == raw_matrix:
                    is_all_integers = (
                        map_blocks(all_integers_chunk, matrix, dtype=bool)
                        .compute()
                        .ravel()
                        .all()
                    )

            min_maxs[matrix_name] = f"{min}-{max}"

            report(f"{matrix_name} min = {min}")
            report(f"{matrix_name} max = {max}")

            if matrix_name == raw_matrix:
                if is_all_integers:
                    report(f"{matrix_name} is all integers", "GOOD")
                else:
                    report(f"{matrix_name} is NOT all integers", "ERROR")

        poss_dups = [
            k for k, v in min_maxs.items() if list(min_maxs.values()).count(v) > 1
        ]

        if poss_dups:
            report(f"possible redundant layers: {poss_dups}", "WARNING")

    def evaluate_all_zero_indices(self, adata: ad.AnnData, worker_type="processes"):
        """
        Function to check if a row/cell contains an all-zero indices array. This can exist
        in Visium datasets with spots that are in_tissue == 0, but should not exist for
        in_tissue == 1 or other single cell raw count data.
        Uses dask to lazily load the raw matrix and to check in parallel
        """
        matrix = adata.raw.X if adata.raw else adata.X
        # check for csr format, need to load small slice to get past dask wrapper
        matrix_slice = matrix[0:1, 0:1].compute()
        assert isinstance(
            matrix_slice, sparse.csr_matrix
        ), f"Matrix not in CSR Format, found {type(matrix_slice)}"

        def find_all_zero_indices_array_chunk(matrix_chunk) -> np.array:
            indices_array = matrix_chunk.indices
            indptr_array = matrix_chunk.indptr

            start, end = 0, matrix_chunk.shape[0]
            chunk_results = []
            while start < end:
                row_indices = indices_array[
                    indptr_array[start] : indptr_array[start + 1]
                ]
                # need size check to exclude visium spots with empty arrays, only want explicit 0 arrays
                if np.all(row_indices == 0) and row_indices.size > 0:
                    chunk_results.append(True)
                else:
                    chunk_results.append(False)
                start += 1

            return np.array([chunk_results]).reshape(-1, 1)

        def create_indices_class(barcode: str, matrix=matrix, adata=adata):
            index = adata.obs.index.get_loc(barcode)
            return IndicesResult(
                index=index,
                barcode=barcode,
                matrix_slice=matrix[index, :],
            )

        def create_final_results_list(barcodes: list[str]):
            final_result = [create_indices_class(barcode) for barcode in barcodes]
            # dask will most efficiently load slices with this compute instead of using loop
            matrix_slices = dask.compute(final_result)[0]
            data_arrays = [matrix.matrix_slice.data for matrix in matrix_slices]
            indices_arrays = [matrix.matrix_slice.indices for matrix in matrix_slices]
            for result_dc, data_array, indices_array in zip(
                final_result, data_arrays, indices_arrays
            ):
                result_dc.data_array = data_array
                result_dc.indices_array = indices_array

            return final_result

        def print_all_zero_results(rows: list[IndicesResult]):
            for result in rows:
                print(f"Row index: {result.index}")
                print(f"Row barcode: {result.barcode}")
                print(f"Row data array: {result.data_array}")
                print(f"Row indices array: {result.indices_array}")
                print("=" * 40)

        with dask.config.set(scheduler=worker_type):
            bool_mask = (
                map_blocks(find_all_zero_indices_array_chunk, matrix, dtype=int)
                .compute()
                .ravel()
            )

        barcodes = adata.obs[bool_mask].index

        # visium in_tissue == 0 can have empty indices array, only report in_tissue == 1
        # wrap in np.array to allow for size attribute check below
        if "in_tissue" in adata.obs.columns:
            barcodes = np.array(
                [
                    barcode
                    for barcode in barcodes
                    if adata.obs.loc[barcode].in_tissue == 1
                ]
            )

        if barcodes.size > 0:
            report(
                "ERROR: All-zero indices array found for the following cell(s):",
                level="ERROR",
            )
            final_results = create_final_results_list(barcodes)
            print_all_zero_results(final_results)
        else:
            report("Indices array per cell are not all-zero", level="GOOD")

    def evaluate_dup_counts(
        self, adata: ad.AnnData, worker_type="processes"
    ) -> pd.DataFrame | None:
        """
        Hash raw counts matrix in parallel with dask
        If there are duplicated rows, return a copy of obs with the duplicate rows' metadata

        This uses single-cell-curation read_h5ad() to lazily load matrices as dask arrays
        Dask will schedule out chunks to workers (processes in this case) and return a
        properly ordered array of hashes

        Key is to reshape the chunking function output to a columnar array with .reshape()
        and then use .ravel() after calling compute to flatten the final desired array.
        This avoids exceptions with broadcasting for the final chunk that != chunk_size

        Joyce did further work for the validator implementation to show that the sha224
        hash seems to be quickest

        Uses presence of obs "in_tissue" column to default to dense array hashing for
        visium datasets; this should prevent false positives for spots at the edge of the
        tissue borders. A further filtration step will remove obs rows with in_tissue==0

        Stupidly fast and efficient; on M1 Max with 64 GB, can hash the largest dataset
        (11.4 million cells) in about 3 1/2 minutes. The r5 16xlarge EC2 will take about a
        minute to 1:15, depending on chunk size. Memory usage is minimal when slicing through
        the data array.

        This may not work in a script, but storing here for the moment to be in the commit
        history

        Actually this does work when imported into a notebook
        """

        matrix = adata.raw.X if adata.raw else adata.X

        def hash_data_array_chunk(matrix_chunk) -> np.array:
            data_array = matrix_chunk.data
            indptr_array = matrix_chunk.indptr

            start, end = 0, matrix_chunk.shape[0]
            chunk_hashes = []
            while start < end:
                val = hashlib.sha224(
                    data_array[indptr_array[start] : indptr_array[start + 1]].tobytes()
                ).hexdigest()
                chunk_hashes.append(val)
                start += 1

            return np.array([chunk_hashes]).reshape(-1, 1)

        def hash_dense_chunk(matrix_chunk) -> np.array:
            chunk_hashes = [
                hashlib.sha224(r.tobytes()).hexdigest() for r in matrix_chunk.toarray()
            ]
            return np.array([chunk_hashes]).reshape(-1, 1)

        hash_chunk_func = (
            hash_dense_chunk
            if "in_tissue" in adata.obs.columns
            else hash_data_array_chunk
        )

        with dask.config.set(scheduler=worker_type):
            hashes = map_blocks(hash_chunk_func, matrix, dtype=int).compute().ravel()

        hash_df = adata.obs.copy()
        hash_df["row_hash"] = hashes

        if "in_tissue" in hash_df.columns:
            obs_to_keep = hash_df[hash_df["in_tissue"] != 0].index
            hash_df = hash_df[hash_df.index.isin(obs_to_keep)]

        dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

        if not dup_df.empty:
            report("duplicated raw counts", "ERROR")
            return dup_df
        report("no duplicated raw counts", "GOOD")


class NonDaskMatrixFunctions(BaseMatrixFunctions):
    def _determine_sparsity(self, x):
        if (
            isinstance(x, sparse.coo_matrix)
            or isinstance(x, sparse.csr_matrix)
            or isinstance(x, sparse.csc_matrix)
        ):
            sparsity = 1 - x.count_nonzero() / float(np.cumprod(x.shape)[-1])
        elif isinstance(x, np.ndarray):
            sparsity = 1 - np.count_nonzero(x) / float(np.cumprod(x.shape)[-1])
        else:
            report(
                f"matrix is of type {type(x)}, sparsity calculation has not been implemented"
            )

        return round(sparsity, 3)

    def evaluate_sparsity(self, adata: ad.AnnData):
        max_sparsity = 0.5

        valid = True
        sparsity = self._determine_sparsity(adata.X)
        report(f"X sparsity: {sparsity}")
        if sparsity > max_sparsity and type(adata.X) != sparse.csr_matrix:
            report("X should be converted to csr sparse", "ERROR")
            valid = False

        if adata.raw:
            sparsity = self._determine_sparsity(adata.raw.X)
            report(f"raw.X sparsity: {sparsity}")
            if sparsity > max_sparsity and type(adata.raw.X) != sparse.csr_matrix:
                report("raw.X should be converted to csr sparse", "ERROR")
                valid = False

        for l in adata.layers:
            sparsity = self._determine_sparsity(adata.layers[l])
            report(f"layers[{l}] sparsity: {sparsity}")
            if sparsity > max_sparsity and type(adata.layers[l]) != sparse.csr_matrix:
                report(f"layers[{l}] should be converted to csr sparse", "ERROR")
                valid = False

        if valid:
            report("all matrices have passed checks", "GOOD")

    def evaluate_data(self, adata: ad.AnnData):
        min_maxs = {}
        if adata.raw:
            raw_min = adata.raw.X.min()
            raw_max = adata.raw.X.max()
            report(f"raw min = {raw_min}")
            report(f"raw max = {raw_max}")
            min_maxs["raw"] = f"{raw_min}-{raw_max}"
            all_integers = np.all(np.round(adata.raw.X.data) == adata.raw.X.data)
        else:
            all_integers = np.all(np.round(adata.X.data) == adata.X.data)

        if all_integers:
            report("raw is all integers", "GOOD")
        else:
            report("raw contains non-integer values", "ERROR")

        X_min = adata.X.min()
        X_max = adata.X.max()
        report(f"X min = {X_min}")
        report(f"X max = {X_max}")
        min_maxs["X"] = f"{X_min}-{X_max}"

        for l in adata.layers:
            min = adata.layers[l].min()
            max = adata.layers[l].max()
            report(f"layers[{l}] min = {min}")
            report(f"layers[{l}] max = {max}")
            min_maxs[l] = f"{min}-{max}"

        poss_dups = [
            k for k, v in min_maxs.items() if list(min_maxs.values()).count(v) > 1
        ]
        if poss_dups:
            report(f"possible redundant layers: {poss_dups}", "WARNING")

    def evaluate_dup_counts(self, adata: ad.AnnData) -> pd.DataFrame | None:
        """
        Hash sparse csr matrix using np.ndarrays that represent sparse matrix data.
        First pass will hash all rows via slicing the data array and append to copy of obs df
        Second pass will hash only duplicate rows in obs copy via the indices array.
        This will keep only true duplicated matrix rows and not rows with an indicental same
        ordering of their data arrays
        """
        if "in_tissue" in adata.obs.columns:
            obs_to_keep = adata.obs[adata.obs["in_tissue"] != 0].index
            adata = adata[obs_to_keep, :]

        matrix = adata.raw.X if adata.raw else adata.X

        if not isinstance(matrix, sparse.csr_matrix):
            print("Matrix not in sparse csr format, please convert before hashing")
            return

        nnz = matrix.nnz

        if not matrix.has_canonical_format:
            print("Csr matrix not in canonical format, converting now...")
            if adata.raw:
                adata.raw.X.sort_indices()
                adata.raw.X.sum_duplicates()
            else:
                adata.X.sort_indices()
                adata.X.sum_duplicates()

        assert matrix.has_canonical_format, "Matrix still in non-canonical format"

        if nnz != matrix.nnz:
            print(f"{nnz - matrix.nnz} duplicates found during canonical conversion")

        data_array = matrix.data
        index_array = matrix.indices
        indptr_array = matrix.indptr

        start, end = 0, matrix.shape[0]
        hashes = []
        while start < end:
            val = hash(
                data_array[indptr_array[start] : indptr_array[start + 1]].tobytes()
            )
            hashes.append(val)
            start += 1

        def index_hash(index):
            obs_loc = adata.obs.index.get_loc(index)
            val = hash(
                index_array[indptr_array[obs_loc] : indptr_array[obs_loc + 1]].tobytes()
            )
            return val

        hash_df = adata.obs.copy()
        hash_df["data_array_hash"] = hashes
        hash_df = hash_df[
            hash_df.duplicated(subset="data_array_hash", keep=False) == True
        ]
        hash_df.sort_values("data_array_hash", inplace=True)

        hash_df["index_array_hash"] = [
            index_hash(row) for row in hash_df.index.to_list()
        ]
        hash_df = hash_df[
            hash_df.duplicated(
                subset=["data_array_hash", "index_array_hash"], keep=False
            )
            == True
        ]

        if not hash_df.empty:
            report("duplicated raw counts", "ERROR")
            return hash_df
        report("no duplicated raw counts", "GOOD")

    def evaluate_all_zero_indices(self, adata: ad.AnnData):
        """
        Function to check if a row/cell contains an all-zero indices array. This can exist
        in Visium datasets with spots that are in_tissue == 0, but should not exist for
        in_tissue == 1 or other single cell raw count data.
        Uses similar strategy as evaluate_dup_counts to slice through csr arrays
        """
        matrix = adata.raw.X if adata.raw else adata.X
        # check for csr format, check not needed for dense matrices
        assert isinstance(
            matrix, sparse.csr_matrix
        ), f"Matrix not in CSR Format, found {type(matrix)}"

        indices_array = matrix.indices
        indptr_array = matrix.indptr

        start, end = 0, matrix.shape[0]
        bool_mask = []
        while start < end:
            row_indices = indices_array[indptr_array[start] : indptr_array[start + 1]]
            # need size check to exclude visium spots with empty arrays, only want explicit 0 arrays
            if np.all(row_indices == 0) and row_indices.size > 0:
                bool_mask.append(True)
            else:
                bool_mask.append(False)
            start += 1

        def create_indices_class(barcode: str, matrix=matrix, adata=adata):
            index = adata.obs.index.get_loc(barcode)
            return IndicesResult(
                index=index,
                barcode=barcode,
                matrix_slice=matrix[index, :],
            )

        def create_final_results_list(barcodes: list[str]):
            result = [create_indices_class(barcode) for barcode in barcodes]
            data_arrays = [matrix.matrix_slice.data for matrix in result]
            indices_arrays = [matrix.matrix_slice.indices for matrix in result]
            for result_dc, data_array, indices_array in zip(
                result, data_arrays, indices_arrays
            ):
                result_dc.data_array = data_array
                result_dc.indices_array = indices_array

            return result

        def print_all_zero_results(rows: list[IndicesResult]):
            for result in rows:
                print(f"Row index: {result.index}")
                print(f"Row barcode: {result.barcode}")
                print(f"Row data array: {result.data_array}")
                print(f"Row indices array: {result.indices_array}")
                print("=" * 40)

        barcodes = adata.obs[bool_mask].index

        # visium in_tissue == 0 can have empty indices array, only report in_tissue == 1
        # wrap in np.array to allow for size attribute check below
        if "in_tissue" in adata.obs.columns:
            barcodes = np.array(
                [
                    barcode
                    for barcode in barcodes
                    if adata.obs.loc[barcode].in_tissue == 1
                ]
            )

        if barcodes.size > 0:
            report(
                "ERROR: All-zero indices array found for the following cell(s):",
                level="ERROR",
            )
            final_results = create_final_results_list(barcodes)
            print_all_zero_results(final_results)
        else:
            report("Indices array per cell are not all-zero", level="GOOD")


def get_adata_memory(
    adata_path: os.PathLike | str,
    print_datasets: bool = False,
    print_results: bool = False,
    sizes: Sizes | None = None,
) -> Sizes:
    """
    Calculate size of AnnData object when fully loaded in memory. Reads header/metadata information
    in h5/h5ad file and returns size of object loaded into RAM and calculated size on disk.
    In-memory sizes also loaded into attr_size_dict[dataset_name, nbytes] so further size calculations can be done
    on specific h5 datasets.

    :param: adata_path: str path to h5/h5ad file. Will only load header/metadata info, not full file
    :param: print_datasets: Default False, set to True to get print out of individual datasets
    :param: print_results: Default False, set to True to get print out of final results

    :returns: Sizes object
    """
    if sizes is None:
        sizes = Sizes()

    print_width = 90

    def dataset_sizes(name, obj):
        if isinstance(obj, h5py.Dataset):
            stor_size = obj.id.get_storage_size()
            if print_datasets:
                header = f"RAM Size for {obj.name}:"
                spaces = (
                    print_width - len(header) - len(f"{obj.nbytes:_}") - 6
                )  # 6 for ' bytes'
                print(f"{header}{' ' * spaces}{obj.nbytes:_} bytes")
            sizes.memory_size += obj.nbytes
            sizes.disk_size += stor_size
            sizes.attr_size_dict[obj.name] = obj.nbytes

    with h5py.File(adata_path, mode="r") as f:
        f.visititems(dataset_sizes)

    if print_results:
        print("-" * print_width)
        for header, result in {
            "Size in RAM:": sizes.memory_size,
            "Size on disk:": sizes.disk_size,
        }.items():
            spaces = print_width - len(f"{header}") - len(f"{result:_}") - 6
            print(f"{header}{' ' * spaces}{result:_} bytes")

    return sizes


def get_read_h5ad_function(file: os.PathLike | str) -> Callable:
    """
    Determine whether to load matrices lazily or fully depending on local memory
    Will use single-cell-curation's read_h5ad() function if adata will take up more
    memory than is locally available.
    """
    total_local_memory = psutil.virtual_memory().total
    total_adata_memory = get_adata_memory(file).memory_size
    if total_adata_memory > total_local_memory:
        print(
            f"Adata will take {total_adata_memory:_} bytes, "
            f"greater than local {total_local_memory:_} bytes"
        )
        print("Lazily loading adata with dask arrays...")
        return read_h5ad

    print(
        f"Adata will take {total_adata_memory:_} bytes, "
        f"less than local {total_local_memory:_} bytes"
    )
    print("Loading adata with ad.read_h5ad()...")
    return ad.read_h5ad
