import anndata as ad
import dask
import hashlib
import h5py
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import scanpy as sc
import squidpy as sq
import subprocess
import sys
from anndata.compat import DaskArray
from dataclasses import dataclass
from dask.array import map_blocks
from pathlib import Path
from scipy import sparse
from cellxgene_ontology_guide.ontology_parser import OntologyParser
import cellxgene_schema.gencode as gencode
import cellxgene_schema.utils as utils
import cellxgene_schema.schema as schema
from typing import Union
from cellxgene_schema.utils import SPARSE_MATRIX_TYPES, read_h5ad


portal_uns_fields = [
    'citation',
    'schema_reference',
    'schema_version',
    'organism'
]

curator_uns_fields = [
    'title',
    'organism_ontology_term_id'
]

portal_var_fields = [
    'feature_name',
    'feature_reference',
    'feature_biotype',
    'feature_type',
    'feature_length'
]

portal_obs_fields = [
    'assay',
    'cell_type',
    'development_stage',
    'disease',
    'self_reported_ethnicity',
    'sex',
    'tissue'
]
non_ontology_fields = ['donor_id','suspension_type','tissue_type','is_primary_data']
curator_obs_fields = [e + '_ontology_term_id' for e in portal_obs_fields] + non_ontology_fields
full_obs_standards = portal_obs_fields + curator_obs_fields
ONTOLOGY_PARSER = OntologyParser(schema_version=schema.get_current_schema_version())

def get_path(search_term: str) -> os.PathLike | str:
    """
    Find path of local repos and API keys regardless of source machine. Use Path objects and
    likely locations instead of glob or rglob to limit search overhead for simple import

    Returns Path when found, otherwise str "Path not found"
    """
    # should start at ./lattice-tools/cellxgene_resources/
    local_path = Path()
    
    likely_locations = [
        local_path.resolve().parent.parent,                 # same level as lattice-tools
        local_path.home(),
        local_path.home() / "CZI",
        local_path.home() / "GitClones",
        local_path.home() / "GitClones" / "CZI",
        local_path.home() / "GitClones" / "Lattice-Data",   # if other local lattice repos beyond lattice-tools
        local_path.home() / "Documents" / "keys",
        local_path.home() / "keys",
        local_path.home() / "Documents",
        local_path.home() / "Desktop" / "Curation",
    ]

    for place in likely_locations:
        if place.exists():
            for item in place.iterdir():
                if search_term in item.name:
                    return item

    return "Path not found"


class CxG_API:
    scc_repo_loc = get_path("single-cell-curation")

    if isinstance(scc_repo_loc, Path):
        api_source = scc_repo_loc.resolve() / "notebooks" / "curation_api" / "python"
        sys.path.append(str(api_source))
    else:
        print("Path not found for single-cell-curation repo")
                

    from src.collection import create_collection,create_revision,get_collection,get_collections,update_collection
    from src.dataset import create_dataset,delete_dataset,get_dataset,get_datasets,upload_datafile_from_link,upload_local_datafile,upload_datafiles_from_manifest,get_dataset_manifest

    def config(env="prod"):
        from src.utils.config import set_api_access_config

        api_key_files = {
            "prod": "cxg-api-key.txt",
            "dev": "cxg-api-key-dev.txt",
            "staging": "cxg-api-key-staging.txt",
        }

        api_key_file_path = get_path(api_key_files[env])
        set_api_access_config(api_key_file_path, env=env)


def report(mess, level=None):
    colors = {
        'GOOD': '\033[32m', #green
        'WARNING': '\033[33m', #yellow
        'ERROR': '\033[31m' #red
    }
    if level:
        c = colors[level]
        print(f'\033[1m{c}{level}: {mess}\033[0m')
    else:
        print(mess)


def revise_cxg(adata):
    for p in portal_uns_fields:
        del adata.uns[p]

    adata.obs.drop(columns=portal_obs_fields, inplace=True)
    adata.obs.drop(columns='observation_joinid', inplace=True)
    adata.var.drop(columns=portal_var_fields, inplace=True)

    if adata.raw:
        adata.raw.var.drop(columns=portal_var_fields, inplace=True)

    return adata


@dataclass
class Sizes:
    memory_size: int = 0
    disk_size: int = 0
    attr_size_dict = {}


def calculate_adata_memory(adata_path: str, print_datasets: bool = False, sizes: Sizes | None = None) -> Sizes:
    """
    Calculate size of AnnData object when fully loaded in memory. Reads header/metadata information
    in h5/h5ad file and returns size of object loaded into RAM and calculated size on disk. 
    In-memory sizes also loaded into attr_size_dict[dataset_name, nbytes] so further size calculations can be done
    on specific h5 datasets.
    
    :param: adata_path: str path to h5/h5ad file. Will only load header/metadata info, not full file
    :param: print_datasets: Default False, set to True to get print out of individual datasets

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
                spaces = print_width - len(header) - len(f"{obj.nbytes:_}") - 6  # 6 for ' bytes'
                print(f"{header}{' ' * spaces}{obj.nbytes:_} bytes")
            sizes.memory_size += obj.nbytes
            sizes.disk_size += stor_size
            sizes.attr_size_dict[obj.name] = obj.nbytes

    with h5py.File(adata_path, mode='r') as f:
        f.visititems(dataset_sizes)

    print("-" * print_width)
    for header, result in {
        "Size in RAM:": sizes.memory_size, 
        "Size on disk:": sizes.disk_size
    }.items():
        spaces = print_width - len(f"{header}") - len(f"{result:_}") - 6
        print(f"{header}{' ' * spaces}{result:_} bytes")

    return sizes


def get_matrix_format(matrix: DaskArray) -> str:
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


def count_matrix_nonzero(matrix: DaskArray) -> int:
    def count_nonzeros(matrix_chunk: Union[np.ndarray, sparse.spmatrix], is_sparse_matrix: bool) -> np.array:
        nnz = matrix_chunk.nnz if is_sparse_matrix else np.count_nonzero(matrix_chunk)
        return np.array([nnz])

    is_sparse_matrix = get_matrix_format(matrix) in SPARSE_MATRIX_TYPES
    # if matrix too small to chunk, then just call function and don't dask distribute
    if len(matrix.chunks[0]) > 1:
        nonzeros = map_blocks(
            count_nonzeros,
            matrix, 
            is_sparse_matrix, 
            drop_axis=1, 
            dtype=int) \
                .compute(scheduler='processes') \
                .sum()
    else:
        nonzeros = count_nonzeros(matrix.compute(), is_sparse_matrix)[0]
    return nonzeros


def determine_sparsity(matrix: DaskArray):
    nnz = count_matrix_nonzero(matrix)
    sparsity = 1 - nnz / np.prod(matrix.shape)
    return round(sparsity, 3)


def get_matrices_to_evaluate(adata: ad.AnnData) -> list[tuple[Union[np.ndarray, sparse.spmatrix], str]]:
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


def evaluate_sparsity(adata_path: os.PathLike | str):
    """
    Evaluate sparsity in qa notebook

    Expects matrices as dask array, so make sure to load Anndata with
    read_h5ad()
    """
    adata = read_h5ad(adata_path)
    max_sparsity = 0.5
    valid = True

    matrices_to_evaluate = get_matrices_to_evaluate(adata)

    for matrix, matrix_name in matrices_to_evaluate:
        format = get_matrix_format(matrix)
        sparsity = determine_sparsity(matrix)
        report(f"{matrix_name} sparsity: {sparsity}")
        if sparsity > max_sparsity and format not in SPARSE_MATRIX_TYPES:
            report(f'{matrix_name} should be converted to csr sparse, found to be {format}', 'ERROR')
            valid = False

    if valid:
        report("All matrices have passed checks", "GOOD")


def evaluate_data(adata_path: os.PathLike | str):
    """
    3 other data checks on matrix: min, max, and if raw matrix, all integer check

    For better efficency, all-integer check only done on raw matrix
    Probably room for improvement to delay chunks, better call compute(), or
    other dask optimizations to not chunk through each matrix 2-3 times
    """
    adata = read_h5ad(adata_path)
    min_maxs = {}
    matrices_to_evaluate = get_matrices_to_evaluate(adata)
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
        with dask.config.set(scheduler='processes'):
            max = map_blocks(get_max_chunk, matrix, dtype=int).compute().ravel().max()
            min = map_blocks(get_min_chunk, matrix, dtype=int).compute().ravel().min()
            if matrix_name == raw_matrix:
                is_all_integers = map_blocks(all_integers_chunk, matrix, dtype=bool).compute().ravel().all()

        min_maxs[matrix_name] = f"{min}-{max}"

        report(f"{matrix_name} min = {min}")
        report(f"{matrix_name} max = {max}")

        if matrix_name == raw_matrix:
            if is_all_integers:
                report(f"{matrix_name} is all integers", "GOOD")
            else:
                report(f"{matrix_name} is NOT all integers", "ERROR")
            

    poss_dups = [k for k,v in min_maxs.items() if list(min_maxs.values()).count(v) > 1]

    if poss_dups:
        report(f'possible redundant layers: {poss_dups}','WARNING')


@dataclass
class IndicesResult:
    index: int
    barcode: str
    matrix_slice: sparse.csr_matrix
    data_array = None
    indices_array = None


def evaluate_all_zero_indices(adata_path: os.PathLike | str, worker_type="processes"):
    """
    Function to check if a row/cell contains an all-zero indices array. This can exist
    in Visium datasets with spots that are in_tissue == 0, but should not exist for
    in_tissue == 1 or other single cell raw count data.
    Uses dask to lazily load the raw matrix and to check in parallel
    """
    adata = read_h5ad(adata_path)
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
            row_indices = indices_array[indptr_array[start] : indptr_array[start + 1]]
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
            matrix_slice=matrix[index, : ],
        )

    def create_final_results_list(barcodes: list[str]):
        final_result = [create_indices_class(barcode) for barcode in barcodes]
        # dask will most efficiently load slices with this compute instead of using loop
        matrix_slices = dask.compute(final_result)[0]
        data_arrays = [matrix.matrix_slice.data for matrix in matrix_slices]
        indices_arrays = [matrix.matrix_slice.indices for matrix in matrix_slices]
        for result_dc, data_array, indices_array in zip(final_result, data_arrays, indices_arrays):
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
            [barcode for barcode in barcodes if adata.obs.loc[barcode].in_tissue == 1]
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


def evaluate_uns_colors(adata):
    numb_types = ['int_', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64','float_', 'float16', 'float32', 'float64']

    colors_keys = [k for k in adata.uns.keys() if k.endswith('_colors')]
    if colors_keys:
        for k in colors_keys:
            colors = len(adata.uns[k])
            obs_field = k[:-(len('_colors'))]

            if obs_field in portal_obs_fields:
                report(f'uns.{k} not allowed, move to uns.{obs_field}_ontology_term_id_colors', 'ERROR')
            elif obs_field not in adata.obs.keys():
                report(f'{obs_field} not found in obs, consider DELETING or RENAMING uns.{k}', 'ERROR')
            else:
                valid = True
                values = len(adata.obs[obs_field].unique())
                if colors < values:
                    report(f'uns.{k} has only {str(colors)} colors but obs.{obs_field} has {str(values)} values', 'ERROR')
                    valid = False
                if adata.obs.dtypes[obs_field].name in numb_types:
                    report(f'uns.{k} is associated with non-categorical {obs_field}', 'ERROR')
                    valid = False
                if valid:
                    report(f'uns.{k} defined appropriately', 'GOOD')
    else:
        report('no _colors keys defined')


def map_filter_gene_ids(adata):
    #map genes
    gene_map = {}
    gene_map_files = ['gene_map_human_v48.json','gene_map_mouse_v37.json','gene_map_fly_v114.json']
    for file in gene_map_files:
        with open(f'../gene_ID_mapping/{file}', 'r') as f:
            data = json.load(f)
            gene_map.update(data)

    approved_file = 'ref_files/genes_approved.csv.gz'
    approved = pd.read_csv(approved_file,dtype='str')

    my_gene_map = {k:v for k,v in gene_map.items() if k in adata.var.index and v not in adata.var.index}
    adata.var.rename(index=my_gene_map, inplace=True)

    #filter out genes
    if not adata.var.index.name:
        adata.var.index.name = 'ensembl_id'
    index_name = adata.var.index.name
    adata.var.reset_index(inplace=True)
    var_to_keep = adata.var[adata.var[index_name].isin(approved['feature_id'])].index #what if it's not called 'gene_ids'
    adata = adata[:, var_to_keep]
    adata.var.set_index(index_name, inplace=True)

    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs) #do we need to define obs?
        my_gene_map = {k:v for k,v in gene_map.items() if k in raw_adata.var.index and v not in raw_adata.var.index}
        raw_adata.var.rename(index=my_gene_map, inplace=True)
        if not raw_adata.var.index.name:
            raw_adata.var.index.name = 'ensembl_id'
        index_name = raw_adata.var.index.name
        raw_adata.var.reset_index(inplace=True)
        var_to_keep = raw_adata.var[raw_adata.var[index_name].isin(approved['feature_id'])].index
        raw_adata = raw_adata[:, var_to_keep]
        raw_adata.var.set_index(index_name, inplace=True)
        
        adata.raw = raw_adata

    return adata


def evaluate_10x_barcodes(obs, visium=False):
    vis_terms = ['EFO:0010961','EFO:0022857','EFO:0022858','EFO:0022859','EFO:0022860']
    if 'assay_ontology_term_id' in obs.columns and [e for e in obs['assay_ontology_term_id'].unique() if e in vis_terms]:
        visium=True

    if visium:
        csv = 'ref_files/visium_barcode_table.csv.gz'
    else:
        csv = 'ref_files/10X_barcode_table.csv.gz'
    ref_df = pd.read_csv(csv, sep=',', header=0, index_col='barcode')
    global barcode_headers
    barcode_headers = ref_df['summary'].unique()

    global no_barcode_v
    no_barcode_v = 'no barcode'

    obs = obs.copy()
    obs['barcode'] = obs.index.str.extract(r'([ACTG]{12,})')[0].tolist()
    if len(set(ref_df.index.to_list()).intersection(set(obs['barcode'].to_list()))) == 0:
        report('Did not find any barcodes in obs index, cannot evaluate barcodes', 'WARNING')
        return
    obs = obs.merge(ref_df[['summary']],on='barcode',how='left').set_index(obs.index)
    obs['summary'] = obs.apply(
        lambda x: no_barcode_v if pd.isna(x['barcode']) else (f"{len(x['barcode'])}nt" if pd.isna(x['summary']) else x['summary']),
        axis=1
    )

    return obs


def parse_barcode_df(df, field):
    if df is None:
        return
    results = {}

    for a in df[field].unique():
        temp = df[df[field] == a]
        results[a] = temp['summary'].value_counts().to_dict()

    df = pd.DataFrame(results).fillna(0).astype(int).transpose()
    for h in list(barcode_headers) + [no_barcode_v]:
        if h not in df.columns:
            df[h] = 0
    df = df[[c for c in df if df[c].sum() > 0 and c not in ['multiple',no_barcode_v] and not c.endswith('nt')]
            + [c for c in df if df[c].sum() > 0 and c.endswith('nt')]
            + [c for c in df if df[c].sum() == 0 and c not in ['multiple',no_barcode_v]]
            + [c for c in ['multiple',no_barcode_v] if c in df]]
    df.sort_values(list(df.columns), ascending=False, inplace=True)

    return df


def evaluate_obs_schema(obs, labels=False):
    if labels:
        for o in portal_obs_fields + non_ontology_fields:
            if o not in obs.keys():
                report(f'{o} not in obs\n', 'ERROR')
            else:
                report(f'{o} {obs[o].unique().tolist()}\n')
    else:
        for o in curator_obs_fields:
            if o not in obs.keys():
                report(f'{o} not in obs\n', 'ERROR')
            else:
                report(f'{o} {obs[o].unique().tolist()}\n')
        for o in portal_obs_fields:
            if o in obs.keys():
                report(f'schema conflict - {o} in obs\n', 'ERROR')
    if 'cell_type_ontology_term_id' in obs.columns and 'unknown' in obs['cell_type_ontology_term_id'].unique():
        if 'in_tissue' in obs.columns:
            num_unknown = obs.loc[(obs['in_tissue']==1) & (obs['cell_type_ontology_term_id']=='unknown')].shape[0]
            perc_unknown = 100*(num_unknown/obs.loc[obs['in_tissue']==1].shape[0])
        else:
            num_unknown = obs[obs['cell_type_ontology_term_id']=='unknown'].shape[0]
            perc_unknown = 100*(num_unknown/obs.shape[0])
        if num_unknown > 20:
            report(f'{num_unknown} ({perc_unknown}%) cells are cell_type:unknown.', 'WARNING')



def evaluate_obs(obs, full_obs_standards):
    long_fields = []
    gradient_fields = []
    uber_dict = {}
    for o in obs.keys():
        vc_dict = obs[o].value_counts(dropna=False).to_dict()
        counts = '_'.join([str(c) for c in vc_dict.values()])
        count_len = len(vc_dict.keys())
        values = [str(i) for i in vc_dict.keys()]
    
        if o.startswith(' ') or o.endswith(' ') or '  ' in o:
            report(f'leading/trailing whitespace: {o}\n')
    
        if o not in full_obs_standards and ' '.join(o.split()).lower() in full_obs_standards:
            report(f'schema conflict: {o}\n')

        numb_types = ['int_', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16',
                      'uint32', 'uint64','float_', 'float16', 'float32', 'float64']
        if obs.dtypes[o].name in numb_types:
            gradient_fields.append(o)
        else:
            #check for long categories as they will not be enabled for coloring
            if count_len > 200 and o != 'observation_joinid':
                long_fields.append(o)
    
            #report value_counts to later look for redundancy
            metadata = {
                'values': values,
                'property': o
            }
            if counts in uber_dict:
                uber_dict[counts].append(metadata)
            else:
                uber_dict[counts] = [metadata]
    for k,v in uber_dict.items():
        if '_' in k and not k.startswith('1_1'):
            props = [e['property'] for e in v]
            if len(v) > 1 and not all(elem in full_obs_standards for elem in props):
                report(f'possible redundancy: {[e["property"] for e in v]}\n')

    if gradient_fields:
        report(f'continuous fields: {gradient_fields}\n')
    if long_fields:
        report(f'long fields: {long_fields}')


def evaluate_dup_counts(adata_path: os.PathLike | str, worker_type="processes") -> pd.DataFrame | None:
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

    adata = read_h5ad(adata_path)
    matrix = adata.raw.X if adata.raw else adata.X

    def hash_data_array_chunk(matrix_chunk) -> np.array:
        data_array = matrix_chunk.data
        indptr_array = matrix_chunk.indptr

        start, end = 0, matrix_chunk.shape[0]
        chunk_hashes = []
        while start < end:
            val = hashlib.sha224(data_array[indptr_array[start]:indptr_array[start + 1]].tobytes()).hexdigest()
            chunk_hashes.append(val)
            start += 1
            
        return np.array([chunk_hashes]).reshape(-1, 1)

        
    def hash_dense_chunk(matrix_chunk) -> np.array:
        chunk_hashes = [hashlib.sha224(r.tobytes()).hexdigest() for r in matrix_chunk.toarray()]
        return np.array([chunk_hashes]).reshape(-1, 1)


    hash_chunk_func = hash_dense_chunk if "in_tissue" in adata.obs.columns else hash_data_array_chunk

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

    
def symbols_to_ids(symbols, var):
    """
    Given a list of gene symbols, look in genes_approved.csv.gz to see if we can map to an Ensembl ID that is found
    in adata.var. If there not a successful mapping of gene symbols to Ensembl ID, will take the lower case version of 
    gene symbol and try mapping again.

    :param symbols: List of upper cased gene symbols that would like to find Ensembl ID mapping for
    :param var: adata.var

    :return ensg_list: List of Ensembl IDs found

    """
    ref_dir = 'ref_files/'
    if not os.path.exists(ref_dir + 'genes_approved.csv.gz'):
        report('There is no genes_approved.csv.gz file present', 'ERROR')
        return

    approved = pd.read_csv(ref_dir + 'genes_approved.csv.gz',dtype='str')
    approved['symbol_only'] = approved['symb'].str.split('_', expand=True)[0]
    
    ensg_list = []
    for s in symbols:
        found = False
        if s in approved['symbol_only'].tolist():
            ensg_ids = approved.loc[approved['symbol_only'] == s, 'feature_id']
            for ensg_id in ensg_ids:
                if ensg_id in var.index:
                    ensg_list.append(ensg_id)
                    report(f'{ensg_id} -- {s}')
                    found = True
        if not found:
            s_lower = s[0] + s[1:].lower()
            if s_lower in approved['symbol_only'].tolist():
                ensg_ids_lower = approved.loc[approved['symbol_only'] == s_lower, 'feature_id']
                for ensg_id_lower in ensg_ids_lower:
                    if ensg_id_lower in var.index:
                        ensg_list.append(ensg_id_lower)
                        report(f'{ensg_id_lower} -- {s_lower}')
                        found = True
        if not found:
            report(f'{s} not found in genes_approved or adata.var')

    return ensg_list


def pick_embed(keys):
    for k in keys:
        if 'umap' in k.lower():
            return k
        elif 'umap' in k.lower():
            return k

    return keys[0]


def plot_vis(adata, cellpop_field):
    ncols = 2
    nrows = 1
    figsize = 4
    wspace = 0.5
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
    )
    plt.subplots_adjust(wspace=wspace)
    lib = [k for k in adata.uns['spatial'].keys() if k != 'is_single'][0]

    sq.pl.spatial_scatter(adata, ax=axs[0], library_id=lib, color=cellpop_field, return_ax=False)
    sq.pl.spatial_scatter(adata, ax=axs[1], library_id=lib)

    if 'fullres' in adata.uns['spatial'][lib]['images']:
        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
        )
        plt.subplots_adjust(wspace=wspace)
        sq.pl.spatial_scatter(adata, ax=axs[0], library_id=lib, img_res_key='fullres', scale_factor=1.0, color=cellpop_field, return_ax=False)
        sq.pl.spatial_scatter(adata, ax=axs[1], library_id=lib, img_res_key='fullres', scale_factor=1.0)
    else:
        report('fullres image is highly recommended', 'WARNING')


def validate(file):
    validate_process = subprocess.run(['cellxgene-schema', 'validate', file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in validate_process.stdout.decode('utf-8').split('\n'):
        report(line)
    for line in validate_process.stderr.decode('utf-8').split('\n'):
        if line.endswith('is_valid=True'):
            report(line, 'GOOD')
            return True
        elif line.endswith('is_valid=False'):
            report(line, 'ERROR')
            return False
        else:
            prefix = line.split(':')[0]
            if prefix in ['ERROR','WARNING']:
                report(line.replace(f'{prefix}:',''), prefix)
            else:
                report(line)


def compare_revision(collection):
    change = False
    if collection.get('revising_in'):
        revision_id = collection['revising_in']
        revision = CxG_API.get_collection(revision_id)
    elif collection.get('revision_of'):
        revision = collection
        collection_id = collection['revision_of']
        collection = CxG_API.get_collection(collection_id)

    should_differ_collection = [
        'collection_id', 'collection_url', 'collection_version_id',
        'created_at', 'revising_in', 'revision_of', 'visibility'
    ]
    should_be_absent = [
        'processing_status'
    ]
    should_differ_dataset = [
        'dataset_version_id','explorer_url','assets','revised_at','citation','processing_status'
    ]
    ont_fields = [
        'assay','cell_type','development_stage','disease',
        'self_reported_ethnicity','sex','tissue','organism'
    ]
    for k,v in revision.items():
        if k not in collection.keys():
            if k not in should_be_absent:
                print('not present: ' + k)
                change = True
        elif collection.get(k) != v and k not in should_differ_collection:
            if k == 'datasets':
                diff_props = set()
                pub_datasets = {d['dataset_id']: d for d in collection[k]}
                rev_datasets = {d['dataset_id']: d for d in revision[k]}
                comp = {}
                new = {}
                removed = {}
                for ds_id,v in rev_datasets.items():
                    if ds_id not in pub_datasets.keys():
                        new[ds_id] = {}
                        for p in ['title','cell_count']:
                            new[ds_id][p] = v[p]
                        for p in ['assay','organism','tissue']:
                            new[ds_id][p] = [a['label'] for a in v[p]]
                        change = True
                    else:
                        comp[ds_id] = {'title': v['title']}
                        for prop,rev_val in v.items():
                            if prop not in should_differ_dataset:
                                pub_val = pub_datasets[ds_id].get(prop)
                                if prop in ont_fields:
                                    rev_val = [t['label'] for t in rev_val]
                                    pub_val = [t['label'] for t in pub_val]
                                if isinstance(rev_val, list) and prop != 'assets':
                                    rev_val.sort()
                                    pub_val.sort()
                                if pub_val != rev_val:
                                    if prop == 'mean_genes_per_cell' and round(rev_val, 5) == round(pub_val, 5):
                                        continue
                                    diff_props.add(prop)
                                    change = True
                                    comp[ds_id][prop + '_REV'] = rev_val
                                    comp[ds_id][prop + '_PUB'] = pub_val
                for ds_id,v in pub_datasets.items():
                    if ds_id not in rev_datasets.keys():
                        removed[ds_id] = {}
                        for p in ['title','cell_count']:
                            removed[ds_id][p] = v[p]
                        for p in ['assay','organism','tissue']:
                            removed[ds_id][p] = [a['label'] for a in v[p]]
                        change = True
            else:
                print('not same: ' + k)
                if k not in ['datasets','publisher_metadata']:
                    if k in ['links']:
                        diff_in_pub = [l for l in collection[k] if l not in v]
                        print('--- published: ', diff_in_pub)
                        diff_in_rev = [l for l in v if l not in collection[k]]
                        print('----- revised: ', diff_in_rev)
                        change = True
                    else:
                        print('--- published: ', str(collection[k]))
                        print('----- revised: ', v)
                        change = True


    comp_df = pd.DataFrame(comp).transpose()
    comp_df = comp_df.dropna(subset=[c for c in comp_df.columns if c != 'title'], how='all')
    if not comp_df.empty:
        print('\033[1mRevised Datasets\033[0m')
        change = True

        cols = list(comp_df)
        cols.insert(0, cols.pop(cols.index('title')))
        comp_df = comp_df.loc[:, cols]

        a = ['title'] + [c[-3:] for c in comp_df.columns if c not in ['title']]
        b = [''] + [c[:-4] for c in comp_df.columns if c not in ['title']]
        comp_df.columns = [b, a]
        display(comp_df.fillna(''))

        print('\033[1mProperty Comparison\033[0m')
        for f in diff_props:
            if f == 'title':
                continue
            temp = comp_df[(comp_df[f]['REV'] != comp_df[f]['PUB']) & (comp_df[f]['PUB'].isna() == False)]
            for i,row in temp.iterrows():
                p = row[f]['PUB']
                if isinstance(p, (int, float)):
                    continue
                r = row[f]['REV']
                only_in_pub = [str(e) for e in p if e not in r]
                only_in_rev = [str(e) for e in r if e not in p]
                print(i + '-' + f)
                if only_in_pub:
                    print('only in pub:' + ','.join(only_in_pub))
                if only_in_rev:
                    print('only in rev:' + ','.join(only_in_rev))
                print('---------')

    if new:
        print('\033[1mNew Datasets\033[0m')
        change = True
        display(pd.DataFrame(new).transpose())

    if removed:
        print('\033[1mRemoved Datasets\033[0m')
        change = True
        display(pd.DataFrame(removed).transpose())

    if not change:
        report('no changes changes detectable based on API response')

    return revision


def generate_fm_dict(female_ids, female_adata, male_ids, male_adata, adata):
    """
    Input: both sets of gene ids, female and male subset anndatas, and original anndata object
    Output: dictionary containing sex-specific gene_ids, subset anndatas and dataframes of summed raw expression counts per gene per donor
    """
    fm_dict = {'female': [female_ids, female_adata], 'male': [male_ids, male_adata]}
    for k,v in fm_dict.items():
        sex_specific_adata = v[1]
        obs_data = pd.DataFrame(sex_specific_adata[:, sex_specific_adata.var.index].X.toarray(), columns=sex_specific_adata.var.index, index=adata.obs.index)
        obs_donor_data = pd.merge(obs_data,adata.obs['donor_id'], how='left',left_index=True, right_index=True)
        df = obs_donor_data.groupby(['donor_id'])[obs_donor_data.columns].sum(numeric_only=True).reset_index()
        missing_genes = [g for g in v[0] if g not in df.columns]
        df[missing_genes] = np.nan
        v.append(df)

    return fm_dict


def assign_sex(x):
    """
    Input: ratio of male to female raw expression counts summed across all genes
    Output: assignment of sex
    """
    if x > 0.35:
        return 'male'
    elif x < 0.05:
        return 'female'
    else:
        return 'unknown'


def check_percent(female_adata,male_adata,female_ids,male_ids):
    """
    Input: subset adatas and gene_id lists
    Output: percent of genes found per sex
    """
    fem = female_adata.shape[1]/len(female_ids)
    male = male_adata.shape[1]/len(male_ids)
    print(f"% Female genes found: {(fem)*100}")
    print(f"% Male genes found: {(male)*100}\n")

    return (fem,male)


def calculate_sex(fm_dict):
    """
    Input: output dictionary from generate_fm_dict()
    Output: dataframe of both male and female donors that have 100+ total counts and stats
    """
    try:
        female_df = fm_dict['female'][2]
        female_df['female_sum'] = female_df.sum(numeric_only=True, axis=1)
        male_df = fm_dict['male'][2]
        male_df['male_sum'] = male_df.sum(numeric_only=True, axis=1)
        male_female_df = pd.merge(male_df, female_df, how='left', on='donor_id')
        male_female_df['total_sum'] = male_female_df[['female_sum','male_sum']].sum(numeric_only=True, axis=1)
        donors_to_remove = male_female_df[male_female_df.total_sum < 100].donor_id.unique()  # Remove donors that have less than 100 total counts
        if len(donors_to_remove) > 0:
            print('Donors with < 100 total counts dropped: ', *(donors_to_remove), sep='\n')
            male_female_df.drop(male_female_df[male_female_df.total_sum < 100].index, inplace=True)
        #Calculate ratio and assign sex
        male_female_df['male_to_female'] = male_female_df['male_sum']/male_female_df['female_sum']
        male_female_df['scRNAseq_sex'] = male_female_df.apply(lambda x: assign_sex(x['male_to_female']), axis=1)

        return male_female_df,donors_to_remove

    except Exception as e:
        print(e)


def evaluate_donors_sex(adata):
    if 'NCBITaxon:9606' != adata.uns['organism_ontology_term_id']:
        print('Cannot calculate sex for non-human data.')
        return None,None
    else:
        genes_file = 'ref_files/sex_analysis_genes.json'
        genes = json.load(open(genes_file))
        female_ids = genes['female'].keys()
        male_ids = genes['male'].keys()
        metadata_list = ['donor_id', 'sex_ontology_term_id','assay_ontology_term_id']
        smart_assay_list = ['EFO:0010184','EFO:0008931','EFO:0008930','EFO:0010022','EFO:0700016','EFO:0022488','EFO:0008442']
        adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
        adata.obs.loc[adata.obs['assay_ontology_term_id'].isin(smart_assay_list) == True, 'donor_id'] += '-smartseq'

        if adata.raw:
            adata = ad.AnnData(sparse.csr_matrix(adata.raw.X), var=adata.raw.var, obs=adata.obs)

        female_adata = adata[:,adata.var.index.isin(female_ids)]
        male_adata = adata[:,adata.var.index.isin(male_ids)]

        genes_found = check_percent(female_adata,male_adata,female_ids,male_ids)
        if genes_found[0] == 0 or genes_found[1] == 0:
            return None,None
        fm_counts_dict = generate_fm_dict(female_ids,female_adata,male_ids,male_adata,adata)
        donor_sex_df,removed_donors = calculate_sex(fm_counts_dict)
        donor_sex_df = donor_sex_df[['donor_id','male_to_female','scRNAseq_sex']]
        donor_sex_df = donor_sex_df.merge(adata.obs[metadata_list].drop_duplicates(), on='donor_id', how='left')
        sex_map = {
            'PATO:0000383':'female',
            'PATO:0000384':'male',
            'unknown':'unknown'
        }
        donor_sex_df['author_annotated_sex'] = donor_sex_df['sex_ontology_term_id'].map(sex_map)
        donor_sex_df.loc[donor_sex_df['assay_ontology_term_id'].isin(smart_assay_list) == True, 'smart_seq'] = True
        donor_sex_df.drop(columns=['sex_ontology_term_id','assay_ontology_term_id'], inplace=True)
        donor_sex_df.sort_values('male_to_female', inplace=True)
        obs_to_keep = []
        ratio_order = []
        smart_seq_donors_rename = {}

        if 'smart_seq' in donor_sex_df.columns:
            donor_sex_df['smart_seq'] = donor_sex_df['smart_seq'].fillna(False).astype('bool')

        if donor_sex_df['smart_seq'].all() or not any(donor_sex_df['smart_seq']):
            adata.obs['donor_id'] = adata.obs['donor_id'].str.split('-smartseq').str[0]
            donor_sex_df['donor_id'] = donor_sex_df['donor_id'].str.split('-smartseq').str[0]
            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin((donor_sex_df[donor_sex_df['donor_id'].isin(removed_donors)!=True]['donor_id']))].index)
            ratio_order.append((donor_sex_df['donor_id'] + ' ' + donor_sex_df['author_annotated_sex'].astype('string')).to_list())

        else:
            for d in pd.Series(donor_sex_df['donor_id'].str.split('-smartseq').str[0]).unique():
                if d in removed_donors:
                    print(f"Donor {d} was removed from analysis, cannot include in plot.")
                else:
                    try:
                        smart_seq_sex = donor_sex_df.loc[(donor_sex_df['donor_id'] == d + '-smartseq') & (donor_sex_df['smart_seq'] == True)]['scRNAseq_sex'].unique()
                        nonsmart_seq_sex = donor_sex_df.loc[(donor_sex_df['donor_id'] == d) & (donor_sex_df['smart_seq'] == False)]['scRNAseq_sex'].unique()

                        if len(smart_seq_sex) > 0 and len(nonsmart_seq_sex) > 0:
                            if smart_seq_sex != nonsmart_seq_sex:
                                print(f'Smart-seq and non-smart-seq scRNAseq_sex for donor ({d}) do not match - both will be included in plot.')
                                d_df = donor_sex_df[(donor_sex_df['donor_id'] == d) | (donor_sex_df['donor_id'] == d + '-smartseq')]
                                obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                                ratio_order.append((d_df['donor_id']  + ' ' + d_df['author_annotated_sex'].astype('string')).to_list())

                            if smart_seq_sex == nonsmart_seq_sex:
                                print(f'Smart-seq and non-smart-seq scRNAseq_sex for donor ({d}) match, dropping Smart-seq from plot.')
                                d_df = donor_sex_df[donor_sex_df['donor_id'] == d]
                                obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                                ratio_order.append((d_df['donor_id']  + ' ' + d_df['author_annotated_sex'].astype('string')).to_list())

                        elif len(smart_seq_sex) > 0 and len(nonsmart_seq_sex) == 0:
                            d_df = donor_sex_df[donor_sex_df['donor_id'] == d + '-smartseq']
                            smart_seq_donors_rename[f'{d}-smartseq'] = d
                            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                            ratio_order.append((d_df['donor_id'].str.split('-smartseq').str[0] + ' ' + d_df['author_annotated_sex'].astype('string')).to_list())


                        elif len(smart_seq_sex) == 0 and len(nonsmart_seq_sex) > 0:
                            d_df = donor_sex_df[donor_sex_df['donor_id'] == d]
                            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                            ratio_order.append((d_df['donor_id'] + ' ' + d_df['author_annotated_sex'].astype('string')).to_list())

                    except Exception as e:
                        print(f"Error: smart-seq and non-smart-seq sex for donor {d} were not calculated. Details: {e}")
                        obs_to_keep, ratio_order = None, None  # Set to None to indicate failure

        flattened_obs_to_keep = [obs for sublist in obs_to_keep for obs in sublist]
        flattened_ratio_order = [ro for sublist1 in ratio_order for ro in sublist1]
        adata_sub = adata[flattened_obs_to_keep, : ].copy()
        adata_sub.obs['donor_id'] = adata_sub.obs['donor_id'].astype('category')
        adata.obs['donor_id'] = adata.obs['donor_id'].str.split('-smartseq').str[0]
        donor_sex_df['donor_id'] = donor_sex_df['donor_id'].str.split('-smartseq').str[0]
        adata_sub.obs['donor_id'] = adata_sub.obs['donor_id'].cat.rename_categories(smart_seq_donors_rename)
        adata_sub.obs['donor_sex'] = adata_sub.obs.apply(lambda x: f"{x['donor_id']} {sex_map[x['sex_ontology_term_id']]}", axis=1).astype('category')
        adata_sub.var.rename(index=genes['female'], inplace=True)
        adata_sub.var.rename(index=genes['male'], inplace=True)
        f_symbs = [g for g in genes['female'].values() if g in adata_sub.var.index]
        m_symbs = [g for g in genes['male'].values() if g in adata_sub.var.index]
        dp = sc.pl.dotplot(
              adata_sub, {'female': f_symbs, 'male': m_symbs}, 'donor_sex',
              use_raw=False, categories_order=flattened_ratio_order, return_fig=True
          )

        return donor_sex_df, dp


def evaluate_var_df(adata):
    """
    Use single-cell-curation classes and fuctions and report warning/error for organism specific minimum number of gene features. Also, this function
    will look that var contains features from only a single organism.

    :param obj adata: AnnData that is being curated

    :return logging: Raises WARNING or ERROR if the number of genes in adata.var and adata.raw.var are fewer than the threshold of 40% or 60%, respectively,
    of 10x preselected biotype of genes. Will also raise ERROR if dataset contains more than a single organism; and raise a WARNING if the dataset has 50% or
    more genes filtered.
    """
    accepted_biotypes = [
        'protein_coding',
        'protein_coding_LoF',
        'lncRNA',
        'IG_C_gene',
        'IG_D_gene',
        'IG_J_gene',
        'IG_LV_gene',
        'IG_V_gene',
        'IG_V_pseudogene',
        'IG_J_pseudogene',
        'IG_C_pseudogene',
        'TR_C_gene',
        'TR_D_gene',
        'TR_J_gene',
        'TR_V_gene',
        'TR_V_pseudogene',
        'TR_J_pseudogene'
    ]

    organisms_with_descendants = [
        'NCBITaxon:9541',
        'NCBITaxon:9544',
        'NCBITaxon:10090',
        'NCBITaxon:9986',
        'NCBITaxon:9598',
        'NCBITaxon:10116',
        'NCBITaxon:9823'
    ]

    # Check that this is single organism both in metadata and var index, exit function if multiple organisms or contains invalid var features
    var_organism_objs = list({gencode.get_organism_from_feature_id(id) for id in adata.var.index.to_list()})
    if None in var_organism_objs:
        report('Features in var.index are gene symbols and/or contain deprecated Ensembl IDs', 'ERROR')
        return
    valid = True
    uns_organism = adata.uns['organism_ontology_term_id']
    var_organisms = [o.value for o in var_organism_objs]

    if 'NCBITaxon:2697049' in var_organisms:
        report('There are covid genes present in var')
        var_organisms.remove('NCBITaxon:2697049')
    if 'NCBITaxon:2697049' == uns_organism:
        report('"Covid is not a supported uns.organism"', 'ERROR')
        valid = False
    if len(var_organisms) > 1:
        report(f'Multiple organisms found in var index: {var_organisms}', 'ERROR')
        valid = False
    if valid:
        if var_organisms[0] in organisms_with_descendants:
            if utils.is_ontological_descendant_of(ONTOLOGY_PARSER,uns_organism,var_organisms[0]):
                report(f'Single organism found: {var_organisms}', 'GOOD')
            else:
                report(f'Uns metadata contains non-descendant of var index organism: {var_organisms[0]}, {uns_organism}', 'ERROR')
                return
        elif uns_organism == var_organisms[0]:
            report(f'Single organism found: {var_organisms}', 'GOOD')
        else:
            report(f'Different organisms found between var index and uns metadata: {var_organisms[0]}, {uns_organism}', 'ERROR')
            return
    else:
        return

    # Check the number of genes threshold base on biotype per specific organism
    org_obj = [i for i in gencode.SupportedOrganisms if i.value==var_organisms[0]][0]
    gene_checker = gencode.GeneChecker(org_obj)
    num_genes_biotype = len([i for i in gene_checker.gene_dict.keys() if gene_checker.gene_dict[i][2] in accepted_biotypes])

    fraction = len(adata.var.index)/num_genes_biotype 
    if fraction < 0.4:
        report(f'{len(adata.var.index)} genes present, compared against {num_genes_biotype} 10x biotype genes; fraction: {fraction} (0.40 threshold)', 'ERROR')
    elif fraction < 0.6:
        report(f'{len(adata.var.index)} genes present, compared against {num_genes_biotype} 10x biotype genes, fraction: {fraction} (0.60 threshold)','WARNING')
    else:
        report(f'{len(adata.var.index)} genes present, compared against {num_genes_biotype} 10x biotype genes; fraction: {fraction}', 'GOOD')

    # Check the number of filtered genes
    if 'feature_is_filtered' in adata.var.columns:
        num_filtered_genes = len(adata.var[adata.var.feature_is_filtered == True])
        if num_filtered_genes/len(adata.var.index) >= 0.5:
            report(f'50% or more genes are filtered', 'WARNING')
        else:
            report(f'Less than 50% of genes are filtered', 'GOOD')
    else:
        report('feature_is_filtered not found in var', 'WARNING')


def return_duplicates_dask(adata: ad.AnnData, worker_type='processes') -> pd.DataFrame:
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
            val = hashlib.sha224(data_array[indptr_array[start]:indptr_array[start + 1]].tobytes()).hexdigest()
            chunk_hashes.append(val)
            start += 1
            
        return np.array([chunk_hashes]).reshape(-1, 1)

        
    def hash_dense_chunk(matrix_chunk) -> np.array:
        chunk_hashes = [hashlib.sha224(r.tobytes()).hexdigest() for r in matrix_chunk.toarray()]
        return np.array([chunk_hashes]).reshape(-1, 1)


    hash_chunk_func = hash_dense_chunk if "in_tissue" in adata.obs.columns else hash_data_array_chunk

    with dask.config.set(scheduler=worker_type):
        hashes = map_blocks(hash_chunk_func, matrix, dtype=int).compute().ravel()

    hash_df = adata.obs.copy()
    hash_df["row_hash"] = hashes

    if "in_tissue" in hash_df.columns:
        obs_to_keep = hash_df[hash_df["in_tissue"] != 0].index
        hash_df = hash_df[hash_df.index.isin(obs_to_keep)]

    dup_df = hash_df[hash_df.duplicated(subset="row_hash", keep=False)].copy()

    return dup_df