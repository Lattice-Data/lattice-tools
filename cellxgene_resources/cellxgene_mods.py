import anndata as ad
import dask.array as da
import h5py
import json
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import os
import pandas as pd
import re
import scanpy as sc
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from scipy import sparse
from cellxgene_ontology_guide.ontology_parser import OntologyParser
import cellxgene_schema.gencode as gencode
import cellxgene_schema.utils as utils
import cellxgene_schema.schema as schema


EXPECTED_BARCODES = {
    'EFO:0009901':                    '3pv1',
    "10x 3' v1":                      '3pv1',
    'EFO:0009899':                    '3pv2_5pv1_5pv2',
    "10x 3' v2":                      '3pv2_5pv1_5pv2',
    'EFO:0009922':                    '3pv3',
    "10x 3' v3":                      '3pv3',
    'EFO:0022604':                    '3pv4',
    "10x 3' v4":                      '3pv4',
    'EFO:0011025':                    '3pv2_5pv1_5pv2',
    "10x 5' v1":                      '3pv2_5pv1_5pv2',
    'EFO:0009900':                    '3pv2_5pv1_5pv2',
    "10x 5' v2":                      '3pv2_5pv1_5pv2',
    'EFO:0030004':                    '3pv2_5pv1_5pv2',
    "10x 5' transcription profiling": '3pv2_5pv1_5pv2',
    'EFO:0022605':                    '5pv3',
    "10x 5' v3":                      '5pv3',
    'EFO:0030059':                    'multiome',
    "10x multiome":                   'multiome',
    'EFO:0920134':                    'multiome',
    "10x GEM-X Epi Multiome":         'multiome',
    'EFO:0920135':                    'multiome',
    "10x Next-GEM Multiome":          'multiome',
    'EFO:0920086':                    'flex_v1',
    "10x gene expression flex v1":    'flex_v1',
    'EFO:0920088':                    'flex_v1',
    "10x GEM-X Flex v1":              'flex_v1',
    'EFO:0920087':                    'flex_v1',
    "10x Next GEM Flex v1":           'flex_v1',
    'EFO:0920089':                    'flex_v2',
    "10x Flex Apex":                  'flex_v2'
}

FLEX_ASSAYS = [
    'EFO:0022606','EFO:0920089','EFO:0920086','EFO:0920088','EFO:0920087'
]

OBS_ONTOLOGY_LABELS_REQUIRED = [
    'assay', 'cell_type', 'development_stage', 'disease',
    'self_reported_ethnicity', 'sex', 'tissue'
]
OBS_ONTOLOGY_LABELS_OPTIONAL = ['experimental_condition']
OBS_ONTOLOGY_LABELS = OBS_ONTOLOGY_LABELS_REQUIRED + OBS_ONTOLOGY_LABELS_OPTIONAL

OBS_ONTOLOGY_IDS_REQUIRED = [
    f"{label}_ontology_term_id" for label in OBS_ONTOLOGY_LABELS_REQUIRED
]
OBS_ONTOLOGY_IDS_OPTIONAL = [
    f"{label}_ontology_term_id" for label in OBS_ONTOLOGY_LABELS_OPTIONAL
]

OBS_NON_ONTOLOGY_CURATOR_REQUIRED = [
    'donor_id', 'suspension_type', 'tissue_type', 'is_primary_data'
]
OBS_NON_ONTOLOGY_CURATOR_OPTIONAL = [
    'genetic_perturbation_id', 'genetic_perturbation_strategy'
]

OBS_NON_ONTOLOGY_PORTAL_REQUIRED = ['observation_joinid']
OBS_NON_ONTOLOGY_PORTAL_OPTIONAL = ['perturbation_types']

OBS_PORTAL_REQUIRED = OBS_ONTOLOGY_LABELS_REQUIRED + OBS_NON_ONTOLOGY_PORTAL_REQUIRED
OBS_PORTAL_OPTIONAL = OBS_ONTOLOGY_LABELS_OPTIONAL + OBS_NON_ONTOLOGY_PORTAL_OPTIONAL
OBS_PORTAL_ALL = OBS_PORTAL_REQUIRED + OBS_PORTAL_OPTIONAL
OBS_CURATOR_REQUIRED = OBS_ONTOLOGY_IDS_REQUIRED + OBS_NON_ONTOLOGY_CURATOR_REQUIRED
OBS_CURATOR_OPTIONAL = OBS_ONTOLOGY_IDS_OPTIONAL + OBS_NON_ONTOLOGY_CURATOR_OPTIONAL
OBS_CURATOR_ALL = OBS_CURATOR_REQUIRED + OBS_CURATOR_OPTIONAL

OBS_FULL_STANDARDS = OBS_PORTAL_ALL + OBS_CURATOR_ALL

UNS_PORTAL_REQUIRED = [
    'citation', 'is_pre_analysis', 'schema_reference', 'schema_version', 'organism'
]

UNS_CURATOR_REQUIRED = ['title', 'organism_ontology_term_id']

VAR_PORTAL_REQUIRED = [
    'feature_name', 'feature_reference', 'feature_biotype', 'feature_type', 'feature_length'
]

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
                

    from src.collection import (
        create_collection,create_revision,delete_collection,get_collection,
        get_collections,get_collection_version,get_collection_versions,update_collection
    )
    from src.dataset import (
        create_dataset,delete_dataset,get_dataset,get_datasets,get_dataset_manifest,
        get_dataset_version,get_dataset_versions,
        upload_local_datafile,upload_datafiles_from_manifest
    )

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
        'ERROR': '\033[31m', #red
        'code': '\033[30;48;5;252m' #grey background
    }
    if level:
        c = colors[level]
        if level not in ['code']:
            mess = f'{level}: {mess}'
        print(f'\033[1m{c}{mess}\033[0m')
    else:
        print(mess)


def revise_genetic_perturbations(uns):
    fields_to_remove = ['derived_genomic_regions', 'derived_features']

    for guide_meta in uns['genetic_perturbations'].values():
        if 'intended_features' in guide_meta:
            guide_meta['intended_features'] = {
                feature_id: '' for feature_id in guide_meta['intended_features'].keys()
            }

        for field in fields_to_remove:
            guide_meta.pop(field, None)

    return uns


def revise_cxg(adata):
    for p in UNS_PORTAL_REQUIRED:
        del adata.uns[p]

    if 'genetic_perturbations' in adata.uns:
        adata.uns = revise_genetic_perturbations(adata.uns)

    adata.obs.drop(columns=[c for c in OBS_PORTAL_ALL if c in adata.obs.columns], inplace=True)
    adata.var.drop(columns=VAR_PORTAL_REQUIRED, inplace=True)

    if adata.raw:
        adata.raw.var.drop(columns=VAR_PORTAL_REQUIRED, inplace=True)

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


def determine_sparsity(x):
    """Calculate sparsity of a matrix."""
    if isinstance(x, (sparse.coo_matrix, sparse.csr_matrix, sparse.csc_matrix)):
        sparsity = 1 - x.count_nonzero() / float(np.prod(x.shape))
    elif isinstance(x, np.ndarray):
        sparsity = 1 - np.count_nonzero(x) / float(np.prod(x.shape))
    else:
        report(f'matrix is of type {type(x)}, sparsity calculation has not been implemented', 'WARNING')
        return None

    return round(sparsity, 3)


def evaluate_sparsity(adata, max_sparsity=0.5):
    """Check sparsity and recommend sparse format conversion if needed."""
    valid = True

    # Check X
    sparsity = determine_sparsity(adata.X)
    report(f'.X sparsity: {sparsity}')
    if sparsity and sparsity > max_sparsity and not isinstance(adata.X, sparse.csr_matrix):
        report('X should be converted to csr sparse', 'ERROR')
        report('adata.X = sparse.csr_matrix(adata.X)', 'code')
        valid = False

    # Check raw.X
    if adata.raw:
        sparsity = determine_sparsity(adata.raw.X)
        report(f'.raw.X sparsity: {sparsity}')
        if sparsity and sparsity > max_sparsity and not isinstance(adata.raw.X, sparse.csr_matrix):
            report('raw.X should be converted to csr sparse', 'ERROR')
            report(
                'adata.raw = ad.AnnData(sparse.csr_matrix(adata.raw.X), var=adata.raw.var, obs=adata.obs)',
                'code'
            )
            valid = False

    # Check layers
    for layer_name in adata.layers:
        sparsity = determine_sparsity(adata.layers[layer_name])
        report(f'layers[{layer_name}] sparsity: {sparsity}')
        if sparsity and sparsity > max_sparsity and not isinstance(adata.layers[layer_name], sparse.csr_matrix):
            report(f'layers[{layer_name}] should be converted to csr sparse', 'ERROR')
            report(f"adata.layers['{layer_name}'] = sparse.csr_matrix(adata.layers['{layer_name}'])", 'code')
            valid = False

    if valid:
        report('all matrices have passed sparsity checks', 'GOOD')


def get_raw_matrix_info(adata):
    """
    Get raw count matrix and its location.
    Returns (matrix, location_string, is_csr_sparse)
    """
    if adata.raw:
        matrix = adata.raw.X
        location = '.raw.X'
    else:
        matrix = adata.X
        location = '.X'

    return matrix, location


def evaluate_raw_matrix(matrix, loc):
    """Validate raw count matrix properties."""
    report(f'raw counts determined to be in {loc}')

    # Check if all values are integers
    # For sparse matrices, only check the data array
    data = matrix.data if hasattr(matrix, 'data') else matrix
    all_integers = np.array_equal(data, np.round(data))

    if all_integers:
        report('raw counts are all integers', 'GOOD')
        if matrix.dtype != np.float32:
            report(f'raw count dtype should be float32, not {matrix.dtype}', 'ERROR')
            if loc == '.raw.X':
                report(
                    "adata.raw = ad.AnnData(sparse.csr_matrix(adata.raw.X.astype('float32')), var=adata.raw.var, obs=adata.obs)",
                    'code'
                )
            else:
                report(
                    "adata.X = sparse.csr_matrix(adata.X.astype('float32'))",
                    'code'
                )
        else:
            report('raw count dtype is float32', 'GOOD')
    else:
        report('raw counts contain non-integer values', 'ERROR')


def get_matrix_range(matrix):
    """Get min and max values from a matrix."""
    return matrix.min(), matrix.max()


def evaluate_data(adata):
    evaluate_sparsity(adata)
    print()
    evaluate_data_range(adata)


def check_matrix_duplicates(matrix_pairs):
    """
    Check if matrices are truly identical.
    Uses fast checks, then goes straight to full comparison.

    Args:
        matrix_pairs: List of (name, matrix) tuples

    Returns:
        List of groups that are true duplicates, or empty list
    """
    if len(matrix_pairs) < 2:
        return []

    # Level 1: Check shapes (instant)
    shapes = [(name, mx.shape) for name, mx in matrix_pairs]
    if len(set(s for _, s in shapes)) > 1:
        return []  # Different shapes, can't be duplicates

    # Level 2: Check sum (very fast)
    sums = [(name, mx.sum()) for name, mx in matrix_pairs]
    sum_groups = {}
    for name, s in sums:
        sum_groups.setdefault(s, []).append(name)

    duplicate_groups = []

    for sum_val, names in sum_groups.items():
        if len(names) < 2:
            continue  # Only one matrix with this sum

        # Get matrices with matching sums
        matching = [(name, mx) for name, mx in matrix_pairs if name in names]

        # Full comparison - they passed the quick checks
        if len(matching) == 2:
            name1, mx1 = matching[0]
            name2, mx2 = matching[1]

            if matrices_equal(mx1, mx2):
                duplicate_groups.append([name1, name2])
        else:
            # For 3+ matrices, compare all pairs
            for i in range(len(matching)):
                for j in range(i + 1, len(matching)):
                    name_i, mx_i = matching[i]
                    name_j, mx_j = matching[j]

                    if matrices_equal(mx_i, mx_j):
                        # Find or create group containing these duplicates
                        found_group = None
                        for group in duplicate_groups:
                            if name_i in group or name_j in group:
                                found_group = group
                                break

                        if found_group:
                            if name_i not in found_group:
                                found_group.append(name_i)
                            if name_j not in found_group:
                                found_group.append(name_j)
                        else:
                            duplicate_groups.append([name_i, name_j])

    return duplicate_groups


def matrices_equal(mx1, mx2):
    """
    Check if two matrices are exactly equal.
    Handles both sparse and dense matrices.
    """
    # Check if both are sparse or both are dense
    mx1_sparse = isinstance(mx1, (sparse.spmatrix, sparse.sparray))
    mx2_sparse = isinstance(mx2, (sparse.spmatrix, sparse.sparray))

    if mx1_sparse != mx2_sparse:
        return False

    if mx1_sparse:
        # For sparse matrices, compare in CSR format
        mx1_csr = mx1.tocsr() if mx1.format != 'csr' else mx1
        mx2_csr = mx2.tocsr() if mx2.format != 'csr' else mx2

        # Compare data, indices, and indptr arrays
        return (np.array_equal(mx1_csr.data, mx2_csr.data) and
                np.array_equal(mx1_csr.indices, mx2_csr.indices) and
                np.array_equal(mx1_csr.indptr, mx2_csr.indptr))
    else:
        # For dense matrices
        return np.array_equal(mx1, mx2)


def _get_matrix_by_name(adata, name):
    """Helper to retrieve matrix by string name."""
    if name == '.X':
        return adata.X
    elif name == '.raw.X':
        return adata.raw.X
    elif name.startswith('layers['):
        layer_name = name[7:-1]
        return adata.layers[layer_name]
    else:
        raise ValueError(f"Unknown matrix name: {name}")


def evaluate_data_range(adata):
    """Check data ranges and detect potential duplicate layers."""
    min_maxs = {}

    # Determine where raw counts are
    if adata.raw:
        raw_min, raw_max = get_matrix_range(adata.raw.X)
        report(f'.raw.X min = {raw_min}')
        report(f'.raw.X max = {raw_max}')
        min_maxs['.raw.X'] = (raw_min, raw_max)
        raw_matrix = adata.raw.X
        raw_loc = '.raw.X'
    else:
        raw_matrix = adata.X
        raw_loc = '.X'

    # Check X
    x_min, x_max = get_matrix_range(adata.X)
    report(f'.X min = {x_min}')
    report(f'.X max = {x_max}')
    min_maxs['.X'] = (x_min, x_max)

    # Check layers
    for layer_name in adata.layers:
        layer_min, layer_max = get_matrix_range(adata.layers[layer_name])
        report(f'layers[{layer_name}] min = {layer_min}')
        report(f'layers[{layer_name}] max = {layer_max}')
        min_maxs[f'layers[{layer_name}]'] = (layer_min, layer_max)

    # Detect potential duplicates based on min/max
    range_groups = {}
    for name, range_val in min_maxs.items():
        range_groups.setdefault(range_val, []).append(name)

    potential_duplicates = [names for names in range_groups.values() if len(names) > 1]
    if potential_duplicates:
        for dup_group in potential_duplicates:
            report(f'possible redundant layers based on min/max: {dup_group}. Checking...', 'WARNING')

            # Get the matrices for this group
            group_matrices = [(name, _get_matrix_by_name(adata,name)) for name in dup_group]

            # Check if they're truly identical
            true_duplicates = check_matrix_duplicates(group_matrices)

            if true_duplicates:
                report(
                    f'Confirmed duplicates: {true_duplicates}\nRemove duplication to reduce object and file size',
                    'ERROR'
                )
            else:
                report('Different matrices (same min/max is coincidental)', 'WARNING')

    print()

    # Validate raw matrix
    evaluate_raw_matrix(raw_matrix, raw_loc)


def evaluate_uns_colors(adata):
    colors_keys = [k for k in adata.uns.keys() if k.endswith('_colors')]
    if colors_keys:
        for k in colors_keys:
            obs_field = k[:-(len('_colors'))]

            if obs_field in OBS_ONTOLOGY_LABELS:
                report(f'uns.{k} not allowed, move to uns.{obs_field}_ontology_term_id_colors', 'ERROR')
            elif obs_field not in adata.obs.columns:
                report(f'{obs_field} not found in obs, consider DELETING or RENAMING uns.{k}', 'ERROR')
            elif adata.obs.dtypes[obs_field].name != 'category':
                report(f'uns.{k} is associated with non-categorical {obs_field}', 'ERROR')
            else:
                colors = len(adata.uns[k])
                values = len(adata.obs[obs_field].cat.categories.values)
                if colors < values:
                    report(f'uns.{k} has only {str(colors)} colors but obs.{obs_field} has {str(values)} values', 'ERROR')
                else:
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

    module_dir = os.path.dirname(os.path.abspath(__file__))
    approved_file = os.path.join(module_dir, 'ref_files', 'genes_approved.csv.gz')
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


def extract_barcodes(index):
    pattern = re.compile(r'[ACTG]{12,}')
    barcodes = []
    affixes = []

    for i in index:
        m = pattern.search(str(i))
        if m:
            barcode = m.group()[:16]
            barcodes.append(barcode)
            affixes.append(i.replace(barcode,''))
        else:
            barcodes.append(None)
            affixes.append(None)

    if not any(barcodes):
        report('No barcodes found in obs.index', 'WARNING')

    return barcodes, affixes


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
    obs[['barcode', 'affix']] = pd.DataFrame(
        zip(*extract_barcodes(obs.index)),
        index=obs.index
    )
    if len(set(ref_df.index.to_list()).intersection(set(obs['barcode'].to_list()))) == 0:
        report('Did not find any barcodes in obs index, cannot evaluate barcodes', 'WARNING')
        return
    obs = obs.merge(ref_df[['summary']],on='barcode',how='left').set_index(obs.index)
    obs['summary'] = obs.apply(
        lambda x: no_barcode_v if pd.isna(x['barcode']) else (f"{len(x['barcode'])}nt" if pd.isna(x['summary']) else x['summary']),
        axis=1
    )

    return obs


def validate_barcode_assignments(df_summary):
    """
    Check for unexpected barcode assignments based on assay type.
    Prints warnings when barcodes don't match expected patterns.
    """
    has_unexpected = False

    for i,row in df_summary.iterrows():
        if i not in EXPECTED_BARCODES:
            continue
        ignore_cols = {'multiple', no_barcode_v, EXPECTED_BARCODES[i]}

        # Check all barcode columns
        for col in df_summary.columns:
            count = row[col]
            if count > 0 and col not in ignore_cols:
                report(f'{col} barcodes marked as {i}','ERROR')
                has_unexpected = True

    if has_unexpected:
        print()


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

    validate_barcode_assignments(df)

    df = df[[c for c in df if df[c].sum() > 0 and c not in ['multiple',no_barcode_v] and not c.endswith('nt')]
            + [c for c in df if df[c].sum() > 0 and c.endswith('nt')]
            + [c for c in df if df[c].sum() == 0 and c not in ['multiple',no_barcode_v]]
            + [c for c in ['multiple',no_barcode_v] if c in df]]
    df.sort_values(list(df.columns), ascending=False, inplace=True)

    return df


def evaluate_obsm(adata, labels=None):
    keys = adata.obsm_keys()

    cellpop_field = 'cell_type' if labels else 'cell_type_ontology_term_id'
    colors_key = f'{cellpop_field}_colors'
    had_colors = colors_key in adata.uns

    plot = False
    sc.set_figure_params(dpi=100)
    for e in keys:
        if e.startswith('X_'):
            sc.pl.embedding(adata, basis=e, color=cellpop_field, legend_loc='on data')
            plot = True
        elif e == 'spatial':
            if np.isnan(adata.obsm['spatial']).any():
                report("obsm[spatial] contains nans", 'ERROR')
            sc.pl.embedding(adata, basis=e, color=cellpop_field, legend_loc='on data')
            plot = True
        else:
            report(f'{e} will not be plotted')

    if not had_colors and colors_key in adata.uns:
        del adata.uns[colors_key]

    if not plot:
        report('No visualizable embeddings in obsm', 'ERROR')

    de = adata.uns.get('default_embedding')
    if de:
        if de not in adata.obsm_keys():
            report(f'uns.default_embedding:{de} not in [{",".join(adata.obsm.keys())}]', 'ERROR')
        else:
            report(f'uns.default_embedding:{de} is in [{",".join(adata.obsm.keys())}]', 'GOOD')


def evaluate_uns_schema(uns, labels=False):
    for f in UNS_CURATOR_REQUIRED:
        if f in uns:
            report(f'{f}: {uns[f]}')
        else:
            report(f'{f} is required', 'ERROR')
    if not labels:
        for f in UNS_PORTAL_REQUIRED:
            if f in uns:
                report(f'{f} should not be present in uns', 'ERROR')


def evaluate_obs_schema(obs, labels=False):
    if labels:
        for o in OBS_PORTAL_REQUIRED:
            if o in ['observation_joinid']:
                continue

            if o in obs.columns:
                report(f'{o} {obs[o].unique().tolist()}\n')
            else:
                report(f'{o} not in obs\n', 'ERROR')
        for o in OBS_PORTAL_OPTIONAL:
            if o in obs.columns:
                report(f'{o} {obs[o].unique().tolist()}\n')
    else:
        for o in OBS_CURATOR_REQUIRED:
            if o in obs.columns:
                report(f'{o} {obs[o].unique().tolist()}\n')
            else:
                report(f'{o} not in obs\n', 'ERROR')
        for o in OBS_ONTOLOGY_IDS_OPTIONAL:
            if o in obs.columns:
                report(f'{o} {obs[o].unique().tolist()}\n')
        for o in OBS_ONTOLOGY_LABELS:
            if o in obs.columns:
                report(f'schema conflict - {o} in obs\n', 'ERROR')
    for o in OBS_NON_ONTOLOGY_CURATOR_OPTIONAL:
        if o in obs.columns:
            if o == 'genetic_perturbation_id':
                uniq_vals = obs[o].unique().tolist()
                report(f'{o} <{len(uniq_vals)} unique values>, {uniq_vals[:5]}...\n')
            else:
                report(f'{o} {obs[o].unique().tolist()}\n')
    if 'cell_type_ontology_term_id' in obs.columns and 'unknown' in obs['cell_type_ontology_term_id'].unique():
        if 'in_tissue' in obs.columns:
            num_unknown = obs.loc[(obs['in_tissue']==1) & (obs['cell_type_ontology_term_id']=='unknown')].shape[0]
            perc_unknown = round(100*(num_unknown/obs.loc[obs['in_tissue']==1].shape[0]), 1)
        else:
            num_unknown = obs[obs['cell_type_ontology_term_id']=='unknown'].shape[0]
            perc_unknown = round(100*(num_unknown/obs.shape[0]), 1)
        if num_unknown > 20:
            report(
                f'{num_unknown} ({perc_unknown}%) cells are cell_type:unknown.\n'
                'Some unknowns are acceptable but confirm there is neither an appropriate CL term nor a term to request',
                'WARNING'
            )

    for o in obs.columns:
        if o not in OBS_FULL_STANDARDS and '_'.join(o.split()).lower() in OBS_FULL_STANDARDS:
            report(f'"close enough" schema conflict: suggest renaming obs.{o}\n', 'ERROR')


def evaluate_obs(obs):
    long_fields = []
    gradient_fields = []
    uber_dict = {}
    for o in obs.columns:
        vc_dict = obs[o].value_counts(dropna=False).to_dict()
        counts = '_'.join([str(c) for c in vc_dict.values()])
        count_len = len(vc_dict.keys())
        values = [str(i) for i in vc_dict.keys()]

        if o.startswith(' ') or o.endswith(' ') or '  ' in o:
            report(f'leading/trailing whitespace: {o}\n')

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
            if len(v) > 1 and not all(elem in OBS_FULL_STANDARDS for elem in props):
                report(f'possible redundancy: {[e["property"] for e in v]}\n')

    if gradient_fields:
        report(f'continuous fields: {gradient_fields}\n')
    if long_fields:
        report(f'long fields: {long_fields}')


def ensure_canonical_csr(matrix, location_desc):
    """Ensure matrix is in canonical CSR format."""
    if not isinstance(matrix, sparse.csr_matrix):
        report(
            f'{location_desc} not in sparse CSR format, conversion required, rerun evaluate_data() for guidance',
            'ERROR'
        )
        return None

    if not matrix.has_canonical_format:
        report(f"{location_desc} not in canonical format, converting now...")
        original_nnz = matrix.nnz
        matrix.sort_indices()
        matrix.sum_duplicates()
        if original_nnz != matrix.nnz:
            report(f"{original_nnz - matrix.nnz} duplicates found during canonical conversion")

    assert matrix.has_canonical_format, f"{location_desc} still in non-canonical format"
    return matrix


def hash_sparse_rows(matrix, obs_df):
    """
    Hash rows of a sparse CSR matrix to detect duplicates.

    Returns DataFrame with only duplicated rows and their hash values.
    """
    data_array = matrix.data
    index_array = matrix.indices
    indptr_array = matrix.indptr

    # First pass: hash data arrays for all rows
    data_hashes = []
    for i in range(matrix.shape[0]):
        row_data = data_array[indptr_array[i]:indptr_array[i + 1]]
        data_hashes.append(hash(row_data.tobytes()))

    # Create working dataframe with data hashes
    hash_df = obs_df.copy()
    hash_df['data_array_hash'] = data_hashes

    # Keep only rows with duplicate data hashes
    hash_df = hash_df[hash_df.duplicated(subset='data_array_hash', keep=False)]

    if hash_df.empty:
        return hash_df

    hash_df.sort_values('data_array_hash', inplace=True)

    # Second pass: hash index arrays for potential duplicates
    def hash_row_indices(obs_index):
        obs_loc = obs_df.index.get_loc(obs_index)
        row_indices = index_array[indptr_array[obs_loc]:indptr_array[obs_loc + 1]]
        return hash(row_indices.tobytes())

    hash_df['index_array_hash'] = hash_df.index.map(hash_row_indices)

    # Keep only true duplicates (both data and indices match)
    hash_df = hash_df[
        hash_df.duplicated(subset=['data_array_hash', 'index_array_hash'], keep=False)
    ]

    return hash_df


def evaluate_dup_counts(adata):
    """
    Detect duplicate raw count rows in the dataset.

    Returns DataFrame of duplicated rows if found, None otherwise.
    """
    # Filter to in-tissue observations for spatial data
    working_adata = adata
    if 'in_tissue' in adata.obs.columns:
        obs_to_keep = adata.obs['in_tissue'] != 0
        working_adata = adata[obs_to_keep, :]
        report(f'Filtered to {obs_to_keep.sum()} in-tissue observations')

    # Get the raw count matrix
    matrix, loc_desc = get_raw_matrix_info(working_adata)

    # Ensure matrix is in canonical CSR format
    matrix = ensure_canonical_csr(matrix, loc_desc)
    if matrix is None:
        return None

    # Hash rows to find duplicates
    dup_df = hash_sparse_rows(matrix, working_adata.obs)
    if not dup_df.empty:
        report(f'Found {len(dup_df)} rows with duplicated raw counts', 'ERROR')
        return dup_df

    report('No duplicated raw counts', 'GOOD')
    return None


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
        found_approved = False
        found_var = False
        if s in approved['symbol_only'].tolist():
            found_approved = True
            ensg_ids = approved.loc[approved['symbol_only'] == s, 'feature_id']
            for ensg_id in ensg_ids:
                if ensg_id in var.index:
                    ensg_list.append(ensg_id)
                    report(f'{ensg_id} -- {s}')
                    found_var = True
        if not found_var:
            s_lower = s[0] + s[1:].lower()
            if s_lower in approved['symbol_only'].tolist():
                found_approved = True
                ensg_ids = approved.loc[approved['symbol_only'] == s_lower, 'feature_id']
                for ensg_id in ensg_ids:
                    if ensg_id in var.index:
                        ensg_list.append(ensg_id)
                        report(f'{ensg_id} -- {s_lower}')
                        found_var = True
        if not found_approved:
            report(f'{s} not found in genes_approved.csv.gz, check for typos', 'WARNING')
        elif not found_var:
            report(f"{s}/{','.join(ensg_ids)} not found in var", 'WARNING')

    return ensg_list


def anndata_to_spatialdata_visium(adata, library_id, cellpop_field):
    '''
    Convert Visium AnnData object to SpatialData object with proper coordinate transformations.
    
    Due to complex package dependencies that prevent installing the below packages on JupyterHub,
    this function will import only in this scope to allow the rest of cellxgene_mods to work
    without issue.
    '''
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=FutureWarning, module='dask.dataframe')
            import geopandas as gpd
            import spatialdata as sd
            import spatialdata_plot
            from shapely.geometry import Point
            from spatialdata.transformations import Identity
    except ImportError as e:
        print(f"Cannot plot spatial data due to import error: {e}")
        print("Please create local conda env according to lattice-tools readme")
        return

    # Extract spatial coordinates and scale factors
    # Extract scalefactors and coordinates
    scalefactors = adata.uns['spatial'][library_id]['scalefactors']
    tissue_hires_scalef = scalefactors['tissue_hires_scalef']
    spot_radius_fullres = scalefactors['spot_diameter_fullres'] / 2
    spot_radius_hires = spot_radius_fullres * tissue_hires_scalef

    coords_fullres = adata.obsm['spatial'].copy()
    coords_hires = coords_fullres * tissue_hires_scalef

    shapes_dfs = {}
    radii_and_coords = [
        # hires first, since will always be present
        (spot_radius_hires, coords_hires, 'hires'),
        (spot_radius_fullres, coords_fullres, 'fullres'),
    ]

    for spot_radius, coords, name in radii_and_coords:
        circles = [Point(x, y).buffer(spot_radius) for x, y in coords]
        shapes_df = gpd.GeoDataFrame({
            'geometry': circles,
            'in_tissue': adata.obs['in_tissue'].values,
            'spot_id': adata.obs.index
        })
        if cellpop_field:
            shapes_df[cellpop_field] = adata.obs[cellpop_field].values
        shapes_dfs[name] = shapes_df

    # Process images
    images = {}
    shapes = {}
    visium_images = adata.uns['spatial'][library_id]['images']
    images_to_process = ['hires', 'fullres'] if 'fullres' in visium_images else ['hires']

    for image in images_to_process:
        img = visium_images[image]
        if len(img.shape) == 3:
            img = np.transpose(img, (2, 0, 1))
        dask_array = da.from_array(img, chunks=img.shape)
        images[image] = sd.models.Image2DModel.parse(
            dask_array, 
            transformations={'global': Identity()}
        )
        shapes_img = sd.models.ShapesModel.parse(
            shapes_dfs[image], 
            transformations={'global': Identity()}
        )
        shapes[f'{library_id}_{image}'] = shapes_img

    # Create table (linked to hires by default)
    table_obs = adata.obs.copy()
    table_obs.drop(columns=[c for c in table_obs.columns if c not in ['in_tissue',cellpop_field]],inplace=True)
    table_obs['region'] = f'{library_id}_hires'
    table_obs['region'] = table_obs['region'].astype('category')
    table_obs['instance_key'] = range(len(adata.obs))

    table_adata = adata.copy()
    table_adata.obs = table_obs

    table = sd.models.TableModel.parse(
        table_adata,
        region=f'{library_id}_hires',
        region_key='region',
        instance_key='instance_key'
    )

    # Create SpatialData object
    sdata = sd.SpatialData(
        images=images,
        shapes=shapes,
        tables=table
    )
    sdata.attrs['scalefactors'] = scalefactors

    return sdata


def visualize_spatial(sdata, library_id, cellpop_field):
    viz_spatial_per_res(sdata, library_id, 'hires', cellpop_field)
    if f'{library_id}_fullres' in sdata.shapes:
        viz_spatial_per_res(sdata, library_id, 'fullres', cellpop_field)
    else:
        report('fullres image is absent - strongly recommended', 'WARNING')


def viz_spatial_per_res(sdata, library_id, res, cellpop_field):
    viz_spatial_per_field(sdata, library_id, res, 'in_tissue')
    if cellpop_field:
        viz_spatial_per_field(sdata, library_id, res, cellpop_field)


def viz_spatial_per_field(sdata, library_id, res, field):
    ncols = 2
    nrows = 1
    figsize = 4
    wspace = 0.5
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
    )
    plt.subplots_adjust(wspace=wspace)

    # Left: Image only
    sdata.pl.render_images(res).pl.show(ax=axes[0])
    axes[0].set_title(f'{res} image', fontsize=12)
    axes[0].axis('off')

    # Right: Image first, then add points
    sdata.pl.render_images(res).pl.show(ax=axes[1])
    sdata.pl.render_shapes(
        f'{library_id}_{res}',
        color=field
    ).pl.show(ax=axes[1])
    axes[1].set_title(f'{res} image + {field}', fontsize=12)
    axes[1].axis('off')

    plt.show()


def plot_vis(adata, cellpop_field=None):
    library_id = [k for k in adata.uns['spatial'].keys() if k != 'is_single'][0]
    sdata = anndata_to_spatialdata_visium(adata, library_id, cellpop_field)
    visualize_spatial(sdata, library_id, cellpop_field)


def evaluate_spatial(adata, cellpop_field):
    if 'spatial' not in adata.uns:
        report('required uns[spatial] is absent', 'ERROR')
        return
    if 'is_single' not in adata.uns['spatial']:
        report('required uns[spatial][is_single] is absent', 'ERROR')
        return
    if adata.uns['spatial']['is_single'] == True:
        if len(adata.uns['spatial']) != 2:
            report(
                'uns[spatial] keys should be is_single + exactly 1 library_id\n'
                f"keys: {', '.join(adata.uns['spatial'].keys())}",
                'ERROR'
            )
            return
        if 'spatial' not in adata.obsm:
            report('required obsm[spatial] is absent', 'ERROR')
            return
        plot_vis(adata, cellpop_field)
    elif len(adata.uns['spatial'].keys()) > 1:
        report(
            'uns[spatial] keys should be only is_single\n'
            f"keys: {', '.join(adata.uns['spatial'].keys())}",
            'ERROR'
        )


def side_by_side_dotplot(adata, gene_list, groupby):
    if not gene_list:
        report('No genes in list', 'ERROR')
        return
    panels = [(False, ".X")] + ([(True, ".raw.X")] if adata.raw else [])
    n = len(panels)

    n_groups = adata.obs[groupby].nunique()
    top_in, gap_in, legend_in, bottom_in = 0.4, 1.3, 1.3, 0.15
    main_in = max(2.0, n_groups * 0.28)
    fig_h = top_in + main_in + gap_in + legend_in + bottom_in
    fig = plt.figure(figsize=(len(gene_list) * 0.3 * n, fig_h))
    outer = fig.add_gridspec(1, n, wspace=0)

    main_y0, main_h = (bottom_in + legend_in + gap_in) / fig_h, main_in / fig_h
    legend_y0, legend_h = bottom_in / fig_h, legend_in / fig_h

    for i, (use_raw, title) in enumerate(panels):
        col = outer[0, i].get_position(fig)
        main_gs = fig.add_gridspec(1, 1, left=col.x0, right=col.x1, bottom=main_y0, top=main_y0 + main_h)
        main_ax = fig.add_subplot(main_gs[0, 0])
        legend_ax = fig.add_axes([col.x0, legend_y0, col.width, legend_h])

        dp = sc.pl.dotplot(adata, gene_list, groupby=groupby, use_raw=use_raw,
                            ax=main_ax, show=False, return_fig=True, title=title)
        dp.legend(show=False)
        dp.make_figure()
        mainplot_ax = dp.get_axes()["mainplot_ax"]
        mainplot_ax.set_position(main_ax.get_position())
        if i > 0:
            mainplot_ax.tick_params(axis="y", left=False, labelleft=False)

        legend_ax.axis("off")
        lp = legend_ax.get_position()
        size_ax = fig.add_axes([lp.x0, lp.y0, lp.width * 0.35, lp.height * 0.7])
        cbar_ax = fig.add_axes([lp.x0 + lp.width * 0.55, lp.y0 + lp.height * 0.15, lp.width * 0.4, lp.height * 0.28])
        dp._plot_size_legend(size_ax)
        norm = Normalize(vmin=dp.dot_color_df.values.min(), vmax=dp.dot_color_df.values.max())
        dp._plot_colorbar(cbar_ax, norm)

    plt.show()


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
            print('Donors with < 100 total counts dropped:', ','.join(donors_to_remove))
            male_female_df.drop(male_female_df[male_female_df.total_sum < 100].index, inplace=True)
        #Calculate ratio and assign sex
        male_female_df['male_to_female'] = male_female_df['male_sum']/male_female_df['female_sum']
        male_female_df['scRNAseq_sex'] = male_female_df.apply(lambda x: assign_sex(x['male_to_female']), axis=1)

        return male_female_df,donors_to_remove

    except Exception as e:
        print(e)

def compare_donor_sex(df):
    inconsistencies = df[df['scRNAseq_sex'] != df['author_annotated_sex']].sort_values('donor_id')
    if inconsistencies.empty:
        report('donor sex metadata is consistent', 'GOOD')
    else:
        curated_unknowns = inconsistencies[inconsistencies['author_annotated_sex'] == 'unknown']
        if not curated_unknowns.empty:
            report(
                'donors with annoted sex:unknown have sex indicated by expression\n'
                'if the contributor is a study author (not reuse), provide the plots and\n'
                'ask if they would like to update the donor sex annotated based on this analysis',
                'WARNING'
            )
            display(curated_unknowns)

        expression_unknowns = inconsistencies[inconsistencies['scRNAseq_sex'] == 'unknown']
        if not expression_unknowns.empty:
            report(f'{len(expression_unknowns)} donors with undetermined sex by expression')
            display(expression_unknowns)

        true_inconsistencies = inconsistencies[
            (
                (inconsistencies['author_annotated_sex'] == 'male') &
                (inconsistencies['scRNAseq_sex'] == 'female')
            ) | (
                (inconsistencies['author_annotated_sex'] == 'female') &
                (inconsistencies['scRNAseq_sex'] == 'male')
            )
        ]
        if not true_inconsistencies.empty:
            report(
                'donor sex metadata inconsistencies\n'
                'the reported sex should be double-checked for these donors in the associated publication\n'
                'if the contributor is a study author (not reuse), ask them to double-check their records',
                'ERROR'
            )
            display(true_inconsistencies)


def evaluate_donors_sex(adata):
    if 'NCBITaxon:9606' != adata.uns['organism_ontology_term_id']:
        report('Cannot calculate sex for non-human data.')
        return None
    else:
        genes_file = 'ref_files/sex_analysis_genes.json'
        genes = json.load(open(genes_file))
        female_ids = genes['female'].keys()
        male_ids = genes['male'].keys()
        metadata_list = ['donor_id', 'sex_ontology_term_id','assay_ontology_term_id']
        smart_assay_list = [
            'EFO:0010184','EFO:0008931','EFO:0008930','EFO:0010022',
            'EFO:0700016','EFO:0022488','EFO:0008442'
        ]
        adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
        mask = adata.obs['assay_ontology_term_id'].isin(smart_assay_list)
        adata.obs.loc[mask, 'donor_id'] += '-smartseq'

        if adata.raw:
            adata = ad.AnnData(sparse.csr_matrix(adata.raw.X), var=adata.raw.var, obs=adata.obs)

        female_adata = adata[:,adata.var.index.isin(female_ids)]
        male_adata = adata[:,adata.var.index.isin(male_ids)]

        genes_found = check_percent(female_adata,male_adata,female_ids,male_ids)
        if genes_found[0] == 0 or genes_found[1] == 0:
            return None
        fm_counts_dict = generate_fm_dict(female_ids,female_adata,male_ids,male_adata,adata)
        donor_sex_df,removed_donors = calculate_sex(fm_counts_dict)
        donor_sex_df = donor_sex_df[['donor_id','male_to_female','scRNAseq_sex']]
        donor_sex_df = donor_sex_df.merge(
            adata.obs[metadata_list].drop_duplicates(),
            on='donor_id',
            how='left'
        )
        sex_map = {
            'PATO:0000383':'female',
            'PATO:0000384':'male',
            'unknown':'unknown'
        }
        donor_sex_df['author_annotated_sex'] = donor_sex_df['sex_ontology_term_id'].map(sex_map)
        mask = donor_sex_df['assay_ontology_term_id'].isin(smart_assay_list)
        donor_sex_df.loc[mask, 'smart_seq'] = True
        donor_sex_df.drop(columns=['sex_ontology_term_id','assay_ontology_term_id'], inplace=True)
        donor_sex_df = donor_sex_df.drop_duplicates().sort_values('male_to_female')
        obs_to_keep = []
        ratio_order = []
        smart_seq_donors_rename = {}

        if 'smart_seq' in donor_sex_df.columns:
            donor_sex_df['smart_seq'] = (
                donor_sex_df['smart_seq']
                .where(donor_sex_df['smart_seq'].notna(), False)
                .astype('bool')
            )

        if donor_sex_df['smart_seq'].all() or not any(donor_sex_df['smart_seq']):
            adata.obs['donor_id'] = adata.obs['donor_id'].str.split('-smartseq').str[0]
            donor_sex_df['donor_id'] = donor_sex_df['donor_id'].str.split('-smartseq').str[0]
            valid_donors = donor_sex_df[~donor_sex_df['donor_id'].isin(removed_donors)]['donor_id']
            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(valid_donors)].index)
            ratio_order.append(
                (donor_sex_df['donor_id'] + ' ' +
                 donor_sex_df['author_annotated_sex'].astype('string')).to_list()
            )

        else:
            for d in pd.Series(donor_sex_df['donor_id'].str.split('-smartseq').str[0]).unique():
                if d not in removed_donors:
                    try:
                        smart_seq_sex = donor_sex_df.loc[
                            (donor_sex_df['donor_id'] == d + '-smartseq') &
                            (donor_sex_df['smart_seq'] == True)
                        ]['scRNAseq_sex'].unique()
                        nonsmart_seq_sex = donor_sex_df.loc[
                            (donor_sex_df['donor_id'] == d) &
                            (donor_sex_df['smart_seq'] == False)
                        ]['scRNAseq_sex'].unique()

                        if len(smart_seq_sex) > 0 and len(nonsmart_seq_sex) > 0:
                            if smart_seq_sex != nonsmart_seq_sex:
                                report(
                                    f'Smart-seq and non-Smart-seq scRNAseq_sex for donor ({d}) '
                                    'do not match - both will be included in plot.',
                                    'WARNING'
                                )
                                d_df = donor_sex_df[
                                    (donor_sex_df['donor_id'] == d) |
                                    (donor_sex_df['donor_id'] == d + '-smartseq')
                                ]
                                obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                                ratio_order.append(
                                    (d_df['donor_id'] + ' ' +
                                     d_df['author_annotated_sex'].astype('string')).to_list()
                                )

                            if smart_seq_sex == nonsmart_seq_sex:
                                report(
                                    f'Smart-seq and non-smart-seq scRNAseq_sex for donor ({d}) '
                                    'match, dropping Smart-seq from plot.'
                                )
                                d_df = donor_sex_df[donor_sex_df['donor_id'] == d]
                                obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                                ratio_order.append(
                                    (d_df['donor_id'].str.split('-smartseq').str[0] + ' ' +
                                     d_df['author_annotated_sex'].astype('string')).to_list()
                                )

                        elif len(smart_seq_sex) > 0 and len(nonsmart_seq_sex) == 0:
                            d_df = donor_sex_df[donor_sex_df['donor_id'] == d + '-smartseq']
                            smart_seq_donors_rename[f'{d}-smartseq'] = d
                            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                            ratio_order.append(
                                (d_df['donor_id'].str.split('-smartseq').str[0] + ' ' +
                                 d_df['author_annotated_sex'].astype('string')).to_list()
                            )


                        elif len(smart_seq_sex) == 0 and len(nonsmart_seq_sex) > 0:
                            d_df = donor_sex_df[donor_sex_df['donor_id'] == d]
                            obs_to_keep.append(adata.obs[adata.obs['donor_id'].isin(d_df['donor_id'])].index)
                            ratio_order.append(
                                (d_df['donor_id'] + ' ' +
                                d_df['author_annotated_sex'].astype('string')).to_list()
                            )

                    except Exception as e:
                        report(
                            f'Error: smart-seq and non-smart-seq sex for donor {d} were not calculated. Details: {e}',
                            'WARNING'
                        )
                        obs_to_keep, ratio_order = None, None  # Set to None to indicate failure

        flattened_obs_to_keep = [obs for sublist in obs_to_keep for obs in sublist]
        flattened_ratio_order = [ro for sublist1 in ratio_order for ro in sublist1]
        adata_sub = adata[flattened_obs_to_keep, : ].copy()
        adata_sub.obs['donor_id'] = adata_sub.obs['donor_id'].astype('category')
        adata.obs['donor_id'] = adata.obs['donor_id'].str.split('-smartseq').str[0]
        donor_sex_df['donor_id'] = donor_sex_df['donor_id'].str.split('-smartseq').str[0]
        adata_sub.obs['donor_id'] = adata_sub.obs['donor_id'].cat.rename_categories(smart_seq_donors_rename)
        adata_sub.obs['donor_sex'] = (
            adata_sub.obs.apply(
                lambda x: f"{x['donor_id']} {sex_map[x['sex_ontology_term_id']]}",
                axis=1
            ).astype('category')
        )
        adata_sub.var.rename(index=genes['female'], inplace=True)
        adata_sub.var.rename(index=genes['male'], inplace=True)
        f_symbs = [g for g in genes['female'].values() if g in adata_sub.var.index]
        m_symbs = [g for g in genes['male'].values() if g in adata_sub.var.index]
        dp = sc.pl.dotplot(
            adata_sub,
            {'female': f_symbs, 'male': m_symbs},
            'donor_sex',
            use_raw=False,
            categories_order=flattened_ratio_order,
            return_fig=True
        )

        if not donor_sex_df.empty:
            compare_donor_sex(donor_sex_df)

        if dp:
            dp.show()


def evaluate_var(adata):
    """
    Use single-cell-curation classes and fuctions and report warning/error for organism specific minimum number of gene features. Also, this function
    will look that var contains features from only a single organism.

    :param obj adata: AnnData that is being curated

    :return logging: Raises WARNING or ERROR if the number of genes in adata.var and adata.raw.var are fewer than the threshold of 40% or 60%, respectively,
    of 10x preselected biotype of genes. Will also raise ERROR if dataset contains more than a single organism; and raise a WARNING if the dataset has 50% or
    more genes filtered.
    """
    accepted_biotypes = [
        'protein_coding','protein_coding_LoF','lncRNA',
        'IG_C_gene','IG_D_gene','IG_J_gene','IG_LV_gene','IG_V_gene',
        'IG_V_pseudogene','IG_J_pseudogene','IG_C_pseudogene',
        'TR_C_gene','TR_D_gene','TR_J_gene','TR_V_gene',
        'TR_V_pseudogene','TR_J_pseudogene'
    ]

    organisms_with_descendants = [
        'NCBITaxon:9541','NCBITaxon:9544','NCBITaxon:10090','NCBITaxon:9986',
        'NCBITaxon:9598','NCBITaxon:10116','NCBITaxon:9823'
    ]

    # Check that this is single organism both in metadata and var index, exit function if multiple organisms or contains invalid var features
    var_organism_objs = list({gencode.get_organism_from_feature_id(id) for id in adata.var.index.to_list()})
    if None in var_organism_objs:
        report(
            'Some features in var.index are not valid gene IDs. index may be gene symbols or contain deprecated IDs',
            'ERROR'
        )
        report('To remove deprecated IDs, run...')
        report('adata = map_filter_gene_ids(adata)', 'code')
        return
    valid = True
    uns_organism = adata.uns['organism_ontology_term_id']
    var_organisms = [o.value for o in var_organism_objs]

    if 'NCBITaxon:2697049' in var_organisms:
        report('There are covid genes present in var')
        var_organisms.remove('NCBITaxon:2697049')
    if 'NCBITaxon:2697049' == uns_organism:
        report('Covid is not a supported uns.organism', 'ERROR')
        valid = False
    if len(var_organisms) > 1:
        report(f'Multiple organisms found in var index: {var_organisms}', 'ERROR')
        valid = False
    if valid:
        if var_organisms[0] in organisms_with_descendants:
            if utils.is_ontological_descendant_of(ONTOLOGY_PARSER,uns_organism,var_organisms[0]):
                report(f'Single organism found: {var_organisms}', 'GOOD')
            else:
                report(f'uns metadata contains non-descendant of var index organism: {var_organisms[0]}, {uns_organism}', 'ERROR')
                return
        elif uns_organism == var_organisms[0]:
            report(f'Single organism found: {var_organisms}', 'GOOD')
        else:
            report(f'Different organisms found between var index ({var_organisms[0]}) and uns metadata ({uns_organism})', 'ERROR')
            return
    else:
        return

    gene_count = len(adata.var)
    # unpaired ATAC have no gene count criteria
    if adata.obs['assay_ontology_term_id'].unique()[0] in ['EFO:0010891','EFO:0030007','EFO:0008925','EFO:0008904','EFO:0022045']:
        return
    elif [a for a in adata.obs['assay_ontology_term_id'].unique() if a in FLEX_ASSAYS]:
        count_type = 'Flex'
        warn_cut = 0.9
        err_cut = 0.7
        if var_organisms[0] == 'NCBITaxon:9606':
            flex_v2_count = 18132
            target_count = 18082
        elif var_organisms[0] == 'NCBITaxon:10090':
            flex_v2_count = 19070
            target_count = 19059
        else:
            report('Update required to support Flex data for non-human/mouse', 'ERROR')
            return

        if gene_count > flex_v2_count:
            report(f'{gene_count} genes present, expecting at most {flex_v2_count} for Flex V2', 'ERROR')
            return
    else:
        warn_cut = 0.6
        err_cut = 0.4
        if adata.uns['organism_ontology_term_id'] in ['NCBITaxon:10090', 'NCBITaxon:9606']:
            target_count = 35_000
            count_type = '10x biotype (CellRanger reference version)'
        else:
            # Check the number of genes threshold base on biotype per specific organism
            org_obj = [i for i in gencode.SupportedOrganisms if i.value==var_organisms[0]][0]
            gene_checker = gencode.GeneChecker(org_obj)
            target_count = len([i for i in gene_checker.gene_dict.keys() if gene_checker.gene_dict[i][2] in accepted_biotypes])
            count_type = '10x biotype'

    fraction = gene_count / target_count
    percent = fraction * 100
    if fraction < err_cut:
        report(
            f'{gene_count} genes present, compared against {target_count} {count_type} genes:'\
            f'{percent:.1f}% ({err_cut} threshold)\n'\
            'Data may not be eligible for submission without a less filtered gene set',
            'ERROR'
        )
    elif fraction < warn_cut:
        report(
            f'{gene_count} genes present, compared against {target_count} {count_type} genes:'\
            f'{percent:.1f}% ({warn_cut} threshold)\n'\
            'A less filtered gene set should be requested',
            'WARNING'
        )
    else:
        report(
            f'{gene_count} genes present, compared against {target_count} {count_type} genes: '\
            f'{percent:.1f}%',
            'GOOD'
        )

    # Check the number of filtered genes
    if 'feature_is_filtered' in adata.var.columns:
        if True in adata.var.feature_is_filtered.unique():
            num_filtered_genes = len(adata.var[adata.var.feature_is_filtered == True])
            frac_filtered = num_filtered_genes / gene_count * 100
            report(f'{num_filtered_genes} ({frac_filtered:.1f}%) genes are filtered from .X')
    else:
        report('feature_is_filtered not found in var', 'ERROR')
