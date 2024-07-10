import anndata as ad
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import squidpy as sq
import subprocess
import sys
from scipy import sparse


portal_uns_fields = [
    'citation',
    'schema_reference',
    'schema_version'
]

curator_uns_fields = [
    'title'
]

portal_var_fields = [
    'feature_name',
    'feature_reference',
    'feature_biotype',
    'feature_length'
]

portal_obs_fields = [
    'assay',
    'cell_type',
    'development_stage',
    'disease',
    'self_reported_ethnicity',
    'organism',
    'sex',
    'tissue'
]
non_ontology_fields = ['donor_id','suspension_type','tissue_type','is_primary_data']
curator_obs_fields = [e + '_ontology_term_id' for e in portal_obs_fields] + non_ontology_fields
full_obs_standards = portal_obs_fields + curator_obs_fields


class CxG_API:
    scc_repo_loc = os.path.expanduser('~/GitClones/CZI/')
    sys.path.append(os.path.abspath(scc_repo_loc + 'single-cell-curation/notebooks/curation_api/python/'))


    from src.collection import create_collection,create_revision,get_collection,get_collections,update_collection
    from src.dataset import create_dataset,delete_dataset,get_dataset,get_datasets,upload_datafile_from_link,upload_local_datafile


    def config(env=None):
        from src.utils.config import set_api_access_config


        if env == 'dev':
            api_key_file_path = os.path.expanduser('~/Documents/keys/cxg-api-key-dev.txt')
            set_api_access_config(api_key_file_path, env='dev')
        elif env == 'staging':
            api_key_file_path = os.path.expanduser('~/Documents/keys/cxg-api-key-staging.txt')
            set_api_access_config(api_key_file_path, env='staging')
        else:
            api_key_file_path = os.path.expanduser('~/Documents/keys/cxg-api-key.txt')
            set_api_access_config(api_key_file_path)


def report(mess, level=None):
    colors = {
        'GOOD': '\033[32m', #green
        'WARNING': '\033[33m', #yellow
        'ERROR': '\033[31m' #red
    }
    if level:
        c = colors[level]
        print(f'\033[1m{c}{level}:{mess}\033[0m')
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


def get_adata_size(adata: ad.AnnData, show_stratified=True) -> int:
    """
    Adapted from the adata.__sizeof__() method. This version will also return sizes of
    adata.raw and the visium image arrays if they exist. Returns int of 
    size and can print out itemized display
    
    :param: show_stratified: print out attribute sizes if set to True
    """
    def get_size(X) -> int:
        if isinstance(X, (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix)):
            return X.data.nbytes + X.indptr.nbytes + X.indices.nbytes
        elif isinstance(X, np.ndarray):
            return X.nbytes
        else:
            return X.__sizeof__()

    def has_visium_uns_images(adata: ad.AnnData) -> bool:
        return (
            "spatial" in adata.uns and
            [k for k in adata.uns if "is_single" not in k]
        )

    size = 0
    attrs = list(["_X", "_obs", "_var"])
    attrs_raw = list(["X", "var"])
    attrs_multi = list(["_uns", "_obsm", "_varm", "varp", "_obsp", "_layers"])
    if adata.raw:
        adata_raw = adata.raw
        attrs.extend(attrs_raw)
        
    for attr in attrs + attrs_multi:
        if attr in attrs_multi:
            keys = getattr(adata, attr).keys()
            s = sum([get_size(getattr(adata, attr)[k]) for k in keys])
        elif attr in attrs_raw:
            s = get_size(getattr(adata_raw, attr))
        else:
            s = get_size(getattr(adata, attr))

        if s > 0 and show_stratified:
            if "_" not in attr:
                str_attr = ".raw." + attr + " " * (13 - len(attr))
            else:
                str_attr = attr.replace("_", ".") + " " * (18 - len(attr))
            print(f"Size of {str_attr}: {'%3.2f' % (s / (1024 ** 2))} MB")

        size += s

    if has_visium_uns_images(adata):
        library_id = [k for k in adata.uns["spatial"].keys() if "is_single" not in k][0]
        print("Visium image arrays:") if show_stratified else None
        for image_name, image_array in adata.uns["spatial"][library_id]["images"].items():
            s = get_size(image_array)
            if show_stratified:
                str_name = image_name + " " * (18 - len(image_name))
                print(f"Size of {str_name}: {'%3.2f' % (s / (1024 ** 2))} MB")
            size += s
            
    return size


def determine_sparsity(x):
    if isinstance(x, sparse.coo_matrix) or isinstance(x, sparse.csr_matrix) or isinstance(x, sparse.csc_matrix):
        sparsity = 1 - x.count_nonzero() / float(np.cumprod(x.shape)[-1])
    elif isinstance(x, np.ndarray):
        sparsity = 1 - np.count_nonzero(x) / float(np.cumprod(x.shape)[-1])
    else:
        report(f'matrix is of type {type(x)}, sparsity calculation has not been implemented')

    return round(sparsity, 3)


def evaluate_sparsity(adata):
    max_sparsity = 0.5

    sparsity = determine_sparsity(adata.X)
    report(f'X sparsity: {sparsity}')
    if sparsity > max_sparsity and type(adata.X) != sparse.csr_matrix:
        report('X should be converted to sparse', 'WARNING')
    
    if adata.raw:
        sparsity = determine_sparsity(adata.raw.X)
        report(f'raw.X sparsity: {sparsity}')
        if sparsity > max_sparsity and type(adata.raw.X) != sparse.csr_matrix:
            report('raw.X should be converted to sparse', 'WARNING')
    
    for l in adata.layers:
        sparsity = determine_sparsity(adata.layers[l])
        report(f'layers[{l}] sparsity: {sparsity}')
        if sparsity > max_sparsity and type(adata.layers[l]) != sparse.csr_matrix:
            report(f'layers[{l}] should be converted to sparse', 'WARNING')


def evaluate_data(adata):
    min_maxs = {}
    if adata.raw:
        raw_min = adata.raw.X.min()
        raw_max = adata.raw.X.max()
        report(f'raw min = {raw_min}')
        report(f'raw max = {raw_max}')
        min_maxs['raw'] = f'{raw_min}-{raw_max}'
        all_integers = np.all(np.round(adata.raw.X.data) == adata.raw.X.data)
    else:
        all_integers = np.all(np.round(adata.X.data) == adata.X.data)

    if all_integers:
        report('raw is all integers', 'GOOD')
    else:
        report('raw contains non-integer values', 'ERROR')

    X_min = adata.X.min()
    X_max = adata.X.max()
    report(f'X min = {X_min}')
    report(f'X max = {X_max}')
    min_maxs['X'] = f'{X_min}-{X_max}'

    for l in adata.layers:
        min = adata.layers[l].min()
        max = adata.layers[l].max()
        report(f'layers[{l}] min = {min}')
        report(f'layers[{l}] max = {max}')
        min_maxs[l] = f'{min}-{max}'

    poss_dups = [k for k,v in min_maxs.items() if list(min_maxs.values()).count(v) > 1]
    if poss_dups:
        report(f'possible redundant layers: {poss_dups}','WARNING')


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
    v44_gene_map = json.load(open('../gene_ID_mapping/gene_map_v44.json'))
    approved_file = 'ref_files/genes_approved.csv'
    approved = pd.read_csv(approved_file,dtype='str')

    my_gene_map = {k:v for k,v in v44_gene_map.items() if k in adata.var.index and v not in adata.var.index}
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
        my_gene_map = {k:v for k,v in v44_gene_map.items() if k in raw_adata.var.index and v not in raw_adata.var.index}
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


def barcode_compare(ref_df, obs_df):
    obs_df_split = obs_df.index.str.split('([ACTG]{16})')
    barcodes = pd.DataFrame([b for l in obs_df_split for b in l if re.match(r".*[ACTG]{16}.*", b)])
    if barcodes.empty:
        return pd.DataFrame({'summary':['no barcode'] * len(obs_df)})
    else:
        barcodes.rename(columns={0:'barcode'},inplace=True)
        barcodes.set_index('barcode', inplace=True)
        barcode_results = barcodes.merge(ref_df,on='barcode',how='left')
        barcode_results.fillna(0, inplace=True)
        barcode_results['summary'].replace(0, None, inplace=True)
        return barcode_results


def evaluate_10x_barcodes(prop, obs):

    results = []
    csv = 'ref_files/10X_barcode_table.csv.gz'
    ref_df = pd.read_csv(csv, sep=',', header=0, index_col='barcode')
    
    for a in obs[prop].unique():
        obs_df = obs[obs[prop] == a]
        r = barcode_compare(ref_df, obs_df)
        r_dict = {'3pv2_5pv1_5pv2': None, '3pv3': None, 'multiome': None,'multiple': None, 'None': None} | r['summary'].value_counts().to_dict()
        r_dict[prop] = a
        results.append(r_dict)
    
    return pd.DataFrame(results).set_index(prop).fillna(0).astype(int)


def evaluate_obs_schema(obs, labels=False):
    if labels:
        for o in portal_obs_fields + non_ontology_fields:
            if o not in obs.keys():
                report(f'{o} not in obs\n', 'ERROR')
            else:
                un = obs[o].unique()
                if un.dtype == 'category':
                    report(f'{o} {un.to_list()}\n')
                else:
                    report(f'{o} {un.tolist()}\n')
    else:
        for o in curator_obs_fields:
            if o not in obs.keys():
                report(f'{o} not in obs\n', 'ERROR')
            else:
                un = obs[o].unique()
                if un.dtype == 'category':
                    report(f'{o} {un.to_list()}\n')
                else:
                    report(f'{o} {un.tolist()}\n')
        for o in portal_obs_fields:
            if o in obs.keys():
                report(f'schema conflict - {o} in obs\n', 'ERROR')


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


def evaluate_dup_counts(adata):
    """
    Hash sparse matrix using np.ndarrays that represent sparse matrix data.
    First pass will hash all rows via slicing the data array and append to copy of obs df
    Second pass will hash only duplicate rows in obs copy via the indices array.
    This will keep only true duplicated matrix rows and not rows with an indicental same
    ordering of their data arrays
    """
    matrix = adata.raw.X if adata.raw else adata.X

    if isinstance(matrix, np.ndarray):
        print("Matrix not in sparse format, please convert before hashing")
        return

    data_array = matrix.data
    index_array = matrix.indices
    indptr_array = matrix.indptr

    start, end = 0, matrix.shape[0]
    hashes = []
    while start < end:
        val = hash(data_array[indptr_array[start]:indptr_array[start + 1]].tobytes())
        hashes.append(val)
        start += 1

    def index_hash(index):
        obs_loc = adata.obs.index.get_loc(index)
        val = hash(index_array[indptr_array[obs_loc]:indptr_array[obs_loc + 1]].tobytes())
        return val
    
    hash_df = adata.obs.copy()
    hash_df['data_array_hash'] = hashes
    hash_df = hash_df[hash_df.duplicated(subset='data_array_hash',keep=False) == True]
    hash_df.sort_values('data_array_hash', inplace=True)

    hash_df['index_array_hash'] = [index_hash(row) for row in hash_df.index.to_list()]
    hash_df = hash_df[hash_df.duplicated(subset=['data_array_hash', 'index_array_hash'], keep=False) == True]

    if not hash_df.empty:
        report('duplicated raw counts', 'ERROR')
        return hash_df
    report('no duplicated raw counts', 'GOOD')

    
def symbols_to_ids(symbols, var):
    
    ref_files = [
        'genes_ercc.csv',
        'genes_homo_sapiens.csv',
        'genes_mus_musculus.csv',
        'genes_sars_cov_2.csv'
    ]
    
    ref_dir = 'ref_files/'
    if not os.path.exists(ref_dir + 'genes_approved.csv'):
        ids = pd.DataFrame()
        for f in ref_files:
            df = pd.read_csv(f, names=['feature_id','symb','num','length'],dtype='str',index_col=False)
            ids = ids.append(df)
            os.remove(f)
        ids.to_csv(ref_dir + 'genes_approved.csv', index=False)
    
    approved = pd.read_csv(ref_dir + 'genes_approved.csv',dtype='str')
    
    ensg_list = []
    for s in symbols:
        if s in approved['symb'].tolist():
            ensg_id = approved.loc[approved['symb'] == s, 'feature_id'].iloc[0]
            if ensg_id in var.index:
                ensg_list.append(ensg_id)
                report(f'{ensg_id} -- {s}')
            else:
                s = s[0] + s[1:].lower()
                if s in approved['symb'].tolist():
                    ensg_id = approved.loc[approved['symb'] == s, 'feature_id'].iloc[0]
                    if ensg_id in var.index:
                        ensg_list.append(ensg_id)
                        report(f'{ensg_id} -- {s}')
                    else:
                        report(f'{s}/{ensg_id} not found in var')    
                else:
                    report(f'{s} not found in gene file')

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
        elif collection.get(k) != v and k not in should_differ_collection:
            if k == 'datasets':
                diff_props = set()
                pub_datasets = {d['dataset_id']: d for d in collection[k]}
                rev_datasets = {d['dataset_id']: d for d in revision[k]}
                comp = {}
                new = {}
                for ds_id,v in rev_datasets.items():
                    if ds_id not in pub_datasets.keys():
                        new[ds_id] = {}
                        for p in ['title','cell_count']:
                            new[ds_id][p] = v[p]
                        for p in ['assay','organism','tissue']:
                            new[ds_id][p] = [a['label'] for a in v[p]]
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
                                    comp[ds_id][prop + '_REV'] = rev_val
                                    comp[ds_id][prop + '_PUB'] = pub_val
            else:
                print('not same: ' + k)
                if k not in ['datasets','publisher_metadata']:
                    if k in ['links']:
                        diff_in_pub = [l for l in collection[k] if l not in v]
                        print('--- published: ', diff_in_pub)
                        diff_in_rev = [l for l in v if l not in collection[k]]
                        print('----- revised: ', diff_in_rev)
                    else:
                        print('--- published: ', str(collection[k]))
                        print('----- revised: ', v)


    comp_df = pd.DataFrame(comp).transpose()
    comp_df = comp_df.dropna(subset=[c for c in comp_df.columns if c != 'title'], how='all')
    if not comp_df.empty:
        print('\033[1mRevised Datasets\033[0m')
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
        display(pd.DataFrame(new).transpose())

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
        return male_female_df
    except Exception as e:
        print(e)


def evaluate_donors_sex(adata):
    genes_file = 'ref_files/sex_analysis_genes.json'
    genes = json.load(open(genes_file))
    female_ids = genes['female']
    male_ids = genes['male']
    metadata_list = ['donor_id', 'sex_ontology_term_id']
    smart_assay_list = ['EFO:0010184','EFO:0008931','EFO:0008930','EFO:0010022','EFO:0700016','EFO:0022488','EFO:0008442']
    adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
    adata.obs.loc[adata.obs['assay_ontology_term_id'].isin(smart_assay_list) == True, 'donor_id'] += '-smartseq'

    if adata.raw:
        female_adata = adata.raw[:,adata.raw.var.index.isin(female_ids)]
        male_adata = adata.raw[:,adata.raw.var.index.isin(male_ids)]
    else:
        female_adata = adata[:,adata.var.index.isin(female_ids)]
        male_adata = adata[:,adata.var.index.isin(male_ids)]

    genes_found = check_percent(female_adata,male_adata,female_ids,male_ids)

    if 0 in genes_found:
        print('Cannot calculate sex - male/female/both gene set(s) not found in dataset.')
        return None
    else:
        fm_counts_dict = generate_fm_dict(female_ids,female_adata,male_ids,male_adata,adata)
        donor_sex_df = calculate_sex(fm_counts_dict)
        donor_sex_df = donor_sex_df[['donor_id','male_to_female','scRNAseq_sex']]
        donor_sex_df = donor_sex_df.merge(adata.obs[metadata_list].drop_duplicates(), on='donor_id', how='left')
        sex_map = {
            'PATO:0000383':'female',
            'PATO:0000384':'male',
            'unknown':'unknown'
        }
        donor_sex_df['author_annotated_sex'] = donor_sex_df['sex_ontology_term_id'].map(sex_map)
        donor_sex_df.drop(columns='sex_ontology_term_id', inplace=True)
        return donor_sex_df
