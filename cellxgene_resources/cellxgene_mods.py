import anndata as ad
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import scanpy as sc
import subprocess
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
        non_integer = np.any(~np.equal(np.mod(adata.raw.X.data, 1), 0))
    else:
        non_integer = np.any(~np.equal(np.mod(adata.X.data, 1), 0))

    if non_integer == False:
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
    if adata.raw:
        hashes = [hash(r.tobytes()) for r in adata.raw.X.toarray()]
    else:
        hashes = [hash(r.tobytes()) for r in adata.X.toarray()]
    
    hash_df = adata.obs.copy()
    hash_df['hashes'] = hashes
    hash_df = hash_df[hash_df.duplicated(subset='hashes',keep=False) == True]
    hash_df.sort_values('hashes', inplace=True)
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

    sc.pl.spatial(adata, ax=axs[0], library_id=lib, color=cellpop_field, show=False)
    sc.pl.spatial(adata, ax=axs[1], library_id=lib)

    if 'fullres' in adata.uns['spatial'][lib]['images']:
        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
        )
        plt.subplots_adjust(wspace=wspace)
        sc.pl.spatial(adata, ax=axs[0], library_id=lib, img_key='fullres', scale_factor=1, color=cellpop_field, show=False)
        sc.pl.spatial(adata, ax=axs[1], library_id=lib, img_key='fullres', scale_factor=1)
    else:
        report('fullres image is highly recommended', 'WARNING')


def validate(file):
    validate_process = subprocess.run(['cellxgene-schema', 'validate', file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in validate_process.stdout.decode('utf-8').split('\n'):
        report(line)
    for line in validate_process.stderr.decode('utf-8').split('\n'):
        if line.endswith('is_valid=True'):
            report(line, 'GOOD')
        elif line.endswith('is_valid=False'):
            report(line, 'ERROR')
        else:
            prefix = line.split(':')[0]
            if prefix in ['ERROR','WARNING']:
                report(line.replace(f'{prefix}:',''), prefix)
            else:
                report(line)
