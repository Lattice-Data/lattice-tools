import anndata as ad
import json
import numpy as np
import os
import pandas as pd
import re
from scipy import sparse


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
    if adata.raw:
        report('raw min = ' + str(adata.raw.X.min()))
        report('raw max = ' + str(adata.raw.X.max()))
        non_integer = np.any(~np.equal(np.mod(adata.raw.X.data, 1), 0))
    else:
        non_integer = np.any(~np.equal(np.mod(adata.X.data, 1), 0))
    
    if non_integer == False:
        report('raw is all integers', 'GOOD')
    else:
        report('raw contains non-integer values', 'ERROR')
    
    report('X min = ' + str(adata.X.min()))
    report('X max = ' + str(adata.X.max()))
    
    for l in adata.layers:
        report(f'layers[{l}] min = ' + str(adata.layers[l].min()))
        report(f'layers[{l}] max = ' + str(adata.layers[l].max()))


def map_filter_gene_ids(adata):
    #map genes
    v44_gene_map = json.load(open('../gene_ID_mapping/gene_map_v44.json'))
    approved_file = 'ref_files/genes_approved.csv'
    approved = pd.read_csv(approved_file,dtype='str')

    my_gene_map = {k:v for k,v in v44_gene_map.items() if k in adata.var.index and v not in adata.var.index}
    adata.var.rename(index=my_gene_map, inplace=True)

    #filter out genes
    index_name = adata.var.index.name
    adata.var.reset_index(inplace=True)
    var_to_keep = adata.var[adata.var[index_name].isin(approved['feature_id'])].index #what if it's not called 'gene_ids'
    adata = adata[:, var_to_keep]
    adata.var.set_index(index_name, inplace=True)

    if adata.raw:
        raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs) #do we need to define obs?
        my_gene_map = {k:v for k,v in v44_gene_map.items() if k in raw_adata.var.index and v not in raw_adata.var.index}
        raw_adata.var.rename(index=my_gene_map, inplace=True)
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
