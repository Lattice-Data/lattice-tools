import argparse
import anndata as ad
import boto3
import lattice
import os
import pandas as pd
import shutil
import sys
import scanpy as sc
import re
import subprocess
import numpy as np
import collections
import logging
import gc
from scipy import sparse
from datetime import datetime
import json
import numbers
import squidpy as sq
import flattener_mods as fm

# Creating empty list for warnings
warning_list = []

EPILOG = '''
Examples:

    python %(prog)s -m local -f LATDF119AAA

For more details:

    python %(prog)s --help
'''

class GlobVals:
	def __init__(self, mfinal_obj, mfinal_adata, cxg_adata, cxg_adata_raw, cxg_obs, cxg_uns, cxg_obsm):
		self.mfinal_obj = mfinal_obj
		self.mfinal_adata = mfinal_adata
		self.cxg_adata = cxg_adata
		self.cxg_adata_raw = cxg_adata_raw
		self.cxg_obs = cxg_obs
		self.cxg_uns = cxg_uns
		self.cxg_obsm = cxg_obsm

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--file', '-f',
                        help='Any identifier for the matrix of interest.')
    parser.add_argument('--mode', '-m',
                        help='The machine to run on.')
    args = parser.parse_args()
    if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit()
    return args



# Download file object from s3
def download_file(download_url, directory, accession=None):
	bucket_name = download_url.split('/')[2]
	file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
	s3client = boto3.client("s3")
	if accession:
		file_ext = download_url.split('.')[-1]
		file_name = accession + '.' + file_ext
	else:
		file_name = download_url.split('/')[-1]
	print(file_name + ' downloading')
	try:
		s3client.download_file(bucket_name, file_path, directory + '/' + file_name)
	except subprocess.CalledProcessError as e:
		logging.error('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
		sys.exit('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
	else:
		print(file_name + ' downloaded')


# Download entire directory contents from S3
def download_directory(download_url, directory):
	bucket_name = download_url.split('/')[2]
	spatial_folder = download_url.replace('s3://{}/'.format(bucket_name), "")
	s3client = boto3.client("s3")
	results = s3client.list_objects_v2(Bucket=bucket_name, Prefix=spatial_folder, Delimiter='/')
	os.mkdir(fm.MTX_DIR+"/spatial")
	for file in results.get('Contents'):
		if file.get('Size') == 0:
			continue
		file_path = file.get('Key')
		file_name = file_path.split('/')[-1]
		try:
			s3client.download_file(bucket_name, file_path, directory + '/spatial/' + file_name)
		except subprocess.CalledProcessError as e:
			logging.error('ERROR: Failed to find file S3://{}/{}'.format(bucket, file_path))
			sys.exit('ERROR: Failed to find file s3://{}/{}'.format(bucket, file_path))
		else:
			print(file_name + ' downloaded')


# Compile all reference annotations for var features into one pandas df
def compile_annotations(files):
	ids = pd.DataFrame()
	urls = 'https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/'
	for key in files:
		filename = fm.MTX_DIR + "/" + files[key] + ".gz"
		if os.path.exists(filename) == False:
			filename = urls + files[key] + '.gz'
		df = pd.read_csv(filename, names=['feature_id','symbol','start','stop'], dtype='str')
		ids = pd.concat([ids, df])
	return ids


# Recreated the final matrix ids, also checking to see if '-1' was removed from original cell identifier
def concatenate_cell_id(mxr_acc, raw_obs_names, glob):
	mfinal_cells = glob.mfinal_adata.obs.index.to_list()
	new_ids = []
	flag_removed = False
	cell_mapping_dct = {}
	cell_location = glob.mfinal_obj['cell_label_location']
	for mapping_dict in glob.mfinal_obj['cell_label_mappings']:
		cell_mapping_dct[mapping_dict['raw_matrix']] = mapping_dict['label']
	for final_id in mfinal_cells:
		if not re.search('[AGCT]+-1', final_id):
			flag_removed = True
	for id in raw_obs_names:
		if flag_removed:
			id = id.replace('-1', '')
		if cell_location == 'prefix':
			new_ids.append(cell_mapping_dct[mxr_acc]+id)
		elif cell_location == 'suffix':
			new_ids.append(id+cell_mapping_dct[mxr_acc])
	return(new_ids)

# Quality check final anndata created for cxg, sync up gene identifiers if necessary
def quality_check(glob):
	if glob.cxg_adata.obs.isnull().values.any():
		warning_list.append("WARNING: There is at least one 'NaN' value in the following cxg anndata obs columns: {}".format\
			(glob.cxg_adata.obs.columns[glob.cxg_adata.obs.isna().any()].tolist()))
	elif 'default_visualization' in glob.cxg_adata.uns:
		if adata.uns['default_visualization'] not in glob.cxg_adata.obs.values:
			logging.error('ERROR: The default_visualization field is not in the cxg anndata obs dataframe.')
			sys.exit("ERROR: The default_visualization field is not in the cxg anndata obs dataframe.")
	elif glob.mfinal_obj['X_normalized'] == True and glob.mfinal_obj['assays'] != ['snATAC-seq']:
		if len(glob.cxg_adata.var.index.tolist()) > len(glob.cxg_adata.raw.var.index.tolist()):
			logging.error('ERROR: There are more genes in normalized genes than in raw matrix.')
			sys.exit("ERROR: There are more genes in normalized genes than in raw matrix.")

# Return value to be stored in disease field based on list of diseases from donor and sample
def clean_list(lst, exp_disease):
	lst = lst.split(',')
	exp_disease_list = [i['term_id'] for i in exp_disease]
	disease_found = [i for i in exp_disease_list if i in lst]
	if disease_found:
		disease = disease_found[0]
		if len(disease_found) > 1:
			warning_list.append("WARNING: There is at least one sample with more than one experimental variable disease:\t{}".format(disease_found))
	else:
		disease = 'PATO:0000461'
	return disease

# Add temporary suffixes to var columns and then concatenate anndata objects in list together
# temp_anndata_list is created so that original anndata_list does not get suffixes
def concat_list(anndata_list, column, uns_merge):
	if column != 'none':
		suffix = 0
		for a in anndata_list:
			a.var.rename(columns={column:column + '-' + str(suffix)}, inplace=True)
			suffix += 1
	if uns_merge == True:
		concat_result = ad.concat(anndata_list, index_unique=None, join='outer', merge='unique',  uns_merge='first')
	else:
		concat_result = ad.concat(anndata_list, index_unique=None, join='outer', merge='unique')
	redundants = [i for i,c in collections.Counter(concat_result.obs.index.to_list()).items() if c>1]
	if len(redundants)>0:
		logging.error('ERROR: cell IDs are found in multiple raw matrix files.\t{}'.format(redundants))
		sys.exit('ERROR: cell IDs are found in multiple raw matrix files.\t{}'.format(redundants))
	drop_columns = [c for c in concat_result.obs.columns if c not in ['raw_matrix_accession', 'in_tissue', 'array_row', 'array_col']]
	if drop_columns:
		concat_result.obs.drop(columns=drop_columns, inplace=True)
	return concat_result

# Determine reported disease as unique of sample and donor diseases, removing unreported value
def report_diseases(mxr_df, exp_disease):
	mxr_df['reported_diseases'] = mxr_df['sample_diseases_term_name'] + ',' + mxr_df['donor_diseases_term_name']
	mxr_df['reported_diseases'] = mxr_df['reported_diseases'].apply(lambda x: '[{}]'.format(','.join([i for i in set(x.split(',')) if i!=fm.UNREPORTED_VALUE])))
	total_reported = mxr_df['reported_diseases'].unique()
	if len(total_reported) == 1:
		if total_reported[0] == '[]':
			mxr_df['reported_diseases'] = '[]'
	elif '[]' in total_reported:
		mxr_df['reported_diseases'].replace({'[]':'none'}, inplace=True)

	if exp_disease == fm.UNREPORTED_VALUE:
		mxr_df['disease_ontology_term_id'] = ['PATO:0000461'] * len(mxr_df.index)
	else:
		mxr_df['disease_ontology_term_id'] = mxr_df['sample_diseases_term_id'] + ',' + mxr_df['donor_diseases_term_id']
		mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].apply(clean_list, exp_disease=exp_disease)
		exp_disease_aslist = ['[{}]'.format(x['term_name']) for x in exp_disease]
		exp_disease_aslist.extend(['none', '[]'])
		if len([x for x in total_reported if x not in exp_disease_aslist])==0:
			mxr_df['reported_diseases'] = '[]'

# Demultiplex experimental metadata by finding demultiplexed suspension 
# Determine overlapping suspension, create library & demultiplexed suspension df
# get cell_metadata from that suspension, merge in library info
# merge with mxr_df on library
def demultiplex(lib_donor_df, library_susp, donor_susp, glob):
	susp_df = pd.DataFrame()
	lattice_donor = {}
	lattice_donor_col = []
	demult_susp_lst = []

	for donor_map in glob.mfinal_obj['donor_mappings']:
		lattice_donor[donor_map['label']] = donor_map['donor']
	for author_don in lib_donor_df['author_donor'].to_list():
		lattice_donor_col.append(lattice_donor[author_don])
	lib_donor_df['author_donor_@id'] = lattice_donor_col
	lib_donor_df['library_donor_@id'] = lib_donor_df['library_@id'] + "," + lib_donor_df['author_donor_@id']

	error = False
	error_list = []
	for lib_donor_unique in lib_donor_df['library_donor_@id'].to_list():
		demult_susp = ''
		lib_uniq = lib_donor_unique.split(',')[0]
		donor_uniq = lib_donor_unique.split(',')[1]

		for susp in donor_susp[donor_uniq]:
			if susp in library_susp[lib_uniq]:
				demult_susp = susp
		if demult_susp == '':
			logging.error('ERROR: Could not find suspension for demultiplexed donor: {}, {}, {}, {}'.format(donor_uniq, lib_uniq, donor_susp[donor_uniq], library_susp[lib_uniq]))
			print('ERROR: Could not find suspension for demultiplexed donor: {}, {}, {}, {}'.format(donor_uniq, lib_uniq, donor_susp[donor_uniq], library_susp[lib_uniq]))
			error = True
			error_list.append(donor_uniq)
		else:
			demult_susp_lst.append(demult_susp)
	if error:
		sys.exit("There are issues with finding common suspension for pooled library for the following donors:\t{}".format(error_list))
	lib_donor_df['suspension_@id'] = demult_susp_lst

	obj_type_subset = ['sample', 'suspension', 'donor']
	for susp in set(lib_donor_df['suspension_@id'].to_list()):
		values_to_add = {}
		susp_obj = lattice.get_object(susp, connection)
		relevant_objects = fm.gather_objects(susp_obj, glob.mfinal_obj, connection, start_type='suspension')
		for obj_type in obj_type_subset:
		 	objs = relevant_objects.get(obj_type, [])
		 	if len(objs) == 1:
		 		values_to_add = fm.gather_metdata(obj_type, fm.CELL_METADATA[obj_type], values_to_add, objs, connection)
		 	else:
		 		logging.error('ERROR: Could not find suspension for demultiplexed donor: {}'.format(obj_type))
		 		sys.exit('ERROR: Could not find suspension for demultiplexed donor: {}'.format(obj_type))
		row_to_add = pd.Series(values_to_add)
		susp_df = pd.concat([susp_df, row_to_add.to_frame().T], ignore_index=True)
	lib_donor_df = lib_donor_df.merge(susp_df, left_on='suspension_@id', right_on='suspension_@id', how='left')
	return(lib_donor_df)


# Ontologize sex from donor.sex enum
def get_sex_ontology(donor_df):
	term_lookup = {
		'female': 'PATO:0000383',
		'male': 'PATO:0000384'
	}
	sexes = donor_df['sex'].unique()
	for sex in sexes:
		if sex in term_lookup:
			donor_df.loc[donor_df['sex'] == sex, 'sex_ontology_term_id'] = term_lookup[sex]
		elif sex == 'unknown' or sex == 'mixed':
			donor_df.loc[donor_df['sex'] == sex, 'sex_ontology_term_id'] = 'unknown'
		else:
			logging.error('ERROR: Unexpected sex: {}'.format(sex))
			sys.exit("ERROR: Unexpected sex: {}".format(sex))


# Make sure cxg_adata and cxg_adata_raw have same number of features
# If not, add implied zeros to csr, and add corresponding 'feature_is_filtered'
def add_zero(glob):
	if glob.cxg_adata_raw.shape[1] > glob.cxg_adata.shape[1]:
		genes_add = [x for x in glob.cxg_adata_raw.var.index.to_list() if x not in glob.cxg_adata.var.index.to_list()]
		new_matrix = sparse.csr_matrix((glob.cxg_adata.X.data, glob.cxg_adata.X.indices, glob.cxg_adata.X.indptr), shape=glob.cxg_adata_raw.shape)
		all_genes = glob.cxg_adata.var.index.to_list()
		all_genes.extend(genes_add)
		new_var = pd.DataFrame(index=all_genes)
		new_var = pd.merge(new_var, glob.cxg_adata_raw.var, left_index=True, right_index=True, how='left')
		new_var['feature_is_filtered'] = False
		new_var.loc[genes_add, 'feature_is_filtered'] = True
		new_adata = ad.AnnData(X=new_matrix, obs=glob.cxg_adata.obs, var=new_var, uns=glob.cxg_adata.uns, obsm=glob.cxg_adata.obsm)
		if glob.cxg_adata.layers:
			for layer in glob.cxg_adata.layers:
				new_layer = sparse.csr_matrix((glob.cxg_adata.layers[layer].data, glob.cxg_adata.layers[layer].indices, 
					glob.cxg_adata.layers[layer].indptr), shape=glob.cxg_adata_raw.shape)
				new_adata.layers[layer] = new_layer
		new_adata = new_adata[:, glob.cxg_adata_raw.var.index.to_list()]
		new_adata.var = new_adata.var.merge(glob.cxg_adata.var, left_index=True, right_index=True, how='left')
		glob.cxg_adata = new_adata
	else:
		glob.cxg_adata.var['feature_is_filtered'] = False


# Use cxg_adata_raw var to map ensembl IDs and use that as index and drop redundants
# Make sure the indices are the same order for both anndata objects & clean up var metadata
# WILL NEED TO ADD NEW BIOTYPE FOR CITE-SEQ
def set_ensembl(redundant, glob):
	raw_cols = glob.cxg_adata_raw.var.columns.to_list()
	glob.cxg_adata_raw.var.drop(columns=[i for i in raw_cols if i!='gene_ids'], inplace=True)
	if glob.mfinal_obj['feature_keys'] == ['gene symbol']:
		if 'gene_ids' in glob.cxg_adata_raw.var.columns.to_list():
			# Drop unmapped genes from normalized matrix
			norm_index = set(glob.cxg_adata.var.index.to_list())
			raw_index = set(glob.cxg_adata_raw.var.index.to_list())
			drop_unmapped = list(norm_index.difference(raw_index))
			glob.cxg_adata = glob.cxg_adata[:, [i for i in glob.cxg_adata.var.index.to_list() if i not in drop_unmapped]]
			if len(drop_unmapped) > 0:
				warning_list.append('WARNING: {} unmapped gene_symbols were dropped. Full list available in logging file. Preview: {}'.format(len(drop_unmapped), drop_unmapped[:10]))
				warning_list.append('WARNING: Full list of the {} unmapped gene_ids dropped:\t{}'.format(len(drop_unmapped), drop_unmapped))
			glob.cxg_adata.var = pd.merge(glob.cxg_adata.var, glob.cxg_adata_raw.var, left_index=True, right_index=True, how='left', copy=True)
			glob.cxg_adata.var = glob.cxg_adata.var.set_index('gene_ids', drop=True)
			glob.cxg_adata_raw.var  = glob.cxg_adata_raw.var.set_index('gene_ids', drop=True)
			glob.cxg_adata.var.index.name = None
			glob.cxg_adata_raw.var.index.name = None

			# Drop redundant by Ensembl ID
			drop_redundant = list(set(redundant).intersection(set(glob.cxg_adata.var.index.to_list())))
			if len(drop_redundant) > 0:
				warning_list.append('WARNING: {} redundant gene_ids dropped:\t{}'.format(len(drop_redundant), drop_redundant))
			glob.cxg_adata = glob.cxg_adata[:, [i for i in glob.cxg_adata.var.index.to_list() if i not in redundant]]

		else:
			warning_list.append("WARNING: raw matrix does not have genes_ids column")
	elif glob.mfinal_obj['feature_keys'] == ['Ensembl gene ID']:
		glob.cxg_adata_raw.var_names_make_unique()
		glob.cxg_adata_raw.var  = glob.cxg_adata_raw.var.set_index('gene_ids', drop=True)
		glob.cxg_adata_raw.var.index.name = None
		unique_to_norm =  set(glob.cxg_adata.var.index.to_list()).difference(set(glob.cxg_adata_raw.var.index.to_list()))
		if len(unique_to_norm) > 0:
			warning_list.append("WARNING: normalized matrix contains {} Ensembl IDs not in raw".format(unique_to_norm))

# Filters the Ensembl IDs based on the compiled list of approved IDs
def filter_ensembl(adata, compiled_annot):
	# Using map file to map old ensembl_ids to new ensembl_ids before filtering
	v44_gene_map = json.load(open('../gene_ID_mapping/gene_map_v44.json'))
	new_gene_map = {k:v for k,v in v44_gene_map.items() if k in adata.var.index and v not in adata.var.index}
	adata.var.rename(index=new_gene_map, inplace=True)
	var_in_approved = adata.var.index[adata.var.index.isin(compiled_annot['feature_id'])]
	adata = adata[:, var_in_approved]
	return adata


# Clean up columns and column names for var, gene name must be gene_id
# WILL NEED TO REVISIT WHEN THERE IS MORE THAN ONE VALUE FOR GENE ID
def clean_var(glob):
	gene_pd = glob.cxg_adata_raw.var[[i for i in glob.cxg_adata_raw.var.columns.values.tolist() if 'gene_ids' in i]]
	gene_pd = gene_pd.replace('nan', np.nan)
	gene_pd = gene_pd.stack().groupby(level=0).apply(lambda x: x.unique()[0]).to_frame(name='gene_ids')
	glob.cxg_adata_raw.var.drop(columns=glob.cxg_adata_raw.var.columns.tolist(), inplace=True)
	glob.cxg_adata_raw.var = glob.cxg_adata_raw.var.merge(gene_pd, left_index=True, right_index=True, how='left')


# Reconcile genes if raw matrices annotated to multiple version by merging raw based on Ensembl ID
# Return raw matrix merged on Ensembl ID and list of gene to drop from normalized matrices
def reconcile_genes(cxg_adata_lst, glob):
	mfinal_adata_genes = glob.mfinal_adata.var.index.to_list()
	redundant = []
	stats = {}

	# Join raw matrices on ensembl, gene symbols stored as metadata
	for cxg_adata in cxg_adata_lst:
		cxg_adata.var['gene_symbols'] = cxg_adata.var.index
		cxg_adata.var = cxg_adata.var.set_index('gene_ids', drop=True)
	cxg_adata_raw_ensembl = concat_list(cxg_adata_lst, 'gene_symbols', False)

	# Join raw matrices on gene symbol, ensembl stored as metadata. Add suffix to make unique, using '.' as to match R default
	count = 0
	for cxg_adata in cxg_adata_lst:
		cxg_adata.var['gene_ids'] = cxg_adata.var.index
		cxg_adata.var =  cxg_adata.var.set_index('gene_symbols' + '-' + str(count), drop=True)
		count += 1
		cxg_adata.var_names_make_unique(join = '.')
	cxg_adata_raw_symbol = concat_list(cxg_adata_lst, 'gene_ids', False)

	# Go through adata indexed on symbol to see which have > 1 Ensembl ID
	gene_pd_symbol = cxg_adata_raw_symbol.var[[i for i in cxg_adata_raw_symbol.var.columns.values.tolist() if 'gene_ids' in i]]
	multiple_ensembl_df = gene_pd_symbol[gene_pd_symbol.stack().groupby(level=0).apply(lambda x: len([i for i in x.unique() if str(i)!='nan'])>1)==True]
	stats['multiple_ensembl'] = list(set(multiple_ensembl_df.index.to_list()))

	# Go through adata indexed on ensembl to see which have > 1 symbol
	gene_pd_ensembl = cxg_adata_raw_ensembl.var[[i for i in cxg_adata_raw_ensembl.var.columns.values.tolist() if 'gene_symbols' in i]]
	multiple_symbols_df = gene_pd_ensembl[gene_pd_ensembl.stack().groupby(level=0).apply(lambda x: len([i for i in x.unique() if str(i)!='nan'])>1)==True]
	stats['multiple_symbols'] = list(set(multiple_symbols_df.stack().groupby(level=0).apply(lambda x: x.unique().tolist()).sum()))

	# Log redundant gene Ensembl IDs from normalized matrix within a single version
	for col in gene_pd_ensembl.columns:
		for gene in [i for i, c in collections.Counter(gene_pd_ensembl[col].dropna().to_list()).items() if c > 1]:
			if True in gene_pd_ensembl[gene_pd_ensembl[col]==gene].index.str.endswith("PAR_Y"):
				redundant.extend([i for i in gene_pd_ensembl[gene_pd_ensembl[col]==gene].index.to_list() if i.endswith('PAR_Y')])
			else:
				redundant.extend(gene_pd_ensembl[gene_pd_ensembl[col] == gene].index.to_list())
	redundant = list(set(redundant))

	# Clean up raw.var in outer join on ensembl and switch to gene symbol for index. Run var_names_make_unique and remove redundants after mapping of Ensembl
	cxg_adata_raw_ensembl.var['gene_ids'] = cxg_adata_raw_ensembl.var.index
	cxg_adata_raw_ensembl.var['gene_symbols'] = gene_pd_ensembl.stack().groupby(level=0).apply(lambda x: x.unique()[0]).to_frame(name='gene_symbols')
	cxg_adata_raw_ensembl.var = cxg_adata_raw_ensembl.var.set_index('gene_symbols', drop=True)
	cxg_adata_raw_ensembl.var.index.name = None
	cxg_adata_raw_ensembl.var_names_make_unique(join='.')

	all_remove = list(set(stats['multiple_ensembl'] + stats['multiple_symbols']))

	for key in stats:
		stats[key] = set(stats[key])
		overlap_norm = set(mfinal_adata_genes).intersection(stats[key])
		warning_list.append("WARNING: Full list of {}\t{}\t{}\t{}".format(key, len(stats[key]), len(overlap_norm), overlap_norm))
		warning_list.append("WARNING: Genes with {}. {} in raw. {} in normalized. Full list available in logging file. Preview of normalized: {}".format(key, 
			len(stats[key]), len(overlap_norm), list(overlap_norm)[:10]))

	return cxg_adata_raw_ensembl, redundant, all_remove
	


# filename will be collectionuuid_datasetuuid_accession_version.h5ad, collectionuuid_accession_version.h5ad, or accession_version.h5ad
# depending on what information is available for the dataset
def get_results_filename(glob):
	results_file = None
	collection_id = None
	dataset = glob.mfinal_obj.get('dataset', [])
	obj_type, filter_url = lattice.parse_ids([dataset])
	dataset_objs = lattice.get_report(obj_type, filter_url, ['cellxgene_urls'], connection)
	if dataset_objs[0].get('cellxgene_urls', []):
		collection_id = dataset_objs[0].get('cellxgene_urls', [])[0]
		collection_id = collection_id.replace("https://cellxgene.cziscience.com/collections/", "")
	if glob.mfinal_obj.get('cellxgene_uuid', []) and collection_id:
		results_file = '{}_{}_{}_v{}.h5ad'.format(collection_id, glob.mfinal_obj['cellxgene_uuid'], glob.mfinal_obj['accession'], fm.SCHEMA_VERSION)
	elif collection_id:
		results_file = '{}_{}_v{}.h5ad'.format(collection_id, glob.mfinal_obj['accession'], fm.SCHEMA_VERSION)
	else:
		results_file = '{}_v{}.h5ad'.format(glob.mfinal_obj['accession'], fm.SCHEMA_VERSION)
	return results_file


# Add antibody metadata to var and raw.var
def map_antibody(glob):
	antibody_meta = pd.DataFrame()
	for anti_mapping in glob.mfinal_obj.get('antibody_mappings'):
		values_to_add = {}
		antibody = anti_mapping.get('antibody')
		values_to_add = fm.gather_metdata('antibody', fm.ANTIBODY_METADATA['antibody'], values_to_add, [antibody], connection)
		values_to_add['host_organism'] = re.sub(r'/organisms/(.*)/', r'\1', values_to_add['host_organism'])
		if not antibody.get('control'):
			target = None
			if len(antibody.get('targets')) > 1:
				for t in antibody.get('targets'):
					name = t.get('organism').get('scientific_name')
					if name == glob.cxg_adata.obs['organism'].unique()[0]:
						target = [t]
			else:
				target = antibody.get('targets')
			values_to_add = fm.gather_metdata('target', fm.ANTIBODY_METADATA['target'], values_to_add, target, connection)
			values_to_add['feature_name'] = values_to_add['target_label']
		else:
			values_to_add['feature_name'] = '{} {} (control)'.format(values_to_add['host_organism'], values_to_add['isotype'])
			for val in fm.ANTIBODY_METADATA['target']:
				latkey = ('target_' + val).replace('.', '_')
				key = prop_map.get(latkey, latkey)
				values_to_add[key] = 'na'
		row_to_add = pd.DataFrame(values_to_add, index=[anti_mapping.get('label')])
		antibody_meta = pd.concat([antibody_meta, row_to_add])
	glob.cxg_adata.var = pd.merge(glob.cxg_adata.var, antibody_meta, left_index=True, right_index=True, how='left')
	glob.cxg_adata_raw.var = pd.merge(glob.cxg_adata_raw.var, antibody_meta, left_index=True, right_index=True, how='left')
	
	glob.cxg_adata.var['author_index'] = glob.cxg_adata.var.index
	glob.cxg_adata_raw.var['author_index'] = glob.cxg_adata.var.index
	glob.cxg_adata_raw.var.drop(columns=['genome'], inplace=True)
	glob.cxg_adata.var['feature_biotype'] = 'antibody-derived tags'
	glob.cxg_adata.var.set_index('feature_name', inplace=True, drop=False)
	glob.cxg_adata_raw.var.set_index('feature_name', inplace=True, drop=False)
	glob.cxg_adata.var_names_make_unique(join='-')
	glob.cxg_adata_raw.var_names_make_unique(join='-')


# Final touches for obs columns, modifying any Lattice fields to fit cxg schema
def clean_obs(glob):
	# For columns in mfinal_obj that contain continuous cell metrics, they are transferred to cxg_obs as float datatype
	# WILL NEED TO REVISIT IF FINAL MATRIX CONTAINS MULTIPLE LAYERS THAT WE ARE WRANGLING
	for author_col in glob.mfinal_obj.get('author_columns', []):
		if author_col in glob.mfinal_adata.obs.columns.to_list():
			glob.cxg_obs = pd.merge(glob.cxg_obs, glob.mfinal_adata.obs[[author_col]], left_index=True, right_index=True, how='left')
		else:
			warning_list.append("WARNING: author_column not in final matrix: {}".format(author_col))
	# if obs category suspension_type does not exist in dataset, create column and fill values with na (for spatial assay)
	if 'suspension_type' not in glob.cxg_obs.columns:
		glob.cxg_obs.insert(len(glob.cxg_obs.columns), 'suspension_type', 'na')
	elif glob.cxg_obs['suspension_type'].isnull().values.any():
		glob.cxg_obs['suspension_type'].fillna(value='na', inplace=True)
	
	add_units = {'tissue_section_thickness': 'tissue_section_thickness_units',
				'suspension_dissociation_time': 'suspension_dissociation_time_units',
				'library_starting_quantity': 'library_starting_quantity_units'}
	for field in add_units.keys():
		if field in glob.cxg_obs.columns:
			glob.cxg_obs[field] = glob.cxg_obs[field].astype(str) + " " + glob.cxg_obs[add_units[field]].astype(str)
			glob.cxg_obs.drop(columns=add_units[field], inplace=True)
			glob.cxg_obs[field].replace({'unknown unknown':'unknown'}, inplace=True)

	make_numeric = ['suspension_percent_cell_viability','donor_BMI_at_collection']
	for field in make_numeric:
		if field in glob.cxg_obs.columns:
			if True in glob.cxg_obs[field].str.contains('[<>-]|'+fm.UNREPORTED_VALUE+'|'+'pooled', regex=True).to_list():
				if True not in glob.cxg_obs[field].str.contains('[<>-]', regex=True).to_list():
					glob.cxg_obs[field].replace({'unknown':np.nan}, inplace=True) 
					glob.cxg_obs[field][np.where(glob.cxg_obs[field].str.contains('pooled') == True)[0].tolist()] = np.nan
					glob.cxg_obs[field]  = glob.cxg_obs[field].astype('float')
			else: 
				glob.cxg_obs[field]  = glob.cxg_obs[field].astype('float')

	change_unreported = ['suspension_enriched_cell_types', 'suspension_depleted_cell_types', 'suspension_enrichment_factors', 
	'suspension_depletion_factors', 'disease_state', 'cell_state']
	for field in change_unreported:
		if field in glob.cxg_obs.columns.to_list():
			glob.cxg_obs[field].replace({fm.UNREPORTED_VALUE: 'na'}, inplace=True)
	valid_tissue_types = ['tissue', 'organoid', 'cellculture']
	glob.cxg_obs['tissue_type'] = glob.cxg_obs['tissue_type'].str.lower()
	for i in glob.cxg_obs['tissue_type'].unique().tolist():
		if i == 'cellculture':
			glob.cxg_obs['tissue_type'].replace({'cellculture':'cell culture'}, inplace=True)
		if i not in valid_tissue_types:
			logging.error('ERROR: not a valid tissue type:\t{}'.format(i))
			print('ERROR: not a valid tissue type:\t{}'.format(i))
	if glob.mfinal_obj['is_primary_data'] == 'mixed':
		primary_portion = glob.mfinal_obj.get('primary_portion')
		glob.cxg_obs['is_primary_data'] = False
		glob.cxg_obs.loc[glob.cxg_obs[primary_portion.get('obs_field')].isin(primary_portion.get('values')), 'is_primary_data'] = True

	glob.cxg_obs[[i for i in glob.cxg_obs.columns.tolist() if i.startswith('family_history_')]] = \
		glob.cxg_obs[[i for i in glob.cxg_obs.columns.tolist() if i.startswith('family_history_')]].fillna(value='unknown')
	
	# map gencode to ensembl version for HCA tier 1
	if 'gene_annotation_version' in glob.cxg_obs.columns:
		glob.cxg_obs['gene_annotation_version'] = glob.cxg_obs['gene_annotation_version'].map(fm.GENCODE_MAP)


# Drop any intermediate or optional fields that are all empty
def drop_cols(celltype_col, glob):
	optional_columns = ['donor_BMI_at_collection', 'donor_family_medical_history', 'reported_diseases', 'donor_times_pregnant', 'sample_preservation_method',\
			'sample_treatment_summary', 'suspension_uuid', 'tissue_section_thickness', 'tissue_section_thickness_units', 'cell_state', 'disease_state',\
			'suspension_enriched_cell_types', 'suspension_enrichment_factors', 'suspension_depletion_factors', 'tyrer_cuzick_lifetime_risk',\
			'donor_living_at_sample_collection', 'donor_menopausal_status', 'donor_smoking_status', 'sample_derivation_process', 'suspension_dissociation_reagent',\
			'suspension_dissociation_time', 'suspension_depleted_cell_types', 'suspension_derivation_process', 'suspension_percent_cell_viability',\
			'library_starting_quantity', 'library_starting_quantity_units', 'tissue_handling_interval', 'suspension_dissociation_time_units', 'alignment_software',\
			'gene_annotation_version', 'reference_genome', 'sequencing_platform', 'sample_source', 'donor_cause_of_death', 'growth_medium', 'genetic_modifications',
			'menstrual_phase_at_collection']
	
	if 'sequencing_platform' in glob.cxg_obs.columns:
		if glob.cxg_obs['sequencing_platform'].isnull().values.any():
			glob.cxg_obs['sequencing_platform'].fillna(fm.UNREPORTED_VALUE, inplace=True)
	for col in optional_columns:
		if col in glob.cxg_obs.columns.to_list():
			col_content = glob.cxg_obs[col].unique()
			if len(col_content) == 1:
				if col_content[0] == fm.UNREPORTED_VALUE or col_content[0] == '[' + fm.UNREPORTED_VALUE + ']' or col_content[0] == '[]':
					glob.cxg_obs.drop(columns=col, inplace=True)

	if len(glob.cxg_obs['donor_age_redundancy'].unique()) == 1:
		if glob.cxg_obs['donor_age_redundancy'].unique():
			glob.cxg_obs.drop(columns='donor_age', inplace=True)
	columns_to_drop = ['raw_matrix_accession', celltype_col, 'sample_diseases_term_id', 'sample_diseases_term_name', 'sample_biosample_ontology_organ_slims',\
			'donor_diseases_term_id', 'donor_diseases_term_name', 'batch', 'library_@id_x', 'library_@id_y', 'author_donor_x', 'author_donor_y',\
			'library_authordonor', 'author_donor_@id', 'library_donor_@id', 'suspension_@id', 'library_@id', 'sex', 'sample_biosample_ontology_cell_slims',\
			'sample_summary_development_ontology_at_collection_development_slims', 'donor_age_redundancy']
	for column_drop in columns_to_drop: 
		if column_drop in glob.cxg_obs.columns.to_list():
			glob.cxg_obs.drop(columns=column_drop, inplace=True)


# Add ontology term names to CXG standardized columns
def add_labels(glob):
	obj = lattice.get_report('OntologyTerm', '', ['term_id', 'term_name'], connection)
	ontology_df = pd.DataFrame(obj)
	id_cols = ['assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id',\
			'tissue_ontology_term_id', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id']
	for col in id_cols:
		name_col = col.replace("_ontology_term_id", "")
		glob.cxg_adata.obs[name_col] = glob.cxg_adata.obs[col]
		for term_id in glob.cxg_adata.obs[name_col].unique():
			if term_id == 'PATO:0000461':
				term_name = 'normal'
			elif term_id == 'PATO:0000384':
				term_name = 'male'
			elif term_id == 'PATO:0000383':
				term_name = 'female'
			elif term_id == 'NCBITaxon:9606':
				term_name = 'Homo sapiens'
			elif term_id == 'unknown':
				term_name = 'unknown'
			elif len(ontology_df.loc[ontology_df['term_id']==term_id, 'term_name'].unique() == 1):
				term_name = ontology_df.loc[ontology_df['term_id']==term_id, 'term_name'].unique()[0]
			else:
				logging.error('ERROR: Found more than single ontology term name for id: {}\t{}'.format(term_id, ontology_df.loc[ontology_df['term_id']==term_id, 'term_name'].unique()))
				sys.exit("ERROR: Found more than single ontology term name for id: {}\t{}".format(term_id, ontology_df.loc[ontology_df['term_id']==term_id, 'term_name'].unique()))
			glob.cxg_adata.obs[name_col].replace(term_id, term_name, inplace=True)


# Look at matrix and only convert to sparse if density is less than 0.5
def check_matrix(m):
	if not sparse.issparse(m):
		if (np.count_nonzero(m)/np.prod(m.shape)) <= 0.5:
			m = sparse.csr_matrix(m)
	elif m.getformat()=='csc':
		m = sparse.csr_matrix(m)
	return m


# Add background spots to raw adata, getting obs metadata from tissue_positions_list.csv
def add_raw_background_spots(adata, glob):
	adata.var_names_make_unique(join='.')
	all_barcodes = pd.read_csv(fm.MTX_DIR+"/spatial/tissue_positions_list.csv", index_col=0, header=None)
	all_barcodes.columns = ['in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
	for c in ['in_tissue', 'array_row','array_col']:
		if all_barcodes[c].dtype != int:
			all_barcodes[c] = all_barcodes[c].astype('int64')
	missing_barcodes = [i for i in all_barcodes.index.to_list() if i not in list(adata.obs.index)]
	empty_matrix = sparse.csr_matrix((len(missing_barcodes), adata.var.shape[0]))
	missing_adata = ad.AnnData(empty_matrix, var=adata.var, obs=pd.DataFrame(index=missing_barcodes))
	comb_adata = ad.concat([adata, missing_adata], uns_merge='first', merge='first')
	comb_adata = comb_adata[all_barcodes.index.to_list()]
	comb_adata.obs = pd.merge(comb_adata.obs, all_barcodes[['in_tissue','array_row','array_col']], left_index=True, right_index=True, how='left')
	comb_adata.obsm['spatial'] = all_barcodes.loc[:,['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()
	return comb_adata


# Add missing obs to mfinal_adata to match raw and also modify cxg_obsm to include missing obs
def add_background_spots(glob):
	if glob.mfinal_adata.obs.shape[0] < 4992:
		missing_barcodes = [i for i in glob.cxg_adata_raw.obs.index.to_list() if i not in list(glob.mfinal_adata.obs.index)]
		empty_matrix = sparse.csr_matrix((len(missing_barcodes), glob.mfinal_adata.var.shape[0]))
		missing_adata = ad.AnnData(empty_matrix, var=glob.mfinal_adata.var, obs=pd.DataFrame(index=missing_barcodes))
		comb_adata = ad.concat([glob.mfinal_adata, missing_adata], uns_merge='first', merge='first')
		comb_adata = comb_adata[glob.cxg_adata_raw.obs.index.to_list(), :]

		for embed in glob.mfinal_adata.obsm.keys():
			new_array = np.empty((comb_adata.obs.shape[0], glob.mfinal_adata.obsm[embed].shape[1]))
			new_array[:] = np.nan
			for orig_row in range(glob.mfinal_adata.obs.shape[0]):
				row = comb_adata.obs.index.get_loc(glob.mfinal_adata.obs.iloc[orig_row].name)
				new_array[row] = glob.mfinal_adata.obsm[embed][orig_row, :]
			comb_adata.obsm[embed] = new_array

		glob.mfinal_adata = comb_adata
		glob.cxg_obsm = glob.mfinal_adata.obsm
		glob.cxg_obs['cell_type_ontology_term_id'] = glob.cxg_obs['cell_type_ontology_term_id'].fillna('unknown')
	glob.cxg_obsm['spatial'] = glob.cxg_adata_raw.obsm['spatial']


def main(mfinal_id):

	glob = GlobVals(None, None, None, None, None, None, None)
	glob.mfinal_obj = lattice.get_object(mfinal_id, connection)

	logging.basicConfig(filename="{}_outfile_flattener.log".format(mfinal_id), filemode='w', level=logging.INFO)
	# Adding date and time to top of logging file
	time_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
	logging.info("Date and time of flattener run: " + time_date)
	# Suppressing specific warnings from anndata
	logging.captureWarnings(True)

	# confirm that the identifier you've provided corresponds to a ProcessedMatrixFile
	mfinal_type = glob.mfinal_obj['@type'][0]
	summary_assay = ''
	if mfinal_type != 'ProcessedMatrixFile':
		logging.error('ERROR: {} is not a ProcessedMatrixFile, but a {}'.format(mfinal_id, mfinal_type))
		sys.exit('ERROR: {} is not a ProcessedMatrixFile, but a {}'.format(mfinal_id, mfinal_type))

	if glob.mfinal_obj['output_types'] == ['gene quantifications']:
		if glob.mfinal_obj['assays'] == ['snATAC-seq']:
			summary_assay = 'ATAC'
		else:
			summary_assay = 'RNA'
	elif glob.mfinal_obj['output_types'] == ['antibody capture quantifications']:
		summary_assay = 'CITE'
	else:
		logging.error('ERROR: Unexpected assay types to generate cxg h5ad: {} {}'.format(glob.mfinal_obj['assays'], glob.mfinal_obj['output_types']))
		sys.exit("ERROR: Unexpected assay types to generate cxg h5ad: {} {}".format(glob.mfinal_obj['assays'], glob.mfinal_obj['output_types']))


	# Dataframe that contains experimental metadata keyed off of raw matrix
	df = pd.DataFrame()

	results = {}
	
	# Checking for presence of matrix_files, and creating if not present
	if os.path.exists(fm.MTX_DIR) == False:
		os.mkdir(fm.MTX_DIR)
		
	# Checking for presence of h5ad, and downloading if not present
	if os.path.exists(fm.MTX_DIR + '/' + glob.mfinal_obj['accession'] + '.h5ad'):
		print(glob.mfinal_obj['accession'] + '.h5ad' + ' was found locally')
	else:
		download_file(glob.mfinal_obj.get('s3_uri'), fm.MTX_DIR, glob.mfinal_obj.get('accession'))

	# Get list of unique final cell identifiers
	file_url = glob.mfinal_obj['s3_uri']
	file_ext = file_url.split('.')[-1]
	mfinal_local_path = '{}/{}.{}'.format(fm.MTX_DIR, glob.mfinal_obj['accession'], file_ext)
	glob.mfinal_adata = sc.read_h5ad(mfinal_local_path)
	mfinal_cell_identifiers = glob.mfinal_adata.obs.index.to_list()

	cxg_adata_lst = []
	redundant = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = fm.gather_rawmatrices(glob.mfinal_obj['derived_from'], connection)
	donor_susp = {}
	library_susp = {}
	mapping_error = False
	error_info = {}

	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		relevant_objects = fm.gather_objects(mxr, glob.mfinal_obj, connection)
		values_to_add = {}

		# Get raw matrix metadata
		values_to_add = fm.gather_metdata('raw_matrix', fm.CELL_METADATA['raw_matrix'], values_to_add, [mxr], connection)

		# If there is a demultiplexed_donor_column, assume it is a demuxlet experiment and demultiplex df metadata
		# Gather library, suspension, and donor associations while iterating through relevant objects
		# Cannot handle multiple pooling events, so will sys.exit
		if 'demultiplexed_donor_column' in glob.mfinal_obj:
			lib_obj = relevant_objects.get('library', [])
			values_to_add = fm.gather_metdata('library', fm.CELL_METADATA['library'], values_to_add, lib_obj, connection)
			for i in range(len(lib_obj)):
				for single_lib_susp in lib_obj[i]['derived_from']:
					if lib_obj[i]['@id'] not in library_susp:
						library_susp[lib_obj[i]['@id']] = [single_lib_susp['@id']]
					else:
						library_susp[lib_obj[i]['@id']].append(single_lib_susp['@id'])
			susp_obj = relevant_objects.get('suspension', [])
			if isinstance(susp_obj, list):
				for single_susp in susp_obj:
					if len(single_susp['donors']) > 1:
						logging.error('ERROR: Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
						sys.exit('ERROR: Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
					if single_susp['donors'][0] not in donor_susp:
						donor_susp[single_susp['donors'][0]] = [single_susp['@id']]
					else:
						donor_susp[single_susp['donors'][0]].append(single_susp['@id'])
			else:
				if len(single_susp['donors']) > 1:
					logging.error('ERROR: Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
					sys.exit('ERROR: Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
				if single_susp['donors'][0] not in donor_susp:
					donor_susp[susp_obj['donors'][0]] = [susp_obj['@id']]
				else:
					donor_susp[susp_obj['donors'][0]].append(susp_obj['@id'])

		# Gather metdata without demultiplexing
		else:
			for obj_type in fm.CELL_METADATA.keys():
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					values_to_add = fm.gather_metdata(obj_type, fm.CELL_METADATA[obj_type], values_to_add, objs, connection)
				elif len(objs) > 1:
					# Check to make sure it is not multimodal before determining that it is pooled
					if obj_type == 'library':
						value = list()
						for obj in objs:
							v = fm.get_value(obj, 'protocol.assay_ontology.term_id')
							value.append(v)
						if set(value) == {'EFO:0030059'}:
							if glob.mfinal_obj.get('assays') == ['snATAC-seq']:
								single_obj = [o for o in objs if o.get('assay')=='snATAC-seq']
							else:
								single_obj = [o for o in objs if o.get('assay')=='snRNA-seq']
							values_to_add = fm.gather_metdata(obj_type, fm.CELL_METADATA[obj_type], values_to_add, single_obj, connection)
						else:
							fm.gather_pooled_metadata(obj_type, fm.CELL_METADATA[obj_type], values_to_add, objs, connection)
					else:
						fm.gather_pooled_metadata(obj_type, fm.CELL_METADATA[obj_type], values_to_add, objs, connection)
		row_to_add = pd.DataFrame(values_to_add, index=[mxr['@id']], dtype=str)
		
		# Add anndata to list of final raw anndatas, only for RNAseq
		if summary_assay in ['RNA', 'CITE']:
			# Checking for presence of mxr file and downloading if not present
			if mxr['s3_uri'].endswith('h5'):
				if os.path.exists(fm.MTX_DIR + '/' + mxr_acc + '.h5'):
					print(mxr_acc + '.h5' + ' was found locally')
				else:
					download_file(mxr.get('s3_uri'), fm.MTX_DIR, mxr.get('accession'))
			elif mxr['s3_uri'].endswith('h5ad'):
				if os.path.exists(fm.MTX_DIR + '/' + mxr_acc + '.h5ad'):
					print(mxr_acc + '.h5ad' + ' was found locally')
				else:
					download_file(mxr.get('s3_uri'), fm.MTX_DIR, mxr.get('accession'))

			mxr_name = '{}.h5'.format(mxr_acc) if mxr['s3_uri'].endswith('h5') else '{}.h5ad'.format(mxr_acc)
			if glob.mfinal_obj.get('spatial_s3_uri', None) and glob.mfinal_obj['assays'] == ['spatial transcriptomics']:
				# Checking for presence of spatial directory and redownloading if present
				if os.path.exists(fm.MTX_DIR + '/spatial'):
					shutil.rmtree(fm.MTX_DIR + '/spatial')
				download_directory(glob.mfinal_obj['spatial_s3_uri'], fm.MTX_DIR)
				# If tissue_positions is present rename to tissue_positions_list and remove header
				if os.path.exists(fm.MTX_DIR + '/spatial/tissue_positions.csv') == True:
					fixed_file = pd.read_csv(fm.MTX_DIR + '/spatial/tissue_positions.csv', skiprows=1, header=None)
					fixed_file.to_csv(fm.MTX_DIR + '/spatial/tissue_positions_list.csv', header=False, index=False)
					os.remove(fm.MTX_DIR + '/spatial/tissue_positions.csv')
				if 'spatial' in glob.mfinal_adata.uns.keys():
					del glob.mfinal_adata.uns['spatial']
				adata_raw = sq.read.visium(fm.MTX_DIR, counts_file=mxr_name)
				if adata_raw.obs.shape[0] < 4992 and len(glob.mfinal_obj.get('libraries'))==1:
					adata_raw = add_raw_background_spots(adata_raw, glob)

			elif mxr['s3_uri'].endswith('h5'):
				adata_raw = sc.read_10x_h5('{}/{}'.format(fm.MTX_DIR, mxr_name), gex_only=False)
			elif mxr['s3_uri'].endswith('h5ad'):
				adata_raw = sc.read_h5ad('{}/{}'.format(fm.MTX_DIR, mxr_name))
			else:
				logging.error('ERROR: Raw matrix file of unknown file extension: {}'.format(mxr['s3_uri']))
				sys.exit('ERROR: Raw matrix file of unknown file extension: {}'.format(mxr['s3_uri']))

			if summary_assay == 'RNA':
				row_to_add['gene_annotation_version'] = mxr['genome_annotation']
				adata_raw = adata_raw[:, adata_raw.var['feature_types']=='Gene Expression']
			else:
				adata_raw = adata_raw[:, adata_raw.var['feature_types']=='Antibody Capture']
			# only make var unique if all raw matrices are same annotation version
			if len(glob.mfinal_obj.get('genome_annotations', [])) == 1:
				for g in [i for i,c in collections.Counter(adata_raw.var.index.to_list()).items() if c > 1]:
					if True in adata_raw.var.loc[g, 'gene_ids'].str.endswith("PAR_Y").to_list():
						redundant.extend([i for i in adata_raw.var.loc[g, 'gene_ids'] if i.endswith('PAR_Y')])
					else:
						redundant.extend(adata_raw.var.loc[g, 'gene_ids'].to_list())
				adata_raw.var_names_make_unique(join='.')
			# Recreate cell_ids and subset raw matrix and add mxr_acc into obs
			if glob.mfinal_obj.get('cell_label_mappings', None):
				concatenated_ids = concatenate_cell_id(mxr['@id'], adata_raw.obs_names, glob)
				adata_raw.obs_names = concatenated_ids
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
			else:
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(adata_raw.obs_names.to_list()))

			# Error check to see that cells in raw matrix match the cell in mfinal_adata
			cell_mapping_dct = {}
			if glob.mfinal_obj.get('cell_label_mappings'):
				for mapping_dict in glob.mfinal_obj.get('cell_label_mappings'):
					cell_mapping_dct[mapping_dict['raw_matrix']] = mapping_dict['label']
				mapping_label = cell_mapping_dct[mxr.get('@id')]
			if glob.mfinal_obj.get('cell_label_location') == 'prefix':
				prefixes = []
				for mapping_dict in glob.mfinal_obj['cell_label_mappings']: # Creating list of all prefixes
					prefixes.append(mapping_dict['label'])
				prefixes.remove(mapping_label) # Removing mapping_label prefix from list of all prefixes
				# Checking to make sure none of the other prefixes contain the mapping_label prefix, if they do, then make sure there's no false match
				if any(prefix.startswith(mapping_label) for prefix in prefixes):
					mfinal_with_label = [i for i in mfinal_cell_identifiers if i.startswith(mapping_label) and not i.startswith(tuple(prefixes))]
				else:
					mfinal_with_label = [i for i in mfinal_cell_identifiers if i.startswith(mapping_label)]
			elif glob.mfinal_obj.get('cell_label_location') == 'suffix':
				mfinal_with_label = [i for i in mfinal_cell_identifiers if i.endswith(mapping_label)]
			else: 
				mfinal_with_label = [i for i in mfinal_cell_identifiers if i in adata_raw.obs_names]
			if len(overlapped_ids) == 0:
				if glob.mfinal_obj['cell_label_location'] == 'prefix':
					if concatenated_ids[0].endswith('-1'):
						concatenated_ids = [re.sub('-1$', '', i) for i in concatenated_ids]
					else:
						concatenated_ids = [i+'-1' for i in concatenated_ids]
					overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
					if len(overlapped_ids) == 0:
						mapping_error = True
						error_info[mxr.get('@id')] = mapping_label
						error_info[mxr.get('@id')] += '; zero_prefix: {}'.format(concatenated_ids[0:5])
					adata_raw.obs_names = concatenated_ids
				else:
					mapping_error = True
					error_info[mxr.get('@id')] = mapping_label
					error_info[mxr.get('@id')] += '; zero: {}'.format(concatenated_ids[0:5])
			elif len(overlapped_ids) != len(mfinal_with_label):
				mapping_error = True
				error_info[mxr.get('@id')] = mapping_label
				error_info[mxr.get('@id')] += '; in_adata: {}, in_raw: {}, overlap: {}'.format(len(mfinal_with_label), len(concatenated_ids), len(overlapped_ids))
			if not mapping_error:
				if not glob.mfinal_obj.get('spatial_s3_uri', None):
					adata_raw = adata_raw[overlapped_ids]
				adata_raw.obs['raw_matrix_accession'] = mxr['@id']
				cxg_adata_lst.append(adata_raw)
       
		df = pd.concat([df, row_to_add])
		redundant = list(set(redundant))
		
	# Removing gene_annotation_version if genome_annotations from ProcMatrixFile is empty
	if not glob.mfinal_obj.get('genome_annotations', None):
		del df['gene_annotation_version']

	if mapping_error:
		logging.error('ERROR: There are {} mapping errors in cell_label_mappings:'.format(len(error_info.keys())))
		print("ERROR: There are {} mapping errors in cell_label_mappings:".format(len(error_info.keys())))
		for er in error_info.keys():
			print("RawMatrixFile: {}, {}".format(er, error_info[er]))
			logging.error("RawMatrixFile: {}, {}".format(er, error_info[er]))
		sys.exit()

	# Get dataset-level metadata and set 'is_primary_data' for obs accordingly as boolean
	ds_results = {}
	ds_results = fm.gather_metdata('matrix', fm.DATASET_METADATA['final_matrix'], ds_results, [glob.mfinal_obj], connection)
	df['is_primary_data'] = ds_results['is_primary_data']
	df['is_primary_data'].replace({'True': True, 'False': False}, inplace=True)

	# Check if default_embedding is unreported_value, and if so remove
	if ds_results['default_embedding'] == fm.UNREPORTED_VALUE:
		del ds_results['default_embedding']

	del ds_results['is_primary_data']

	# Should add error checking to make sure all matrices have the same number of vars
	feature_lengths = []
	for adata in cxg_adata_lst:
		feature_lengths.append(adata.shape[1])
	feature_lengths = list(set(feature_lengths))

	# Set up dataframe for cell annotations keyed off of author_cell_type
	annot_df = pd.DataFrame()
	for annot_obj in glob.mfinal_obj['cell_annotations']:
		annot_lst = []
		annot_lst.append(annot_obj)
		annot_metadata = {}
		annot_metadata = fm.gather_metdata('cell_annotation', fm.ANNOT_FIELDS, annot_metadata, annot_lst, connection)
		annot_row = pd.DataFrame(annot_metadata, index=[annot_obj['author_cell_type']])
		annot_df = pd.concat([annot_df, annot_row])

	# For CITE and RNA datasets, concatenate all anndata objects in list, but no reconciling genes for CITE
	# For ATAC datasets, assumption is that there is no scale.data, and raw count is taken from mfinal_adata.raw.X
	raw_matrix_mapping = []
	cell_mapping_rev_dct = {}
	if summary_assay == 'RNA':
		# If raw matrices are annotated to multiple gencode versions, concatenate on ensembl ID and remove ambiguous symbols
		if len(glob.mfinal_obj.get('genome_annotations', [])) > 1:
			glob.cxg_adata_raw, redundant, all_remove  = reconcile_genes(cxg_adata_lst, glob)
			drop_removes = set(glob.mfinal_adata.var.index.to_list()).intersection(set(all_remove))
			logging.info('drop_all_removes:\t{}\t{}'.format(len(drop_removes), drop_removes))
			glob.mfinal_adata = glob.mfinal_adata[:, [i for i in glob.mfinal_adata.var.index.to_list() if i not in all_remove]]
		elif len(feature_lengths) > 1:
			glob.cxg_adata_raw = concat_list(cxg_adata_lst, 'gene_ids', True)
		else:
			glob.cxg_adata_raw = concat_list(cxg_adata_lst, 'none', True)
			if len(feature_lengths) == 1:
				if glob.cxg_adata_raw.var.shape[0] != feature_lengths[0]:
					logging.error('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
					sys.exit('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
		if not glob.mfinal_obj.get('spatial_s3_uri', None):
			glob.cxg_adata_raw = glob.cxg_adata_raw[mfinal_cell_identifiers]
			if glob.cxg_adata_raw.shape[0] != glob.mfinal_adata.shape[0]:
				logging.error('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
				sys.exit('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
	elif summary_assay == 'CITE':
		glob.cxg_adata_raw = concat_list(cxg_adata_lst, 'none', False)
		if len(feature_lengths) == 1:
			if glob.cxg_adata_raw.var.shape[0] != feature_lengths[0]:
				logging.error('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
				sys.exit('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
		glob.cxg_adata_raw = glob.cxg_adata_raw[mfinal_cell_identifiers]
		if glob.cxg_adata_raw.shape[0] != glob.mfinal_adata.shape[0]:
			logging.error('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
			sys.exit('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
	elif summary_assay == 'ATAC':
		for mapping_dict in glob.mfinal_obj['cell_label_mappings']:
			cell_mapping_rev_dct[mapping_dict['label']] = mapping_dict['raw_matrix']
		flag_removed = False
		for final_id in mfinal_cell_identifiers:
			if not re.search('[AGCT]+-1', final_id):
				flag_removed = True
		for cell_id in mfinal_cell_identifiers:
			if glob.mfinal_obj['cell_label_location'] == 'prefix':
				label = re.search('^(.*)[AGCT]{16}.*$', cell_id).group(1)
			elif glob.mfinal_obj['cell_label_location'] == 'suffix':
				if flag_removed:
					label = re.search(r'^[AGCT]+(.*)$', cell_id).group(1)
				else:
					label = re.search(r'^[AGCT]+-1(.*)$', cell_id).group(1)
			raw_matrix_mapping.append(cell_mapping_rev_dct[label])
		atac_obs = pd.DataFrame({'raw_matrix_accession': raw_matrix_mapping}, index=mfinal_cell_identifiers)
		if glob.mfinal_adata.raw == None:
			glob.cxg_adata_raw = ad.AnnData(sparse.csr_matrix(glob.mfinal_adata.X.shape), var=glob.mfinal_adata.var, obs=atac_obs)
		else:
			glob.cxg_adata_raw = ad.AnnData(glob.mfinal_adata.raw.X, var=glob.mfinal_adata.var, obs=atac_obs)

	# Set uns and obsm parameters (obsm is originally set from mfinal_adata, since there may be other embeddings)
	glob.cxg_uns = ds_results
	glob.cxg_obsm = glob.mfinal_adata.obsm.copy()
	drop_obsm = []
	for embed in glob.cxg_obsm.keys():
		if type(glob.cxg_obsm[embed]) == np.ndarray:
			if glob.cxg_obsm[embed].ndim >= 2:
				if embed.startswith('X_'):
					if glob.cxg_obsm[embed].shape[1] < 2:
						warning_list.append("WARNING: Embedding starting with 'X_' that has length < 2 for second dimension is dropped: {}\t{}".format(glob.cxg_obsm[embed].shape, embed))
						drop_obsm.append(embed)
			else:
				warning_list.append("WARNING: Embedding that is not 2 dimensions is dropped: {}\t{}".format(glob.cxg_obsm[embed].shape, embed))
				drop_obsm.append(embed)
		else:
			warning_list.append("WARNING: Embedding that is not a numpy array is dropped: {}\t{}".format(type(glob.cxg_obsm[embed]), embed))
			drop_obsm.append(embed)
	for k in drop_obsm:
		del glob.cxg_obsm[k]
	if len([i for i in glob.cxg_obsm.keys() if i.startswith('X_')]) < 1:
		if glob.mfinal_obj['assays'] != ['spatial transcriptomics']:
			logging.error("ERROR: At least one embedding that starts with 'X_' is required")
			sys.exit("ERROR: At least one embedding that starts with 'X_' is required")

	# Merge df with raw_obs according to raw_matrix_accession, and add additional cell metadata from mfinal_adata if available
	# Also add calculated fields to df 
	celltype_col = glob.mfinal_obj['author_cell_type_column']
	if glob.mfinal_adata.obs[celltype_col].dtype != 'object':
		glob.mfinal_adata.obs[celltype_col] = glob.mfinal_adata.obs[celltype_col].astype('string')
	glob.cxg_obs = pd.merge(glob.cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='left')
	glob.cxg_obs = pd.merge(glob.cxg_obs, glob.mfinal_adata.obs[[celltype_col]], left_index=True, right_index=True, how='left')
	glob.cxg_obs = pd.merge(glob.cxg_obs, annot_df, left_on=celltype_col, right_index=True, how='left')
	glob.cxg_obs['cell_state'].fillna(fm.UNREPORTED_VALUE, inplace=True)

	if not glob.mfinal_obj.get('spatial_s3_uri', None):
		if glob.cxg_obs['cell_type_ontology_term_id'].isnull().values.any():
			warning_list.append("WARNING: Cells did not sucessfully map to CellAnnotations with author cell type and counts: {}".\
				format(glob.cxg_obs.loc[glob.cxg_obs['cell_type_ontology_term_id'].isnull()==True, celltype_col].value_counts().to_dict()))
		if glob.cxg_obs[celltype_col].isna().any():
			logging.error("ERROR: author_cell_type column contains 'NA' values, unable to perform CellAnnotation mapping.")
			sys.exit("ERROR: author_cell_type column contains 'NA' values, unable to perform CellAnnotation mapping.")
		if len([i for i in annot_df.index.to_list() if i not in glob.cxg_obs[celltype_col].unique().tolist()]) > 0:
			warning_list.append("WARNING: CellAnnotation that is unmapped: {}\n".format([i for i in annot_df.index.to_list() if i not in glob.cxg_obs[celltype_col].unique().tolist()]))

	if 'author_cluster_column' in glob.mfinal_obj:
		cluster_col = glob.mfinal_obj['author_cluster_column']
		glob.cxg_obs = pd.merge(glob.cxg_obs, glob.mfinal_adata.obs[[cluster_col]], left_index=True, right_index=True, how='left')
		glob.cxg_obs.rename(columns={cluster_col: 'author_cluster'}, inplace=True)
		glob.cxg_obs['author_cluster'] = glob.cxg_obs['author_cluster'].astype('category')

	# After getting experimental metadata keyed off of mxr, if there is demultiplexed_donor_column, run demultiplex
	if 'demultiplexed_donor_column' in glob.mfinal_obj:
		donor_col = glob.mfinal_obj['demultiplexed_donor_column']
		glob.cxg_obs = pd.merge(glob.cxg_obs, glob.mfinal_adata.obs[[donor_col]], left_index=True, right_index=True, how='left')
		glob.cxg_obs.rename(columns={donor_col: 'author_donor'}, inplace=True)
		glob.cxg_obs['library_@id'] = glob.cxg_obs['library_@id'].astype(str)
		glob.cxg_obs['author_donor'] = glob.cxg_obs['author_donor'].astype(str)
		glob.cxg_obs['library_authordonor'] = glob.cxg_obs['library_@id'] + ',' + glob.cxg_obs['author_donor']

		lib_donor_df = glob.cxg_obs[['library_@id', 'author_donor', 'library_authordonor']].drop_duplicates().reset_index(drop=True)
		donor_df = demultiplex(lib_donor_df, library_susp, donor_susp, glob)

		report_diseases(donor_df, glob.mfinal_obj.get('experimental_variable_disease', fm.UNREPORTED_VALUE))
		get_sex_ontology(donor_df)

		# Retain cell identifiers as index
		glob.cxg_obs = glob.cxg_obs.reset_index().merge(donor_df, how='left', on='library_authordonor').set_index('index')
		if glob.mfinal_adata.X.shape[0] != glob.cxg_obs.shape[0]:
			logging.error('ERROR: cxg_obs does not contain the same number of rows as final matrix: {} vs {}'.format(glob.mfinal_adata.X.shape[0], glob.cxg_obs.shape[0]))
			sys.exit('ERROR: cxg_obs does not contain the same number of rows as final matrix: {} vs {}'.format(glob.mfinal_adata.X.shape[0], glob.cxg_obs.shape[0]))
	else:
		# Go through donor and biosample diseases and calculate cxg field accordingly
		report_diseases(df, glob.mfinal_obj.get('experimental_variable_disease', fm.UNREPORTED_VALUE))
		get_sex_ontology(df)
		glob.cxg_obs = pd.merge(glob.cxg_obs, df[['disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id']], left_on="raw_matrix_accession", right_index=True, how="left")

	# Clean up columns in obs to follow cxg schema and drop any unnecessary fields
	drop_cols(celltype_col, glob)
	clean_obs(glob)
	del cxg_adata_lst
	gc.collect()

	# Check that primary_portion.obs_field of ProcessedMatrixFile is present in cxg_obs
	if glob.mfinal_obj.get('primary_portion', None): # Checking for presence of 'primary_portion'
		primary_portion = glob.mfinal_obj.get('primary_portion')
		if primary_portion.get('obs_field') not in glob.cxg_obs.columns:
			logging.error("ERROR: 'obs_field' value '{}' not found in cxg_obs columns".format(primary_portion.get('obs_field')))
			sys.exit("ERROR: 'obs_field' value '{}' not found in cxg_obs columns".format(primary_portion.get('obs_field')))

		# Check that all primary_portion.values of ProcessedMatrixFile are found in the 'obs_field' column of cxg_obs
		missing = [f for f in primary_portion.get('values') if f not in glob.cxg_obs[primary_portion.get('obs_field')].tolist()]
		if missing:
			logging.error("ERROR: cxg_obs column '{}' doesn't contain values present in 'primary_portion.obs_field' of ProcessedMatrixFile: {}".format(primary_portion.get('obs_field'),missing))
			sys.exit("ERROR: cxg_obs column '{}' doesn't contain values present in 'primary_portion.obs_field' of ProcessedMatrixFile: {}".format(primary_portion.get('obs_field'),missing))

	# If final matrix file is h5ad, take expression matrix from .X to create cxg anndata
	results_file  = get_results_filename(glob)
	glob.mfinal_adata.var_names_make_unique()
	cxg_var = pd.DataFrame(index=glob.mfinal_adata.var.index.to_list())
	keep_types = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
	if summary_assay == 'CITE':
		keep_types.append('object')
	var_meta = glob.mfinal_adata.var.select_dtypes(include=keep_types)

	# Add spatial information to adata.uns, which is assay dependent. Assumption is that the spatial dataset is from a single assay
	if glob.mfinal_obj['assays'] == ['spatial transcriptomics']:
		warnings = fm.process_spatial(glob)
		if warnings:
			warning_list.append(warnings)
	# Check to see if need to add background spots
	if len(glob.mfinal_obj.get('libraries'))==1 and glob.mfinal_obj.get('spatial_s3_uri', None):
		add_background_spots(glob)
	glob.cxg_adata = ad.AnnData(glob.mfinal_adata.X, obs=glob.cxg_obs, obsm=glob.cxg_obsm, var=cxg_var, uns=glob.cxg_uns)
	glob.cxg_adata.var = glob.cxg_adata.var.merge(var_meta, left_index=True, right_index=True, how='left')
	
	# Removing feature_length column from var if present
	if 'feature_length' in glob.cxg_adata.var.columns:
		adata.var.drop(columns=['feature_length'], inplace=True)

	# Check matrix density
	glob.cxg_adata.X = check_matrix(glob.cxg_adata.X)

	# Check that cxg_adata_raw.X is correct datatype
	if not glob.cxg_adata_raw.X.dtype == 'float32':
		glob.cxg_adata_raw.X = glob.cxg_adata_raw.X.astype(np.float32) 
		
	# Adding layers from 'layers_to_keep' to cxg_adata.layers	
	if 'layers_to_keep' in glob.mfinal_obj:
		for k in glob.mfinal_obj['layers_to_keep']:
			glob.cxg_adata.layers[k] = glob.mfinal_adata.layers[k]
			glob.cxg_adata.layers[k] = check_matrix(glob.cxg_adata.layers[k])
				
	# Convert gene symbols to ensembl and filter to approved set
	if len(feature_lengths) > 1 and len(glob.mfinal_obj['genome_annotations'])==1:
		clean_var(glob)

	# For ATAC gene activity matrices, it is assumed there are no genes that are filtered
	# For CITE, standardize antibody index and metadata and no filtering
	if summary_assay == 'RNA':
		compiled_annot = compile_annotations(fm.REF_FILES)
		set_ensembl(redundant, glob)
		add_zero(glob)
		glob.cxg_adata_raw = filter_ensembl(glob.cxg_adata_raw, compiled_annot)
		glob.cxg_adata = filter_ensembl(glob.cxg_adata, compiled_annot)
	elif summary_assay == 'ATAC':
		compiled_annot = compile_annotations(fm.REF_FILES)
		glob.cxg_adata_raw = filter_ensembl(glob.cxg_adata_raw, compiled_annot)
		glob.cxg_adata = filter_ensembl(glob.cxg_adata, compiled_annot)
		glob.cxg_adata.var['feature_is_filtered'] = False
	elif summary_assay == 'CITE':
		add_labels(glob)
		map_antibody(glob)
		add_zero(glob)

	# Copy over any additional data from mfinal_adata to cxg_adata
	reserved_uns = ['schema_version', 'title', 'default_embedding', 'X_approximate_distribution', 'schema_reference', 'citation']
	warnings = fm.copy_over_uns(glob, reserved_uns)
	if warnings:
		warning_list.append(warnings)
	

	if glob.mfinal_adata.obsp:
		glob.cxg_adata.obsp = glob.mfinal_adata.obsp

	# Check if mfinal_obj matrix is normalized,if so set cxg_adata.raw to raw, if not then place raw in adata.X
	if glob.mfinal_obj['X_normalized']:
		if summary_assay != 'ATAC' or glob.mfinal_adata.raw != None:
			glob.cxg_adata.raw = glob.cxg_adata_raw
	else:
		glob.cxg_adata.var['feature_is_filtered'] = False
		glob.cxg_adata = ad.AnnData(glob.cxg_adata_raw.X, obs=glob.cxg_adata.obs, obsm=glob.cxg_adata.obsm, var=glob.cxg_adata.var, uns=glob.cxg_adata.uns)
	quality_check(glob)
	glob.cxg_adata.write_h5ad(results_file, compression='gzip')

	# Printing out list of warnings
	for n in warning_list:
		if 'WARNING: Full list' not in n:
			print(n, end='\n')
		if 'Full list available in logging file.' not in n:
			logging.warning(n)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
