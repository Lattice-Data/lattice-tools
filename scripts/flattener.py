import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
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


cell_metadata = {
	'donor': [
		'uuid',
		'age_display',
		'sex',
		'ethnicity.term_name',
		'ethnicity.term_id',
		'life_stage',
		'life_stage_term_id',
		'diseases.term_id',
		'diseases.term_name',
		'body_mass_index',
		'times_pregnant',
		'family_history_breast_cancer'
		],
	'sample': [
		'uuid',
		'preservation_method',
		'biosample_ontology.term_name',
		'biosample_ontology.term_id',
		'diseases.term_id',
		'diseases.term_name'
		],
	'suspension': [
		'uuid',
		'suspension_type',
		'@id'
		],
	'library': [
		'uuid',
		'protocol.title',
		'protocol.term_id',
		'@id'
		]
	}

dataset_metadata = {
	'final_matrix': [
		'description',
		'genome_annotation',
		'default_visualization',
		'default_embedding'
		]
	}

annot_fields = [
	'cell_ontology.term_name',
	'cell_ontology.term_id',
	'author_cell_type',
	'cell_ontology.cell_slims'
]

prop_map = {
	'sample_biosample_ontology_term_name': 'tissue',
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'library_protocol_title': 'assay',
	'library_protocol_term_id': 'assay_ontology_term_id',
	'donor_body_mass_index': 'donor_BMI',
	'donor_sex': 'sex',
	'donor_ethnicity_term_name': 'ethnicity',
	'donor_ethnicity_term_id': 'ethnicity_ontology_term_id',
	'donor_life_stage': 'development_stage',
	'donor_life_stage_term_id': 'development_stage_ontology_term_id',
	'donor_age_display': 'donor_age',
	'donor_family_history_breast_cancer': 'family_history_breast_cancer',
	'matrix_genome_annotation': 'reference_annotation_version',
	'matrix_description': 'title',
	'cell_annotation_author_cell_type': 'author_cell_type',
	'cell_annotation_cell_ontology_term_id': 'cell_type_ontology_term_id',
	'cell_annotation_cell_ontology_term_name': 'cell_type',
	'matrix_default_visualization': 'default_field',
	'matrix_default_embedding': 'default_embedding',
	'cell_annotation_cell_ontology_cell_slims': 'cell_type_category',
	'suspension_suspension_type': 'suspension_type'
}

unreported_value = ''
corpora_schema_version = '1.1.0'
corpora_encoding_version = '0.1.0'
flat_version = '2'

EPILOG = '''
Examples:

    python %(prog)s -m local -f LATDF119AAA

For more details:

    python %(prog)s --help
'''

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


def gather_rawmatrices(derived_from):
	my_raw_matrices = []
	df_ids = []
	for identifier in derived_from:
		obj = lattice.get_object(identifier, connection)
		if obj['@type'][0] == 'MatrixFile' and obj['layers'][0]['normalized'] != True and \
				obj['layers'][0]['value_scale'] == 'linear' and len(obj['layers']) == 1 and \
				'cell calling' in obj['derivation_process']:
			my_raw_matrices.append(obj)
		else:
			# grab the derived_from in case we need to go a layer deeper
			for i in obj['derived_from']:
				df_ids.append(i)
	if not my_raw_matrices:
		for identifier in df_ids:
			obj = lattice.get_object(identifier, connection)
			if obj['@type'][0] == 'MatrixFile' and obj['normalized'] != True and \
				obj['value_scale'] == 'linear' and obj['file_format'] == 'hdf5' and \
				'cell calling' in obj['derivation_process']:
				my_raw_matrices.append(obj)
	return my_raw_matrices


def gather_objects(input_object, start_type=None):
	if start_type == None:
		lib_ids = input_object['libraries']
	libraries = []
	susp_ids = []
	suspensions = []
	prep_susp_ids = []
	prepooled_susps = []
	sample_ids = []
	samples = []
	donor_ids = []
	donors = []

	if start_type == None:
		for i in lib_ids:
			obj = lattice.get_object(i, connection)
			libraries.append(obj)
			for o in obj['derived_from']:
				if o.get('uuid') not in susp_ids:
					suspensions.append(o)
					susp_ids.append(o.get('uuid'))
			for o in obj['donors']:
				if o.get('uuid') not in donor_ids:
					donors.append(o)
					donor_ids.append(o.get('uuid'))
	elif start_type == 'suspension':
		susp_ids = [input_object['uuid']]
		suspensions = [input_object]
		for susp in suspensions:
			for o in susp['donors']:
				if o.get('uuid') not in donor_ids:
					donors.append(o)
					donor_ids.append(o.get('uuid'))

	for o in suspensions:
		for i in o['derived_from']:
			sample_ids.append(i)
	remaining = set(sample_ids)
	seen = set()
	while remaining:
		seen.update(remaining)
		next_remaining = set()
		for i in remaining:
			obj = lattice.get_object(i, connection)
			if 'Biosample' in obj['@type']:
				samples.append(obj)
			else:
				if 'Suspension' in obj['@type'] and obj['uuid'] not in prep_susp_ids:
					prepooled_susps.append(obj)
					next_remaining.update(obj['derived_from'])
					prep_susp_ids.append(obj['uuid'])
		remaining = next_remaining - seen

	objs = {
		'donor': donors,
		'sample': samples,
		'suspension': suspensions
		}
	if start_type == None:
		objs['library'] = libraries
	if prepooled_susps:
		objs['prepooled_suspension'] = prepooled_susps
		objs['pooled_suspension'] = objs.pop('suspension')
	return objs


def get_value(obj, prop):
	path = prop.split('.')
	if len(path) == 1:
		return obj.get(prop, unreported_value)
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2, unreported_value) for i in obj[key1]]
			return list(set(values))
		elif obj.get(key1):
			value = obj[key1].get(key2, unreported_value)
			if key1 == 'biosample_ontology' and 'Culture' in obj['@type']:
				obj_type = obj['@type'][0]
				if obj_type == 'Organoid':
					obj_type_conv = 'organoid'
				elif obj_type == 'CellCulture':
					obj_type_conv = 'cell culture'
				return  '{} ({})'.format(value, obj_type_conv)
			else:
				return value
		else:
			return obj.get(key1,unreported_value)
	else:
		return 'unable to traverse more than 1 embedding'


def gather_metdata(obj_type, properties, values_to_add, objs):
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = value


def gather_pooled_metadata(obj_type, properties, values_to_add, objs):
	for prop in properties:
		value = list()
		for obj in objs:
			v = get_value(obj, prop)
			if isinstance(v, list):
				value.extend(v)
			else:
				value.append(v)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		value_str = [str(i) for i in value]
		value_set = set(value_str)
		if len(value_set) > 1:
			value_str = [re.sub(r'^$', 'unknown', i) for i in value_str]
			values_to_add[key] = 'pooled samples: [{}]'.format(','.join(value_str))
		else:
			values_to_add[key] = next(iter(value_set))

def report_dataset(donor_objs, matrix, dataset):
	ds_results = {}
	ds_obj = lattice.get_object(dataset, connection)
	for prop in dataset_metadata['final_matrix']:
		value = get_value(matrix, prop)
		if isinstance(value, list):
			value = ','.join(value)
		if value != unreported_value:
			latkey = 'matrix_' + prop.replace('.','_')
			key = prop_map.get(latkey, latkey)
			ds_results[key] = value
	layer_descs = {
		'raw.X': 'raw'
	}
	derived_by = matrix.get('derivation_process')
	for layer in matrix.get('layers'):
		units = layer.get('value_units', unreported_value)
		scale = layer.get('value_scale', unreported_value)
		norm_meth = layer.get('normalization_method', unreported_value)
		desc = '{} counts; {} scaling; normalized using {}; derived by {}'.format(units, scale, norm_meth, ', '.join(derived_by))
		layer_descs['X'] = desc
	ds_results['layer_descriptions'] = layer_descs
	org_id = set()
	org_name = set()
	for obj in donor_objs:
		org_id.add(obj['organism']['taxon_id'])
		org_name.add(obj['organism']['scientific_name'])
	ds_results['organism_ontology_term_id'] = ','.join(org_id)
	ds_results['organism'] = ','.join(org_name)
	return ds_results


def download_file(file_obj, directory):
	if file_obj.get('s3_uri'):
		download_url = file_obj.get('s3_uri')
		bucket_name = download_url.split('/')[2]
		file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
		file_ext = download_url.split('.')[-1]
		s3client = boto3.client("s3")
		file_name = file_obj.get('accession') + '.' + file_ext
		print(file_name + ' downloading')
		try:
			s3client.download_file(bucket_name, file_path, directory + '/' + file_name)
		except subprocess.CalledProcessError as e:
		 	sys.exit('Failed to find file {} on s3'.format(file_obj.get('@id')))
		else:
			print(file_name + ' downloaded')
	elif file_obj.get('external_uri'):
		download_url = file_obj.get('external_uri')
		ftp_server = download_url.split('/')[2]
		ftp = FTP(ftp_server)
		ftp.login(user='anonymous', passwd = 'password')
		file_path = download_url.replace('ftp://{}/'.format(ftp_server), '')
		file_ext = download_url.split('.')[-1]
		file_name = file_obj.get('accession') + '.' + file_ext
		print(file_name + ' downloading')
		try:
			ftp.retrbinary('RETR ' + file_path, open(directory + '/' + file_name, 'wb').write)
		except error_perm as e:
			os.remove(file_name)
			sys.exit(e)
		else:
			ftp.quit()
			print(file_name + ' downloaded')
	else:
		sys.exit('File {} has no uri defined'.format(file_obj['@id']))


# Recreated the final matrix ids, also checking to see if '-1' was removed from original cell identifier
def concatenate_cell_id(mfinal_obj, mxr_acc, raw_obs_names, mfinal_cells):
	new_ids = []
	flag_removed = False
	cell_mapping_dct = {}
	cell_location = mfinal_obj['cell_label_location']
	for mapping_dict in mfinal_obj['cell_label_mappings']:
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


# Get cell embeddings from final matrix object
# Skip pca or harmony embeddings, and only transfer umap and tsne embeddings
def get_embeddings(mfinal_adata):
	if len(mfinal_adata.obsm_keys()) == 0:
		sys.exit('At least 1 set of cell embeddings is required in final matrix')
	final_embeddings = mfinal_adata.obsm.copy()
	all_embedding_keys = mfinal_adata.obsm_keys()
	for embedding in all_embedding_keys:
		if embedding == 'X_pca' or embedding == 'X_harmony':
			final_embeddings.pop(embedding)
		elif embedding != 'X_umap' and embedding != 'X_tsne':
			sys.exit('There is an unrecognized embedding in final matrix: {}'.format(embedding))
	return final_embeddings


# From R object, create and return h5ad in temporary drive
# Returns list of 
def convert_from_rds(path_rds, assays, temp_dir, cell_col):
	converted_h5ad = []
	utils = rpackages.importr('utils')
	#base = rpackages.importr('base')
	utils.chooseCRANmirror(ind=1)
	packnames = ('Seurat', 'reticulate', 'BiocManager', 'devtools')
	names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
	if len(names_to_install) > 0:
		utils.install_packages(StrVector(names_to_install))

	devtools = rpackages.importr('devtools')
	github_packages = ('cli', 'crayon', 'hdf5r', 'Matrix', 'R6', 'rlang', 'stringi', 'withr')
	github_install = [x for x in github_packages if not rpackages.isinstalled(x)]
	if len(github_install) > 0:
		devtools.install_github(StrVector(bioc_to_install))
	seurat = rpackages.importr('Seurat')
	seuratdisk = rpackages.importr('SeuratDisk')
	
	# Load Seurat object into a h5Seurat file, and create h5ad for each assay
	# If scaled.data is present, data is in raw.X, else data is in X
	robjects.r('print(sessionInfo())')
	h5s_file = '{}/r_obj.h5Seurat'.format(temp_dir)
	robjects.r('robj <- readRDS("{}")'.format(path_rds))
	robjects.r('updated_robj <- Seurat::UpdateSeuratObject(object = robj)')
	for assay in assays:
		robjects.r('robj_assay <- Seurat::GetAssay(updated_robj, assay="{}")'.format(assay))
		robjects.r('robj_new <- Seurat::CreateSeuratObject(robj_assay, assay = "{}", meta.data = updated_robj@meta.data)'.format(assay))
		robjects.r('reducs <- names(updated_robj@reductions)')
		reductions = robjects.r('print(reducs)')
		for reduc in reductions:
			robjects.r('robj_reduc <- Seurat::CreateDimReducObject(embeddings = Embeddings(updated_robj, reduction="{}"), loadings = Loadings(updated_robj, reduction="umap"), global=TRUE, assay = "{}")'.format(reduc, assay))
			robjects.r('robj_new@reductions${} <- robj_reduc'.format(reduc))
		robjects.r('robj_new@meta.data${} <- as.character(robj_new@meta.data${})'.format(cell_col, cell_col))
		robjects.r('SeuratDisk::SaveH5Seurat(robj_new, filename="{}")'.format(h5s_file))
		scaled_matrix = robjects.r('dim(Seurat::GetAssayData(object = robj_new, assay="{}", slot="scale.data"))'.format(assay))
		h5ad_file = '{}/{}.h5ad'.format(temp_dir, assay)
		seuratdisk.Convert(h5s_file, dest="h5ad", assay=assay, overwrite = 'FALSE')
		os.rename(h5s_file.replace('h5Seurat', 'h5ad'), h5ad_file)
		print("Converting to h5ad: {}".format(h5ad_file))
		if scaled_matrix[0] == 0 and scaled_matrix[0] == 0:
			converted_h5ad.append((h5ad_file, 'X', assay))
		else:
			converted_h5ad.append((h5ad_file,'raw.X', assay))
	return converted_h5ad


# If cell slims is a list, narrow down to one ontology
def trim_cell_slims(df_annot):
	cell_term_list = df_annot['cell_type_category'].tolist()
	for i in range(len(cell_term_list)):
		cell_terms = cell_term_list[i].split(',')
		if len(cell_terms) > 1:
			if df_annot.iloc[i]['cell_type'] in cell_terms:
				for index in range(len(cell_terms)):
					if df_annot.iloc[i]['cell_type'] == cell_terms[index]:
						if index == len(cell_terms) - 1:
							print("WARNING, there is a cell_slims that is more than a single ontology: {} changed to {}".format(cell_term_list[i], cell_terms[index]))
							cell_term_list[i] = cell_terms[index]
						else:
							print("WARNING, there is a cell_slims that is more than a single ontology: {} changed to {}".format(cell_term_list[i], cell_terms[index+1]))
							cell_term_list[i] = cell_terms[index+1]
			else:
				print("WARNING, there is a cell_slims that is more than a single ontology: {} changed to {}".format(cell_term_list[i], cell_terms[0]))
				cell_term_list[i] = cell_terms[0]
	df_annot['cell_type_category'] = cell_term_list
	df_annot['cell_type'] = df_annot['cell_type'].astype(str)
	df_annot['cell_type_category'] = df_annot['cell_type_category'].astype(str)
	if (len(df_annot.loc[df_annot['cell_type_category']=='']) > 1):
		print("WARNING, there are cells that do not have a cell slim, so will just use cell_type")
		df_annot.loc[df_annot['cell_type_category']=='', 'cell_type_category'] = df_annot.loc[df_annot['cell_type_category']=='', 'cell_type']
	return df_annot


# Quality check final anndata created for cxg, sync up gene identifiers if necessary
def quality_check(adata):
	if adata.obs.isnull().values.any():
		sys.exit("There is at least one 'NaN' value in the cxg anndata obs dataframe.")
	elif 'default_visualization' in adata.uns:
		if adata.uns['default_visualization'] not in adata.obs.values:
			sys.exit("The default_visualization field is not in the cxg anndata obs dataframe.")
	elif len(adata.var.index.tolist()) > len(adata.raw.var.index.tolist()):
		sys.exit("There are more genes in normalized genes than in raw matrix.")


# Return uniqued list separated by ', '
def clean_list(lst):
	lst = lst.split(',')
	if '' in lst:
		lst.remove('')
	lst = list(set(lst))
	diseases_str = '[{}]'.format(', '.join(lst))
	return diseases_str


# Remove unused disease fields, and summarize with 'disease', 'disease_ontology_term_id', and 'reported_diseases'
# List in pandas are still considered strings, so need to join and split to create a list of all diseases
def report_diseases(mxr_df, exp_disease):
	if not exp_disease:
		mxr_df['disease'] = ['normal'] * len(mxr_df.index)
		mxr_df['disease_ontology_term_id'] = ['PATO:0000461'] * len(mxr_df.index)
	elif exp_disease['term_id'] in ','.join(mxr_df['sample_diseases_term_id'].unique()).split(','):
		if len(','.join(mxr_df['sample_diseases_term_id'].unique()).split(',')) <= 2:
			mxr_df['disease'] = mxr_df['sample_diseases_term_name'].str.contains(exp_disease['term_name']).astype('string')
			mxr_df['disease'] = mxr_df['disease'].str.replace('^False$', 'normal', regex=True)
			mxr_df['disease'] = mxr_df['disease'].str.replace('^True$', exp_disease['term_name'], regex=True)
			mxr_df['disease_ontology_term_id'] = mxr_df['sample_diseases_term_id'].str.contains(exp_disease['term_id']).astype('string')
			mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].str.replace('^False$', 'PATO:0000461', regex=True)
			mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].str.replace('^True$', exp_disease['term_id'], regex=True)
		else:
			sys.exit("There is unexpected extra disease states in biosamples: {}".format(mxr_df['sample_diseases_term_id'].unique()))
	elif exp_disease['term_id'] in ','.join(mxr_df['donor_diseases_term_id'].unique()).split(','):
		mxr_df['disease'] = mxr_df['donor_diseases_term_name'].str.contains(exp_disease['term_name']).astype('string')
		mxr_df['disease'] = mxr_df['disease'].str.replace('^False$', 'normal', regex=True)
		mxr_df['disease'] = mxr_df['disease'].str.replace('^True$', exp_disease['term_name'], regex=True)
		mxr_df['disease_ontology_term_id'] = mxr_df['donor_diseases_term_id'].str.contains(exp_disease['term_id']).astype('string')
		mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].str.replace('^False$', 'PATO:0000461', regex=True)
		mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].str.replace('^True$', exp_disease['term_id'], regex=True)
	else:
		sys.exit("Cannot find the experimental_variable_disease in donor or sample diseases: {}; {}; {}".format(exp_disease['term_id'], \
			','.join(mxr_df['sample_diseases_term_id'].unique()).split(','), ','.join(mxr_df['donor_diseases_term_id'].unique()).split(',')))
	# 'reported_diseases' is a list of all unique diseases from donor and samples
	mxr_df['reported_diseases'] = mxr_df['sample_diseases_term_name'] + ',' + mxr_df['donor_diseases_term_name']
	mxr_df['reported_diseases'] = mxr_df['reported_diseases'].apply(clean_list)
	#return mxr_df


# Demultiplex experimental metadata by finding demultiplexed suspension
# Determine overlapping suspension, create library & demultiplexed suspension df
# get cell_metadata from that suspension, merge in library info
# merge with mxr_df on library
def demultiplex(lib_donor_df, library_susp, donor_susp, mfinal_obj):
	susp_df = pd.DataFrame()
	lattice_donor = {}
	lattice_donor_col = []
	demult_susp_lst = []

	for donor_map in mfinal_obj['donor_mappings']:
		lattice_donor[donor_map['label']] = donor_map['donor']
	for author_don in lib_donor_df['author_donor'].to_list():
		lattice_donor_col.append(lattice_donor[author_don])
	lib_donor_df['author_donor_@id'] = lattice_donor_col
	lib_donor_df['library_donor_@id'] = lib_donor_df['library_@id'] + "," + lib_donor_df['author_donor_@id']

	for lib_donor_unique in lib_donor_df['library_donor_@id'].to_list():
		demult_susp = ''
		lib_uniq = lib_donor_unique.split(',')[0]
		donor_uniq = lib_donor_unique.split(',')[1]
		for susp in donor_susp[donor_uniq]:
			if susp in library_susp[lib_uniq]:
				demult_susp = susp
		demult_susp_lst.append(demult_susp)
		if demult_susp == '':
	 		print('ERROR: Could not find suspension for demultiplexed donor: {}, {}, {}'.format(donor, donor_susp[donor], library_susp[assoc_lib]))
	lib_donor_df['suspension_@id'] = demult_susp_lst

	obj_type_subset = ['sample', 'suspension', 'donor']
	for susp in lib_donor_df['suspension_@id'].to_list():
		values_to_add = {}
		susp_obj = lattice.get_object(susp, connection)
		relevant_objects = gather_objects(susp_obj, start_type='suspension')
		for obj_type in obj_type_subset:
		 	objs = relevant_objects.get(obj_type, [])
		 	if len(objs) == 1:
		 		gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, objs)
		 	else: 
		 		print('ERROR: Could not find suspension for demultiplexed donor: {}'.format(obj_type))
		row_to_add = pd.Series(values_to_add)
		susp_df = susp_df.append(row_to_add, ignore_index=True)
	lib_donor_df = lib_donor_df.merge(susp_df, left_on='suspension_@id', right_on='suspension_@id', how='left')
	return(lib_donor_df)

def main(mfinal_id):
	mfinal_obj = lattice.get_object(mfinal_id, connection)

	# confirm that the identifier you've provided corresponds to a MatrixFile
	mfinal_type = mfinal_obj['@type'][0]
	assay10x = ''
	if mfinal_type != 'MatrixFile':
		sys.exit('{} is not a MatrixFile, but a {}'.format(mfinal_id, mfinal_type))
	if mfinal_obj['assays'] == ['snATAC-seq']:
		assay10x = 'ATAC'
	elif mfinal_obj['assays'] == ['snRNA-seq'] or mfinal_obj['assays'] == ['scRNA-seq'] or\
			mfinal_obj['assays'] == ['snRNA-seq', 'scRNA-seq']:
		assay10x = 'RNA'
	else:
		sys.exit("Unexpected assay types to generate cxg h5ad: {}".format(mfinal_obj['assays']))

	# set the metadata keys based on defined metadata fields
	headers = []
	for obj_type in cell_metadata.keys():
		for prop in cell_metadata[obj_type]:
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = prop_map.get(latkey, latkey)
			headers.append(key)

	# Dataframe that contains experimental metadata keyed off of raw matrix
	df = pd.DataFrame()

	results = {}
	tmp_dir = 'matrix_files'
	os.mkdir(tmp_dir)
	download_file(mfinal_obj, tmp_dir)

	# Get list of unique final cell identifiers
	file_url = mfinal_obj['s3_uri']
	file_ext = file_url.split('.')[-1]
	mfinal_local_path = '{}/{}.{}'.format(tmp_dir, mfinal_obj['accession'], file_ext)
	mfinal_adata = None
	assays = []
	converted_h5ad = []
	for layer in mfinal_obj['layers']:
		if 'assay' in layer:
			assays.append(layer['assay'])
	if mfinal_obj['file_format'] == 'hdf5' and re.search('h5ad$', mfinal_local_path):
		mfinal_adata = sc.read_h5ad(mfinal_local_path)
	elif mfinal_obj['file_format'] == 'rds':
		converted_h5ad = convert_from_rds(mfinal_local_path, assays, tmp_dir, mfinal_obj['author_cell_type_column'])
		mfinal_adata = sc.read_h5ad(converted_h5ad[0][0])
	else:
		sys.exit('Do not recognize file format or exention {} {}'.format(mfinal_obj['file_format'], mfinal_local_path))
	mfinal_cell_identifiers = mfinal_adata.obs.index.to_list()

	cxg_adata_lst = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	donor_susp = {}
	library_susp = {}
	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		relevant_objects = gather_objects(mxr)
		values_to_add = {}

		# If there is a author_donor_column, assume it is a demuxlet experiment and demultiplex df metadata
		# Gather library, suspension, and donor associations while iterating through relevant objects
		# Cannot handle multiple pooling events, so will sys.exit
		if 'author_donor_column' in mfinal_obj:
			lib_obj = relevant_objects.get('library', [])
			gather_metdata('library', cell_metadata['library'], values_to_add, lib_obj)
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
						sys.exit('Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
					if single_susp['donors'][0] not in donor_susp:
						donor_susp[single_susp['donors'][0]] = [single_susp['@id']]
					else:
						donor_susp[single_susp['donors'][0]].append(single_susp['@id'])
			else:
				if len(single_susp['donors']) > 1:
					sys.exit('Not currently able to handle 2 pooling events: {}, {}'.format(single_susp['@id'], single_susp['donors']))
				if single_susp['donors'][0] not in donor_susp:
					donor_susp[susp_obj['donors'][0]] = [susp_obj['@id']]
				else:
					donor_susp[susp_obj['donors'][0]].append(susp_obj['@id'])

		# Gather metdata without demultiplexing
		else:
			for obj_type in cell_metadata.keys():
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, objs)
				elif len(objs) > 1:
					gather_pooled_metadata(obj_type, cell_metadata[obj_type], values_to_add, objs)
		row_to_add = pd.Series(values_to_add, name=mxr['@id'], dtype=str)
		df = df.append(row_to_add)
		
		# Add anndata to list of final raw anndatas, only for RNAseq
		if assay10x == 'RNA':
			download_file(mxr, tmp_dir)
			local_path = '{}/{}.h5'.format(tmp_dir,mxr_acc)
			adata_raw = sc.read_10x_h5(local_path)
			adata_raw.var_names_make_unique(join = '.')
			# Recreate cell_ids and subset raw matrix and add mxr_acc into obs
			concatenated_ids = concatenate_cell_id(mfinal_obj, mxr['@id'], adata_raw.obs_names, mfinal_cell_identifiers)
			adata_raw.obs_names = concatenated_ids
			overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
			adata_raw = adata_raw[overlapped_ids]
			adata_raw.obs['raw_matrix_accession'] = [mxr['@id']]*len(overlapped_ids)
			cxg_adata_lst.append(adata_raw)

	# get dataset-level metadata
	ds_results = report_dataset(relevant_objects['donor'], mfinal_obj, mfinal_obj['dataset'])

	# Should add error checking to make sure all matrices have the same number of vars
	feature_lengths = []
	for adata in cxg_adata_lst:
		feature_lengths.append(adata.shape[1])
	if len(set(feature_lengths)) > 1:
		sys.exit('The number of genes in all raw matrices need to match.')

	# Set up dataframe for cell annotations keyed off of author_cell_type
	annot_df = pd.DataFrame()
	for annot_obj in mfinal_obj['cell_annotations']:
		annot_lst = []
		annot_lst.append(annot_obj)
		annot_metadata = {}
		gather_metdata('cell_annotation', annot_fields, annot_metadata, annot_lst)
		annot_row = pd.Series(annot_metadata, name=annot_obj['author_cell_type'])
		annot_df = annot_df.append(annot_row)
	annot_df = trim_cell_slims(annot_df)

	# For RNA datasets, concatenate all anndata objects in list,
	# For ATAC datasets, assumption is that there is no scale.data, and raw count is taken from mfinal_adata.raw.X
	cell_mapping_rev_dct = {}
	raw_matrix_mapping = []
	for mapping_dict in mfinal_obj['cell_label_mappings']:
		cell_mapping_rev_dct[mapping_dict['label']] = mapping_dict['raw_matrix']
	if assay10x == 'RNA':
		cxg_adata_raw = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None)
		cxg_adata_raw = cxg_adata_raw[mfinal_cell_identifiers]
		if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
			sys.exit('The number of cells do not match between final matrix and cxg h5ad.')
	else:
		flag_removed = False
		for final_id in mfinal_cell_identifiers:
			if not re.search('[AGCT]+-1', final_id):
				flag_removed = True
		for cell_id in mfinal_cell_identifiers:
			if mfinal_obj['cell_label_location'] == 'prefix':
				label = re.search('^(.*)[AGCT]{16}.*$', cell_id).group(1)
			elif mfinal_obj['cell_label_location'] == 'suffix':
				if flag_removed:
					label = re.search(r'^[AGCT]+(.*)$', cell_id).group(1)
				else:
					label = re.search(r'^[AGCT]+-1(.*)$', cell_id).group(1)
			raw_matrix_mapping.append(cell_mapping_rev_dct[label])
		atac_obs = pd.DataFrame({'raw_matrix_accession': raw_matrix_mapping}, index = mfinal_cell_identifiers)
		cxg_adata_raw = ad.AnnData(mfinal_adata.raw.X, var = mfinal_adata.var, obs = atac_obs)

	# Set uns and obsm parameters
	cxg_uns = ds_results
	cxg_uns['version'] = {}
	cxg_uns['version']['corpora_schema_version'] = corpora_schema_version
	cxg_uns['version']['corpora_encoding_version'] = corpora_encoding_version
	cxg_obsm = get_embeddings(mfinal_adata)

	# Merge df with raw_obs according to raw_matrix_accession, and add additional cell metadata from mfinal_adata if available
	celltype_col = mfinal_obj['author_cell_type_column']
	cxg_obs = pd.merge(cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[celltype_col]], left_index=True, right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, annot_df, left_on=celltype_col, right_index=True, how='left')

	if 'author_cluster_column' in mfinal_obj:
		cluster_col = mfinal_obj['author_cluster_column']
		cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[cluster_col]], left_index=True, right_index=True, how='left')
		cxg_obs.rename(columns={cluster_col: 'author_cluster'}, inplace=True)
		cxg_obs['author_cluster'] = cxg_obs['author_cluster'].astype('category')

	# After getting experimental metadata keyed off of mxr, if there is author_donor_column, run demultiplex
	if 'author_donor_column' in mfinal_obj:
		donor_col = mfinal_obj['author_donor_column']
		cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[donor_col]], left_index=True, right_index=True, how='left')
		cxg_obs.rename(columns={donor_col: 'author_donor'}, inplace=True)
		cxg_obs['library_@id'] = cxg_obs['library_@id'].astype(str)
		cxg_obs['author_donor'] = cxg_obs['author_donor'].astype(str)
		cxg_obs['library_authordonor'] = cxg_obs['library_@id'] + ',' + cxg_obs['author_donor']

		lib_donor_df = cxg_obs[['library_@id', 'author_donor', 'library_authordonor']].drop_duplicates().reset_index(drop=True)
		donor_df = demultiplex(lib_donor_df, library_susp, donor_susp, mfinal_obj)
		report_diseases(donor_df, mfinal_obj.get('experimental_variable_disease'))
		# Retain cell identifiers as index
		cxg_obs = cxg_obs.reset_index().merge(donor_df, how='left', on='library_authordonor').set_index('index')
	else:
		# Go through donor and biosample diseases and calculate cxg field accordingly
		report_diseases(df, mfinal_obj.get('experimental_variable_disease'))
		cxg_obs = pd.merge(cxg_obs, df[['disease', 'disease_ontology_term_id', 'reported_diseases']], left_on="raw_matrix_accession", right_index=True, how="left" )

	# Drop columns that were used as intermediate calculations
	# Also check to see if optional columns are all empty, then drop those columns as well
	columns_to_drop = ['raw_matrix_accession', celltype_col, 'sample_diseases_term_id', 'sample_diseases_term_name',\
			'donor_diseases_term_id', 'donor_diseases_term_name', 'batch', 'library_@id_x', 'library_@id_y', 'author_donor_x',\
			'author_donor_y', 'library_authordonor', 'author_donor_@id', 'library_donor_@id', 'suspension_@id', 'library_@id']
	for column_drop in  columns_to_drop: 
		if column_drop in cxg_obs.columns.to_list():
			cxg_obs.drop(columns=column_drop, inplace=True)
	optional_columns = ['donor_BMI', 'family_history_breast_cancer', 'reported_diseases', 'donor_times_pregnant']
	for col in optional_columns:
		if col in cxg_obs.columns.to_list():
			col_content = cxg_obs[col].unique()
			if len(col_content) == 1:
				if col_content == ['[]'] or col_content == unreported_value:
					cxg_obs.drop(columns=col, inplace=True)

	if cxg_uns['organism'] == 'Homo sapiens':
		cxg_obs['ethnicity'] = cxg_obs['ethnicity'].str.replace('^$', 'unknown', regex=True)
	else:
		cxg_obs['ethnicity'] = cxg_obs['ethnicity'].str.replace('^$', 'na', regex=True)

	# Make sure gene ids match before using mfinal_data.var for cxg_adata
	for gene in list(mfinal_adata.var_names):
		if gene not in list(cxg_adata_raw.var_names):
			if re.search(r'^[A-Za-z]\S+-[0-9]$', gene):
				modified_gene = re.sub(r'(^[A-Z]\S+)-([0-9])$', r'\1.\2', gene)
				mfinal_adata.var.rename(index={gene: modified_gene}, inplace=True)
			else:
				sys.exit('There is a genes in the final matrix that is not in the raw matrix: {}'.format(gene))

	# If final matrix file is h5ad, take expression matrix from .X to create cxg anndata
	if converted_h5ad == []:
		cxg_var = cxg_adata_raw.var.loc[list(mfinal_adata.var_names),]
		cxg_adata = ad.AnnData(mfinal_adata.X, obs=cxg_obs, obsm=cxg_obsm, var=cxg_var, uns=cxg_uns)
		cxg_adata.raw = cxg_adata_raw
		quality_check(cxg_adata)
		results_file = '{}_v{}.h5ad'.format(mfinal_obj['accession'], flat_version)
		cxg_adata.write(results_file)
	else:
		# For seurat objects, create an anndata object for each assay; append assay name if more than 1 file
		for i in range(len(converted_h5ad)):
			if  converted_h5ad[i][1] == 'X':
				if i != 0:
					mfinal_adata = sc.read_h5ad(converted_h5ad[i][0])
				matrix_loc = mfinal_adata.X
				final_var = cxg_adata_raw.var.loc[list(mfinal_adata.var_names),]
			else:
				if i != 0:
					mfinal_adata = sc.read_h5ad(converted_h5ad[i][0])
				matrix_loc = mfinal_adata.raw.X
				final_var = cxg_adata_raw.var.loc[list(mfinal_adata.raw.var.index),]
			cxg_adata = ad.AnnData(matrix_loc, obs=cxg_obs, obsm=cxg_obsm, var=final_var, uns=cxg_uns)
			cxg_adata.raw = cxg_adata_raw
			quality_check(cxg_adata)
			if len(converted_h5ad) == 1:
				results_file = '{}_v{}.h5ad'.format(mfinal_obj['accession'], flat_version)
			else:
				results_file = '{}_{}_v{}.h5ad'.format(mfinal_obj['accession'], converted_h5ad[i][2], flat_version)
			cxg_adata.write(results_file)

	shutil.rmtree(tmp_dir)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
