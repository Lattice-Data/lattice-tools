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
from urllib.parse import urljoin
import requests
import numpy as np
import collections
import logging
import gc
from scipy import sparse

# Reference files by which the flattener will filter var features
ref_files = {
	'ercc':'genes_ercc.csv',
	'human':'genes_homo_sapiens.csv',
	'mouse':'genes_mus_musculus.csv',
	'sars':'genes_sars_cov_2.csv'
}

# Metadata to be gathered for each object type
cell_metadata = {
	'donor': [
		'donor_id',
		'age_display',
		'sex',
		'ethnicity.term_id',
		'diseases.term_id',
		'diseases.term_name',
		'body_mass_index',
		'times_pregnant',
		'family_history_breast_cancer',
		'organism.taxon_id',
		'risk_score_tyrer_cuzick_lifetime'
		],
	'sample': [
		'age_development_stage_redundancy',
		'uuid',
		'preservation_method',
		'biosample_ontology.term_id',
		'biosample_ontology.organ_slims',
		'biosample_ontology.cell_slims',
		'summary_development_ontology_at_collection.development_slims',
		'summary_development_ontology_at_collection.term_id',
		'diseases.term_id',
		'diseases.term_name',
		'disease_state',
		'treatment_summary'
		],
	'tissue_section': [
		'uuid',
		'thickness',
		'thickness_units'
	],
	'suspension': [
		'cell_depletion_factors',
		'enriched_cell_types.term_name',
		'enrichment_factors',
		'uuid',
		'suspension_type',
		'@id'
		],
	'library': [
		'uuid',
		'protocol.assay_ontology.term_id',
		'@id'
		]
	}

dataset_metadata = {
	'final_matrix': [
		'description',
		'default_embedding',
		'is_primary_data'
		]
	}

annot_fields = [
	'cell_ontology.term_id',
	'author_cell_type',
	'cell_state'
]

# Mapping of field name (object_type + "_" + property) and what needs to be in the final cxg h5ad
prop_map = {
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'sample_summary_development_ontology_at_collection_term_id': 'development_stage_ontology_term_id',
	'sample_age_development_stage_redundancy': 'donor_age_redundancy',
	'sample_disease_state': 'disease_state',
	'library_protocol_assay_ontology_term_id': 'assay_ontology_term_id',
	'donor_body_mass_index': 'donor_BMI',
	'donor_sex': 'sex',
	'donor_donor_id': 'donor_id',
	'donor_organism_taxon_id': 'organism_ontology_term_id',
	'donor_ethnicity_term_id': 'self_reported_ethnicity_ontology_term_id',
	'donor_age_display': 'donor_age',
	'donor_family_history_breast_cancer': 'family_history_breast_cancer',
	'donor_risk_score_tyrer_cuzick_lifetime': 'tyrer_cuzick_lifetime_risk',
	'matrix_description': 'title',
	'matrix_default_embedding': 'default_embedding',
	'matrix_is_primary_data': 'is_primary_data',
	'cell_annotation_author_cell_type': 'author_cell_type',
	'cell_annotation_cell_ontology_term_id': 'cell_type_ontology_term_id',
	'cell_annotation_cell_state': 'cell_state',
	'suspension_suspension_type': 'suspension_type',
	'suspension_enriched_cell_types_term_name': 'suspension_enriched_cell_types',
	'suspension_cell_depletion_factors': 'suspension_depletion_factors'
}

# Global variables
unreported_value = 'unknown'
schema_version = '3.0.0'
flat_version = '5'
tmp_dir = 'matrix_files'
mfinal_obj = None
mfinal_adata = None
cxg_adata = None
cxg_adata_raw = None
cxg_obs = None

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


# Gather raw matrices by object type and 'background_barcodes_included' to select for filtered matrix from CR output
def gather_rawmatrices(derived_from):
	my_raw_matrices = []
	df_ids = []
	for identifier in derived_from:
		obj = lattice.get_object(identifier, connection)
		if obj['@type'][0] == 'RawMatrixFile':
			my_raw_matrices.append(obj)
		else:
			# grab the derived_from in case we need to go a layer deeper
			for i in obj['derived_from']:
				df_ids.append(i)
	if not my_raw_matrices:
		for identifier in df_ids:
			obj = lattice.get_object(identifier, connection)
			if obj['@type'][0] == 'RawMatrixFile':
				my_raw_matrices.append(obj)
	return my_raw_matrices


# Gather all objects up the experimental graph, assuming that there is either a suspension or tissue section
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
	tissue_section_ids = []
	tissue_sections = []

	if start_type == None:
		for i in lib_ids:
			obj = lattice.get_object(i, connection)
			libraries.append(obj)
			for o in obj['derived_from']:
				if o.get('uuid') not in susp_ids:
					if 'Suspension' in o['@type']:
						suspensions.append(o)
						susp_ids.append(o.get('uuid'))
					elif 'TissueSection' in o['@type']:
						tissue_sections.append(o)
						tissue_section_ids.append(o.get('uuid'))

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

	if len(suspensions) > 0:
		for o in suspensions:
			for i in o['derived_from']:
				sample_ids.append(i)
	else:
		for o in tissue_sections:
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
		'suspension': suspensions,
		'tissue_section': tissue_sections
		}
	if start_type == None:
		objs['library'] = libraries
	if prepooled_susps:
		objs['prepooled_suspension'] = prepooled_susps
		objs['pooled_suspension'] = objs['suspension']

	return objs


# Get property value for given object, can only traverse embedded objects that are embedded 2 levels in
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
	elif len(path) == 3:
		key1 = path[0]
		key2 = path[1]
		key3 = path[2]
		if isinstance(obj.get(key1), list):
			embed_objs = obj.get(key1, unreported_value)
			values = []
			for embed_obj in embed_objs:
				if isinstance(embed_obj.get(key2), list):
					values += [k.get(key3, unreported_value) for k in embed_obj[key2]]
				else:
					values += embed_obj[key2].get(key3, unreported_value)
			return list(set(values))
		# Will need to revisit cell culture and organoid values 
		elif obj.get(key1):
			embed_obj = obj.get(key1, unreported_value)
			if isinstance(embed_obj.get(key2, unreported_value), list):
				return [v.get(key3, unreported_value) for v in embed_obj[key2]]
			elif embed_obj.get(key2, unreported_value) == unreported_value:
				return unreported_value
			else:
				return embed_obj[key2].get(key3, unreported_value)
		else:
			return obj.get(key1, unreported_value)
	else:
		return 'unable to traverse more than 2 embeddings'


# Gather object metadata, convert property name to cxg required field names
def gather_metdata(obj_type, properties, values_to_add, objs):
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = value


# Gather metadata for pooled objects
# For required cxg fields, these need to be a single value, so development_stage_ontology_term_id needs to be a commnon slim
def gather_pooled_metadata(obj_type, properties, values_to_add, objs):
	dev_list = []
	for prop in properties:
		value = list()
		for obj in objs:
			v = get_value(obj, prop)
			if prop == 'summary_development_ontology_at_collection.development_slims':
				dev_list.append(v)
			if isinstance(v, list):
				value.extend(v)
			else:
				value.append(v)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		value_str = [str(i) for i in value]
		value_set = set(value_str)
		cxg_fields = ['disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id',\
						 'sex', 'tissue_ontology_term_id', 'development_stage_ontology_term_id']
		if len(value_set) > 1:
			if key in cxg_fields:
				if key == 'development_stage_ontology_term_id':
					dev_in_all = list(set.intersection(*map(set, dev_list)))
					if dev_in_all == []:
						sys.exit("There is no common development_slims that can be used for development_stage_ontology_term_id")
					else:
						query_url = urljoin(server, 'search/?type=OntologyTerm&term_name=' + dev_in_all[0] + '&format=json')
						r = requests.get(query_url, auth=connection.auth)
						try:
							r.raise_for_status()
						except requests.HTTPError:
							sys.exit("Error in getting development_slims as development_stage ontology: {}".format(query_url))
						else:
							if r.json()['total']==1:
								values_to_add[key] = r.json()['@graph'][0]['term_id']
							else:
								sys.exit("Error in getting development_slims as development_stage ontology: {}".format(query_url))
				elif key == 'sex':
					values_to_add[key] = 'unknown'
				else:
					sys.exit("Cxg field is a list")
			else:		
				values_to_add[key] = 'pooled [{}]'.format(','.join(value_str))
		else:
			values_to_add[key] = next(iter(value_set))


# Gather dataset metadata for adata.uns
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
	return ds_results


# Download file object from s3
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


# Compile all reference annotations for var features into one pandas df
def compile_annotations(files):
	ids = pd.DataFrame()
	client = boto3.client('s3')
	bucket_name = 'submissions-lattice'
	for key in files:
		filename = tmp_dir + "/" + files[key]
		try:
			client.download_file(bucket_name, 'cxg_migration/var_refs/' + files[key], filename)
		except subprocess.CalledProcessError as e:
			sys.exit('Failed to find file {} on s3'.format(file_obj.get('@id')))
		else:
			print("Downloading reference: {}".format(files[key]))
		df = pd.read_csv(filename, names=['feature_id','symbol','num'], dtype='str')
		ids  = pd.concat([ids,df])
	return ids


# Recreated the final matrix ids, also checking to see if '-1' was removed from original cell identifier
def concatenate_cell_id(mxr_acc, raw_obs_names, mfinal_cells):
	global mfinal_obj
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

# Quality check final anndata created for cxg, sync up gene identifiers if necessary
def quality_check(adata):
	if adata.obs.isnull().values.any():
		print("WARNING: There is at least one 'NaN' value in the cxg anndata obs dataframe.")
	elif 'default_visualization' in adata.uns:
		if adata.uns['default_visualization'] not in adata.obs.values:
			sys.exit("The default_visualization field is not in the cxg anndata obs dataframe.")
	elif len(adata.var.index.tolist()) > len(adata.raw.var.index.tolist()):
		sys.exit("There are more genes in normalized genes than in raw matrix.")


# Return value to be stored in disease field based on list of diseases from donor and sample
def clean_list(lst, exp_disease):
	lst = lst.split(',')
	disease = exp_disease['term_id'] if exp_disease['term_id'] in lst else 'PATO:0000461'
	return disease


# Determine reported disease as unique of sample and donor diseases, removing unreported value
def report_diseases(mxr_df, exp_disease):
	mxr_df['reported_diseases'] = mxr_df[['sample_diseases_term_name','donor_diseases_term_name']].stack().groupby(level=0).apply(lambda x: [i for i in x.unique() if i != unreported_value])
	mxr_df['reported_diseases'] = mxr_df['reported_diseases'].astype(dtype='string')
	mxr_df['reported_diseases'] = mxr_df['reported_diseases'].apply(lambda x: x.replace("'", ""))
	total_reported = mxr_df['reported_diseases'].unique()
	if len(total_reported) == 1:
		if total_reported[0] == '[]':
			mxr_df['reported_diseases'] = '[]'
	elif '[]' in total_reported:
		mxr_df['reported_diseases'].replace({'[]':'none'}, inplace=True)

	if exp_disease == unreported_value:
		mxr_df['disease_ontology_term_id'] = ['PATO:0000461'] * len(mxr_df.index)
	else:
		mxr_df['disease_ontology_term_id'] = mxr_df['sample_diseases_term_id'] + ',' + mxr_df['donor_diseases_term_id']
		mxr_df['disease_ontology_term_id'] = mxr_df['disease_ontology_term_id'].apply(clean_list, exp_disease=exp_disease)
		exp_disease_aslist = '[{}]'.format(exp_disease['term_name'])
		if len([x for x in total_reported if x not in ['none', exp_disease_aslist,'[]']])==0:
			mxr_df['reported_diseases'] = '[]'

# Demultiplex experimental metadata by finding demultiplexed suspension 
# Determine overlapping suspension, create library & demultiplexed suspension df
# get cell_metadata from that suspension, merge in library info
# merge with mxr_df on library
def demultiplex(lib_donor_df, library_susp, donor_susp):
	global mfinal_obj
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
		if demult_susp == '':
			print('ERROR: Could not find suspension for demultiplexed donor: {}, {}, {}'.format(donor, donor_susp[donor], library_susp[assoc_lib]))
		else:
			demult_susp_lst.append(demult_susp)
	lib_donor_df['suspension_@id'] = demult_susp_lst

	obj_type_subset = ['sample', 'suspension', 'donor']
	for susp in set(lib_donor_df['suspension_@id'].to_list()):
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


# For cell culture, tissue is not UBERON, use cell slims to get CL
def get_cell_slim(df_series, suffix):
	cell = df_series['sample_biosample_ontology_cell_slims'].split("'")[1].replace(" ", "+")
	df_series.drop(labels='sample_biosample_ontology_cell_slims', inplace=True)
	query_url = urljoin(server, 'search/?type=OntologyTerm&term_name=' + cell + '&format=json')
	r = requests.get(query_url, auth=connection.auth)
	try:
		r.raise_for_status()
	except requests.HTTPError:
		sys.exit("Error in getting cell slim as tissue ontology: {}".format(query_url))
	else:
		if r.json()['total']==1:
			df_series['tissue_ontology_term_id'] = r.json()['@graph'][0]['term_id'] + suffix
		else:
			sys.exit("Error in getting organ slim as tissue ontology: {}".format(query_url))


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
			sys.exit("Unexpected sex: {}".format(sex))


# Make sure cxg_adata and cxg_adata_raw have same number of features
# If not, add implied zeros to csr, and add corresponding 'feature_is_filtered'
def add_zero():
	global cxg_adata
	global cxg_adata_raw
	if cxg_adata_raw.shape[1] > cxg_adata.shape[1]:
		genes_add = [x for x in cxg_adata_raw.var.index.to_list() if x not in cxg_adata.var.index.to_list()]
		new_matrix = sparse.csr_matrix((cxg_adata.X.data, cxg_adata.X.indices, cxg_adata.X.indptr), shape = cxg_adata_raw.shape)
		all_genes = cxg_adata.var.index.to_list()
		all_genes.extend(genes_add)
		new_var = pd.DataFrame(index=all_genes)
		new_var = pd.merge(new_var, cxg_adata_raw.var, left_index=True, right_index=True, how='left')
		new_var['feature_is_filtered'] = False
		new_var.loc[genes_add, 'feature_is_filtered'] = True
		new_adata = ad.AnnData(X=new_matrix, obs=cxg_adata.obs, var=new_var, uns=cxg_adata.uns, obsm=cxg_adata.obsm)
		new_adata = new_adata[:,cxg_adata_raw.var.index.to_list()]
		new_adata.var = new_adata.var.merge(cxg_adata.var, left_index=True, right_index=True, how='left')
		cxg_adata = new_adata
	else:
		cxg_adata.var['feature_is_filtered'] = False


# Use cxg_adata_raw var to map ensembl IDs and use that as index and drop redundants
# Make sure the indices are the same order for both anndata objects & clean up var metadata
# WILL NEED TO ADD NEW BIOTYPE FOR CITE-SEQ
def set_ensembl(redundant, feature_keys):
	global cxg_adata
	global cxg_adata_raw
	raw_cols = cxg_adata_raw.var.columns.to_list()
	cxg_adata_raw.var.drop(columns=[i for i in raw_cols if i!='gene_ids'], inplace=True)
	if feature_keys == ['gene symbol']:
		if 'gene_ids' in cxg_adata_raw.var.columns.to_list():
			# Drop unmapped genes from normalized matrix
			norm_index = set(cxg_adata.var.index.to_list())
			raw_index = set(cxg_adata_raw.var.index.to_list())
			drop_unmapped = list(norm_index.difference(raw_index))
			cxg_adata = cxg_adata[:, [i for i in cxg_adata.var.index.to_list() if i not in drop_unmapped]]
			logging.info('drop_unmapped:\t{}\t{}'.format(len(drop_unmapped),drop_unmapped))

			cxg_adata.var = pd.merge(cxg_adata.var, cxg_adata_raw.var, left_index=True, right_index=True, how='left', copy = True)
			cxg_adata.var = cxg_adata.var.set_index('gene_ids', drop=True)
			cxg_adata_raw.var  = cxg_adata_raw.var.set_index('gene_ids', drop=True)
			cxg_adata.var.index.name = None
			cxg_adata_raw.var.index.name = None

			# Drop redundant by Ensembl ID
			drop_redundant = list(set(redundant).intersection(set(cxg_adata.var.index.to_list())))
			logging.info('drop_redundant:\t{}\t{}'.format(len(drop_redundant),drop_redundant))
			cxg_adata = cxg_adata[:, [i for i in cxg_adata.var.index.to_list() if i not in redundant]]

		else:
			print("WARNING: raw matrix does not have genes_ids column")
	elif feature_keys == ['Ensembl gene ID']:
		cxg_adata_raw.var_names_make_unique()
		cxg_adata_raw.var  = cxg_adata_raw.var.set_index('gene_ids', drop=True)
		cxg_adata_raw.var.index.name = None
		unique_to_norm =  set(cxg_adata.var.index.to_list()).difference(set(cxg_adata_raw.var.index.to_list()))
		if len(unique_to_norm) > 0:
			print("WARNING: normalized matrix contains Ensembl IDs not in raw: {}".format(unique_to_norm))


def filter_ensembl(adata, compiled_annot):
	var_in_approved = adata.var.index[adata.var.index.isin(compiled_annot['feature_id'])]
	adata = adata[:, var_in_approved]
	return adata


# Clean up columns and column names for var, gene name must be gene_id
# WILL NEED TO REVISIT WHEN THERE IS MORE THAN ONE VALUE FOR GENE ID
def clean_var():
	global cxg_adata
	global cxg_adata_raw
	global mfinal_adata
	gene_pd = cxg_adata_raw.var[[i for i in cxg_adata_raw.var.columns.values.tolist() if 'gene_ids' in i]]
	gene_pd = gene_pd.replace('nan', np.nan)
	gene_pd = gene_pd.stack().groupby(level=0).apply(lambda x: x.unique()[0]).to_frame(name='gene_ids')
	cxg_adata_raw.var.drop(columns = cxg_adata_raw.var.columns.tolist(), inplace=True)
	cxg_adata_raw.var = cxg_adata_raw.var.merge(gene_pd, left_index = True, right_index=True, how = 'left')


# Reconcile genes if raw matrices annotated to multiple version by merging raw based on Ensembl ID
# Return raw matrix merged on Ensembl ID and list of gene to drop from normalized matrices
def reconcile_genes(cxg_adata_lst):
	global mfinal_obj
	global mfinal_adata
	mfinal_adata_genes = mfinal_adata.var.index.to_list()
	redundant = []
	stats = {}

	# Join raw matrices on ensembl, gene symbols stored as metadata
	for cxg_adata in cxg_adata_lst:
		cxg_adata.var['gene_symbols'] = cxg_adata.var.index
		cxg_adata.var = cxg_adata.var.set_index('gene_ids', drop=True)
	cxg_adata_raw_ensembl = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None, join='outer')

	# Join raw matrices on gene symbol, ensembl stored as metadata. Add suffix to make unique, using '.' as to match R default
	for cxg_adata in cxg_adata_lst:
		cxg_adata.var['gene_ids'] = cxg_adata.var.index
		cxg_adata.var =  cxg_adata.var.set_index('gene_symbols', drop=True)
		cxg_adata.var_names_make_unique(join = '.')
	cxg_adata_raw_symbol = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None, join='outer')

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
	        redundant.extend(gene_pd_ensembl[gene_pd_ensembl[col] == gene].index.to_list())
	redundant = list(set(redundant))
	stats['redundant'] = redundant

	# Clean up raw.var in outer join on ensembl and switch to gene symbol for index. Run var_names_make_unique and remove redundants after mapping of Ensembl
	cxg_adata_raw_ensembl.var['gene_ids'] = cxg_adata_raw_ensembl.var.index
	cxg_adata_raw_ensembl.var['gene_symbols'] = gene_pd_ensembl.stack().groupby(level=0).apply(lambda x: x.unique()[0]).to_frame(name='gene_symbols')
	cxg_adata_raw_ensembl.var = cxg_adata_raw_ensembl.var.set_index('gene_symbols', drop=True)
	cxg_adata_raw_ensembl.var.index.name = None
	cxg_adata_raw_ensembl.var_names_make_unique(join = '.')

	all_remove = list(set(stats['multiple_ensembl'] + stats['multiple_symbols']))

	for key in stats:
		stats[key] = set(stats[key])
		overlap_norm = set(mfinal_adata_genes).intersection(stats[key])
		logging.info("{}\t{}\t{}\t{}\t{}".format(key, len(stats[key]), len(overlap_norm), overlap_norm, stats[key]))

	return cxg_adata_raw_ensembl, redundant, all_remove


# filename will be collectionuuid_datasetuuid_accession_version.h5ad, collectionuuid_accession_version.h5ad, or accession_version.h5ad
# depending on what information is available for the dataset
def get_results_filename(mfinal_obj):
	results_file = None
	dataset = mfinal_obj.get('dataset',[])
	dataset_obj = lattice.get_object(dataset, connection)
	collection_id = None
	if dataset_obj.get('cellxgene_urls',[]):
		collection_id = dataset_obj.get('cellxgene_urls',[])[0]
		collection_id = collection_id.replace("https://cellxgene.cziscience.com/collections/","")
	if mfinal_obj.get('cellxgene_uuid',[]) and collection_id:
		results_file = '{}_{}_{}_v{}.h5ad'.format(collection_id, mfinal_obj['cellxgene_uuid'], mfinal_obj['accession'], flat_version)
	elif collection_id:
		results_file = '{}_{}_v{}.h5ad'.format(collection_id, mfinal_obj['accession'], flat_version)
	else:
		results_file = '{}_v{}.h5ad'.format(mfinal_obj['accession'], flat_version)
	return results_file


# Final touches for obs columns, dropping unnecessary column and modifying any Lattice fields to fit cxg schema
def clean_obs(celltype_col):
	global cxg_obs
	global mfinal_obj
	# For columns in mfinal_obj that contain continuous cell metrics, they are transferred to cxg_obs as float datatype
	# WILL NEED TO REVISIT IF FINAL MATRIX CONTAINS MULTIPLE LAYERS THAT WE ARE WRANGLING
	for author_col in mfinal_obj.get('author_columns',[]):
		if author_col in mfinal_adata.obs.columns.to_list():
			cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[author_col]], left_index=True, right_index=True, how='left')
		else:
			print("WARNING: author_column not in final matrix: {}".format(author_col))

	if 'NCIT:C17998' in cxg_obs['self_reported_ethnicity_ontology_term_id'].unique():
		cxg_obs.loc[cxg_obs['organism_ontology_term_id'] == 'NCBITaxon:9606', 'self_reported_ethnicity_ontology_term_id'] = cxg_obs['self_reported_ethnicity_ontology_term_id'].str.replace('NCIT:C17998', 'unknown')

	# if the donor has multiple ethnicities, self_reported_ethnicity_ontology_term_id is a list, set ontology term to multiethnic
	# need to complete test on this section.
	if len([i for i in cxg_obs['self_reported_ethnicity_ontology_term_id'].unique() if ',' in i]) > 0:
		for multi in [i for i in cxg_obs['self_reported_ethnicity_ontology_term_id'].unique() if ',' in i]:
			cxg_obs['self_reported_ethnicity_ontology_term_id'].replace({multi:'multiethnic'}, inplace=True)
	
	# if obs category suspension_type does not exist in dataset, create column and fill values with na (for spatial assay)
	if 'suspension_type' not in cxg_obs.columns:
		cxg_obs.insert(len(cxg_obs.columns),'suspension_type', 'na')
	elif cxg_obs['suspension_type'].isnull().values.any():
		cxg_obs['suspension_type'].fillna(value='na', inplace=True)
	
	# Drop columns that were used as intermediate calculations
	# Also check to see if optional columns are all empty, then drop those columns as well
	optional_columns = ['donor_BMI', 'family_history_breast_cancer', 'reported_diseases', 'donor_times_pregnant', 'sample_preservation_method',\
			'sample_treatment_summary', 'suspension_uuid', 'tissue_section_thickness', 'tissue_section_thickness_units','cell_state',\
			'suspension_enriched_cell_types', 'suspension_enrichment_factors', 'suspension_depletion_factors', 'tyrer_cuzick_lifetime_risk', 'disease_state']
	for col in optional_columns:
		if col in cxg_obs.columns.to_list():
			col_content = cxg_obs[col].unique()
			if len(col_content) == 1:
				if col_content[0] == unreported_value or col_content[0] == '[' + unreported_value + ']' or col_content[0] == '[]':
					cxg_obs.drop(columns=col, inplace=True)

	if len(cxg_obs['donor_age_redundancy'].unique()) == 1:
		if cxg_obs['donor_age_redundancy'].unique():
			cxg_obs.drop(columns='donor_age', inplace=True)
	columns_to_drop = ['raw_matrix_accession', celltype_col, 'sample_diseases_term_id', 'sample_diseases_term_name',\
			'donor_diseases_term_id', 'donor_diseases_term_name', 'batch', 'library_@id_x', 'library_@id_y', 'author_donor_x',\
			'author_donor_y', 'library_authordonor', 'author_donor_@id', 'library_donor_@id', 'suspension_@id', 'library_@id', 'sex',
			'sample_biosample_ontology_cell_slims', 'sample_summary_development_ontology_at_collection_development_slims','donor_age_redundancy',\
			'sample_biosample_ontology_organ_slims']
	for column_drop in  columns_to_drop: 
		if column_drop in cxg_obs.columns.to_list():
			cxg_obs.drop(columns=column_drop, inplace=True)
	if 'tissue_section_thickness' in cxg_obs.columns.to_list() and 'tissue_section_thickness_units' in cxg_obs.columns.to_list():
		cxg_obs['tissue_section_thickness'] = cxg_obs['tissue_section_thickness'].astype(str) + cxg_obs['tissue_section_thickness_units'].astype(str)
		cxg_obs.drop(columns='tissue_section_thickness_units', inplace=True)
	change_unreported = ['suspension_enriched_cell_types', 'suspension_enrichment_factors', 'suspension_depletion_factors', 'disease_state', 'cell_state']
	for field in change_unreported:
		if field in cxg_obs.columns.to_list():
			cxg_obs[field].replace({unreported_value: 'na'}, inplace=True)


def main(mfinal_id):
	global mfinal_obj
	global mfinal_adata
	global cxg_adata
	global cxg_adata_raw
	global cxg_obs
	mfinal_obj = lattice.get_object(mfinal_id, connection)
	logging.basicConfig(filename='outfile_flattener.log', level=logging.INFO)

	# confirm that the identifier you've provided corresponds to a ProcessedMatrixFile
	mfinal_type = mfinal_obj['@type'][0]
	summary_assay = ''
	if mfinal_type != 'ProcessedMatrixFile':
		sys.exit('{} is not a ProcessedMatrixFile, but a {}'.format(mfinal_id, mfinal_type))
	if mfinal_obj['assays'] == ['snATAC-seq']:
		summary_assay = 'ATAC'
	elif mfinal_obj['assays'] == ['snRNA-seq'] or mfinal_obj['assays'] == ['scRNA-seq'] or\
			mfinal_obj['assays'] == ['snRNA-seq', 'scRNA-seq'] or mfinal_obj['assays'] == ['spatial transcriptomics'] or\
			mfinal_obj['assays'] == ['scRNA-seq', 'snRNA-seq']:
		summary_assay = 'RNA'
	else:
		sys.exit("Unexpected assay types to generate cxg h5ad: {}".format(mfinal_obj['assays']))

	# Dataframe that contains experimental metadata keyed off of raw matrix
	df = pd.DataFrame()

	results = {}
	os.mkdir(tmp_dir)
	download_file(mfinal_obj, tmp_dir)

	# Get list of unique final cell identifiers
	file_url = mfinal_obj['s3_uri']
	file_ext = file_url.split('.')[-1]
	mfinal_local_path = '{}/{}.{}'.format(tmp_dir, mfinal_obj['accession'], file_ext)
	mfinal_adata = sc.read_h5ad(mfinal_local_path)
	mfinal_cell_identifiers = mfinal_adata.obs.index.to_list()
	if 'counts' in mfinal_adata.layers:
		del(mfinal_adata.layers['counts'])

	cxg_adata_lst = []
	redundant = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	donor_susp = {}
	library_susp = {}


	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		relevant_objects = gather_objects(mxr)
		values_to_add = {}

		# If there is a demultiplexed_donor_column, assume it is a demuxlet experiment and demultiplex df metadata
		# Gather library, suspension, and donor associations while iterating through relevant objects
		# Cannot handle multiple pooling events, so will sys.exit
		if 'demultiplexed_donor_column' in mfinal_obj:
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
		row_to_add = pd.DataFrame(values_to_add, index=[mxr['@id']], dtype=str)

		# make sure donor_df contains UBERON for tissue, may need to revisit 'if' statement
		if 'demultiplexed_donor_column' not in mfinal_obj:
			if not row_to_add.loc[mxr['@id'],'tissue_ontology_term_id'].startswith('UBERON'):
				if row_to_add.loc[mxr['@id'],'tissue_ontology_term_id'].endswith('(cell culture)'):
					get_cell_slim(row_to_add, ' (cell culture)')
				else:
					sys.exit('Tissue should have an UBERON ontology term: {}'.format(row_to_add['tissue_ontology_term_id']))
		
		# Add anndata to list of final raw anndatas, only for RNAseq
		if summary_assay == 'RNA':
			download_file(mxr, tmp_dir)
			row_to_add['mapped_reference_annotation'] = mxr['genome_annotation']
			if mxr['s3_uri'].endswith('h5'):
				local_path = '{}/{}.h5'.format(tmp_dir, mxr_acc)
				adata_raw = sc.read_10x_h5(local_path)
			elif mxr['s3_uri'].endswith('h5ad'):
				local_path = '{}/{}.h5ad'.format(tmp_dir, mxr_acc)
				adata_raw = sc.read_h5ad(local_path)
			else:
				sys.exit('Raw matrix file of unknown file extension: {}'.format(mxr['s3_uri']))
			# only make var unique if all raw matrices are same annotation version
			if len(mfinal_obj.get('genome_annotations', [])) == 1:
				for gene in [i for i, c in collections.Counter(adata_raw.var['gene_ids'].dropna().to_list()).items() if c > 1]:
					redundant.extend(adata_raw.var[adata_raw.var['gene_ids'] == gene].index.to_list())
				adata_raw.var_names_make_unique(join = '.')
			# Recreate cell_ids and subset raw matrix and add mxr_acc into obs
			if mfinal_obj.get('cell_label_mappings', None):
				concatenated_ids = concatenate_cell_id(mxr['@id'], adata_raw.obs_names, mfinal_cell_identifiers)
				adata_raw.obs_names = concatenated_ids
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
			else:
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(adata_raw.obs_names.to_list()))
			if len(overlapped_ids) == 0:
				if mfinal_obj['cell_label_location'] == 'prefix':
					if concatenated_ids[0].endswith('-1'):
						concatenated_ids = [re.sub('-1$', '', i) for i in concatenated_ids]
					else:
						concatenated_ids = [i+'-1' for i in concatenated_ids]
					overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
					if len(overlapped_ids) == 0:
						sys.exit("Could not find any matching cell identifiers: {}".format(concatenated_ids[0:5]))
					adata_raw.obs_names = concatenated_ids
				else:
					sys.exit("Could not find any matching cell identifiers: {}".format(concatenated_ids[0:5]))
			adata_raw = adata_raw[overlapped_ids]
			adata_raw.obs['raw_matrix_accession'] = mxr['@id']
			cxg_adata_lst.append(adata_raw)

		df = pd.concat([df, row_to_add])
		redundant = list(set(redundant))

	# get dataset-level metadata and set 'is_primary_data' for obs accordingly as boolean
	ds_results = report_dataset(relevant_objects['donor'], mfinal_obj, mfinal_obj['dataset'])
	df['is_primary_data'] = ds_results['is_primary_data']
	df['is_primary_data'].replace({'True': True, 'False': False}, inplace=True)

	del ds_results['is_primary_data']

	# Should add error checking to make sure all matrices have the same number of vars
	feature_lengths = []
	for adata in cxg_adata_lst:
		feature_lengths.append(adata.shape[1])
	feature_lengths = list(set(feature_lengths))

	# Set up dataframe for cell annotations keyed off of author_cell_type
	annot_df = pd.DataFrame()
	for annot_obj in mfinal_obj['cell_annotations']:
		annot_lst = []
		annot_lst.append(annot_obj)
		annot_metadata = {}
		gather_metdata('cell_annotation', annot_fields, annot_metadata, annot_lst)
		annot_row = pd.DataFrame(annot_metadata, index=[annot_obj['author_cell_type']])
		annot_df = pd.concat([annot_df, annot_row])

	# For RNA datasets, concatenate all anndata objects in list,
	# For ATAC datasets, assumption is that there is no scale.data, and raw count is taken from mfinal_adata.raw.X
	raw_matrix_mapping = []
	cell_mapping_rev_dct = {}
	if summary_assay == 'RNA':
		# If raw matrices are annotated to multiple gencode versions, concatenate on ensembl ID and remove ambiguous symbols
		if len(mfinal_obj.get('genome_annotations',[])) > 1:
			cxg_adata_raw, redundant, all_remove  = reconcile_genes(cxg_adata_lst)
			drop_removes = set(mfinal_adata.var.index.to_list()).intersection(set(all_remove))
			logging.info('drop_all_removes:\t{}\t{}'.format(len(drop_removes), drop_removes))
			mfinal_adata = mfinal_adata[:, [i for i in mfinal_adata.var.index.to_list() if i not in all_remove]]
		else:
			cxg_adata_raw = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None, join='outer')
			if len(feature_lengths) == 1:
				if cxg_adata_raw.var.shape[0] != feature_lengths[0]:
					sys.exit('There should be the same genes for raw matrices if only a single genome annotation')
		cxg_adata_raw = cxg_adata_raw[mfinal_cell_identifiers]
		if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
			sys.exit('The number of cells do not match between final matrix and cxg h5ad.')
	elif summary_assay == 'ATAC':
		for mapping_dict in mfinal_obj['cell_label_mappings']:
			cell_mapping_rev_dct[mapping_dict['label']] = mapping_dict['raw_matrix']
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
	cxg_uns['schema_version'] = schema_version
	cxg_obsm = mfinal_adata.obsm.copy()
	if len([i for i in cxg_obsm.keys() if i.startswith('X_')]) < 1:
		sys.exit("At least one embedding that starts with 'X_' is required")

	# Merge df with raw_obs according to raw_matrix_accession, and add additional cell metadata from mfinal_adata if available
	# Also add calculated fields to df 
	celltype_col = mfinal_obj['author_cell_type_column']
	if mfinal_adata.obs[celltype_col].dtype != 'object':
		mfinal_adata.obs[celltype_col] = mfinal_adata.obs[celltype_col].astype('string')
	cxg_obs = pd.merge(cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[celltype_col]], left_index=True, right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, annot_df, left_on=celltype_col, right_index=True, how='left')

	if 'author_cluster_column' in mfinal_obj:
		cluster_col = mfinal_obj['author_cluster_column']
		cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[cluster_col]], left_index=True, right_index=True, how='left')
		cxg_obs.rename(columns={cluster_col: 'author_cluster'}, inplace=True)
		cxg_obs['author_cluster'] = cxg_obs['author_cluster'].astype('category')

	# After getting experimental metadata keyed off of mxr, if there is demultiplexed_donor_column, run demultiplex
	if 'demultiplexed_donor_column' in mfinal_obj:
		donor_col = mfinal_obj['demultiplexed_donor_column']
		cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[donor_col]], left_index=True, right_index=True, how='left')
		cxg_obs.rename(columns={donor_col: 'author_donor'}, inplace=True)
		cxg_obs['library_@id'] = cxg_obs['library_@id'].astype(str)
		cxg_obs['author_donor'] = cxg_obs['author_donor'].astype(str)
		cxg_obs['library_authordonor'] = cxg_obs['library_@id'] + ',' + cxg_obs['author_donor']

		lib_donor_df = cxg_obs[['library_@id', 'author_donor', 'library_authordonor']].drop_duplicates().reset_index(drop=True)
		donor_df = demultiplex(lib_donor_df, library_susp, donor_susp)

		report_diseases(donor_df, mfinal_obj.get('experimental_variable_disease', unreported_value))
		get_sex_ontology(donor_df)

		# Retain cell identifiers as index
		cxg_obs = cxg_obs.reset_index().merge(donor_df, how='left', on='library_authordonor').set_index('index')
		if mfinal_adata.X.shape[0] != cxg_obs.shape[0]:
			sys.exit('WARNING: cxg_obs does not contain the same number of rows as final matrix: {} vs {}'.format(mfinal_adata.X.shape[0], cxg_obs.shape[0]))
	else:
		# Go through donor and biosample diseases and calculate cxg field accordingly
		report_diseases(df, mfinal_obj.get('experimental_variable_disease', unreported_value))
		get_sex_ontology(df)
		cxg_obs = pd.merge(cxg_obs, df[['disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id']], left_on="raw_matrix_accession", right_index=True, how="left" )

	# Clean up columns in obs to follow cxg schema and drop any unnecessary fields
	clean_obs(celltype_col)

	# If final matrix file is h5ad, take expression matrix from .X to create cxg anndata
	results_file  = get_results_filename(mfinal_obj)
	cxg_var = pd.DataFrame(index = mfinal_adata.var.index.to_list())
	keep_types = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
	var_meta = mfinal_adata.var.select_dtypes(include=keep_types)
	cxg_adata = ad.AnnData(mfinal_adata.X, obs=cxg_obs, obsm=cxg_obsm, var=cxg_var, uns=cxg_uns)
	cxg_adata.var = cxg_adata.var.merge(var_meta, left_index=True, right_index=True, how='left')
	if not sparse.issparse(cxg_adata.X):
		cxg_adata.X = sparse.csr_matrix(cxg_adata.X)
	elif cxg_adata.X.getformat()=='csc':
		cxg_adata.X = sparse.csr_matrix(cxg_adata.X)


	# Convert gene symbols to ensembl and filter to approved set
	if len(feature_lengths) > 1 and len(mfinal_obj['genome_annotations'])==1:
		clean_var()

	compiled_annot = compile_annotations(ref_files)
	# For ATAC gene activity matrices, it is assumed there are no genes that are filtered
	if summary_assay != 'ATAC':
		set_ensembl(redundant, mfinal_obj['feature_keys'])
		cxg_adata_raw = filter_ensembl(cxg_adata_raw, compiled_annot)
		cxg_adata = filter_ensembl(cxg_adata, compiled_annot)
		add_zero()
	else:
		cxg_adata_raw = filter_ensembl(cxg_adata_raw, compiled_annot)
		cxg_adata = filter_ensembl(cxg_adata, compiled_annot)
		cxg_adata.var['feature_is_filtered'] = False
	
	if not sparse.issparse(cxg_adata_raw.X):
		cxg_adata_raw = ad.AnnData(X = sparse.csr_matrix(cxg_adata_raw.X), obs = cxg_adata_raw.obs, var = cxg_adata_raw.var)
	elif cxg_adata.X.getformat()=='csc':
		cxg_adata.X = sparse.csr_matrix(cxg_adata.X)

	# Copy over any additional data from mfinal_adata to cxg_adatda
	reserved_uns = ['schema_version', 'title', 'batch_condition', 'default_embedding', 'X_approximate_distribution']
	for i in mfinal_adata.uns.keys():
		if i not in reserved_uns:
			cxg_adata.uns[i] = mfinal_adata.uns[i]
	if mfinal_adata.obsp:
		cxg_adata.obsp = mfinal_adata.obsp

	cxg_adata.raw = cxg_adata_raw
	quality_check(cxg_adata)
	cxg_adata.write(results_file, compression = 'gzip')

	shutil.rmtree(tmp_dir)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
