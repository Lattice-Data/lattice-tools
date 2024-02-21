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
from datetime import datetime
import matplotlib.colors as mcolors

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
		'ethnicity',
		'causes_of_death.term_name',
		'diseases.term_id',
		'diseases.term_name',
		'family_medical_history',
		'living_at_sample_collection',
		'menopausal_status',
		'organism.taxon_id',
		'risk_score_tyrer_cuzick_lifetime',
		'smoker',
		'times_pregnant'
		],
	'sample': [
		'age_development_stage_redundancy',
		'uuid',
		'preservation_method',
		'biosample_ontology.term_id',
		'biosample_ontology.organ_slims',
		'summary_development_ontology_at_collection.development_slims',
		'summary_development_ontology_at_collection.term_id',
		'derivation_process',
		'diseases.term_id',
		'diseases.term_name',
		'disease_state',
		'menstrual_phase_at_collection',
		'source',
		'summary_body_mass_index_at_collection',
		'treatment_summary',
		'growth_medium',
		'genetic_modifications',
		'@type'
		],
	'tissue_section': [
		'uuid',
		'thickness',
		'thickness_units'
	],
	'suspension': [
		'cell_depletion_factors',
		'depleted_cell_types.term_name',
		'derivation_process',
		'dissociation_reagent',
		'dissociation_time',
		'dissociation_time_units',
		'enriched_cell_types.term_name',
		'enrichment_factors',
		'percent_cell_viability',
		'uuid',
		'suspension_type',
		'tissue_handling_interval',
		'@id'
		],
	'library': [
		'uuid',
		'protocol.assay_ontology.term_id',
		'starting_quantity',
		'starting_quantity_units',
		'@id'
	],
	'raw_matrix': [
		'assembly',
		'genome_annotation',
		'software'
	],
	'seq_run': [
		'platform'
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

antibody_metadata = {
	'antibody': [
		'oligo_sequence',
		'host_organism',
		'source',
		'product_ids',
		'clone_id',
		'control',
		'isotype'
	],
	'target': [
		'label',
		'organism.scientific_name'
	]
}


# Mapping of field name (object_type + "_" + property) and what needs to be in the final cxg h5ad
prop_map = {
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'sample_summary_development_ontology_at_collection_term_id': 'development_stage_ontology_term_id',
	'sample_age_development_stage_redundancy': 'donor_age_redundancy',
	'sample_disease_state': 'disease_state',
	'sample_summary_body_mass_index_at_collection': 'donor_BMI_at_collection',
	'sample_growth_medium': 'growth_medium',
	'sample_genetic_modifications': 'genetic_modifications',
	'sample_menstrual_phase_at_collection': 'menstrual_phase_at_collection',
	'library_protocol_assay_ontology_term_id': 'assay_ontology_term_id',
	'donor_sex': 'sex',
	'sample_@type': 'tissue_type',
	'donor_donor_id': 'donor_id',
	'donor_organism_taxon_id': 'organism_ontology_term_id',
	'donor_ethnicity': 'self_reported_ethnicity_ontology_term_id',
	'donor_age_display': 'donor_age',
	'donor_risk_score_tyrer_cuzick_lifetime': 'tyrer_cuzick_lifetime_risk',
	'donor_smoker': 'donor_smoking_status',
	'donor_causes_of_death_term_name': 'donor_cause_of_death',
	'matrix_description': 'title',
	'matrix_default_embedding': 'default_embedding',
	'matrix_is_primary_data': 'is_primary_data',
	'cell_annotation_author_cell_type': 'author_cell_type',
	'cell_annotation_cell_ontology_term_id': 'cell_type_ontology_term_id',
	'cell_annotation_cell_state': 'cell_state',
	'suspension_suspension_type': 'suspension_type',
	'suspension_enriched_cell_types_term_name': 'suspension_enriched_cell_types',
	'suspension_depleted_cell_types_term_name': 'suspension_depleted_cell_types',
	'suspension_cell_depletion_factors': 'suspension_depletion_factors',
	'suspension_tissue_handling_interval': 'tissue_handling_interval',
	'antibody_oligo_sequence': 'barcode',
	'antibody_source': 'vendor',
	'antibody_product_ids': 'vender_product_ids',
	'antibody_clone_id': 'clone_id',
	'antibody_isotype': 'isotype',
	'antibody_host_organism': 'host_organism',
	'target_organism_scientific_name': 'target_organism',
	'raw_matrix_software': 'alignment_software',
	'raw_matrix_genome_annotation': 'mapped_reference_annotation',
	'raw_matrix_assembly': 'mapped_reference_assembly',
	'seq_run_platform': 'sequencing_platform'
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
warning_list = []

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
	global mfinal_obj
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
	sequencing_runs = []
	raw_seq_ids = []
	raw_seqs = []
	seq_run_ids = []
	seq_runs = []

	if start_type == None:
		for i in lib_ids:
			obj = lattice.get_object(i, connection)
			if mfinal_obj.get('output_types') == ['gene quantifications']:
				if obj.get('assay') in ['scRNA-seq','snRNA-seq','spatial transcriptomics','bulk RNA-seq', 'snATAC-seq']:
					libraries.append(obj)
			elif mfinal_obj.get('output_types') == ['antibody capture quantifications']:
				if obj.get('assay') == 'CITE-seq':
					libraries.append(obj)
			elif len(mfinal_obj.get('output_types')>1):
				logging.error("ERROR: The flattener cannot flatten multimodal ProcessedMatrixFile")
				sys.exit("ERROR: The flattener cannot flatten multimodal ProcessedMatrixFile")
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
		for d in input_object['derived_from']:
			obj = lattice.get_object(d, connection)
			if 'RawSequenceFile' in obj.get('@type'):
				for run in obj['derived_from']:
					if run.get('@id') not in seq_run_ids:
						seq_runs.append(run)
						seq_run_ids.append(run.get('@id'))
			else:
				break

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
		'tissue_section': tissue_sections,
		'seq_run':  seq_runs
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
		if path[0] == '@type':
			value = obj.get('@type')[0]
			return value
		else:
			return obj.get(prop, unreported_value)
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2, unreported_value) for i in obj[key1]]
			return list(set(values))
		elif obj.get(key1):
			value = obj[key1].get(key2, unreported_value)
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
		if prop == 'family_medical_history':
			history_list = get_value(obj, prop)
			if history_list != 'unknown':
				for history in history_list:
					ontology = lattice.get_object(history.get('diagnosis'), connection)
					key = 'family_history_' + str(ontology.get('term_name')).replace(' ','_')
					values_to_add[key] = history.get('present')
		elif prop == 'ethnicity':
			ethnicity_list = []
			ethnicity_dict_list = get_value(obj,prop)
			if ethnicity_dict_list != None:
				for ethnicity_dict in ethnicity_dict_list:
					if ethnicity_dict.get('term_id') == 'NCIT:C17998':
						ethnicity_list.append('unknown')
					else:
						ethnicity_list.append(ethnicity_dict.get('term_id'))
				ethnicity_list.sort()
				value = ','.join(ethnicity_list)
				latkey = (obj_type + '_' + prop).replace('.','_')
				key = prop_map.get(latkey, latkey)
				values_to_add[key] = value
		else:
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
		if prop == 'family_medical_history':
			values_df = pd.DataFrame()
			unknowns = []
			for obj in objs:
				ident = obj.get('@id')
				history_list = get_value(obj, prop)
				if history_list != 'unknown':
					for history in history_list:
						ontology = lattice.get_object(history.get('diagnosis'), connection)
						key = 'family_history_' + str(ontology.get('term_name')).replace(' ','_')
						value = str(history.get('present'))
						values_df.loc[key,ident] = value
				else:
					values_df[ident] = np.nan
					unknowns.append(ident)
			for i in unknowns:
				values_df[i] = 'unknown'
			for index, row in values_df.iterrows():
				values_to_add[index] = 'pooled [{}]'.format(','.join(row.to_list()))
		elif prop == 'ethnicity':
			values_df = pd.DataFrame()
			latkey = (obj_type + '_' + prop).replace('.','_')
			key = prop_map.get(latkey, latkey)
			for obj in objs:
				ident = obj.get('@id')
				ethnicity_list = []
				ethnicity_dict_list = get_value(obj,prop)
				if ethnicity_dict_list != None:
					for ethnicity_dict in ethnicity_dict_list:
						if ethnicity_dict.get('term_id') == 'NCIT:C17998':
							ethnicity_list.append('unknown')
						else:
							ethnicity_list.append(ethnicity_dict.get('term_id'))
				if len(ethnicity_list) == 1 and ethnicity_list[0] == 'unknown':
					values_df.loc[key,ident] = 'unknown'
				elif 'unknown' in ethnicity_list:
					values_df.loc[key,ident] = 'unknown'
				elif len(set(ethnicity_list)) == len(ethnicity_list):
					value = ethnicity_list[0]
					values_df.loc[key,ident] = value
				elif 'unknown' in ethnicity_list:
					values_df.loc[key,ident] = 'unknown'
			for index, row in values_df.iterrows():
				values_to_add[index] = str(row[0])
		else:
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
			cxg_fields = ['disease_ontology_term_id', 'organism_ontology_term_id',\
							 'sex', 'tissue_ontology_term_id', 'development_stage_ontology_term_id']
			if len(value_set) > 1:
				if key in cxg_fields:
					if key == 'development_stage_ontology_term_id':
						dev_in_all = list(set.intersection(*map(set, dev_list)))
						if dev_in_all == []:
							logging.error('ERROR: There is no common development_slims that can be used for development_stage_ontology_term_id')
							sys.exit("ERROR: There is no common development_slims that can be used for development_stage_ontology_term_id")
						else:
							query_url = urljoin(server, 'search/?type=OntologyTerm&term_name=' + dev_in_all[0] + '&format=json')
							r = requests.get(query_url, auth=connection.auth)
							try:
								r.raise_for_status()
							except requests.HTTPError:
								logging.error('ERROR: Unable to get development_slims as development_stage ontology: {}'.format(query_url))
								sys.exit("ERROR: Unable to get development_slims as development_stage ontology: {}".format(query_url))
							else:
								if r.json()['total']==1:
									values_to_add[key] = r.json()['@graph'][0]['term_id']
								else:
									logging.error('ERROR: Unable to get development_slims as development_stage ontology: {}'.format(query_url))
									sys.exit("ERROR: Unable to get development_slims as development_stage ontology: {}".format(query_url))
					elif key == 'sex':
						values_to_add[key] = 'unknown'
					else:
						logging.error('ERROR: Cxg field is a list')
						sys.exit("ERROR: Cxg field is a list")
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
			logging.error('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
			sys.exit('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
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
			logging.error('ERROR: The following error occured while downloading file {}: \n {}'.format(file_name,e))
			sys.exit('ERROR: The following error occured while downloading file{}: \n {}'.format(file_name,e))
		else:
			ftp.quit()
			print(file_name + ' downloaded')
	else:
		logging.error('ERROR: File {} has no uri defined'.format(file_obj['@id']))
		sys.exit('ERROR: File {} has no uri defined'.format(file_obj['@id']))


# Download entire directory contents from S3
def download_directory(download_url, directory):
	bucket_name = download_url.split('/')[2]
	spatial_folder = download_url.replace('s3://{}/'.format(bucket_name),"")
	s3client = boto3.client("s3")
	results = s3client.list_objects_v2(Bucket=bucket_name, Prefix=spatial_folder, Delimiter='/')
	os.mkdir(tmp_dir+"/spatial")
	for file in results.get('Contents'):
		if file.get('Size') == 0:
			continue
		file_path = file.get('Key')
		file_name = file_path.split('/')[-1]
		try:
			s3client.download_file(bucket_name, file_path, directory + '/spatial/' + file_name)
		except subprocess.CalledProcessError as e:
			logging.error('ERROR: Failed to find file S3://{}/{}'.format(bucket,file_path))
			sys.exit('ERROR: Failed to find file s3://{}/{}'.format(bucket,file_path))
		else:
			print(file_name + ' downloaded')


# Compile all reference annotations for var features into one pandas df
def compile_annotations(files):
	ids = pd.DataFrame()
	client = boto3.client('s3')
	bucket_name = 'submissions-lattice'
	for key in files:
		filename = tmp_dir + "/" + files[key]
		if os.path.exists(filename) == False:
			try:
				client.download_file(bucket_name, 'cxg_migration/var_refs/' + files[key], filename)
			except subprocess.CalledProcessError as e:
				logging.error('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
				sys.exit('ERROR: Failed to find file {} on s3'.format(file_obj.get('@id')))
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
	global mfinal_obj
	if adata.obs.isnull().values.any():
		warning_list.append("WARNING: There is at least one 'NaN' value in the following cxg anndata obs columns: {}".format(adata.obs.columns[adata.obs.isna().any()].tolist()))
	elif 'default_visualization' in adata.uns:
		if adata.uns['default_visualization'] not in adata.obs.values:
			logging.error('ERROR: The default_visualization field is not in the cxg anndata obs dataframe.')
			sys.exit("ERROR: The default_visualization field is not in the cxg anndata obs dataframe.")
	elif mfinal_obj['X_normalized'] == True and mfinal_obj['assays'] != ['snATAC-seq']:
		if len(adata.var.index.tolist()) > len(adata.raw.var.index.tolist()):
			logging.error('ERROR: There are more genes in normalized genes than in raw matrix.')
			sys.exit("ERROR: There are more genes in normalized genes than in raw matrix.")


# Check validity of colors before adding to cxg_adata.uns
def colors_check(adata, color_column, column_name):
	column_name = column_name.replace('_colors','')
	# Check that obs column exists
	if column_name not in adata.obs.columns:
		error = 'the corresponding column is not present in obs.'
		return False, error
	# Check that the corresponding column is the right datatype
	if column_name in adata.obs.columns:
		if adata.obs[column_name].dtype.name != 'category':
			error = 'the corresponding column in obs. is the wrong datatype ({})'.format(adata.obs[column_name].dtype.name)
			return False, error
	# Verify color_column is a numpy array
	if color_column is None or not isinstance(color_column, np.ndarray):
		error = 'the column is not a numpy array.'
		return False, error
	# Verify that the numpy array contains strings
	if not all(isinstance(color,str) for color in color_column):
		error = 'the column does not contain strings.'
		return False, error
	# Verify that we have atleast as many colors as unique values in the obs column
	if len(color_column) < len(adata.obs[column_name].unique()):
		error = 'the column has less colors than unique values in the corresponding obs. column.'
		return False, error
	# Verify that either all colors are hex OR all colors are CSS4 named colors strings
	all_hex_colors = all(re.match(r"^#([0-9a-fA-F]{6})$", color) for color in color_column)
	all_css4_colors = all(color in mcolors.CSS4_COLORS for color in color_column)
	if not (all_hex_colors or all_css4_colors):
		error = 'the colors are not all hex or CSS4 named color strings.'
		return False, error
	else:
		error = 'none'
		return True, error


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
def concat_list(anndata_list,column,uns_merge):
	if column != 'none':
		suffix = 0
		for a in anndata_list:
			a.var.rename(columns={column:column + '-' + str(suffix)}, inplace=True)
			suffix += 1
	if uns_merge == True:
		concat_result = ad.concat(anndata_list,index_unique=None, join='outer', merge = 'unique',  uns_merge='first')
	else:
		concat_result = ad.concat(anndata_list,index_unique=None, join='outer', merge = 'unique')
	redundants = [i for i,c in collections.Counter(concat_result.obs.index.to_list()).items() if c>1]
	if len(redundants)>0:
		logging.error('ERROR: cell IDs are found in multiple raw matrix files.\t{}'.format(redundants))
		sys.exit('ERROR: cell IDs are found in multiple raw matrix files.\t{}'.format(redundants))
	drop_columns = [c for c in concat_result.obs.columns if c!='raw_matrix_accession']
	if drop_columns:
		concat_result.obs.drop(columns=drop_columns, inplace=True)
	return concat_result

# Determine reported disease as unique of sample and donor diseases, removing unreported value
def report_diseases(mxr_df, exp_disease):
	mxr_df['reported_diseases'] = mxr_df['sample_diseases_term_name'] + ',' + mxr_df['donor_diseases_term_name']
	mxr_df['reported_diseases'] = mxr_df['reported_diseases'].apply(lambda x: '[{}]'.format(','.join([i for i in set(x.split(',')) if i!=unreported_value])))
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
		exp_disease_aslist = ['[{}]'.format(x['term_name']) for x in exp_disease]
		exp_disease_aslist.extend(['none','[]'])
		if len([x for x in total_reported if x not in exp_disease_aslist])==0:
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
		relevant_objects = gather_objects(susp_obj, start_type='suspension')
		for obj_type in obj_type_subset:
		 	objs = relevant_objects.get(obj_type, [])
		 	if len(objs) == 1:
		 		gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, objs)
		 	else:
		 		logging.error('ERROR: Could not find suspension for demultiplexed donor: {}'.format(obj_type))
		 		sys.exit('ERROR: Could not find suspension for demultiplexed donor: {}'.format(obj_type))
		row_to_add = pd.Series(values_to_add)
		susp_df = pd.concat([susp_df,row_to_add.to_frame().T], ignore_index=True)
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
		if cxg_adata.layers:
			for layer in cxg_adata.layers:
				new_layer = sparse.csr_matrix((cxg_adata.layers[layer].data, cxg_adata.layers[layer].indices, cxg_adata.layers[layer].indptr), shape = cxg_adata_raw.shape)
				new_adata.layers[layer] = new_layer
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
			if len(drop_unmapped) > 0:
				warning_list.append('WARNING: {} unmapped gene_symbols were dropped. Full list available in logging file. Preview: {}'.format(len(drop_unmapped),drop_unmapped[:10]))
				warning_list.append('WARNING: Full list of the {} unmapped gene_ids dropped:\t{}'.format(len(drop_unmapped),drop_unmapped))
			cxg_adata.var = pd.merge(cxg_adata.var, cxg_adata_raw.var, left_index=True, right_index=True, how='left', copy = True)
			cxg_adata.var = cxg_adata.var.set_index('gene_ids', drop=True)
			cxg_adata_raw.var  = cxg_adata_raw.var.set_index('gene_ids', drop=True)
			cxg_adata.var.index.name = None
			cxg_adata_raw.var.index.name = None

			# Drop redundant by Ensembl ID
			drop_redundant = list(set(redundant).intersection(set(cxg_adata.var.index.to_list())))
			if len(drop_redundant) > 0:
				warning_list.append('WARNING: {} redundant gene_ids dropped:\t{}'.format(len(drop_redundant),drop_redundant))
			cxg_adata = cxg_adata[:, [i for i in cxg_adata.var.index.to_list() if i not in redundant]]

		else:
			warning_list.append("WARNING: raw matrix does not have genes_ids column")
	elif feature_keys == ['Ensembl gene ID']:
		cxg_adata_raw.var_names_make_unique()
		cxg_adata_raw.var  = cxg_adata_raw.var.set_index('gene_ids', drop=True)
		cxg_adata_raw.var.index.name = None
		unique_to_norm =  set(cxg_adata.var.index.to_list()).difference(set(cxg_adata_raw.var.index.to_list()))
		if len(unique_to_norm) > 0:
			warning_list.append("WARNING: normalized matrix contains {} Ensembl IDs not in raw".format(unique_to_norm))


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
	cxg_adata_raw_ensembl = concat_list(cxg_adata_lst,'gene_symbols',False)

	# Join raw matrices on gene symbol, ensembl stored as metadata. Add suffix to make unique, using '.' as to match R default
	count = 0
	for cxg_adata in cxg_adata_lst:
		cxg_adata.var['gene_ids'] = cxg_adata.var.index
		cxg_adata.var =  cxg_adata.var.set_index('gene_symbols' + '-' + str(count), drop=True)
		count += 1
		cxg_adata.var_names_make_unique(join = '.')
	cxg_adata_raw_symbol = concat_list(cxg_adata_lst,'gene_ids',False)

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
	cxg_adata_raw_ensembl.var_names_make_unique(join = '.')

	all_remove = list(set(stats['multiple_ensembl'] + stats['multiple_symbols']))

	for key in stats:
		stats[key] = set(stats[key])
		overlap_norm = set(mfinal_adata_genes).intersection(stats[key])
		warning_list.append("WARNING: Full list of {}\t{}\t{}\t{}".format(key, len(stats[key]), len(overlap_norm), overlap_norm))
		warning_list.append("WARNING: Genes with {}. {} in raw. {} in normalized. Full list available in logging file. Preview of normalized: {}".format(key, len(stats[key]), len(overlap_norm), list(overlap_norm)[:10]))

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


# Add antibody metadata to var and raw.var
def map_antibody():
	global cxg_adata
	global cxg_adata_raw
	global mfinal_obj
	antibody_meta = pd.DataFrame()
	for anti_mapping in mfinal_obj.get('antibody_mappings'):
		values_to_add = {}
		antibody = anti_mapping.get('antibody')
		gather_metdata('antibody', antibody_metadata['antibody'], values_to_add, [antibody])
		values_to_add['host_organism'] = re.sub(r'/organisms/(.*)/', r'\1', values_to_add['host_organism'])
		if not antibody.get('control'):
			target = None
			if len(antibody.get('targets')) > 1:
				for t in antibody.get('targets'):
					name = t.get('organism').get('scientific_name')
					if name == cxg_adata.obs['organism'].unique()[0]:
						target = [t]
			else:
				target = antibody.get('targets')
			gather_metdata('target', antibody_metadata['target'], values_to_add, target)
			values_to_add['feature_name'] = values_to_add['target_label']
		else:
			values_to_add['feature_name'] = '{} {} (control)'.format(values_to_add['host_organism'], values_to_add['isotype'])
			for val in antibody_metadata['target']:
				latkey = ('target_' + val).replace('.', '_')
				key = prop_map.get(latkey, latkey)
				values_to_add[key] = 'na'
		row_to_add = pd.DataFrame(values_to_add, index=[anti_mapping.get('label')])
		antibody_meta = pd.concat([antibody_meta, row_to_add])
	cxg_adata.var = pd.merge(cxg_adata.var, antibody_meta, left_index=True, right_index=True, how='left')
	cxg_adata_raw.var = pd.merge(cxg_adata_raw.var, antibody_meta, left_index=True, right_index=True, how='left')
	
	cxg_adata.var['author_index'] = cxg_adata.var.index
	cxg_adata_raw.var['author_index'] = cxg_adata.var.index
	cxg_adata_raw.var.drop(columns=['genome'], inplace=True)
	cxg_adata.var['feature_biotype'] = 'antibody-derived tags'
	cxg_adata.var.set_index('feature_name', inplace=True, drop=False)
	cxg_adata_raw.var.set_index('feature_name', inplace=True, drop=False)
	cxg_adata.var_names_make_unique(join='-')
	cxg_adata_raw.var_names_make_unique(join='-')


# Final touches for obs columns, modifying any Lattice fields to fit cxg schema
def clean_obs():
	global cxg_obs
	global mfinal_obj
	# For columns in mfinal_obj that contain continuous cell metrics, they are transferred to cxg_obs as float datatype
	# WILL NEED TO REVISIT IF FINAL MATRIX CONTAINS MULTIPLE LAYERS THAT WE ARE WRANGLING
	for author_col in mfinal_obj.get('author_columns',[]):
		if author_col in mfinal_adata.obs.columns.to_list():
			cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[author_col]], left_index=True, right_index=True, how='left')
		else:
			warning_list.append("WARNING: author_column not in final matrix: {}".format(author_col))
	# if obs category suspension_type does not exist in dataset, create column and fill values with na (for spatial assay)
	if 'suspension_type' not in cxg_obs.columns:
		cxg_obs.insert(len(cxg_obs.columns),'suspension_type', 'na')
	elif cxg_obs['suspension_type'].isnull().values.any():
		cxg_obs['suspension_type'].fillna(value='na', inplace=True)
	
	add_units = {'tissue_section_thickness': 'tissue_section_thickness_units',
				'suspension_dissociation_time': 'suspension_dissociation_time_units',
				'library_starting_quantity': 'library_starting_quantity_units'}
	for field in add_units.keys():
		if field in cxg_obs.columns:
			cxg_obs[field] = cxg_obs[field].astype(str) + " " + cxg_obs[add_units[field]].astype(str)
			cxg_obs.drop(columns=add_units[field], inplace=True)
			cxg_obs[field].replace({'unknown unknown': 'unknown'}, inplace=True)

	make_numeric = ['suspension_percent_cell_viability','donor_BMI_at_collection']
	for field in make_numeric:
		if field in cxg_obs.columns:
			if True in cxg_obs[field].str.contains('[<>-]|'+unreported_value+'|'+'pooled', regex=True).to_list():
				if True not in cxg_obs[field].str.contains('[<>-]', regex=True).to_list():
					cxg_obs[field].replace({'unknown':np.nan}, inplace=True) 
					cxg_obs[field][np.where(cxg_obs[field].str.contains('pooled') == True)[0].tolist()] = np.nan
					cxg_obs[field]  = cxg_obs[field].astype('float')
			else: 
				cxg_obs[field]  = cxg_obs[field].astype('float')

	change_unreported = ['suspension_enriched_cell_types','suspension_depleted_cell_types','suspension_enrichment_factors','suspension_depletion_factors','disease_state','cell_state']
	for field in change_unreported:
		if field in cxg_obs.columns.to_list():
			cxg_obs[field].replace({unreported_value: 'na'}, inplace=True)
	valid_tissue_types = ['tissue', 'organoid', 'cellculture']
	cxg_obs['tissue_type'] = cxg_obs['tissue_type'].str.lower()
	for i in cxg_obs['tissue_type'].unique().tolist():
		if i == 'cellculture':
			cxg_obs['tissue_type'].replace({'cellculture':'cell culture'}, inplace=True)
		if i not in valid_tissue_types:
			logging.error('ERROR: not a valid tissue type:\t{}'.format(i))
			print('ERROR: not a valid tissue type:\t{}'.format(i))
	if mfinal_obj['is_primary_data'] == 'mixed':
		primary_portion = mfinal_obj.get('primary_portion')
		cxg_obs['is_primary_data'] = False
		cxg_obs.loc[cxg_obs[primary_portion.get('obs_field')].isin(primary_portion.get('values')),'is_primary_data'] = True

	cxg_obs[[i for i in cxg_obs.columns.tolist() if i.startswith('family_history_')]] = \
		cxg_obs[[i for i in cxg_obs.columns.tolist() if i.startswith('family_history_')]].fillna(value='unknown')


# Drop any intermediate or optional fields that are all empty
def drop_cols(celltype_col):
	global cxg_obs
	optional_columns = ['donor_BMI_at_collection', 'donor_family_medical_history', 'reported_diseases', 'donor_times_pregnant', 'sample_preservation_method',\
			'sample_treatment_summary', 'suspension_uuid', 'tissue_section_thickness', 'tissue_section_thickness_units','cell_state', 'disease_state',\
			'suspension_enriched_cell_types', 'suspension_enrichment_factors', 'suspension_depletion_factors', 'tyrer_cuzick_lifetime_risk',\
			'donor_living_at_sample_collection','donor_menopausal_status','donor_smoking_status','sample_derivation_process','suspension_dissociation_reagent',\
			'suspension_dissociation_time','suspension_depleted_cell_types','suspension_derivation_process','suspension_percent_cell_viability',\
			'library_starting_quantity','library_starting_quantity_units','tissue_handling_interval','suspension_dissociation_time_units','alignment_software',\
			'mapped_reference_annotation','mapped_reference_assembly','sequencing_platform','sample_source','donor_cause_of_death', 'growth_medium','genetic_modifications',
			'menstrual_phase_at_collection']
	
	if 'sequencing_platform' in cxg_obs.columns:
		if cxg_obs['sequencing_platform'].isnull().values.any():
			cxg_obs['sequencing_platform'].fillna(unreported_value, inplace=True)
	for col in optional_columns:
		if col in cxg_obs.columns.to_list():
			col_content = cxg_obs[col].unique()
			if len(col_content) == 1:
				if col_content[0] == unreported_value or col_content[0] == '[' + unreported_value + ']' or col_content[0] == '[]':
					cxg_obs.drop(columns=col, inplace=True)

	if len(cxg_obs['donor_age_redundancy'].unique()) == 1:
		if cxg_obs['donor_age_redundancy'].unique():
			cxg_obs.drop(columns='donor_age', inplace=True)
	columns_to_drop = ['raw_matrix_accession', celltype_col, 'sample_diseases_term_id', 'sample_diseases_term_name', 'sample_biosample_ontology_organ_slims',\
			'donor_diseases_term_id', 'donor_diseases_term_name', 'batch', 'library_@id_x', 'library_@id_y', 'author_donor_x', 'author_donor_y',\
			'library_authordonor', 'author_donor_@id', 'library_donor_@id', 'suspension_@id', 'library_@id', 'sex', 'sample_biosample_ontology_cell_slims',\
			'sample_summary_development_ontology_at_collection_development_slims','donor_age_redundancy','in_tissue','array_row','array_col']
	for column_drop in  columns_to_drop: 
		if column_drop in cxg_obs.columns.to_list():
			cxg_obs.drop(columns=column_drop, inplace=True)


# Add ontology term names to CXG standardized columns
def add_labels():
	global cxg_adata
	global cxg_adata_raw
	query_url = urljoin(server, 'search/?type=OntologyTerm&field=term_id&field=term_name&limit=all')
	r = requests.get(query_url, auth=connection.auth)
	try:
		r.raise_for_status()
	except requests.HTTPError:
		logging.error('ERROR: Error in ontology term ids and names: {}'.format(query_url))
		sys.exit("ERROR: Error in ontology term ids and names: {}".format(query_url))
	else:
		ontology_df = pd.DataFrame(r.json()['@graph'])
	id_cols = ['assay_ontology_term_id','disease_ontology_term_id','cell_type_ontology_term_id','development_stage_ontology_term_id','sex_ontology_term_id',\
			'tissue_ontology_term_id','organism_ontology_term_id','self_reported_ethnicity_ontology_term_id']
	for col in id_cols:
		name_col = col.replace("_ontology_term_id","")
		cxg_adata.obs[name_col] = cxg_adata.obs[col]
		for term_id in cxg_adata.obs[name_col].unique():
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
			elif len(ontology_df.loc[ontology_df['term_id']==term_id,'term_name'].unique() == 1):
				term_name = ontology_df.loc[ontology_df['term_id']==term_id,'term_name'].unique()[0]
			else:
				logging.error('ERROR: Found more than single ontology term name for id: {}\t{}'.format(term_id, ontology_df.loc[ontology_df['term_id']==term_id,'term_name'].unique()))
				sys.exit("ERROR: Found more than single ontology term name for id: {}\t{}".format(term_id, ontology_df.loc[ontology_df['term_id']==term_id,'term_name'].unique()))
			cxg_adata.obs[name_col].replace(term_id, term_name, inplace=True)


# Look at matrix and only convert to sparse if density is less than 0.5
def check_matrix(m):
	if not sparse.issparse(m):
		if (np.count_nonzero(m)/np.prod(m.shape)) <= 0.5:
			m = sparse.csr_matrix(m)
	elif m.getformat()=='csc':
		m = sparse.csr_matrix(m)
	return m


def main(mfinal_id):
	global mfinal_obj
	global mfinal_adata
	global cxg_adata
	global cxg_adata_raw
	global cxg_obs
	mfinal_obj = lattice.get_object(mfinal_id, connection)
	logging.basicConfig(filename= mfinal_id + '_outfile_flattener.log', filemode='w', level=logging.INFO)
	# Adding date and time to top of logging file
	time_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
	logging.info("Date and time of flattener run: " + time_date)
	# Suppressing specific warnings from anndata
	logging.captureWarnings(True)

	# confirm that the identifier you've provided corresponds to a ProcessedMatrixFile
	mfinal_type = mfinal_obj['@type'][0]
	summary_assay = ''
	if mfinal_type != 'ProcessedMatrixFile':
		logging.error('ERROR: {} is not a ProcessedMatrixFile, but a {}'.format(mfinal_id, mfinal_type))
		sys.exit('ERROR: {} is not a ProcessedMatrixFile, but a {}'.format(mfinal_id, mfinal_type))

	if mfinal_obj['output_types'] == ['gene quantifications']:
		if mfinal_obj['assays'] == ['snATAC-seq']:
			summary_assay = 'ATAC'
		else:
			summary_assay = 'RNA'
	elif mfinal_obj['output_types'] == ['antibody capture quantifications']:
		summary_assay = 'CITE'
	else:
		logging.error('ERROR: Unexpected assay types to generate cxg h5ad: {} {}'.format(mfinal_obj['assays'], mfinal_obj['output_types']))
		sys.exit("ERROR: Unexpected assay types to generate cxg h5ad: {} {}".format(mfinal_obj['assays'], mfinal_obj['output_types']))


	# Dataframe that contains experimental metadata keyed off of raw matrix
	df = pd.DataFrame()

	results = {}
	
	# Checking for presence of matrix_files, and creating if not present
	if os.path.exists(tmp_dir) == False:
		os.mkdir(tmp_dir)
		
	# Checking for presence of h5ad, and downloading if not present
	if os.path.exists(tmp_dir + '/' + mfinal_obj['accession'] + '.h5ad'):
		print(mfinal_obj['accession'] + '.h5ad' + ' was found locally')
	else:
		download_file(mfinal_obj, tmp_dir)

	# Get list of unique final cell identifiers
	file_url = mfinal_obj['s3_uri']
	file_ext = file_url.split('.')[-1]
	mfinal_local_path = '{}/{}.{}'.format(tmp_dir, mfinal_obj['accession'], file_ext)
	mfinal_adata = sc.read_h5ad(mfinal_local_path)
	mfinal_cell_identifiers = mfinal_adata.obs.index.to_list()

	cxg_adata_lst = []
	redundant = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	donor_susp = {}
	library_susp = {}
	mapping_error = False
	error_info = {}

	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		relevant_objects = gather_objects(mxr)
		values_to_add = {}

		# Get raw matrix metadata
		gather_metdata('raw_matrix', cell_metadata['raw_matrix'], values_to_add, [mxr])

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
			for obj_type in cell_metadata.keys():
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, objs)
				elif len(objs) > 1:
					gather_pooled_metadata(obj_type, cell_metadata[obj_type], values_to_add, objs)
		row_to_add = pd.DataFrame(values_to_add, index=[mxr['@id']], dtype=str)
		
		# Add anndata to list of final raw anndatas, only for RNAseq
		if summary_assay in ['RNA','CITE']:
			# Checking for presence of mxr file and downloading if not present
			if mxr['s3_uri'].endswith('h5'):
				if os.path.exists(tmp_dir + '/' + mxr_acc + '.h5'):
					print(mxr_acc + '.h5' + ' was found locally')
				else:
					download_file(mxr, tmp_dir)
			elif mxr['s3_uri'].endswith('h5ad'):
				if os.path.exists(tmp_dir + '/' + mxr_acc + '.h5ad'):
					print(mxr_acc + '.h5ad' + ' was found locally')
				else:
					download_file(mxr, tmp_dir)
			if mfinal_obj.get('spatial_s3_uri', None) and mfinal_obj['assays'] == ['spatial transcriptomics']:
				if mxr['s3_uri'].endswith('h5'):
					mxr_name = '{}.h5'.format(mxr_acc)
				elif mxr['s3_uri'].endswith('h5ad'):
					mxr_name = '{}.h5ad'.format(mxr_acc)
				# Checking for presence of spatial directory and redownloading if present
				if os.path.exists(tmp_dir + '/spatial'):
					shutil.rmtree(tmp_dir + '/spatial')
				download_directory(mfinal_obj['spatial_s3_uri'], tmp_dir)
				# If tissue_positions is present rename to tissue_positions_list and remove header
				if os.path.exists(tmp_dir + '/spatial/tissue_positions.csv') == True:
					fixed_file = pd.read_csv(tmp_dir + '/spatial/tissue_positions.csv', skiprows = 1, header = None)
					fixed_file.to_csv(tmp_dir + '/spatial/tissue_positions_list.csv', header = False, index = False)
					os.remove(tmp_dir + '/spatial/tissue_positions.csv')
				if 'spatial' in mfinal_adata.uns.keys():
					del mfinal_adata.uns['spatial']
				adata_raw = sc.read_visium(tmp_dir, count_file=mxr_name)
			elif mxr['s3_uri'].endswith('h5'):
				mxr_name = '{}.h5'.format(mxr_acc)
				adata_raw = sc.read_10x_h5('{}/{}'.format(tmp_dir,mxr_name), gex_only = False)
			elif mxr['s3_uri'].endswith('h5ad'):
				mxr_name = '{}.h5ad'.format(mxr_acc)
				adata_raw = sc.read_h5ad('{}/{}'.format(tmp_dir,mxr_name))
			else:
				logging.error('ERROR: Raw matrix file of unknown file extension: {}'.format(mxr['s3_uri']))
				sys.exit('ERROR: Raw matrix file of unknown file extension: {}'.format(mxr['s3_uri']))

			if summary_assay == 'RNA':
				row_to_add['mapped_reference_annotation'] = mxr['genome_annotation']
				adata_raw = adata_raw[:,adata_raw.var['feature_types']=='Gene Expression']
			else:
				adata_raw = adata_raw[:,adata_raw.var['feature_types']=='Antibody Capture']
			# only make var unique if all raw matrices are same annotation version
			if len(mfinal_obj.get('genome_annotations', [])) == 1:
				for g in [i for i,c in collections.Counter(adata_raw.var.index.to_list()).items() if c > 1]:
					if True in adata_raw.var.loc[g, 'gene_ids'].str.endswith("PAR_Y").to_list():
						redundant.extend([i for i in adata_raw.var.loc[g,'gene_ids'] if i.endswith('PAR_Y')])
					else:
						redundant.extend(adata_raw.var.loc[g,'gene_ids'].to_list())
				adata_raw.var_names_make_unique(join = '.')
			# Recreate cell_ids and subset raw matrix and add mxr_acc into obs
			if mfinal_obj.get('cell_label_mappings', None):
				concatenated_ids = concatenate_cell_id(mxr['@id'], adata_raw.obs_names, mfinal_cell_identifiers)
				adata_raw.obs_names = concatenated_ids
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
			else:
				overlapped_ids = list(set(mfinal_cell_identifiers).intersection(adata_raw.obs_names.to_list()))
			# Error check to see that cells in raw matrix match the cell in mfinal_adata
			cell_mapping_dct = {}
			if mfinal_obj.get('cell_label_mappings'):
				for mapping_dict in mfinal_obj.get('cell_label_mappings'):
					cell_mapping_dct[mapping_dict['raw_matrix']] = mapping_dict['label']
				mapping_label = cell_mapping_dct[mxr.get('@id')]
			if mfinal_obj.get('cell_label_location') == 'prefix':
				prefixes = []
				for mapping_dict in mfinal_obj['cell_label_mappings']: # Creating list of all prefixes
					prefixes.append(mapping_dict['label'])
				prefixes.remove(mapping_label) # Removing mapping_label prefix from list of all prefixes
				# Checking to make sure none of the other prefixes contain the mapping_label prefix, if they do, then make sure there's no false match
				if any(prefix.startswith(mapping_label) for prefix in prefixes):
					mfinal_with_label = [i for i in mfinal_cell_identifiers if i.startswith(mapping_label) and not i.startswith(tuple(prefixes))]
				else:
					mfinal_with_label = [i for i in mfinal_cell_identifiers if i.startswith(mapping_label)]
			elif mfinal_obj.get('cell_label_location') == 'suffix':
				mfinal_with_label = [i for i in mfinal_cell_identifiers if i.endswith(mapping_label)]
			else: 
				mfinal_with_label = [i for i in mfinal_cell_identifiers if i in adata_raw.obs_names]
			if len(overlapped_ids) == 0:
				if mfinal_obj['cell_label_location'] == 'prefix':
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
				adata_raw = adata_raw[overlapped_ids]
				adata_raw.obs['raw_matrix_accession'] = mxr['@id']
				cxg_adata_lst.append(adata_raw)
       
		df = pd.concat([df, row_to_add])
		redundant = list(set(redundant))
		
	# Removing mapped_reference_annotation if genome_annotations from ProcMatrixFile is empty
	if not mfinal_obj.get('genome_annotations', None):
		del df['mapped_reference_annotation']

	if mapping_error:
		logging.error('ERROR: There are {} mapping errors in cell_label_mappings:'.format(len(error_info.keys())))
		print("ERROR: There are {} mapping errors in cell_label_mappings:".format(len(error_info.keys())))
		for er in error_info.keys():
			print("RawMatrixFile: {}, {}".format(er, error_info[er]))
			logging.error("RawMatrixFile: {}, {}".format(er, error_info[er]))
		sys.exit()

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

	# For CITE and RNA datasets, concatenate all anndata objects in list, but no reconciling genes for CITE
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
		elif len(feature_lengths) > 1:
			cxg_adata_raw = concat_list(cxg_adata_lst,'gene_ids',True)
		else:
			cxg_adata_raw = concat_list(cxg_adata_lst,'none',True)
			if len(feature_lengths) == 1:
				if cxg_adata_raw.var.shape[0] != feature_lengths[0]:
					logging.error('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
					sys.exit('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
		cxg_adata_raw = cxg_adata_raw[mfinal_cell_identifiers]
		if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
			logging.error('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
			sys.exit('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
	elif summary_assay == 'CITE':
		cxg_adata_raw = concat_list(cxg_adata_lst,'none',False)
		if len(feature_lengths) == 1:
			if cxg_adata_raw.var.shape[0] != feature_lengths[0]:
				logging.error('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
				sys.exit('ERROR: There should be the same genes for raw matrices if only a single genome annotation')
		cxg_adata_raw = cxg_adata_raw[mfinal_cell_identifiers]
		if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
			logging.error('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
			sys.exit('ERROR: The number of cells do not match between final matrix and cxg h5ad.')
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
		if mfinal_adata.raw == None:
			cxg_adata_raw = ad.AnnData(sparse.csr_matrix(mfinal_adata.X.shape), var = mfinal_adata.var, obs = atac_obs)
		else:
			cxg_adata_raw = ad.AnnData(mfinal_adata.raw.X, var = mfinal_adata.var, obs = atac_obs)

	# Set uns and obsm parameters, moving over spatial information if applicable
	cxg_uns = ds_results
	if 'spatial' in cxg_adata_raw.uns.keys():
		cxg_uns['spatial'] = cxg_adata_raw.uns['spatial']
		spatial_lib = list(cxg_uns['spatial'].keys())[0]
		cxg_uns['image'] = cxg_uns['spatial'][spatial_lib]['images']['hires']
	cxg_obsm = mfinal_adata.obsm.copy()
	if mfinal_obj['assays'] == ['spatial transcriptomics']:
		if 'spatial' in cxg_adata_lst[0].obsm.keys():
			cxg_obsm['spatial'] = cxg_adata_lst[0][mfinal_adata.obs.index.to_list()].obsm['spatial']
			if 'X_spatial' not in mfinal_adata.obsm.keys():
				spatial_lib = list(cxg_uns['spatial'].keys())[0]
				cxg_obsm['X_spatial'] = cxg_obsm['spatial'] * cxg_uns['spatial'][spatial_lib]['scalefactors']['tissue_hires_scalef']
	if len([i for i in cxg_obsm.keys() if i.startswith('X_')]) < 1:
		logging.error("ERROR: At least one embedding that starts with 'X_' is required")
		sys.exit("ERROR: At least one embedding that starts with 'X_' is required")


	# Merge df with raw_obs according to raw_matrix_accession, and add additional cell metadata from mfinal_adata if available
	# Also add calculated fields to df 
	celltype_col = mfinal_obj['author_cell_type_column']
	if mfinal_adata.obs[celltype_col].dtype != 'object':
		mfinal_adata.obs[celltype_col] = mfinal_adata.obs[celltype_col].astype('string')
	cxg_obs = pd.merge(cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[celltype_col]], left_index=True, right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, annot_df, left_on=celltype_col, right_index=True, how='left')
	if cxg_obs['cell_type_ontology_term_id'].isnull().values.any():
		warning_list.append("WARNING: Cells did not sucessfully map to CellAnnotations with author cell type and counts: {}".\
			format(cxg_obs.loc[cxg_obs['cell_type_ontology_term_id'].isnull()==True, celltype_col].value_counts().to_dict()))
	if cxg_obs[celltype_col].isna().any():
		logging.error("ERROR: author_cell_type column contains 'NA' values, unable to perform CellAnnotation mapping.")
		sys.exit("ERROR: author_cell_type column contains 'NA' values, unable to perform CellAnnotation mapping.")
	if len([i for i in annot_df.index.to_list() if i not in cxg_obs[celltype_col].unique().tolist()]) > 0:
		warning_list.append("WARNING: CellAnnotation that is unmapped: {}\n".format([i for i in annot_df.index.to_list() if i not in cxg_obs[celltype_col].unique().tolist()]))

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
			logging.error('ERROR: cxg_obs does not contain the same number of rows as final matrix: {} vs {}'.format(mfinal_adata.X.shape[0], cxg_obs.shape[0]))
			sys.exit('ERROR: cxg_obs does not contain the same number of rows as final matrix: {} vs {}'.format(mfinal_adata.X.shape[0], cxg_obs.shape[0]))
	else:
		# Go through donor and biosample diseases and calculate cxg field accordingly
		report_diseases(df, mfinal_obj.get('experimental_variable_disease', unreported_value))
		get_sex_ontology(df)
		cxg_obs = pd.merge(cxg_obs, df[['disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id']], left_on="raw_matrix_accession", right_index=True, how="left" )

	# Clean up columns in obs to follow cxg schema and drop any unnecessary fields
	drop_cols(celltype_col)
	clean_obs()
	del cxg_adata_lst
	gc.collect()

	# Check that primary_portion.obs_field of ProcessedMatrixFile is present in cxg_obs
	if mfinal_obj.get('primary_portion', None): # Checking for presence of 'primary_portion'
		primary_portion = mfinal_obj.get('primary_portion')
		if primary_portion.get('obs_field') not in cxg_obs.columns:
			logging.error("ERROR: 'obs_field' value '{}' not found in cxg_obs columns".format(primary_portion.get('obs_field')))
			sys.exit("ERROR: 'obs_field' value '{}' not found in cxg_obs columns".format(primary_portion.get('obs_field')))

		# Check that all primary_portion.values of ProcessedMatrixFile are found in the 'obs_field' column of cxg_obs
		missing = [f for f in primary_portion.get('values') if f not in cxg_obs[primary_portion.get('obs_field')].tolist()]
		if missing:
			logging.error("ERROR: cxg_obs column '{}' doesn't contain values present in 'primary_portion.obs_field' of ProcessedMatrixFile: {}".format(primary_portion.get('obs_field'),missing))
			sys.exit("ERROR: cxg_obs column '{}' doesn't contain values present in 'primary_portion.obs_field' of ProcessedMatrixFile: {}".format(primary_portion.get('obs_field'),missing))

	# If final matrix file is h5ad, take expression matrix from .X to create cxg anndata
	results_file  = get_results_filename(mfinal_obj)
	mfinal_adata.var_names_make_unique()
	cxg_var = pd.DataFrame(index = mfinal_adata.var.index.to_list())
	keep_types = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
	if summary_assay == 'CITE':
		keep_types.append('object')
	var_meta = mfinal_adata.var.select_dtypes(include=keep_types)
	cxg_adata = ad.AnnData(mfinal_adata.X, obs=cxg_obs, obsm=cxg_obsm, var=cxg_var, uns=cxg_uns)
	cxg_adata.var = cxg_adata.var.merge(var_meta, left_index=True, right_index=True, how='left')
	
	# Removing feature_length column from var if present
	if 'feature_length' in cxg_adata.var.columns:
		adata.var.drop(columns=['feature_length'], inplace=True)

	# Check matrix density
	cxg_adata.X = check_matrix(cxg_adata.X)
	cxg_adata_raw.X = check_matrix(cxg_adata_raw.X)

	# Check that cxg_adata_raw.X is correct datatype
	if not cxg_adata_raw.X.dtype == 'float32':
		cxg_adata_raw.X = cxg_adata_raw.X.astype(np.float32) 
		
	# Adding layers from 'layers_to_keep' to cxg_adata.layers	
	if 'layers_to_keep' in mfinal_obj:
		for k in mfinal_obj['layers_to_keep']:
			cxg_adata.layers[k] = mfinal_adata.layers[k]
			cxg_adata.layers[k] = check_matrix(cxg_adata.layers[k])
				
	# Convert gene symbols to ensembl and filter to approved set
	if len(feature_lengths) > 1 and len(mfinal_obj['genome_annotations'])==1:
		clean_var()

	# For ATAC gene activity matrices, it is assumed there are no genes that are filtered
	# For CITE, standardize antibody index and metadata and no filtering

	if summary_assay == 'RNA':
		compiled_annot = compile_annotations(ref_files)
		set_ensembl(redundant, mfinal_obj['feature_keys'])
		cxg_adata_raw = filter_ensembl(cxg_adata_raw, compiled_annot)
		cxg_adata = filter_ensembl(cxg_adata, compiled_annot)
		add_zero()
	elif summary_assay == 'ATAC':
		compiled_annot = compile_annotations(ref_files)
		cxg_adata_raw = filter_ensembl(cxg_adata_raw, compiled_annot)
		cxg_adata = filter_ensembl(cxg_adata, compiled_annot)
		cxg_adata.var['feature_is_filtered'] = False
	elif summary_assay == 'CITE':
		add_labels()
		map_antibody()
		add_zero()

	# Copy over any additional data from mfinal_adata to cxg_adata
	reserved_uns = ['schema_version', 'title', 'default_embedding', 'X_approximate_distribution','schema_reference','citation']
	for i in mfinal_adata.uns.keys():
		if i == 'batch_condition':
			if not isinstance(mfinal_adata.uns['batch_condition'], list) and not isinstance(mfinal_adata.uns['batch_condition'], np.ndarray) :
				warning_list.append("WARNING: adata.uns['batch_condition'] is not a list and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
			else:
				if len([x for x in mfinal_adata.uns['batch_condition'] if x not in cxg_adata.obs.columns]) > 0:
					warning_list.append("WARNING: adata.uns['batch_condition'] contains column names not found and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
				elif len(set(mfinal_adata.uns['batch_condition'])) != len(mfinal_adata.uns['batch_condition']):
					warning_list.append("WARNING: adata.uns['batch_condition'] contains redundant column names and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
				else:
					cxg_adata.uns['batch_condition'] = mfinal_adata.uns['batch_condition']
		elif i.endswith('_colors'):
			colors_result = colors_check(cxg_adata, mfinal_adata.uns[i], i)
			if colors_result[0]:
				cxg_adata.uns[i] = mfinal_adata.uns[i]
			else:
				warning_list.append("WARNING: '{}' has been dropped from uns dict due to being invalid because '{}' \n".format(i,colors_result[1]))
		elif i not in reserved_uns:
			cxg_adata.uns[i] = mfinal_adata.uns[i]
		else:
			warning_list.append("WARNING: The key '{}' has been dropped from uns dict due to being reserved \n".format(i))

	if mfinal_adata.obsp:
		cxg_adata.obsp = mfinal_adata.obsp

	# Check if mfinal_obj matrix is normalized,if so set cxg_adata.raw to raw, if not then place raw in adata.X
	if mfinal_obj['X_normalized']:
		if summary_assay != 'ATAC' or mfinal_adata.raw != None:
			cxg_adata.raw = cxg_adata_raw
	else:
		cxg_adata.var['feature_is_filtered'] = False
		cxg_adata = ad.AnnData(cxg_adata_raw.X, obs=cxg_adata.obs, obsm=cxg_adata.obsm, var=cxg_adata.var, uns=cxg_adata.uns)
	quality_check(cxg_adata)
	cxg_adata.write_h5ad(results_file, compression = 'gzip')

	# Printing out list of warnings
	for n in warning_list:
		if 'WARNING: Full list' not in n:
			print(n, end = '\n')
		if 'Full list available in logging file.' not in n:
			logging.warning(n)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
