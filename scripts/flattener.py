import argparse
import boto3
import h5py
import lattice
import os
import pandas as pd
import shutil
import sys
import scanpy as sc


cell_metadata = {
	'donor': [
		'uuid',
		'age_display',
		'sex',
		'ethnicity.term_name',
		'ethnicity.term_id',
		'life_stage',
		'life_stage_term_id',
		'diseases.term_name',
		'diseases.term_id'
		],
	'sample': [
		'uuid',
		'preservation_method',
		'biosample_ontology.term_name',
		'biosample_ontology.term_id',
		'diseases.term_name',
		'diseases.term_id'
		],
	'suspension': [
		'uuid',
		'suspension_type'
		],
	'library': [
		'uuid',
		'protocol.title',
		'protocol.term_id',
		'cell_mapping_identifier'
		]
	}

dataset_metadata = {
	'dataset': [
		'award.title',
		'uuid'
		],
	'final_matrix': [
		'genome_annotation',
		'value_scale',
		'value_units',
		'normalized',
		'normalization_method',
		'mapping_column',
		'cell_mapping'
		]
	}


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
    return args


def gather_rawmatrices(derived_from):
	my_raw_matrices = []
	for identifier in derived_from:
		obj = lattice.get_object(identifier, connection)
		if obj['@type'][0] == 'MatrixFile' and obj['normalized'] != True and \
			obj['value_scale'] == 'linear' and 'cell calling' in obj['derivation_process']:
			my_raw_matrices.append(obj)
		else: 
			# grab the derived_from in case we need to go a layer deeper
			for i in obj['derived_from']:
				df_ids.append(i)
	if not my_raw_matrices:
		for identifier in df_ids:
			obj = lattice.get_object(identifier, connection)
			if obj['@type'][0] == 'MatrixFile' and obj['normalized'] != True and \
				obj['value_scale'] == 'linear' and 'cell calling' in obj['derivation_process']:
				my_raw_matrices.append(obj)
	return my_raw_matrices


def gather_objects(raw_matrix_file):
	libraries = []
	suspensions = []
	samples = []
	donors = []

	lib_ids = raw_matrix_file['libraries']
	sample_ids = []

	for i in lib_ids:
		obj = lattice.get_object(i, connection)
		libraries.append(obj)
		for o in obj['derived_from']:
			suspensions.append(o)
		for o in obj['donors']:
			donors.append(o)
	for o in suspensions:
		for i in o['derived_from']:
			sample_ids.append(i)
	samples = [lattice.get_object(i, connection) for i in sample_ids]

	return {
		'donor': donors,
		'sample': samples,
		'suspension': suspensions,
		'library': libraries
		}


def get_value(obj, prop):
	path = prop.split('.')
	if len(path) == 1:
		return obj.get(prop,' ')
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2,' ') for i in obj[key1]]
			return list(set(values))
		elif obj.get(key1):
			return obj[key1].get(key2)
		else:
			return obj.get(key1,' ')
	else:
		return 'unable to traverse more than 1 embedding'


def gather_metdata(obj_type, values_to_add, mxr_acc, objs):
	obj = objs[0]
	for prop in cell_metadata[obj_type]:
		value = get_value(obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		key = (obj_type + '_' + prop).replace('.', '_')
		values_to_add[key] = value


def gather_poooled_metadata(obj_type, values_to_add, mxr_acc, objs):
	for prop in cell_metadata[obj_type]:
		key = (obj_type + '_' + prop).replace('.', '_')
		# NEED TO DETERMINE HOW TO REPORT POOLED DATA
		values_to_add[key] = 'POOLED'


def report_dataset(matrix, dataset):
	ds_results = {}
	ds_obj = lattice.get_object(dataset, connection)
	for prop in dataset_metadata['dataset']:
		value = get_value(ds_obj, prop)
		ds_results['dataset' + '_' + prop.replace('.','_')] = value
	for prop in dataset_metadata['final_matrix']:
		value = get_value(matrix, prop)
		ds_results['matrix' + '_' + prop.replace('.','_')] = value
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


# Takes ds_results and converts information to cxg fields
# WILL NEED CORPORA VERSION INFORMATION AS INPUT, CURRENTLY HARDCODED
def organize_uns_data(ds_results):
	cxg_uns = {}
	cxg_uns['version'] = {}
	cxg_uns['version']['corpora_schema_version'] = '1.1.0'
	cxg_uns['version']['corpora_encoding_version'] = '0.2.0'
	cxg_uns['title'] = ds_results['dataset_award_title']
	cxg_uns['layer_descriptions'] = {}

	# Summarize normalization
	if ds_results['matrix_normalized']:
		normalization = "normalized using " + ds_results['matrix_normalization_method']
	else:
		normalization = "not normalized"
	cxg_uns['layer_descriptions']['X'] = ds_results['matrix_value_units'] + " counts; " + ds_results['matrix_value_scale'] + " scaling; " + normalization

	return cxg_uns


def main(mfinal_id):
	mfinal_obj = lattice.get_object(mfinal_id, connection)

	# confirm that the identifier you've provided corresponds to a MatrixFile
	mfinal_type = mfinal_obj['@type'][0]
	if mfinal_type != 'MatrixFile':
		sys.exit('{} is not a MatrixFile, but a {}'.format(mfinal_id, mfinal_type))

	# set the metadata keys based on defined metadata fields
	headers = []
	for obj_type in cell_metadata.keys():
		for prop in cell_metadata[obj_type]:
			key = (obj_type + '_' + prop).replace('.', '_')
			headers.append(key)

	tmp_dir = 'matrix_files'
	os.mkdir(tmp_dir)
	download_file(mfinal_obj, tmp_dir)

	# Create initial cell metadata will desired cellIDs and library cell_mapping_identifiers to allow for merging of metadata
	# Rename column name to match library_cell_mapping_identifiers
	final_local_path = '{}/{}.h5ad'.format(tmp_dir,mfinal_obj['accession'])
	final_matrix_adata = sc.read_h5ad(final_local_path)
	cell_metadata_df = final_matrix_adata.obs[[mfinal_obj['mapping_column']]]
	cell_metadata_df.columns = ['library_cell_mapping_identifier']

	# Create library metadata dataframe, which will eventually be merged back to cell_metadata_df
	library_metadata_df = pd.DataFrame(columns=headers)

	# get dataset-level metadata
	# reorganize cell_mapping_lst
	ds_results = report_dataset(mfinal_obj, mfinal_obj['dataset'])
	all_cell_mapping_lst = ds_results['matrix_cell_mapping']
	all_cell_mapping_dct = {}
	for mapping_dict in all_cell_mapping_lst:
		k = mapping_dict['cell_mapping_identifier']
		v = [mapping_dict['cell_mapping_string'], mapping_dict['cell_mapping_location'], mapping_dict['cell_mapping_connector']]
		all_cell_mapping_dct[k] = v

	# create list of anndata objects that need to be concatenated
	final_adata_lst = []
	final_batch_categories = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	for mxr in mxraws:
		mxr_acc = mxr['accession']
		# get all of the objects necessary to pull the desired metadata
		relevant_objects = gather_objects(mxr)
		values_to_add = {}
		for obj_type in cell_metadata.keys():
			if len(relevant_objects[obj_type]) == 1:
				gather_metdata(obj_type, values_to_add, mxr_acc, relevant_objects[obj_type])
			else:
				gather_pooled_metadata(obj_type, values_to_add, mxr_acc, relevant_objects[obj_type])
		download_file(mxr, tmp_dir)
		local_path = '{}/{}.h5'.format(tmp_dir,mxr_acc)

		# Add library information intoo library dataframe
		row_to_add = pd.Series(values_to_add)
		library_metadata_df = library_metadata_df.append(row_to_add, ignore_index=True)

		# Add new anndata to list of final anndatas
		# Should add error checking to make sure all matrices have the same number of vars
		adata_raw = sc.read_10x_h5(local_path)
		adata_raw.var_names_make_unique()
		final_adata_lst.append(adata_raw)
		final_batch_categories.append(all_cell_mapping_dct[values_to_add['library_cell_mapping_identifier']][0])

	# Merge library metadata and cell meta on cell_mapping_identifier
	merged_cell_df = cell_metadata_df.merge(library_metadata_df, on='library_cell_mapping_identifier', how='left')
	merged_cell_df.index = cell_metadata_df.index

	# Concatenate all anndata objects in list
	# Subset according to cellIDs in merged cell metadata
	if len(final_adata_lst) > 1:
		final_adata = final_adata_lst[0].concatenate(final_adata_lst[1:], batch_categories = final_batch_categories)
		final_adata = final_adata[merged_cell_df.index]
	else:
		final_adata = final_adata_lst[0]
		final_adata = final_adata[merged_cell_df.index]

	# NEED TO ADD A FIELD NAME CONVERSION FOR A HANDFUL OF FIELDS TO MEET CXG REQS
	# Load uns data info final anndata object
	final_adata.uns = organize_uns_data(ds_results)
	final_adata.obs = final_adata.obs[[]].merge(merged_cell_df, left_index=True, right_index=True, how='inner')

	# NEED TO CHANGE ENSEMBL IDS TO GENE_SYMBOLS
	# Move matrices into appropriate layer
	final_adata.raw = final_adata
	final_adata.X = final_matrix_adata.X

	# print the fields into a report
	results_file = "final_cxg.h5ad"
	final_adata.write(results_file)
	merged_cell_df.to_csv('cell_temp.csv', index=True)

	shutil.rmtree(tmp_dir)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
