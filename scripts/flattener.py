import argparse
import anndata as ad
import boto3
import lattice
import os
import pandas as pd
import shutil
import sys


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
		'protocol.term_id'
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
		'value_units'
		]
	}


EPILOG = '''
Examples:

    python %(prog)s -m prod -f LATDF119AAA

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

	df = pd.DataFrame(columns=headers)

	results = {}
	tmp_dir = 'matrix_files'
	os.mkdir(tmp_dir)
	download_file(mfinal_obj, tmp_dir)

	# get dataset-level metadata
	ds_results = report_dataset(mfinal_obj, mfinal_obj['dataset'])

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
		row_to_add = pd.Series(values_to_add, name=mxr_acc)
		df = df.append(row_to_add)
		download_file(mxr, tmp_dir)

		# NEED TO ADD CELL METADATA & VALUES TO THE NEW ANNDATA OBJ

	# NEED TO ADD A FIELD NAME CONVERSION FOR A HANDFUL OF FIELDS TO MEET CXG REQS
	# NEED TO CHANGE ENSEMBL IDS TO GENE_SYMBOLS
	# NEED TO ADD FINAL MX VALUES TO THE NEW ANNDATA OBJ
	# NEED TO ADD DATASET-LEVEL METADATAA TO THE NEW ANNDATA OBJ

	# print the fields into a report
	df.to_csv('temp.csv')
	print(ds_results)

	shutil.rmtree(tmp_dir)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
