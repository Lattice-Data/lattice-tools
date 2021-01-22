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
		'suspension_type',
		'percent_cell_viability'
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
		'references.preprint_doi',
		'references.publication_doi',
		'urls'
		],
	'final_matrix': [
		'genome_annotation',
		'value_scale',
		'value_units',
		'normalized',
		'normalization_method'
		]
	}

prop_map = {
	'sample_biosample_ontology_term_name': 'tissue',
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'library_protocol_title': 'assay',
	'library_protocol_term_id': 'assay_ontology_term_id',
	'donor_diseases_term_name': 'disease',
	'donor_diseases_term_id': 'disease_ontology_term_id',
	'donor_sex': 'sex',
	'donor_ethnicity_term_name': 'ethnicity',
	'donor_ethnicity_term_id': 'ethnicity_ontology_term_id',
	'donor_life_stage': 'development_stage',
	'donor_life_stage_term_id': 'development_stage_ontology_term_id',
	'matrix_genome_annotation': 'reference_annotation_version',
	'dataset_references_publication_doi': 'publication_doi',
	'dataset_references_preprint_doi': 'preprint_doi'
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
	lib_ids = raw_matrix_file['libraries']
	libraries = []
	susp_ids = []
	suspensions = []
	prep_susp_ids = []
	prepooled_susps = []
	sample_ids = []
	samples = []
	donor_ids = []
	donors = []

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
		'suspension': suspensions,
		'library': libraries
		}
	if prepooled_susps:
		objs['prepooled_suspension'] = prepooled_susps
		objs['pooled_suspension'] = objs.pop('suspension')
	return objs


def get_value(obj, prop):
	unreported_value = 'unknown'
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


def gather_metdata(obj_type, properties, values_to_add, mxr_acc, objs):
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = value


def gather_pooled_metadata(obj_type, properties, values_to_add, mxr_acc, objs):
	for prop in properties:
		value = set()
		for obj in objs:
			v = get_value(obj, prop)
			value.add(v)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = 'multiple {}s ({})'.format(obj_type, ','.join(value))


def report_dataset(donor_objs, matrix, dataset):
	ds_results = {}
	ds_obj = lattice.get_object(dataset, connection)
	for prop in dataset_metadata['dataset']:
		value = get_value(ds_obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		if value != 'unknown':
			latkey = 'dataset_' + prop.replace('.','_')
			key = prop_map.get(latkey, latkey)
			ds_results[key] = value
	for prop in dataset_metadata['final_matrix']:
		value = get_value(matrix, prop)
		if isinstance(value, list):
			value = ','.join(value)
		if value != 'unknown':
			latkey = 'matrix_' + prop.replace('.','_')
			key = prop_map.get(latkey, latkey)
			ds_results[key] = value

	#org_id = set()
	#org_name = set()
	#for obj in donor_objs:
	#	org_id.add(obj['organism']['taxon_id'])
	#	org_name.add(obj['organism']['scientific_name'])
	#ds_results['organism_ontology_term_id'] = ','.join(org_id)
	#ds_results['organism'] = ','.join(org_name)

	if ds_results.get('publication_doi') and ds_results.get('preprint_doi'):
		del ds_results['preprint_doi']
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


def remove_consistent(df, ds_results):
	# if all values are equal for any metadata field, move those from the cell metadata to the dataset metadata
	to_drop = []
	for label, content in df.items():
		if content.nunique() == 1:
			ds_results[label] = content[0]
			to_drop.append(label)
	df = df.drop(columns=to_drop)
	return df


# Takes ds_results and converts information to cxg fields
# WILL NEED CORPORA VERSION INFORMATION AS INPUT, CURRENTLY HARDCODED
def organize_uns_data(ds_results):
	cxg_uns = {}
	cxg_uns['version'] = {}
	cxg_uns['version']['corpora_schema_version'] = '1.1.0'
	cxg_uns['version']['corpora_encoding_version'] = '0.2.0'
	cxg_uns['title'] = ds_results['dataset_award_title']
	cxg_uns['layer_descriptions'] = {}
	# Summarize normalization, making assumption that .X is normalized and .raw.X is raw
	# ->-> NEED TO CHANGE LOGIC BASED ON WHAT CXG SAYS ON WHAT TO PUT WHERE
	cxg_uns['layer_descriptions']['X'] = '{} counts; {} scaling; normalized using {}'.\
		format(ds_results['matrix_value_units'], ds_results['matrix_value_scale'], ds_results['matrix_normalization_method'])
	cxg_uns['layer_descriptions']['.raw.X'] = 'raw'
	return cxg_uns


# Recreated the final matrix ids, also checking to see if '-1' was removed from original cell identifier
def concatenate_cell_id(cell_mapping, mxr_acc, raw_obs_names, mfinal_cells):
	new_ids = []
	flag_removed = False
	cell_mapping_dct = {}
	for mapping_dict in cell_mapping['label_mappings']:
		cell_mapping_dct[mapping_dict['raw_matrix']] = mapping_dict['string']

	for final_id in mfinal_cells:
		if not re.search('[AGCT]+-1', final_id):
			flag_removed = True
	for id in raw_obs_names:
		if flag_removed:
			id = id.replace('-1', '')
		if cell_mapping['location'] == 'prefix':
			new_ids.append(cell_mapping_dct[mxr_acc]+id)
		elif cell_mapping['location'] == 'suffix':
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
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = prop_map.get(latkey, latkey)
			headers.append(key)

	df = pd.DataFrame(columns=headers)

	results = {}
	tmp_dir = 'matrix_files'
	os.mkdir(tmp_dir)
	download_file(mfinal_obj, tmp_dir)

	# Create initial cell metadata will desired cellIDs and library cell_mapping_identifiers to allow for merging of metadata
	# Get list of cell identifiers
	file_url = mfinal_obj['s3_uri']
	file_ext = file_url.split('.')[-1]
	mfinal_local_path = '{}/{}.{}'.format(tmp_dir, mfinal_obj['accession'], file_ext)
	mfinal_adata = None
	if mfinal_obj['file_format'] == 'hdf5' and re.search('h5ad$', mfinal_local_path):
		mfinal_adata = sc.read_h5ad(mfinal_local_path)
	elif mfinal_obj['file_format'] == 'rds':
		sys.exit('Cannot read from rds {}'.format(final_local_path))
	else:
		sys.exit('Do not recognize file format or exention {} {}'.format(mfinal_obj['file_format'], mfinal_local_path))
	mfinal_cell_identifiers = list(mfinal_adata.obs_names)

	cxg_adata_lst = []

	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	for mxr in mxraws:
		mxr_acc = mxr['accession']
		# get all of the objects necessary to pull the desired metadata
		relevant_objects = gather_objects(mxr)
		values_to_add = {}
		for obj_type in cell_metadata.keys():
			objs = relevant_objects.get(obj_type, [])
			if len(objs) == 1:
				gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, mxr_acc, objs)
			elif len(objs) > 1:
				gather_pooled_metadata(obj_type, cell_metadata[obj_type], values_to_add, mxr_acc, objs)
		if relevant_objects.get('prepooled_suspension'):
			for obj_type in ['prepooled_suspension', 'pooled_suspension']:
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					gather_metdata(obj_type, cell_metadata['suspension'], values_to_add, mxr_acc, objs)
				elif len(objs) > 1:
					gather_pooled_metadata(obj_type, cell_metadata['suspension'], values_to_add, mxr_acc, objs)
		row_to_add = pd.Series(values_to_add, name=mxr['@id'])
		df = df.append(row_to_add)
		download_file(mxr, tmp_dir)

		# Add new anndata to list of final anndatas
		# Should add error checking to make sure all matrices have the same number of vars
		local_path = '{}/{}.h5'.format(tmp_dir,mxr_acc)
		adata_raw = sc.read_10x_h5(local_path)
		adata_raw.var_names_make_unique()

		# Recreate cell_ids and subset raw matrix and add mxr_acc into obs
		concatenated_ids = concatenate_cell_id(mfinal_obj['cell_mapping'], mxr['@id'], adata_raw.obs_names, mfinal_cell_identifiers)
		adata_raw.obs_names = concatenated_ids
		overlapped_ids = list(set(mfinal_cell_identifiers).intersection(concatenated_ids))
		adata_raw = adata_raw[overlapped_ids]
		adata_raw.obs['raw_matrix_accession'] = [mxr['@id']]*len(overlapped_ids)
		cxg_adata_lst.append(adata_raw)

	# get dataset-level metadata
	ds_results = report_dataset(relevant_objects['donor'], mfinal_obj, mfinal_obj['dataset'])

	# Concatenate all anndata objects in list and load parameters
	# NEED TO ADD A FIELD NAME CONVERSION FOR A HANDFUL OF FIELDS TO MEET CXG REQS
	cxg_adata_raw = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None)
	if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
		sys.exit('The number of cells do not match between final matrix and cxg h5ad.')
	cxg_uns = organize_uns_data(ds_results)
	cxg_obs = pd.merge(cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='inner')
	cxg_obsm = get_embeddings(mfinal_adata)
	cxg_adata = ad.AnnData(mfinal_adata.X, obs=cxg_obs, obsm=cxg_obsm, var=mfinal_adata.var, uns=cxg_uns)
	cxg_adata.raw = cxg_adata_raw

	# NEED TO CHANGE ENSEMBL IDS TO GENE_SYMBOLS

	# print the fields into a report
	results_file = "final_cxg.h5ad"
	cxg_adata.write(results_file)
	print(ds_results)
	print(cxg_adata.uns)

	shutil.rmtree(tmp_dir)

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.file)
