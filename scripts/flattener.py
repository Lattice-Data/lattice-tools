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
		'life_stage_term_id'
		],
	'sample': [
		'uuid',
		'preservation_method',
		'biosample_ontology.term_name',
		'biosample_ontology.term_id'
		],
	'suspension': [
		'uuid',
		'suspension_type',
		'percent_cell_viability'
		],
	'library': [
		'uuid',
		'protocol.title',
		'protocol.term_id'
		]
	}

dataset_metadata = {
	'dataset': [
		'references.preprint_doi',
		'references.publication_doi',
		'urls'
		],
	'final_matrix': [
		'description',
		'genome_annotation'
		]
	}

annot_fields = [
	'cell_ontology.term_name',
	'cell_ontology.term_id',
	'author_cell_type'
]

prop_map = {
	'sample_biosample_ontology_term_name': 'tissue',
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'library_protocol_title': 'assay',
	'library_protocol_term_id': 'assay_ontology_term_id',
	'donor_sex': 'sex',
	'donor_ethnicity_term_name': 'ethnicity',
	'donor_ethnicity_term_id': 'ethnicity_ontology_term_id',
	'donor_life_stage': 'development_stage',
	'donor_life_stage_term_id': 'development_stage_ontology_term_id',
	'matrix_genome_annotation': 'reference_annotation_version',
	'matrix_description': 'title',
	'dataset_references_publication_doi': 'publication_doi',
	'dataset_references_preprint_doi': 'preprint_doi',
	'cell_annotation_author_cell_type': 'author_cell_type',
	'cell_annotation_cell_ontology_term_id': 'cell_type_ontology_term_id',
	'cell_annotation_cell_ontology_term_name': 'cell_type'
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

	layer_descs = {
		'.raw.X': 'raw'
	}
	for layer in matrix.get('layers'):
		label = layer.get('label')
		units = layer.get('value_units', 'unknown')
		scale = layer.get('value_scale', 'unknown')
		norm_meth = layer.get('normalization_method', 'unknown')
		desc = '{} counts; {} scaling; normalized using {}'.format(units, scale, norm_meth)
		layer_descs[label] = desc
	ds_results['layer_descriptions'] = layer_descs

	pub_doi = set()
	for pub in ds_obj['references']:
		for i in pub['identifiers']:
			if i.startswith('doi:'):
				pub_doi.add(i.replace('doi:', 'https://doi.org/'))
	if pub_doi:
		ds_results['doi'] = ','.join(pub_doi)

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


def remove_consistent(df, ds_results):
	# if all values are equal for any metadata field, move those from the cell metadata to the dataset metadata
	to_drop = []
	for label, content in df.items():
		if content.nunique() == 1:
			ds_results[label] = content[0]
			to_drop.append(label)
	df = df.drop(columns=to_drop)

	return df


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


def report_diseases(values_to_add, donor_objs, sample_objs):
	names = set()
	ids = set()
	my_donors = []
	for o in sample_objs:
		dis_objs = o.get('diseases')
		for donor in donor_objs:
			if donor['@id'] in o.get('donors'):
				my_donors.append(donor)
		if dis_objs:
			for do in dis_objs:
				ids.add(do.get('term_id'))
				names.add(do.get('term_name'))
		for o in my_donors:
			dis_objs = o.get('diseases')
			if dis_objs:
				for do in dis_objs:
					ids.add(do.get('term_id'))
					names.add(do.get('term_name'))
		if not ids:
			ids.add('unknown')
		if not names:
			names.add('unknown')

	if len(sample_objs) > 1:
		values_to_add['disease'] = 'multiple samples ({})'.format(','.join(names))
		values_to_add['disease_ontology_term_id'] = 'multiple samples ({})'.format(','.join(ids))
	else:
		values_to_add['disease'] = ','.join(names)
		values_to_add['disease_ontology_term_id'] = ','.join(ids)


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
	if mfinal_obj['file_format'] == 'hdf5' and re.search('h5ad$', mfinal_local_path):
		mfinal_adata = sc.read_h5ad(mfinal_local_path)
	elif mfinal_obj['file_format'] == 'rds':
		sys.exit('Cannot read from rds {}'.format(mfinal_local_path))
	else:
		sys.exit('Do not recognize file format or exention {} {}'.format(mfinal_obj['file_format'], mfinal_local_path))
	mfinal_cell_identifiers = list(mfinal_adata.obs_names)

	cxg_adata_lst = []


	# get the list of matrix files that hold the raw counts corresponding to our Final Matrix
	mxraws = gather_rawmatrices(mfinal_obj['derived_from'])
	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		relevant_objects = gather_objects(mxr)
		values_to_add = {}
		for obj_type in cell_metadata.keys():
			objs = relevant_objects.get(obj_type, [])
			if len(objs) == 1:
				gather_metdata(obj_type, cell_metadata[obj_type], values_to_add, objs)
			elif len(objs) > 1:
				gather_pooled_metadata(obj_type, cell_metadata[obj_type], values_to_add, objs)
		if relevant_objects.get('prepooled_suspension'):
			for obj_type in ['prepooled_suspension', 'pooled_suspension']:
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					gather_metdata(obj_type, cell_metadata['suspension'], values_to_add, objs)
				elif len(objs) > 1:
					gather_pooled_metadata(obj_type, cell_metadata['suspension'], values_to_add, objs)
		report_diseases(values_to_add, relevant_objects['donor'], relevant_objects['sample'])
		row_to_add = pd.Series(values_to_add, name=mxr['@id'])
		df = df.append(row_to_add)
		download_file(mxr, tmp_dir)

		# Add new anndata to list of final anndatas
		local_path = '{}/{}.h5'.format(tmp_dir,mxr_acc)
		adata_raw = sc.read_10x_h5(local_path)
		adata_raw.var_names_make_unique()

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

	print(annot_df)
	print(ds_results)
	# Concatenate all anndata objects in list and load parameters
	# ->->STILL NEED TO ADD CXG SCHEMA VERSION IN UNS
	cxg_adata_raw = cxg_adata_lst[0].concatenate(cxg_adata_lst[1:], index_unique=None)
	cxg_adata_raw = cxg_adata_raw[mfinal_cell_identifiers]
	if cxg_adata_raw.shape[0] != mfinal_adata.shape[0]:
		sys.exit('The number of cells do not match between final matrix and cxg h5ad.')
	cxg_uns = ds_results
	cxg_obsm = get_embeddings(mfinal_adata)

	# Prep obs dataframe
	celltype_col = mfinal_obj['author_cell_type_column']
	cxg_obs = pd.merge(cxg_adata_raw.obs, df, left_on='raw_matrix_accession', right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, mfinal_adata.obs[[celltype_col]], left_index=True, right_index=True, how='left')
	cxg_obs = pd.merge(cxg_obs, annot_df, left_on=celltype_col, right_index=True, how='left')
	cxg_obs.drop(columns=['raw_matrix_accession','batch', celltype_col], inplace=True)

	cxg_adata = ad.AnnData(mfinal_adata.X, obs=cxg_obs, obsm=cxg_obsm, var=mfinal_adata.var, uns=cxg_uns)
	cxg_adata.raw = cxg_adata_raw

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
