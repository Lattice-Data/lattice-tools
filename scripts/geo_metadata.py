import argparse
import anndata as ad
import lattice
import os
import pandas as pd
import shutil
import sys
import scanpy as sc
import re
import subprocess
from urllib.parse import urljoin
import requests
import gspread


EPILOG = '''
Examples:

    python %(prog)s -m local -d LATDS977NKW

For more details:

    python %(prog)s --help
'''

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--dataset', '-d',
                        help='Any identifier for the dataset of interest.')
    parser.add_argument('--mode', '-m',
                        help='The machine to run on.')
    args = parser.parse_args()
    if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit()
    return args

library_metadata = {
	'sample': [
		'biosample_ontology.term_name',
		'diseases.term_name',
		'treatment_summary'
		],
	'library': [
		'aliases',
		'protocol.assay_ontology.term_name',
		'protocol.biological_macromolecule',
		'donors.sex',
		'donors.ethnicity.term_name',
		'donors.age_display',
		'donors.living_at_sample_collection'
		],
	'suspension': [
		'cell_depletion_factors',
		'enrichment_factors'
		],
	'sequencing_run': [
		'aliases',
		'read_1N_file',
		'read_1_file',
		'read_2N_file',
		'read_2_file',
		'i5_index_file',
		'i7_index_file',
		'platform'
	],
	'sequence_file': [
		's3_uri',
		'md5sum'
	]
}



# Mapping of field name (object_type + "_" + property) and what needs to be in the final cxg h5ad
prop_map = {
	'library_donors_age_display': 'age',
	'library_donors_sex': 'sex',
	'library_protocol_assay_ontology_term_name': 'assay',
	'library_protocol_biological_macromolecule': 'molecule',
	'library_donors_ethnicity_term_name': 'ethnicity',
	'sample_biosample_ontology_term_name': 'tissue',
	'sample_diseases_term_name': 'disease',
	'sequencing_run_read_1N_file': 'read_1N',
	'sequencing_run_read_1_file': 'read_1',
	'sequencing_run_read_2N_file': 'read_2N',
	'sequencing_run_read_2_file': 'read_2',
	'sequencing_run_i5_index_file': 'index_2',
	'sequencing_run_i7_index_file': 'index_1',
	'sequencing_run_platform': 'platform',
	'sequence_file_s3_uri': 'filename'
}


# Global variables
unreported_value = ''
geo_study = pd.DataFrame()
geo_samples = pd.DataFrame()
geo_protocols = pd.DataFrame()
geo_sequences = pd.DataFrame()

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
	sequencing_run_ids = []
	sequencing_runs = []
	fastq_ids = []
	fastqs =[]

	if start_type == None:
		for i in lib_ids:
			obj = lattice.get_object(i, connection)
			libraries.append(obj)
			for o in obj['derived_from']:
				if o.get('uuid') not in susp_ids:
					if 'Suspension' in o['@type']:
						suspensions.append(o)
						susp_ids.append(o.get('uuid'))
		for s in input_object['derived_from']:
			s_obj = lattice.get_object(s, connection)
			if s_obj.get('uuid') not in fastq_ids:
				fastqs.append(s_obj)
				fastq_ids.append(s_obj)
			for r in s_obj['derived_from']:
				if r.get('uuid') not in sequencing_run_ids:
					sequencing_runs.append(r)
					sequencing_run_ids.append(r.get('uuid'))

	elif start_type == 'suspension':
		susp_ids = [input_object['uuid']]
		suspensions = [input_object]

	if len(suspensions) > 0:
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
		'sample': samples,
		'suspension': suspensions,
		'sequencing_run': sequencing_runs,
		'sequence_file': fastqs
		}
	if start_type == None:
		objs['library'] = libraries
	if prepooled_susps:
		objs['prepooled_suspension'] = prepooled_susps
		objs['pooled_suspension'] = objs['suspension']

	return objs


# Gather raw matrices by object type and 'background_barcodes_included' to select for filtered matrix from CR output
def gather_rawmatrices(dataset):
	my_raw_matrices = []
	unfiltered_matrices = []
	raw_ids = []
	query_url = urljoin(server, 'search/?type=RawMatrixFile&dataset=/datasets/' + dataset + '/&format=json' + '&limit=all')
	r = requests.get(query_url, auth=connection.auth)
	try:
		r.raise_for_status()
	except requests.HTTPError:
		sys.exit("Error in getting raw matrix files: {}".format(query_url))
	else:
		for raw in r.json()['@graph']:
			raw_ids.append(raw['@id'])
	for identifier in raw_ids:
		obj = lattice.get_object(identifier, connection)
		if obj['@type'][0] == 'RawMatrixFile':
			if not obj.get('background_barcodes_included'):
				my_raw_matrices.append(obj)
			else:
				unfiltered_matrices.append(obj)
		else:
			continue
	return my_raw_matrices, unfiltered_matrices


def get_dataset(dataset):
	global geo_study
	results = lattice.get_object(dataset, connection)
	title = results.get('dataset_title','')
	desc = results.get('description','')

	all_tissues = []
	for lib in results.get('libraries'):
		for tissue in lib.get('biosample_summary'):
			if tissue not in all_tissues:
				all_tissues.append(tissue)
	design = "{} of {} from {} human donors".format(results.get('libraries')[0].get('assay'), ','.join(all_tissues), results.get('donor_count'))
	keys = ['title','summary (abstract)','experimental design']
	vals = [title,desc,design]
	
	if results.get('corresponding_contributors',None):
		con1 = results.get('corresponding_contributors')[0]['title']
		keys.append('contributor')
		vals.append(con1)
		
	if results.get('internal_contact',None):
		obj = lattice.get_object(results.get('internal_contact'), connection)
		con2 = obj.get('first_name','')+' '+obj.get('last_name','')
		keys.append('contributor')
		vals.append(con2)

	geo_study[0] = keys
	geo_study[1] = vals
	return(results.get('files'))


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
def gather_metadata(obj_type, properties, values_to_add, objs):
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			if len(value)==1:
				value = value[0]
			else:
				value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = value
	if obj_type=='sequencing_run':
		if values_to_add['read_1'] == unreported_value and values_to_add['read_1N'] != unreported_value:
			values_to_add['read_1'] = values_to_add['read_1N']
		if values_to_add['read_2'] == unreported_value and values_to_add['read_2N'] != unreported_value:
			values_to_add['read_2'] =  values_to_add['read_2N']
		del values_to_add['read_1N']
		del values_to_add['read_2N']

# Return last element if split by '/'
def get_filename(uri):
	l = uri.split('/')
	return l[len(l)-1]


def main(dataset):
	global geo_study
	global geo_samples
	global geo_protocols
	global geo_sequences
	fastq_meta = pd.DataFrame()
	all_s3_uri = []

	# Get dataset metadata
	all_files = get_dataset(dataset)
	raw_seq_files = [i.split('/')[2] for i in all_files if i.startswith("/raw-sequence-files")]
	raw_matrix_files = [i.split('/')[2] for i in all_files if i.startswith("/raw-matrix-files")]

	derived_seq_files = []
	derived_raw_matrix = []
	derived_unfiltered_matrix = []
	mxraws_output = gather_rawmatrices(dataset)
	mxraws = mxraws_output[0]
	for unfiltered in mxraws_output[1]:
		derived_unfiltered_matrix.append(unfiltered.get('accession'))

	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata
		mxr_acc = mxr['accession']
		derived_raw_matrix.append(mxr_acc)
		relevant_objects = gather_objects(mxr)
		values_to_add = {} #library metadata for each raw matrix
		fastq_to_add = {} #fastq metadata
		runs_to_add = {} #sequencing run metadata
		all_runs = {} #all fastq files for library

		for obj_type in library_metadata.keys():
			objs = relevant_objects.get(obj_type, [])
			if obj_type == 'sequencing_run':
				i = 1
				for o in objs:
					gather_metadata(obj_type, library_metadata[obj_type], runs_to_add, [o])
					alias = runs_to_add.get('sequencing_run_aliases')
					platform = runs_to_add.get('platform')
					runs_to_add.pop('sequencing_run_aliases')
					runs_to_add.pop('platform')
					runs_to_add.pop('instrument model', None)
					geo_sequences = pd.concat([geo_sequences, pd.DataFrame(runs_to_add, index=[alias])])
					runs_to_add['instrument model'] = platform
					for k,v in runs_to_add.items():
						new_k = k+"_run"+str(i)
						all_runs[new_k] = v
					i+=1
			elif obj_type == 'sequence_file':
				for o in objs:
					gather_metadata(obj_type, library_metadata[obj_type], fastq_to_add, [o])
					fastq_to_add['filename'] = get_filename(fastq_to_add['filename'])
					fastq_meta = pd.concat([fastq_meta, pd.DataFrame(fastq_to_add, index=[o.get('accession')])])
					derived_seq_files.append(o.get('accession'))
					all_s3_uri.append(o.get('s3_uri'))
			else:
				if len(objs)>1:
					sys.exit("Cannot handle multiplexed libraries")
				gather_metadata(obj_type, library_metadata[obj_type], values_to_add, objs)
		values_to_add['processed data file'] = get_filename(mxr.get('s3_uri'))
		all_s3_uri.append(mxr.get('s3_uri'))
		values_to_add.update(all_runs)
		alias = values_to_add.get('library_aliases').split(':')[1]
		geo_samples = pd.concat([geo_samples, pd.DataFrame(values_to_add, index=[alias])])

	for col in geo_sequences.columns.to_list():
		if len(geo_sequences[col].unique())==1 and geo_sequences[col].unique()[0]==unreported_value:
			geo_sequences.drop(columns=[col], inplace=True)
		elif col == 'sequencing_run_aliases':
			geo_sequences.drop(columns=[col], inplace=True)
		else:
			geo_sequences = geo_sequences.merge(fastq_meta, left_on=col, right_index=True, how='left')
			geo_sequences.rename(columns={'filename': col+'_filename', 'sequence_file_md5sum': col+'_md5sum'}, inplace=True)
			geo_sequences.drop(columns=[col], inplace=True)

	for col in geo_samples.columns.to_list():
		if len(geo_samples[col].unique())==1 and geo_samples[col].unique()[0]==unreported_value:
			geo_samples.drop(columns=[col], inplace=True)
		elif col=='library_aliases':
			geo_samples.drop(columns=[col], inplace=True)
		elif col.startswith(('read','index')):
			geo_samples = geo_samples.merge(fastq_meta, left_on=col, right_index=True, how='left')
			geo_samples.drop(columns=[col,'sequence_file_md5sum'], inplace=True)
			geo_samples.rename(columns={'filename':col}, inplace=True)

	# Check files with files field in dataset object
	dataset_seq_uniq = [i for i in raw_seq_files if i not in derived_seq_files]
	derived_seq_uniq = [i for i in derived_seq_files if i not in raw_seq_files]
	dataset_matrix_uniq = [i for i in raw_matrix_files if i not in derived_raw_matrix]
	derived_matrix_uniq = [i for i in derived_raw_matrix if i not in raw_matrix_files]
	if len(dataset_seq_uniq)>1:
		print("raw seq only in dataset: {}".format(dataset_seq_uniq))
	elif len(derived_seq_uniq)>1:
		print("raw seq only in derived: {}".format(derived_seq_uniq))
	elif len(dataset_matrix_uniq)>1:
		filtered_uniq = [i for i in dataset_matrix_uniq if i not in derived_unfiltered_matrix]
		if len(filtered_uniq)>1:
			print("matrix only in dataset filtered: {}".format(filtered_uniq))
	elif len(derived_matrix_uniq)>1:
		print("matrix only in derived: {}".format(derived_matrix_uniq))

	# Write to files
	# all_df = [geo_study,geo_samples,geo_sequences]
	with open(dataset+"_metadata.csv",'a') as f:
		geo_study.to_csv(f, header=False, index=False)
		f.write('\n')
		geo_samples.to_csv(f)
		f.write('\n')
		geo_sequences.to_csv(f)
	with open(dataset+"_s3_uri.csv", "w") as f:
		f.write('\n'.join(all_s3_uri))	


args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.dataset)







