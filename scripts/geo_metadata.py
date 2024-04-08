import argparse
import lattice
import pandas as pd
import sys
from urllib.parse import urljoin
import requests
import numpy as np
import time
import gc

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

LIBRARY_METADATA = {
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
		'donors.living_at_sample_collection',
		'donors.donor_id'
		],
	'suspension': [
		'cell_depletion_factors',
		'enrichment_factors',
		'treatment_summary',
		'preservation_method'
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

RAW_MATRIX = [
	'normalized',
	'value_scale',
	'feature_keys',
	'cellranger_assay_chemistry',
	'background_barcodes_included',
	'assembly',
	'genome_annotation',
	'file_format',
	'value_units',
	'assays'
]


# Mapping of field name (object_type + "_" + property) and what needs to be in the final cxg h5ad
PROP_MAP = {
	'library_donors_age_display': 'age',
	'library_donors_sex': 'sex',
	'library_donors_donor_id': 'donor_id',
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
UNREPORTED_VALUE = ''
output_dir = lattice.create_subdirectory('geo_metadata')


def gather_objects(input_object, start_type=None):
	"""
	Gather all objects up the experimental graph, assuming that there is either a suspension or tissue section
	"""
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
			field_lst = ['aliases','protocol','donors','derived_from']
			obj = lattice.get_report('Library', '&@id={}'.format(i), field_lst, connection)[0]
			time.sleep(2)
			libraries.append(obj)
			for o in obj['derived_from']:
				if o.get('uuid') not in susp_ids:
					if 'Suspension' in o['@type']:
						suspensions.append(o)
						susp_ids.append(o.get('uuid'))

		# Will be assuming that raw matrix files are derived from raw sequence files, or else should not run geo script
		for s in input_object['derived_from']:
			obj_type, filter_url = lattice.parse_ids([s])
			field_lst = LIBRARY_METADATA['sequence_file']+['uuid','derived_from', 'accession']
			s_obj = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
			time.sleep(2)
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
			time.sleep(2)
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


def gather_rawmatrices(dataset):
	"""
	Gather raw matrices by object type associated with dataset with following filtering:
	 - 'background_barcodes_included' to select for filtered matrix from CR output
	 - 'status' = 'in progress' to select for valid objects
	"""
	my_raw_matrices = []
	unfiltered_matrices = []

	field_lst = RAW_MATRIX+['@id','accession', 'libraries', 'derived_from', 's3_uri', 'status', 'md5sum']
	filter_url = '&dataset=/datasets/{}/&background_barcodes_included=0&status=in+progress'.format(dataset)
	obj_type = 'RawMatrixFile'
	all_raw_matrices = lattice.get_report(obj_type,filter_url,field_lst,connection)
	for mtx in all_raw_matrices:
		if not mtx.get('background_barcodes_included'):
			my_raw_matrices.append(mtx)
		else:
			unfiltered_matrices.append(mtx)
	return my_raw_matrices


def get_dataset(dataset, geo_study):
	"""
	Get dataset object to get dataset metadata for geo study metadata. Also gather all files for a quick sanity check of file count
	"""
	dataset_id = '/datasets/{}/'.format(dataset)
	field_lst = ['dataset_title', 'description', 'libraries', 'donor_count', 'corresponding_contributors', 'internal_contact', 'original_files']
	obj_type, filter_url = lattice.parse_ids([dataset_id])
	results = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
	title = results.get('dataset_title','')
	desc = results.get('description','')

	all_tissues = []
	all_assays = []
	for lib in results.get('libraries'):
		for tissue in lib.get('biosample_summary'):
			if tissue not in all_tissues:
				all_tissues.append(tissue)
		all_assays.append(lib.get('assay'))
	assays_uniq = list(set(all_assays))
	if len(assays_uniq) == 1:
		design = "{} of {} from {} human donors".format(assays_uniq[0], ', '.join(all_tissues), results.get('donor_count'))
	else:
		design = "{} of {} from {} human donors".format(' and '.join(assays_uniq), ', '.join(all_tissues), results.get('donor_count'))
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
	return(results.get('original_files'))


def get_value(obj, prop):
	"""
	Get property value for given object, can only traverse embedded objects that are embedded 2 levels in
	"""
	path = prop.split('.')
	if len(path) == 1:
		return obj.get(prop, UNREPORTED_VALUE)
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2, UNREPORTED_VALUE) for i in obj[key1]]
			return list(set(values))

		elif obj.get(key1):
			value = obj[key1].get(key2, UNREPORTED_VALUE)
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
			return obj.get(key1,UNREPORTED_VALUE)
	elif len(path) == 3:
		key1 = path[0]
		key2 = path[1]
		key3 = path[2]
		if isinstance(obj.get(key1), list):
			embed_objs = obj.get(key1, UNREPORTED_VALUE)
			values = []
			for embed_obj in embed_objs:
				if isinstance(embed_obj.get(key2), list):
					values += [k.get(key3, UNREPORTED_VALUE) for k in embed_obj[key2]]
				else:
					values += embed_obj[key2].get(key3, UNREPORTED_VALUE)
			return list(set(values))
		# Will need to revisit cell culture and organoid values 
		elif obj.get(key1):
			embed_obj = obj.get(key1, UNREPORTED_VALUE)
			if isinstance(embed_obj.get(key2, UNREPORTED_VALUE), list):
				return [v.get(key3, UNREPORTED_VALUE) for v in embed_obj[key2]]
			elif embed_obj.get(key2, UNREPORTED_VALUE) == UNREPORTED_VALUE:
				return UNREPORTED_VALUE
			else:
				return embed_obj[key2].get(key3, UNREPORTED_VALUE)
		else:
			return obj.get(key1, UNREPORTED_VALUE)
	else:
		return 'unable to traverse more than 2 embeddings'


def gather_metadata(obj_type, properties, values_to_add, objs):
	"""
	Gather object metadata, convert property name to desired field names for geo metadata spreadsheet
	"""
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			if len(value)==1:
				value = value[0]
			else:
				value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = PROP_MAP.get(latkey, latkey)
		values_to_add[key] = value
	if obj_type=='sequencing_run':
		if values_to_add['read_1'] == UNREPORTED_VALUE and values_to_add['read_1N'] != UNREPORTED_VALUE:
			values_to_add['read_1'] = values_to_add['read_1N']
		if values_to_add['read_2'] == UNREPORTED_VALUE and values_to_add['read_2N'] != UNREPORTED_VALUE:
			values_to_add['read_2'] =  values_to_add['read_2N']
		del values_to_add['read_1N']
		del values_to_add['read_2N']


def get_filename(uri):
	"""
	Return last element if split by '/'
	"""
	l = uri.split('/')
	return l[len(l)-1]


def format_protocols(field, key, value, protocols_results):
	"""
	Format protocol metadata to desired geo spreadsheet format
	"""
	if field in protocols_results.keys():
		protocols_results[field] += "; "+key+':'+str(value)
	else:
		protocols_results[field] = key+':'+str(value)
	return protocols_results


def calculate_protocols(geo_protocols, geo_samples):
	"""
	For single values in columns in 'raw_matrix' fields, move over to protocols dataframe
	"""
	protocols_input = {}
	protocols_results = {}
	for c in RAW_MATRIX:
		col = "raw_matrix_"+c
		if col == 'raw_matrix_cellranger_assay_chemistry':
			continue
		if len(geo_samples[col].unique())==1:
			protocols_input[col]=geo_samples[col].unique()[0]
			geo_samples.drop(columns=[col], inplace=True)
	for k,v in protocols_input.items():
		if k in ['raw_matrix_assembly','raw_matrix_genome_annotation']:
			protocols_results = format_protocols('genome build/assembly',k,v,protocols_results)
		elif k in ['raw_matrix_normalized','raw_matrix_value_scale','raw_matrix_background_barcodes_included',\
				'raw_matrix_feature_keys','raw_matrix_cellranger_assay_chemistry']:
			protocols_results = format_protocols('data processing step',k,v,protocols_results)
		elif k in ['raw_matrix_file_format','raw_matrix_value_units']:
			protocols_results = format_protocols('processed data files format and content',k,v,protocols_results)
		elif k == 'raw_matrix_assays':
			protocols_results = format_protocols('library strategy',k,v,protocols_results)
	geo_protocols = pd.DataFrame(protocols_results.items())
	return geo_protocols


def gather_pooled_metadata(obj_type, properties, values_to_add, objs):
	"""
	Gather metadata for pooled objects, which will be a list within double quotes
	"""
	dev_list = []
	for prop in properties:
		if prop == 'ethnicity':
			values_df = pd.DataFrame()
			latkey = (obj_type + '_' + prop).replace('.','_')
			key = PROP_MAP.get(latkey, latkey)
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
				elif len(set(ethnicity_list)) == len(ethnicity_list):
					value = ethnicity_list[0]
					values_df.loc[key,ident] = value
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
			key = PROP_MAP.get(latkey, latkey)
			value_str = [str(i) for i in value]
			value_set = set(value_str)
			if len(value_set) > 1:
				values_to_add[key] = '{}'.format(','.join(value_str))
			else:
				values_to_add[key] = next(iter(value_set))


def main(dataset):
	geo_study = pd.DataFrame()
	geo_samples = pd.DataFrame()
	geo_protocols = pd.DataFrame()
	geo_sequences_list = []
	geo_md5 = pd.DataFrame()
	fastq_meta_list = []
	fastq_meta = pd.DataFrame()
	all_s3_uri = []
	matrix_meta = pd.DataFrame()

	# Get dataset metadata
	all_files = get_dataset(dataset, geo_study)
	print("Total files: {}".format(len(all_files)))

	derived_seq_files = []
	derived_raw_matrix = []
	mxraws_output = gather_rawmatrices(dataset)
	mxraws = mxraws_output
	time.sleep(3)
	print("Total raw matrix files: {}".format(len(mxraws)))

	for mxr in mxraws:
		# get all of the objects necessary to pull the desired metadata

		mxr_acc = mxr['accession']
		print("Processing {}".format(mxr_acc))
		derived_raw_matrix.append(mxr_acc)
		relevant_objects = gather_objects(mxr)
		values_to_add = {} #library metadata for each raw matrix
		all_runs = {} #all fastq files for library
		matrix_to_add = {}
		mxr_to_add = {}

		gather_metadata('raw_matrix', RAW_MATRIX, matrix_to_add, [mxr])
		time.sleep(3)

		for obj_type in LIBRARY_METADATA.keys():
			objs = relevant_objects.get(obj_type, [])
			# For sequencing run objects, add related sequences to geo_sequences and additional metadata to runs_to_add which are added to all_runs
			if obj_type == 'sequencing_run':
				i = 1
				for o in objs:
					runs_to_add = {}
					gather_metadata(obj_type, LIBRARY_METADATA[obj_type], runs_to_add, [o])
					alias = runs_to_add.get('sequencing_run_aliases')
					platform = runs_to_add.get('platform')
					runs_to_add.pop('sequencing_run_aliases')
					runs_to_add.pop('platform')
					runs_to_add.pop('instrument model', None)
					geo_sequences_list.append(runs_to_add)
					runs_to_add['instrument model'] = platform
					for k,v in runs_to_add.items():
						new_k = k+"_run"+str(i)
						all_runs[new_k] = v
					i+=1
			# Go through raw sequence file objects to store metadata in fastq_meta
			elif obj_type == 'sequence_file':
				for o in objs:
					fastq_to_add = {}
					gather_metadata(obj_type, LIBRARY_METADATA[obj_type], fastq_to_add, [o])
					fastq_to_add['filename'] = get_filename(fastq_to_add['filename'])
					fastq_to_add['new_index'] = o.get('accession')
					fastq_meta_list.append(fastq_to_add)
					derived_seq_files.append(o.get('accession'))
					all_s3_uri.append(o.get('s3_uri'))
			else:
				if len(objs)>1:
					gather_pooled_metadata(obj_type, LIBRARY_METADATA[obj_type], values_to_add, objs)
				gather_metadata(obj_type, LIBRARY_METADATA[obj_type], values_to_add, objs)
		# Get mxr metadata
		mxr_to_add['processed file name'] = get_filename(mxr.get('s3_uri'))
		mxr_to_add['processed checksum'] = mxr.get('md5sum')
		geo_md5 = pd.concat([geo_md5, pd.DataFrame(mxr_to_add, index=[geo_md5.shape[0]])])
		values_to_add['processed data file'] = get_filename(mxr.get('s3_uri'))
		all_s3_uri.append(mxr.get('s3_uri'))
		values_to_add.update(matrix_to_add)
		values_to_add.update(all_runs)
		alias = values_to_add.get('library_aliases').split(':')[1]
		geo_samples = pd.concat([geo_samples, pd.DataFrame(values_to_add, index=[alias])])

	gc.collect()

	# Merge fastq metadata to raw sequence files
	# Add sequences to geo_md5
	print("MERGING FASTQ METADATA")
	fastq_meta = pd.DataFrame(fastq_meta_list)
	fastq_meta.set_index('new_index', inplace=True)
	geo_sequences = pd.DataFrame(geo_sequences_list)
	geo_seq_md5 = pd.DataFrame()
	for col in geo_sequences.columns.to_list():
		if len(geo_sequences[col].unique())==1 and geo_sequences[col].unique()[0]==UNREPORTED_VALUE:
			geo_sequences.drop(columns=[col], inplace=True)
		elif col == 'sequencing_run_aliases':
			geo_sequences.drop(columns=[col], inplace=True)
		else:
			geo_sequences = geo_sequences.merge(fastq_meta, left_on=col, right_index=True, how='left')
			geo_subset = pd.DataFrame(geo_sequences[['filename', 'sequence_file_md5sum']])
			geo_subset['filename'].replace('', np.nan, inplace=True)
			geo_subset.dropna(subset=['filename'], inplace=True)
			geo_seq_md5 = pd.concat([geo_seq_md5, geo_subset], ignore_index=True)
			geo_sequences.rename(columns={'filename': col+'_filename'}, inplace=True)
			geo_sequences.drop(columns=[col,'sequence_file_md5sum'], inplace=True)
	geo_sequences.rename(columns={'read_1_filename':'filename 1', 'read_2_filename':'filename 2', 'index_1_filename':'filename 3', 'index_2_filename':'filename 4'}, inplace=True)
	if 'filename 4' in geo_sequences.columns:
		geo_sequences = geo_sequences[['filename 1','filename 2', 'filename 3', 'filename 4']]
	elif 'filename 3' in geo_sequences.columns:
		geo_sequences = geo_sequences[['filename 1','filename 2', 'filename 3']]
	else:
		geo_sequences = geo_sequences[['filename 1','filename 2']]
	geo_md5 = pd.concat([geo_seq_md5,geo_md5], axis=1)
	geo_md5.rename(columns={'processed file name':'filename','processed checksum':'checksum'})
	
	# Replace raw sequence file accession with fastq file name
	print("REPLACING SEQUENCE FILE ACCESSION WITH FILE")
	for col in geo_samples.columns.to_list():
		if len(geo_samples[col].unique())==1 and geo_samples[col].unique()[0]==UNREPORTED_VALUE:
			geo_samples.drop(columns=[col], inplace=True)
		elif col=='library_aliases':
			geo_samples.drop(columns=[col], inplace=True)
	geo_samples.replace(fastq_meta['filename'].to_dict(), inplace=True)

	print("CALCULATING PROTOCOLS")
	# For single values in columns in 'raw_matrix' fields, move over to protocols dataframe
	geo_protocols = calculate_protocols(geo_protocols, geo_samples)

	# Collapse platform
	collapse = []
	for c in geo_samples.columns:	
		if c.startswith('instrument model'):
			collapse.append(c)
	geo_samples['instrument model'] = geo_samples[collapse].stack().groupby(level=0).apply(lambda x: [i for i in x.unique() if i != UNREPORTED_VALUE])
	geo_samples.drop(columns=collapse, inplace=True)
	ordered_cols = [c for c in geo_samples.columns if not c.startswith(('read_','index_'))] + [c for c in geo_samples.columns if c.startswith(('read_','index_'))]
	geo_samples = geo_samples[ordered_cols]

	# Write to files
	# all_df = [geo_study,geo_samples,geo_sequences]
	with open(output_dir + '/' + dataset+"_metadata.csv",'w') as f:
		f.write("STUDY\n")
		geo_study.to_csv(f, header=False, index=False)
		f.write('\nSAMPLES\n')
		geo_samples.to_csv(f)
		f.write('\nPROTOCOLS\n')
		geo_protocols.to_csv(f, header=False, index=False)
		f.write('\nPAIRED-END EXPERIMENTS\n')
		geo_sequences.to_csv(f, index=False)

	with open(output_dir + '/' + dataset+"_s3_uri.csv", "w") as f:
		f.write('\n'.join(all_s3_uri))	

	with open(output_dir + '/' + dataset+"_md5sum.csv",'w') as f:
		f.write('RAW FILES,,PROCESSED DATA FILES\n')
		geo_md5.to_csv(f, index=False)


args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.dataset)







