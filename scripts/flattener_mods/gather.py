import sys
import pandas as pd
import numpy as np
import logging
import flattener_mods.constants as constants

# Backtracking to scripts folder to import lattice.py
sys.path.insert(0, '../')
import lattice
sys.path.pop(0)

# Attaching logger to Flattener logger
logger = logging.getLogger(__name__)

def gather_rawmatrices(derived_from, connection):
	'''
	Utilize lattice parse_ids and get_report to get all the raw matrix objects at once

	:param List[str] derived_from: Raw matrices up the experimental graph which the final matrix object is derived from
	:param obj connection: Information for connecting to lattice database

	:returns List[dict] my_raw_matrices: List of raw matrices which the final matrix object is derived from
	'''
	new_derived_from = []
	field_lst = ['@id','accession','s3_uri','genome_annotation','libraries','derived_from']
	
	obj_type, filter_lst = lattice.parse_ids(derived_from)
	
	# If object type of the original derived_from ids is not raw matrix file, go another layer down
	if obj_type != 'RawMatrixFile':
		objs = lattice.get_report(obj_type,filter_list,['derived_from'],connection)
		for obj in objs:
			new_derived_from.append(obj['derived_from'])
		obj_type, filter_lst = lattice.parse_ids(new_derived_from)
		
	my_raw_matrices = lattice.get_report(obj_type,filter_lst,field_lst,connection)

	return my_raw_matrices



def gather_objects(input_object, mfinal_obj, connection, start_type=None):
	'''
	Gather all objects up the experimental graph, assuming that there is either a suspension or tissue section

	:param obj input_object: Object for which objects up the experimental graph will be gathered
	:param obj mfinal_obj: Processed Matrix File object
	:param obj connection: Information for connecting to lattice database
	:param str start_type: Starting object type of input_object, if not specified defaults to None

	:returns dict objs: Dictionary with metadata information of all objects up the experimental graph from input_object

	'''
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
				logger.error("ERROR: The flattener cannot flatten multimodal ProcessedMatrixFile")
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



def gather_metdata(obj_type, properties, values_to_add, objs, connection):
	'''
	Gather object metadata, convert property name to cxg required field names

	:prop str obj_type: Object type of object in objs
	:prop List[str] properties: List of properties for which metadata needs to be gathered
	:prop dict values_to_add: Dictionary containing properties as keys and metadata as values
	:prop List[obj] objs: List containing singular obj for which metadata is being gathered
	:prop obj connection: Information for connecting to lattice database

	:returns dict values_to_add: Modified values_to_add dictionary with metadata to be returned
	'''
	obj = objs[0]
	for prop in properties:
		value = get_value(obj,prop)
		if prop == 'family_medical_history':
			if value != 'unknown':
				for history in value:
					ontology = lattice.get_object(history.get('diagnosis'), connection)
					key = 'family_history_' + str(ontology.get('term_name')).replace(' ','_')
					values_to_add[key] = history.get('present')
		elif prop == 'ethnicity':
			ethnicity_list = []
			if value != None:
				for ethnicity_dict in value:
					if ethnicity_dict.get('term_id') == 'NCIT:C17998':
						ethnicity_list.append('unknown')
					else:
						ethnicity_list.append(ethnicity_dict.get('term_id'))
				ethnicity_list.sort()
				value = ','.join(ethnicity_list)
				latkey = (obj_type + '_' + prop).replace('.','_')
				key = constants.PROP_MAP.get(latkey, latkey)
				values_to_add[key] = value
		elif prop == 'cell_ontology.term_id':
			if value == 'NCIT:C17998':
				value = 'unknown'
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = constants.PROP_MAP.get(latkey, latkey)
			values_to_add[key] = value
		else:
			if isinstance(value, list):
				value = ','.join(value)
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = constants.PROP_MAP.get(latkey, latkey)
			values_to_add[key] = value

	return values_to_add



def gather_pooled_metadata(obj_type, properties, values_to_add, objs, connection):
	'''
	Gather metadata for pooled objects
	For required cxg fields, these need to be a single value, so development_stage_ontology_term_id needs to be a commnon slim

	:prop str obj_type: Object type of object in objs
	:prop List[str] properties: List of properties for which metadata needs to be gathered
	:prop dict values_to_add: Dictionary containing properties as keys and metadata as values
	:prop List[obj] objs: List containing singular obj for which metadata is being gathered
	:prop obj connection: Information for connecting to lattice database

	:returns dict values_to_add: Modified values_to_add dictionary with metadata to be returned

	'''
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
			key = constants.PROP_MAP.get(latkey, latkey)
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
			for index, row in values_df.iterrows():
				values_to_add[index] = str(row[0])
		else:
			value = list()
			for obj in objs:
				v = get_value(obj, prop)
				if prop == 'summary_development_ontology_at_collection.development_slims':
					dev_list.append(v)
				if prop == 'cell_ontology.term_id':
					if v == 'NCIT:C17998':
						v = 'unknown'
					value.append(v)
				if isinstance(v, list):
					value.extend(v)
				else:
					value.append(v)
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = constants.PROP_MAP.get(latkey, latkey)
			value_str = [str(i) for i in value]
			value_set = set(value_str)
			cxg_fields = ['disease_ontology_term_id', 'organism_ontology_term_id',\
							 'sex', 'tissue_ontology_term_id', 'development_stage_ontology_term_id']
			if len(value_set) > 1:
				if key in cxg_fields:
					if key == 'development_stage_ontology_term_id':
						dev_in_all = list(set.intersection(*map(set, dev_list)))
						if dev_in_all == []:
							logger.error('ERROR: There is no common development_slims that can be used for development_stage_ontology_term_id')
							sys.exit("ERROR: There is no common development_slims that can be used for development_stage_ontology_term_id")
						else:
							obj = lattice.get_report('OntologyTerm','&term_name='+dev_in_all[0], ['term_id'], connection)
							values_to_add[key] = obj[0].get('term_id')
					elif key == 'sex':
						values_to_add[key] = 'unknown'
					else:
						logger.error('ERROR: Cxg field is a list')
						sys.exit("ERROR: Cxg field is a list")
				else:
					values_to_add[key] = 'pooled [{}]'.format(','.join(value_str))
			else:
				values_to_add[key] = next(iter(value_set))

	return values_to_add



def get_value(obj, prop):
	'''
	Get property value for given object, can only traverse embedded objects that are embedded 2 levels in

	:prop obj obj: Object for which the property value is desired
	:prop str prop: The property value that is desired

	:return varys depending on what is found when traversing the embeddings
	'''
	path = prop.split('.')
	if len(path) == 1:
		if path[0] == '@type':
			value = obj.get('@type')[0]
			return value
		else:
			return obj.get(prop, constants.UNREPORTED_VALUE)
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2, constants.UNREPORTED_VALUE) for i in obj[key1]]
			return list(set(values))
		elif obj.get(key1):
			value = obj[key1].get(key2, constants.UNREPORTED_VALUE)
			return value
		else:
			return obj.get(key1,constants.UNREPORTED_VALUE)
	elif len(path) == 3:
		key1 = path[0]
		key2 = path[1]
		key3 = path[2]
		if isinstance(obj.get(key1), list):
			embed_objs = obj.get(key1, constants.UNREPORTED_VALUE)
			values = []
			for embed_obj in embed_objs:
				if isinstance(embed_obj.get(key2), list):
					values += [k.get(key3, constants.UNREPORTED_VALUE) for k in embed_obj[key2]]
				else:
					values += embed_obj[key2].get(key3, constants.UNREPORTED_VALUE)
			return list(set(values))
		# Will need to revisit cell culture and organoid values 
		elif obj.get(key1):
			embed_obj = obj.get(key1, constants.UNREPORTED_VALUE)
			if isinstance(embed_obj.get(key2, constants.UNREPORTED_VALUE), list):
				return [v.get(key3, constants.UNREPORTED_VALUE) for v in embed_obj[key2]]
			elif embed_obj.get(key2, constants.UNREPORTED_VALUE) == constants.UNREPORTED_VALUE:
				return constants.UNREPORTED_VALUE
			else:
				return embed_obj[key2].get(key3, constants.UNREPORTED_VALUE)
		else:
			return obj.get(key1, constants.UNREPORTED_VALUE)
	else:
		return 'unable to traverse more than 2 embeddings'