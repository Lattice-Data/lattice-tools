import argparse
import hashlib
import crcmod
import json
import lattice
import logging
import os
import requests
import sys
import zlib
from datetime import datetime, timezone
from pint import UnitRegistry
from urllib.parse import urljoin
import DCP_mods.request_to_gcp as request_to_gcp
from DCP_mods.DCP_flatten import tsv_report
from DCP_mods.property_mapping import (
	dcp_versions,
	lattice_to_dcp
)
from DCP_mods.staging_area_validator import (
	dcp_validation
)


EPILOG = '''
Get a Dataset from the Lattice database and put into DCP schema and exported in json files.

Examples:

    python %(prog)s --mode prod --dataset LATDS940AQG

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
    parser.add_argument('--data-only',
                        default=False,
                        action='store_true',
                        help='If true and --update, only data files will be transferred to DCP, no metadata files.')
    parser.add_argument('--metadata-only',
                        default=False,
                        action='store_true',
                        help='If true and --update, only metadata files will be transferred to DCP, no data files.')
    parser.add_argument('--no-validate',
                        default=False,
                        action='store_true',
                        help='If true, nothing will be validated against the DCP schema.')
    parser.add_argument('--validate-only',
                        default=False,
                        action='store_true',
                        help='If true, existing metadata files will be validated against the DCP schema.')
    parser.add_argument('--revision',
                        default=False,
                        action='store_true',
                        help='If true, treated as updating an existing DCP project.  Default is False')
    parser.add_argument('--update',
                        default=False,
                        action='store_true',
                        help='Let the script proceed with the changes.  Default is False')
    args = parser.parse_args()
    return args


dcp_types = {
	'project': 'project',
	'donor_organism': 'biomaterial',
	'specimen_from_organism': 'biomaterial',
	'imaged_specimen': 'biomaterial',
	'cell_line': 'biomaterial',
	'organoid': 'biomaterial',
	'cell_suspension': 'biomaterial',
	'sequence_file': 'file',
	'library_preparation_protocol': 'protocol/sequencing',
	'sequencing_protocol': 'protocol/sequencing',
	'dissociation_protocol': 'protocol/biomaterial_collection',
	'collection_protocol': 'protocol/biomaterial_collection',
	'enrichment_protocol': 'protocol/biomaterial_collection',
	'differentiation_protocol': 'protocol/biomaterial_collection',
	'protocol': 'protocol',
	'supplementary_file': 'file',
	'process': 'process'
}


def test_mapping():
	error_count = 0
	schema_url = urljoin(server, 'profiles/?format=json')
	schemas = requests.get(schema_url, timeout=60).json()
	for k,v in lattice_to_dcp.items():
		schema_props = schemas[k]['properties']
		for prop in v.keys():
			flag = False
			if prop == 'class':
				flag = True
			else:
				lat_prop = v[prop]['lattice']
				if '.' in lat_prop:
					path = lat_prop.split('.')
					if path[0] in schema_props:
						flag = True
					elif schema_props[path[0]].get('linkTo'):
						link_obj = schema_props[path[0]]['linkTo']
						if path[1] in schemas[link_obj]['properties']:
							flag = True
					elif schema_props[path[0]].get('properties'):
						if path[1] in schema_props[path[0]]['properties']:
							flag = True
					elif schema_props[path[0]].get('items'):
						if schema_props[path[0]]['items'].get('properties'):
							if path[1] in schema_props[path[0]]['items']['properties']:
								flag = True
						elif schema_props[path[0]]['items'].get('linkTo'):
							link_obj = schema_props[path[0]]['items']['linkTo']
							if path[1] in schemas[link_obj]['properties']:
								flag = True
					else:
						logging.info(k)
						logging.info(path)
					if len(path) > 2:
						logging.info('-----long path----')
				elif lat_prop in schema_props:
					flag = True
			if flag != True:
				error_count += 1
				logging.info('{} not found in schema for {}'.format(lat_prop,k))
	return error_count


def split_dbxrefs(my_obj, dbxrefs, dbxref_map):
	for i in dbxrefs:
		path = i.split(':')
		prefix = path[0]
		v = path[1]
		if dbxref_map.get(prefix):
			prop = dbxref_map[prefix]
			if '.' in prop:
				path = prop.split('.')
				if len(path) == 2:
					if my_obj.get(path[0]):
						my_obj[path[0]][path[1]] = v
					else:
						my_obj[path[0]] = {}
						my_obj[path[0]][path[1]] = v
			elif my_obj.get(prop):
				my_obj[prop].append(v)
			else:
				my_obj[prop] = [v]


def add_value(my_obj, temp_obj, prop_map, prop):
	lat_prop = prop_map['lattice']
	v = temp_obj[lat_prop]
	if prop_map.get('value_map'):
		v = prop_map['value_map'].get(v)
	if not v:
		add_to_not_incl(lat_prop, v)
	else:
		# lists of objects need to be mapped each item at a time
		if prop_map.get('subprop_map'):
			if isinstance(v, list):
				new_v = []
				for i in v:
					new_i = {}
					for subprop in prop_map['subprop_map'].keys():
						lat_subprop = prop_map['subprop_map'][subprop]['lattice']
						if i.get(lat_subprop):
							add_value(new_i, i, prop_map['subprop_map'][subprop], subprop)
					new_v.append(new_i)
					v = new_v
			elif isinstance(v, dict):
				new_v = {}
				for subprop in prop_map['subprop_map'].keys():
					lat_subprop = prop_map['subprop_map'][subprop]['lattice']
					if v.get(lat_subprop):
						add_value(new_v, v, prop_map['subprop_map'][subprop], subprop)
				v = new_v
		# map and create subobjects
		if '.' in prop:
			path = prop.split('.')
			if len(path) == 2:
				if my_obj.get(path[0]):
					my_obj[path[0]][path[1]] = v
				else:
					my_obj[path[0]] = {}
					my_obj[path[0]][path[1]] = v
			elif len(path) == 3:
				if my_obj.get(path[0]):
					if my_obj[path[0]].get(path[1]):
						my_obj[path[0]][path[1]][path[2]] = v
					else:
						my_obj[path[0]][path[1]] = {}
						my_obj[path[0]][path[1]][path[2]] = v
				else:
					my_obj[path[0]] = {}
					my_obj[path[0]][path[1]] = {}
					my_obj[path[0]][path[1]][path[2]] = v
		# Lattice dbxrefs are split into various DCP fields
		elif prop == 'dbxrefs':
			split_dbxrefs(my_obj, v, prop_map)
		# add a mapped property and its value
		else:
			my_obj[prop] = v


def add_to_not_incl(prop, value):
	ok_to_drop = ['@context','@id','@type','age_display','award','biosample_classification','biosample_ontologies','biosample_ontology',
		'biosample_summary','children','content_md5sum','dataset','demultiplexed_type','derivation_process','description','derived_from',
		'donors','ethnicity','href','i5_index_file','i7_index_file','internal_contact','lab','no_file_available','notes','organism',
		'output_types','protocol','read_1_file','read_1N_file','read_2_file','read_2N_file','read_length_units','reference_annotation',
		'reference_assembly','references','schema_version','sequence_elements','software','title','url','uuid','validated']

	base_prop = prop.split('.')[0]
	if base_prop not in ok_to_drop:
		if isinstance(value, list):
			if not_incl.get(prop):
				for i in value:
					if str(i) not in not_incl[prop]:
						not_incl[prop].append(str(i))
			else:
				not_incl[prop] = []
				for i in value:
					if str(i) not in not_incl[prop]:
						not_incl[prop].append(str(i))			
		elif isinstance(value, dict):
			not_incl[prop] = ['dictionary values']
		else:
			if not_incl.get(prop):
				if str(value) not in not_incl[prop]:
					not_incl[prop].append(str(value))
			else:
				not_incl[prop] = [str(value)]


def file_stats(local_path, obj):
	# get file size
    file_stat = os.stat(local_path)
    obj['file_size'] = file_stat.st_size

    # get the sha256
    sha256_hash = hashlib.sha256()
    with open(local_path,"rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            sha256_hash.update(byte_block)
        obj['sha256'] = sha256_hash.hexdigest()

    # get the crc32c
    with open(local_path, 'rb') as f:
        crc32c_func = crcmod.predefined.Crc('crc-32c')
        while True:
            s = f.read(65536)
            if not s:
                break
            crc32c_func.update(s)
        obj['crc32c'] = crc32c_func.hexdigest().lower()


def get_object(temp_obj, variable_age=False):
	obj_type = temp_obj['@type'][0]

	# if there are ERROR-level audits, flag it to stop the script
	if (temp_obj.get('audit') and temp_obj['audit'].get('ERROR')):
		freq = {}
		for a in temp_obj['audit']['ERROR']:
			if a['category'] in freq:
				freq[a['category']] += 1
			else:
				freq[a['category']] = 1
		for k,v in freq.items():
			print('ERROR audit:{}x {} on {} {}'.format(
				str(v),
				k,
				obj_type,
				temp_obj.get('uuid')))
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to one or more ERROR audits')

	# drop unused properties, flatten subobjects
	temp_obj = flatten_obj(temp_obj)

	remove = set()
	my_obj = {}

	if temp_obj.get('documents'):
		[handle_doc(d) for d in temp_obj['documents'] if d not in handled_docs]

	for prop in lattice_to_dcp[obj_type].keys():
		if prop != 'class':
			lat_prop = lattice_to_dcp[obj_type][prop]['lattice']
			if temp_obj.get(lat_prop) or temp_obj.get(lat_prop) == False:
				remove.add(lat_prop)
				add_value(my_obj, temp_obj, lattice_to_dcp[obj_type][prop], prop)

	dcp_obj_type = lattice_to_dcp[obj_type]['class']

	# remove mapped properties from the object and stash the rest to report later
	for p in remove:
		del temp_obj[p]
	for k,v in temp_obj.items():
		add_to_not_incl(k, v)

	if variable_age:
		new_id = uuid_make(my_obj['provenance']['document_id'] + ''.join(variable_age))
		my_obj['provenance']['document_id'] = new_id
		my_obj['timecourse'] = {
			'value': variable_age[0],
			'unit': {'text': variable_age[1]},
			'relevance': 'donor age at collection'
		}
		#remove the non-HsapDv terms, but leave the development_stage.text
		del my_obj['development_stage']['ontology']
		del my_obj['development_stage']['ontology_label']

	# add the object of mapped properties
	if whole_dict.get(dcp_obj_type):
		whole_dict[dcp_obj_type].append(my_obj)
	else:
		whole_dict[dcp_obj_type] = [my_obj]


def flatten_obj(obj):
	ignore = ['accession','actions','aliases','audit',
		'date_created','original_files','status','submitted_by',
		'superseded_by','supersedes']
	new_obj = {}
	for k,v in obj.items():
		if k not in ignore:
			if isinstance(v, dict) and k not in ['test_results']:
				for x,y in v.items():
					new_obj[k + '.' + x] = y
			elif (isinstance(v, list) and all(isinstance(i, dict) for i in v)):
				new_v = []
				for i in v:
					new_i = flatten_obj(i)
					new_v.append(new_i)
				new_obj[k] = new_v
			else:
				new_obj[k] = v
	return new_obj


def get_derived_from(temp_obj, next_remaining, links):
	uuid = temp_obj['uuid']
	# get identifiers for any object that referenced by this one
	if temp_obj.get('derived_from'):
		if temp_obj.get('age_at_collection'): #these will need their donors split
			variable_age = (temp_obj['age_at_collection'], temp_obj['age_at_collection_units'])
		else:
			variable_age = False

		if isinstance(temp_obj['derived_from'][0], str): #object is not embedded
			if variable_age:
				identifier = [(temp_obj['derived_from'][0], variable_age)]
				next_remaining.update(identifier)
				add_links(temp_obj, tuple(identifier), links)
			else:
				next_remaining.update(temp_obj['derived_from'])
				add_links(temp_obj, tuple(temp_obj['derived_from']), links)
		elif isinstance(temp_obj['derived_from'][0], dict): #object is embedded
			if variable_age:
				identifier = [(temp_obj['derived_from'][0]['@id'], variable_age)]
				next_remaining.update(identifier)
				add_links(temp_obj, tuple(identifier), links)

			else:
				der_fr = ()
				for o in temp_obj['derived_from']:
					next_remaining.add(o['@id'])
					der_fr = der_fr + (o['@id'],)
				add_links(temp_obj, der_fr, links)


def get_links(temp_obj, der_fr, links_dict):
	lat_type = temp_obj['@type'][0]
	out_type = lattice_to_dcp[lat_type]['class']
	out_id = temp_obj['uuid']
	outs = {'output_type': out_type, 'output_id': out_id}
	if links_dict.get(der_fr):
		links_dict[der_fr]['outputs'].append(outs)
	else:
		links_dict[der_fr] = {'outputs': [outs]}


def uuid_make(value):
	s = str(value)
	md5 = hashlib.md5(s.encode('utf-8')).hexdigest()
	return '-'.join([
		md5[0:8],
		md5[8:12],
		md5[12:16],
		md5[16:20],
		md5[20:]
		])


def seq_to_susp(links_dict):
	links = []
	lib_cell_counts = {}
	all_susps = set()
	for seqruns in list(links_dict.keys()):
		l = {}
		ins = []
		susps = []
		protocols = []
		for sr in seqruns:
			# handle the sequencing metadata
			field_lst = ['derived_from','demultiplexed_link','platform']
			sr_obj = lattice.get_report('SequencingRun', f'&@id={sr}', field_lst, connection)[0]
			lib = sr_obj['derived_from'][0]
			field_lst = ['observation_count', 'derived_from'] + \
				[v['lattice'] for v in lattice_to_dcp['Library'].values() if isinstance(v, dict)]
			lib_obj = lattice.get_report('Library', f'&@id={lib}', field_lst, connection)[0]
			lib_cell_counts[lib] = lib_obj.get('observation_count') #grab to sum the project counts

			# see if we need to skip back over a pooling step for demultiplexed sequence data
			if sr_obj.get('demultiplexed_link'):
				obj_type, filter_url = lattice.parse_ids([sr_obj['demultiplexed_link']])
				field_lst = ['uuid']
				prepooled_obj = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
				if obj_type == 'Suspension':
					susps.extend([prepooled_obj['@id']])
					in_type = lattice_to_dcp[obj_type]['class']
					in_id = prepooled_obj['uuid']
					ins.append({'input_type': in_type, 'input_id': in_id})
					seq_method = 'high throughput sequencing, demultiplexing'
				# cannot yet handle multiple Tissues pooled into a single Suspension
				else:
					sys.exit('ERROR: not yet equiped to handle pooling non-Suspensions directly pooled into a Suspension')
			else:
				susps.extend([i['@id'] for i in lib_obj['derived_from']])
				for obj in lib_obj['derived_from']:
					lat_type = obj['@type'][0]
					in_type = lattice_to_dcp[lat_type]['class']
					in_id = obj['uuid']
					ins.append({'input_type': in_type, 'input_id': in_id})
					seq_method = 'high throughput sequencing'

			seq_prot = {
				'instrument_manufacturer_model': {
					'text': ','.join(sr_obj.get('platform',[]))
				},
				'method': {
					'text': seq_method
				},
				'paired_end': False
			}
			seq_prot_uuid = uuid_make(seq_prot)
			seq_prot['protocol_core'] = {'protocol_id': seq_prot_uuid}
			seq_prot['provenance'] = {'document_id': seq_prot_uuid}
			if not whole_dict.get('sequencing_protocol'):
				whole_dict['sequencing_protocol'] = [seq_prot]
			else:
				seq_prots = [i['provenance']['document_id'] for i in whole_dict['sequencing_protocol']]
				if seq_prot_uuid not in seq_prots:
					whole_dict['sequencing_protocol'].append(seq_prot)

			seq_prot_short = {'protocol_type': 'sequencing_protocol', 'protocol_id': seq_prot_uuid}
			protocols.append(seq_prot_short)

			# handle the library metadata
			if not whole_dict.get('library_preparation_protocol'):
				get_object(lib_obj)
			else:
				lib_prots = [i['provenance']['document_id'] for i in whole_dict['library_preparation_protocol']]
				if lib_obj['protocol']['uuid'] not in lib_prots:
					get_object(lib_obj)
			lib_type = lib_obj['@type'][0]
			dcp_type = lattice_to_dcp[lib_type]['class']
			lib_prot = {'protocol_type': dcp_type, 'protocol_id': lib_obj['protocol']['uuid']}
			protocols.append(lib_prot)

			all_susps.update(susps)
		l = {
			'protocols': protocols,
			'inputs': ins,
			'outputs': links_dict[seqruns]['outputs']
			}
		link_hash = uuid_make(l)
		l['link_type'] = 'process_link'
		l['process_type'] = 'process'
		l['process_id'] = link_hash
		link = {'links': [l]}
		links.append(link)

		add_process(link_hash)

	return all_susps, links


def handle_doc(doc_id):
	field_lst = ['attachment'] + \
		[v['lattice'] for v in lattice_to_dcp['Document'].values() if isinstance(v, dict)]
	doc_obj = lattice.get_report('Document', f'&@id={doc_id}', field_lst, connection)[0]

	link_info = {'file_type': 'supplementary_file', 'file_id': doc_obj['uuid']}
	doc_files.append(link_info)

	download_url = urljoin(server, doc_id + doc_obj['attachment']['href'])
	r = requests.get(download_url, auth=connection.auth, timeout=60)
	file_name = doc_obj['attachment']['download']
	open(file_name, 'wb').write(r.content)

	get_object(doc_obj)
	handled_docs.append(doc_id)


def create_protocol(in_type, out_type, out_obj):
	prots = []
	der_process = out_obj['derivation_process']
	my_obj = {
		'provenance': {},
		'protocol_core': {},
		'method': {
			'text': ','.join(der_process)
		}
	}

	if in_type == 'donor_organism':
		pr_type = 'collection_protocol'
		if out_obj.get('spatial_information'):
			my_obj['protocol_core']['protocol_description'] = out_obj['spatial_information']
	elif out_type == 'cell_suspension' and in_type != 'cell_suspension':
		pr_type = 'dissociation_protocol'
		desc = []
		if out_obj.get('dissociation_time'):
			desc.append('dissociation time:' + str(out_obj['dissociation_time']) + ' ' + out_obj['dissociation_time_units'])
		if out_obj.get('red_blood_cell_lysis') == True:
			desc.append('included red blood cell lysis')
		if out_obj.get('dissociation_reagent'):
			desc.append('used {}'.format(out_obj['dissociation_reagent']))
		if desc:
			my_obj['protocol_core']['protocol_description'] = ', '.join(desc)
		if out_obj.get('enrichment_factors') or out_obj.get('enriched_cell_types'):
			enrich_derivs = []
			my_derivs = []
			for d in der_process:
				if d in ['detergent solubilization','enzymatic dissociation','mechanical dissociation']:
					my_derivs.append(d)
				else:
					enrich_derivs.append(d)
			if my_derivs:
				my_obj['method']['text'] = ','.join(my_derivs)
			else:
				my_obj['method']['text'] = '{} dissociation'.format(out_obj['suspension_type'])
	elif out_type == 'organoid' and in_type == 'cell_line':
		pr_type = 'differentiation_protocol'
		my_obj['method'] = my_obj['method']['text']
	else:
		pr_type = 'protocol'
		del my_obj['method']

	for u in out_obj.get('urls', []):
		if 'protocols.io' in u and 'doi.org' in u:
			protocols_doi = u.split('doi.org/')[1]
			my_obj['protocol_core']['protocols_io_doi'] = protocols_doi

	pr_id = uuid_make([pr_type + in_type + out_type + str(my_obj)])
	prots.append({'protocol_type': pr_type, 'protocol_id': pr_id})
	my_obj['protocol_core']['protocol_id'] = pr_id
	my_obj['provenance']['document_id'] = pr_id
	if whole_dict.get(pr_type):
		whole_dict[pr_type].append(my_obj)
	else:
		whole_dict[pr_type] = [my_obj]

	if out_obj.get('treatment_summary'):
		pr_type = 'protocol'
		tr_obj = {
			'provenance': {},
			'protocol_core': {},
			'type': {
				'text': 'treated with ' + out_obj['treatment_summary']
				}
			}
		pr_id = uuid_make([pr_type + in_type + out_type + str(tr_obj)])
		prots.append({'protocol_type': pr_type, 'protocol_id': pr_id})
		tr_obj['protocol_core']['protocol_id'] = pr_id
		tr_obj['provenance']['document_id'] = pr_id
		if whole_dict.get(pr_type):
			whole_dict[pr_type].append(tr_obj)
		else:
			whole_dict[pr_type] = [tr_obj]

	if out_obj.get('enrichment_factors') or out_obj.get('enriched_cell_types'):
		pr_type = 'enrichment_protocol'
		enr_obj = {
			'provenance': {},
			'protocol_core': {},
			'method': {}
			}
		meth_txt = ['enrichment']
		if out_obj.get('enriched_cell_types'):
			enr_cells = [ct['term_name'] for ct in out_obj['enriched_cell_types']]
			meth_txt.append('for {}'.format(','.join(enr_cells)))
		try:
			enrich_derivs
		except NameError:
			enrich_derivs = []
		if enrich_derivs:
			meth_txt.append('by {}'.format(','.join(enrich_derivs)))
		enr_obj['method']['text'] = ' '.join(meth_txt)
		if out_obj.get('enrichment_factors'):
			enr_obj['markers'] = ','.join(out_obj['enrichment_factors'])
		enpr_id = uuid_make([pr_type + in_type + out_type + str(enr_obj)])
		prots.append({'protocol_type': pr_type,'protocol_id': enpr_id})	
		enr_obj['protocol_core']['protocol_id'] = enpr_id
		enr_obj['provenance']['document_id'] = enpr_id
		if whole_dict.get(pr_type):
			whole_dict[pr_type].append(enr_obj)
		else:
			whole_dict[pr_type] = [enr_obj]

	return prots


def add_links(temp_obj, der_fr, links):
	ins = []
	for i in der_fr:
		if isinstance(i, tuple): #donor ID + variable_age
			i = i[0]
		obj_type, filter_url = lattice.parse_ids([i])
		field_lst = ['uuid']
		obj = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
		in_type = lattice_to_dcp[obj_type]['class']
		in_id = obj['uuid']
		ins.append({'input_type': in_type, 'input_id': in_id})
	lat_type = temp_obj['@type'][0]
	out_type = lattice_to_dcp[lat_type]['class']
	out_id = temp_obj['uuid']
	outs = [{'output_type': out_type, 'output_id': out_id}]
	prots = create_protocol(in_type, out_type, temp_obj)
	link = {
		'outputs': outs,
		'inputs': ins,
		'protocols': prots
		}
	link_hash = uuid_make(link)
	link['link_type'] = 'process_link'
	link['process_type'] = 'process'
	link['process_id'] = link_hash

	add_process(link_hash)

	for i in links.copy():
		index = links.index(i)
		for i2 in i['links']:
			for i3 in i2['inputs']:
				if temp_obj['uuid'] == i3['input_id']:
					links[index]['links'].append(link)


def add_process(process_id):
	process = {
		'process_core': {'process_id': process_id},
		'provenance': {'document_id': process_id}
		}

	if whole_dict.get('process'):
		whole_dict['process'].append(process)
	else:
		whole_dict['process'] = [process]


def number_conversion(value):
	range_indicators = ['-','<','>']
	if True not in [c in value for c in range_indicators]:
		if '.' in value:
			return float(value)
		else:
			return int(value)
	else:
		return None


def unit_conversion(value, current_u, desired_u):
	ureg = UnitRegistry()
	curr = value * ureg(current_u)
	curr.ito(desired_u)
	return curr.magnitude


def customize_fields(obj, obj_type, dataset_id, docid_updates):
	if obj.get('provenance'):
		obj['provenance']['submission_date'] = dt
	else:
		obj['provenance'] = {'submission_date': dt}

	if obj_type in no_reuse_types:
		old_docid = obj['provenance']['document_id']
		new_docid = uuid_make(old_docid + dataset_id)
		obj['provenance']['document_id'] = uuid_make(old_docid + dataset_id)
		docid_updates[old_docid] = new_docid

	if obj.get('biomaterial_core'):
		if obj['biomaterial_core'].get('ncbi_taxon_id'):
			if obj_type == 'donor_organism':
				obj['biomaterial_core']['ncbi_taxon_id'] = [int(obj['biomaterial_core']['ncbi_taxon_id'])]
			else:
				temp = set()
				for d in obj['biomaterial_core']['ncbi_taxon_id']:
					path = d.get('organism.taxon_id').split(':')
					temp.add(int(path[1]))
				obj['biomaterial_core']['ncbi_taxon_id'] = list(temp)

	if obj.get('treatment_summary'):
		del obj['treatment_summary']

	if obj_type == 'project':
		if obj.get('project_core'):
			if obj['project_core'].get('project_title'):
				short = obj['project_core']['project_title'].replace(' ','')[:50]
				obj['project_core']['project_short_name'] = short
		if obj.get('funders'):
			if obj['funders'].get('organization'):
				funds = []
				for o in obj['funders']['organization']:
					funds.append({'organization': o, 'grant_id': 'unspecified'})
				obj['funders'] = funds
		if obj.get('corresponding_contributors'):
			if not obj.get('contributors'):
				obj['contributors'] = []
			for cc in obj['corresponding_contributors']:
				cc['corresponding_contributor'] = True
				obj['contributors'].append(cc)
			del obj['corresponding_contributors']
		if obj.get('publications'):
			for p in obj['publications']:
				if p.get('pmid'):
					p['pmid'] = int(p['pmid'])
				if p.get('authors'):
					p['authors'] = [a.lstrip() for a in p['authors'].split(',')]
				if p.get('doi'):
					p['url'] = 'https://doi.org/' + p['doi']
				if 'official_hca_publication' not in p:
					p['official_hca_publication'] = False
		if obj.get('supplementary_links'):
			for l in obj['supplementary_links']:
				if l == f'https://explore.data.humancellatlas.org/projects/{dataset_id}':
					obj['supplementary_links'].remove(l)

	elif obj_type == 'donor_organism':
		if obj.get('genus_species'):
			obj['genus_species'] = [obj['genus_species']]
		if obj.get('organism_age') in ['unknown','variable']:
			del obj['organism_age']
		elif '>' in obj.get('organism_age',''):
			obj['organism_age'] = obj['organism_age'].replace('>','') + '-122'
			del obj['organism_age_unit']
		if not obj.get('is_living'):
			if obj.get('gestational_age'):
				obj['is_living'] = 'not applicable'
			else:
				obj['is_living'] = 'unknown'
		if obj.get('gestational_age'):
			if obj['gestational_age'] == 'unknown':
				del obj['gestational_age']
			elif obj['gestational_age_unit']['text'] == 'month' or '-' in obj['gestational_age']:
				del obj['gestational_age']
				del obj['gestational_age_unit']
			else:
				age = int(obj['gestational_age'])
				if obj['gestational_age_unit']['text'] == 'week':
					obj['gestational_age'] = str(age + 2)
				elif obj['gestational_age_unit']['text'] == 'day':
					obj['gestational_age'] = str(age + 14)
		if obj.get('human_specific'):
			if obj['human_specific'].get('ethnicity'):
				new_eths = []
				for e in obj['human_specific']['ethnicity']:
					new_eths.append({
						'ontology': e['term_id'],
						'ontology_label': e['term_name'],
						'text': e['term_name']
						})
				obj['human_specific']['ethnicity'] = new_eths
			if obj['human_specific'].get('body_mass_index'):
				if '-' in obj['human_specific']['body_mass_index'] or obj['human_specific']['body_mass_index'] == 'variable':
					del obj['human_specific']['body_mass_index']
				else:
					obj['human_specific']['body_mass_index'] = float(obj['human_specific']['body_mass_index'])
		if obj.get('mouse_specific'):
			if obj['mouse_specific'].get('strain'):
				obj['mouse_specific']['strain'] = [obj['mouse_specific']['strain']]
		if obj.get('medical_history'):
			if obj['medical_history'].get('test_results'):
				tr = [k + ':' + v for k,v in obj['medical_history']['test_results'].items()]
				obj['medical_history']['test_results'] = ','.join(tr)
		if obj.get('death'):
			causes_of_death = [f"{cd['term_name']} [{cd['term_id']}]" for cd in obj['death']['cause_of_death']]
			obj['death'] = {'cause_of_death': ','.join(causes_of_death)}

	elif obj_type == 'specimen_from_organism':
		if obj.get('spatial_information'):
			del obj['spatial_information']
		if obj.get('purchased_specimen'):
			if not obj['purchased_specimen'].get('catalog_number'):
				del obj['purchased_specimen']
		if obj.get('organ_parts'):
			obj['organ_parts'] = [obj['organ_parts']]
		if obj.get('organ'):
			if obj['organ'].get('text'):
				obj['organ']['text'] = obj['organ']['text'][0]
		else:
			obj['organ'] = {'text': obj['organ_parts'][0]['ontology_label']}
		if obj.get('preservation_storage'):
			if obj['preservation_storage'].get('storage_time'):
				new_v = number_conversion(obj['preservation_storage']['storage_time'])
				if new_v:
					obj['preservation_storage']['storage_time'] = new_v
				else:
					del obj['preservation_storage']['storage_time']
					del obj['preservation_storage']['storage_time_unit']
		if obj.get('state_of_specimen'):
			if obj['state_of_specimen'].get('ischemic_time'):
				new_v = number_conversion(obj['state_of_specimen']['ischemic_time'])
				if new_v:
					obj['state_of_specimen']['ischemic_time'] = new_v
				else:
					del obj['state_of_specimen']['ischemic_time']
					del obj['state_of_specimen']['ischemic_time_units']
			if obj['state_of_specimen'].get('ischemic_time_units'):
				if obj['state_of_specimen']['ischemic_time_units'] != 'second':
					new_v = unit_conversion(obj['state_of_specimen']['ischemic_time'],obj['state_of_specimen']['ischemic_time_units'], 'second')
					obj['state_of_specimen']['ischemic_time'] = new_v
					del obj['state_of_specimen']['ischemic_time_units']
			if obj['state_of_specimen'].get('postmortem_interval'):
				new_v = number_conversion(obj['state_of_specimen']['postmortem_interval'])
				if new_v:
					obj['state_of_specimen']['postmortem_interval'] = new_v
				else:
					del obj['state_of_specimen']['postmortem_interval']
					del obj['state_of_specimen']['postmortem_interval_units']
			if obj['state_of_specimen'].get('postmortem_interval_units'):
				if obj['state_of_specimen']['postmortem_interval_units'] != 'second':
					new_v = unit_conversion(obj['state_of_specimen']['postmortem_interval'],obj['state_of_specimen']['postmortem_interval_units'], 'second')
					obj['state_of_specimen']['postmortem_interval'] = new_v
					del obj['state_of_specimen']['postmortem_interval_units']

	elif obj_type == 'imaged_specimen':
		if obj.get('slice_thickness'):
			obj['slice_thickness'] = float(obj['slice_thickness'])

	elif obj_type == 'cell_line':
		if obj.get('cell_line_type'):
			celltypes = obj['cell_line_type']
			if 'induced pluripotent stem cell' in celltypes:
				obj['cell_line_type'] = 'induced pluripotent'
			elif 'stem cell derived cell line' in celltypes:
				obj['cell_line_type'] = 'stem cell-derived'
			elif 'stem cell' in celltypes:
				obj['cell_line_type'] = 'stem cell'
			else:
				obj['cell_line_type'] = 'primary'
		else:
			obj['cell_line_type'] = 'primary'
		if obj['biomaterial_core'].get('insdc_sample_accession'):
			obj['biomaterial_core']['insdc_biomaterial'] = obj['biomaterial_core']['insdc_sample_accession']
			del obj['biomaterial_core']['insdc_sample_accession']
		if obj['biomaterial_core'].get('biosamples_accession'):
			obj['biomaterial_core']['biosd_biomaterial'] = obj['biomaterial_core']['biosamples_accession']
			del obj['biomaterial_core']['biosamples_accession']
		if obj.get('tissue'):
			if obj['tissue'].get('text'):
				obj['tissue']['text'] = obj['tissue']['text'][0]
		if obj.get('supplier') and not obj.get('catalog_number'):
			del obj['supplier']
	elif obj_type == 'organoid':
		if obj.get('age_unit') and not obj.get('age'):
			obj['age'] = 0
		if obj.get('model_organ'):
			if obj['model_organ'].get('text'):
				obj['model_organ']['text'] = obj['model_organ']['text'][0]
		else:
			obj['model_organ'] = {'text': obj['model_organ_part']['ontology_label']}

	elif obj_type == 'cell_suspension':
		if obj.get('dissociation_time'):
			del obj['dissociation_time']
			del obj['dissociation_time_units']
		if obj.get('dissociation_reagent'):
			del obj['dissociation_reagent']
		if obj.get('red_blood_cell_lysis'):
			del obj['red_blood_cell_lysis']
		if obj.get('enrichment_factors'):
			del obj['enrichment_factors']
		if obj.get('suspension_type'):
			del obj['suspension_type']
		if obj.get('cell_morphology'):
			if obj['cell_morphology'].get('cell_size'):
				obj['cell_morphology']['cell_size'] = str(obj['cell_morphology']['cell_size'])
			if obj['cell_morphology'].get('percent_cell_viability'):
				new_v = number_conversion(obj['cell_morphology']['percent_cell_viability'])
				if new_v:
					obj['cell_morphology']['percent_cell_viability'] = new_v
				else:
					del obj['cell_morphology']['percent_cell_viability']
		if obj.get('estimated_count_units'):
			if obj['estimated_count_units'] not in ['cells', 'nuclei']:
				del obj['estimated_cell_count']
			del obj['estimated_count_units']
		if obj.get('estimated_cell_count'):
			new_v = number_conversion(obj['estimated_cell_count'])
			if new_v:
				obj['estimated_cell_count'] = new_v
			else:
				del obj['estimated_cell_count']

	elif obj_type == 'library_preparation_protocol':
		if not obj.get('strand'):
			obj['strand'] = 'not provided'
		if not obj.get('end_bias'):
			obj['end_bias'] = 'full length'
		if obj.get('polyA_selection'):
			obj['input_nucleic_acid_molecule'] = {'text': 'polyA ' +obj['input_nucleic_acid_molecule']['text']}
			del obj['polyA_selection']
		if obj.get('input_nucleic_acid_molecule') and obj['input_nucleic_acid_molecule']['text'] in input_onts:
			obj['input_nucleic_acid_molecule']['ontology'] = input_onts[obj['input_nucleic_acid_molecule']['text']]
		if obj.get('read_structure'):
			for rs in obj['read_structure']:
				length = rs['end'] - rs['start'] + 1
				read_type = rs['located_in_read_type'].strip('N').replace('index','Index')
				if rs['sequence_element'] == 'cell barcode':
					obj['cell_barcode'] = {
						'barcode_read': read_type,
						'barcode_offset': rs['start'] - 1, #1-based to 0-based
						'barcode_length': length
					}
					if obj.get('cell_barcode_whitelist'):
						obj['cell_barcode']['white_list_file'] = obj['cell_barcode_whitelist']
						del obj['cell_barcode_whitelist']
				if rs['sequence_element'] == 'spatial barcode':
					obj['spatial_barcode'] = {
						'barcode_read': read_type,
						'barcode_offset': rs['start'] - 1, #1-based to 0-based
						'barcode_length': length
					}
				if rs['sequence_element'] == 'UMI':
					obj['umi_barcode'] = {
						'barcode_read': read_type,
						'barcode_offset': rs['start'] - 1, #1-based to 0-based
						'barcode_length': length
					}
			del obj['read_structure']

	elif obj_type == 'sequence_file':
		if obj.get('insdc_run_accessions'):
			v = []
			for i in obj['insdc_run_accessions']:
				if i.get('dbxrefs'):
					for ref in i['dbxrefs']:
						no_pre = ref.split(':')[1]
						v.append(no_pre)
		if v:
			obj['insdc_run_accessions'] = v
		else:
			del obj['insdc_run_accessions']

	elif obj_type == 'supplementary_file':
		file_format = obj['file_core']['file_name'].split('.')[-1]
		obj['file_core']['format'] = file_format


def file_descript(obj, obj_type, dataset):
	if obj.get('content_type'):
		content_type = obj['content_type']
		del obj['content_type']
	elif obj_type == 'sequence_file':
		content_type = 'application/gzip'

	file_descriptor = {
		'describedBy': 'https://schema.humancellatlas.org/system/{}/file_descriptor'.format(dcp_versions['file_descriptor']),
		'schema_type': 'file_descriptor',
		'file_name': obj['file_core']['file_name'],
		'file_id': obj['provenance']['document_id'],
		'file_version': dt,
		'content_type': content_type,
		'size': obj['file_size'],
		'sha256': obj['sha256'],
		'crc32c': obj['crc32c']
	}
	del obj['file_size']
	del obj['sha256']
	del obj['crc32c']
	with open(f'{output_dir}/{dataset}/descriptors/{obj_type}/{file_descriptor["file_id"]}_{file_descriptor["file_version"]}.json', 'w') as outfile:
		json.dump(file_descriptor, outfile, indent=4)


def main():
	logging.info('GETTING THE DATASET')
	field_lst = ['status','original_files','audit'] + \
		[v['lattice'] for v in lattice_to_dcp['Dataset'].values() if isinstance(v, dict)]
	ds_obj = lattice.get_report('Dataset', f'&accession={args.dataset}', field_lst, connection)[0]

	# check status of the dataset
	if ds_obj.get('status') not in ['in progress', 'released']:
		print('WARNING: Dataset status is {}'.format(ds_obj.get('status')))
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to one or more ERROR audits')

	# in case a project is already in the DCP, we may need to match the project ID to that
	lattice_dataset_id = ds_obj['uuid']
	dcp_projects = {
		'2a13992b-c518-4f12-8536-09c92d51d707': '16e99159-78bc-44aa-b479-55a5e903bf50', # van Zyl et al 2022 (Sanes)
		'49104d6a-0e30-4180-94b9-bf552b110686': 'dbd836cf-bfc2-41f0-9834-41cc6c0b235a', # Lavaert et al 2020
		'e62da5d2-33fb-4a3d-bdd7-4daca3a042cd': '9c20a245-f2c0-43ae-82c9-2232ec6b594f', # Chen retina
		'320f295d-9667-4d64-a58c-e4e9bd051546': 'abe1a013-af7a-45ed-8c26-f3793c24a1f4', # Stewart et al 2019
		'35c1d480-ceac-452e-b08a-f710029f73d9': '577c946d-6de5-4b55-a854-cd3fde40bff2', # Wilson et al 2019
		'5efd8bd0-2125-4c61-bc8c-f03dbd8370da': '4bec484d-ca7a-47b4-8d48-8830e06ad6db', # Voigt et al 2019
		'a4aee8b9-8d32-4e14-91e8-42961671edc1': 'c1810dbc-16d2-45c3-b45e-3e675f88d87b', # Park et al 2020
		'a52f6a97-a674-4f46-a32d-77ce7b03092c': '8185730f-4113-40d3-9cc3-929271784c2b', # Lukowski et al 2019
		'bdbe72d3-76f7-4dd7-bb5d-80cd02382399': 'ad98d3cd-26fb-4ee3-99c9-8a2ab085e737', # Litviňuková et al 2020
		'c9aa45d0-cacd-4c0e-999c-f10df0609a31': '425c2759-db66-4c93-a358-a562c069b1f1'  # Guilliams et al 2022
		}

	dataset_id = dcp_projects.get(lattice_dataset_id, lattice_dataset_id)
	ds_obj['uuid'] = dataset_id

	if args.validate_only:
		dcp_errors = dcp_validation(dataset_id, output_dir)
		if dcp_errors != 0:
			logging.info('WARNING: {} files with DCP schema errors'.format(str(dcp_errors)))
		sys.exit('Exiting after DCP schema validation, as per validate-only option')

	logging.info('WILL WRITE FILES TO THE {} DIRECTORY'.format(dataset_id))

	# convert the dataset to DCP schema
	get_object(ds_obj)

	# get the validated raw sequence files from that dataset
	logging.info('GETTING RAW SEQUENCE FILES')
	links_dict = {}

	files = [f for f in ds_obj['original_files'] if 'raw-sequence-files' in f]
	field_lst = ['validated'] + \
		[v['lattice'] for v in lattice_to_dcp['RawSequenceFile'].values() if isinstance(v, dict)]
	file_objs = lattice.get_report(
		'RawSequenceFile',
		f'&dataset=/datasets/{args.dataset}/&status!=archived&status!=deleted',
		field_lst, connection
		)

	for temp_obj in file_objs:
		if temp_obj.get('validated') == False:
			logging.info('{} has not been validated, will be excluded'.format(temp_obj['@id']))
			not_valid.append(temp_obj['@id'])
		else:
			# convert each to DCP schema
			get_object(temp_obj)
			# pull the derived_from to store for later formation to links
			der_from = [i['@id'] for i in temp_obj['derived_from']]
			get_links(temp_obj, tuple(der_from), links_dict)

	if not links_dict:
		sys.exit('No validated RawSequenceFiles associated with this dataset')

	# special walkback of graph until Suspension object
	# gather all the Suspension objects to traverse next
	# set up links between sequence_file and suspension as the start of each subgraph
	logging.info('GETTING THE GRAPH FROM RAW SEQUENCE FILES BACK TO SUSPENSIONS')
	susps, links = seq_to_susp(links_dict)

	# walkback graph the rest of the way
	logging.info('GETTING THE GRAPH FROM SUSPENSIONS BACK TO DONORS')
	seen = set()
	remaining = susps
	while remaining:
		seen.update(remaining)
		next_remaining = set()
		for identifier in remaining:
			if isinstance(identifier, tuple):
				i = identifier[0]
				obj_type, filter_url = lattice.parse_ids([i])
				field_lst = [v['lattice'] for v in lattice_to_dcp[obj_type].values() if isinstance(v, dict)]
				field_lst.extend(['derived_from','derivation_process'])
				temp_obj = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
				get_object(temp_obj, identifier[1])
			else:
				obj_type, filter_url = lattice.parse_ids([identifier])
				field_lst = [v['lattice'] for v in lattice_to_dcp[obj_type].values() if isinstance(v, dict)]
				field_lst.extend(['derived_from','derivation_process'])
				temp_obj = lattice.get_report(obj_type, filter_url, field_lst, connection)[0]
				get_object(temp_obj)
			get_derived_from(temp_obj, next_remaining, links)
		remaining = next_remaining - seen

	dir_to_make = ['', 'metadata', 'links', 'data', 'descriptors', 'descriptors/sequence_file']

	if doc_files:
		doc_link = {
			'link_type': 'supplementary_file_link',
			'entity': {'entity_type': 'project', 'entity_id': dataset_id},
			'files': doc_files
		}
		for i in links:
			i['links'].append(doc_link)
		dir_to_make.append('descriptors/supplementary_file')

	# make directory named after dataset
	logging.info('WRITING THE JSON FILES')
	for d in dir_to_make:
		os.mkdir(f'{output_dir}/{dataset_id}/{d}')

	s3_uris = []
	ftp_uris = []

	# write a json file for each object
	docid_updates = {}
	for k in whole_dict.keys():
		os.mkdir(f'{output_dir}/{dataset_id}/metadata/{k}')
		for o in whole_dict[k]:
			if k == 'sequence_file':
				if o.get('s3_uri'):
					s3_uris.append(o['s3_uri'])
					o['file_core']['file_name'] = o['s3_uri'].split('/')[-1]
					del o['s3_uri']
				elif o.get('external_uri'):
					ftp_uris.append(o['external_uri'])
					o['file_core']['file_name'] = o['external_uri'].split('/')[-1]
					del o['external_uri']
				file_descript(o, k, dataset_id)
			elif k == 'supplementary_file':
				file_name = o['file_core']['file_name']
				file_stats(file_name, o)
				os.rename(file_name, f'{output_dir}/{dataset_id}/data/{file_name}')
				file_descript(o, k, dataset_id)
			customize_fields(o, k, dataset_id, docid_updates)
			o['schema_type'] = dcp_types[k].split('/')[0]
			o['schema_version'] = dcp_versions[k]
			o['describedBy'] = 'https://schema.humancellatlas.org/type/{}/{}/{}'.format(dcp_types[k], dcp_versions[k], k)
			with open(f'{output_dir}/{dataset_id}/metadata/{k}/{o["provenance"]["document_id"]}_{dt}.json', 'w', encoding='utf8') as outfile:
				json.dump(o, outfile, indent=4, ensure_ascii=False)

	# reformat links to use uuids and convert to DCP schema
	for i in links:
		i['schema_type'] = 'links'
		i['schema_version'] = dcp_versions['links']
		i['describedBy'] = 'https://schema.humancellatlas.org/system/{}/links'.format(dcp_versions['links'])
		first_id = i['links'][0]['process_id']

		for l in i['links']:
			for inp in l['inputs']:
				inp['input_id'] = docid_updates.get(inp['input_id'],inp['input_id'])
			for outp in l['outputs']:
				outp['output_id'] = docid_updates.get(outp['output_id'],outp['output_id'])
			for prot in l['protocols']:
				prot['protocol_id'] = docid_updates.get(prot['protocol_id'],prot['protocol_id'])

		with open(f'{output_dir}/{dataset_id}/links/{first_id}_{dt}_{dataset_id}.json', 'w') as outfile:
			json.dump(i, outfile, indent=4)
			outfile.close()

	# make a staging_area.json to meet DCP requirements
	text = {
		"is_delta": args.revision
	}
	with open(f'{output_dir}/{dataset_id}/staging_area.json', 'w', encoding='utf8') as outfile:
		json.dump(text, outfile, indent=4, ensure_ascii=False)

	# logging.info a close approximation of the DCP metadata.tsv
	logging.info('PREPARING MOCK DCP METADATA TSV')
	tsv_report(dataset_id, output_dir)

	# report metadata not mapped to DCP schema
	for k,v in not_incl.items():
		not_incl[k] = list(v)
	with open(f'{output_dir}/not_included.json', 'w') as outfile:
		json.dump(not_incl, outfile, indent=4)

	# report files that are not validated, and thus, not included in the mapping
	if not_valid:
		with open(f'{output_dir}/not_validated.txt', 'w') as outfile:
			outfile.write('\n'.join(not_valid))

	if not args.no_validate:
		logging.info('VALIDATING AGAINST DCP SCHEMA')
		dcp_errors = dcp_validation(dataset_id, output_dir)
		if dcp_errors != 0:
			print('WARNING: {} files with DCP schema errors'.format(str(dcp_errors)))
			if args.update:
				i = input('Continue? y/n: ')
				if i.lower() not in ['y','yes']:
					sys.exit('Stopped due to one or more DCP schema errors')

	if args.update:
		if not args.data_only:
			# transfer the metadata directory to the DCP Google Cloud project
			logging.info('TRANSFERRING METADATA DIRECTORIES')
			request_to_gcp.local_dir_transfer(f'{output_dir}/{dataset_id}', output_dir)

			if args.metadata_only:
				sys.exit('Metadata directories transferred, data not transferred due to --metadata-only')

		# transfer the data files from S3 to the DCP Google Cloud project
		if s3_uris:
			logging.info('TRANSFERRING {} FILES FROM S3'.format(str(len(s3_uris))))
			request_to_gcp.aws_file_transfer(dataset_id, s3_uris)

		# transfer the data files from external FTPs to the DCP Google Cloud project
		if ftp_uris:
			logging.info('TRANSFERRING {} EXTERNAL FILES'.format(str(len(ftp_uris))))
			request_to_gcp.ftp_file_transfer(dataset_id, ftp_uris)

	else:
		logging.info('Metadata directories written locally, but not transferring without the --update option')
		quit()

if __name__ == '__main__':
	output_dir = lattice.create_subdirectory('DCPmap')

	logging.basicConfig(filename=f'{output_dir}/mapper.log', level=logging.INFO)
	logging.info('STARTED')
	# set the current date time, used to version throughout
	d_now = datetime.now(tz=timezone.utc).isoformat(timespec='auto')
	dt = str(d_now).replace('+00:00', 'Z')

	# used to track metdata not mapped to DCP schema
	if os.path.exists('not_included.json'):
		not_incl = json.load(open('not_included.json'))
	else:
		not_incl = {}

	whole_dict = {}
	not_valid = []
	doc_files = []
	handled_docs = []

	input_onts = {
		'DNA': 'CHEBI:16991',
		'RNA': 'CHEBI:33697',
		'protein': 'CHEBI:36080'
	}

	no_reuse_types = [
		'cell_suspension',
		'collection_protocol',
		'dissociation_protocol',
		'donor_organism',
		'enrichment_protocol',
		'imaged_specimen',
		'library_preparation_protocol',
		'protocol',
		'sequencing_protocol',
		'specimen_from_organism'
	]

	args = getArgs()
	if not args.dataset:
		sys.exit('ERROR: --dataset is required')
	if not args.mode:
		sys.exit('ERROR: --mode is required')

	connection = lattice.Connection(args.mode)
	server = connection.server

	# check for Lattice schema errors, possibly outdated fields
	schema_errors = test_mapping()
	if schema_errors > 0:
		print(str(schema_errors) + ' schema errors found')
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to schema errors')

	main()
