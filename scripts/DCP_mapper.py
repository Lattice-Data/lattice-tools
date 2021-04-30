import argparse
import hashlib
import json
import lattice
import os
import request_to_gcp
import requests
import sys
import zlib
from datetime import datetime, timezone
from pint import UnitRegistry
from urllib.parse import urljoin
from property_mapping import (
    lattice_to_dcp
)


EPILOG = '''
Get a Dataset from the Lattice database and put into DCP schema and exported in json files.

Examples:

    python %(prog)s --mode prod --dcp /Users/jason/GitClones/HCA/ --dataset LATDS940AQG

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
    parser.add_argument('--dcp',
                        help='The directory to find the DCP metadata-schema.')
    args = parser.parse_args()
    return args


def test_mapping():
	error_count = 0
	schema_url = urljoin(server, 'profiles/?format=json')
	schemas = requests.get(schema_url).json()
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
						print(k)
						print(path)
					if len(path) > 2:
						print('-----long path----')
				elif lat_prop in schema_props:
					flag = True
			if flag != True:
				error_count += 1
				print('{} not found in schema for {}'.format(lat_prop,k))
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
			elif my_obj.get('prop'):
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
		'donors','ethnicity','files','href','i5_index_file','i7_index_file','internal_contact','lab','no_file_available','notes','organism',
		'output_types','protocol','read_1_file','read_1N_file','read_2_file','read_2N_file','read_length_units','reference_annotation',
		'reference_assembly','references','revoked_files','schema_version','sequence_elements','software','title','url','uuid','validated']

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
        hash = 0
        while True:
            s = f.read(65536)
            if not s:
                break
            hash = zlib.crc32(s, hash)
        obj['crc32c'] = "%08X" % (hash & 0xFFFFFFFF)


def get_object(temp_obj):
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

	# drop unneeded properties, flatten subobjects
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

	# add the object of mapped properties
	if whole_dict.get(dcp_obj_type):
		whole_dict[dcp_obj_type].append(my_obj)
	else:
		whole_dict[dcp_obj_type] = [my_obj]


def flatten_obj(obj):
	ignore = ['accession','actions','aliases','audit','contributing_files',
		'date_created','libraries','original_files','status','submitted_by',
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
		if isinstance(temp_obj['derived_from'], str):
			next_remaining.add(temp_obj['derived_from'])
			add_links(temp_obj, tuple(temp_obj['derived_from']), links)
		elif isinstance(temp_obj['derived_from'], dict):
			der_fr = temp_obj['derived_from']['@id']
			next_remaining.add(der_fr)
			add_links(temp_obj, tuple(der_fr), links)
		elif isinstance(temp_obj['derived_from'], list):
			if isinstance(temp_obj['derived_from'][0], str):
				next_remaining.update(temp_obj['derived_from'])
				add_links(temp_obj, tuple(temp_obj['derived_from']), links)
			elif isinstance(temp_obj['derived_from'][0], dict):
				der_fr = ()
				for o in temp_obj['derived_from']:
					next_remaining.add(o['@id'])
					der_fr = der_fr + (o['@id'],)
				add_links(temp_obj, der_fr, links)


def get_links(temp_obj, der_fr, links_dict):
	lat_type = temp_obj['@type'][0]
	out_type = lattice_to_dcp[lat_type]['class']
	outs = {'output_type': out_type, 'output_id': temp_obj['uuid']}
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
	all_susps = set()
	for seqruns in list(links_dict.keys()):
		l = {}
		ins = []
		susps = []
		protocols = []
		for sr in seqruns:
			url = urljoin(server, sr + '/?format=json')
			sr_obj = requests.get(url, auth=connection.auth).json()
			get_object(sr_obj)
			lib = sr_obj['derived_from'][0]
			lib_url = urljoin(server, lib + '/?format=json')
			lib_obj = requests.get(lib_url, auth=connection.auth).json()
			if not whole_dict.get('library_preparation_protocol'):
				get_object(lib_obj)
			else:
				lib_prots = [i['provenance']['document_id'] for i in whole_dict['library_preparation_protocol']]
				if lib_obj['protocol']['uuid'] not in lib_prots:
					get_object(lib_obj)
			susps.extend([i['uuid'] for i in lib_obj['derived_from']])
			all_susps.update(susps)
			for obj in lib_obj['derived_from']:
				lat_type = obj['@type'][0]
				in_type = lattice_to_dcp[lat_type]['class']
				ins.append({'input_type': in_type, 'input_id': obj['uuid']})
			lib_type = lib_obj['@type'][0]
			dcp_type = lattice_to_dcp[lib_type]['class']
			lib_prot = {'protocol_type': dcp_type, 'protocol_id': lib_obj['protocol']['uuid']}
			protocols.append(lib_prot)
			sr_type = sr_obj['@type'][0]
			dcp_type = lattice_to_dcp[sr_type]['class']
			seq_prot = {'protocol_type': dcp_type, 'protocol_id': sr_obj['uuid']}
			protocols.append(seq_prot)
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
	return all_susps, links


def handle_doc(doc_id):
	doc_url = urljoin(server, doc_id + '/?format=json')
	doc_obj = requests.get(doc_url, auth=connection.auth).json()

	link_info = {'file_type': 'supplementary_file', 'file_id': doc_obj['uuid']}
	doc_files.append(link_info)

	download_url = urljoin(server, doc_id + doc_obj['attachment']['href'])
	r = requests.get(download_url, auth=connection.auth)
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
	elif out_type == 'organoid' and in_type == 'cell_line':
		pr_type = 'differentiation_protocol'
		my_obj['method'] = my_obj['method']['text']
	else:
		pr_type = 'protocol'
		del my_obj['method']

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
			'method': {
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
			'method': {
				'text': 'enrichment'
				}
			}
		if out_obj.get('enriched_cell_types'):
			enr_cells = [ct['term_name'] for ct in out_obj['enriched_cell_types']]
			enr_obj['method']['text'] = 'enrichment for {}'.format(','.join(enr_cells))
		if out_obj.get('enrichment_factors'):
			enr_obj['markers'] = out_obj['enrichment_factors']
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
		url = urljoin(server, i + '/?format=json')
		obj = requests.get(url, auth=connection.auth).json()
		lat_type = obj['@type'][0]
		in_type = lattice_to_dcp[lat_type]['class']
		ins.append({'input_type': in_type, 'input_id': obj['uuid']})
	lat_type = temp_obj['@type'][0]
	out_type = lattice_to_dcp[lat_type]['class']
	outs = [{'output_type': out_type, 'output_id': temp_obj['uuid']}]
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
	for i in links.copy():
		index = links.index(i)
		for i2 in i['links']:
			for i3 in i2['inputs']:
				if temp_obj['uuid'] == i3['input_id']:
					links[index]['links'].append(link)


def get_dcp_schema_ver(directory):
	vers = {}
	if not directory.endswith('/'):
		directory = directory + '/'
	v_file = directory + 'json_schema/versions.json'
	versions = json.load(open(v_file))
	for v in versions['version_numbers'].values():
		for k2,v2 in v.items():
			if isinstance(v2, dict):
				for k3,v3 in v2.items():
					if isinstance(v3, dict):
						for k4,v4 in v3.items():
							vers[k4] = v4
					else:
						vers[k3] = v3
			else:
				vers[k2] = v2
	return vers


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


def customize_fields(obj, obj_type):
	if obj.get('provenance'):
		obj['provenance']['submission_date'] = dt
	else:
		obj['provenance'] = {'submission_date': dt}
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
					funds.append({'organization': o})
				obj['funders']['organization'] = 'Chan Zuckerberg Initiative'
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
		if obj.get('supplementary_links'):
			for l in obj['supplementary_links']:
				if l.startswith('https://data.humancellatlas.org'):
					obj['supplementary_links'].remove(l)
	elif obj_type == 'donor_organism':
		if not obj.get('is_living'):
			if obj['development_stage']['ontology_label'] in ['embryonic','fetal']:
				obj['is_living'] = 'not applicable'
			else:
				obj['is_living'] = 'unknown'
		if obj.get('human_specific'):
			if obj['human_specific'].get('ethnicity'):
				obj['human_specific']['ethnicity'] = [obj['human_specific']['ethnicity']]
		if obj.get('mouse_specific'):
			if obj['mouse_specific'].get('strain'):
				obj['mouse_specific']['strain'] = [obj['mouse_specific']['strain']]
		if obj.get('medical_history'):
			if obj['medical_history'].get('test_results'):
				tr = [k + ':' + v for k,v in obj['medical_history']['test_results'].items()]
				obj['medical_history']['test_results'] = ','.join(tr)
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
	elif obj_type == 'cell_line':
		if obj.get('type'):
			celltypes = obj['type']
			if 'induced pluripotent stem cell' in celltypes:
				obj['type'] = 'induced pluripotent'
			elif 'stem cell derived cell line' in celltypes:
				obj['type'] = 'stem cell-derived'
			elif 'stem cell' in celltypes:
				obj['type'] = 'stem cell'
			else:
				obj['type'] = 'primary'
		else:
			obj['type'] = 'primary'
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
	elif obj_type == 'sequencing_protocol':
		if 'paired_end' not in obj:
			obj['paired_end'] = False
		if not obj.get('method'):
			obj['method'] = {'text': 'high throughput sequencing'}


def remove_cell_lines(links, whole_dict):
	keep_cell_line = set()
	consolidated_links = []
	for l in links:
		consolidate = {}
		sublinks = l['links']
		for link in sublinks:
			for i in link['inputs']:
				if i['input_type'] == 'cell_line':
					ins = sorted([x['input_id'] for x in link['inputs']])
					if '#'.join(ins) not in consolidate:
						consolidate['#'.join(ins)] = {}
					consolidate['#'.join(ins)]['ins'] = link
			for o in link['outputs']:
				if o['output_type'] == 'cell_line':
					outs = sorted([x['output_id'] for x in link['outputs']])
					if '#'.join(outs) not in consolidate:
						consolidate['#'.join(outs)] = {}
					consolidate['#'.join(outs)]['outs'] = link
		if consolidate:
			for k,v in consolidate.items():
				# we don't want to collapse any cell_line-to-suspension links
				if v['ins']['outputs'][0]['output_type'] == 'cell_suspension':
					keep_cell_line.update(k.split('#'))
				else:
					new_link = {
						'outputs': v['ins']['outputs'],
						'inputs': v['outs']['inputs'],
						'protocols': v['ins']['protocols'] + v['outs']['protocols'],
						'link_type': 'process_link',
						'process_type': 'process',
						'process_id': v['outs']['process_id']
					}
					sublinks.remove(v['ins'])
					sublinks.remove(v['outs'])
					sublinks.append(new_link)
			consolidated_links.append({'links': sublinks})
		else:
			consolidated_links.append(l)

	# we need to remove cell_line objects if they don't feed directly to suspension
	if whole_dict.get('cell_line'):
		for i in whole_dict['cell_line']:
			if i['biomaterial_core']['biomaterial_id'] not in keep_cell_line:
				whole_dict['cell_line'].remove(i)

	return consolidated_links


def file_descript(obj, obj_type, dataset):
	if obj.get('content_type'):
		content_type = obj['content_type']
		del obj['content_type']
	elif obj_type == 'sequence_file':
		content_type = 'application/gzip'

	file_descriptor = {
		'describedBy': 'https://schema.humancellatlas.org/system/{}/file_descriptor'.format(dcp_vs['file_descriptor']),
		'schema_type': 'file_descriptor',
		'file_name': obj['file_core']['file_name'],
		'file_id': obj['provenance']['document_id'],
		'file_version': dt,
		'content_type': content_type,
		'size': obj['file_size'],
		'sha256': '', # obj['sha256'] UPDATE AFTER RELEASE
		'crc32c': '' # obj['crc32c'] UPDATE AFTER RELEASE
	}
	del obj['file_size']
	# del obj['sha256'] UPDATE AFTER RELEASE
	# del obj['crc32c'] UPDATE AFTER RELEASE
	with open(dataset + '/descriptors/' + obj_type + '/' + file_descriptor['file_id'] + '_' + file_descriptor['file_version'] + '.json', 'w') as outfile:
		json.dump(file_descriptor, outfile, indent=4)


def main():
	# get the dataset, and convert it to a project object
	print('GETTING THE DATASET')
	url = urljoin(server, args.dataset + '/?format=json')
	ds_obj = requests.get(url, auth=connection.auth).json()

	# check status of the dataset
	if ds_obj.get('status') not in ['in progress', 'released']:
		print('WARNING: Dataset status is {}'.format(ds_obj.get('status')))
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to one or more ERROR audits')

	dataset_id = ds_obj['uuid']
	get_object(ds_obj)

	# get the raw sequence files from that dataset
	print('GETTING RAW SEQUENCE FILES')
	links_dict = {}
	files = [i for i in ds_obj['files']]
	for f in files:
		url = urljoin(server, f + '/?format=json')
		temp_obj = requests.get(url, auth=connection.auth).json()
		obj_type = temp_obj['@type'][0]
		if obj_type == 'RawSequenceFile':
			if temp_obj.get('validated') == False:
				print('{} has not been validated, will be excluded'.format(temp_obj['@id']))
				not_valid.append(temp_obj['@id'])
			else:
				# convert each to a sequence_file
				get_object(temp_obj)
				# pull the derived_from to store for later formation to links
				der_from = [i['@id'] for i in temp_obj['derived_from']]
				get_links(temp_obj, tuple(der_from), links_dict)

	if not links_dict:
		sys.exit('No RawSequenceFiles associated with this dataset')

	# special walkback of graph until Suspension object
	# gather all the Suspension objects to traverse next
	# set up links between sequence_file and suspension as the start of each subgraph
	print('GETTING THE GRAPH FROM RAW SEQUENCE FILES BACK TO SUSPENSIONS')
	susps, links = seq_to_susp(links_dict)

	# walkback graph the rest of the way
	print('GETTING THE GRAPH FROM SUSPENSIONS BACK TO DONORS')
	seen = set()
	remaining = susps
	while remaining:
		seen.update(remaining)
		next_remaining = set()
		for identifier in remaining:
			url = urljoin(server, identifier + '/?format=json')
			temp_obj = requests.get(url, auth=connection.auth).json()
			get_object(temp_obj)
			get_derived_from(temp_obj, next_remaining, links)
		remaining = next_remaining - seen

	dir_to_make = ['', 'metadata', 'links', 'data', 'descriptors', 'descriptors/sequence_file']

	# the dcp does not capture cell_lines used just to grow organoids
	links = remove_cell_lines(links, whole_dict)
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
	print('WRITING THE JSON FILES')
	for d in dir_to_make:
		os.mkdir(dataset_id + '/' + d)

	# reformat links to use uuids and put in required schema
	for i in links:
		i['schema_type'] = 'links'
		i['schema_version'] = dcp_vs['links']
		i['describedBy'] = 'https://schema.humancellatlas.org/system/{}/links'.format(dcp_vs['links'])
		first_id = i['links'][0]['process_id']
		with open(dataset_id + '/links/' + first_id + '_' + dt + '_' + dataset_id + '.json', 'w') as outfile:
			json.dump(i, outfile, indent=4)
			outfile.close()

	dcp_types = {
		'project': 'project',
		'donor_organism': 'biomaterial',
		'specimen_from_organism': 'biomaterial',
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
		'supplementary_file': 'file'
	}

	s3_uris = []
	ftp_uris = []

	# write a json file for each object
	for k in whole_dict.keys():
		os.mkdir(dataset_id + '/metadata/' + k)
		for o in whole_dict[k]:
			if k == 'sequence_file':
				file_descript(o, k, dataset_id)
				if o.get('s3_uri'):
					s3_uris.append(o['s3_uri'])
					del o['s3_uri']
				elif o.get('external_uri'):
					ftp_uris.append(o['external_uri'])
					del o['external_uri']
			elif k == 'supplementary_file':
				file_name = o['file_core']['file_name']
				file_stats(file_name, o)
				os.rename(file_name, dataset_id + '/data/' + file_name)
				file_descript(o, k, dataset_id)
			customize_fields(o, k)
			o['schema_type'] = dcp_types[k].split('/')[0]
			o['schema_version'] = dcp_vs[k]
			o['describedBy'] = 'https://schema.humancellatlas.org/type/{}/{}/{}'.format(dcp_types[k], dcp_vs[k], k)
			with open(dataset_id + '/metadata/' + k + '/' + o['provenance']['document_id'] + '_' + dt + '.json', 'w') as outfile:
				json.dump(o, outfile, indent=4)

	# transfer the metadata directory to the DCP Google Cloud project
	request_to_gcp.local_dir_transfer(dataset_id)

	# transfer the data files from S3 to the DCP Google Cloud project
	if s3_uris:
		request_to_gcp.aws_file_transfer(dataset_id, s3_uris)

	# transfer the data files from external FTPs to the DCP Google Cloud project
	if ftp_uris:
		request_to_gcp.ftp_file_transfer(dataset_id, ftp_uris)

	for k,v in not_incl.items():
		not_incl[k] = list(v)
	with open('not_included.json', 'w') as outfile:
		json.dump(not_incl, outfile, indent=4)

	if not_valid:
		with open('not_validated.json', 'w') as outfile:
			outfile.write('\n'.join(not_valid))
			outfile.write('\n')

if __name__ == '__main__':
	d_now = datetime.now(tz=timezone.utc).isoformat(timespec='auto')
	dt = str(d_now).replace('+00:00', 'Z')

	if os.path.exists('not_included.json'):
		not_incl = json.load(open('not_included.json'))
	else:
		not_incl = {}

	whole_dict = {}
	not_valid = []
	doc_files = []
	handled_docs = []

	args = getArgs()
	if not args.dataset:
		sys.exit('ERROR: --dataset is required')
	if not args.mode:
		sys.exit('ERROR: --mode is required')
	if not args.dcp:
		sys.exit('ERROR: --dcp is required')

	connection = lattice.Connection(args.mode)
	server = connection.server
	dcp_vs = get_dcp_schema_ver(args.dcp)

	schema_errors = test_mapping()
	if schema_errors > 0:
		print(str(schema_errors) + ' found')
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to schema errors')

	main()
