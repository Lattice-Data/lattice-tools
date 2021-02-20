import argparse
import hashlib
import json
import lattice
import os
import requests
import sys
from datetime import datetime, timezone
from urllib.parse import urljoin
from property_mapping import (
    lattice_to_dcp
)


EPILOG = '''
Get a Dataset from the Lattice database and put into json files.

Examples:

    python %(prog)s --dataset LATDS444YHX

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
	schema_url = urljoin(server, 'profiles/?format=json')
	schemas = requests.get(schema_url).json()
	for k,v in lattice_to_dcp.items():
		schema_props = schemas[k]['properties']
		for prop in v.keys():
			if prop != 'class':
				lat_prop = v[prop]['lattice']
				if '.' in lat_prop:
					path = lat_prop.split('.')
					if path[0] not in schema_props:
						print(k)
						print(lat_prop)
						print('--------')
				elif lat_prop not in schema_props:
					print(k)
					print(lat_prop)
					print('--------')
	sys.exit()


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
		if isinstance(v, list) and prop_map.get('subprop_map'):
			new_v = []
			for i in v:
				new_i = {}
				for subprop in prop_map['subprop_map'].keys():
					lat_subprop = prop_map['subprop_map'][subprop]['lattice']
					if i.get(lat_subprop):
						add_value(new_i, i, prop_map['subprop_map'][subprop], subprop)
				new_v.append(new_i)
				v = new_v
		elif isinstance(v, dict) and prop_map.get('subprop_map'):
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
	if not_incl.get(prop):
		not_incl[prop].add(str(value))
	else:
		not_incl[prop] = {str(value)}


def get_object(temp_obj):
	# drop unneeded properties, flatten subobjects
	temp_obj = flatten_obj(temp_obj)

	remove = set()
	my_obj = {}
	obj_type = temp_obj['@type'][0]

	for prop in lattice_to_dcp[obj_type].keys():
		if prop != 'class':
			lat_prop = lattice_to_dcp[obj_type][prop]['lattice']
			if temp_obj.get(lat_prop):
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
			if isinstance(v, dict) and k not in ['test_results', 'award']:
				for x,y in v.items():
					new_obj[k + '.' + x] = y
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
			lib = sr_obj['derived_from']
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
			'outputs': links_dict[seqruns]
			}
		link_hash = hashlib.md5(str(l).encode('utf-8')).hexdigest()
		l['link_type'] = 'process_link'
		l['process_type'] = 'process'
		l['process_id'] = link_hash
		link = {'links': [l]}
		links.append(link)
	return all_susps, links


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
	elif out_type == 'cell_suspension' and in_type != 'cell_suspension':
		pr_type = 'dissociation_protocol'
	elif out_type == 'organoid' and in_type == 'cell_line':
		pr_type = 'differentiation_protocol'
	else:
		pr_type = 'protocol'
	pr_id = hashlib.md5(str([pr_type + in_type + out_type]).encode('utf-8')).hexdigest()

	prots.append({'protocol_type': pr_type, 'protocol_id': pr_id})
	
	my_obj['protocol_core']['protocol_id'] = pr_id
	my_obj['provenance']['document_id'] = pr_id
	if whole_dict.get(pr_type):
		whole_dict[pr_type].append(my_obj)
	else:
		whole_dict[pr_type] = [my_obj]
	
	if 'enrichment_factors' in out_obj:
		pr_type = 'enrichment_protocol'
		enpr_id = hashlib.md5(str([pr_type + in_type + out_type]).encode('utf-8')).hexdigest()
		prots.apppend({
			'protocol_type': pr_type,
			'protocol_id': enpr_id
			})
		enr_obj = {
			'provenance': {
				'document_id': enpr_id
			},
			'protocol_core': {
				'protocol_id': enpr_id
				},
			'method': {
				'text': ','.join(der_process)
				},
			'enirchment_markers': out_type[pr_type]
			}
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
	link_hash = hashlib.md5(str(link).encode('utf-8')).hexdigest()
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
	v_file = directory + 'metadata-schema/json_schema/versions.json'
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


def customize_fields(obj, obj_type):
	if obj_type == 'project':
		if obj.get('project_core'):
			if obj['project_core'].get('project_title'):
				short = obj['project_core']['project_title'].replace(' ','')[:50]
				obj['project_core']['project_short_name'] = short
		if obj.get('publications'):
			for i in obj['publications']:
				index = obj['publications'].index(i)
				if i.get('doi'):
					obj['publications'][index]['doi'] = i['doi'][0]
				if i.get('pmid'):
					obj['publications'][index]['pmid'] = int(i['pmid'][0])
		if obj.get('funders'):
			obj['funders']['organization'] = 'Chan Zuckerberg Initiative'
			obj['funders'] = [obj['funders']]
	elif obj_type == 'donor_organism':
		if not obj.get('is_living'):
			if obj['development_stage']['ontology_label'] in ['embryonic','fetal']:
				obj['is_living'] = 'not applicable'
			else:
				obj['is_living'] = 'unknown'
		if obj.get('human_specific'):
			if obj['human_specific'].get('ethnicity'):
				obj['human_specific']['ethnicity'] = [obj['human_specific']['ethnicity']]
		if obj.get('medical_history'):
			if obj['medical_history'].get('test_results'):
				tr = [k + ':' + v for k,v in obj['medical_history']['test_results'].items()]
				obj['medical_history']['test_results'] = ','.join(tr)
	elif obj_type == 'specimen_from_organism':
		if obj.get('purchased_specimen'):
			if not obj['purchased_specimen'].get('catalog_number'):
				del obj['purchased_specimen']
	elif obj_type == 'cell_line':
		if obj.get('supplier') and not obj.get('catalog_number'):
			del obj['supplier']
	elif obj_type == 'cell_suspension':
		if obj.get('estimated_count_units'):
			if obj['estimated_count_units'] not in ['cells', 'nuclei']:
				del obj['estimated_cell_count']
			del obj['estimated_count_units']
	elif obj_type == 'library_preparation_protocol':
		if not obj.get('strand'):
			obj['strand'] = 'not provided'


def main():
	# get the dataset, and convert it to a project object
	url = urljoin(server, args.dataset + '/?format=json')
	ds_obj = requests.get(url, auth=connection.auth).json()
	dataset_id = ds_obj['uuid']
	get_object(ds_obj)

	# get the raw sequence files from that dataset
	links_dict = {}
	files = [i for i in ds_obj['files']]
	for f in files:
		url = urljoin(server, f + '/?format=json')
		temp_obj = requests.get(url, auth=connection.auth).json()
		obj_type = temp_obj['@type'][0]
		if obj_type == 'RawSequenceFile':
			# convert each to a sequence_file
			get_object(temp_obj)
			# pull the derived_from to store for later formation to links
			get_links(temp_obj, tuple(temp_obj['derived_from']), links_dict)

	# special walkback of graph until Suspension object
	# gather all the Suspension objects to traverse next
	# set up links between sequence_file and suspension as the start of each subgraph
	susps, links = seq_to_susp(links_dict)

	# walkback graph the rest of the way
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

	d_now = datetime.now(tz=timezone.utc).isoformat(timespec='auto')
	dt = str(d_now).replace('+00:00', 'Z')

	# make directory named after dataset
	os.mkdir(dataset_id)
	os.mkdir(dataset_id + '/metadata')
	os.mkdir(dataset_id + '/links/')

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
		'enrichment_protcol': 'protocol/biomaterial_collection',
		'differentiation_protocol': 'protocol/biomaterial_collection',
		'protocol': 'protocol'
	}

	# write a json file for each object
	for k in whole_dict.keys():
		os.mkdir(dataset_id + '/metadata/' + k)
		for o in whole_dict[k]:
			customize_fields(o, k)
			o['schema_type'] = dcp_types[k].split('/')[0]
			o['schema_version'] = dcp_vs[k]
			o['describedBy'] = 'https://schema.humancellatlas.org/type/{}/{}/{}'.format(dcp_types[k], dcp_vs[k], k)
			with open(dataset_id + '/metadata/' + k + '/' + o['provenance']['document_id'] + '_' + dt + '.json', 'w') as outfile:
				json.dump(o, outfile, indent=4)
				outfile.close()

	for k,v in not_incl.items():
		not_incl[k] = list(v)
	with open(dataset_id + '/not_included.json', 'w') as outfile:
		json.dump(not_incl, outfile, indent=4)

if __name__ == '__main__':
	whole_dict = {}
	not_incl = {}
	args = getArgs()
	if not args.dataset:
		sys.exit('ERROR: --dataset is required')
	if not args.mode:
		sys.exit('ERROR: --mode is required')
	connection = lattice.Connection(args.mode)
	server = connection.server
	dcp_vs = get_dcp_schema_ver(args.dcp)
	main()
