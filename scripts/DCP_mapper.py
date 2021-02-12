import argparse
import json
import lattice
import os
import requests
import sys
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
    args = parser.parse_args()
    return args


def get_object(temp_obj):
	obj_type = temp_obj['@type'][0]

	# remove fields with empty values or unnecessary fields
	my_obj = {}
	for k,v in temp_obj.items():
		if lattice_to_dcp[obj_type].get(k):
			dcp_prop = lattice_to_dcp[obj_type][k]
			if '.' in dcp_prop:
				path = dcp_prop.split('.')
				if my_obj.get(path[0]):
					my_obj[path[0]][path[1]] = v
				else:
					my_obj[path[0]] = {}
					my_obj[path[0]][path[1]] = v
			else:
				my_obj[dcp_prop] = v

	dcp_obj_type = lattice_to_dcp[obj_type]['class']

	# add the object if we haven't already done so
	if whole_dict.get(dcp_obj_type):
		whole_dict[dcp_obj_type].append(my_obj)
	else:
		whole_dict[dcp_obj_type] = [my_obj]


def get_derived_from(temp_obj):
	uuid = temp_obj['uuid']
	# get identifiers for any object that referenced by this one
	if temp_obj.get('derived_from'):
		if isinstance(temp_obj['derived_from'], str):
			next_remaining.add(temp_obj['derived_from'])
		elif isinstance(temp_obj['derived_from'], dict):
			next_remaining.add(temp_obj['derived_from']['@id'])
		elif isinstance(temp_obj['derived_from'], list):
			if isinstance(temp_obj['derived_from'][0], str):
				next_remaining.update(temp_obj['derived_from'])
			elif isinstance(temp_obj['derived_from'][0], dict):
				for o in temp_obj['derived_from']:
					next_remaining.add(o['@id'])

args = getArgs()
if not args.dataset:
	sys.exit('ERROR: --dataset is required')
if not args.mode:
	sys.exit('ERROR: --mode is required')

connection = lattice.Connection(args.mode)
server = connection.server

schema_url = urljoin(server, 'profiles/?format=json')
schemas = requests.get(schema_url).json()

whole_dict = {}
url = urljoin(server, args.dataset + '/?format=json')
ds_obj = requests.get(url, auth=connection.auth).json()
dataset_id = ds_obj['uuid']
get_object(ds_obj)

files = [i for i in ds_obj['files']]
libs = set()
for f in files:
	url = urljoin(server, f + '/?format=json')
	temp_obj = requests.get(url, auth=connection.auth).json()
	obj_type = temp_obj['@type'][0]
	if obj_type == 'RawSequenceFile':
		get_object(temp_obj)
		libs.update(temp_obj['libraries'])

seen = set()
remaining = libs # seed the dataset identifier to start
while remaining:
	seen.update(remaining)
	next_remaining = set()
	for identifier in remaining:
		url = urljoin(server, identifier + '/?format=json')
		temp_obj = requests.get(url, auth=connection.auth).json()
		get_object(temp_obj)
		get_derived_from(temp_obj)
	remaining = next_remaining - seen

# make directory named after dataset
os.mkdir(dataset_id)
os.mkdir(dataset_id + '/metadata')

# write json files for each object type in the dataset directory
for k in whole_dict.keys():
	os.mkdir(dataset_id + '/metadata/' + k)
	for o in whole_dict[k]:	
	    with open(dataset_id + '/metadata/' + k + '/' + o['provenance']['document_id'] + '.json', 'w') as outfile:
	        json.dump(o, outfile, indent=4)
	        outfile.close()
