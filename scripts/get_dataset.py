import argparse
import json
import lattice
import os
import requests
from urllib.parse import urljoin


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

def get_object_and_links(obj_id):
    url = urljoin(server, obj_id + '/?format=json&frame=object')
    temp_obj = requests.get(url, auth=connection.auth).json()

    # remove fields with empty values or unnecessary fields
    my_obj = {k: v for k, v in temp_obj.items() if v and k not in ignore_props}

    obj_type = my_obj['@type'][0]

    # add the object if we haven't already done so
    if whole_dict.get(obj_type) and my_obj.get('uuid') not in added_uuids:
        whole_dict[obj_type].append(my_obj)
        added_uuids.add(my_obj['uuid'])
    elif my_obj.get('uuid') not in added_uuids:
        whole_dict[obj_type] = [my_obj]
        added_uuids.add(my_obj['uuid'])

    # get identifiers for any object that referenced by this one
    for prop in my_obj.keys():
        if schemas[obj_type]['properties'][prop].get('linkTo'):
            next_remaining.add(my_obj[prop])
        elif schemas[obj_type]['properties'][prop].get('items'):
            if schemas[obj_type]['properties'][prop]['items'].get('linkTo'):
                next_remaining.update(my_obj[prop])

ignore_props = [
    'original_files', # these are all captured in files
    'fileset', # we don't need to grab the ReferenceFileSet objects
    'submitted_by' # our internal users aren't relevant
]

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

schema_url = urljoin(server, 'profiles/?format=json')
schemas = requests.get(schema_url).json()

dataset_id = args.dataset

whole_dict = {}
added_uuids = set() # we want to keep track to avoid duplicating objects
seen = set()
remaining = set([dataset_id])

while remaining:
    seen.update(remaining)
    next_remaining = set()
    for identifier in remaining:
        get_object_and_links(identifier)
    remaining = next_remaining - seen

# make directory named after dataset
os.mkdir(dataset_id)

# write json files for each object type in the dataset directory
for k in whole_dict.keys():
    with open(dataset_id + '/' + k + '.json', 'w') as outfile:
        json.dump(whole_dict[k], outfile, indent=4)
        outfile.close()
