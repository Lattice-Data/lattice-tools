import argparse
import gspread
import json
import lattice
import requests
import string
import sys
from collections import OrderedDict
from gspread_formatting import *
from oauth2client.service_account import ServiceAccountCredentials
from urllib.parse import urljoin


def getArgs():
    parser = argparse.ArgumentParser(
    	formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    parser.add_argument('-t','--type',
                        help="the object type to return a template for")
    parser.add_argument('-m','--mode',
                        help="the server to look-up schema, if not local")
    parser.add_argument('-c','--creds',
                        help="the location of google drive client_secret.json file")
    parser.add_argument('-s','--sheet',
                        help="the key for the google sheet")
    args = parser.parse_args()
    return args

args = getArgs()
if not args.type:
	sys.exit('ERROR: --type is required')
if not args.creds:
	sys.exit('ERROR: --creds is required')
if not args.sheet:
	sys.exit('ERROR: --sheet is required')
if not args.mode:
	sys.exit('ERROR: --mode is required')

system_props = [
	'status',
	'uuid',
	'schema_version',
	'date_created',
	'submitted_by',
	'accession',
	'alternate_accessions',
	'supersedes',
	'superseded_by',
	'notes',
	'aliases',
	'@id',
	'@type'
]

expand_me = {
	'LibraryProtocol': 'library_protocol',
	'OntologyTerm': 'ontology_term',
	'SequencingRun': 'sequencing_run'
}

all_schemas = [
	'HumanPostnatalDonor',
	'HumanPrenatalDonor',
	'MousePostnatalDonor',
	'MousePrenatalDonor',
	'Tissue',
	'CellCulture',
	'Organoid',
	'Suspension',
	'Library',
	'SequencingRun',
	'RawSequenceFile',
	'MatrixFile',
	'ReferenceFile',
	'CellAnnotation'
]

if args.type == 'all':
	schemas = all_schemas
else:
	schemas = [args.type]

# follow instructions here to enable API & generate credentials
# https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html
creds = ServiceAccountCredentials.from_json_keyfile_name(args.creds, 'https://www.googleapis.com/auth/drive')
client = gspread.authorize(creds)
sheet = client.open_by_key(args.sheet)
for schema in schemas:
	for tab in sheet.worksheets():
		if tab.title.lower() == schema.lower():
			sheet.del_worksheet(tab)
	tab = sheet.add_worksheet(title=schema,rows='100',cols='52')

	connection = lattice.Connection(args.mode)
	server = connection.server

	# grab all of the submittable properties
	props = {}
	schema_url = urljoin(server, 'profiles/' + schema + '/?format=json')
	schema = requests.get(schema_url).json()
	for p in schema['properties'].keys():
		if p not in system_props:
			props[p] = schema['properties'][p]

	ordered_props = OrderedDict(props)

	# grab all of the properties of subobjects
	subprops = {}
	non_submit = [] # collect the base property so we can grey it out in favor of the subproperties
	for p in props.keys():
		if props[p]['type'] == 'object' or \
		(props[p]['type'] == 'array' and props[p]['items']['type'] == 'object'):
			subprops[p] = props[p]
			ordered_props.pop(p)
			non_submit.append(p)
			if props[p]['type'] == 'array':
				for sp in props[p]['items']['properties'].keys():
					if props[p]['items']['properties'][sp]['type'] == 'object' or \
					(props[p]['items']['properties'][sp]['type'] == 'array' and props[p]['items']['properties'][sp]['items']['type'] == 'object'):
						subprops[p + '.' + sp] = props[p]['items']['properties'][sp]
						non_submit.append(p + '.' + sp)
						if props[p]['items']['properties'][sp]['type'] == 'array':
							for ssp in props[p]['items']['properties'][sp]['items']['properties'].keys():
								subprops[p + '.' + sp + '.' + ssp] = props[p]['items']['properties'][sp]['items']['properties'][ssp]
						else:
							for ssp in props[p]['items']['properties'][sp]['items']['properties'].keys():
								subprops[p + '.' + sp + '.' + ssp] = props[p]['items']['properties'][sp]['properties'][ssp]
					else:
						subprops[p + '.' + sp] = props[p]['items']['properties'][sp]
			else:
				my_props = props[p]['properties']
				for sp in my_props.keys():
					subprops[p + '.' + sp] = my_props[sp]
		elif isinstance(props[p].get('linkTo'), str):
			if props[p].get('linkTo','') in expand_me.keys():
				subprops[p] = props[p]
				ordered_props.pop(p)
				non_submit.append(p)
				embed_obj = expand_me[props[p].get('linkTo','')]
				embed_url = urljoin(server, 'profiles/' + embed_obj + '/?format=json')
				embed_schema = requests.get(embed_url).json()
				for ep in embed_schema['properties'].keys():
					if ep not in system_props:
						subprops[p + '.' + ep] = embed_schema['properties'][ep]
		elif props[p]['type'] == 'array':
			if isinstance(props[p]['items'].get('linkTo'), str):
				if props[p]['items'].get('linkTo','') in expand_me.keys():
					subprops[p] = props[p]
					ordered_props.pop(p)
					non_submit.append(p)
					embed_obj = expand_me[props[p]['items'].get('linkTo','')]
					embed_url = urljoin(server, 'profiles/' + embed_obj + '/?format=json')
					embed_schema = requests.get(embed_url).json()
					for ep in embed_schema['properties'].keys():
						if ep not in system_props:
							subprops[p + '.' + ep] = embed_schema['properties'][ep]
	ordered_props.update(subprops)

	non_submit_row = []
	for p in non_submit:
		non_submit_row.append(list(ordered_props.keys()).index(p) + 2)

	# list the attributes we want to know about each property
	descriptor_list = [
		'title',
		'description',
		'comment',
		'type',
		'linkTo',
		'enum',
		'pattern'
	]

	uber_list = [['property'] + descriptor_list]

	# gather the attributes of each property
	for p in ordered_props.keys():
		this_list = [p]
		for descriptor in descriptor_list:
			if ordered_props[p]['type'] == 'array' and descriptor in ['type','enum','linkTo']:
				if ordered_props[p]['items'].get(descriptor):
					this_list.append('array of ' + str(ordered_props[p]['items'].get(descriptor,'')))
				else:
					this_list.append('')
			else:
				this_list.append(str(ordered_props[p].get(descriptor,'')))
		uber_list.append(this_list)

	# write the whole thing to the google sheet
	tab.update('A1',uber_list)

	# bold the first column
	tab.format('A:A', {'textFormat': {'bold': True}})

	# set the whole sheet to clip text
	tab.format('A1:AZ100',{'wrapStrategy': 'CLIP'})

	set_frozen(tab, rows=1)

	# for the properties with embedded objects, shade the non-submittable property
	for row in non_submit_row:
		cells = 'A' + str(row) + ':H' + str(row)
		grey = color(0.85, 0.85, 0.85)
		format_cell_range(tab, cells, cellFormat(backgroundColor=grey))
