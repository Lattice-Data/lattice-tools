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

schema_name = args.type

# follow instructions here to enable API & generate credentials
# https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html
creds = ServiceAccountCredentials.from_json_keyfile_name(args.creds, 'https://www.googleapis.com/auth/drive')
client = gspread.authorize(creds)
sheet = client.open_by_key(args.sheet)
for tab in sheet.worksheets():
	if tab.title == schema_name:
		sheet.del_worksheet(tab)
tab = sheet.add_worksheet(title=schema_name,rows='100',cols='52')

abcs = string.ascii_uppercase
cell_grid = list(abcs) + ['A' + i for i in abcs]

connection = lattice.Connection(args.mode)
server = connection.server

# grab the OntologyTerm term_name & term_id schemas to put in places that linkTo OntologyTerm
ont_schema_url = urljoin(server, 'profiles/ontology_term/?format=json')
ont_schema = requests.get(ont_schema_url).json()
term_id_props = ont_schema['properties']['term_id']
term_name_props = ont_schema['properties']['term_name']

# grab all of the submittable properties
props = {}
schema_url = urljoin(server, 'profiles/' + schema_name + '/?format=json')
schema = requests.get(schema_url).json()
for p in schema['properties'].keys():
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
ordered_props.update(subprops)

remove_props = []
ont_props = []
for p in ordered_props.keys():
	if str(ordered_props[p].get('comment')).startswith('Do not submit') \
	or ordered_props[p].get('notSubmittable') == True:
		remove_props.append(p)
		if p in non_submit:
			non_submit.remove(p)
	elif ordered_props[p].get('linkTo') == 'OntologyTerm':
		remove_props.append(p)
		ont_props.append(p)
for p in remove_props:
	del ordered_props[p]
for p in ont_props:
	ordered_props[p + '.term_id'] = term_id_props
	ordered_props[p + '.term_name'] = term_name_props

non_submit_col = []
for p in non_submit:
	non_submit_col.append(cell_grid[list(ordered_props.keys()).index(p) + 1])

# collect required fields & move fields to the front
req_props = []
if schema.get('required'):
	req_count = 0
	req_props = schema['required']
	for i in req_props:
		if i in ordered_props:
			ordered_props.move_to_end(i, False)
			req_count += 1
		else:
			ordered_props.move_to_end(i + '.term_id', False)
			ordered_props.move_to_end(i + '.term_name', False)
			req_count += 2

# get the required field columns so we can color them later
req_columns = []
if req_props:
	if 'aliases' in ordered_props.keys():
		ordered_props.move_to_end('aliases', False)
		req_start_col = 'C'
		req_stop_col = cell_grid[req_count + 1]
	else:
		req_start_col = 'B'
		req_stop_col = cell_grid[req_count]
	req_columns = ':'.join([req_start_col, req_stop_col])

# list the attributes we want to know about each property
descriptor_list = [
	'title',
	'description',
	'comment',
	'type',
	'linkTo',
	'enum'
]

uber_list = []

# gather the top row list of schema_version followed by the property names
schema_version = schema['properties']['schema_version']['default']
prop_list = ['schema_version=' + schema_version]
for p in ordered_props.keys():
	prop_list.append(p)
uber_list.append(prop_list)

# gather the attributes of each property
for descriptor in descriptor_list:
	this_list = ['#' + descriptor]
	for p in ordered_props.keys():
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

# set cell validation in the first input row for all boolean fields or fields with an enum list
count = 0
for p in ordered_props.keys():			
	count += 1
	if ordered_props[p].get('enum') or ordered_props[p].get('type') == 'boolean':
		col = cell_grid[count] 
		cell_to_format = col + str(len(descriptor_list) + 2) + ':' + col + '100'
		validation_rule = DataValidationRule(BooleanCondition('ONE_OF_LIST',
							ordered_props[p].get('enum', ['TRUE','FALSE'])),
							showCustomUi=True)
		set_data_validation_for_cell_range(tab, cell_to_format, validation_rule)

# aliases should be the first property listed, so freeze that column and the descriptor column
if ordered_props.get('aliases'):
	set_frozen(tab, rows=len(descriptor_list) + 1, cols=2)
else: #if no aliases propertry, then just freeze the descriptor column
	set_frozen(tab, rows=len(descriptor_list) + 1, cols=1)

# shade all of the columns with required properties
if req_columns:	
	green = color(0.58, 0.77, 0.49)
	format_cell_range(tab, req_columns, cellFormat(backgroundColor=green))

# for the properties with embedded objects, shade the non-submittable property
for column in non_submit_col:
	grey = color(0.85, 0.85, 0.85)
	format_cell_range(tab, column, cellFormat(backgroundColor=grey))
