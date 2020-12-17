import argparse
import gspread
import json
import lattice
import os
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

# grab all of the submittable properties
props = {}
schema_url = urljoin(server, 'profiles/' + schema_name + '/?format=json')
schema = requests.get(schema_url).json()
for prop in schema['properties'].keys():
	if not str(schema['properties'][prop].get('comment')).startswith('Do not submit') \
	and schema['properties'][prop].get('notSubmittable') != True:
		props[prop] = {}
		for i in schema['properties'][prop].keys():
			props[prop][i] = schema['properties'][prop][i]

ordered_props = OrderedDict(props)

# grab all of the properties of subobjects
subprops = {}
non_submit = [] # collect the base property so we can grey it out in favor of the subproperties
for prop in props.keys():
	if props[prop]['type'] == 'object' or \
	(props[prop]['type'] == 'array' and props[prop]['items']['type'] == 'object'):
		if props[prop]['type'] == 'array':
			my_props = props[prop]['items']['properties'] 
		else:
			my_props = props[prop]['properties']
		ordered_props.pop(prop)
		subprops[prop] = props[prop]
		non_submit.append(prop)
		for sp in my_props.keys():
			subprops[prop + '.' + sp] = {}
			for i in my_props[sp].keys():
				subprops[prop + '.' + sp][i] = my_props[sp][i]
ordered_props.update(subprops)

non_submit_col = []
for prop in non_submit:
	non_submit_col.append(cell_grid[list(ordered_props.keys()).index(prop) + 1])

# collect required fields & move fields to the front
req_props = []
if schema.get('required'):
	req_props = schema['required']
	for i in req_props:
		ordered_props.move_to_end(i, False)

# get the required field columns so we can color them later
req_columns = []
if req_props:
	if 'aliases' in ordered_props.keys():
		ordered_props.move_to_end('aliases', False)
		req_start_col = 'C'
		req_stop_col = cell_grid[len(req_props) + 1]
	else:
		req_start_col = 'B'
		req_stop_col = cell_grid[len(req_props)]
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
for prop in ordered_props.keys():
	prop_list.append(prop)
uber_list.append(prop_list)

# gather the attributes of each property
for descriptor in descriptor_list:
	this_list = ['#' + descriptor]
	for prop in ordered_props.keys():
		if ordered_props[prop]['type'] == 'array' and descriptor in ['type','enum','linkTo']:
			if ordered_props[prop]['items'].get(descriptor):
				this_list.append('array of ' + str(ordered_props[prop]['items'].get(descriptor,'')))
			else:
				this_list.append('')
		else:
			this_list.append(str(ordered_props[prop].get(descriptor,'')))
	uber_list.append(this_list)

# write the whole thing to the google sheet
tab.update('A1',uber_list)

# bold the first column
tab.format('A:A', {'textFormat': {'bold': True}})

# set the whole sheet to clip text
tab.format('A1:AZ100',{'wrapStrategy': 'CLIP'})

# set cell validation in the first input row for all boolean fields or fields with an enum list
count = 0
for prop in ordered_props.keys():			
	count += 1
	if ordered_props[prop].get('enum') or ordered_props[prop].get('type') == 'boolean':
		col = cell_grid[count] 
		cell_to_format = col + str(len(descriptor_list) + 2) + ':' + col + '100'
		validation_rule = DataValidationRule(
						    BooleanCondition('ONE_OF_LIST',
						    ordered_props[prop].get('enum', ['TRUE','FALSE'])),
						    showCustomUi=True
							)
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
