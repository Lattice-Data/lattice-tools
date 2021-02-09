import argparse
import gspread
import json
import os
from gspread_formatting import *
from oauth2client.service_account import ServiceAccountCredentials


def getArgs():
    parser = argparse.ArgumentParser(
    	formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    parser.add_argument('-c','--creds',
                        help="the location of google drive client_secret.json file")
    parser.add_argument('-s','--sheet',
                        help="the key for the google sheet")
    args = parser.parse_args()
    return args

args = getArgs()
if not args.creds:
	sys.exit('ERROR: --creds is required')
if not args.sheet:
	sys.exit('ERROR: --sheet is required')

schemas_to_pull = [
	'donor_organism',
	'specimen_from_organism',
	'cell_line',
	'organoid',
	'cell_suspension',
	'library_preparation_protocol',
	'sequence_file'
]

# follow instructions here to enable API & generate credentials
# https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html
creds = ServiceAccountCredentials.from_json_keyfile_name(args.creds, 'https://www.googleapis.com/auth/drive')
client = gspread.authorize(creds)
sheet = client.open_by_key(args.sheet)

schema_dir = '/Users/jason/GitClones/HCA/metadata-schema/json_schema/'

main_obj_dir = 'type/'

versions = json.load(open(os.path.join(schema_dir, 'versions.json')))

descriptor_list = [
	'property',
	'user_friendly',
	'description',
	'comment',
	'type',
	'pattern',
	'$ref',
	'enum',
	'example',
	'items',
	'guidelines'
]

for path, subdirs, files in os.walk(schema_dir + main_obj_dir):
	for name in files:
		schema_name = name.split('.')[0]
		if schema_name in schemas_to_pull:
			non_submit = []
			all_props = {}
			schema = json.load(open(os.path.join(path, name)))
			for tab in sheet.worksheets():
				if tab.title == schema_name:
					sheet.del_worksheet(tab)
			tab = sheet.add_worksheet(title=schema_name,rows='100',cols='52')
			for prop in schema['properties'].keys():
				all_props[prop] = {}
				all_props[prop]['property'] = prop
				for k in schema['properties'][prop].keys():
					all_props[prop][k] = str(schema['properties'][prop][k])
				ref = None
				if '$ref' in schema['properties'][prop].keys():
					ref = schema['properties'][prop]['$ref']
				elif schema['properties'][prop].get('items'):
					if '$ref' in schema['properties'][prop]['items'].keys():
						ref = schema['properties'][prop]['items']['$ref']
				if ref:
					non_submit.append(prop)
					ref_path = ref.split('#')
					ref_file = ref_path[0]
					subschema = json.load(open(os.path.join(schema_dir, ref_file)))
					if len(ref_path) == 1:
						for subprop in subschema['properties'].keys():
							up = '{}.{}'.format(prop, subprop)
							all_props[up] = {}
							all_props[up]['property'] = up
							for k in subschema['properties'][subprop].keys():
								all_props[up][k] = str(subschema['properties'][subprop][k])
							ref2 = None
							if '$ref' in subschema['properties'][subprop].keys():
								ref2 = subschema['properties'][subprop]['$ref']
							elif subschema['properties'][subprop].get('items'):
								if '$ref' in subschema['properties'][subprop]['items'].keys():
									ref2 = subschema['properties'][subprop]['items']['$ref']
							if ref2:
								non_submit.append('{}.{}'.format(prop, subprop))
								ref2_path = ref2.split('#')
								ref2_file = ref2_path[0]
								subsubschema = json.load(open(os.path.join(schema_dir, ref2_file)))
								if len(ref2_path) == 1:
									for subsubprop in subsubschema['properties'].keys():
										up = '{}.{}.{}'.format(prop, subprop, subsubprop)
										all_props[up] = {}
										all_props[up]['property'] = up
										for k in subsubschema['properties'][subsubprop].keys():
											all_props[up][k] = str(subsubschema['properties'][subsubprop][k])
					else:
						subpath = ref_path[1].split('/')
						if len(subpath) == 3 and subpath[0] == '':
							for subprop in subschema[subpath[1]][subpath[2]]['properties'].keys():
								up = '{}.{}'.format(prop, subprop)
								all_props[up] = {}
								all_props[up]['property'] = up
								for k in subschema[subpath[1]][subpath[2]]['properties'][subprop].keys():
									all_props[up][k] = str(subschema[subpath[1]][subpath[2]]['properties'][subprop][k])
			uber_list = []
			uber_list.append(descriptor_list)
			for k in all_props.keys():
				next_list = []
				for i in descriptor_list:
					next_list.append(all_props[k].get(i))
				uber_list.append(next_list)

			tab.update('A1',uber_list)

			# bold the first column
			tab.format('A:A', {'textFormat': {'bold': True}})

			# set the whole sheet to clip text
			tab.format('A1:AZ100',{'wrapStrategy': 'CLIP'})

			set_frozen(tab, rows=1)

			non_submit_row = []
			for p in non_submit:
				non_submit_row.append(list(all_props.keys()).index(p) + 2)

			# for the properties with embedded objects, shade the non-submittable property
			for row in non_submit_row:
				cells = 'A' + str(row) + ':K' + str(row)
				grey = color(0.85, 0.85, 0.85)
				format_cell_range(tab, cells, cellFormat(backgroundColor=grey))
