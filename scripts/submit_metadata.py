import argparse
import ast
import os
import json
import lattice
import mimetypes
import re
import requests
import sys
import magic  # install me with 'pip install python-magic'
import numpy as np
import pandas as pd
from PIL import Image  # install me with 'pip install Pillow'
from urllib.parse import urljoin, quote, unquote
from urllib.request import Request,urlopen
from base64 import b64encode
from bs4 import BeautifulSoup


EPILOG = '''
Thorough instructions have been documented in lattice-tools/docs/submit_metadata.md

For more details:

		%(prog)s --help
'''


# update this from encoded/src/encoded/loadxl.py as needed
ORDER = [
	'user',
	'award',
	'lab',
	'organism',
	'gene',
	'publication',
	'document',
	'ontology_term',
	'library_protocol',
	'target',
	'antibody',
	'treatment',
	'human_postnatal_donor',
	'human_prenatal_donor',
	'mouse_postnatal_donor',
	'mouse_prenatal_donor',
	'tissue',
	'cell_culture',
	'organoid',
	'suspension',
	'tissue_section',
	'dataset',
	'library',
	'sequencing_run',
	'raw_sequence_file',
	'sequence_alignment_file',
	'raw_matrix_file',
	'processed_matrix_file',
	'rna_metrics',
	'antibody_capture_metrics',
	'rna_aggregate_metrics',
	'atac_metrics',
	'atac_aggregate_metrics',
	'multiome_metrics',
	'spatial_metrics',
	'cell_annotation',
	'image',
	'page',
	'access_key'
]


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	parser.add_argument('sheet_id',
						help='the ID of the google sheet containing object data to import')
	parser.add_argument('--justtype', '-jt',
						help='the type of the objects to import if only one type is desired')
	parser.add_argument('--starttype', '-st',
						help='the type of the objects to start importing and continue in sequence')
	parser.add_argument('--mode', '-m',
						help='The machine to run on.',
						required=True)
	parser.add_argument('--debug', '-d',
						default=False,
						action='store_true',
						help='Print debug messages.  Default is False')
	parser.add_argument('--update',
						default=False,
						action='store_true',
						help='Let the script proceed with the changes.  Default is False'),
	parser.add_argument('--patchall', '-p',
						default=False,
						action='store_true',
						help='PATCH existing objects.  Default is False \
						and will only PATCH with user override')
	parser.add_argument('--remove',
						default=False,
						action='store_true',
						help='Will remove all values in the provided properties.  Default is False'),
	args = parser.parse_args()
	return args


def report_schema_error(prop, schema):
	print('ERROR: Property "{}" not found in schema {}'.format(prop, schema))


def properties_validator(keys, schema_name, schema, dup_cols, remove, patchall):
	flag = False
	prop_types = {}
	dup_keys = {}
	base_props = []
	linkTos = []
	ont_props = ['term_name','term_id']
	schema_properties = schema['properties']

	for key in keys:
		props = key.split('.')
		if '.'.join(props[:-1]) in dup_keys and props[-1] == str(dup_keys['.'.join(props[:-1])]):
			flag = True #this is the second or more instance of a duplicated header
			dup_keys['.'.join(props[:-1])] += 1
		else:
			if key in dup_cols: # check for duplicated headers
				print('Property {} found multiple times in {} sheet headers'.format(key, schema_name))
				dup_keys[key] = 1
				flag = True
			propA = props[0].split('-')[0]
			base_props.append(propA)
			if propA not in schema_properties.keys():
				report_schema_error(propA, schema_name)
				flag = True
			elif len(props) == 1:
				immediate_schema = schema_properties[propA]
				prop_types[key] = immediate_schema['type']
				if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
					linkTos.append(key)
			elif len(props) > 1:
				propB = props[1].split('-')[0]
				if schema_properties[propA]['type'] == 'array':
					if propB not in schema_properties[propA]['items']['properties'].keys():
						report_schema_error(propB, schema_name + '.' + propA)
						flag = True
					elif len(props) == 2:
						immediate_schema = schema_properties[propA]['items']['properties'][propB]
						prop_types[key] = immediate_schema['type']
						if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
							linkTos.append(key)
					elif len(props) > 2:
						propC = props[2].split('-')[0]
						if schema_properties[propA]['items']['properties'][propB]['type'] == 'array':
							if propC not in schema_properties[propA]['items']['properties'][propB]['items']['properties'].keys():
								report_schema_error(propC, schema_name + '.' + propA + '.' + propB)
								flag = True
							else:
								immediate_schema = schema_properties[propA]['items']['properties'][propB]['items']['properties'][propC]
								prop_types[key] = immediate_schema['type']
								if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
									linkTos.append(key)
						elif schema_properties[propA]['items']['properties'][propB]['type'] == 'object':
							if propC not in schema_properties[propA]['items']['properties'][propB]['properties'].keys():
								report_schema_error(propC, schema_name + '.' + propA + '.' + propB)
								flag = True
							else:
								immediate_schema = schema_properties[propA]['items']['properties'][propB]['properties'][propC]
								prop_types[key] = immediate_schema['type']
								if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
									linkTos.append(key)
						elif schema_properties[propA]['items']['properties'][propB].get('linkTo') == "OntologyTerm":
							if propC not in ont_props:
								report_schema_error(propC, 'OntologyTerm via ' + schema_name + '.' + propA + '.' + propB)
								flag = True
							elif propC == 'term_id':
								prop_types[key] = 'ontology.term_id' #DO WE NEED LINKTOS HERE?
							else:
								prop_types[key] = 'string' #because term_name is strings #DO WE NEED LINKTOS HERE?
						else:
							print('Not expecting subproperties for property {}.{} in schema {} ({} given) '.format(propA, propB, schema_name, propC))
							flag = True
				elif schema_properties[propA]['type'] == 'object':
					if propB not in schema_properties[propA]['properties'].keys() and schema_properties[propA].get('additionalProperties') != True:
						report_schema_error(propB, schema_name + '.' + propA)
						flag = True
					elif len(props) == 2:
						if propB in schema_properties[propA]['properties'].keys():
							immediate_schema = schema_properties[propA]['properties'][propB]
							prop_types[key] = immediate_schema['type']
							if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
								linkTos.append(key)
						elif schema_properties[propA].get('additionalProperties') == True:
							prop_types[key] = 'string'
					elif len(props) > 2:
						propC = props[2].split('-')[0]
						if schema_properties[propA]['properties'][propB]['type'] == 'array':
							if propC not in schema_properties[propA]['properties'][propB]['items']['properties'].keys():
								report_schema_error(propC, schema_name + '.' + propA + '.' + propB)
								flag = True
							else:
								immediate_schema = schema_properties[propA]['properties'][propB]['items']['properties'][propC]
								prop_types[key] = immediate_schema['type']
								if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
									linkTos.append(key)
						elif schema_properties[propA]['properties'][propB]['type'] == 'object':
							if propC not in schema_properties[propA]['properties'][propB]['properties'].keys():
								report_schema_error(propC, schema_name + '.' + propA + '.' + propB)
								flag = True
							else:
								immediate_schema = schema_properties[propA]['properties'][propB]['properties'][propC]
								prop_types[key] = immediate_schema['type']
								if immediate_schema.get('linkTo') or immediate_schema.get('items',{}).get('linkTo'):
									linkTos.append(key)
						else:
							print('Not expecting subproperties for property {} in schema {}.{} ({} given) '.format(propB, schema_name, propA, propC))
							flag = True
				elif schema_properties[propA].get('linkTo') == "OntologyTerm":
					if propB not in ont_props:
						report_schema_error(propB, 'OntologyTerm via ' + schema_name + '.' + propA)
						flag = True
					else:
						prop_types[key] = 'string' #because term_id and term_name are both strings #DO WE NEED LINKTOS HERE?
				else:
					print('Not expecting subproperties for property {} in schema {} ({} given) '.format(propA, schema_name, propB))
					flag = True

	if not remove and not patchall:
		req_missing = [p for p in schema['required'] if p not in base_props]
		if req_missing:
			print('Required properties {} missing in schema {} '.format(','.join(req_missing), schema_name))
			flag = True

	return flag, prop_types, linkTos


def dict_patcher(old_dict, schema_properties):
	new_dict = {}
	array_o_objs_dict = {}
	onts = {}
	for key, value in old_dict.items():
		path = key.split('.')
		if len(path) == 1:
			if schema_properties[key].get('linkTo') == 'OntologyTerm' or \
				(schema_properties[key].get('items') and schema_properties[key]['items'].get('linkTo') == 'OntologyTerm'):
				if isinstance(value, list):
					value = [unquote(v).replace(':','_') for v in value]
				else:
					value = unquote(value).replace(':','_')
			elif key == 'term_id':
				value = value.replace('_',':')
			new_dict[key] = value
		elif len(path) == 2: # embedded object, need to build the mini dictionary to put this in
			if '-' in path[0]: # this has a number next to it and we expect an array of objects
				group = path[0].split('-')[1]
			else:
				group = '1'
			propA = path[0].split('-')[0]
			propB = path[1]
			if schema_properties[propA].get('linkTo') == 'OntologyTerm':
				if onts.get(propA):
					onts[propA][propB] = value
				else:
					onts[propA] = {}
					onts[propA][propB] = value
				if propB == 'term_id':
					new_dict[propA] = value.replace(':','_') # replace the term_name/term_id with the linkTo
			elif schema_properties[propA]['type'] == 'array': # this is a single object that needs to end as a 1 item array of objects
				if array_o_objs_dict.get(propA): # I have already added the embedded object to the new dictionary
					if array_o_objs_dict[propA].get(group):
						# this should be that we have started putting in the new object
						array_o_objs_dict[propA][group].update({propB: value})
					else:
						# this means we have not added any part of new item to the list
						array_o_objs_dict[propA][group] = {propB: value}
				else: # I need to make new item in dictionary
					array_o_objs_dict[propA] = {}
					array_o_objs_dict[propA][group] = {propB: value}
			else: # there is no number next to it, this is our only object to handle
				if new_dict.get(propA): # I have already added the embedded object to the new dictionary
					new_dict[propA].update({propB: value})
				else: # I need to make new item in dictionary
					temp_dict = {propB: value}
					new_dict[propA] = temp_dict
		elif len(path) == 3: # we have an object embedded within an embedded object
			if '-' in path[0]:
				group = path[0].split('-')[1]
			else:
				group = '1'
			if '-' in path[1]:
				subgroup = path[1].split('-')[1]
			else:
				subgroup = '1'
			propA = path[0].split('-')[0]
			propB = path[1].split('-')[0]
			propC = path[2]
			if schema_properties[propA]['type'] == 'array': # this is a single object that needs to end as a 1 item array of objects
				if schema_properties[propA]['items']['properties'][propB]['type'] == 'array': # we have an array of objects within an array of objects
					if array_o_objs_dict.get(propA):
						if array_o_objs_dict[propA].get(group):
							if array_o_objs_dict[propA][group].get(propB):
								if array_o_objs_dict[propA][group][propB].get(subgroup):
									array_o_objs_dict[propA][group][propB][subgroup].update({propC: value})
								else:
									array_o_objs_dict[propA][group][propB][subgroup] = {propC: value}
							else:
								array_o_objs_dict[propA][group][propB] = {}
								array_o_objs_dict[propA][group][propB][subgroup] = {propC: value}
						else:
							array_o_objs_dict[propA][group] = {}
							array_o_objs_dict[propA][group][propB] = {}
							array_o_objs_dict[propA][group][propB][subgroup] = {propC: value}
					else:
						array_o_objs_dict[propA] = {}
						array_o_objs_dict[propA][group] = {}
						array_o_objs_dict[propA][group][propB] = {}
						array_o_objs_dict[propA][group][propB][subgroup] = {propC: value}
				else: # this is an object within an array of objects
					if schema_properties[propA]['items']['properties'][propB].get('linkTo') == 'OntologyTerm':
						if onts.get(propB):
							onts[propB][propC] = value
						else:
							onts[propB] = {}
							onts[propB][propC] = value
						if propC == 'term_id':
							if array_o_objs_dict.get(propA):
								if array_o_objs_dict[propA].get(group):
									array_o_objs_dict[propA][group][propB] = value.replace(':','_') # replace the term_name/term_id with the linkTo
								else:
									array_o_objs_dict[propA][group] = {}
									array_o_objs_dict[propA][group][propB] = value.replace(':','_') # replace the term_name/term_id with the linkTo
							else:
								array_o_objs_dict[propA] = {}
								array_o_objs_dict[propA][group] = {}
								array_o_objs_dict[propA][group][propB] = value.replace(':','_') # replace the term_name/term_id with the linkTo
					else:
						print('ATTENTION IS NEEDED HERE')
			else:
				if schema_properties[propA][propB]['type'] == 'array': # this is an array of objects within an object
					print('ATTENTION IS NEEDED HERE')
				else: # this is an object within an object
					if new_dict.get(propA): # I have already added the embedded object to the new dictionary
						if new_dict[propA].get(propB): # I have already added the double-embedded object to the new dictionary
							new_dict[propA][propB].update({propC: value})
						else: # I need to add the double-embedded object to the new dictionary
							temp_dict = {propC: value}
							new_dict[propA].update({propB: temp_dict})
					else: # I need to make new item in dictionary
						temp_dict = {propC: value}
						new_dict[propA] = {propB: temp_dict}
		elif len(path) > 3:
			print('ERROR:' + key + ': Not prepared to do triple-embedded properties, check for errant period')

	for prop in array_o_objs_dict.keys():
		new_list = []
		for group_id in array_o_objs_dict[prop].keys():
			for i in array_o_objs_dict[prop][group_id].keys():
				if isinstance(array_o_objs_dict[prop][group_id][i], dict):
					nn_list = []
					for g_id in array_o_objs_dict[prop][group_id][i].keys():
						nn_list.append(array_o_objs_dict[prop][group_id][i][g_id])
					array_o_objs_dict[prop][group_id][i] = nn_list
			new_list.append(array_o_objs_dict[prop][group_id])
		new_dict[prop] = new_list
	return new_dict, onts


def attachment(path):
	''' Create an attachment upload object from a filename
	Embeds the attachment as a data url.
	'''
	if not os.path.isfile(path):
		r = requests.get(path)
		path = path.split('/')[-1]
		with open(path, 'wb') as outfile:
			outfile.write(r.content)
	filename = os.path.basename(path)
	mime_type, encoding = mimetypes.guess_type(path)
	major, minor = mime_type.split('/')
	detected_type = magic.from_file(path, mime=True)

	# XXX This validation logic should move server-side.
	if not (detected_type == mime_type or
			detected_type == 'text/plain' and major == 'text'):
		raise ValueError('Wrong extension for %s: %s' %
						 (detected_type, filename))

	with open(path, 'rb') as stream:
		attach = {
			'download': filename,
			'type': mime_type,
			'href': 'data:%s;base64,%s' % (mime_type, b64encode(stream.read()).decode('ascii'))
		}

		if mime_type in ('application/pdf', 'text/plain', 'text/tab-separated-values', 'text/html'):
			# XXX Should use chardet to detect charset for text files here.
			return attach

		if major == 'image' and minor in ('png', 'jpeg', 'gif', 'tiff'):
			# XXX we should just convert our tiffs to pngs
			stream.seek(0, 0)
			im = Image.open(stream)
			im.verify()
			if im.format != minor.upper():
				msg = 'Image file format %r does not match extension for %s'
				raise ValueError(msg % (im.format, filename))

			attach['width'], attach['height'] = im.size
			return attach

	raise ValueError('Unknown file type for %s' % filename)


def get_tab_ids(soup):
	tab_ids = {}
	pattern = re.compile('var bootstrapData = (.*?)};')
	for s in soup.find_all('script'):
		if pattern.search(str(s)):
			d = pattern.search(str(s)).group()[20:-1]
			data = json.loads(d)
			for t in data['changes']['topsnapshot']:
				u = t[1].split('"')
				if len(u) > 5:
					tab_ids[u[5]] = u[1]
	return tab_ids


def booleanify(s):
	if str(s).lower() in ['true', '1']:
		return True
	elif str(s).lower() in ['false', '0']:
		return False
	else:
		return s


def stringify(s):
	s = str(s).strip()
	if s.lower() == 'true':
		return 'True'
	elif s.lower() == 'false':
		return 'False'
	else:
		return s


def set_value_types(df, prop_types, linkTos):
	for c in list(df.columns):
		val_type = prop_types[c]
		if val_type == 'string':
			if c in linkTos:
				df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else quote(str(x[c]).strip()), axis=1)
			else:
				df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else stringify(x[c]), axis=1)
		elif val_type == 'ontology.term_id':
			df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else str(x[c]).replace('_',':'), axis=1)
		elif val_type == 'array':
			if c in linkTos:
				df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else [quote(v.strip()) for v in str(x[c]).split(',')], axis=1)
			else:
				df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else [v.strip() for v in str(x[c]).split(',')], axis=1)
				# if array of objects, 	old_value = old_value.replace('},{', '}${') THEN [ast.literal_eval(x) for x in old_value.split('$')]
		elif val_type == 'boolean':
			df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else booleanify(x[c]), axis=1)
		elif val_type == 'integer':
			df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else int(float(str(x[c]).replace(',',''))), axis=1).astype('Int64')
		elif val_type == 'number':
			df[c] = df.apply(lambda x: np.nan if pd.isnull(x[c]) else float(str(x[c]).replace(',','')), axis=1)


def check_existing_obj(post_json, schema, connection):
	if post_json.get('uuid'):
		temp_identifier = post_json['uuid']
	elif post_json.get('accession'):
		temp_identifier = post_json['accession']
	elif post_json.get('external_accession'):
		temp_identifier = post_json['external_accession']
	elif post_json.get('name'):
		temp_identifier = schema + '/' + post_json['name']
	elif post_json.get('term_id'):
		temp_identifier = schema + '/' + post_json['term_id'].replace(':','_')
	elif post_json.get('aliases'):
		temp_identifier = quote(post_json['aliases'][0])
	try:
		temp_identifier
	except NameError:
		temp = {}
		temp_identifier = ''
	else:
		url = urljoin(connection.server, temp_identifier + '/?format=json')
		temp = requests.get(url, auth=connection.auth).json()
	return temp_identifier, temp


def main():
	args = getArgs()

	connection = lattice.Connection(args.mode)
	server = connection.server
	print(f'Running on {server}')

	if args.remove and args.patchall:
		sys.exit('ERROR: cannot specify both --remove and --patch')

	if args.update:
		if args.patchall:
			print('Running to UPDATE, data will be PATCHED')
		elif args.remove:
			print('Running to UPDATE, data will be REMOVED')
		else:
			print('Running to UPDATE, data will be POSTED')
	else:
		print('Running in dry-run mode. No updates will be made')

	sheet_url = f'https://docs.google.com/spreadsheets/d/{args.sheet_id}'
	req = Request(sheet_url, headers={'User-Agent' : "Magic Browser"})
	s = urlopen(req)
	soup = BeautifulSoup(s, 'html.parser')
	h1 = soup.find('h1')
	if (h1 and h1.text == 'Sign in'):
		sys.exit('Google Sheet is not accessible. Ensure anyone with link can View.')
	tabs = [x.text.lower() for x in soup.find_all('div', {'class':'goog-inline-block docs-sheet-tab-caption'})]

	#establish what objects are being loaded and in what order
	load_order = ORDER
	if args.starttype:
		if args.justtype:
			sys.exit('ERROR: cannot specify both --justtype and --starttype')
		else:
			st_index = load_order.index(args.starttype)
			load_order = load_order[st_index:]
	if args.justtype:
		if args.justtype in load_order:
			load_order = [args.justtype]
		else:
			sys.exit('ERROR: --justtype {jt} not valid object type'.format(jt=args.justtype))

	not_schema = [t for t in tabs if t not in load_order]
	print('Sheets will NOT be considered: ' + ','.join(not_schema))

	to_load = [t for t in load_order if t in tabs]
	print('Sheets WILL be considered: ' + ','.join(to_load))

	#get the tab ids in order to fill in csv access links
	tab_ids = get_tab_ids(soup)

	all_posts = {}
	summary_report = []

	for schema_to_load in to_load:
		#read the sheet from csv into pandas
		g_id = tab_ids[schema_to_load]
		g_url = f'https://docs.google.com/spreadsheets/d/{args.sheet_id}/export?gid={g_id}&format=csv'
		df = pd.read_csv(g_url, dtype=str)

		#remove columns that are all NaNs
		if not args.remove:
			df = df.dropna(how='all',axis='columns')

		#get rid of rows that start with "#"
		df = df[(df.iloc[:,0].str.startswith('#') == False) | (df.iloc[:,0].isna())]

		#get rid of columns with a header that start with "#"
		drop = [h for h in df.columns if h.startswith('#')]
		df = df.drop(columns=drop)

		#remove columns that are all NaNs now that possible notes have been cleared
		if not args.remove:
			df = df.dropna(how='all',axis='columns')

		if df.empty:
			print('skipping {}: no objects to submit'.format(schema_to_load))
			tabs.remove(schema_to_load)
			continue

		#check for duplicated headers
		dup_cols = []
		for k,v in pd.read_csv(g_url, header=None, nrows=1).iloc[0,:].value_counts().items():
			if v > 1:
				dup_cols.append(k)

		headers = list(df.columns)
		schema_url = urljoin(server, 'profiles/' + schema_to_load + '/?format=json')
		schema = requests.get(schema_url).json()

		invalid_flag, prop_types, linkTos = properties_validator(headers, schema_to_load, schema, dup_cols, args.remove, args.patchall)

		if invalid_flag == True:
			print('{}: invalid schema, check the headers'.format(schema_to_load))
			summary_report.append('{}: invalid schema, check the headers'.format(schema_to_load))
			continue

		set_value_types(df, prop_types, linkTos)

		obj_posts = []

		for row_count, row in df.iterrows():
			row_count += 2 #1 for the header row, 1 to go from 0-based to 1-based
			if not args.remove:
				row = row.dropna()
			post_json = json.loads(row.to_json())

			#split out the OntologyTerms to load
			if not args.remove:
				post_json, post_ont = dict_patcher(post_json, schema['properties'])
				for k, v in post_ont.items():
					if 'ontology_term_' + schema_to_load not in load_order:
						load_order.insert(0, 'ontology_term_' + schema_to_load)
					all_posts.setdefault('ontology_term_' + schema_to_load, []).append(('', v))

				# add attchments here
				if post_json.get('attachment'):
					attach = attachment(post_json['attachment'])
					post_json['attachment'] = attach

			obj_posts.append((row_count, post_json))

		all_posts[schema_to_load] = obj_posts

	failed_postings = []
	idprops = ['uuid', 'accession', 'aliases', 'term_id', 'name']
	for schema_to_report in load_order:
		schema = 'ontology_term' if schema_to_report.startswith('ontology_term_') else schema_to_report
		if all_posts.get(schema_to_report):
			total = 0
			error = 0
			success = 0
			patch = 0
			obj_failed = []
			for row_count, post_json in all_posts[schema_to_report]:
				total += 1

				#check for an existing object based on any possible identifier
				temp_id, temp = check_existing_obj(post_json, schema, connection)
				temp_id = unquote(temp_id)

				if temp.get('uuid'): # if there is an existing corresponding object
					if args.remove:
						if schema_to_report.startswith('ontology_term_'): #the secondary ontology_terms should be ignored for remove
							total -= 1
						else:
							existing_json = lattice.get_object(temp['uuid'], connection, frame="edit")
							for k in post_json.keys():
								if k not in existing_json.keys():
									print('ERROR: {} ROW {}: Cannot remove {}, may be calculated property, or is not submitted'.format(schema_to_report.upper(), row_count, k))
								elif k not in idprops:
									existing_json.pop(k)
									print('{} ROW {}:Removing value: {}'.format(schema_to_report.upper(), row_count, k))
							if args.update:
								e = lattice.replace_object(temp['uuid'], connection, existing_json)
								if e['status'] == 'error':
									error += 1
								elif e['status'] == 'success':
									new_patched_object = e['@graph'][0]
									print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':identifier: {}'.format((new_patched_object.get(
										'accession', new_patched_object.get('uuid')))))
									patch += 1
					elif args.patchall:
						if schema_to_report.startswith('ontology_term_'):
							if post_json.get('term_name') and post_json['term_name'] != temp['term_name']:
								print('ERROR: {}: sumbission term_name {} of {} does not match existing term_name {}'.format(schema_to_report.upper(), post_json['term_name'], post_json['term_id'], temp['term_name']))
								error += 1
							else: #secondary ontology_terms that exist can be ignored if the term_name matches existing object
								total -= 1
						elif args.update:
							e = lattice.patch_object(temp['uuid'], connection, post_json)
							if e['status'] == 'error':
								error += 1
							elif e['status'] == 'success':
								new_patched_object = e['@graph'][0]
								print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':identifier: {}'.format((new_patched_object.get(
									'accession', new_patched_object.get('uuid')))))
								patch += 1
					elif schema_to_report.startswith('ontology_term_'):
						if post_json.get('term_name') and post_json['term_name'] != temp['term_name']:
							print('ERROR: {}: sumbission term_name {} of {} does not match existing term_name {}'.format(schema_to_report.upper(), post_json['term_name'], post_json['term_id'], temp['term_name']))
							error += 1
						else: #secondary ontology_terms that exist can be ignored if the term_name matches existing object
							total -= 1
					else:
						print('ERROR: {} ROW {}: Object {} already exists.'.format(schema_to_report.upper(), row_count, temp_id))
						error += 1

				else: # there is no existing object found
					if args.remove:
						print('ERROR: {} ROW {}: Object does not exist.'.format(schema_to_report.upper(), row_count))
						error += 1
					elif args.patchall:
						if schema_to_report.startswith('ontology_term_'): #possible we need to post new ontology_terms to patch existing objects
							if args.update:
								print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':POSTing data!')
								e = lattice.post_object(schema, connection, post_json)
								if e['status'] == 'error':
									error += 1
									obj_failed.append(schema_to_report.upper() + ' ROW ' + str(row_count) + ':' + 'term_id')
								elif e['status'] == 'success':
									new_object = e['@graph'][0]
									print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':New accession/UUID: {}'.format((new_object.get(
										'accession', new_object.get('uuid')))))
									success += 1
						else:
							print('ERROR: {} ROW {}: Object does not exist.'.format(schema_to_report.upper(), row_count))
							error += 1
					else: # post new object
						if args.update:
							print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':POSTing data!')
							e = lattice.post_object(schema, connection, post_json)
							if e['status'] == 'error':
								error += 1
								obj_failed.append(schema_to_report.upper() + ' ROW ' + str(row_count) + ':' + post_json.get(
									'aliases', ['alias not specified'])[0])
							elif e['status'] == 'success':
								new_object = e['@graph'][0]
								print(schema_to_report.upper() + ' ROW ' + str(row_count) + ':New accession/UUID: {}'.format((new_object.get(
									'accession', new_object.get('uuid')))))
								success += 1

			if total > 0:
				print('{sheet}: {success} posted, {patch} patched, {error} errors out of {total} total'.format(
					sheet=schema_to_report.upper(), success=success, total=total, error=error, patch=patch))
				summary_report.append('{sheet}: {success} posted, {patch} patched, {error} errors out of {total} total'.format(
					sheet=schema_to_report.upper(), success=success, total=total, error=error, patch=patch))

			if obj_failed:
				print('Posting failed for {} object(s):'.format(len(obj_failed)))
				for a in obj_failed:
					print(a)
				failed_postings.extend(obj_failed)

	if failed_postings:
		print('-------Summary of failed objects-------')
		for a in failed_postings:
			print(a)
	print('-------Summary of all objects-------')
	print('\n'.join(summary_report))


if __name__ == '__main__':
	main()
