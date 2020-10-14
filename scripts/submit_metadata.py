import argparse
import os
import json
import lattice
import sys
import logging
import mimetypes
import requests
import magic  # install me with 'pip install python-magic'
from PIL import Image  # install me with 'pip install Pillow'
from openpyxl import load_workbook
from urllib.parse import urljoin, quote
from base64 import b64encode


EPILOG = '''
This script takes in an Excel file with the data
This is a dryrun-default script, run with --update or --patchall to work

This relies on local variables to be defined based on the --mode you provide
to direct the updates to a server and to provide permissions
For example, if specifying --mode local, to make the changes on a local instance,
the following variables need to be defined...
LOCAL_KEY, LOCAL_SECRET, LOCAL_SERVER

By DEFAULT:
If there is a uuid, alias, accession, name, or term_id (for OntologyTerms) in the document,
and an existing object is found using that identifier, it will ask if you want to PATCH that object
Use '--patchall' if you want to patch ALL objects in your document and ignore that message

DEFINING OBJECT TYPES:
	Name each sheet of the excel file the name of the object type you are using
	with the format used on https://www.lattice-data.org/profiles/,
	capitalization and underscores do not matter
	Ex: Dataset, Tissue, Document

	The objects will be loaded in the order specified in encoded/src/loadxl.py ORDER

	Use the --justtype argument to only submit one of the object types, even if your file contains more sheets
	Ex: %(prog)s mydata.xsls --justtype Experiment

	Use the --starttype argument to start at an object type and only submit the sequential objects

SPREADSHEET FORMATTING:
	The top row of each sheet should be the names of the fields specified in the schema.
	The validity of these fields will be check during a 'dry-run' (i.e. if '--update' is not used)

	If the first value in a row (column A) begins with '#', the entire row will be ignored
	This is useful if you want rows with property descriptors (e.g. description, enum, linkTo)
	These rows to ignore must be above metadata-to-submit rows or the Row number print-outs will be inaccurate

	If the first cell (A1) is either empty or begins with 'schema_version=', the entire first column (A) will be ignored
	This action happens after the above avoidance of rows starting with '#'

	Any sheet with the name 'Cover Sheet' will be ignored, capitalization does not matter, but the space does

To upload objects with attachments, have a column titled 'attachment'
containing the name of the file you wish to attach

FOR EMBEDDED SUBOBJECTS:
	Embedded objects are considered to be things like:
	 - flowcell_details in RawSequenceFile
	 - filtering_cutoffs in MatrixFile
	 They are assumed to be of the format of dictionary objects or a list of dictionary objects

	 If you are submitting just one dictionary object...
	 Ex:
	 'filtering_cutoffs': [
			{
				'cutoff_value': 30,
				'cutoff_units': 'genes/cell',
				'cutoff_type': 'minimum'
			}
		]
	Formatting in the document should be as follows for the above example:
	filtering_cutoffs.cutoff_value	filtering_cutoffs.cutoff_units
	30								genes/cell

	If you are submitting a list of multiple dictionary objects...
	Ex:
	'filtering_cutoffs': [
			{
				'cutoff_value': 30,
				'cutoff_units': 'genes/cell',
				'cutoff_type': 'minimum'
			},
			{
				'cutoff_value': 25,
				'cutoff_units': 'cells/gene',
				'cutoff_type': 'minimum'
			}
		]

	An identifier (number or letter) should be appended to the property names w/ '-' in order to group them appropriately
	filtering_cutoffs.cutoff_value-1	filtering_cutoffs.cutoff_units-1	filtering_cutoffs.cutoff_type-1	filtering_cutoffs.cutoff_value-2	filtering_cutoffs.cutoff_units-2	filtering_cutoffs.cutoff_type-2
	30									genes/cell							minimum							25					  				cells/genes							minimum

For more details:

		%(prog)s --help
'''


ORDER = [
    'user',
    'award',
    'lab',
    'organism',
    'gene',
    'publication',
    'document',
    'library_protocol',
    'antibody_lot',
    'ontology_term',
    'treatment',
    'human_postnatal_donor',
    'human_prenatal_donor',
    'mouse_postnatal_donor',
    'mouse_prenatal_donor',
    'tissue',
    'cell_culture',
    'organoid',
    'suspension',
    'dataset',
    'reference_file_set',
    'library',
    'reference_file',
    'sequencing_run',
    'raw_sequence_file',
    'sequence_alignment_file',
    'matrix_file',
    'rna_metrics',
    'antibody_capture_metrics',
    'rna_aggregate_metrics',
    'atac_metrics',
    'atac_aggregate_metrics',
    'image',
    'page',
    'access_key'
]


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	parser.add_argument('infile',
						help='the datafile containing object data to import')
	parser.add_argument('--justtype', '-jt',
						help='the type of the objects to import if only one type is desired')
	parser.add_argument('--starttype', '-st',
						help='the type of the objects to start importing and continue in sequence')
	parser.add_argument('--mode', '-m',
						help='The machine to run on.')
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
	args = parser.parse_args()
	return args


def reader(book, sheetname):
	sheet = book[sheetname]
	row_count = 1 # start with 1 for the headers
	rows = []
	for row in sheet.iter_rows():
		row = ['' if v.value is None else cell_value(v) for v in row]
		if len(row) != row.count(''):
			if row[0] != '' and row[0][0] == '#': # avoid rows starting with comments
				row_count += 1 # but add one row for each, assuming they are above submission rows
			else:
				rows.append(list(row))
	cell_A1 = rows[0][0]
	if cell_A1 is None or cell_A1.startswith('schema_version='): # may need to delete values in the first column
		for row in rows:
			del row[0]
	# now remove the values in empty columns
	bad_elements = []
	for i in range(0,len(rows[0])):
		temp = []
		for row in rows:
			temp.append(row[i])
		if len(temp) == temp.count(''):
			bad_elements.append(i)
	for ele in sorted(bad_elements, reverse = True):
		for row in rows:
			del row[ele]
	return(row_count, rows)


def cell_value(cell):
	ctype = cell.data_type
	cvalue = cell.value

	if ctype == 'e': # TYPE_ERROR
		raise ValueError(repr(cell), 'cell error')

	elif ctype == 'b': # TYPE_BOOL
		return str(cvalue).upper()

	elif ctype == 'n': # TYPE_NUMERIC
		return int(cvalue)

	elif ctype == 's': # TYPE_STRING
		return cvalue

	elif ctype == 'd': # datetime
		return cvalue.date()

	elif ctype == 'f': # TYPE_FORMULA
		raise ValueError(repr(cell), 'formula in cell')

	else:
		raise ValueError(repr(cell), '-'.join(['unknown cell type',ctype,cvalue]))


def properties_validator(keys, schema, schema_properties):
	flag = False
	dup_keys = []
	for key in keys:
		if key is not None:
			if keys.count(key) >1 and key not in dup_keys: # check for duplicated headers
				print('Property {prop} found {count} times in {schema} sheet headers'.format(prop=key, count=keys.count(key), schema=schema))
				dup_keys.append(key)
				flag = True
			key_no_iter = key.split('-')[0] # remove the '-n' notation to check that each header is a property in the specified schema
			if '.' in key_no_iter:
				if key_no_iter.split('.')[0] not in schema_properties.keys():
					print('Property {prop} not found in schema {schema}'.format(prop=key_no_iter.split('.')[0], schema=schema))
					flag = True
				elif schema_properties[key_no_iter.split('.')[0]]['type'] == 'array':
					if schema_properties[key_no_iter.split('.')[0]]['items'].get('additionalProperties','false') == 'false':
						if key_no_iter.split('.')[1] not in schema_properties[key_no_iter.split('.')[0]]['items']['properties'].keys():
							print('Property {prop} not found in schema {schema}'.format(prop=key_no_iter, schema=schema))
							flag = True
				else:
					if schema_properties[key_no_iter.split('.')[0]].get('additionalProperties','false') == 'false':
						if key_no_iter.split('.')[1] not in schema_properties[key_no_iter.split('.')[0]]['properties'].keys():
							print('Property {prop} not found in schema {schema}'.format(prop=key_no_iter, schema=schema))
							flag = True
			elif key_no_iter not in schema_properties.keys():
				print('Property {prop} not found in schema {schema}'.format(prop=key_no_iter, schema=schema))
				flag = True
	return flag


def type_formatter(old_value, schema_properties, key1, key2=None):
	# determine the type specified in the schema
	linkTo_flag = False
	desired_type = ''
	if not key2:
		desired_type = schema_properties[key1]['type']
		if schema_properties[key1].get('linkTo'):
			linkTo_flag = True
		if schema_properties[key1].get('items'):
			if schema_properties[key1]['items'].get('linkTo'):
				linkTo_flag = True
	else:
		if schema_properties[key1].get('items'):
			if schema_properties[key1]['items']['properties'].get(key2):
				desired_type = schema_properties[key1]['items']['properties'][key2]['type']
				if schema_properties[key1]['items']['properties'][key2].get('linkTo'):
					linkTo_flag = True
		else:
			if schema_properties[key1]['properties'].get(key2):
				desired_type = schema_properties[key1]['properties'][key2]['type']
				if schema_properties[key1]['items']['properties'][key2].get('linkTo'):
					linkTo_flag = True
	# adjust the value to the specified type
	if desired_type == 'array' and linkTo_flag == True:
		return [quote(x.strip()) for x in old_value.split(',')]
	if desired_type == 'array':
		return [x.strip() for x in old_value.split(',')]
	elif desired_type == 'boolean':
		if str(old_value).lower() in ['true', '1']:
			return True
		elif str(old_value).lower() in ['false', '0']:
			return False
		else:
			print('Boolean was expected for {prop} but value is {value}'.format(prop=key1, value=old_value))
	elif desired_type == 'number':
		return float(old_value)
	elif desired_type == 'integer':
		return int(old_value)
	elif desired_type == 'string' and linkTo_flag == True:
		return quote(str(old_value))
	elif desired_type == 'string':
		return str(old_value)
	else:
		return old_value


def dict_patcher(old_dict, schema_properties):
	new_dict = {}
	array_o_objs_dict = {}
	for key in sorted(old_dict.keys()):
		if old_dict[key] != '':  # this removes empty cells
			path = key.split('.')
			if len(path) == 1:
				new_dict[key] = type_formatter(old_dict[key], schema_properties, key)
			elif len(path) == 2: # embedded object, need to build the mini dictionary to put this in
				if '-' in key: # this has a number next to it and we expect an array of objects
					value = path[1].split('-')
					if array_o_objs_dict.get(path[0]): # I have already added the embedded object to the new dictionary
						if value[1] not in array_o_objs_dict[path[0]].keys():
							# this means we have not added any part of new item to the list
							array_o_objs_dict[path[0]][value[1]] = {}
							array_o_objs_dict[path[0]][value[1]] = {value[0]: type_formatter(old_dict[key], schema_properties, path[0], value[0])}
						else:
							# this should be that we have started putting in the new object
							array_o_objs_dict[path[0]][value[1]].update({value[0]: type_formatter(old_dict[key], schema_properties, path[0], value[0])})
					else: # I need to make new item in dictionary
						array_o_objs_dict[path[0]] = {}
						array_o_objs_dict[path[0]][value[1]] = {value[0]: type_formatter(old_dict[key], schema_properties, path[0], value[0])}
				else: # there is no number next to it, this is our only object to handle
					if new_dict.get(path[0]): # I have already added the embedded object to the new dictionary
						new_dict[path[0]].update({path[1]: type_formatter(old_dict[key], schema_properties, path[0], path[1])})
					else: # I need to make new item in dictionary
						temp_dict = {path[1]: type_formatter(old_dict[key], schema_properties, path[0], path[1])}
						new_dict[path[0]] = temp_dict
			elif len(path) > 2:
				print('ERROR:' + key + ': Not prepared to do double-embedded properties')
	new_list = []
	for prop in array_o_objs_dict.keys():
		for group_id in array_o_objs_dict[prop].keys():
			new_list.append(array_o_objs_dict[prop][group_id])
		new_dict[prop] = new_list
	return new_dict


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
	detected_type = magic.from_file(path, mime=True).decode('ascii')

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


def main():
	summary_report = []
	args = getArgs()
	if not args.mode:
		sys.exit('ERROR: --mode is required')
	connection = lattice.Connection(args.mode)
	server = connection.server
	print('Running on {server}'.format(server=server))
	if not os.path.isfile(args.infile):
		sys.exit('ERROR: file {filename} not found!'.format(filename=args.infile))
	book = load_workbook(args.infile)
	names = {}
	if args.justtype:
		if args.starttype:
			sys.exit('ERROR: cannot specify both --justtype and --starttype')
		else:
			for sheet in book.sheetnames:
				if sheet.lower().replace('_','') == args.justtype.lower().replace('_',''):
					names[sheet.lower().replace('_','')] = sheet
	else:
		for sheet in book.sheetnames:
			names[sheet.lower().replace('_','')] = sheet
	profiles = requests.get(server + 'profiles/?format=json').json()
	supported_collections = [s.lower() for s in list(profiles.keys())] # get accepted object types
	supported_collections.append('cover sheet')
	for n in names.keys():
		if n not in supported_collections: # check that each sheet name corresponds to an object type
			print('ERROR: Sheet name {name} not part of supported object types!'.format(
				name=n), file=sys.stderr)
	load_order = ORDER # pull in the order used to load test inserts on a local instance
	if args.starttype:
		st_index = load_order.index(args.starttype)
		load_order = load_order[st_index:]
	for schema_to_load in load_order: # go in order to try and get objects posted before they are referenced by another object
		obj_type = schema_to_load.replace('_','')
		if obj_type in names.keys():
			row_count, rows = reader(book, names[obj_type])
			headers = rows.pop(0)
			schema_url = urljoin(server, 'profiles/' + schema_to_load + '/?format=json')
			schema_properties = requests.get(schema_url).json()['properties']
			invalid_flag = properties_validator(headers, schema_to_load, schema_properties)
			if invalid_flag == True:
				print('{sheet}: invalid schema, check the headers'.format(sheet=schema_to_load))
				summary_report.append('{sheet}: invalid schema, check the headers'.format(sheet=schema_to_load))
				continue
			total = 0
			error = 0
			success = 0
			patch = 0
			new_accessions_aliases = []
			failed_postings = []
			for row in rows:
				row_count += 1
				total += 1
				post_json = dict(zip(headers, row))
				# convert values to the type specified in the schema, including embedded json objects
				post_json = dict_patcher(post_json,schema_properties)
				# add attchments here
				if post_json.get('attachment'):
					attach = attachment(post_json['attachment'])
					post_json['attachment'] = attach
				#check for an existing object based on any possible identifier
				if post_json.get('uuid'):
					url = urljoin(server, post_json['uuid'] + '/?format=json')
				elif post_json.get('aliases'):
					url = urljoin(server, quote(post_json['aliases'][0]) + '/?format=json')
				elif post_json.get('accession'):
					url = urljoin(server, post_json['accession'] + '/?format=json')
				elif post_json.get('name'):
					url = urljoin(server, schema_to_load + '/' + post_json['name'] + '/?format=json')
				elif post_json.get('term_id'):
					url = urljoin(server, schema_to_load + '/' + post_json['term_id'].replace(':','_') + '/?format=json')
				patch_flag = False

				try:
					url
				except NameError:
					temp = {}
				else:
					temp = requests.get(url, auth=connection.auth).json()

				if temp.get('uuid'): # if there is an existing corresponding object
					if args.patchall:
						patch_flag = True
					else: # patch wasn't specified, see if the user wants to patch
						print(schema_to_load.upper() + ' ROW ' + str(row_count) + ':Object {} already exists.  Would you like to patch it instead?'.format(
							temp['uuid']))
						i = input('PATCH? y/n ')
						if i.lower() == 'y':
							patch_flag = True
					if patch_flag == True and args.update:
						e = lattice.patch_object(temp['uuid'], connection, post_json)
						if e['status'] == 'error':
							error += 1
						elif e['status'] == 'success':
							new_patched_object = e['@graph'][0]
							# Print now and later
							print(schema_to_load.upper() + ' ROW ' + str(row_count) + ':identifier: {}'.format((new_patched_object.get(
								'accession', new_patched_object.get('uuid')))))
							patch += 1
				else: # we have new object to post
					if args.patchall:
						print(schema_to_load.upper() + ' ROW ' + str(row_count) + ':Object not found. Check identifier or consider removing --patchall to post a new object')
						error += 1
					elif args.update:
						print(schema_to_load.upper() + ' ROW ' + str(row_count) + ':POSTing data!')
						e = lattice.post_object(schema_to_load, connection, post_json)
						if e['status'] == 'error':
							error += 1
							failed_postings.append(schema_to_load.upper() + ' ROW ' + str(row_count) + ':' + str(post_json.get(
								'aliases', 'alias not specified')))
						elif e['status'] == 'success':
							new_object = e['@graph'][0]
							# Print now and later
							print(schema_to_load.upper() + ' ROW ' + str(row_count) + ':New accession/UUID: {}'.format((new_object.get(
								'accession', new_object.get('uuid')))))
							new_accessions_aliases.append(('ROW ' + str(row_count), new_object.get(
								'accession', new_object.get('uuid')), new_object.get('aliases', new_object.get('name'))))
							success += 1
			# Print now and later
			print('{sheet}: {success} posted, {patch} patched, {error} errors out of {total} total'.format(
				sheet=schema_to_load.upper(), success=success, total=total, error=error, patch=patch))
			summary_report.append('{sheet}: {success} posted, {patch} patched, {error} errors out of {total} total'.format(
				sheet=schema_to_load.upper(), success=success, total=total, error=error, patch=patch))
			if new_accessions_aliases:
				print('New accessions/UUIDs and aliases:')
				for (row, accession, alias) in new_accessions_aliases:
					if alias == None:
						alias = 'alias not specified'
					else:
						alias = ', '.join(alias) if isinstance(alias, list) else alias
					print(row, accession, alias)
			if failed_postings:
				print('Posting failed for {} object(s):'.format(len(failed_postings)))
				for alias in failed_postings:
					print(', '.join(alias) if isinstance(alias, list) else alias)
	print('-------Summary of all objects-------')
	print('\n'.join(summary_report))


if __name__ == '__main__':
	main()
