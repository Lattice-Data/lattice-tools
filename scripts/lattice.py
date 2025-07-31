import json
import logging
import os
import requests
import sys
from urllib.parse import urljoin


class Connection(object):
    def __init__(self, mode):
        if not (os.environ.get(mode.upper() + '_KEY') 
            and os.environ.get(mode.upper() + '_SECRET') 
            and os.environ.get(mode.upper() + '_SERVER')):
            sys.exit('ERROR: ' + mode.upper() + '_KEY, ' + mode.upper() + '_SECRET, ' + mode.upper() + "_SERVER not all defined. Try 'conda env config vars list' to list existing variables")
        self.authid = os.environ.get(mode.upper() + '_KEY')
        self.authpw = os.environ.get(mode.upper() + '_SECRET')
        self.server = os.environ.get(mode.upper() + '_SERVER')
        if not self.server.endswith('/'):
            self.server += '/'
        self.headers = {'content-type': 'application/json',
                        'accept': 'application/json'}
        self.auth = (self.authid, self.authpw)


def parse_ids(ids_lst):
	"""
	Takes a list of @ids and returns the object type and filters that designate all the ids as a complete string
	"""
	at_id = ids_lst[0].split('/')[1]
	if at_id == 'libraries':
		obj_type = 'Library'
	else:
		obj_type = ''.join([x.capitalize() for x in at_id.split('-')]).rstrip('s')
	filter_url = ''.join(["&@id="+i for i in ids_lst])
	return obj_type, filter_url


def get_report(obj_type, filter_url, field_lst, connection):
	"""
	Constructs a report url of fields in field_list for the objects determined by the filter and returns list of dictionaries from @graph
	Will split url into two requests if url is > 8000 characters, and still return a single list of dictionaries from @graph
	"""
	field_url = ''.join(["&field="+i for i in field_lst])
	url1 = urljoin(connection.server, "report/?type={}{}{}&format=json&limit=all".format(obj_type, filter_url, field_url))
	urls = []
	if len(url1) > 8000:
		filter_lst = filter_url.split('&@id=')
		filter_url1 = ''.join(["&@id="+i for i in filter_lst[:len(filter_lst)//2]])
		filter_url2 = ''.join(["&@id="+i for i in filter_lst[len(filter_lst)//2:]])
		url1 = urljoin(connection.server, "report/?type={}{}{}&format=json&limit=all".format(obj_type, filter_url1, field_url))
		url2 = urljoin(connection.server, "report/?type={}{}{}&format=json&limit=all".format(obj_type, filter_url2, field_url))
		urls.append(url2)
	urls.append(url1)
	graph = []
	for url in urls:
		try:
			obj = requests.get(url, auth=connection.auth)
			obj.raise_for_status()
		except requests.exceptions.HTTPError as err:
			print("HTTP Error: ", err)
			sys.exit()
		except requests.exceptions.Timeout as err:
			print ("Timeout Error: ",err)
			sys.exit()
		except requests.exceptions.RequestException as err:
			print ("Requests error: ",err)
			sys.exit()
		else:
			graph.extend(obj.json().get('@graph'))	
	return graph


def get_object(obj_id, connection, frame=None):
	"""
	Constructs a url for a single object and returns the object with option of having frame=object
	"""
	obj_id = obj_id.replace(':','%3A')
	if frame is None:
		url = urljoin(connection.server, obj_id)
	else:
		url = urljoin(connection.server, obj_id + '/?frame=' + frame)
	try:
		obj = requests.get(url, auth=connection.auth)
		obj.raise_for_status()
	except requests.exceptions.HTTPError as err:
		print("HTTP Error: ", err)
		sys.exit()
	except requests.exceptions.Timeout as err:
		print ("Timeout Error: ",err)
		sys.exit()
	except requests.exceptions.RequestException as err:
		print ("Requests error: ",err)
		sys.exit()
	else:
		return obj.json()


def post_object(schema, connection, post_json):
	if isinstance(post_json, dict):
		json_payload = json.dumps(post_json)
	elif isinstance(post_json, str):
		json_payload = post_json
	else:
		print('Datatype to POST is not string or dict.', file=sys.stderr)
	url = urljoin(connection.server, schema)
	logging.debug('POST URL : %s' % (url))
	logging.debug('POST data: %s' % (json.dumps(post_json,
												sort_keys=True, indent=4,
												separators=(',', ': '))))
	response = requests.post(url, auth=connection.auth,
							 headers=connection.headers, data=json_payload)
	logging.debug('POST RESPONSE: %s' % (json.dumps(response.json(),
													indent=4, separators=(',', ': '))))
	if not response.status_code == 201:
		logging.warning('POST failure. Response = %s' % (response.text))
	logging.debug('Return object: %s' % (json.dumps(response.json(),
													sort_keys=True, indent=4,
													separators=(',', ': '))))
	return response.json()


def patch_object(obj_id, connection, patch_input):
	if isinstance(patch_input, dict):
		json_payload = json.dumps(patch_input)
	elif isinstance(patch_input, str):
		json_payload = patch_input
	else:
		print('Datatype to PATCH is not string or dict.', file=sys.stderr)
	url = urljoin(connection.server, obj_id)
	logging.debug('PATCH URL : %s' % (url))
	logging.debug('PATCH data: %s' % (json_payload))
	response = requests.patch(url, auth=connection.auth, data=json_payload,
							  headers=connection.headers)
	logging.debug('PATCH RESPONSE: %s' % (json.dumps(response.json(), indent=4,
													 separators=(',', ': '))))
	if not response.status_code == 200:
		logging.warning('PATCH failure.  Response = %s' % (response.text))
	return response.json()


def replace_object(obj_id, connection, post_json):
    if isinstance(post_json, dict):
        json_payload = json.dumps(post_json)
    elif isinstance(post_json, str):
        json_payload = post_json
    else:
        logging.warning('Datatype to PUT is not string or dict.')
    url = urljoin(connection.server, obj_id)
    logging.debug('PUT URL : %s' % (url))
    logging.debug('PUT data: %s' % (json_payload))
    response = requests.put(url, auth=connection.auth, data=json_payload,
                            headers=connection.headers)
    logging.debug('PUT RESPONSE: %s' % (json.dumps(response.json(), indent=4,
                                                   separators=(',', ': '))))
    if not response.status_code == 200:
        logging.warning('PUT failure.  Response = %s' % (response.text))
    return response.json()


def create_subdirectory(subdirectory):
    if os.path.exists('./outputs') == False:
        os.mkdir('./outputs')
    if os.path.exists('./outputs/' + subdirectory) == False:
        os.mkdir('./outputs/' + subdirectory)
    full_filepath = './outputs/' + subdirectory
    return full_filepath


def check_audit(obj):
	obj_type = obj['@type'][0]
	if 'accession' in obj:
		i = obj['accession']
	else:
		i = obj.get('uuid')
    # if there are ERROR-level audits, flag it to stop the script
	if (obj.get('audit') and obj['audit'].get('ERROR')):
		freq = {}
		for a in obj['audit']['ERROR']:
			if a['category'] in freq:
				freq[a['category']] += 1
			else:
				freq[a['category']] = 1
		for k,v in freq.items():
			print('ERROR audit:{}x {} on {} {}'.format(
				str(v),
				k,
				obj_type,
				i
			))
		i = input('Continue? y/n: ')
		if i.lower() not in ['y','yes']:
			sys.exit('Stopped due to one or more ERROR audits')
