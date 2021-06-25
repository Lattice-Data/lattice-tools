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


def get_object(obj_id, connection, frame=None):
	obj_id = obj_id.replace(':','%3A')
	if frame is None:
		url = urljoin(connection.server, obj_id)
	else:
		url = urljoin(connection.server, obj_id + '/?frame=' + frame)
	obj = requests.get(url, auth=connection.auth).json()
	return obj


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
