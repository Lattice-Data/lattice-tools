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
			if obj.status_code == 404 and '@graph' in obj.json():
				graph.extend(obj.json()['@graph'])
			else:
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