import argparse
import lattice
import requests
import sys
from urllib.parse import urljoin, quote


EPILOG = '''
Query the entire Lattice Database for duplicate files, based on checksums.

Examples:

    python %(prog)s --mode prod

This relies on local variables to be defined based on the --mode you provide
to direct the updates to a server and to provide permissions
For example, if specifying --mode prod, to make the changes on a local instance,
the following variables need to be defined...
PROD_KEY, PROD_SECRET, PROD_SERVER

For more details:

        python %(prog)s --help
'''

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--mode', '-m',
                        help='The machine to run on.')
    args = parser.parse_args()
    return args


args = getArgs()
if not args.mode:
    sys.exit('ERROR: --mode is required')
connection = lattice.Connection(args.mode)
server = connection.server

checksums = [
    'md5sum',
    'content_md5sum',
    'crc32c',
    'sha256'
]

obj_types = [
	'RawSequenceFile',
	'RawMatrixFile'
	]

for obj_type in obj_types:
	print(obj_type)
	url = urljoin(server, 'search/?limit=all&format=json&type={}&status!=deleted'.format(obj_type))

	for c in checksums:
	    url += '&field={}'.format(c)
	results = requests.get(url, auth=connection.auth).json()

	for c in checksums:
	    full = [r.get(c) for r in results['@graph']]
	    dups = [e for e in full if e and full.count(e) >1]
	    print(c + ':' + str(len(full)) + ' checked, ' + str(len(dups)) + ' duplicated')
	    if dups:
	        for d in set(dups):
	            print(server + f'report/?type={obj_type}&{c}={d}')
	print('-----')
