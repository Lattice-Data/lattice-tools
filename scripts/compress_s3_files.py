import boto3
import scanpy as sc
import argparse
import sys
import os



EPILOG = '''
This script will take a list of s3 uri of uncompressed h5ad files as input, and it will create a compressed h5ad in the same S3 directory, 
with '_curated' added to file name.

Examples:

    python %(prog)s -f list_of_uri.txt

For more details:

    python %(prog)s --help
'''


TEMP_DIR = 'temp_dir'

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--file', '-f',
                        help='A text list of s3 uri of files to be compressed')
    args = parser.parse_args()
    if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit()
    return args


def download_file(download_url):
	"""
	Given an S3 uri, the file will be downloaded locally
	"""
	bucket_name = download_url.split('/')[2]
	file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
	s3client = boto3.client("s3")
	file_name = download_url.split('/')[-1]
	new_filename = file_name.replace('.h5ad', '_curated.h5ad')
	print(file_name + ' downloading')
	try:
		s3client.download_file(bucket_name, file_path, TEMP_DIR + '/' + file_name)
	except subprocess.CalledProcessError as e:
		logging.error('ERROR: Failed to find file {} on s3'.format(download_url))
		sys.exit('ERROR: Failed to find file {} on s3'.format(download_url))
	else:
		print(file_name + ' downloaded')


def main(s3_uri_file):
	if os.path.exists(TEMP_DIR) == False:
		os.mkdir(TEMP_DIR)
	print('test')



args = getArgs()
if __name__ == '__main__':
    main(args.file)




