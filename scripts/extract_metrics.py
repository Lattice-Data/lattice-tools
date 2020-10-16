import argparse
import boto3
import csv
import os
import requests
import subprocess
from urllib.parse import urljoin


EPILOG = '''
Run sanity checks on files.

Examples:

    python %(prog)s --s3-file s3://submissions-czi012eye/chen_2020/D001-12_NeuNM_outs/filtered_feature_bc_matrix.h5 --alias rui-chen:D001-12_NeuNM_matrix_filtered

For more details:

        python %(prog)s --help
'''


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--s3-file',
                        help="path to the matrix file at s3 to check")
    parser.add_argument('--alias',
                        help="an identifier to use for the file object")
    parser.add_argument('--assay',
                        help="use if pulling metrics from rna data")
    args = parser.parse_args()
    return args

args = getArgs()

if args.assay == 'atac':
	file_name = 'summary.csv'
	obj_name = 'atac_metrics'
elif args.assay == 'rna':
	file_name = 'metrics_summary.csv'
	obj_name = 'rna_metrics'
else:
	print('ERROR: must use one of --assay rna or atac')
	quit()

server = 'https://www.lattice-data.org/'
schema_url = urljoin(server, 'profiles/' + obj_name + '/?format=json')
schema_props = requests.get(schema_url).json()['properties']

matrix_uri = args.s3_file
bucket_name = matrix_uri.split('/')[2]
matrix_path = matrix_uri.replace('s3://{}/'.format(bucket_name), '')
outs_dir_path = '/'.join(matrix_path.split('/')[:-1])

s3client = boto3.client("s3")
try:
    s3client.download_file(bucket_name, outs_dir_path + '/' + file_name, file_name)
except subprocess.CalledProcessError as e:
    errors['file not found'] = 'Failed to find file on s3'
else:
    print(file_name + ' downloaded')

with open(file_name, newline='') as csvfile:
	spamreader = csv.reader(csvfile)
	rows = list(spamreader)
	headers = rows[0]
	new_headers = [header.lower().replace(' ','_') for header in headers]
	values = rows[1]
	new_values = [value.strip('%') for value in values]
	post_json = dict(zip(new_headers, new_values))
csvfile.close()
os.remove(file_name)
print(file_name + ' removed')

post_json['quality_metric_of'] = args.alias

final_headers = []
final_values = []
for key in post_json.keys():
	if key in schema_props:
		final_headers.append(key)
		final_values.append(post_json[key])
	else:
		print('{} not in schema'.format(key))

out_file = open('metrics.tsv', 'a')
out_file.write('\t'.join(final_headers) + '\n')
out_file.write('\t'.join(final_values) + '\n')
out_file.close()

print('metrics recorded')

