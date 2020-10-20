import argparse
import boto3
import csv
import lattice
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


schema_mapping = {
	"valid_barcodes": "percent_valid_barcodes",
	"sequencing_saturation": "percent_sequencing_saturation",
	"q30_bases_in_barcode": "percent_q30_bases_in_barcode",
	"q30_bases_in_rna_read": "percent_q30_bases_in_rna_read",
	"q30_bases_in_sample_index": "percent_q30_bases_in_sample_index",
	"q30_bases_in_umi": "percent_q30_bases_in_umi",
	"reads_mapped_to_genome": "percent_reads_mapped_to_genome",
	"reads_mapped_confidently_to_genome": "percent_reads_mapped_confidently_to_genome",
	"reads_mapped_confidently_to_intergenic_regions": "percent_reads_mapped_confidently_to_intergenic_regions",
	"reads_mapped_confidently_to_intronic_regions": "percent_reads_mapped_confidently_to_intronic_regions",
	"reads_mapped_confidently_to_exonic_regions": "percent_reads_mapped_confidently_to_exonic_regions",
	"reads_mapped_confidently_to_transcriptome": "percent_reads_mapped_confidently_to_transcriptome",
	"reads_mapped_antisense_to_gene": "percent_reads_mapped_antisense_to_gene",
	"fraction_reads_in_cells": "percent_reads_in_cells",
}


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
	parser.add_argument('--mode', '-m',
						help='The machine to pull schema from.')
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
	sys.exit('ERROR: must use one of --assay rna or atac')

if not args.mode:
	sys.exit('ERROR: --mode is required')

connection = lattice.Connection(args.mode)
server = connection.server
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
	headers = [header.lower().replace(' ','_') for header in rows[0]]
	new_headers = [schema_mapping.get(header, header) for header in headers]
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
		print('{}:{} not in schema'.format(key, post_json[key]))

out_file = open('metrics.tsv', 'a')
out_file.write('\t'.join(final_headers) + '\n')
out_file.write('\t'.join(final_values) + '\n')
out_file.close()

print('metrics recorded')

