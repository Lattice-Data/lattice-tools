import argparse
import boto3
import botocore
import csv
import json
import lattice
import os
import pandas as pd
import re
import requests
import subprocess
import sys
from urllib.parse import urljoin
from bs4 import BeautifulSoup


EPILOG = '''
Extract summary info and QC metrics from a cellranger pipeline run.

Examples:

    python %(prog)s -m production -a atac -d submissions-czi009kid/muto_humphreys_2020/Control_5/outs
    python %(prog)s -m local -a rna -d submissions-czi012eye/chen_2020/19D013_foveaR_outs

For more details:

        python %(prog)s --help
'''

schema_mapping = {
	"valid_barcodes": "percent_valid_barcodes",
	"sequencing_saturation": "percent_sequencing_saturation",
	"q30_bases_in_barcode": "percent_q30_bases_in_barcode",
	"q30_bases_in_rna_read": "percent_q30_bases_in_rna_read",
	"q30_bases_in_rna_read_2": "percent_q30_bases_in_rna_read_2",
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
	"frac_valid_barcode": "frac_valid_barcodes",
	"frac_waste_duplicate": "frac_waste_dup",
	"estimated_fraction_cells_annotated": "estimated_frac_cells_annotated",
	"fraction_cell_calling_noise": "frac_cell_calling_noise",
	"fraction_gelbead_doublets_cells": "frac_gelbead_doublets_cells",
	"fraction_of_genome_within_2000bp_of_peaks": "frac_of_genome_within_2000bp_of_peaks",
	"bc_q30_bases_fract": "frac_q30_bases_in_barcode",
	"r1_q30_bases_fract": "frac_q30_bases_in_read1",
	"r2_q30_bases_fract": "frac_q30_bases_in_read2",
	"si_q30_bases_fract": "frac_q30_bases_in_sample_index",
	"median_frags_overlapping_peaks_per_cell": "median_fragments_overlapping_peaks_per_cell",
	"frac_waste_mitochondrial": "frac_waste_mito",
	"waste_non_cell_barcode_fragments": "waste_noncell_barcode_fragments",
	"frac_waste_non_cell_barcode": "frac_waste_noncell_barcode",
	"frac_valid_noncell": "frac_valid_noncells",
	"num_fragments": "number_fragments",
	"num_reads": "number_reads",
	"r1_tot_bases": "total_bases_in_read1",
	"r2_tot_bases": "total_bases_in_read2",
	"si_tot_bases": "total_bases_in_sample_index",
	"bc_tot_bases": "total_bases_in_barcode",
	"si_q30_bases": "q30_bases_in_sample_index",
	"r2_q30_bases": "q30_bases_in_read2",
	"r1_q30_bases": "q30_bases_in_read1",
	"bc_q30_bases": "q30_bases_in_barcode",
	"waste_duplicate_fragments": "waste_dup_fragments",
	"waste_mitochondrial_fragments": "waste_mito_fragments"
}

value_mapping = {
	"Single Cell 3' v1": "SC3Pv1",
	"Single Cell 3' v2": "SC3Pv2",
	"Single Cell 3' v3": "SC3Pv3",
	"Single Cell 5' PE": "SC5P-PE",
	"Single Cell 5' R2-only": "SC5P-R2",
	"Single Cell 5' R1-only": "SC5P-R1"
}

should_match = {
	"annotated_cells": "cells_detected",
	"waste_cell_dup_fragments": "waste_dup_fragments",
	"waste_cell_lowmapq_fragments": "waste_lowmapq_fragments",
	"waste_cell_mito_fragments": "waste_mito_fragments",
	"waste_cell_unmapped_fragments": "waste_unmapped_fragments"
}

def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('--dir', '-d',
						help="s3 path to the cellranger outs directory or local path to a file that lists directories")
	parser.add_argument('--assay', '-a',
						help="specify atac or rna")
	parser.add_argument('--mode', '-m',
						help='The machine to pull schema from.')
	args = parser.parse_args()
	return args

args = getArgs()

if not args.mode:
	sys.exit('ERROR: --mode is required')

s3client = boto3.client("s3")

connection = lattice.Connection(args.mode)
server = connection.server

if args.assay == 'rna':
	obj_name = 'rna_metrics'
	files_to_check = [
		'metrics_summary.csv',
		'web_summary.html'
	]

elif args.assay == 'atac':
	obj_name = 'atac_metrics'
	files_to_check = [
		'summary.json',
		'web_summary.html'
	]

else:
	sys.exit('must specify rna or atac for --assay')

schema_url = urljoin(server, 'profiles/{}/?format=json'.format(obj_name))
full_schema = requests.get(schema_url).json()
schema_props = list(full_schema['properties'].keys())
schema_version = 'schema_version=' + (full_schema['properties']['schema_version']['default'])

dir_list = args.dir

if os.path.isfile(dir_list):
    directories = [line.rstrip('\n') for line in open(dir_list)]
else:
    directories = dir_list.split(',')

in_schema = {}
out_schema = {}

for direct in directories:
	direct = direct.replace('s3://', '')
	full_path = direct.rstrip('/')
	bucket_name = full_path.split('/')[0]
	outs_dir_path = full_path.replace(bucket_name + '/', '')

	report_json = {'quality_metric_of': '<linkTo RawMatrixFile - filtered matrix .h5>'}

	summary_file = 'web_summary.html'
	try:
	    s3client.download_file(bucket_name, outs_dir_path + '/' + summary_file, summary_file)
	except botocore.exceptions.ClientError:
		print('Failed to find {} on s3'.format(summary_file))
	else:
		print(summary_file + ' downloaded')
		with open(summary_file) as html_doc:
			match_flag = False
			soup = BeautifulSoup(html_doc, 'html.parser')
			for x in soup.find_all('script'):
				match = re.search("const data = ", x.string)
				if match:
					match_flag = True
					end = match.end()
					data = json.loads(x.string[end:])
					if data.get('pipeline_info_table'):
						pipeline_info_table = data.get('pipeline_info_table')
					else:
						pipeline_info_table = data['summary']['summary_tab']['pipeline_info_table']
					info_list = pipeline_info_table['rows']
					for pair in info_list:
						report_json[pair[0]] = value_mapping.get(pair[1], pair[1])
			if match_flag == False:
				for x in soup.find_all('table', id='sample_table'):
					for row in x.find_all('tr'):
						col_count = 1
						columns = row.find_all('td')
						for column in columns:
							if col_count == 1:
								field = column.get_text().strip()
							else:
								value = column.get_text().strip()
							col_count += 1
						report_json[field] = value_mapping.get(value, value)

		os.remove(summary_file)
		print(summary_file + ' removed')

	if args.assay == 'atac':
		metrics_file = 'summary.json'
		try:
		    s3client.download_file(bucket_name, outs_dir_path + '/' + metrics_file, metrics_file)
		except botocore.exceptions.ClientError:
			print('Failed to find {} on s3'.format(metrics_file))
		else:
			with open(metrics_file) as summary_json:
				post_json = json.load(summary_json)
				my_props = list(post_json.keys())
				for prop in my_props:
					if prop in schema_mapping.keys():
						post_json[schema_mapping[prop]] = post_json[prop]
						del post_json[prop]
	else:
		metrics_file = 'metrics_summary.csv'
		try:
		    s3client.download_file(bucket_name, outs_dir_path + '/' + metrics_file, metrics_file)
		except botocore.exceptions.ClientError:
			print('Failed to find {} on s3'.format(metrics_file))
		else:
			with open(metrics_file, newline='') as csvfile:
				spamreader = csv.reader(csvfile)
				rows = list(spamreader)
				headers = [header.lower().replace(' ','_') for header in rows[0]]
				new_headers = [schema_mapping.get(header, header) for header in headers]
				values = rows[1]
				new_values = [value.strip('%') for value in values]
				post_json = dict(zip(new_headers, new_values))

	os.remove(metrics_file)
	print(metrics_file + ' removed')

	report_json.update(post_json)

	final_values = {}
	extra_values = {}

	for prop, value in report_json.items():
		if prop in schema_props:
			if (full_schema['properties'][prop]['type'] == 'integer') and (str(value).endswith('.0') == True):
				value = str(value).strip('.0')
			elif value == None:
				value = ''
			else:
				value = str(value)
			final_values[prop] = value
		else:
			extra_values[prop] = value
			if prop in should_match.keys():
				if report_json[prop] != report_json[should_match[prop]]:
					print('WARNING: {} does not match {}'.format(should_match[prop], prop))
				else:
					print('all good: {} does match {}'.format(should_match[prop], prop))
	in_schema[direct] = final_values
	out_schema[direct] = extra_values

df = pd.DataFrame(in_schema).transpose()
df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
df.index.name = schema_version
df.to_csv(args.assay + '_metrics.tsv', sep='\t')

df = pd.DataFrame(out_schema).transpose()
df.index.name = schema_version
df.to_csv(args.assay + '_metrics_not_in_schema.tsv', sep='\t')
