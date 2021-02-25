import argparse
import boto3
import botocore
import csv
import json
import lattice
import os
import re
import requests
import subprocess
import sys
import gzip
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

system_props = [
	"status",
	"uuid",
	"schema_version",
	"date_created",
	"submitted_by",
	"notes",
	"documents",
	"@id",
	"@type"
]


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('--dir', '-d',
						help="s3 path to the cellranger outs directory")
	parser.add_argument('--assay', '-a',
						help="use if pulling metrics from rna data")
	parser.add_argument('--mode', '-m',
						help='The machine to pull schema from.')
	args = parser.parse_args()
	return args

args = getArgs()

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

if not args.mode:
	sys.exit('ERROR: --mode is required')

connection = lattice.Connection(args.mode)
server = connection.server
schema_url = urljoin(server, 'profiles/' + obj_name + '/?format=json')
full_schema = requests.get(schema_url).json()
schema_props = list(full_schema['properties'].keys())
for a in system_props:
	schema_props.remove(a)

full_path = args.dir.rstrip('/')
bucket_name = full_path.split('/')[0]
outs_dir_path = full_path.replace(bucket_name + '/', '')

s3client = boto3.client("s3")

report_json = {}

for file_name in files_to_check:
	try:
	    s3client.download_file(bucket_name, outs_dir_path + '/' + file_name, file_name)
	except botocore.exceptions.ClientError:
		print('Failed to find {} on s3'.format(file_name))
	else:
		print(file_name + ' downloaded')

		if file_name == 'web_summary.html':
			with open(file_name) as html_doc:
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

			os.remove(file_name)
			print(file_name + ' removed')

		else:
			if file_name == 'metrics_summary.csv':
				with open(file_name, newline='') as csvfile:
					spamreader = csv.reader(csvfile)
					rows = list(spamreader)
					headers = [header.lower().replace(' ','_') for header in rows[0]]
					new_headers = [schema_mapping.get(header, header) for header in headers]
					values = rows[1]
					new_values = [value.strip('%') for value in values]
					post_json = dict(zip(new_headers, new_values))
					post_json['quality_metric_of'] = '<linkTo filtered_feature_bc_matrix.h5>'
				csvfile.close()
				os.remove(file_name)
				print(file_name + ' removed')

			elif file_name == 'summary.json':
				with open(file_name) as summary_json:
					post_json = json.load(summary_json)
					my_props = list(post_json.keys())
					for prop in my_props:
						if prop in schema_mapping.keys():
							post_json[schema_mapping[prop]] = post_json[prop]
							del post_json[prop]
					post_json['quality_metric_of'] = '<linkTo filtered_peak_bc_matrix.h5>'
				os.remove(file_name)
				print(file_name + ' removed')

			report_json.update(post_json)

report_json['directory'] = full_path

extra_headers = []
extra_values = []
final_headers = ['directory'] + schema_props
final_values = []

full_schema['properties']['directory'] = {'type':'string'}

for prop in final_headers:
	if (full_schema['properties'][prop]['type'] == 'integer') and (str(report_json.get(prop, '')).endswith('.0') == True):
		final_values.append(str(report_json.get(prop, '')).strip('.0'))
	elif str(report_json.get(prop, '')) == 'None':
		final_values.append('')
	else:
		final_values.append(str(report_json.get(prop, '')))

for key in report_json.keys():
	if key not in schema_props:
		extra_headers.append(key)
		extra_values.append(str(report_json[key]))
		if key in should_match.keys():
			if report_json[key] != report_json[should_match[key]]:
				print('WARNING: {} does not match {}'.format(should_match[key], key))
			else:
				print('all good: {} does match {}'.format(should_match[key], key))

if len(extra_headers) > 1:
	none_schema_outfile = open('metrics_not_in_schema.tsv', 'a')
	none_schema_outfile.write('\t'.join(extra_headers) + '\n')
	none_schema_outfile.write('\t'.join(extra_values) + '\n')
	none_schema_outfile.close()

if len(final_headers) > 1:
	outfile = open('metrics.tsv', 'a')
	outfile.write('\t'.join(final_headers) + '\n')
	outfile.write('\t'.join(final_values) + '\n')
	outfile.close()

print('metrics recorded')
