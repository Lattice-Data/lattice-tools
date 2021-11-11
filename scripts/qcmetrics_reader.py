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
import qcmetrics_mapper
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


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('--dir', '-d',
						help="s3 path to the cellranger outs directory or dragen html, or local path to a file that lists those")
	parser.add_argument('--assay', '-a',
						help="specify atac or rna")
	parser.add_argument('--mode', '-m',
						help='The machine to pull schema from.')
	parser.add_argument('--pipeline', '-p',
						help='The pipeline that generated the metrics.')
	args = parser.parse_args()
	return args


def schemify(value, prop_type):
	if (prop_type == 'integer') and (str(value).endswith('.0') == True):
		return str(value).strip('.0')
	elif value == None:
		return ''
	else:
		return str(value)

args = getArgs()

if not args.mode:
	sys.exit('ERROR: --mode is required')
if not args.pipeline:
	sys.exit('ERROR: --pipeline is required')

s3client = boto3.client("s3")

connection = lattice.Connection(args.mode)
server = connection.server

dir_list = args.dir

if os.path.isfile(dir_list):
    directories = [line.rstrip('\n') for line in open(dir_list)]
else:
    directories = dir_list.split(',')

in_schema = {}
out_schema = {}
genotypemetrics = []

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

if args.pipeline.lower() in ['cr','cellranger']:
	value_mapping = qcmetrics_mapper.cellranger['value_mapping']
	schema_mapping = qcmetrics_mapper.cellranger['schema_mapping']
	should_match = qcmetrics_mapper.cellranger['should_match']

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
									if not value:
										value = ''
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
				final_values[prop] = schemify(value, full_schema['properties'][prop]['type'])
			else:
				extra_values[prop] = value
				# CHANGE TO for k,v in should_match.items():
				if prop in should_match.keys():
					if report_json[prop] != report_json[should_match[prop]]:
						print('WARNING: {} does not match {}'.format(should_match[prop], prop))
					else:
						print('all good: {} does match {}'.format(should_match[prop], prop))
		in_schema[direct] = final_values
		out_schema[direct] = extra_values

elif args.pipeline.lower() == 'dragen':
	schema_mapping = qcmetrics_mapper.dragen['schema_mapping']
	value_mapping = qcmetrics_mapper.dragen['value_mapping']
	for file in directories:
		file = file.replace('s3://', '')
		full_path = file.rstrip('/')
		bucket_name = full_path.split('/')[0]
		outs_dir_path = '/'.join(full_path.split('/')[1:-1])

		report_json = {'quality_metric_of': '<linkTo RawMatrixFile - filtered matrix .h5>'}

		summary_file = full_path.split('/')[-1]
		try:
		    s3client.download_file(bucket_name, outs_dir_path + '/' + summary_file, summary_file)
		except botocore.exceptions.ClientError:
			print('Failed to find {} on s3'.format(summary_file))
		else:
			print(summary_file + ' downloaded')
			with open(summary_file) as html_doc:
				soup = BeautifulSoup(html_doc, 'html.parser')

				sc_metrics = soup.find('main', id='scrnaseq-metrics-page')
				x = sc_metrics.find('table', {"class":"table table-striped table-hover"})
				for row in x.find_all('tr'):
					col_count = 1
					columns = row.find_all('td')
					for column in columns:
						if col_count == 1:
							field = column.get_text().strip()
							field = schema_mapping.get(field, field)
						else:
							value = column.get_text().strip()
							if field in value_mapping:
								factor = value_mapping[field]['factor']
								action = value_mapping[field]['action']
								if action == 'multiply':
									value = str(float(value) * factor)
							report_json[field] = value
						col_count += 1

				map_metrics = soup.find('main', id='mapping-metrics-page')
				x = map_metrics.find('table', {"class":"table table-striped"})
				for row in x.find_all('tr'):
					col_count = 1
					columns = row.find_all('td')
					for column in columns:
						if col_count == 1:
							field = column.get_text().strip()
							field = schema_mapping.get(field, field)
						elif col_count == 2:
							value = column.get_text().strip()
							report_json[field] = value
						else:
							per = column.get_text().strip()
							if per:
								field = field + ' %'
								field = schema_mapping.get(field, field)
								report_json[field] = per
						col_count += 1

				genotype_metrics = soup.find( 'main', id='genotype-demultiplexing-page')
				tables = genotype_metrics.find_all('table', {"class":"table table-striped table-hover"})
				x = tables[1]
				for row in x.find_all('tr'):
					col_count = 1
					columns = row.find_all('td')
					for column in columns:
						if col_count == 1:
							field = column.get_text().strip()
							field = schema_mapping.get(field, field)
						else:
							value = column.get_text().strip()
							report_json[field] = value
						col_count += 1

				x = tables[0]
				rows = x.find_all('tr')
				headers = rows[0].find_all('th', scope='col')
				headers = [h.text for h in headers]
				for row in rows[1:]:
					sample_id = row.find('th')
					values = [sample_id] + row.find_all('td')
					values = [v.text for v in values]
					d = {'file': summary_file}
					d.update(dict(zip(headers, values)))
					genotypemetrics.append(d)

			os.remove(summary_file)
			print(summary_file + ' removed')

		final_values = {}
		extra_values = {}

		for prop, value in report_json.items():
			if prop in schema_props:
				final_values[prop] = schemify(value, full_schema['properties'][prop]['type'])
			else:
				extra_values[prop] = value
		in_schema[file] = final_values
		out_schema[file] = extra_values

else:
	sys.exit('ERROR: --pipline not recognized, should be cellranger or dragen')

if genotypemetrics:
	df = pd.DataFrame(genotypemetrics)
	df.to_csv('genotype_metrics.tsv', sep='\t')

df = pd.DataFrame(in_schema).transpose()
df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
df.index.name = schema_version
df.to_csv(args.assay + '_metrics.tsv', sep='\t')

df = pd.DataFrame(out_schema).transpose()
df.to_csv(args.assay + '_metrics_not_in_schema.tsv', sep='\t')
