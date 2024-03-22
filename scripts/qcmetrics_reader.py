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
from qcmetrics_mapper import mappings
from urllib.parse import urljoin
from bs4 import BeautifulSoup


EPILOG = '''
Extract summary info and QC metrics from a cellranger pipeline run.

Examples:

    python %(prog)s -m production -a atac -p cr -d submissions-czi009kid/muto_humphreys_2020/Control_5/outs
    python %(prog)s -m local -a rna -p cr -d submissions-czi012eye/chen_2020/19D013_foveaR_outs

For more details:

        python %(prog)s --help
'''


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('--dir', '-d',
						required=True,
						help="s3 path to the cellranger outs directory or dragen html")
	parser.add_argument('--assay', '-a',
						required=True,
						help="specify rna, atac, multiome, or spatial")
	parser.add_argument('--mode', '-m',
						required=True,
						help='The machine to pull schema from.')
	parser.add_argument('--pipeline', '-p',
						required=True,
						help='specify cr or cellranger for CellRanger, dragen, or star')
	args = parser.parse_args()
	return args


def read_cr_csv(file):
	s3client.download_file(bucket_name, outs_dir_path + '/' + file, file)
	print(file + ' downloaded')
	with open(file, newline='') as csvfile:
		spamreader = csv.reader(csvfile)
		rows = list(spamreader)
		headers = rows[0]
		values = rows[1]
		temp_json = dict(zip(headers, values))

	os.remove(file)
	print(file + ' removed')

	return temp_json


def read_star_csv(file):
	s3client.download_file(bucket_name, outs_dir_path + '/' + file, file)
	print(file + ' downloaded')
	temp_json = {}
	with open(file, newline='') as csvfile:
		spamreader = csv.reader(csvfile)
		rows = list(spamreader)
		for row in rows:
			temp_json[row[0]] = row[1]

	os.remove(file)
	print(file + ' removed')

	return temp_json


def read_cr_json(file):
	s3client.download_file(bucket_name, outs_dir_path + '/' + file, file)
	print(file + ' downloaded')
	with open(file) as summary_json:
		temp_json = json.load(summary_json)

	os.remove(file)
	print(file + ' removed')

	return temp_json


def read_cr_html(file):
	s3client.download_file(bucket_name, outs_dir_path + '/' + file, file)
	print(file + ' downloaded')
	temp_json = {}
	with open(file) as html_doc:
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
				elif data.get('joint_pipeline_info_table'):
					pipeline_info_table = data.get('joint_pipeline_info_table')
				else:
					pipeline_info_table = data['summary']['summary_tab']['pipeline_info_table']
				info_list = pipeline_info_table['rows']
				for pair in info_list:
					temp_json[pair[0]] = pair[1]
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
					temp_json[field] = value

	os.remove(file)
	print(file + ' removed')

	return temp_json


def fractionize(v):
	dig = len(v.replace('.',''))
	v = float(v)
	if v < 10: #WHAT IF V IS < 1?
		v = round(v/100,dig+1)
	else:
		v = round(v/100,dig)
	return v


def schemify(value, prop_type):
	if value in [None, '', 'nan']:
		return ''
	elif (prop_type == 'integer') and (str(value).endswith('.0') == True):
		return str(value)[:-2]
	else:
		return str(value)

pipe_map = {
	'cr': 'cellranger',
	'starsolo': 'star'
}

output_dir = 'outputs'
qcmtrics_dir = 'qcmetrics'

# Checking for presence / creating output folder and associated sub-folders
if os.path.exists(output_dir) == False:
	os.mkdir(output_dir)
if os.path.exists(output_dir + '/' + qcmtrics_dir) == False:
	os.mkdir(output_dir + '/' + qcmtrics_dir)

args = getArgs()

assay = args.assay

pipeline = args.pipeline.lower()
pipeline = pipe_map.get(pipeline, pipeline)

s3client = boto3.client("s3")

connection = lattice.Connection(args.mode)
server = connection.server

dir_list = args.dir
if os.path.isfile(dir_list):
	directories = [line.rstrip('\n') for line in open(dir_list)]
else:
    directories = dir_list.split(',')

assays = [
	'rna',
	'atac',
	'antibody_capture',
	'multiome',
	'spatial'
]

schemas = {}
for s in assays:
	url = urljoin(server, f'profiles/{s}_metrics/?format=json')
	schema = requests.get(url).json()
	props = schema['properties']
	schemas[s] = props


in_schema = {}
in_rna_schema = {}
in_atac_schema = {}
in_ac_schema = {}
in_mu_schema = {}
in_sp_schema = {}
out_schema = {}
genotypemetrics = []

if pipeline in ['cellranger', 'star']:
	value_mapping = mappings[pipeline].get('value_mapping',{})
	schema_mapping = mappings[pipeline][assay].get('schema_mapping',{})
	should_match = mappings[pipeline][assay].get('should_match',{})
	perc_to_frac = mappings[pipeline][assay].get('perc_to_frac',[])

	for direct in directories:
		print('starting: {}'.format(direct))
		direct = direct.replace('s3://', '')
		full_path = direct.rstrip('/')
		bucket_name = full_path.split('/')[0]
		outs_dir_path = full_path.replace(bucket_name + '/', '')

		report_json = {}

		summary_file = 'web_summary.html'
		objects = s3client.list_objects_v2(Bucket=bucket_name,Prefix=outs_dir_path + '/')
		summaries = [o['Key'] for o in objects['Contents'] if o['Key'].lower().endswith(summary_file)]
		if len(summaries) == 1:
			summary_file = summaries[0].split('/')[-1]
			report_json.update(read_cr_html(summary_file))
		elif len(summaries) > 1:
			print('WARNING: multiple {} files found on s3'.format(summary_file))
		elif len(summaries) == 0:
			print('WARNING: no {} files found on s3'.format(summary_file))

		metrics_json_file = 'summary.json'
		jsons = [o['Key'] for o in objects['Contents'] if o['Key'].lower().endswith(metrics_json_file)]
		if len(jsons) == 1:
			metrics_json_file = jsons[0].split('/')[-1]
			report_json.update(read_cr_json(metrics_json_file))
		elif len(jsons) > 1:
			print('WARNING: multiple {} files found on s3'.format(metrics_json_file))
		elif len(jsons) == 0 and assay == 'atac':
				print('WARNING: no {} files found on s3'.format(metrics_json_file))

		metrics_csv_file = 'summary.csv'
		csvs = [o['Key'] for o in objects['Contents'] if o['Key'].lower().endswith(metrics_csv_file)]
		if len(csvs) == 1:
			metrics_csv_file = csvs[0].split('/')[-1]
			if pipeline == 'star':
				report_json.update(read_star_csv(metrics_csv_file))
			else:
				report_json.update(read_cr_csv(metrics_csv_file))
		elif len(csvs) > 1:
			print(csvs)
			print('WARNING: multiple {} files found on s3'.format(metrics_csv_file))
		elif len(csvs) == 0:
			print('WARNING: no {} files found on s3'.format(metrics_csv_file))

		rna_values = {}
		atac_values = {}
		ac_values = {}
		mu_values = {}
		extra_values = {}
		final_values = {}


		for k,v in report_json.items():
			if v is None:
				try_v = ''
			else:
				try_v = str(v)
				try_v = try_v.strip('%')

			try_k = k.lower().replace(' ','_')
			try_k = schema_mapping.get(try_k, try_k)

			if try_k in schemas[assay]:
				if try_k in perc_to_frac:
					try_v = fractionize(try_v)
				final_values[try_k] = schemify(try_v, schemas[assay][try_k]['type'])

			elif assay == 'rna' and try_k.startswith('antibody:_'):
				try_k = '_'.join(try_k.split('_')[1:])
				try_k = mappings[pipeline]['antibody_capture']['schema_mapping'].get(try_k, try_k)
				if try_k in schemas['antibody_capture']:
					if try_k in mappings[pipeline]['antibody_capture']['perc_to_frac']:
						try_v = fractionize(try_v)
					ac_values[try_k] = schemify(try_v, schemas['antibody_capture'][try_k]['type'])
				else:
					extra_values[k] = v

			elif assay in 'multiome' and try_k.startswith('gex_'):
				try_k = '_'.join(try_k.split('_')[1:])
				try_k = mappings[pipeline]['rna']['schema_mapping'].get(try_k, try_k)
				if try_k in schemas['rna']:
					rna_values[try_k] = schemify(try_v, schemas['rna'][try_k]['type'])
				else:
					extra_values[k] = v

			elif assay == 'multiome' and try_k.startswith('atac_'):
				try_k = '_'.join(try_k.split('_')[1:])
				try_k = mappings[pipeline]['atac']['schema_mapping'].get(try_k, try_k)
				if try_k in schemas['atac']:
					atac_values[try_k] = schemify(try_v, schemas['atac'][try_k]['type'])
				else:
					extra_values[k] = v

			elif assay == 'multiome' and try_k == 'total_cells_detected':
				rna_values[try_k] = schemify(try_v, schemas['rna'][try_k]['type'])
				atac_values[try_k] = schemify(try_v, schemas['atac'][try_k]['type'])

			elif k in schemas[assay]:
				if k in perc_to_frac:
					try_v = fractionize(try_v)
				final_values[try_k] = schemify(try_v, schemas[assay][k]['type'])

			elif assay in 'spatial':
				try_k = mappings[pipeline]['rna']['schema_mapping'].get(try_k, try_k)
				if try_k in schemas['rna']:
					rna_values[try_k] = schemify(try_v, schemas['rna'][try_k]['type'])
				else:
					extra_values[k] = v

			elif type(v) in [dict]:
				print('WARNING: {} is a {}, not included in report'.format(k, type(v).__name__))

			else:
				if isinstance(v, list):
					v = ','.join(v)
				extra_values[k] = value_mapping.get(v, v)

		for k,v in should_match.items():
			if k in extra_values and v in final_values:
				if float(extra_values[k]) != float(final_values[v]):
					print('ERROR: {}:{} does not match {}:{}'.format(k, str(extra_values[k]),v,str(final_values[v])))
				else:
					del extra_values[k]
 
		if final_values:
			final_values['quality_metric_of'] = '<linkTo RawMatrixFile - filtered matrix .h5>'
			in_schema[direct] = final_values
		if rna_values:
			rna_values['quality_metric_of'] = '<linkTo RawMatrixFile - filtered matrix .h5>'
			in_rna_schema[direct] = rna_values
		if atac_values:
			atac_values['quality_metric_of'] = '<linkTo RawMatrixFile - filtered matrix .h5>'
			in_atac_schema[direct] = atac_values
		if mu_values:
			mu_values['quality_metric_of'] = '<linkTo RawMatrixFile - filtered matrix .h5>'
			in_mu_schema[direct] = mu_values
		if ac_values:
			if (final_values and 'total_cells_detected' in final_values):
				ac_values['total_cells_detected'] = final_values['total_cells_detected']
			ac_values['quality_metric_of'] = '<linkTo RawMatrixFile - filtered matrix .h5>'
			in_ac_schema[direct] = ac_values
		out_schema[direct] = extra_values


elif pipeline == 'dragen':
	schema_mapping = mappings[pipeline][assay]['schema_mapping']
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
							report_json[field] = column.get_text().strip()
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
			if prop in mappings[pipeline][assay]['perc_to_frac']:
				dig = len(value.replace('.',''))
				value = float(value)
				if value < 10:
					value = round(value/100,dig+1)
				else:
					value = round(value/100,dig)
			if prop in schemas['rna']:
				final_values[prop] = schemify(value, schemas['rna'][prop]['type'])
			else:
				extra_values[prop] = value
		in_rna_schema[file] = final_values
		out_schema[file] = extra_values

else:
	sys.exit('ERROR: --pipeline not recognized, should be cellranger, star or dragen')

if in_schema:
	df = pd.DataFrame(in_schema).transpose()
	df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
	df.index.name = 'schema_version=' + str(schemas[assay]['schema_version']['default'])
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + assay + '_metrics.tsv', sep='\t')

if genotypemetrics:
	df = pd.DataFrame(genotypemetrics)
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'genotype_metrics.tsv', sep='\t', index=False)

if in_atac_schema:
	df = pd.DataFrame(in_atac_schema).transpose()
	df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
	df.index.name = 'schema_version=' + str(schemas['atac']['schema_version']['default'])
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'atac_metrics.tsv', sep='\t')

if in_rna_schema:
	df = pd.DataFrame(in_rna_schema).transpose()
	df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
	df.index.name = 'schema_version=' + str(schemas['rna']['schema_version']['default'])
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'rna_metrics.tsv', sep='\t')

if in_mu_schema:
	df = pd.DataFrame(in_mu_schema).transpose()
	df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
	df.index.name = 'schema_version=' + str(schemas['multiome']['schema_version']['default'])
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'multiome_metrics.tsv', sep='\t')

if in_ac_schema:
	df = pd.DataFrame(in_ac_schema).transpose()
	df = df[['quality_metric_of'] + [col for col in df.columns if col != 'quality_metric_of']]
	df.index.name = 'schema_version=' + str(schemas['antibody_capture']['schema_version']['default'])
	df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'antibody_capture_metrics.tsv', sep='\t')

df = pd.DataFrame(out_schema).transpose()
df.to_csv(output_dir + '/' + qcmtrics_dir + '/' + 'metrics_not_in_schema.tsv', sep='\t')
