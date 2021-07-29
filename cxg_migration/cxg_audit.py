import argparse
import boto3
import gc
import json
import os
import re
import requests
import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
from random import randint
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('-f', '--file', default=None)
	args = parser.parse_args()
	return args

s3client = boto3.client("s3")

ref_files = {
	'ont': 's3://latticed-build/ontology/ontology-2021-06-25.json',
	'v2_barcodes': 's3://submissions-lattice/cellranger-whitelist/737K-august-2016.txt',
	'v3_barcodes': 's3://submissions-lattice/cellranger-whitelist/3M-february-2018.txt'
}

for s3_uri in ref_files.values():
    bucket_name = s3_uri.split('/')[2]
    file_path = s3_uri.replace('s3://{}/'.format(bucket_name), '')
    file_name = s3_uri.split('/')[-1]
    if not os.path.exists(file_name):
	    try:
	        s3client.download_file(bucket_name, file_path, file_name)
	    except subprocess.CalledProcessError as e:
	        sys.exit('ERROR: {} not found, check uri'.format(file_name))

ont_file = ref_files['ont'].split('/')[-1]
ont = json.load(open(ont_file))
ont['PATO:0000461'] = {'name': 'normal'}
ont['NCBITaxon:9606'] = {'name': 'Homo sapiens'}
ont['NCBITaxon:10090'] = {'name': 'Mus musculus'}
ont['NCBITaxon:9539'] = {'name': 'Macaca'}
ont['NCBITaxon:9483'] = {'name': 'Callithrix jacchus'}

v2_file = ref_files['v2_barcodes'].split('/')[-1]
v3_file = ref_files['v3_barcodes'].split('/')[-1]

v2_barcode_list = [line.strip() for line in open(v2_file, 'r')]
v3_barcode_list = [line.strip() for line in open(v3_file, 'r')]

barcode_pattern = '[ACTG]{16}'

obs_ont_standards = [
	'assay',
	'cell_type',
	'development_stage',
	'disease',
	'ethnicity',
	'tissue',
]

full_standards = ['sex']
for a in obs_ont_standards:
	full_standards.append(a)
	full_standards.append(a + '_ontology_term_id')

ont_dbs = {
	'organism': 'NCBITaxon',
	'assay': 'EFO',
	'cell_type': 'CL',
	'development_stage': 'HsapDv',
	'disease': 'MONDO',
	'ethnicity': 'HANCESTRO',
	'tissue': 'UBERON'
}


def upload_file(file_name):
    object_name = 'cxg_migration/original/' + file_name
    bucket = 'submissions-lattice'

    s3_client = boto3.client('s3')
    try:
        response = s3_client.upload_file(file_name, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True


def id_label_comp(ds, prop, label, term_id, org_id=None):
	if '(' in label:
		label = label.split('(')[0].strip()
		term_id = term_id.split('(')[0].strip()
		if label.split('(')[1:] != term_id.split('(')[1:]:
			report_error(ds, 'mismatched appendings between {} and {}'.format(label, term_id))
	if label == 'unknown' and prop in ['ethnicity', 'sex', 'development_stage']:
		if term_id != '':
			report_error(ds, prop + ' ' + term_id + ' not in ontology, associated with label:' + label)
	elif label == 'na' and org_id != 'NCBITaxon:9606' and prop == 'ethnicity':
		if term_id != '':
			report_error(ds, prop + ' ' + term_id + ' not in ontology, associated with label:' + label)
	elif not ont.get(term_id):
		report_error(ds, prop + ' ' + term_id + ' not in ontology, associated with label:' + label)
	elif ont[term_id]['name'] != label:
		report_error(ds, prop + ' ' + label + ' does not match ' + 'label:' + ont[term_id]['name'] + ' of ' + term_id)
	if ':' in term_id and term_id.split(':')[0] != ont_dbs[prop]:
		if prop != 'disease' and term_id != 'PATO:0000461':
			if prop == 'development_stage':
				if org_id == 'NCBITaxon:9606':
					report_error(ds, prop + ' ' + term_id + ' should be from db ' + ont_dbs[prop])
				elif term_id.split(':')[0] != 'EFO':
					report_error(ds, prop + ' ' + term_id + ' should be from db EFO')
			else:
				report_error(ds, prop + ' ' + term_id + ' should be from db ' + ont_dbs[prop])


def report_error(ds, err):
	with open('cxg_errors.txt', 'a') as out:
		out.write(ds + '\t' + err + '\n')


def check_raw(ds, raw):
	if type(raw).__name__ != 'ndarray':
		vals = list(np.unique(raw.toarray()))[:1000]
	else:
		vals = list(np.unique(raw))[:1000]
	while len(vals):
		v = vals.pop(0)
		if v - int(v) != 0:
			report_error(ds, 'non-whole numbers found in raw.X')
			break


def matrix_info(local_path):
	dataset = local_path

	# report max values of raw.X and .X
	# Look for whole numbers in raw layer then clear memory
	adata_full = sc.read_h5ad(local_path)
	if adata_full.raw:
		raw_layer = adata_full.raw.X
		raw_max = 'max is ' + str(adata_full.raw.X.max())
		raw_min = 'min is ' + str(adata_full.raw.X.min())
		raw_var_count = 'var count is ' + str(adata_full.raw.shape[1])
	else:
		raw_layer = adata_full.X
		raw_max = 'is not present'
		raw_min = 'is not present'
		raw_var_count = 'is not present'
	report_error(dataset, '.raw.X ' + raw_max + ',.X max is ' + str(adata_full.X.max()))
	report_error(dataset, '.raw.X ' + raw_min + ',.X min is ' + str(adata_full.X.min()))
	report_error(dataset, '.raw.X ' + raw_var_count + ',.X var count is ' + str(adata_full.X.shape[1]))
	del adata_full
	gc.collect()
	check_raw(dataset, raw_layer)
	del raw_layer
	gc.collect()

	adata = sc.read_h5ad(local_path, backed='r')

	# check 1000 random barcodes against 10x lists
	if re.search(barcode_pattern, adata.obs.index[5,]):
		obs_count = adata.obs.index.shape[0]
		random_indices = [randint(0, obs_count - 1) for p in range(0, 1000)]
		barcodes = {'v2': 0,'v3': 0,'both': 0,'neither': 0}
		for i in random_indices:
			if re.search(barcode_pattern, adata.obs.index[i,]):
				barcode = re.search(barcode_pattern, adata.obs.index[i,]).group(0)
				if barcode in v2_barcode_list and barcode in v3_barcode_list:
					barcodes['both'] += 1
				elif barcode in v2_barcode_list:
					barcodes['v2'] += 1
				elif barcode in v3_barcode_list:
					barcodes['v3'] += 1
				else:
					barcodes['neither'] += 1
		report_list = ['10x barcode check']
		for k,v in barcodes.items():
			if v > 0:
				report_list.append(str(v) + ' ' + k)
		report_error(dataset, ','.join(report_list))

	# check for ensembl ids
	ensembl_flag = False
	for k in adata.var_keys():
		sample_val = str(adata.var.iloc[5,][k])
		if sample_val.startswith('ENSG'):
			ensembl_flag = True
			report_error(dataset, 'SUCCESS:Ensembl IDs in var.' + k)
		elif sample_val.startswith('ENSMUSG'):
			ensembl_flag = True
			report_error(dataset, 'SUCCESS:Ensembl IDs in var.' + k)
	if ensembl_flag == False:
		report_error(dataset, 'needs Ensembl IDs')

	# check for default_embedding value in obsm_keys()
	if 'default_embedding' in adata.uns:
		de = adata.uns['default_embedding']
		if 'X_' + de not in adata.obsm_keys():
			report_error(dataset, 'uns.default_embedding:{} not in obsm keys'.format(de))

	# check title field appears exactly once
	if 'title' not in adata.uns:
		report_error(dataset, 'missing uns.title')
	elif adata.uns_keys().count('title') > 1:
		report_error(dataset, 'title appears {} times in uns keys'.format(adata.uns_keys().count('title')))

	# check for organism fields exactly once, validate labe/ID pair
	if 'organism' not in adata.uns or 'organism_ontology_term_id' not in adata.uns:
		if 'organism' not in adata.uns:
			report_error(dataset, 'missing uns.organism')
		if 'organism_ontology_term_id' not in adata.uns:
			report_error(dataset, 'missing uns.organism_ontology_term_id')
	else:
		id_label_comp(dataset, 'organism', adata.uns['organism'], adata.uns['organism_ontology_term_id'])

	# provide all EFO dev ontology terms so we know what will need to be migrated
	if 'organism_ontology_term_id' in adata.uns and 'development_stage_ontology_term_id' in adata.obs:
		org = adata.uns['organism_ontology_term_id']
		for i in adata.obs['development_stage_ontology_term_id'].value_counts().to_dict().keys():
			if i.startswith('EFO:'):
				with open('cxg_efo_dev.txt', 'a') as out:
					out.write(dataset + '\t' + org + '\t' + i + '\n')
	elif 'organism_ontology_term_id' in adata.uns and 'development_stage' in adata.obs:
		org = adata.uns['organism_ontology_term_id']
		for i in adata.obs['development_stage'].value_counts().to_dict().keys():
			if i.startswith('EFO:'):
				with open('cxg_efo_dev.txt', 'a') as out:
					out.write(dataset + '\t' + org + '\t' + i + '\n')

	# check sex field appears exactly once and contains valid values
	if 'sex' not in adata.obs:
		report_error(dataset, 'missing obs.sex')
	else:
		if adata.obs.dtypes['sex'].name != 'category':
			report_error(dataset, 'sex dtype is {}, not category'.format(adata.obs.dtypes['sex'].name))
		if adata.obs_keys().count('sex') > 1:
			report_error(dataset, 'sex appears {} times in obs keys'.format(adata.obs_keys().count('sex')))
		for k in adata.obs['sex'].value_counts().to_dict().keys():
			if k not in ['male', 'female', 'unknown']:
				report_error(dataset, 'obs.sex:{} invalid'.format(k))

	# check remaining schema fields
	for o in obs_ont_standards:
		ont_field = o + '_ontology_term_id'
		if o in adata.obs and ont_field in adata.obs:
			if adata.obs_keys().count(o) > 1:
				report_error(dataset, '{} appears {} times in obs keys'.format(o, adata.obs_keys().count(o)))
			if adata.obs_keys().count(ont_field) > 1:
				report_error(dataset, '{} appears {} times in obs keys'.format(ont_field, adata.obs_keys().count(ont_field)))
			for k in adata.obs[[o,ont_field]].value_counts().to_dict().keys():
				id_label_comp(dataset, o, k[0],k[1], adata.uns['organism_ontology_term_id'])
		else:
			if o not in adata.obs and ont_field not in adata.obs:
				report_error(dataset, 'obs.{} missing'.format(o))
				report_error(dataset, 'obs.{} missing'.format(ont_field))
			elif o not in adata.obs:
				report_error(dataset, 'obs.{} missing'.format(o))
				for k in adata.obs[ont_field].value_counts().to_dict().keys():
					id_label_comp(dataset, o, 'missing', k, adata.uns['organism_ontology_term_id'])
			else:
				report_error(dataset, 'obs.{} missing'.format(ont_field))
				for k in adata.obs[o].value_counts().to_dict().keys():
					id_label_comp(dataset, o, k, 'missing', adata.uns['organism_ontology_term_id'])

	uber_dict = {}
	heatmap_fields = []
	obj_fields = []
	long_fields = []
	for o in adata.obs.keys():
		value_counts_dict = adata.obs[o].value_counts().to_dict()
		counts = '_'.join([str(c) for c in value_counts_dict.values()])
		values = [str(i) for i in value_counts_dict.keys()]

		# look for fields with only 1 value of 'null' or equivalent
		if len(value_counts_dict.keys()) == 1 and o not in full_standards:
			report_error(dataset, 'all values for ' + o + ' are ' + str(list(value_counts_dict.keys())[0]))

		# look for a field with age values, that can inform development ontology term
		age_flag = False
		if 'development_stage' not in o:
			for i in value_counts_dict.keys():
				if 'year' in str(i).lower():
					age_flag = True
			if 'age' in o.lower() or 'year' in o.lower():
				age_flag = True
		if age_flag == True:
			report_error(dataset, 'possible age field - ' + o + ' - ' + ','.join(values))

		# check for heatmap fields to ensure they 'make sense' (not cluster ID)
		if adata.obs.dtypes[o].name in ['float32', 'float64', 'int32']:
			heatmap_fields.append(o)
		else:
			# check for long categories as they will not be enabled for coloring
			if len(value_counts_dict.keys()) > 200:
				long_fields.append(o + ':' + str(len(value_counts_dict.keys())))
			# check for dtype object fields as they will not be enabled for coloring
			if adata.obs.dtypes[o].name in ['object']:
				obj_fields.append(o)

			# look for a field with suspension type values - possible optional schema
			susp_flag = False
			for i in value_counts_dict.keys():
				if 'nucle' in str(i).lower() or 'cell' in str(i).lower():
					susp_flag = True
			if susp_flag == True:
				report_error(dataset, 'possible suspension field - ' + o + ' - ' + ','.join(values))

			# report value_counts to later look for redundancy
			metadata = {
				'values': values,
				'property': o
			}
			if counts in uber_dict:
				uber_dict[counts].append(metadata)
			else:
				uber_dict[counts] = [metadata]

	# comb value_counts to report possible redundancy
	for k,v in uber_dict.items():
		if '_' in k and not k.startswith('1_1'):
			props = [e['property'] for e in v]
			if len(v) > 1 and not all(elem in full_standards for elem in props):
				for e in v:
					with open('cxg_redundant.txt', 'a') as out:
						out.write(dataset + '\t' + k + '\t' + e['property'] + '\t' + ','.join(e['values']) + '\n')

	# report various collected fields
	if heatmap_fields:
		report_error(dataset, 'heatmap fields - ' + ','.join(heatmap_fields))
	if obj_fields:
		report_error(dataset, 'object fields - ' + ','.join(obj_fields))
	if long_fields:
		report_error(dataset, 'long fields - ' + ','.join(long_fields))

args = getArgs()
if args.file:
	matrix_info(args.file)
else:
	s3 = boto3.resource('s3')
	my_bucket = s3.Bucket('submissions-lattice')
	for s3_object in my_bucket.objects.filter(Prefix="cxg_migration/original"):
		filename = os.path.split(s3_object.key)[1]
		if filename:
			my_bucket.download_file(s3_object.key, filename)
			matrix_info(filename)
			os.remove(filename)
