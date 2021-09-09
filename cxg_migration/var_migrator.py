import boto3
import gc
import gzip
import json
import os
import requests
import subprocess
import anndata as ad
import pandas as pd
import scanpy as sc
import pandas as pd


def report(counts):
	props = [
		'dataset',
		'starting',
		'not_mapped',
		'not_mapped_samples',
		'multiple_ids',
		'multiple_id_samples',
		'duplicate_ids',
		'duplicate_ids_samples',
		'not_approved',
		'not_approved_samples',
		'final'
	]
	report_out = [counts[e] for e in props]
	outfile = 'var_report.tsv'
	with open(outfile, 'a') as f:
		f.write('\t'.join(report_out) + '\n')


def compile_annotations():
	ref_files = [
		'genes_ercc.csv',
		'genes_homo_sapiens.csv',
		'genes_mus_musculus.csv',
		'genes_sars_cov_2.csv'
	]
	
	ids = pd.DataFrame()
	for f in ref_files:
		client.download_file(bucket_name, 'cxg_migration/var_refs/' + f, f)
		df = pd.read_csv(f, names=['feature_id','symb','num'])
		ids = ids.append(df)
		os.remove(f)

	ids.to_csv('approved_ids.csv',index=False)


def fixup_var(var, strategy):
	counts = {}
	var_len_orig = len(var)
	counts['starting'] = str(var_len_orig)
	
	var_to_keep = var.index.tolist()

	# feature_id, will become index
	# if IDs in var already
	if strategy.startswith('var.'):
		field = '.'.join(strategy.split('.')[1:])
	# if mapping IDs from data source
	else:
		if 'gene_ids' in var.keys():
			var.rename(columns={'gene_ids':'gene_ids_original'},inplace=True)
		file_path = 'cxg_migration/var_refs/' + strategy
		if not os.path.exists(strategy):
			try:
				client.download_file(bucket_name, file_path, strategy)
			except subprocess.CalledProcessError as e:
				sys.exit('ERROR: {} not found, check uri'.format(strategy))
		map_df = pd.read_csv(strategy,sep='\t')
		var = var.merge(map_df,left_index=True,right_on='gene_symbols',how='left').set_index(var.index)
		field='gene_ids'
	no_gene_id = var.index[var[field].isnull()]
	var_to_keep = list(set(var_to_keep) - set(no_gene_id))
	counts['not_mapped'] = str(len(no_gene_id))
	no_gene_id_samples = [str(i) for i in no_gene_id[0:10]]
	counts['not_mapped_samples'] = '_'.join(no_gene_id_samples)
	var = var.rename(columns={field:'feature_id'})

	# filter if symbol mapped to multiple genes
	before = len(var_to_keep)
	multi_mapping = var.index[var['feature_id'] == 'multiple'].tolist()
	var_to_keep = list(set(var_to_keep) - set(multi_mapping))
	counts['multiple_ids'] = str(before - len(var_to_keep))
	multiple_id_samples = [str(i) for i in multi_mapping[0:10]]
	counts['multiple_id_samples'] = '_'.join(multiple_id_samples)

	# filter if ID is duplicated within df
	before = len(var_to_keep)
	dups = var.index[var.duplicated(subset='feature_id',keep=False)].tolist()
	dups_still_kept = [str(i) for i in dups if i in var_to_keep]
	var_to_keep = list(set(var_to_keep) - set(dups))
	counts['duplicate_ids'] = str(before - len(var_to_keep))
	counts['duplicate_ids_samples'] = '_'.join(dups_still_kept)

	# filter on approved annotation references
	before = len(var_to_keep)
	approved = pd.read_csv('approved_ids.csv',dtype='str')['feature_id']
	var_in_approved = var.index[var['feature_id'].isin(approved)].tolist()
	not_approved_samples = [str(i) for i in var_to_keep if i not in var_in_approved][0:10]
	var_to_keep = [e for e in var_to_keep if e in var_in_approved]
	counts['not_approved'] = str(before - len(var_to_keep))
	counts['not_approved_samples'] = '_'.join(not_approved_samples)

	counts['final'] = str(len(var_to_keep))

	# feature_biotype
	# add 'gene' to all
	var['feature_biotype'] = 'gene'
	# merge with a 'spike-in' dataframe to update just those to 'spike-in'
	spikein_id_file = 'cms_095047_identifiers.txt'
	file_path = 'cxg_migration/var_refs/' + spikein_id_file
	if not os.path.exists(spikein_id_file):
		try:
			client.download_file(bucket_name, file_path, spikein_id_file)
		except subprocess.CalledProcessError as e:
			sys.exit('ERROR: {} not found, check uri'.format(spikein_id_file))
	ercc_df = pd.read_csv(spikein_id_file)
	ercc_df['feature_biotype'] = 'spike-in'
	var.loc[var.feature_id.isin(ercc_df.feature_id), ['feature_biotype']] = ercc_df[['feature_biotype']]

	# remove columns from var
	portal_props = ['feature_reference','feature_name','gene_symbols']
	redundant_props = ['feature_types','feature_type','feature_types-Harvard-Nuclei','feature_types-Sanger-CD45','feature_types-Sanger-Cells','feature_types-Sanger-Nuclei',
		'features','ensemblid','gene_ids','gene_ids-Harvard-Nuclei','gene_ids-Sanger-CD45','gene_ids-Sanger-Cells','gene_ids-Sanger-Nuclei','hgnc_gene_symbol','type']
	remove_var = []
	for k in var.keys():
		if k in portal_props + redundant_props:
			remove_var.append(k)
	var = var.drop(columns=remove_var)

	return var, var_to_keep, counts


def main(ds, strategy):
	adata = sc.read_h5ad(ds + '.h5ad')

	if adata.raw:
		raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
		var, to_keep, counts = fixup_var(raw_adata.var, strategy)
		counts['dataset'] = ds + '-raw'
		report(counts)
		raw_adata.var = var
		raw_adata = raw_adata[:, to_keep]
		raw_adata.var.set_index('feature_id', inplace=True)
		adata.raw = raw_adata
		del raw_adata
		gc.collect()

	var, to_keep, counts = fixup_var(adata.var, strategy)
	var['feature_is_filtered'] = False # feature_is_filtered is default False
	counts['dataset'] = ds + '-X'
	report(counts)
	adata.var = var
	adata = adata[:, to_keep]
	adata.var.set_index('feature_id', inplace=True)

	# write the new object to the file
	adata.write(filename=ds + '.h5ad')
	del adata
	gc.collect()


# avoid these datasets for now
attn_needed = [
	'7edef704-f63a-462c-8636-4bc86a9472bd_b83559d1-156f-4ba9-9f6a-b165f83ef43f', # no raw counts
	'f70ebd97-b3bc-44fe-849d-c18e08fe773d_e0ed3c55-aff6-4bb7-b6ff-98a2d90b890c', # no raw counts, 2 non-raw layers
	'a238e9fa-2bdf-41df-8522-69046f99baff_66d15835-5dc8-4e96-b0eb-f48971cb65e8', # cell don't group by cluster/cell_type
	'9b02383a-9358-4f0f-9795-a891ec523bcc_13a027de-ea3e-432b-9a5e-6bc7048498fc', # Lattice dataset, not yet in working/
	'00109df5-7810-4542-8db5-2288c46e0424_fe2479fd-daff-41a8-97dc-a50457ab1871', # adata.write(filename=ds + '.h5ad') numpy.core._exceptions.MemoryError: Unable to allocate 59.9 GiB for an array with shape (292010, 55050) and data type float32
	'c114c20f-1ef4-49a5-9c2e-d965787fb90c_f7c1c579-2dc0-47e2-ba19-8165c5a0e353', # adata.write(filename=ds + '.h5ad') numpy.core._exceptions.MemoryError: Unable to allocate 17.4 GiB for an array with shape (2332585730,) and data type int64
	'0a839c4b-10d0-4d64-9272-684c49a2c8ba_9dbab10c-118d-496b-966a-67f1763a6b7d', # adata.write(filename=ds + '.h5ad') numpy.core._exceptions.MemoryError: Unable to allocate 17.3 GiB for an array with shape (2315747958,) and data type int64
	'e5f58829-1a66-40b5-a624-9046778e74f5_5a11f879-d1ef-458a-910c-9b0bdfca5ebf', # raw.var is just index (numerical, no metadata to migrate)
	'e5f58829-1a66-40b5-a624-9046778e74f5_97a17473-e2b1-4f31-a544-44a60773e2dd', # raw.var is just index (numerical, no metadata to migrate)
	'e5f58829-1a66-40b5-a624-9046778e74f5_c5d88abe-f23a-45fa-a534-788985e93dad', # raw.var is just index (numerical, no metadata to migrate)
	'e5f58829-1a66-40b5-a624-9046778e74f5_a68b64d8-aee3-4947-81b7-36b8fe5a44d2', # raw.var is just index (numerical, no metadata to migrate)
	'e5f58829-1a66-40b5-a624-9046778e74f5_53d208b0-2cfd-4366-9866-c3c6114081bc', # raw.var is just index (numerical, no metadata to migrate)
	'38833785-fac5-48fd-944a-0f62a4c23ed1_2adb1f8a-a6b1-4909-8ee8-484814e2d4bf' # adata.raw = raw_adata numpy.core._exceptions.MemoryError: Unable to allocate 32.5 GiB for an array with shape (14533, 599926) and data type float32
	]

# get the specified mapping strategy/file for each dataset
sheet_id = '18e5PG2wCaN8kf9-KVm_yomgEx8TYka0Ldd7_swVxiJk'
sheet_name = 'datasets'
url = 'https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}'.format(sheet_id, sheet_name)
ds_df = pd.read_csv(url)[['coll_ds','var_mapping']]
ds_df = ds_df.loc[ds_df['var_mapping'] != 'Lattice dataset'].dropna()

already_run = []
client = boto3.client('s3')
resource = boto3.resource('s3')
bucket_name = 'submissions-lattice'
your_bucket = resource.Bucket(bucket_name)
for s3_file in your_bucket.objects.all():
	if 'cxg_migration/final' in s3_file.key:
		already_run.append(s3_file.key.split('/')[-1].split('.')[0])

# compile approved feature_ids if not already local
if not os.path.exists('approved_ids.csv'):
	compile_annotations()

outfile = 'var_report.tsv'
props = [
	'dataset',
	'starting',
	'not_mapped',
	'not_mapped_samples',
	'multiple_ids',
	'multiple_id_samples',
	'duplicate_ids',
	'duplicate_ids_samples',
	'not_approved',
	'not_approved_samples',
	'final'
]
with open(outfile, 'a') as f:
	f.write('\t'.join(props) + '\n')

for index,row in ds_df.iterrows():
	ds = row['coll_ds']
	if ds not in attn_needed and ds not in already_run:
		print('PROCESSING:' + ds)
		file = ds + '.h5ad'
		client.download_file('submissions-lattice', 'cxg_migration/working/' + file, file)
		main(ds, row['var_mapping'])
		client.upload_file(file, bucket_name, 'cxg_migration/final/' + file)
		os.remove(file)
