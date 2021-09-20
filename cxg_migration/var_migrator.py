import boto3
import gc
import gzip
import json
import os
import requests
import subprocess
import gene_symbol_custom
import anndata as ad
import pandas as pd
import scanpy as sc
import pandas as pd


def fixup(df):
	hgnc = 'hgnc_complete_set.txt'
	file_path = 'cxg_migration/var_refs/' + hgnc
	if not os.path.exists(hgnc):
		try:
			client.download_file(bucket_name, file_path, hgnc)
		except subprocess.CalledProcessError as e:
			sys.exit('ERROR: {} not found, check uri'.format(hgnc))
	df = gene_symbol_custom.get_upgraded_var_index(df, hgnc)

	return df


def validate_cxg(ds):
	file = ds + '.h5ad'
	validate_process = subprocess.run(['cellxgene-schema', 'validate', file], stdout=subprocess.PIPE)
	with open('validate_logs.txt', 'a') as f:
		for line in validate_process.stdout.decode('utf-8').split('\n'):
			if line:
				f.write(ds + '\t' + line + '\n')


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
	# https://github.com/chanzuckerberg/single-cell-curation/tree/main/cellxgene_schema_cli/cellxgene_schema/ontology_files
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


def curate_var(var, strategy, fixedup):
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
		if fixedup:
			map_df = fixup(map_df)
		var = var.merge(map_df,left_index=True,right_on='gene_symbols',how='left').set_index(var.index)
		field='gene_ids'
	no_gene_id = var.index[var[field].isnull()]
	var_to_keep = list(set(var_to_keep) - set(no_gene_id))
	counts['not_mapped'] = str(len(no_gene_id))
	no_gene_id_samples = [str(i) for i in no_gene_id]
	counts['not_mapped_samples'] = ','.join(no_gene_id_samples)
	var = var.rename(columns={field:'feature_id'})

	# filter if symbol mapped to multiple genes
	before = len(var_to_keep)
	multi_mapping = var.index[var['feature_id'] == 'multiple'].tolist()
	var_to_keep = list(set(var_to_keep) - set(multi_mapping))
	counts['multiple_ids'] = str(before - len(var_to_keep))
	multiple_id_samples = [str(i) for i in multi_mapping]
	counts['multiple_id_samples'] = ','.join(multiple_id_samples)

	# filter if ID is duplicated within df
	before = len(var_to_keep)
	dups = var.index[var.duplicated(subset='feature_id',keep=False)].tolist()
	dups_still_kept = [str(i) for i in dups if i in var_to_keep]
	var_to_keep = list(set(var_to_keep) - set(dups))
	counts['duplicate_ids'] = str(before - len(var_to_keep))
	counts['duplicate_ids_samples'] = ','.join(dups_still_kept)

	# filter on approved annotation references
	before = len(var_to_keep)
	approved = pd.read_csv('approved_ids.csv',dtype='str')['feature_id']
	var_in_approved = var.index[var['feature_id'].isin(approved)].tolist()
	not_approved_samples = [str(i) for i in var_to_keep if i not in var_in_approved]
	var_to_keep = [e for e in var_to_keep if e in var_in_approved]
	counts['not_approved'] = str(before - len(var_to_keep))
	counts['not_approved_samples'] = ','.join(not_approved_samples)

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
	redundant_props = ['features','ensemblid','gene_ids','gene_ids_original','feature_types','feature_type','type',
		'feature_types-Harvard-Nuclei','feature_types-Sanger-CD45','feature_types-Sanger-Cells','feature_types-Sanger-Nuclei',
		'gene_ids-Harvard-Nuclei','gene_ids-Sanger-CD45','gene_ids-Sanger-Cells','gene_ids-Sanger-Nuclei','hgnc_gene_symbol']
	remove_var = []
	for k in var.keys():
		if k in portal_props + redundant_props:
			remove_var.append(k)
	var = var.drop(columns=remove_var)

	return var, var_to_keep, counts


def main(ds, strategy, fixup_genes):
	adata = sc.read_h5ad(ds + '.h5ad')

	if fixup_genes == 'Y':
		fixedup = True
	else:
		fixedup = False

	if adata.raw:
		raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
		var, to_keep, counts = curate_var(raw_adata.var, strategy, fixedup)
		counts['dataset'] = ds + '-raw'
		report(counts)
		raw_adata.var = var
		raw_adata = raw_adata[:, to_keep]
		raw_adata.var.set_index('feature_id', inplace=True)
		adata.raw = raw_adata
		del raw_adata
		gc.collect()

	var, to_keep, counts = curate_var(adata.var, strategy, fixedup)
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
	'7edef704-f63a-462c-8636-4bc86a9472bd_b83559d1-156f-4ba9-9f6a-b165f83ef43f', # Voigt/Scheetz retina, no raw counts
	'a238e9fa-2bdf-41df-8522-69046f99baff_66d15835-5dc8-4e96-b0eb-f48971cb65e8', # Enge pancreas, cell don't group by cluster/cell_type
	'9b02383a-9358-4f0f-9795-a891ec523bcc_13a027de-ea3e-432b-9a5e-6bc7048498fc' # Lattice dataset, not yet in working/
	]

# get the specified mapping strategy/file for each dataset
sheet_id = '18e5PG2wCaN8kf9-KVm_yomgEx8TYka0Ldd7_swVxiJk'
sheet_name = 'datasets'
url = 'https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}'.format(sheet_id, sheet_name)
ds_df = pd.read_csv(url)[['coll_ds','var_mapping','fixup_genes']]
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
		main(ds, row['var_mapping'], row['fixup_genes'])
		validate_cxg(ds)
		client.upload_file(file, bucket_name, 'cxg_migration/final/' + file, ExtraArgs={'ACL':'public-read'})
		os.remove(file)
