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


def curate_var(var, strategy):
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
		if 'hgnc_gene_symbol' in var.keys():
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


def main(ds, strategy):
	adata = sc.read_h5ad(ds + '.h5ad')

	if adata.raw:
		raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
		var, to_keep, counts = curate_var(raw_adata.var, strategy)
		counts['dataset'] = ds + '-raw'
		report(counts)
		raw_adata.var = var
		raw_adata = raw_adata[:, to_keep]
		raw_adata.var.set_index('feature_id', inplace=True)
		adata.raw = raw_adata
		del raw_adata
		gc.collect()

	var, to_keep, counts = curate_var(adata.var, strategy)
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

next_up = [
	'24d42e5e-ce6d-45ff-a66b-a3b3b715deaf_01209dce-3575-4bed-b1df-129f57fbc031',
	'cdfb9ead-cb58-4a53-879d-5e4ed5329e73_077b0429-0f47-48e0-879a-39eaae531d42',
	'9132fae8-bdfe-480f-9e45-45bc77f320b3_257adc73-8152-414b-a2c7-73861b8e0c0a',
	'2a79d190-a41e-4408-88c8-ac5c4d03c0fc_030faa69-ff79-4d85-8630-7c874a114c19',
	'f70ebd97-b3bc-44fe-849d-c18e08fe773d_e0ed3c55-aff6-4bb7-b6ff-98a2d90b890c',
	'9c8808ce-1138-4dbe-818c-171cff10e650_cfa3c355-ee77-4fc8-9a00-78e61d23024c',
	'9c8808ce-1138-4dbe-818c-171cff10e650_26ae14da-9e5f-4d18-abae-18a5a328feef',
	'ed9185e3-5b82-40c7-9824-b2141590c7f0_30cd5311-6c09-46c9-94f1-71fe4b91813c',
	'ed9185e3-5b82-40c7-9824-b2141590c7f0_21d3e683-80a4-4d9b-bc89-ebb2df513dde',
	'b953c942-f5d8-434f-9da7-e726ba7c1481_9813a1d4-d107-459e-9b2e-7687be935f69',
	'b953c942-f5d8-434f-9da7-e726ba7c1481_85c60876-7f35-40c5-a256-7808d84c6ba5',
	'e2a4a67f-6a18-431a-ab9c-6e77dd31cc80_ea426edb-4e86-4c53-ab17-5b952d94a31e',
	'e2a4a67f-6a18-431a-ab9c-6e77dd31cc80_e2a3c32d-71e2-4f38-b19c-dfcb8729cf46',
	'e2a4a67f-6a18-431a-ab9c-6e77dd31cc80_c3fe3c1e-5bf8-4678-b74a-79899243ad41',
	'367d95c0-0eb0-4dae-8276-9407239421ee_9b686bb6-1427-4e13-b451-7ee961115cf9',
	'367d95c0-0eb0-4dae-8276-9407239421ee_b6203114-e133-458a-aed5-eed1028378b4',
	'367d95c0-0eb0-4dae-8276-9407239421ee_f7a068f1-0fdb-48e8-8029-db870ff11d9e',
	'367d95c0-0eb0-4dae-8276-9407239421ee_6acb6637-ac08-4a65-b2d1-581e51dc7ccf',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_df287f8d-f50d-4620-ab96-489d559e6adc',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_a810e511-c18b-4b2a-8fdf-98a6a0d433a7',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_88c483bf-477d-4be5-90d3-4fb101dd601f',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_8b2e5453-faf7-46ea-9073-aea69b283cb7',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_20634fa3-f3cf-44b5-8bc3-b825610bfe8c',
	'fbc5881f-1ee3-4ffe-8095-35e15e1a08fc_2aef80da-acb4-4e15-8f7d-6c0322b86b2f',
	'4b54248f-2165-477c-a027-dd55082e8818_43770b51-4b0e-4180-a83b-e77664c60736',
	'4b54248f-2165-477c-a027-dd55082e8818_5ba85070-a41c-4184-9c18-cf34c3fd0f62',
	'4b54248f-2165-477c-a027-dd55082e8818_d3a83885-5198-4b04-8314-b753b66ef9a8',
	'4b54248f-2165-477c-a027-dd55082e8818_dd018fc0-8da7-4033-a2ba-6b47de8ebb4f',
	'4b54248f-2165-477c-a027-dd55082e8818_03c0e874-f984-4e6c-9d2a-26ef8507dbbc',
	'4b54248f-2165-477c-a027-dd55082e8818_2a262b59-7936-4ecd-b656-248247a0559f',
	'4b54248f-2165-477c-a027-dd55082e8818_574e9f9e-f8b4-41ef-bf19-89a9964fd9c7',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_1009f384-b12d-448e-ba9f-1b7d2ecfbb4e',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_f75f2ff4-2884-4c2d-b375-70de37a34507',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_84f1a631-910b-4fbb-9f76-d915a07316d2',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_ed852810-a003-4386-9846-1638362cee39',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_9d584fcb-a28a-4b91-a886-ceb66a88ef81',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_572f3f3e-d3e4-4d13-8e2b-88215e508481',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_78fd69d2-75e4-4207-819a-563139f273c6',
	'b52eb423-5d0d-4645-b217-e1c6d38b2e72_d4e69e01-3ba2-4d6b-a15d-e7048f78f22e',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_bdacc907-7c26-419f-8808-969eab3ca2e8',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_9f1049ac-f8b7-45ad-8e31-6e96c3e5058f',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_06b91002-4d3d-4d2e-8484-20c3b31e232c',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_b94e3bdf-a385-49cc-b312-7a63cc28b77a',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_f9ad5649-f372-43e1-a3a8-423383e5a8a2',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_75a881cf-5d88-46e2-bf9b-97e5cbc1bd56',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_24066994-8183-488d-b037-ef6bb524af39',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_1492eb6b-7d50-4c4d-94ac-c801a7d5555c',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_cd77258f-b08b-4c89-b93f-6e6f146b1a4d',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_873ff933-4fda-4936-9a70-67df11af90ac',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_2727d83a-0af0-443a-bff8-58dc7028289a',
	'180bff9c-c8a5-4539-b13b-ddbc00d643e6_6c600df6-ddca-4628-a8bb-1d6de1e3f9b4',
	'60358420-6055-411d-ba4f-e8ac80682a2e_08e94873-c2a6-4f7d-ab72-aeaff3e3f929',
	'60358420-6055-411d-ba4f-e8ac80682a2e_aa0b5adb-957d-4f15-ab83-2c5cc2843f77',
	'60358420-6055-411d-ba4f-e8ac80682a2e_4269074c-f2c1-4d88-b2c3-0946f59d5449',
	'60358420-6055-411d-ba4f-e8ac80682a2e_b9b4cf27-9c22-410d-8bd8-5d43e379485b',
	'60358420-6055-411d-ba4f-e8ac80682a2e_2d66790a-6621-4a49-8f0d-4002db5cc98d',
	'60358420-6055-411d-ba4f-e8ac80682a2e_58679288-9ecc-4647-9781-12a3a8f8c6fd',
	'60358420-6055-411d-ba4f-e8ac80682a2e_774de9c6-9752-4e39-89a9-2a88c869d52a',
	'60358420-6055-411d-ba4f-e8ac80682a2e_abd889c6-f60a-4fbd-924e-ee1e9dcf909b',
	'60358420-6055-411d-ba4f-e8ac80682a2e_04b0eb97-d816-44bb-93a5-8b2968791aa0',
	'60358420-6055-411d-ba4f-e8ac80682a2e_bbd16004-09e8-4b6c-b465-73ff83a52837',
	'60358420-6055-411d-ba4f-e8ac80682a2e_e006d4e3-35fa-44b4-9981-09a66c4322e5',
	'60358420-6055-411d-ba4f-e8ac80682a2e_c42c8ad3-9761-49e5-b9bf-ee8ebd50416f',
	'60358420-6055-411d-ba4f-e8ac80682a2e_4506d9e3-4543-4464-aeae-b0b04eee1cea',
	'60358420-6055-411d-ba4f-e8ac80682a2e_9dfd2243-74d6-4924-86bd-c206ca9287b1',
	'60358420-6055-411d-ba4f-e8ac80682a2e_4d2e0563-cf4a-48bd-aa7f-efc26025b53a',
	'60358420-6055-411d-ba4f-e8ac80682a2e_fd89be61-2869-4342-a86e-e1fce3a8f269',
	'60358420-6055-411d-ba4f-e8ac80682a2e_9d5df009-eb76-43a3-b6cd-22017cc53700'
	]

# get the specified mapping strategy/file for each dataset
sheet_id = '18e5PG2wCaN8kf9-KVm_yomgEx8TYka0Ldd7_swVxiJk'
sheet_name = 'datasets'
url = 'https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}'.format(sheet_id, sheet_name)
ds_df = pd.read_csv(url)[['coll_ds','var_mapping']]
ds_df = ds_df.loc[ds_df['var_mapping'] != 'Lattice dataset'].dropna()


client = boto3.client('s3')
resource = boto3.resource('s3')
bucket_name = 'submissions-lattice'
your_bucket = resource.Bucket(bucket_name)

already_run = []
for s3_file in your_bucket.objects.all():
	if 'cxg_migration/final' in s3_file.key:
		already_run.append(s3_file.key.split('/')[-1].split('.')[0])

in_working = []
for s3_file in your_bucket.objects.all():
	if 'cxg_migration/working' in s3_file.key:
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
	if ds not in attn_needed and ds not in already_run and ds in in_working:
		print('PROCESSING:' + ds)
		file = ds + '.h5ad'
		client.download_file('submissions-lattice', 'cxg_migration/working/' + file, file)
		main(ds, row['var_mapping'])
		validate_cxg(ds)
		client.upload_file(file, bucket_name, 'cxg_migration/final/' + file, ExtraArgs={'ACL':'public-read'})
		os.remove(file)
