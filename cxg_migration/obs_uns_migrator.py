import boto3
import json
import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def main(ds):
	guide = json.load(open('guides/{}.json'.format(ds)))

	adata = sc.read_h5ad(ds + '.h5ad')

	obs = adata.obs
	uns = adata.uns

	remove_obs = guide['remove_obs']

	# rename obs fields
	obs = obs.rename(columns=guide['prop_rename'])

	# map values to new properties and within existing properties
	prop_map = guide['prop_map']
	if 'sex_ontology_term_id' not in prop_map and 'sex_ontology_term_id' not in obs:
		prop_map.update({
			'sex_ontology_term_id': {
				'update_from': 'sex',
				'value_map': {
					'female': 'PATO:0000383',
					'male': 'PATO:0000384',
					'mixed': 'unknown',
					'unknown': 'unknown'
				}
			}
		})
	for k,v in prop_map.items():
		old = v.get('update_from', k)
		new = k
		old_dtype = obs[old].dtype
		if v.get('update_from'):
			map_df = pd.DataFrame(columns=[old, new])
			for k2, v2 in v['value_map'].items():
				if old_dtype in ['<f8','float64']:
					k2 = float(k2)
				elif old_dtype in ['int64','int32']:
					k2 = int(k2)
				elif old_dtype != 'object':
					if old_dtype == 'category':
						if obs[old].cat.categories.dtype in ['int64','int32']:
							k2 = int(k2)
				map_df.loc[map_df.shape[0]] = [k2, v2]
			if new in obs.columns:
				map_df = map_df.append(obs[obs[old].isin(map_df[old]) == False][[old,new]].drop_duplicates())
			obs = obs.reset_index(drop=True).merge(map_df,how='left',on=old,suffixes=('_x',None)).set_index(obs.index)
			obs[old] = obs[old].astype(old_dtype)
			if new == 'is_primary_data':
				obs[new] = obs[new].astype('bool')
			remove_obs.append(new + '_x')
		else:
			obs[k].cat = obs[k].cat.add_categories(list(set(v['value_map'].values())))
			obs[k] = obs[k].replace(v['value_map'])
			obs[k] = obs[k].astype('category')

	# fill in is_primary_data
	if 'is_primary_data' not in obs:
		obs['is_primary_data'] = guide['is_primary_data']

	# migration organism from uns to obs if not already present
	if 'organism_ontology_term_id' not in obs.columns:
		obs['organism_ontology_term_id'] = uns['organism_ontology_term_id']

	# update dtypes
	obs = obs.astype(guide['update_dtypes'])

	# add required schema_version and other specified uns fields
	uns['schema_version'] = '2.0.0'
	uns.update(guide['add_uns'])

	# drop raw.X if specified (if it is duplicate of .X)
	if guide.get('drop_raw_x') == True:
		del adata.raw

	# remove fields from uns
	deprecated_uns = ['organism', 'organism_ontology_term_id', 'default_field', 'tags', 'version']
	remove_uns = deprecated_uns + guide['remove_uns']
	for k in uns.keys():
		if k.endswith('_original'):
			remove_uns.append(k)
	for e in remove_uns:
		if e in uns:
			del uns[e]

	# remove columns from obs
	portal_props = ['assay','tissue','cell_type','sex','development_stage','ethnicity','disease','organism']
	remove = []
	for k in obs.keys():
		if k in portal_props + remove_obs:
			remove.append(k)
	obs = obs.drop(columns=remove)

	# swap raw.X and .X
	if guide.get('swap_layers'):
		orig_raw = adata.raw
		adata.raw = adata
		adata.X = orig_raw.X

	if adata.raw:
		if type(adata.raw.X) != sparse.csr.csr_matrix:
			raw_adata = ad.AnnData(adata.raw.X, var=adata.raw.var, obs=adata.obs)
			raw_adata.X = sparse.csr_matrix(raw_adata.X)
			adata.raw = raw_adata
	if type(adata.X) != sparse.csr.csr_matrix:
		adata.X = sparse.csr_matrix(adata.X)
	for l in adata.layers:
		if type(adata.layers[l]) != sparse.csr.csr_matrix:
			adata.layers[l] = sparse.csr_matrix(adata.layers[l])
	
	# write the new object to the file
	adata.obs = obs
	adata.uns = uns
	adata.write(filename=ds + '.h5ad')

attn_needed = [
	'7edef704-f63a-462c-8636-4bc86a9472bd_b83559d1-156f-4ba9-9f6a-b165f83ef43f', # Voigt/Scheetz retina, no raw counts
	'a238e9fa-2bdf-41df-8522-69046f99baff_66d15835-5dc8-4e96-b0eb-f48971cb65e8' # Enge pancreas, cell don't group by cluster/cell_type
	]

client = boto3.client('s3')
guides = os.listdir('guides')

already_run = []
resource = boto3.resource('s3')
your_bucket = resource.Bucket('submissions-lattice')
for s3_file in your_bucket.objects.all():
	if 'cxg_migration/working' in s3_file.key:
		already_run.append(s3_file.key.split('/')[-1].split('.')[0])

for g in guides:
	ds = g.split('.')[0]
	if ds not in attn_needed and ds not in already_run:
		print('PROCESSING:' + ds)
		file = ds + '.h5ad'
		client.download_file('submissions-lattice', 'cxg_migration/original/' + file, file)
		main(ds)
		client.upload_file(file, 'submissions-lattice', 'cxg_migration/working/' + file)
		os.remove(file)
