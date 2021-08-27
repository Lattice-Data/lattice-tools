import boto3
import json
import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def main(ds):
	guide = json.load(open('guides/{}.json'.format(ds)))
	adata = sc.read_h5ad(ds + '.h5ad')
	obs = adata.obs
	uns = adata.uns
	obs_len = len(obs)

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
		old_dtype = obs.dtypes[old]
		if new in obs.columns:
			for k2, v2 in v['value_map'].items():
				if old_dtype == '<f8':
					k2 = float(k2)
				elif old_dtype == 'int64':
					k2 = int(k2)
				obs[new] = np.where((obs[old] == k2), v2,obs[new])
		else:
			data = {old: [], new: []}
			for k2, v2 in v['value_map'].items():
				if old_dtype == '<f8':
					k2 = float(k2)
				elif old_dtype == 'int64':
					k2 = int(k2)
				data[old].append(k2)
				data[new].append(v2)
			map_df = pd.DataFrame.from_dict(data)
			obs = obs.merge(map_df, on=old)

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
		if uns.get(e):
			del uns[e]

	# remove columns from obs
	portal_fields = ['assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'sex', 'tissue', 'organism']
	remove_obs = []
	for k in obs.keys():
		if k in portal_fields + guide['remove_obs']:
			remove_obs.append(k)
	obs = obs.drop(columns=remove_obs)

	# swap raw.X and .X
	if guide.get('swap_layers'):
		orig_raw = adata.raw
		adata.raw = adata
		adata.X = orig_raw.X

	# write the new object to the file
	adata.obs = obs
	adata.uns = uns
	adata.write(filename=ds + '.h5ad')

attn_needed = [
	'7edef704-f63a-462c-8636-4bc86a9472bd_b83559d1-156f-4ba9-9f6a-b165f83ef43f',
	'f70ebd97-b3bc-44fe-849d-c18e08fe773d_e0ed3c55-aff6-4bb7-b6ff-98a2d90b890c'
	]

s3 = boto3.client('s3')
guides = os.listdir('guides')
for g in guides:
	ds = g.split('.')[0]
	if ds not in attn_needed:
		file = ds + '.h5ad'
		s3.download_file('submissions-lattice', 'cxg_migration/original/' + file, file)
		main(ds)
		s3.upload_file(file, 'submissions-lattice', 'cxg_migration/working/' + file)
		os.remove(file)
