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

	obs = obs.rename(columns=guide.get('prop_rename', {}))

	prop_map = {
		'sex': {
			'update_in': 'sex_ontology_term_id',
			'value_map': {
				'female': 'PATO:0000383',
				'male': 'PATO:0000384',
				'mixed': 'unknown',
				'unknown': 'unknown'
			}
		}
	}
	prop_map.update(guide.get('prop_map', {}))
	for k,v in prop_map.items():
		old = k
		new = v.get('update_in', k)
		if new in obs.columns:
			for k2, v2 in v['value_map'].items():
				obs[new] = np.where((obs[old] == k2), v2,obs[new])
		else:
			data = {old: [], new: []}
			for k2, v2 in v['value_map'].items():
				data[old].append(k2)
				data[new].append(v2)
			map_df = pd.DataFrame.from_dict(data)
			obs = obs.merge(map_df, on=old)

	obs['is_primary_data'] = guide.get('is_primary_data', True)

	if 'organism_ontology_term_id' not in obs.columns:
		obs['organism_ontology_term_id'] = uns['organism_ontology_term_id']

	print(obs.dtypes)
	obs = obs.astype(guide.get('update_dtypes', {}))
	print(obs.dtypes)
	quit()

	uns['schema_version'] = '2.0.0'
	uns.update(guide.get('add_uns', {}))

	# UPDATE OBSM

	if guide.get('drop_raw_x') == True:
		del adata.raw

	# FILL IN NP.NAN?

	# REMOVE FIELDS FROM VAR (portal_fields = ['feature_name'])
	# OR DEAL WITH DURING GENE MIGRATION?

	# remove fields from uns
	deprecated_uns = ['organism', 'organism_ontology_term_id', 'default_field', 'tags', 'version']
	remove_uns = deprecated_uns + guide.get('remove_uns', [])
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
		if k.endswith('_original') or k in portal_fields + guide.get('remove_obs', []):
			remove_obs.append(k)
	obs = obs.drop(columns=remove_obs)

	# SWAP X AND RAW.X
	#if ds_guide.get('swap_layers'):

	adata.obs = obs
	adata.uns = uns
	adata.write(filename=ds + '.h5ad')

s3 = boto3.client('s3')
guides = os.listdir('guides')
for g in guides:
	ds = g.split('.')[0]
	file = ds + '.h5ad'
	s3.download_file('submissions-lattice', 'cxg_migration/original/' + file, file)
	main(ds)
	s3.upload_file(file, 'submissions-lattice', 'cxg_migration/working/' + file)
	os.remove(file)
