import json
import pandas as pd


sheet_id = '18e5PG2wCaN8kf9-KVm_yomgEx8TYka0Ldd7_swVxiJk'

sheet_name = 'datasets'
url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'
ds_df = pd.read_csv(url).fillna('') #fillna needed before will fill in values
ds_df = ds_df[['coll_ds','X_normalization','is_primary_data']]
ds_df = ds_df.loc[ds_df['X_normalization'] != 'Lattice dataset']

guides = {}
no_pri = 0
no_nor = 0
for i,r in ds_df.iterrows():
	if r['is_primary_data'] not in [True, False]:
		no_pri += 1
	if not r['X_normalization']:
		no_nor += 1
	guides[r['coll_ds']] = {
		'is_primary_data': r['is_primary_data'],
		'prop_rename': {},
		'prop_map': {},
		'update_dtypes': {},
		'add_uns': {
			'X_normalization': r['X_normalization']
		},
		'remove_obs': [],
		'remove_uns': []
	}
print(str(no_pri) + ' datasets without is_primary_data')
print(str(no_nor) + ' datasets without X_normalization')

v = lambda x: x.split(' [')[0]

# Need to convert numbers to string - Select 'notes' column, Format-->Number-->Plain Text
# 8 'notes' with FALSE/TRUE, replace with false/true
sheet_name = 'errors'
url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'
err_df = pd.read_csv(url, converters={'category': v, 'note': v, 'update': v}).fillna('')

# check for rows with value in 'update' AND 'reason dismissed'
both_input = len(err_df.loc[(err_df['update'] != '') & (err_df['reason dismissed'] != '')])
print(str(both_input) + ' rows with input in both update & dismissed')

# check for rows without value in either
no_input = len(err_df.loc[(err_df['update'] == '') & (err_df['reason dismissed'] == '')])
print(str(no_input) + ' rows with input in neither update & dismissed')

# filter down to only rows with specified updates
err_df = err_df.loc[err_df['update'] != '']

for i, r in err_df.iterrows():
	ds = r['coll_ds']
	cat = r['category']
	up = r['update']
	field = r['impacted field']
	if cat in ['remove_field', 'update_field'] and up == 'remove_field':
		if field not in guides[ds]['remove_obs']:
			guides[ds]['remove_obs'].append(field)
	elif cat == 'update_field' and up == 'strip_whitespaces':
		if field in guides[ds]['prop_rename']:
			print(str(i) + ':MULTIPLE RENAME FIELD')
		guides[ds]['prop_rename'][field] = ' '.join(field.split())
	elif cat == 'update_field':
		if field in guides[ds]['prop_rename']:
			print(str(i) + ':MULTIPLE RENAME FIELD')
		guides[ds]['prop_rename'][field] = up
	elif cat == 'data' and up == 'remove raw.X':
		guides[ds]['drop_raw_x'] = True
	elif cat == 'data' and up == 'swap raw.X and .X':
		guides[ds]['swap_layers'] = True
	elif cat == 'update_value':
		if field == 'is_primary_data':
			up = up.lower() == 'true'
		if r['impacted field'] not in guides[ds]['prop_map']:
			prop_map = {'value_map': {r['note']: up}}
			if r['map-from field']:
				prop_map['update_from'] = r['map-from field']
			guides[ds]['prop_map'][field] = prop_map
		else:
			if r['map-from field'] and guides[ds]['prop_map'][field].get('update_from') != r['map-from field']:
				print(str(i) + ':VARIABLE MAP-FROM FIELD')
			if r['note'] in guides[ds]['prop_map'][field]['value_map']:
				print(str(i) + ':REDUNDANT VALUE_MAP')
			guides[ds]['prop_map'][field]['value_map'][r['note']] = up
	elif cat == 'add_uns_field' and up == 'add_field':
		if field in guides[ds]['add_uns']:
			print(str(i) + ':MULTIPLE UNS ADD')
		guides[ds]['add_uns'][field] = r['note']
	elif cat == 'add_uns_field_list' and up == 'add_field':
		if field in guides[ds]['add_uns']:
			print(str(i) + ':MULTIPLE UNS ADD')
		guides[ds]['add_uns'][field] = r['note'].split(',')
	elif cat == 'remove_uns_field' and up == 'remove_field':
		if field not in guides[ds]['remove_uns']:
			guides[ds]['remove_uns'].append(field)
	elif cat == 'update_dtype':
		if field in guides[ds]['update_dtypes']:
			print(str(i) + ':MULTIPLE DTYPE UPDATE')
		guides[ds]['update_dtypes'][field] = up
	else:
		print(str(i) + ':CATEGORY NOT MAPPED TO GUIDE')

for ds,g in guides.items():
	with open('guides/{}.json'.format(ds), 'w') as f:
		json.dump(g, f, indent=4)
