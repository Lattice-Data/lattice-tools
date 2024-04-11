import numpy as np



def colors_check(adata, color_column, column_name):
	'''
	Check validity of colors before adding to cxg_adata.uns
	Check that obs column exists
	Check that the corresponding column is the right datatype
	Verify color_column is a numpy array
	Verify that the numpy array contains strings
	Verify that we have atleast as many colors as unique values in the obs column
	Verify that either all colors are hex OR all colors are CSS4 named colors strings
	'''
	column_name = column_name.replace('_colors','')
	if column_name not in adata.obs.columns:
		error = 'the corresponding column is not present in obs.'
		return False, error
	if column_name in adata.obs.columns:
		if adata.obs[column_name].dtype.name != 'category':
			error = 'the corresponding column in obs. is the wrong datatype ({})'.format(adata.obs[column_name].dtype.name)
			return False, error
	if color_column is None or not isinstance(color_column, np.ndarray):
		error = 'the column is not a numpy array.'
		return False, error
	if not all(isinstance(color,str) for color in color_column):
		error = 'the column does not contain strings.'
		return False, error
	if len(color_column) < len(adata.obs[column_name].unique()):
		error = 'the column has less colors than unique values in the corresponding obs. column.'
		return False, error
	all_hex_colors = all(re.match(r"^#([0-9a-fA-F]{6})$", color) for color in color_column)
	all_css4_colors = all(color in mcolors.CSS4_COLORS for color in color_column)
	if not (all_hex_colors or all_css4_colors):
		error = 'the colors are not all hex or CSS4 named color strings.'
		return False, error
	else:
		error = 'none'
		return True, error



def check_not_empty(value):
	'''
	Return False if value is considered empty
	'''
	if any([
		isinstance(value, sparse_class)
		for sparse_class in (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix)
	]):
		if value.nnz == 0:
			return False
	elif (
		value is not None
		and not isinstance(value, numbers.Number)
		and type(value) is not bool
		and not (isinstance(value, (np.bool_, np.bool)))
		and len(value) == 0
	):
		return False
	else:
		return True



def copy_over_uns(reserved_uns, mfinal_adata, cxg_adata, warning_list)
'''
Copy over any additional data from mfinal_adata to cxg_adata
'''
	for k,v in mfinal_adata.uns.items():
		if k == 'batch_condition':
			if not isinstance(mfinal_adata.uns['batch_condition'], list) and not isinstance(mfinal_adata.uns['batch_condition'], np.ndarray) :
				warning_list.append("WARNING: adata.uns['batch_condition'] is not a list and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
			else:
				if len([x for x in mfinal_adata.uns['batch_condition'] if x not in cxg_adata.obs.columns]) > 0:
					warning_list.append("WARNING: adata.uns['batch_condition'] contains column names not found and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
				elif len(set(mfinal_adata.uns['batch_condition'])) != len(mfinal_adata.uns['batch_condition']):
					warning_list.append("WARNING: adata.uns['batch_condition'] contains redundant column names and did not get copied over to flattened h5ad: {}".format(mfinal_adata.uns['batch_condition']))
				else:
					cxg_adata.uns['batch_condition'] = mfinal_adata.uns['batch_condition']
		elif k.endswith('_colors'):
			colors_result = colors_check(cxg_adata, v, k)
			if colors_result[0]:
				cxg_adata.uns[k] = v
			else:
				warning_list.append("WARNING: '{}' has been dropped from uns dict due to being invalid because '{}' \n".format(k,colors_result[1]))
		elif k not in reserved_uns:
			if check_not_empty(v):
				cxg_adata.uns[k] = v
			else:
				warning_list.append("WARNING: The key '{}' has been dropped from uns due to having an empty value\n".format(k))
		else:
			warning_list.append("WARNING: The key '{}' has been dropped from uns dict due to being reserved \n".format(k))
	return cxg_adata, warning_list

