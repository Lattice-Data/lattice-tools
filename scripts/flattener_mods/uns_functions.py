import sys
import pandas as pd
import numpy as np
from scipy import sparse
import os
import numbers
import logging
import re
import matplotlib.colors as mcolors
import flattener_mods.constants as constants
import flattener_mods as fm
from PIL import Image
Image.MAX_IMAGE_PIXELS = 933120000


# Attaching logger to Flattener logger
logger = logging.getLogger(__name__)


def colors_check(glob, color_column, column_name):
	'''
	Check validity of colors before adding to cxg_uns

	:param Class glob: Class object containing certain commonly used variables
	:param df[str] color_column: Anndata dataframe column containing color code strings
	:param str column_name: Name of obs column that color_column indicates colors for

	:returns Boolean, str: Returns whether colors are correct, and if not then what the issue is.
	'''
	column_name = column_name.replace('_colors', '')
	# Check that obs column exists
	if column_name not in glob.cxg_obs.columns:
		error = 'the corresponding column is not present in obs.'
		return False, error
	# Check that the corresponding column is the right datatype
	if column_name in glob.cxg_obs.columns:
		if glob.cxg_adata.obs[column_name].dtype.name != 'category':
			error = 'the corresponding column in obs. is the wrong datatype ({})'.format(glob.cxg_adata.obs[column_name].dtype.name)
			return False, error
	# Verify color_column is a numpy array
	if color_column is None or not isinstance(color_column, np.ndarray):
		error = 'the column is not a numpy array.'
		return False, error
	# Verify that the numpy array contains strings
	if not all(isinstance(color, str) for color in color_column):
		error = 'the column does not contain strings.'
		return False, error
	# Verify that we have atleast as many colors as unique values in the obs column
	if len(color_column) < len(glob.cxg_obs[column_name].unique()):
		error = 'the column has less colors than unique values in the corresponding obs. column.'
		return False, error
	# Verify that either all colors are hex OR all colors are CSS4 named colors strings
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

	:param obj value: Dictionary key object that is being checked

	:return Boolean: Whether or not the value is empty

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
		and not (isinstance(value, (np.bool_, bool)))
		and len(value) == 0
	):
		return False
	else:
		return True

def copy_over_uns(glob, reserved_uns):
	'''
	Copy over uns information from mfinal_adata to cxg_uns

	:param Class glob: Class object containing certain commonly used variables
	:param List[str] reserved_uns: List containing uns column names that are reserved by schema

	:returns List[str] warnings: Any warnings that occur during the copy over process are returned to append to total warning list
	'''
	warnings = []
	for k,v in glob.mfinal_adata.uns.items():
		if k == 'batch_condition':
			if not isinstance(glob.mfinal_adata.uns['batch_condition'], list) and not isinstance(glob.mfinal_adata.uns['batch_condition'], np.ndarray) :
				warnings.append("WARNING: adata.uns['batch_condition'] is not a list and did not get copied over to flattened h5ad: {}".format(glob.mfinal_adata.uns['batch_condition']))
			else:
				if len([x for x in glob.mfinal_adata.uns['batch_condition'] if x not in glob.cxg_obs.columns]) > 0:
					warnings.append("WARNING: adata.uns['batch_condition'] contains column names not found and did not get copied over to flattened h5ad: {}".format(glob.mfinal_adata.uns['batch_condition']))
				elif len(set(glob.mfinal_adata.uns['batch_condition'])) != len(glob.mfinal_adata.uns['batch_condition']):
					warnings.append("WARNING: adata.uns['batch_condition'] contains redundant column names and did not get copied over to flattened h5ad: {}".format(glob.mfinal_adata.uns['batch_condition']))
				else:
					glob.cxg_uns['batch_condition'] = glob.mfinal_adata.uns['batch_condition']
		elif k.endswith('_colors'):
			colors_result = colors_check(glob, v, k)
			if colors_result[0]:
				glob.cxg_uns[k] = v
			else:
				warnings.append("WARNING: '{}' has been dropped from uns dict due to being invalid because '{}' \n".format(k, colors_result[1]))
		elif k not in reserved_uns:
			if check_not_empty(v):
				glob.cxg_uns[k] = v
			else:
				warnings.append("WARNING: The key '{}' has been dropped from uns due to having an empty value\n".format(k))
		else:
			warnings.append("WARNING: The key '{}' has been dropped from uns dict due to being reserved \n".format(k))
	return warnings



def process_spatial(glob):
	'''
	If spatial transcriptomics, add adata.uns['spatial'] fields and eventually add full res

	:param Class glob: Class object containing certain commonly used variables

	:returns List[str] warnings: Any warnings that occur during the process are returned to append to total warning list
	'''
	warnings = []
	if len(glob.mfinal_obj.get('libraries'))==1:
		if glob.cxg_obs['assay_ontology_term_id'].unique()[0] == 'EFO:0010961':
			glob.cxg_uns['spatial'] = glob.cxg_adata_raw.uns['spatial']
			glob.cxg_uns['spatial']['is_single'] = True
			spatial_lib = list(glob.cxg_uns['spatial'].keys())[0]
			# Moving spacial metadata from cxg_uns['spatial']
			if 'metadata' in glob.cxg_uns['spatial'][spatial_lib]:
				glob.cxg_uns['spatial_metadata'] = glob.cxg_uns['spatial'][spatial_lib]['metadata']
				del glob.cxg_uns['spatial'][spatial_lib]['metadata']
			# Deleting unwanted spacial information from spatial
			for i in list(glob.cxg_uns['spatial'][spatial_lib]['images']):
				if i not in ['hires','fullres']:
					del glob.cxg_uns['spatial'][spatial_lib]['images'][i]
			for i in list(glob.cxg_uns['spatial'][spatial_lib]['scalefactors']):
				if i not in ['spot_diameter_fullres','tissue_hires_scalef']:
					del glob.cxg_uns['spatial'][spatial_lib]['scalefactors'][i]

			if glob.mfinal_obj.get('fullres_s3_uri', None):
				filename = glob.mfinal_obj.get('fullres_s3_uri').split('/')[-1]
				if os.path.exists(constants.MTX_DIR+"/"+filename):
					print("{} was found locally".format(filename))
				else:
					fm.download_file(glob.mfinal_obj.get('fullres_s3_uri'), constants.MTX_DIR)
				if filename.endswith(('tif', 'tiff', 'jpg')):
					fullres_np = np.asarray(Image.open(constants.MTX_DIR+"/"+filename))
					glob.cxg_uns['spatial'][spatial_lib]['images']['fullres'] = fullres_np
				else:
					warnings.append("WARNING: Did not recognize fullres file format:\t{}".format(glob.mfinal_obj.get('fullres_s3_uri')))
		else:
			glob.cxg_uns['spatial'] = {}
			glob.cxg_uns['spatial']['is_single'] = True
			if 'X_spatial' not in glob.cxg_obsm:
				logging.error('ERROR: X_spatial embedding is required for Slide-seqV2')
				sys.exit('ERROR: X_spatial embedding is required for Slide-seqV2')
			else:
				glob.cxg_obsm['spatial'] = glob.cxg_obsm['X_spatial']
				del glob.cxg_obsm['X_spatial']
	else:
		glob.cxg_uns['spatial'] = {}
		glob.cxg_uns['spatial']['is_single'] = False
	return warnings

