import os
import io
import pandas as pd
import shutil
import sys
import scanpy as sc
import argparse
import anndata as ad
import gzip
import re
import numpy as np
from scipy import sparse


EPILOG = '''
Example:
	
		python %(prog)s -f GSM3988006_SAM24362284.txt -d tab

For more details:
		python %(prog)s --help

'''

temp_dir = 'unzipped_files'

def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('--indir', '-i', help='The directory of txt file that needs to be converted.', type=str, required=True)
	parser.add_argument('--delimiter', '-d', help="The delimiter used in the files being converted: , or '	'.", default='tab', type=str, required=True)
	parser.add_argument('--remove', '-r', help='Remove first observation entry', type=bool, default=False)
	parser.add_argument('--cell', '-c', help='The file containing cell metadata, assuming the same delimiter used.', type=str)
	parser.add_argument('--gtf', '-g', help='The gtf file that contains gene metadata', type=str, default=None)
	parser.add_argument('--symbol', '-n', help='If matrix is keyed off of gene symbol and need to map ensembl', type=bool, default=False)
	parser.add_argument('--outdir', '-o', help='The out directory where h5ad files will be stored.', type=str, required=True)
	parser.add_argument('--transpose', '-t', help='True/False if matrix is cell x gene', type=bool, required=True, default=False)
	parser.add_argument('--genome_version', '-v', help='Version of genome used in alignment: default GRCh38.', type=str, default='GRCh38')
	parser.add_argument('--feature_type', '-f', help='feature_types for expression matrix: default Gene Expression', type=str, default='Gene Expression')
	parser.add_argument('--layer', '-l', help='Directory coontaining scaled layer matrices', type=str, default=None)
	parser.add_argument('--beadloc', '-b', help='SlideSeq BeadLocations file directory', type=str, default=None)

	args = parser.parse_args()
	return args

# create dataframe of relevant gtf information
def calc_gene(args):
	df_gtf = pd.DataFrame()
	ensembl = ""
	symbol = ""
	#with gzip.open(args.gtf, mode='rt') as file:
	with open(args.gtf, mode='r') as file:
		for line in file:
			line_strip = line.strip()
			if args.genome_version == 'GRCH38':
				if re.search(r"\tgene\t",line_strip) and not re.search(r"_PAR_Y",line_strip):
					line_info = line_strip.split("\t")
					gene_info = line_info[8].split(";")
					ensembl = re.search(r'(ENSG\d+)', gene_info[0])[1]
					symbol = re.search(r'"(.*)"', gene_info[2])[1]
					df_gtf = df_gtf.append(pd.DataFrame({'symbol': symbol, 'feature_types': args.feature_type, 'genome': args.genome_version}, index=[ensembl]))
			else:
				if not re.search(r"pseudogene", line_strip) and not re.search(r"^H", line_strip):
					gene_info = line_strip.split(";")
					ensembl = re.search(r'(ENSG\d+)', gene_info[1])[1]
					symbol = re.search(r'"(.*)"', gene_info[2])[1]
					df_gtf = df_gtf.append(pd.DataFrame({'symbol': symbol}, index=[ensembl]))
	print(df_gtf)
	return df_gtf


# add gene information to anndata, including gene metadata, switch to symbol as index
def add_gene(adata, args, df_gtf):
	adata.var = pd.merge(adata.var, df_gtf, left_index=True, right_index=True, how='left')
	adata.var['gene_identifier'] = adata.var.index
	adata.var = adata.var.set_index('symbol', drop=True)
	col_order = adata.var.columns.to_list()
	col_order = col_order[-1:] + col_order[:-1]
	adata.var = adata.var[col_order]
	adata.var.index.name = None


# map ensembl identifier from gtf dataframe, errors when there are multiple mappings
def get_ensembl(adata, args, df_gtf):
	var = pd.merge(adata.var, df_gtf, left_index=True, right_on="symbol", how="left")
	if var.shape[0] == adata.var.shape[0]:
		adata.var['gene_identifier'] = var.index
		adata.var['feature_types'] = args.feature_type
		adata.var['genome'] = args.genome_version
	else:
		sys.exit("Cannot merge 2 matrices:{} {}".format(var.shape, adata.var.shape))
			

# add bead location as embedding and celltype in obs
def add_beadloc(filename, adata, args):
	beadfile = filename.replace("MappedDGEForR_sct_cts_extensive", "BeadLocationsForR")
	bead_df = pd.read_csv(os.path.join(args.beadloc, beadfile), index_col='barcodes')
	bead_df.index.name = None
	bead_df = bead_df.loc[adata.obs.index.to_list(),]
	bead_np = bead_df[['xcoord', 'ycoord']].to_numpy(copy=True)
	adata.obsm['X_spatial'] = bead_np
	adata.obs = pd.merge(adata.obs, bead_df[['cell_type']], left_index=True, right_index=True, how='left', validate="1:1")
	adata.obs = adata.obs.rename(columns={'cell_type': 'celltype'})


# add scaled layer to AnnData.layers
def add_layer(filename, adata, args):
	layerfile = filename.replace("cts_extensive", "scaled")
	layer_adata = ad.read_csv(os.path.join(args.layer, layerfile), delimiter = args.delimiter, first_column_names=True)
	layer_adata = layer_adata.transpose()
	layer_adata = layer_adata[adata.obs_names.to_list()]
	adata.layers['scaled'] = layer_adata.X


def main(args):
	all_files = os.listdir(args.indir)
	os.mkdir(temp_dir)
	if args.gtf:
		df_gtf = calc_gene(args)
	for filename in all_files:
		if filename.endswith('gz'):
			with gzip.open(os.path.join(args.indir, filename), mode='rb') as f_in:
				with open(os.path.join(temp_dir, filename.replace('.gz', '')), 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			adata = ad.read_csv(os.path.join(temp_dir, filename.replace('.gz', '')), delimiter = args.delimiter)
			print("converting {}".format(filename))
		elif filename.endswith('csv'):
			adata = ad.read_csv(os.path.join(args.indir, filename), delimiter = args.delimiter, first_column_names=True)
			print("converting {}".format(filename))
		else:
			continue

		if type(adata.X) == np.ndarray:
			adata.X = sparse.csr_matrix(adata.X)
		if args.transpose:
			print("transposed the data")
			adata = adata.transpose()
		if args.remove:
			if re.search(r"Unnamed", adata.var.index[0]):
				adata = adata[:,1:adata.n_vars]
		if args.gtf is not None:
			if args.symbol:
				get_ensembl(adata, args, df_gtf)
			else:
				add_gene(adata, args, df_gtf)
		if args.beadloc:
			add_beadloc(filename, adata, args)
		if args.layer:
			add_layer(filename, adata, args)
		adata.write(os.path.join(args.outdir, filename.replace('sct_cts_extensive.csv', 'final.h5ad')))

	shutil.rmtree(temp_dir)


args = getArgs()

if __name__ == '__main__':
	main(args)
