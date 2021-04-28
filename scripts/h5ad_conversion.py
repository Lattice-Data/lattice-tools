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
	parser.add_argument('--remove', '-r', help='Remove the line corresponding to this index from matrix', type=int, default=None)
	parser.add_argument('--cell', '-c', help='The file containing cell metadata, assuming the same delimiter used.', type=str)
	parser.add_argument('--gtf', '-g', help='The gtf file that contains gene metadata', type=str)
	parser.add_argument('--outdir', '-o', help='The out directory where h5ad files will be stored.', type=str, required=True)
	parser.add_argument('--transpose', '-t', help='True/False if matrix is cell x gene', type=bool, required=True)
	parser.add_argument('--genome_version', '-v', help='Version of genome used in alignment: default GRCh38.', type=str, default='GRCh38')
	parser.add_argument('--feature_type', '-f', help='feature_types for expression matrix: default Gene Expression', type=str, default='Gene Expression')
	args = parser.parse_args()
	return args


def calc_gene(args):
	df_gtf = pd.DataFrame()
	ensembl = ""
	symbol = ""
	with gzip.open(args.gtf, mode='rt') as file:
		for line in file:
			line_strip = line.strip()
			if re.search(r"\tgene\t",line_strip) and not re.search(r"_PAR_Y",line_strip):
				line_info = line_strip.split("\t")
				gene_info = line_info[8].split(";")
				ensembl = re.search(r'(ENSG\d+)', gene_info[0])[1]
				symbol = re.search(r'"(.*)"', gene_info[2])[1]
				df_gtf = df_gtf.append(pd.DataFrame({'symbol': symbol, 'feature_types': args.feature_type, 'genome': args.genome_version}, index=[ensembl]))
	return df_gtf


def add_gene(adata, args, df_gtf):
	adata.var = pd.merge(adata.var, df_gtf, left_index=True, right_index=True, how='left')
	adata.var['gene_ids'] = adata.var.index
	adata.var = adata.var.set_index('symbol', drop=True)
	col_order = adata.var.columns.tolist()
	col_order = col_order[-1:] + col_order[:-1]
	adata.var = adata.var[col_order]
	adata.var.index.name = None


def main(args):
	all_files = os.listdir(args.indir)
	os.mkdir(temp_dir)
	if args.gtf is not None:
		df_gtf = calc_gene(args)
	for filename in all_files:
		if filename.endswith('gz'):
			with gzip.open(os.path.join(args.indir, filename), mode='rb') as f_in:
				with open(os.path.join(temp_dir, filename.replace('.gz', '')), 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			adata = ad.read_csv(os.path.join(temp_dir, filename.replace('.gz', '')), delimiter = args.delimiter)
		else:
			adata = ad.read_csv(os.path.join(args.indir, filename), delimiter = args.delimiter)

		if args.transpose:
			adata = adata.transpose()
		if args.gtf is not None:
			add_gene(adata, args, df_gtf)
		if type(adata.X) == np.ndarray:
			adata.X = sparse.csr_matrix(adata.X)
		adata.write(os.path.join(args.outdir, filename.replace('.txt.gz', '.h5ad')))

	shutil.rmtree(temp_dir)


args = getArgs()

if __name__ == '__main__':
	main(args)
