import argparse
import json
import os
import re
import sys
import tarfile
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from urllib.parse import ParseResult, urlparse
from urllib.request import Request, urlopen

import anndata as ad
import boto3
import fsspec
import numpy as np
import pandas as pd
import requests
import scanpy as sc
from bs4 import BeautifulSoup
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from cellxgene_schema.write_labels import AnnDataLabelAppender

sys.path.append(os.path.dirname(os.path.abspath('../cellxgene_resources')))
from cellxgene_resources.cellxgene_mods import map_filter_gene_ids

EPILOG = """
Script will take a single GroupID from a project on CZI S3 as input
and create curated matrices for all subsamples of this GroupID. This will need to be
run on JupyterHub to allow access to CZI S3 buckets and reference indices for Guidescan.
The current assumption is that this is run on 10x Flex data, and will need to
update anticipated directory structure if it is another 10x assay.

Can also pass txt file with args in the following format within the txt file:
One arg per line, with '=' between arg and value like follows:
--bucket=czi-psomagen
--project=marson-mapping-grns-perturb-seq

Example:
    python curate_matrices.py --bucket czi-psomagen --sheet your_sheet_id --project marson-mapping-grns-perturb-seq --group CD4i_R1L01 --csvofguidescan guidescan_out.csv
    python curate_matrices.py @args.txt

For more details:
    python %(prog)s --help

"""


TEMP_DIR = Path("temp_cellranger/")


@dataclass
class URIPath:
    """
    Dataclass to deal with S3 URIs and common metadata associated with them.
    full_uri should be in the format {s3://bucket/key/...}
    """
    full_uri: str
    local_dir: Path = TEMP_DIR

    def __post_init__(self):
        self._parsed: ParseResult = urlparse(self.full_uri)
        self.info: dict = FS.info(self.full_uri)
        self.type: str = self.info["type"]
        self.size: int = self.info["size"]

    @property
    def bucket(self) -> str:
        return self._parsed.netloc
    
    @property
    def key(self) -> str:
        return self._parsed.path.lstrip("/")
    
    @property
    def file_name(self) -> str | None:
        if self.type == "file":
            return self.full_uri.split("/")[-1]
        return None

    @property
    def local_path(self) -> Path:
        if self.file_name:
            return self.local_dir / self.file_name
        return self.local_dir


@dataclass
class LatticeMetadata:
    '''
    Dataclass to hold and map ontologies to a Lattice metadata spreadsheet.
    The functions allow for this object to read in a specific tab for a google sheet into a dataframe
    Subset should be a dict with key as column to subset, and value to filter column on
    '''
    
    sheet_id: str
    tab_name: str
    subset: dict[str, str] | None = None

    def get_gid(self):
        '''
        Given sheet id and tab name, return gid
        '''
        sheet_url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}'
        req = Request(sheet_url, headers={'User-Agent' : "Magic Browser"})
        s = urlopen(req)
        soup = BeautifulSoup(s, 'html.parser')
        tab_ids = {}
        pattern = re.compile('var bootstrapData = (.*?)};')
        for s in soup.find_all('script'):
            if pattern.search(str(s)):
                d = pattern.search(str(s)).group()[20:-1]
                data = json.loads(d)
                for t in data['changes']['topsnapshot']:
                    u = t[1].split('"')
                    if len(u) > 5:
                        tab_ids[u[5]] = u[1]
        return tab_ids[self.tab_name]

    def get_metadata_df(self):
        '''
        Given sheet id and gid, return lattice metadata in a dataframe
        Subset if subset dictionary is present
        '''
        url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}/export?format=csv&gid={self.gid}'
        response = requests.get(url)
        sample_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str).dropna(axis=1,how='all')

        if self.subset:
            sample_df = sample_df[sample_df[self.subset["column"]] == self.subset["filter_value"]]

        return sample_df

    def __post_init__(self):
        self.gid = self.get_gid()
        self.metadata_df = self.get_metadata_df()


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars="@",
    )
    parser.add_argument(
        "--bucket", 
        "-b",
        help="bucket where h5 files live",
        default="czi-psomagen",
        required=True
    )
    parser.add_argument(
        "--sheet", 
        "-s",
        help="google sheet id of metadata",
        required=True
    )
    parser.add_argument(
        "--project", 
        "-p",
        help="Name of the project directory on S3 bucket",
        required=True
    )
    parser.add_argument(
        "--group",
        "-g",
        help="the GroupID of the sample that for curated matrices generation",
        required=True
    )
    parser.add_argument(
        "--csvofguidescan",
        "-c",
        help="Guidescan output CSV",
        required=True
    )
    args = parser.parse_args()

    match sys.argv:
        case ["%(prog)s", txt_input] if txt_input.startswith("@"):
            print(f"Using args from {txt_input[1:]}")
        case ["%(prog)s"]: 
            parser.print_help()
            print("No arguments provided")
            sys.exit()
        case ["%(prog)s", *other_args] if len(other_args) > 10:   # 5 required args + 5 required values
            print("Too many arguments provided")
            parser.print_help()
            sys.exit()

    return args


def download_files(s_dir, bucket, lib_samp):
    """
    Download needed cellranger files from S3. For cellranger v10, there is no "count" subdirectory
    and crispr analysis is no longer tarred.
    """
    uri_stem = s_dir
    mx_h5_file = "sample_filtered_feature_bc_matrix.h5"
    cri_file = "crispr_analysis/protospacer_calls_per_cell.csv"
    metrics_csv = "metrics_summary.csv"

    if FS.isdir(f"{bucket}/{uri_stem}/count/"):
        uri_stem = f"{uri_stem}/count"
    if FS.isfile(f'{bucket}/{uri_stem}/crispr_analysis.tar.gz'):
        cri_file = "crispr_analysis.tar.gz"

    for file_path in [mx_h5_file, metrics_csv, cri_file]:
        f = file_path.split('/')[-1]
        file_path = f"{uri_stem}/{file_path}"
        full_path = os.path.join("temp_cellranger/",f"{lib_samp}_{f}")
        if not os.path.isfile(full_path):
            s3client.download_file(bucket, file_path, full_path)


def custom_var_to_obs(adata):
    moved = []
    for gene_index in adata.var[
        (adata.var['feature_types'] == 'Gene Expression') &
        (adata.var['gene_ids'].str.startswith('ENS') == False)
    ].index:
        gene_id = adata.var.loc[gene_index]['gene_ids']
        gene_values = adata.X[:, adata.var.index.get_loc(gene_index)]
        adata.obs[gene_id] = gene_values.A.flatten() if hasattr(gene_values, 'A') else gene_values.toarray().flatten().tolist()
        moved.append(gene_id)

    adata.var.set_index('gene_ids', inplace=True)
    var_to_keep = [i for i in adata.var.index if i not in moved]
    adata = adata[:, var_to_keep]
    adata.var = adata.var.replace('', np.nan).dropna(axis=1, how='all')

    return adata


def check_standard_presence(sample_df):
    '''
    Checks for all required columns in sample sheet dataframe

    :param dataframe sample_df: the sample metadata from given google sheet

    :returns None
    '''
    required_columns = [
    'sample_name','sample_probe_barcode','is_pilot_data',
    'donor_id', #'donor_living_at_sample_collection', 'donor_body_mass_index',
    'organism', 'sex', 'self_reported_ethnicity', 'disease', 'tissue', 'preservation_method',
    'development_stage', 'assay', 'tissue_type',
    'suspension_type']
    chk_exit = False
    for col in (col for col in required_columns if col not in sample_df.columns):
        if sample_df['organism'].unique()[0] == 'Homo sapiens':
            print(f"ERROR: Column '{col}' not present in sample sheet")
            chk_exit=True
        elif col not in ['donor_living_at_sample_collection', 'donor_body_mass_index']:
            print(f"ERROR: Column '{col}' not present in sample sheet")
            chk_exit=True
    if chk_exit:
        sys.exit()


def gather_crispr(samp):
    '''
    obtain guide assignement from 10x cellranger output files

    return df: dataframe containing crispr guide calling
    '''
    df = pd.DataFrame()
    if os.path.isfile(f'temp_cellranger/{samp}_protospacer_calls_per_cell.csv'):
        df = pd.read_csv(f'temp_cellranger/{samp}_protospacer_calls_per_cell.csv').rename(columns={'num_umis':'num_umis_guide_id'})
    else:
        with tarfile.open(f'temp_cellranger/{samp}_crispr_analysis.tar.gz', 'r:gz') as tar:
            f = tar.extractfile('protospacer_calls_per_cell.csv')
            df = pd.read_csv(f).rename(columns={'num_umis':'num_umis_guide_id'})
    df['genetic_perturbation_id'] = df['feature_call'].apply(lambda x: x.replace('|',' || '))
    df['num_umis_guide_id'] = df['num_umis_guide_id'].apply(lambda x: x.replace('|',' || '))
    df = df[['cell_barcode','genetic_perturbation_id','num_umis_guide_id']].set_index('cell_barcode')

    return df


def gather_metrics(samp, lib_samp):
    df = pd.read_csv(f'temp_cellranger/{lib_samp}_metrics_summary.csv')
    df['Metric'] = df.apply(lambda x: f"{x['Metric Name']}, {x['Library Type']}, {x['Category']}", axis=1)

    probe_barcodes = ' || '.join(df[
        (df['Metric Name'] == 'Sample ID') &
        (df['Metric Value'] == samp)
    ]['Group Name'].unique())

    df = df[
        (df['Grouped By'].isin(['Fastq ID','Probe barcode ID']) == False) &
        (df['Category'].isin(['Library','Cells']))
    ]
    df = df[['Metric','Metric Value']].set_index('Metric').transpose()

    for c in df.columns:
        v = df[c].iloc[0]
        if v.endswith('%'):
            df[c] = df[c].apply(lambda x: float(x.rstrip('%')) / 100)
        elif '%' in v:
            df.drop(columns=c, inplace=True)
        elif v.startswith('0.'):
            df[c] = df[c].apply(lambda x: float(x.replace(',','')))
        else:
            df[c] = df[c].apply(lambda x: int(x.replace(',','')))

    df['Probe barcode IDs'] = probe_barcodes

    keep = {
        'Confidently mapped reads in cells, Gene Expression, Cells': 'fraction_mapped_reads_in_cells_gex_sample',
        'Confidently mapped reads in cells, Gene Expression, Library': 'fraction_mapped_reads_in_cells_gex_library',
        'Median genes per cell, Gene Expression, Cells': 'median_genes_per_cell_sample',
        'Cells, CRISPR Guide Capture, Cells': 'cell_count_cri_sample',
        'Cells, Gene Expression, Cells': 'cell_count_gex_sample',
        'Cells, CRISPR Guide Capture, Library': 'cell_count_cri_library',
        'Cells, Gene Expression, Library': 'cell_count_gex_library',
        'Valid barcodes, CRISPR Guide Capture, Library': 'fraction_valid_barcodes_cri_library',
        'Valid barcodes, Gene Expression, Library': 'fraction_valid_barcodes_gex_library',
        'Reads mapped to probe set, Gene Expression, Cells': 'fraction_reads_mapped_gex_sample',
        'Reads mapped to probe set, Gene Expression, Library': 'fraction_reads_mapped_gex_library',
        'Mapped to genome, Gene Expression, Library': 'fraction_reads_mapped_gex_library',
        'Mean reads per cell, Gene Expression, Cells': 'mean_reads_per_cell_gex_sample',
        'Mean reads per cell, CRISPR Guide Capture, Library': 'mean_reads_per_cell_cri_library',
        'Mean reads per cell, Gene Expression, Library': 'mean_reads_per_cell_gex_library',
        'Sequencing saturation, CRISPR Guide Capture, Library': 'sequencing_saturation_cri_library',
        'Sequencing saturation, Gene Expression, Library': 'sequencing_saturation_gex_library',
        'Median UMI counts per cell, CRISPR Guide Capture, Cells': 'median_umi_counts_per_cell_cri_sample',
        'Median UMI counts per cell, Gene Expression, Cells': 'median_umi_counts_per_cell_gex_sample',
        'Probe barcode IDs': 'sample_probe_barcode'
    }

    df = df[[f for f in keep.keys() if f in df.columns]].rename(columns=keep)

    return df


def cxg_add_labels(adata):
    adata.obs['cell_type_ontology_term_id'] = 'unknown'
    labeler = AnnDataLabelAppender(adata)
    labeler._add_labels()
    adata.obs.drop(columns=['cell_type_ontology_term_id','cell_type'],inplace=True)

    schema_v = labeler.schema_version
    adata.uns['schema_version'] = schema_v
    adata.uns['schema_reference'] = labeler._build_schema_reference_url(schema_v)


def add_guide_metadata(adata, sheet, guide_gid, guidescan_output):
    '''
    Add guide metadata into adata.uns from guidescan_pipeline.py output

    :param obj adata: the anndata object that is being transformed into the curated matrix
    :param str guidescan_output: File containing output from guidescan_pipeline.py

    :returns obj adata: modified adata to contain guide metadata

    '''

    url = f'https://docs.google.com/spreadsheets/d/{sheet}/export?format=csv&gid={guide_gid}'
    response = requests.get(url)
    guide_sheet_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str)
    guide_df = pd.read_csv(guidescan_output, dtype=str)
    genetic_perturbations = {}

    for row in guide_df.itertuples():
        if row.id not in genetic_perturbations.keys():
            genetic_perturbations[row.id] = {}
            genetic_perturbations[row.id]['sequence'] = row.sequence
            if not pd.isna([row.start,row.end,row.sense]).all():
                chr_loc = str(row.chromosome).replace("chr","") + ":" + str(row.start) + "-" + str(row.end) + "(" + str(row.sense) + ")"
                genetic_perturbations[row.id]['derived_genomic_regions'] = [chr_loc]
        if not pd.isna(row.gene_id):
            if 'derived_features' not in genetic_perturbations[row.id].keys():
                genetic_perturbations[row.id]['derived_features'] = {}
            genetic_perturbations[row.id]['derived_features'][row.gene_id.split(".")[0]] = row.gene_name

    for row in guide_sheet_df.itertuples():
        genetic_perturbations[row.guide_id]['role'] = 'targeting' if row.guide_role == 'Targeting a Gene' else 'control'
        genetic_perturbations[row.guide_id]['protospacer_adjacent_motif'] = row.guide_PAM

    adata.uns['genetic_perturbations'] = genetic_perturbations
    return adata


def determine_perturbation_strategy(adata: ad.AnnData) -> ad.AnnData:
    '''
    Assess feature_call from protospacer_calls_per_cell.csv, where if all guides
    assigned to a single cell are all control, then 'control'. Otherwise, it is "no perturbations"
    if no guids or one of the following if targeting:
        - "CRISPR activation screen"
        - "CRISPR interference screen"
        - "CRISPR knockout mutant"
        - "CRISPR knockout screen"
    
    :param obj adata: the anndata object that is being transformed into the curated matrix

    :returns obj adata: modified adata to contain perturbation_strategy as cell metadata
    '''
    calculated_col = "strategy_calculated"
    adata.obs[calculated_col] = adata.obs['genetic_perturbation_id']

    lambdas = [
        lambda x: x.split(' || ') if pd.notna(x) else 'no perturbations',
        lambda x: [adata.uns['genetic_perturbations'][i]['role'] for i in x] if isinstance(x, list) else x,
        lambda x: 'control' if isinstance(x, list) and 'targeting' not in set(x) else x
    ]

    for rule in lambdas:
        adata.obs[calculated_col] = adata.obs[calculated_col].apply(rule)

    for value in ["control", "no perturbations"]:
        adata.obs.loc[adata.obs[calculated_col]==value, 'genetic_perturbation_strategy'] = value

    adata.obs.drop(columns=[calculated_col], inplace=True)
    
    return adata


def map_ontologies(sample_df):
    '''
    Takes the sample metadata dataframe and standardizes ontologies
    Also checks that standard fields are only filled out for appropriate organism

    :param dataframe sample_df: the sample metadata from given google sheet

    :returns dataframe sample_df: sample metadata with ontologies added
    '''
    col_ont_map = {
        'organism':'NCBITaxon',
        'sex':'PATO',
        'self_reported_ethnicity':{'NCBITaxon:9606':'HANCESTRO',
                                   'other':'none'},
        'disease':'MONDO',
        'assay':'EFO',
        'development_stage':{'NCBITaxon:6239':'WBls', # C. Elegans
                             'NCBITaxon:7227':'FBdv', # Drosophila
                             'NCBITaxon:10090':'MmusDv', # Mouse
                             'NCBITaxon:7955':'ZFS', # Zebrafish
                             'other':'HsapDv' # For all other organisms, use HsapDv
                            },
        'tissue':{'NCBITaxon:6239':'WBbt', # C. Elegans
                  'NCBITaxon:7227':'FBbt', # Drosophila
                  'NCBITaxon:7955':'ZFA', # Zebrafish
                  'other':'UBERON' # For all other organisms, use UBERON
                 }
    }
    ontology_parser = OntologyParser()
    ont_err_lst = []
    
    for col in col_ont_map:
        map_dict = {}
        for label in sample_df[col].unique():
            term_id = None
            if col == 'disease' and label == 'normal': # Normal is not in MONDO ontology
                term_id = 'PATO:0000461'
            elif label in ['unknown','na']: # Unknown and na won't be in ontologies, pass along
                map_dict[label] = label
                continue
            elif col in ['tissue','development_stage','self_reported_ethnicity']:
                if col == 'tissue':
                    # Find what tissue type is at label row
                    if sample_df.loc[sample_df[col] == label, 'tissue_type'].tolist()[0] != 'tissue':
                        map_dict[label] = label # Don't map cell type in tissue
                        continue
                # Find what organism term id is at label row
                org_term_id = sample_df.loc[sample_df[col] == label, 'organism_ontology_term_id'].tolist()[0]
                if org_term_id in col_ont_map[col]:
                    # Get ontology of specific organism and map label
                    species_ont = col_ont_map[col][org_term_id]
                    term_id = ontology_parser.get_term_id_by_label(label, species_ont)
                else:
                    if col_ont_map[col]['other'] == 'none':
                        map_dict[label] = label
                        continue
                    else:
                        term_id = ontology_parser.get_term_id_by_label(label, col_ont_map[col]['other'])
            else:
                term_id = ontology_parser.get_term_id_by_label(label, col_ont_map[col])
            if term_id == None:
                if org_term_id:
                    if org_term_id in col_ont_map[col]:
                        ont_err_lst.append(f"Error: Matching '{col_ont_map[col][org_term_id]}' term id not found for label '{label}' in column '{col}'")
                    else:
                        ont_err_lst.append(f"Error: Matching '{col_ont_map[col]['other']}' term id not found for label '{label}' in column '{col}'")
                else:
                    ont_err_lst.append(f"Error: Matching '{col_ont_map[col]}' term id not found for label '{label}' in column '{col}'")
                map_dict[label] = label
                continue
            map_dict[label] = term_id
        sample_df[col + '_ontology_term_id'] = sample_df[col].map(map_dict)
        del sample_df[col]
    
    ### Print out any errors from ontologizing
    if ont_err_lst:
        for e in ont_err_lst:
            print(e)
        sys.exit()

    ### Convert string to boolean for is_pilot_data and donor_living_at_sample_collection
    ### Check that donor_living_at_sample_collection is not filled out for non-human
    b_type = ['is_pilot_data','donor_living_at_sample_collection']
    for c in b_type:
        if c in sample_df.columns:
            if c == 'donor_living_at_sample_collection':
                for val in sample_df[c].unique():
                    if val != 'na' and sample_df.loc[sample_df[c] == val, 
                    'organism_ontology_term_id'].tolist()[0] != 'NCBITaxon:9606':
                        print(f"ERROR: donor_living_at_sample_collection for non-human data should be 'na' but '{val}' is present")
                        sys.exit()
            sample_df[c] == sample_df[c].replace({'FALSE':False, 'TRUE':True})
    
    ### Blank fields in worksheet result in NaN values in dataframe, replacing these with na?
    ### Could also replace with unknown for certain columns using fillna options?
    sample_df.fillna('na', inplace=True)
    sample_df.drop(columns=[c for c in sample_df.columns if c.startswith('!')], inplace=True)

    return sample_df


def get_gid(sheet, tab_name):
    '''
    Given sheet id and tab name, return gid
    '''
    sheet_url = f'https://docs.google.com/spreadsheets/d/{sheet}'
    req = Request(sheet_url, headers={'User-Agent' : "Magic Browser"})
    s = urlopen(req)
    soup = BeautifulSoup(s, 'html.parser')
    tab_ids = {}
    pattern = re.compile('var bootstrapData = (.*?)};')
    for s in soup.find_all('script'):
        if pattern.search(str(s)):
            d = pattern.search(str(s)).group()[20:-1]
            data = json.loads(d)
            for t in data['changes']['topsnapshot']:
                u = t[1].split('"')
                if len(u) > 5:
                    tab_ids[u[5]] = u[1]
    return tab_ids[tab_name]


args = getArgs()
s3client = boto3.client('s3')
FS = fsspec.filesystem("s3")

if __name__ == '__main__':
    
    ### Read in metadata and create ontologized sample_df
    tab_name = 'sample template'
    gid = get_gid(args.sheet, tab_name)
    url = f'https://docs.google.com/spreadsheets/d/{args.sheet}/export?format=csv&gid={gid}'
    response = requests.get(url)
    sample_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str).dropna(axis=1,how='all')
    check_standard_presence(sample_df)
    sample_df = map_ontologies(sample_df)

    ### Obtain subsamples associated with library that matches args.group
    my_dir = f'{args.project}/'
    r = s3client.list_objects(Bucket=args.bucket, Prefix=my_dir, Delimiter="/")
    orders = [o['Prefix'].replace(my_dir,'') for o in r['CommonPrefixes']]
    samples = []
    for o in [i for i in orders if i.startswith(('AN','NV'))]:
        r = s3client.list_objects(Bucket=args.bucket, Prefix=f'{my_dir}{o}', Delimiter='/')
        libs = [l['Prefix'].replace(f'{my_dir}{o}','') for l in r['CommonPrefixes']]
        for l in libs:
            if l == args.group+"/":
                r = s3client.list_objects(Bucket=args.bucket, Prefix=f'{my_dir}{o}{l}processed/cellranger/', Delimiter='/')
                dates = [d['Prefix'].replace(f'{my_dir}{o}{l}processed/cellranger/','') for d in r['CommonPrefixes']]
                for d in dates:
                    r = s3client.list_objects(Bucket=args.bucket, Prefix=f'{my_dir}{o}{l}processed/cellranger/{d}outs/per_sample_outs/', Delimiter='/')
                    subs = [s['Prefix'].replace(f'{my_dir}{o}{l}processed/cellranger/{d}outs/per_sample_outs/','') for s in r['CommonPrefixes']]
                    for s in subs:
                        samples.append({
                            'order': o,
                            'library': l.rstrip('/'),
                            'date': d,
                            'sample': s.rstrip('/')
                        })

    
    ### for all samples associated with GroupID, download files and run final curated matrix generation
    for s in samples:
        order = s['order']
        lib = s['library']
        run_date = s['date']
        samp = s['sample']
        print(f"Processing {s} from {lib}")

        ### Download and read in h5 file
        for d in ['temp_cellranger','curated_matrices']:
            if os.path.exists(d) == False:
                os.mkdir(d)
        s_dir = f'{args.project}/{order}{lib}/processed/cellranger/{run_date}outs/per_sample_outs/{samp}'
        lib_samp = f'{lib}_{samp}'
        download_files(s_dir, args.bucket, lib_samp)
        h5_file = f'temp_cellranger/{lib_samp}_sample_filtered_feature_bc_matrix.h5'
        adata = sc.read_10x_h5(h5_file, gex_only=True)

        ### Track additional tab gids from Lattice wrangling sheet
        guide_gid = get_gid(args.sheet, 'guide template')

        adata = custom_var_to_obs(adata)
        
        adata.obs['library_name'] = lib
        adata.obs['lab'] = args.project.split('-')[0]
        adata.obs['project'] = args.project
        adata.obs['sample_name'] = samp
        
        adata.obs = adata.obs.merge(sample_df, on='sample_name', how='left').set_index(adata.obs.index)
        adata.uns['organism_ontology_term_id'] = adata.obs['organism_ontology_term_id'].unique()[0]
        adata.uns['title'] = lib_samp
        adata.obs.drop(columns=['organism_ontology_term_id'],inplace=True)
        adata.obs['is_primary_data'] = True
        
        crispr_df = gather_crispr(lib_samp)
        adata.obs = adata.obs.merge(
            crispr_df, left_index=True, right_index=True, how='left'
        ).set_index(adata.obs.index)
        
        metrics_df = gather_metrics(samp, lib_samp)
        for c in metrics_df.columns:
            adata.obs[c] = metrics_df[c].values[0]
        
        adata = map_filter_gene_ids(adata)
        cxg_add_labels(adata)
        adata.var.drop(columns=['feature_types','genome'], inplace=True)
        adata.var['feature_is_filtered'] = False
        
        ### Startswith('mt-') works for human & mouse, more attn needed for other organisms
        adata.var['mt'] = adata.var['feature_name'].str.lower().str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        
        ### Add guide schema metadata to adata.uns and adata.obs
        adata = add_guide_metadata(adata, args.sheet, guide_gid, args.csvofguidescan)
        adata = determine_perturbation_strategy(adata)
        adata.obs['genetic_perturbation_id'] = adata.obs['genetic_perturbation_id'].astype('category')
        adata.obs['genetic_perturbation_id'] = adata.obs['genetic_perturbation_id'].cat.add_categories(['na'])
        adata.obs['genetic_perturbation_id'].fillna('na', inplace=True)

        ### Write to directories: curated and temp
        order = order.rstrip('/')
        adata.write(filename=f'curated_matrices/{lib}__{samp}__{order}__curated.h5ad', compression='gzip')
        for file in ['sample_filtered_feature_bc_matrix.h5', 'crispr_analysis.tar.gz', 'metrics_summary.csv', 'protospacer_calls_per_cell.csv']:
            if os.path.isfile(f"temp_cellranger/{lib_samp}_{file}"):
                os.remove(f"temp_cellranger/{lib_samp}_{file}")



