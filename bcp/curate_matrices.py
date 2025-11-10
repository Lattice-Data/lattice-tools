import boto3
import numpy as np
import os
import pandas as pd
import scanpy as sc
import tarfile
from cellxgene_schema.write_labels import AnnDataLabelAppender
from urllib.parse import quote
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from cellxgene_ontology_guide.supported_versions import CXGSchema, load_supported_versions
import requests
from io import BytesIO
from urllib.request import Request,urlopen
from bs4 import BeautifulSoup
import re
import json
import argparse
import sys
sys.path.append(os.path.dirname(os.path.abspath('../cellxgene_resources')))
from cellxgene_resources.cellxgene_mods import map_filter_gene_ids



EPILOG = f"""
Script will take a single GroupID from a project on CZI S3 as input
and create curated matrices for all subsamples of this GroupID. This will need to be
run on JupyterHub to allow access to CZI S3 buckets and reference indices for Guidescan.
The current assumption is that this is run on 10x Flex data, and will need to
update anticipated directory structure if it is another 10x assay.

Example:
    python curate_matrices.py --bucket czi-psomagen --sheet your_sheet_id --project marson-mapping-grns-perturb-seq --group CD4i_R1L01

For more details:
    python %(prog)s --help

"""


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    args = parser.parse_args()
    if len(sys.argv) < 4:
    	parser.print_help()
    	sys.exit()

    return args


def download_files(s_dir, bucket, samp):
    mx_h5 = f'{s_dir}/count/sample_filtered_feature_bc_matrix.h5'
    metrics_csv = f'{s_dir}/metrics_summary.csv'
    cri_tar = f'{s_dir}/count/crispr_analysis.tar.gz'
    
    for file_path in [mx_h5, metrics_csv, cri_tar]:
        f = file_path.split('/')[-1]
        full_path = os.path.join("temp_cellranger/",f"{samp}_{f}")
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



def gather_crispr(samp):
    '''
    obtain guide assignement from 10x cellranger output files

    return df: dataframe containing crispr guide calling
    '''
    df = pd.DataFrame()
    with tarfile.open(f'temp_cellranger/{samp}_crispr_analysis.tar.gz', 'r:gz') as tar:
        f = tar.extractfile('protospacer_calls_per_cell.csv')
        df = pd.read_csv(f).rename(columns={'num_umis':'num_umis_guide_id'})
    df['genetic_perturbation_id'] = df['feature_call'].apply(lambda x: x.replace('|',' || '))
    df['num_umis_guide_id'] = df['num_umis_guide_id'].apply(lambda x: x.replace('|',' || '))
    df = df[['cell_barcode','genetic_perturbation_id','num_umis_guide_id']].set_index('cell_barcode')

    return df


def gather_metrics(samp):
    df = pd.read_csv(f'temp_cellranger/{samp}_metrics_summary.csv')
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


def add_guide_metadata(adata, sheet, guide_gid):
    '''
    Add guide metadata into adata.uns from Lattice wrangling sheet
    
    :param obj adata: the anndata object that is being transformed into the curated matrix
    :param obj guide_df: the dataframe containing guide metadata from wrangling sheet
    
    :returns obj adata: modified adata to contain guide metadata
    '''
    url = f'https://docs.google.com/spreadsheets/d/{sheet}/export?format=csv&gid={guide_gid}'
    response = requests.get(url)
    guide_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str)
    genetic_perturbations = {}
    
    for row in guide_df.itertuples():
        genetic_perturbations[row.guide_id] = {}
        genetic_perturbations[row.guide_id]['role'] = 'targeting' if row.guide_role == 'Targeting a Gene' else 'control'
        genetic_perturbations[row.guide_id]['protospacer_sequence'] = row.guide_protospacer
        genetic_perturbations[row.guide_id]['protospacer_adjacent_motif'] = row.guide_PAM
        if not pd.isna([row.start,row.end,row.strand]).all():
            chr_loc = str(row.chromosome).replace("chr","") + ":" + str(row.start) + "-" + str(row.end) + "(" + str(row.strand) + ")"
            genetic_perturbations[row.guide_id]['target_genomic_regions'] = [chr_loc]
        if not pd.isna(row.overlapping_gene_ids):
            genetic_perturbations[row.guide_id]['target_features'] = {}
            for i in range(len(row.overlapping_gene_ids.split(";"))):
                genetic_perturbations[row.guide_id]['target_features'][row.overlapping_gene_ids.split(";")[i]] = row.overlapping_gene_names.split(";")[i]
                                                                             
            
    adata.uns['genetic_perturbations'] = genetic_perturbations
    
    return adata


def determine_perturbation_strategy(adata):
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
    adata.obs['genetic_perturbation_strategy_calculated'] = adata.obs['genetic_perturbation_id']
    adata.obs['genetic_perturbation_strategy_calculated'] = adata.obs['genetic_perturbation_strategy_calculated'].apply(
        lambda x: x.split(' || ') if pd.notna(x) else 'no perturbations'
    )
    adata.obs['genetic_perturbation_strategy_calculated'] = adata.obs['genetic_perturbation_strategy_calculated'].apply(
        lambda x: [adata.uns['genetic_perturbations'][i]['role'] for i in x] if isinstance(x, list)
            else x        
    )
    adata.obs['genetic_perturbation_strategy_calculated'] = adata.obs['genetic_perturbation_strategy_calculated'].apply(
         lambda x: 'control' if isinstance(x, list) and 'targeting' not in set(x)
            else x
    )
    adata.obs.loc[adata.obs['genetic_perturbation_strategy_calculated']=='control', 'genetic_perturbation_strategy'] = 'control'
    adata.obs.loc[adata.obs['genetic_perturbation_strategy_calculated']=='no perturbations', 'genetic_perturbation_strategy'] = 'no perturbations'
    adata.obs.drop(columns=['genetic_perturbation_strategy_calculated'], inplace=True)
    
    return adata


def map_ontologies(sample_df):
    '''
    Takes the sample metadata dataframe and standardizes ontologies

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
                print(f"Matching '{col_ont_map[col]}' term id not found for label '{label}' in column '{col}'")
                map_dict[label] = label
                continue
            map_dict[label] = term_id
        sample_df[col + '_ontology_term_id'] = sample_df[col].map(map_dict)
        del sample_df[col]
    
    ### Convert string to boolean for is_pilot_data and donor_living_at_sample_collection
    ### Will need to look further into this for next iteration
    b_type = ['is_pilot_data','donor_living_at_sample_collection']
    for c in b_type:
        sample_df[c] = sample_df[c].replace({'FALSE':False, 'TRUE':True})
    
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

if __name__ == '__main__':
    
    ### Read in metadata and create ontologized sample_df
    tab_name = 'sample template'
    gid = get_gid(args.sheet, tab_name)
    url = f'https://docs.google.com/spreadsheets/d/{args.sheet}/export?format=csv&gid={gid}'
    response = requests.get(url)
    sample_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str).dropna(axis=1,how='all')
    sample_df = map_ontologies(sample_df)

    ### Obtain subsamples associated with library that matches args.group
    my_dir = f'{args.project}/'
    r = s3client.list_objects(Bucket=args.bucket, Prefix=my_dir, Delimiter="/")
    orders = [o['Prefix'].replace(my_dir,'') for o in r['CommonPrefixes']]
    samples = []
    for o in orders:
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
        download_files(s_dir, args.bucket, samp)
        h5_file = f'temp_cellranger/{samp}_sample_filtered_feature_bc_matrix.h5'
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
        adata.uns['title'] = samp
        adata.obs.drop(columns=['organism_ontology_term_id'],inplace=True)
        adata.obs['is_primary_data'] = True
        
        crispr_df = gather_crispr(samp)
        adata.obs = adata.obs.merge(
            crispr_df, left_index=True, right_index=True, how='left'
        ).set_index(adata.obs.index)
        
        metrics_df = gather_metrics(samp)
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
        adata = add_guide_metadata(adata, args.sheet, guide_gid)
        adata = determine_perturbation_strategy(adata)
        adata.obs['genetic_perturbation_id'] = adata.obs['genetic_perturbation_id'].astype('category')
        adata.obs['genetic_perturbation_id'] = adata.obs['genetic_perturbation_id'].cat.add_categories(['na'])
        adata.obs['genetic_perturbation_id'].fillna('na', inplace=True)

        ### Write to directories: curated and temp
        order = order.rstrip('/')
        adata.write(filename=f'curated_matrices/{lib}__{samp}__{order}__curated.h5ad', compression='gzip')
        for file in [f'sample_filtered_feature_bc_matrix.h5', f'crispr_analysis.tar.gz', f'metrics_summary.csv']:
            os.remove(f"temp_cellranger/{samp}_{file}")



