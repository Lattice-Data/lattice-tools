import boto3
import json
import os
import pandas as pd
import re
from bs4 import BeautifulSoup


s3client = boto3.client('s3')


chemistries = {
    "Single Cell 5' R2-only v3": "5p",
    "Single Cell 5' R2-only": "5p",
    "Single Cell 3' v4 (polyA)": "3p",
    "Single Cell 3' v3": "3p",
    "Flex Gene Expression": "flex"
}

valid_assays = [
    'CRI','GEX','ATAC','viral_ORF',
    'GEX_hash_oligo','hash_oligo'
]

#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex
#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-flex-outputs-frp
#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview
cellranger_expected = {
    'cellranger-9.0.1':{
        'nonflex': {
            'outs': [
                'config.csv',
                'multi/count/feature_reference.csv',
                'multi/count/raw_cloupe.cloupe',
                'multi/count/raw_feature_bc_matrix.h5',
                'multi/count/raw_feature_bc_matrix.tar.gz',
                'multi/count/raw_molecule_info.h5',
                'multi/count/unassigned_alignments.bam',
                'multi/count/unassigned_alignments.bam.bai'
            ],
            'per_sample': [
                'count/analysis.tar.gz',
                'count/feature_reference.csv',
                'count/sample_cloupe.cloupe',
                'count/sample_alignments.bam',
                'count/sample_alignments.bam.bai',
                'count/sample_filtered_barcodes.csv',
                'count/sample_filtered_feature_bc_matrix.h5',
                'count/sample_filtered_feature_bc_matrix.tar.gz',
                'count/sample_molecule_info.h5',
                'metrics_summary.csv',
                'web_summary.html'
            ]
        },
        'flex': {
            'outs': [
                'config.csv',
                'multi/count/raw_cloupe.cloupe',
                'multi/count/raw_feature_bc_matrix.h5',
                'multi/count/raw_feature_bc_matrix.tar.gz',
                'multi/count/raw_molecule_info.h5',
                'multi/count/raw_probe_bc_matrix.h5'
            ],
            'per_sample': [
                'count/analysis.tar.gz',
                'count/probe_set.csv',
                'count/sample_cloupe.cloupe',
                'count/sample_filtered_barcodes.csv',
                'count/sample_filtered_feature_bc_matrix.h5',
                'count/sample_filtered_feature_bc_matrix.tar.gz',
                'count/sample_molecule_info.h5',
                'count/sample_raw_feature_bc_matrix.h5',
                'count/sample_raw_feature_bc_matrix.tar.gz',
                'count/sample_raw_probe_bc_matrix.h5',
                'metrics_summary.csv',
                'web_summary.html'
            ]
        }
    },
    'cellranger-10.0.0': {
        'nonflex': {
            'outs': [
                'config.csv',
                'filtered_feature_bc_matrix/barcodes.tsv.gz',
                'filtered_feature_bc_matrix/features.tsv.gz',
                'filtered_feature_bc_matrix/matrix.mtx.gz',
                'filtered_feature_bc_matrix.h5',
                'multiplexing_analysis/cells_per_tag.json',
                'qc_library_metrics.csv',
                'qc_report.html',
                'qc_sample_metrics.csv',
                'raw_cloupe.cloupe',
                'raw_feature_bc_matrix/barcodes.tsv.gz',
                'raw_feature_bc_matrix/features.tsv.gz',
                'raw_feature_bc_matrix/matrix.mtx.gz',
                'raw_feature_bc_matrix.h5',
                'raw_molecule_info.h5'
            ],
            'per_sample_outs': [
                'sample_filtered_feature_bc_matrix/barcodes.tsv.gz',
                'sample_filtered_feature_bc_matrix/features.tsv.gz',
                'sample_filtered_feature_bc_matrix/matrix.mtx.gz',
                'sample_raw_feature_bc_matrix/barcodes.tsv.gz',
                'sample_raw_feature_bc_matrix/features.tsv.gz',
                'sample_raw_feature_bc_matrix/matrix.mtx.gz',
                'metrics_summary.csv',
                'sample_cloupe.cloupe',
                'sample_filtered_barcodes.csv',
                'sample_filtered_feature_bc_matrix.h5',
                'sample_molecule_info.h5',
                'sample_raw_feature_bc_matrix.h5',
                'web_summary.html'
            ]
        },
        'flex': {
            'outs':[
                'config.csv',
                'feature_reference.csv',
                'filtered_feature_bc_matrix/barcodes.tsv.gz',
                'filtered_feature_bc_matrix/features.tsv.gz',
                'filtered_feature_bc_matrix/matrix.mtx.gz',
                'filtered_feature_bc_matrix.h5',
                'multiplexing_analysis/cells_per_tag.json',
                'multiplexing_analysis/frp_gem_barcode_overlap.csv',
                'probe_set.csv',
                'qc_library_metrics.csv',
                'qc_report.html',
                'qc_sample_metrics.csv',
                'raw_cloupe.cloupe',
                'raw_feature_bc_matrix/barcodes.tsv.gz',
                'raw_feature_bc_matrix/features.tsv.gz',
                'raw_feature_bc_matrix/matrix.mtx.gz',
                'raw_feature_bc_matrix.h5',
                'raw_molecule_info.h5',
                'raw_probe_bc_matrix.h5'
            ],
            'per_sample_outs': [
                'crispr_analysis/protospacer_calls_per_cell.csv',
                'metrics_summary.csv',
                'sample_cloupe.cloupe',
                'sample_filtered_barcodes.csv',
                'sample_filtered_feature_bc_matrix.h5',
                'sample_filtered_feature_bc_matrix/barcodes.tsv.gz',
                'sample_filtered_feature_bc_matrix/features.tsv.gz',
                'sample_filtered_feature_bc_matrix/matrix.mtx.gz',
                'sample_molecule_info.h5',
                'sample_raw_feature_bc_matrix.h5',
                'sample_raw_feature_bc_matrix/barcodes.tsv.gz',
                'sample_raw_feature_bc_matrix/features.tsv.gz',
                'sample_raw_feature_bc_matrix/matrix.mtx.gz',
                'sample_raw_probe_bc_matrix.h5',
                'web_summary.html'
            ]
        }
    },
    'count': {
        'outs': [
            'analysis/clustering/gene_expression_graphclust/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_3_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_4_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_5_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_6_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_7_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_8_clusters/clusters.csv',
            'analysis/clustering/gene_expression_kmeans_9_clusters/clusters.csv',
            'analysis/diffexp/gene_expression_graphclust/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_10_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_2_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_3_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_4_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_5_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_6_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_7_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_8_clusters/differential_expression.csv',
            'analysis/diffexp/gene_expression_kmeans_9_clusters/differential_expression.csv',
            'analysis/pca/gene_expression_10_components/components.csv',
            'analysis/pca/gene_expression_10_components/dispersion.csv',
            'analysis/pca/gene_expression_10_components/features_selected.csv',
            'analysis/pca/gene_expression_10_components/projection.csv',
            'analysis/pca/gene_expression_10_components/variance.csv',
            'analysis/tsne/gene_expression_2_components/projection.csv',
            'analysis/umap/gene_expression_2_components/projection.csv',
            'cloupe.cloupe',
            'filtered_feature_bc_matrix/barcodes.tsv.gz',
            'filtered_feature_bc_matrix/features.tsv.gz',
            'filtered_feature_bc_matrix/matrix.mtx.gz',
            'filtered_feature_bc_matrix.h5',
            'metrics_summary.csv',
            'molecule_info.h5',
            'possorted_genome_bam.bam',
            'possorted_genome_bam.bam.bai',
            'raw_feature_bc_matrix/barcodes.tsv.gz',
            'raw_feature_bc_matrix/features.tsv.gz',
            'raw_feature_bc_matrix/matrix.mtx.gz',
            'raw_feature_bc_matrix.h5',
            'web_summary.html'
        ]
    }
}

raw_expected = {
    'sci_jumbo': [
        '.cram',
        '.cram-metadata.json',
        '.csv',
        '.json',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv',
        '_FlowQ.metric',
        '_SNVQ.metric'
    ],
    'sci_plex': [
        '.cram',
        '.cram-metadata.json',
        '.csv',
        '.json',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv'
    ],
    '10x': [
        '.csv',
        '.json',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv',
        '_unmatched.cram',
        '_unmatched.cram-metadata.json',
        '_unmatched.csv',
        '_unmatched.json',
        '_S1_L001_R1_001.csv',
        '_S1_L001_R1_001.fastq.gz',
        '_S1_L001_R1_001.fastq.gz-metadata.json',
        '_S1_L001_R1_001.json',
        '_S1_L001_R1_001_sample.fastq.gz',
        '_S1_L001_R1_001_sample.fastq.gz-metadata.json',
        '_S1_L001_R2_001.csv',
        '_S1_L001_R2_001.fastq.gz',
        '_S1_L001_R2_001.fastq.gz-metadata.json',
        '_S1_L001_R2_001.json',
        '_S1_L001_R2_001_sample.fastq.gz',
        '_S1_L001_R2_001_sample.fastq.gz-metadata.json'
    ],
    '10x_viral_ORF': [
        '.csv',
        '.json',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv',
        '.cram',
        '.cram-metadata.json',
        '_FlowQ.metric',
        '_SNVQ.metric'
    ]
}

raw_optional = {
    '10x': [
        '.scRNA.applicationQC.h5',
        '.scRNA.applicationQC.html',
        '_Log.final.out',
        '_Log.out',
        '_Log.progress.out',
        '_ReadsPerGene.out.tab',
        '_SJ.out.tab'
    ]
}


def parse_met_summ(f):
    df = pd.read_csv(f)
    if len(df) == 1:
        report = {
            'GEX_reads': int(df['Number of Reads'].iloc[0].replace(',',''))
        }

        return report

    lib_reads = df[
        (df['Metric Name'].isin(['Number of reads','Number of short reads skipped'])) &
        (df['Grouped By'] == 'Fastq ID')
    ]
    lib_reads.loc[:, 'Metric Value'] = lib_reads['Metric Value'].str.replace(',', '').astype(int)

    gex_reads = lib_reads[lib_reads['Library Type'] == 'Gene Expression']
    report = {'GEX_reads': gex_reads['Metric Value'].sum()}

    if 'CRISPR Guide Capture' in lib_reads['Library Type'].unique():
        cri_reads = lib_reads[lib_reads['Library Type'] == 'CRISPR Guide Capture']
        report['CRI_reads'] = cri_reads['Metric Value'].sum()

    return report


def parse_web_summ(f):
    report = {
        'extra': [],
        'alerts': []
    }
    with open(f) as html_doc:
        soup = BeautifulSoup(html_doc, 'html.parser')
    for x in soup.find_all('script'):
        match = re.search("const data = ", x.string)
        if match:
            end = match.end()
            data = json.loads(x.string[end:])
    
    if 'summary' in data:
        report['sub'] = data['summary']['sample']['subcommand']
        gex_tab = {row[0].lower():row[1] for row in data['summary']['summary_tab']['pipeline_info_table']['rows']}
    else:
        gex_tab = {row[0].lower():row[1] for row in data['library']['data']['gex_tab']['content']['parameters_table']['rows']}
        if 'sample' in data:
            report['sub'] = data['sample']['subcommand']


    chem = gex_tab['chemistry']
    report['chem'] = chemistries.get(chem, chem)

    report['extra'] = []
    if 'library' in data:
        if data['library']['data']['crispr_tab']:
            report['extra'].append('CRISPR')
            crispr_tab = {f'crispr {row[0]}':row[1] for row in data['library']['data']['crispr_tab']['content']['parameters_table']['rows']}
        if data['library']['data']['antibody_tab']:
            report['extra'].append('Antibody')
    if 'experimental_design' in data:
        for line in data['experimental_design']['csv'].split('\n'):
            if line.startswith('['):
                cat = line
                if cat == '[samples]':
                    report['multiplex'] = True
            if ',' in line:
                path = line.strip().split(',')
                if path[0] == 'skip-cell-annotation' and path[1] == 'false':
                    report['extra'].append('CellAnnotate')
                elif path[0] == 'min-crispr-umi':
                    report['min-crispr-umi'] = path[1]
                elif path[0] == 'create-bam':
                    report['create-bam'] = path[1]
                elif path[0] == 'reference' and cat == '[gene-expression]':
                    report['ref'] = path[1]

    #location of some additional info to QA
    report['ref'] = gex_tab['transcriptome']
    if chem != 'Flex Gene Expression':
        report['incl_int'] = gex_tab['include introns'].lower()

    if 'pipeline version' in gex_tab:
        report['software'] = gex_tab['pipeline version']
    elif 'pipeline_version' in data:
        report['software'] = 'cellranger-'+data['pipeline_version']
    elif 'library' in data:
        report['software'] = data['library']['data']['header_info']['Pipeline Version']
        report['gex_alerts'] = data['library']['data']['gex_tab']['alerts']
        if data['library']['data']['crispr_tab']:
            report['crispr_alerts'] = data['library']['data']['crispr_tab']['alerts']

    if 'summary' in data:
        report['gex_alerts'] = data['summary']['alarms']['alarms']

    return report


def grab_trimmer_stats(trimmer_failure_stats, rf, bucket):
    exp = '/'.join(rf.split('/')[1:3])
    if exp not in trimmer_failure_stats:
        trimmer_failure_stats[exp] = {'rsq': [], 'trimmer_fail': []}
    s3client.download_file(bucket, rf, 'trimmer-failure_codes.csv')
    trimmer_fail = 0
    stats_df = pd.read_csv('trimmer-failure_codes.csv')
    stats_df.columns = stats_df.columns.str.replace(' ', '_')
    for row in stats_df.itertuples():
        total_reads = row.total_read_count
        if row.reason == 'rsq file':
            trimmer_failure_stats[exp]['rsq'].append(100*row.failed_read_count/total_reads)
        else:
            trimmer_fail += row.failed_read_count
    trimmer_fail_pct = trimmer_fail/total_reads
    trimmer_failure_stats[exp]['trimmer_fail'].append(100*trimmer_fail_pct)
    os.remove('trimmer-failure_codes.csv')


def parse_raw_filename(f):
    path = f.split('/')[-1].split('-')
    if len(path) < 3:
        return None

    run = path[0]
    group_assay = path[1]
    if group_assay.endswith('GEX_hash_oligo'):
        assay = 'GEX_hash_oligo'
    else:
        match = False
        for v_a in valid_assays:
            if group_assay.endswith(v_a):
                assay = v_a
                match = True
        if not match:
            assay = group_assay.split('_')[-1]
    group = group_assay.replace(f'_{assay}','')
    ug = path[2]
    barcode = path[3].split('_')[0].split('.')[0]

    return run,group,assay,ug,barcode


def load_files_from_manifest(
    manifest_path: str,
    delimiter: str,
    s3_column: int,
    has_header: bool = False
):
    """
    Load S3 file paths from a CSV/TSV manifest file.

    Args:
        manifest_path: Path to the manifest file
        delimiter: Field delimiter (',' for CSV, '\t' for TSV)
        s3_column: Column index (0-based) containing S3 URIs
        has_header: Whether the file has a header row to skip

    Returns:
        Tuple of:
        - all_raw_files: List of raw file S3 keys (paths containing '/raw/')
        - all_proc_files: Dict of {group: [processed file keys]} (paths containing '/processed/')
    """
    df = pd.read_csv(manifest_path, sep=delimiter, header=0 if has_header else None)
    s3_uris = df.iloc[:, s3_column].tolist()

    all_raw_files = []
    all_proc_files = {}

    for uri in s3_uris:
        # Strip 's3://bucket-name/' prefix to get just the key
        if uri.startswith('s3://'):
            # s3://bucket-name/path/to/file -> path/to/file
            parts = uri[5:].split('/', 1)  # Remove 's3://', split on first '/'
            if len(parts) > 1:
                key = parts[1]  # The key (path after bucket)
            else:
                continue  # Invalid URI, skip
        else:
            key = uri  # Assume it's already a key

        # Separate into raw vs processed
        if '/raw/' in key:
            all_raw_files.append(key)
        elif '/processed/' in key:
            # Extract group name from path: .../GROUP/processed/...
            path_before_processed = key.split('/processed/')[0]
            group = path_before_processed.split('/')[-1]
            if group not in all_proc_files:
                all_proc_files[group] = []
            all_proc_files[group].append(key)

    return all_raw_files, all_proc_files
