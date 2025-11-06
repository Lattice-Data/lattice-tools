import json
import pandas as pd
import re
from bs4 import BeautifulSoup


chemistries = {
    "Single Cell 5' R2-only v3": "5p",
    "Single Cell 5' R2-only": "5p",
    "Single Cell 3' v4 (polyA)": "3p",
    "Single Cell 3' v3": "3p",
    "Flex Gene Expression": "flex"
}

#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex
#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-flex-outputs-frp
#https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview
cellranger_expected = {
    'nonflex': {
        'outs': [
            'config.csv',
            'multi/count/feature_reference.csv',
            'multi/count/raw_cloupe.cloupe',
            'multi/count/raw_feature_bc_matrix.h5',
            'multi/count/raw_feature_bc_matrix.tar.gz',
            'multi/count/raw_molecule_info.h5',
            'multi/count/unassigned_alignments.bam',
            'multi/count/unassigned_alignments.bam.bai',
            'multi/multiplexing_analysis.tar.gz'
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
            'multi/count/raw_probe_bc_matrix.h5',
            'multi/multiplexing_analysis.tar.gz'
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
    },
    'count': {
        'outs': [
            'analysis.tar.gz',
            'cloupe.cloupe',
            'filtered_feature_bc_matrix.tar.gz',
            'filtered_feature_bc_matrix.h5',
            'metrics_summary.csv',
            'molecule_info.h5',
            'possorted_genome_bam.bam',
            'possorted_genome_bam.bam.bai',
            'raw_feature_bc_matrix.tar.gz',
            'raw_feature_bc_matrix.h5',
            'web_summary.html'
        ]
    }
}

raw_expected = {
    'standard': [
        '.csv',
        '.json',
        '.scRNA.applicationQC.h5',
        '.scRNA.applicationQC.html',
        '_Log.final.out',
        '_Log.out',
        '_Log.progress.out',
        '_ReadsPerGene.out.tab',
        '_SJ.out.tab',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv',
        '_unmatched.cram',
        '_unmatched.csv',
        '_unmatched.json',
        '_S1_L001_R1_001.csv',
        '_S1_L001_R1_001.fastq.gz',
        '_S1_L001_R1_001.json',
        '_S1_L001_R1_001_sample.fastq.gz',
        '_S1_L001_R2_001.csv',
        '_S1_L001_R2_001.fastq.gz',
        '_S1_L001_R2_001.json',
        '_S1_L001_R2_001_sample.fastq.gz'
    ],
    'viralORF': [
        '.csv',
        '.json',
        '_trimmer-failure_codes.csv',
        '_trimmer-stats.csv',
        '.cram',
        '_FlowQ.metric',
        '_SNVQ.metric'
    ]
}


def parse_met_summ(f):
    df = pd.read_csv(f)
    lib_reads = df[
        (df['Metric Name'].isin(['Number of reads','Number of short reads skipped'])) &
        (df['Grouped By'] == 'Fastq ID')
    ]
    lib_reads.loc[:, 'Metric Value'] = lib_reads['Metric Value'].str.replace(',', '').astype(int)

    gex_reads = lib_reads[lib_reads['Library Type'] == 'Gene Expression']
    print('GEX', gex_reads['Metric Value'].sum())

    if 'CRISPR Guide Capture' in lib_reads['Library Type'].unique():
        cri_reads = lib_reads[lib_reads['Library Type'] == 'CRISPR Guide Capture']
        print('CRI', cri_reads['Metric Value'].sum())


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

    report['sub'] = data['sample']['subcommand']
    gex_tab = {row[0]:row[1] for row in data['library']['data']['gex_tab']['content']['parameters_table']['rows']}
    chem = gex_tab['Chemistry']
    report['chem'] = chemistries.get(chem, chem)
    
    report['extra'] = []
    if data['library']['data']['crispr_tab']:
        report['extra'].append('CRISPR')
        crispr_tab = {f'crispr {row[0]}':row[1] for row in data['library']['data']['crispr_tab']['content']['parameters_table']['rows']}
    if data['library']['data']['antibody_tab']:
        report['extra'].append('Antibody')
    for line in data['experimental_design']['csv'].split('\n'):
        if ',' in line:
            path = line.strip().split(',')
            if path[0] == 'skip-cell-annotation' and path[1] == 'false':
                report['extra'].append('CellAnnotate')
            elif path[0] == 'min-crispr-umi':
                report['min-crispr-umi'] = path[1]
            elif path[0] == 'create-bam':
                report['create-bam'] = path[1]

    #location of some additional info to QA
    report['ref'] = gex_tab['Transcriptome']
    if chem != 'Flex Gene Expression':
        report['incl_int'] = gex_tab['Include Introns']

    report['software'] = data['library']['data']['header_info']['Pipeline Version']

    for m in data['library']['metrics']:
        if m['key'] == 'number_of_reads' and m['grouping_key'] not in ['GEX_1','CGC_1']:
            report[m['grouping_key']] = m['value']

    report['gex_alerts'] = data['library']['data']['gex_tab']['alerts']

    if data['library']['data']['crispr_tab']:
        report['crispr_alerts'] = data['library']['data']['crispr_tab']['alerts']

    return report