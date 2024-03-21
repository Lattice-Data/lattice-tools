import argparse
import json
import os
import pandas as pd
import subprocess
import sys


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--version', '-v',
                        help='The GENCODE version to map previous versions to.')
    args = parser.parse_args()
    return args


#https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver
#hg19ToHg38.over.chain.gz from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
def lift_over(df):
    temp_file = 'preLift.tsv'

    coord_cols = ['seqname','start','stop']
    og_columns = df.columns

    slim_df = df.reset_index()[coord_cols + ['index']]
    slim_df['seqname'] = slim_df['seqname'].apply(lambda x: 'chr' + str(x) if str(x).isnumeric() else x)
    slim_df.to_csv(temp_file, header=False, index=False, sep='\t')

    chain_file = 'hg19ToHg38.over.chain.gz'
    converted = 'conversions.bed'
    unmapped = 'unMapped'

    lift_process = subprocess.run(['./liftOver', temp_file, chain_file, converted, unmapped], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in lift_process.stdout.decode('utf-8').split('\n'):
        print(line)
    for line in lift_process.stderr.decode('utf-8').split('\n'):
        print(line)

    lifted_df = pd.read_csv(converted, names=slim_df.columns, sep='\t')
    
    df.drop(columns=coord_cols, inplace=True)
    df.reset_index(inplace=True)
    df = df.merge(lifted_df, on='index', how='inner')
    df = df[og_columns]

    os.remove(temp_file)
    os.remove(converted)
    os.remove(unmapped)

    return df


def create_df(file, hg19=False):
    column_names = ['seqname','source','feature','start','stop','score','strand','frame','expand_me']
    gtf_df = pd.read_table(file, names=column_names, comment='#')

    #parse out transcripts
    transcript_df = gtf_df[gtf_df['feature'] == 'transcript']
    transcript_df = transcript_df['expand_me'].str.split(';',expand=True)
    for c in transcript_df.columns:
        if str(transcript_df[c].iloc[0]).startswith(' transcript_id'):
            transcript_col = c
    transcript_df = transcript_df[[0,transcript_col]]
    transcript_df.rename(columns={0: 'gene_id', transcript_col: 'transcript_id'}, inplace=True)
    transcript_df['gene_id'] = transcript_df['gene_id'].apply(lambda x: x.split('"')[1].split('.')[0])
    transcript_df['transcript_id'] = transcript_df['transcript_id'].apply(lambda x: x.split('"')[1].split('.')[0])
    transcript_df = transcript_df.groupby('gene_id').agg(lambda x: ','.join(sorted(set(x)))).reset_index()

    #parse out genes
    gene_df = gtf_df[gtf_df['feature'] == 'gene']
    gene_df[['gene_id','gene_name']] = gene_df['expand_me'].str.split(';',expand=True)[[0,2]]
    gene_df['gene_name'] = gene_df['gene_name'].apply(lambda x: x.split('"')[1])
    gene_df['gene_id'] = gene_df['gene_id'].apply(lambda x: x.split('"')[1].split('.')[0])

    #merge transcript lists to the gene table
    gene_df = gene_df.merge(transcript_df, on='gene_id', how='left')

    if hg19:
        gene_df = lift_over(gene_df)
    gene_df['coordinates'] = gene_df.apply(
        lambda x: f"{x['seqname']}:{x['start']}-{x['stop']}{x['strand']}",
        axis=1)
    gene_df = gene_df[['gene_id','gene_name','transcript_id','coordinates']]

    return gene_df


def find_map(x):
    if x['gene_id_new_fromCoordinates'] == x['gene_id_new_fromTranscripts']:
        return x['gene_id_new_fromTranscripts']
    elif x['gene_id_new_fromTranscripts'] == '':
        return x['gene_id_new_fromCoordinates']
    elif x['gene_id_new_fromCoordinates'] == '':
        return x['gene_id_new_fromTranscripts']
    else:
        return 'mismatch mappings'


def compare(schema_df, input_df):
    not_in_schema = input_df[input_df['gene_id'].isin(schema_df['gene_id']) == False].drop(columns='gene_name')

    not_in_schema = not_in_schema.merge(
        schema_df[['transcript_id','gene_id']], how='left',
        on='transcript_id', suffixes=(None,'_new_fromTranscripts'))
    not_in_schema = not_in_schema.merge(
        schema_df[['coordinates','gene_id']], how='left',
        on='coordinates', suffixes=(None,'_new_fromCoordinates'))
    not_in_schema.drop(columns=['transcript_id','coordinates'], inplace=True)
    not_in_schema.fillna('', inplace=True)
    not_in_schema['new_ENSG'] = not_in_schema.apply(lambda x: find_map(x), axis=1)
    not_in_schema = not_in_schema[not_in_schema['new_ENSG'] != '']

    not_in_schema = not_in_schema[not_in_schema['new_ENSG'].isin(input_df['gene_id']) == False]

    not_in_schema = not_in_schema.merge(schema_df[['gene_id','gene_name']],
                        how='left', left_on='new_ENSG', right_on='gene_id',
                        suffixes=(None,'_y')).drop(columns=['gene_id_y'])

    return not_in_schema


args = getArgs()
if not args.version:
    sys.exit('ERROR: --version is required')
schema_v = args.version

schema_ref = f'gencode.v{schema_v}.primary_assembly.annotation.gtf.gz'
if not os.path.isfile(schema_ref):
    sys.exit(f'ERROR: {schema_ref} missing')

refs = []
for n in range(19,22):
    refs.append(f'gencode.v{n}.chr_patch_hapl_scaff.annotation.gtf.gz')
for n in range(22, int(schema_v)):
    refs.append(f'gencode.v{n}.primary_assembly.annotation.gtf.gz')

missing_files = [r for r in refs if r not in os.listdir()]
if missing_files:
    print(f'ERROR: {len(missing_files)} reference files missing')
    print('\n'.join(missing_files))
    sys.exit()

schema_df = create_df(schema_ref)

all_inputs = pd.DataFrame()
ref_comp = pd.DataFrame()

#compare each version to the target version
for r in refs:
    print(r)
    if 'v19' in r:
        input_df = create_df(r, hg19=True)
    else:
        input_df = create_df(r)
    input_df['input'] = r
    all_inputs = pd.concat([all_inputs, input_df])

    comparison_df = compare(schema_df, input_df)
    ref_comp = pd.concat([ref_comp, comparison_df])

min_report = ref_comp[['gene_id','input','new_ENSG']]

min_report['gene_id-input'] = min_report.apply(lambda x: x['gene_id'] + '-' + x['input'], axis=1)
min_report['new_ENSG-input'] = min_report.apply(lambda x: x['new_ENSG'] + '-' + x['input'], axis=1)
min_report['gene_id-new_ENSG'] = min_report.apply(lambda x: x['gene_id'] + '-' + x['new_ENSG'], axis=1)

min_report['input'] = min_report['input'].apply(
    lambda x: x.replace('gencode.','').replace('.annotation.gtf.gz','')
        .replace('.chr_patch_hapl_scaff','').replace('.primary_assembly','')
    )
#remove any where two IDs map to the same target ID in the same old version
invalid_mappings = min_report[min_report.duplicated(subset='new_ENSG-input',keep=False)]['gene_id-new_ENSG'].unique()
min_report = min_report[min_report['gene_id-new_ENSG'].isin(invalid_mappings) == False]

#remove any where two IDs map to the same target ID amongst all old versions
invalid_new_ENSG = []
all_new = [i.split('-')[-1] for i in min_report['gene_id-new_ENSG'].unique()]
dup_new = set([i for i in all_new if all_new.count(i) > 1])
for i in dup_new:
    gene_ids = min_report[min_report['new_ENSG'] == i]['gene_id'].unique()
    input_for_gene_ids = all_inputs[all_inputs['gene_id'].isin(gene_ids)]
    if len(input_for_gene_ids['input'].unique()) < len(input_for_gene_ids['input']):
        invalid_new_ENSG.append(i)
min_report = min_report[min_report['new_ENSG'].isin(invalid_new_ENSG) == False]

min_all_inputs = all_inputs[all_inputs['gene_id'].isin(min_report['gene_id'].unique())][['gene_id','input']]
min_all_inputs['gene_id-input'] = min_all_inputs.apply(lambda x: x['gene_id'] + '-' + x['input'], axis=1)
min_all_inputs = min_all_inputs[min_all_inputs['gene_id-input'].isin(min_report['gene_id-input']) == False]

min_report = min_report.groupby('gene_id-new_ENSG').agg(','.join).reset_index()
min_report[['gene_id','new_ENSG']] = min_report['gene_id-new_ENSG'].str.split('-',expand=True)[[0,1]]
min_report = min_report[min_report.duplicated(subset='gene_id', keep=False) == False]

for i,row in min_all_inputs.iterrows():
    new_ENSG = min_report.loc[min_report['gene_id'] == row['gene_id'], 'new_ENSG'].iloc[0]
    full_df = all_inputs[all_inputs['input'] == row['input']]
    if new_ENSG in full_df['gene_id'].unique():
        min_all_inputs.loc[i, 'new_ENSG in version'] = row['input']
        
min_all_inputs['input'] = min_all_inputs['input'].apply(
    lambda x: x.replace('gencode.','').replace('.annotation.gtf.gz','')
        .replace('.chr_patch_hapl_scaff','').replace('.primary_assembly','')
    )

min_all_inputs.fillna('', inplace=True)
min_all_inputs['new_ENSG in version'] = min_all_inputs['new_ENSG in version'].apply(
    lambda x: x.replace('gencode.','').replace('.annotation.gtf.gz','')
        .replace('.chr_patch_hapl_scaff','').replace('.primary_assembly','')
    )

min_all_inputs = min_all_inputs.groupby('gene_id').agg(lambda x: ','.join([e for e in x if e != ''])).reset_index()
min_all_inputs['new_ENSG in version'] = min_all_inputs['new_ENSG in version'].apply(lambda x: '' if x and x[0] == ',' else x)
min_all_inputs.rename(columns={'input':'present, not mapped'}, inplace=True)

min_report.rename(columns={'input':'mapped in'}, inplace=True)
min_report = min_report[['gene_id','mapped in','new_ENSG']]

min_report = min_report.merge(schema_df[['gene_id','gene_name']],
                    how='left', left_on='new_ENSG', right_on='gene_id',
                    suffixes=(None,'_y'))
min_report = min_report.merge(min_all_inputs, how='left', on='gene_id')
min_report.fillna('', inplace=True)

min_report = min_report[min_report['new_ENSG in version'] == '']

with open(f'gene_map_v{schema_v}.json', 'w') as f:
    json.dump(
        min_report.set_index('gene_id').to_dict()['new_ENSG'],
        f,
        indent=4
    )
