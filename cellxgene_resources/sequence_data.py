import ftplib
import json
import os
import pandas as pd
import requests
import tarfile
import xml.etree.ElementTree as ET


base_urls = {
    'explore.data.humancellatlas.org/projects': 'hca',
    'ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GS': 'geo',
    'ncbi.nlm.nih.gov/projects/gap': 'dbgap',
    'ncbi.nlm.nih.gov/bioproject/?term=PRJ': 'bioproj',
    'ega-archive.org': 'ega',
    'ebi.ac.uk/ena/browser/view/': 'ena',
    'ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB': 'arrex',
    'assets.nemoarchive.org/': 'nemo'
}


def parse_url(url):
    for k,v in base_urls.items():
        if k in url:
            return v
    return 'other'


def validate_raw_ega(url):
    #https://metadata.ega-archive.org/spec
    raw_data_formats = ['fastq.gz','bam','cram']

    acc = url.split('/')[-1]
    obj_type = url.split('/')[-2]
    api_base = 'https://metadata.ega-archive.org'

    if obj_type == 'studies':
        ds_query = f'{api_base}/studies/{acc}/datasets?limit=100000'
        r = requests.get(ds_query).json()
        datasets = [d['accession_id'] for d in r]
    else:
        datasets = [acc]

    for d in datasets:
        files_query = f'{api_base}/datasets/{d}/files?limit=100000'
        r = requests.get(files_query).json()
        for f in r:
            if f['extension'] in raw_data_formats:
                return True

    return False


def validate_raw_nemo(url):
    raw_data_formats = ['fastq.tar']

    for df in pd.read_html(url):
        if 'Dataset Collection URL' in df['Field'].unique():
            coll_url = df.loc[df['Field'] == 'Dataset Collection URL']['Value'].iloc[0]
            if str(coll_url).endswith('.tgz'):
                r = requests.get(coll_url)
                with open('temp.tgz','wb') as f:
                    f.write(r.content)
                file_list = []
                tar = tarfile.open('temp.tgz', 'r:gz')
                for item in tar:
                    if item.name.endswith('fetch.txt'):
                        tar.extract(item)
                        with open(item.name, 'r') as f:
                            for line in (f.read().split('\n')):
                                file_list.append(line.split('\t')[0])
                        os.remove(item.name)
                        os.rmdir(item.name.split('/')[0])
                        os.remove('temp.tgz')
                raw_files = [f for f in file_list if f.endswith(tuple(raw_data_formats))]
                if raw_files:
                    return True
        else:
            i = df.loc[df['Field'] == 'Identifier']['Value'].iloc[0].split(':')[1]
            #ATTN - can we just return validate_raw_nemo('https://assets.nemoarchive.org/' + i)?
            raw_present = validate_raw_nemo('https://assets.nemoarchive.org/' + i)
            if raw_present == True:
                return True

    return False


def validate_raw_arrex(url):
    acc = url.split('/')[-1]
    api_base = 'https://www.ebi.ac.uk/biostudies/api/v1'
    q_url = f'{api_base}/studies/{acc}/info'
    r = requests.get(q_url).json()
    ftp_link = r['ftpLink']

    ftp = ftplib.FTP('ftp.ebi.ac.uk', 'anonymous', 'anonymous@')
    ftp.cwd(ftp_link.replace('ftp://ftp.ebi.ac.uk','') + '/Files')

    filename = f'{acc}.sdrf.txt'
    with open(filename, 'wb') as file:
        ftp.retrbinary(f'RETR {filename}', file.write)

    erxs = []
    df = pd.read_csv(filename, sep='\t')
    for c in df.columns:
        if '[ENA_EXPERIMENT]' in c:
            erxs.extend(df[c].dropna().unique())
    os.remove(filename)
    ftp.quit()

    if erxs:
        return True

    return False


def validate_raw_hca(url):
    raw_data_formats = ['fastq.gz','fastq','fq.gz']

    api_base = 'https://service.azul.data.humancellatlas.org' 
    pj_id = url.split('/')[-1]
    query = {
        'projectId': {'is': [pj_id]}
    }
    q_url = f'{api_base}/index/files/?filters={json.dumps(query)}&size=250'
    r = requests.get(q_url, headers={'Content-Type': 'application/json'}).json()
    hits = r['hits']
    while r['pagination']['next']:
        next_endpoint = r['pagination']['next']
        r = requests.get(next_endpoint).json()
        hits.extend(r['hits'])
    for h in hits:
        for f in h['files']:
            if f['format'] in raw_data_formats:
                return True

    return False


def validate_raw_ncbi(url):
    raw_data_formats = ['fastq','TenX','bam']

    acc = url.split('/')[-1]
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    esearch_base = f'{eutils_base}esearch.fcgi'
    efetch_base = f'{eutils_base}efetch.fcgi'
    prj_flag = False
    if acc.startswith('GS'):
        url1 = f'{esearch_base}?db=bioproject&term={acc}[Project Accession]&retmode=json'
        r1 = requests.get(url1).json()
        if r1['esearchresult']['idlist']:
            i = r1['esearchresult']['idlist'][0] #list of ids, ideally - only search entry type:Series
            url2 = f'{efetch_base}?db=bioproject&id={i}'
            r2 = requests.get(url2)
            rXml = ET.fromstring(r2.text)
            for a in rXml.iter('ArchiveID'):
                prj = a.attrib['accession']
                prj_flag = True
    else:
        prj = acc
        prj_flag = True

    if not prj_flag:
        url5 = f'{esearch_base}?db=gds&term={acc}[GEO Accession]&retmode=json'
        r5 = requests.get(url5).json()
        if r5['esearchresult']['idlist']:
            i = r5['esearchresult']['idlist'][0] #list of ids, ideally - only search entry type:Series
            url6 = f'{efetch_base}?db=gds&id={i}'
            r6 = requests.get(url6)
            for line in r6.text.split('\n'):
                if line.startswith('SRA Run Selector:'):
                    prj_flag = True
                    prj = line.split('acc=')[-1]

    if prj_flag:
        url3 = f'{esearch_base}?db=sra&term={prj}&retmode=json&retmax=100000'
        r3 = requests.get(url3).json()
        idlist = r3['esearchresult']['idlist']
        sublists = [idlist[i:i+500] for i in range(0, len(idlist), 500)]
        for sub in sublists:
            ids = ','.join(sub)
            url4 = f'{efetch_base}?db=sra&id={ids}'
            r4 = requests.get(url4)
            rXml = ET.fromstring(r4.text)
            for ep in rXml.iter('EXPERIMENT_PACKAGE'):
                for run in ep.iter('RUN'):
                    for cf in run.iter('CloudFile'):
                        if cf.attrib['filetype'] in raw_data_formats:
                            return True

    return False


def detect_sequence_data(url):
    resource = parse_url(url)

    if resource in ['geo','bioproj','ena','dbgap']:
        raw_present = validate_raw_ncbi(url)
    elif resource == 'hca':
        raw_present = validate_raw_hca(url)
    elif resource == 'arrex':
        raw_present = validate_raw_arrex(url)
    elif resource == 'nemo':
        raw_present = validate_raw_nemo(url)
    elif resource == 'ega':
        raw_present = validate_raw_ega(url)
    else: #ATTN - what to do w/ undetermined?
        raw_present = 'undetermined'

    return raw_present
