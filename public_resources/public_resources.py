import ftplib
import json
import os
import pandas as pd
import re
import requests
import tarfile
import xml.etree.ElementTree as ET


def cleanhtml(raw_html):
    cleantext = re.sub(re.compile('<.*?>'), '', raw_html).replace('\n','')
    cleantext = ' '.join(cleantext.split()).replace('‚Äì','-')

    return cleantext


def query_preprints(doi, jour):
    r = requests.get(f'https://api.biorxiv.org/details/{jour}/{doi}/na/json')
    r = r.json()
    for c in r['collection']:
        if c['published'] != 'NA':
            return c['published']


def get_journal(r):
    if r['message']['container-title']:
        return r['message']['container-title'][0]
    elif 'institution' in r['message']:
        return r['message']['institution'][0]['name']
    else:
        return r['message'].get('group-title')


def doi_checker(doi):
    out = {'DOI': doi}

    r = requests.get('https://api.crossref.org/works/' + str(doi))
    if r.status_code != 200:
        out['invalid DOI'] = 'yes'
        return out
    r = r.json()

    if r['message']['DOI'].lower() != doi.lower():
        out['updated_doi'] = r['message']['DOI']
        out['invalid DOI'] = 'format'

    if r['message']['relation']:
        if 'is-preprint-of' in r['message']['relation']:
            doi = r['message']['relation']['is-preprint-of'][0]['id']
            out['updated_doi'] = doi
            out['invalid DOI'] = 'outdated'
            r = requests.get('https://api.crossref.org/works/' + str(doi)).json()

    out['Journal'] = get_journal(r)

    #in case bioRxiv API has a published DOI that crossref does not
    if out.get('Journal') in ['bioRxiv','medRxiv']:
        jour = out['Journal'].lower()
        if doi := query_preprints(doi, jour):
            out['updated_doi'] = doi
            out['invalid DOI'] = 'outdated'
            r = requests.get('https://api.crossref.org/works/' + str(doi)).json()
            out['Journal'] = get_journal(r)

    out['Title'] = cleanhtml(r['message']['title'][0])

    if 'published' in r['message']:
        out['Year'] = r['message']['published']['date-parts'][0][0]

    first_auths = []
    for a in r['message'].get('author',[]):
        if a['sequence'] == 'first':
            if 'name' in a:
                first_auths.append((a['name']))
            elif 'given' in a:
                first_auths.append((a['given'] + ' ' + a['family']))
            else:
                first_auths.append((a['family']))
    out['First authors'] = ','.join(first_auths)
    
    return out


def pubtator_search(search_term):
    url = f'https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text="{search_term}"'
    r = requests.get(url).json()
    dois = [h.get('doi','') for h in r['results']]
    for i in range(2,r['total_pages']+1):
        time.sleep(0.5)
        url = f'https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text="{search_term}"&page={i}'
        r = requests.get(url).json()
        dois.update([h.get('doi','') for h in r['results']])

    return dois


def pmid_to_doi(pmid):
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    url = f'{eutils_base}efetch.fcgi?db=pubmed&id={pmid}'
    r = requests.get(url)
    responseXml = ET.fromstring(r.text)
    for pmdata in responseXml.iter('PubmedData'):
        for ai in pmdata.find('ArticleIdList').iter('ArticleId'):
            if ai.attrib['IdType'] == 'doi':
                return ai.text


def pubmed_search(search_term):
    dois = []
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    url = f'{eutils_base}esearch.fcgi?db=pubmed&term={search_term}&retmode=json&retmax=2000'
    r = requests.get(url).json()
    for i in r['esearchresult']['idlist']:
        doi = pmid_to_doi(i)
        dois.append(doi)

    return dois


def publication_search(search_term):
    dois = pubtator_search(search_term)
    dois.extend(pubmed_search(search_term))

    return dois


def validate_raw_ega(url):
    #https://metadata.ega-archive.org/spec
    formats = set()
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
        fs = [f['extension'] for f in r if f['extension'] in raw_data_formats]
        formats.update(fs)

    return list(formats)


def validate_raw_nemo(url):
    raw_data_formats = ['fastq.tar']

    for df in pd.read_html(url):
        if 'Field' in df.columns:
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
                raw_present = validate_raw_nemo('https://assets.nemoarchive.org/' + i)
                if raw_present == True:
                    return True
        elif 'Name' in df.columns:
            raw_files = [f for f in df['Name'].unique() if str(f).endswith(tuple(raw_data_formats))]
            if raw_files:
                return True

    return False


def validate_raw_hca(url):
    formats = set()
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
        fs = [f['format'] for f in h['files'] if f['format'] in raw_data_formats]
        formats.update(fs)

    return list(formats)


def validate_raw_insdc(url, arrex=False):
    raw_data_formats = ['fastq','TenX','bam','cram','10X Genomics bam file']

    acc = url.split('/')[-1].split('=')[-1]
    insdc_attrs = insdc_meta(acc, arrex)
    df = pd.DataFrame(insdc_attrs)
    formats = [f for c in df.columns for f in df[c].unique() if 'semantic' in c]

    return list(set([f for f in formats if f in raw_data_formats]))


def validate_raw_dbgap(url):
    formats = set()
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    esearch_base = f'{eutils_base}esearch.fcgi'
    raw_data_formats = ['fastq','TenX','bam','cram','10X Genomics bam file']
    acc = url.split('/')[-1].split('=')[-1]
    if '.' in acc:
        new_acc = acc.split('.')[0]
        errors = f'remove version from dbGaP URL: {url.replace(acc, new_acc)}'
        acc = new_acc

    url3 = f'{esearch_base}?db=sra&term={acc}&retmode=json&retmax=100'
    r3 = requests.get(url3).json()
    idlist = r3['esearchresult']['idlist']
    sublists = [idlist[i:i+500] for i in range(0, len(idlist), 500)]
    for sub in sublists:
        ids = ','.join(sub)
        url4 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={ids}'
        r4 = requests.get(url4)
        responseXml = ET.fromstring(r4.text)
        for ep in responseXml.iter('EXPERIMENT_PACKAGE'):
            formats.update([f.attrib['filetype'] for f in ep.iter('CloudFile')])

    raw_formats = list(set([f for f in formats if f in raw_data_formats]))

    if errors:
        return raw_formats, errors

    return raw_formats


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

    raw_files = []
    df = pd.read_csv(filename, sep='\t')
    if 'Comment[FASTQ_URI]' in df.columns:
        raw_files = [f for f in df['Comment[FASTQ_URI]'] if f.endswith(tuple(raw_data_formats))]
    os.remove(filename)
    ftp.quit()

    if raw_files:
        return True

    return False


def insdc_meta(acc, arrex=False):
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    esearch_base = f'{eutils_base}esearch.fcgi'
    efetch_base = f'{eutils_base}efetch.fcgi'
    prj = []
    if acc.startswith('GS'):
        url1 = f'{esearch_base}?db=bioproject&term={acc}[Project Accession]&retmode=json'
        r1 = requests.get(url1).json()
        if r1['esearchresult']['idlist']:
            i = r1['esearchresult']['idlist'][0]
            url2 = f'{efetch_base}?db=bioproject&id={i}'
            r2 = requests.get(url2)
            rXml = ET.fromstring(r2.text)
            for a in rXml.iter('ArchiveID'):
                prj = [a.attrib['accession']]
        if not prj:
            url5 = f'{esearch_base}?db=gds&term={acc}[Accn]&retmode=json&retmax=100'
            r5 = requests.get(url5).json()
            if r5['esearchresult']['idlist']:
                ids = ','.join(r5['esearchresult']['idlist'])
                url6 = f'{efetch_base}?db=gds&id={ids}'
                r6 = requests.get(url6)
                for line in r6.text.split('\n'):
                    if line.startswith('SRA Run Selector:'):
                        p = line.split('acc=')[-1]
                        prj.append(p)
    else:
        prj = [acc]

    if prj:
        attributes = []
        idlist = set()
        for p in prj:
            url3 = f'{esearch_base}?db=sra&term={p}&retmode=json&retmax=100'
            r3 = requests.get(url3).json()
            idlist.update(r3['esearchresult']['idlist'])
        idlist = list(idlist)
        sublists = [idlist[i:i+500] for i in range(0, len(idlist), 500)]
        for sub in sublists:
            ids = ','.join(sub)
            url4 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={ids}'
            r4 = requests.get(url4)
            #parse the records for needed information & write report
            responseXml = ET.fromstring(r4.text)
            for ep in responseXml.iter('EXPERIMENT_PACKAGE'):
                if arrex:
                    s = ep.find('STUDY')
                    if s.attrib['alias'] != prj[0]:
                        continue
                e = ep.find('EXPERIMENT')
                attr1 = {'Experiment:title': e.find('TITLE').text}
                for k,v in e.attrib.items():
                    attr1['Experiment:' + k] = v
                for p in ['LIBRARY_STRATEGY','LIBRARY_SOURCE','LIBRARY_SELECTION']:
                    attr1['Experiment:' + p.lower()] = e.find('DESIGN').find('LIBRARY_DESCRIPTOR').find(p).text
                for s in ep.iter('SAMPLE'):
                    attr1['Sample:primary_id'] = s.find('IDENTIFIERS').find('PRIMARY_ID').text
                    if s.find('title'):
                        attr1['Sample:title'] = s.find('TITLE').text
                    for ei in s.iter('EXTERNAL_ID'):
                        attr1['Sample:' + ei.attrib['namespace'] + '_id'] = ei.text
                    attr1['Sample:taxon_id'] = s.find('SAMPLE_NAME').find('TAXON_ID').text
                    attr1['Sample:scientific_name'] = s.find('SAMPLE_NAME').find('SCIENTIFIC_NAME').text
                    for sa in s.iter('SAMPLE_ATTRIBUTE'):
                        attr1['Sample:' + sa.find('TAG').text] = sa.find('VALUE').text
                for run in ep.iter('RUN'):
                    attr2 = attr1.copy()
                    for k,v in run.attrib.items():
                        attr2['Run:' + k] = v
                    file_count = 0
                    for f in run.iter('SRAFile'):
                        if 'SRA' not in f.attrib['semantic_name']:
                            file_count += 1
                            for k,v in f.attrib.items():
                                attr2['Run:file_' + str(file_count) + ' ' + k] = v
                            for alt in f.iter('Alternatives'):
                                attr2['Run:file_' + str(file_count) + ' url'] = alt.attrib['url']
                    attributes.append(attr2)

        return attributes


data_repo_bases = {
    'explore.data.humancellatlas.org': 'hca',
    'ncbi.nlm.nih.gov/geo': 'geo',
    'ncbi.nlm.nih.gov/projects/gap': 'dbgap',
    'ncbi.nlm.nih.gov/bioproject': 'bioproj',
    'ega-archive.org': 'ega',
    'ebi.ac.uk/ena/browser/view': 'ena',
    'ebi.ac.uk/biostudies/arrayexpress': 'arrex',
    'nemoarchive.org': 'nemo',
    'ngdc.cncb.ac.cn': 'ngdc',
    'synapse.org': 'synapse',
    'zenodo.org': 'zenodo',
    'registry.opendata.aws': 'aws',
    'duos': 'duos'
}


def parse_data_repo_url(url):
    for k,v in data_repo_bases.items():
        if k in url:
            return v
    return 'other'


def detect_sequence_data(url):
    resource = parse_data_repo_url(url)

    if resource in ['geo','bioproj','ena']:
        raw_present = validate_raw_insdc(url)
    elif resource == 'arrex':
        raw_present = validate_raw_insdc(url, arrex=True)
        if not raw_present:
            raw_present = validate_raw_arrex(url)
    elif resource == 'dbgap':
        raw_present = validate_raw_dbgap(url)
    elif resource == 'hca':
        raw_present = validate_raw_hca(url)
    elif resource == 'nemo':
        raw_present = validate_raw_nemo(url)
    elif resource == 'ega':
        raw_present = validate_raw_ega(url)
    else:
        raw_present = 'undetermined'

    return raw_present
