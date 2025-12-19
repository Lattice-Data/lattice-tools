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
        print('dbGaP URL with version included')
        acc = acc.split('.')[0]

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

    return list(set([f for f in formats if f in raw_data_formats]))


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


def figshare_pick_acc(x):
    if pd.notnull(x['url_public_html']):
        return x['url_public_html']
    if x['doi'] != x['resource_doi_y']:
        return 'https://doi.org/' + x['doi']


def figshare_pick_doi(x):
    if x['doi'] != x['resource_doi_y']:
        return x['resource_doi_y']
    for p in ['doi','resource_doi_x','resource_doi_y']:
        if x[p]:
            return x[p]


def figshare_query(search_term):
    api_base = 'https://api.figshare.com/v2'
    url = f'{api_base}/collections'
    query = {'search_for': search_term, 'limit': 1000}
    coll_r = requests.get(url, params=query).json()
    if coll_r:
        for c in coll_r:
            full_c = requests.get(c['url']).json()
            c['resource_doi'] = full_c['resource_doi']
        fs_coll_df = pd.DataFrame(coll_r)[['title','doi','resource_doi']]

    url = f'{api_base}/item_types'
    r = requests.get(url).json()
    for it in r:
        if it['name'] == 'dataset':
            ds_item = it['id']
    query['item_type'] = ds_item
    url = f'{api_base}/articles'
    dataset_r = requests.get(url, params=query).json()
    if dataset_r:
        for ds in dataset_r:
            file_r = requests.get(ds['url']).json()
            if 'files' in file_r:
                formats = set()
                for f in file_r['files']:
                    path = f['name'].split('.')
                    form = path[-1] if path[-1] != 'gz' else path[-2]
                    formats.add(form.upper())
                ds['file formats'] = ','.join(formats)
        fs_ds_df = pd.DataFrame(dataset_r)[['resource_doi','group_id','url_public_html','file formats','title']]

        no_group_df = fs_ds_df[fs_ds_df['group_id'].isna()]
        fs_ds_df = fs_ds_df[fs_ds_df['group_id'].isna() == False]
        fs_ds_df['groupby_field'] = fs_ds_df.apply(lambda x: f'{x["group_id"]}_{x["resource_doi"]}', axis=1)
        fs_ds_df = fs_ds_df.groupby('groupby_field', as_index=False).agg(lambda x: ','.join(set(x)))
        fs_ds_df = pd.concat([fs_ds_df,no_group_df]).drop(columns=['groupby_field','group_id'])
        fshare_df = fs_ds_df
        if coll_r:
            fshare_df = fshare_df.merge(fs_coll_df, left_on='resource_doi', right_on='doi', how='outer')
        else:
            print(search_term,'no coll but yes ds')
    elif coll_r:
        print(search_term,'no ds but yes coll')
        fshare_df = fs_coll_df
    else:
        print(search_term,'none of either')
        return pd.DataFrame()

    fshare_df['accession'] = fshare_df.apply(lambda x: figshare_pick_acc(x), axis=1)
    fshare_df['doi'] = fshare_df.apply(lambda x: figshare_pick_doi(x), axis=1)
    fshare_df['title'] = fshare_df.apply(lambda x: x['title_x'] if pd.notnull(x['title_x']) else x['title_y'], axis=1)
    fshare_df.fillna('', inplace=True)
    fshare_df['file formats'] = fshare_df['file formats'].apply(lambda x: x.split(','))
    fshare_df = fshare_df[['accession','title','doi','file formats']]

    return fshare_df


def list_recursive_ftp(ftp, files):
    for f in ftp.nlst():
        try:
            ftp.cwd(f)
        except ftplib.error_perm:
            files.append(f)
        else:
            list_recursive_ftp(ftp, files)
            ftp.cwd('..')


def biostudies_query(search_term):
    reports = []
    api_base = 'https://www.ebi.ac.uk/biostudies/api/v1'
    queries = [
        f'query="{search_term}"'
    ]
    url = f'{api_base}/search?{"&".join(queries)}&pageSize=100'
    search_r = requests.get(url).json()
    if search_r['totalHits'] == 0:
        return pd.DataFrame()
    accessions = [h['accession'] for h in search_r['hits']]

    for acc in accessions:
        rep = {'accession': acc, 'doi': [], 'pmid': []}
        standard_files = [f'{acc}.{ext}' for ext in ['idf.txt','sdrf.txt','json','tsv','xml']]

        url = f'{api_base}/studies/{acc}'
        r = requests.get(url).json()
        if not r.get('errorMessage') == 'Study not found':
            for ss in r['section']['subsections']:
                if isinstance(ss, dict):
                    if ss['type'] == 'Publication':
                        for a in ss['attributes']:
                            if a['name'] == 'DOI':
                                rep['doi'].append(a['value'])
                            elif a['name'] == 'Pubmed ID':
                                rep['pmid'].append(a['value'])
                    for l in ss.get('links', []):
                        if isinstance(l, list):
                            for ll in l:
                                if ll['attributes'][0]['value'] == 'DOI':
                                    rep['doi'].append(ll['url'])
                                if ll['attributes'][0]['value'] == 'PMID':
                                    rep['pmid'].append(ll['url'])
                        else:
                            if l['attributes'][0]['value'] == 'DOI':
                                rep['doi'].append(l['url'])
                            if l['attributes'][0]['value'] == 'PMID':
                                rep['pmid'].append(l['url'])
            for a in r['attributes']:
                if a['name'] == 'Title':
                    rep['title'] = a['value']

            url = f'{api_base}/studies/{acc}/info'
            r = requests.get(url).json()
            ftp_link = r['ftpLink']

            ftp = ftplib.FTP('ftp.ebi.ac.uk', 'anonymous', 'anonymous@')
            ftp.cwd(ftp_link.replace('ftp://ftp.ebi.ac.uk',''))

            files = []
            list_recursive_ftp(ftp, files)
            files = [f for f in files if f not in standard_files]

            formats = set()
            for f in files:
                path = f.split('.')
                form = path[-1] if path[-1] != 'gz' else path[-2]
                formats.add(form.upper())
            rep['file formats'] = list(formats)

            if 'Files' in ftp.nlst():
                ftp.cwd('Files')
                filename = f'{acc}.sdrf.txt'
                if filename in ftp.nlst():
                    with open(filename, 'wb') as file:
                        ftp.retrbinary(f'RETR {filename}', file.write)
                    temp_df = pd.read_csv(filename, sep='\t')
                    for c in ['developmental stage','organism part','disease','organism']:
                        if f'Characteristics[{c}]' in temp_df.columns:
                            rep[c] = ','.join([str(e) for e in temp_df[f'Characteristics[{c}]'].unique()])
                    os.remove(filename)

            reports.append(rep)

    ftp.quit()
    arrex_df = pd.DataFrame(reports)
    for c in ['doi','pmid']:
        arrex_df[c] = arrex_df[c].apply(lambda x: ','.join(set(x)))

    return arrex_df


def geo_query(search_term):
    pubmeds = []
    reports = []

    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

    url = f'{eutils_base}esearch.fcgi?db=gds&term="{search_term}"&retmode=json&retmax=2000'
    r = requests.get(url).json()
    if 'errorlist' in r['esearchresult']:
        if search_term in r['esearchresult']['errorlist']['phrasesnotfound']:
            return pd.DataFrame()

    ids = r['esearchresult']['idlist']
    for i in range(0, len(ids) , 250):
        chunk_ids = ids[i:i+250]
        url = f'{eutils_base}efetch.fcgi?db=gds&id={chunk_ids}'
        r = requests.get(url)
        for record in r.text.split('\n\n'):
            rep = {}
            for line in record.split('\n'):
                if line.startswith('Organism:'):
                    orgs = line.split('\t')[1].replace('; ',',')
                    rep['organism'] = orgs
                elif line.startswith('Series'):
                    rep['accession'] = line.split('Accession: ')[-1].split('\t')[0].strip(' ')
                elif line.startswith('Platform:') and 'Series:' in line:
                    rep['accession'] = line.split('Series: ')[-1].strip(' ')
                elif line.startswith('FTP download:'):
                    if '(' in line:
                        fformats = line.split('(')[1].split(')')[0]
                        rep['file formats'] = fformats.replace(' ','')
            reports.append(rep)

        url = f'{eutils_base}esummary.fcgi?db=gds&id={",".join(chunk_ids)}&retmode=json'
        r = requests.get(url).json()
        for k,v in r['result'].items():
            if k != 'uids' and v['pubmedids']:
                pubmeds.append({'accession': ','.join([v['accession']]), 'pmid': ','.join(v['pubmedids'])})

    geo_df = pd.DataFrame(reports)
    geo_df = geo_df.merge(pd.DataFrame(pubmeds), on='accession', how='left')
    geo_df.drop_duplicates(inplace=True)
    geo_df.fillna('', inplace=True)

    new_rows = []
    for i, row in geo_df.iterrows():
        for a in row['accession'].split(' '):
            new_row = {'accession': a}
            for c in ['organism','file formats','pmid']:
                new_row[c] = row[c]
            new_rows.append(new_row)
    geo_df = pd.DataFrame(new_rows)
    geo_df = geo_df.groupby('accession', as_index=False).agg(','.join)
    geo_df['pmid'] = geo_df['pmid'].apply(lambda x: ','.join([e for e in set(x.split(',')) if e]))

    dups = geo_df[(geo_df.duplicated(subset='pmid', keep=False) == True) & (geo_df['pmid'] != '')]
    dups = dups.groupby('pmid', as_index=False).agg(','.join)
    geo_df = geo_df[(geo_df.duplicated(subset='pmid', keep=False) == False) | (geo_df['pmid'] == '')]
    geo_df = pd.concat([geo_df, dups])

    for c in ['organism']:
        geo_df[c] = geo_df[c].apply(lambda x: ','.join([e for e in set(x.split(',')) if e]))
    geo_df['file formats'] = geo_df['file formats'].apply(lambda x: [e for e in set(x.split(',')) if e])

    return geo_df


def zenodo_query(search_term):
    base_url = 'https://zenodo.org/api/'
    url = f'{base_url}/records'
    query = {'q': f'"{search_term}"', 'type': 'dataset', 'size': 1000}
    r = requests.get(url, params=query).json()
    reports = []
    for h in r['hits']['hits']:
        rep = {'accession': h['links']['self_html']}
        if 'zenodo' not in h['doi_url']:
            rep['doi'] = h['doi_url']
        rep['title'] = h['metadata']['title']
        formats = set()
        for f in h['files']:
            path = f['key'].split('.')
            form = path[-1] if path[-1] != 'gz' else path[-2]
            formats.add(form.upper())
        rep['file formats'] = list(formats)
        reports.append(rep)
    z_df = pd.DataFrame(reports)

    return z_df
