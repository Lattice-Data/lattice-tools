import re
import requests


CLEANR = re.compile('<.*?>') 


def cleanhtml(raw_html):
    cleantext = re.sub(CLEANR, '', raw_html).replace('\n','')
    cleantext = ' '.join(cleantext.split()).replace('‚Äì','-')
    return cleantext


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

    out['Title'] = cleanhtml(r['message']['title'][0])

    if r['message']['container-title']:
        out['Journal'] = r['message']['container-title'][0]
    elif 'institution' in r['message']:
        out['Journal'] = r['message']['institution'][0]['name']
    else:
        out['Journal'] = r['message'].get('group-title')

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
        pmid_dois.update([h.get('doi','') for h in r['results']])

    return dois
