import io
import json
import random
import requests
import zipfile
from collections import defaultdict


configs = [
    {
        'field': 'assay',
        'ontology_name': 'EFO',
        'old_version': '3.82.0',
        'new_version': '3.90.0',
        'url_template': 'https://github.com/EBISPOT/efo/releases/download/v{version}/efo.json',
        'target_terms': [
            'http://www.ebi.ac.uk/efo/EFO_0002772',
            'http://www.ebi.ac.uk/efo/EFO_0010183'
        ]
    },
    {
        'field': 'tissue',
        'ontology_name': 'UBERON',
        'old_version': '2025-08-15',
        'new_version': '2026-04-01',
        'url_template': 'https://github.com/obophenotype/uberon/releases/download/v{version}/uberon-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/UBERON_0001062'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'UBERON',
        'old_version': '2025-08-15',
        'new_version': '2026-04-01',
        'url_template': 'https://github.com/obophenotype/uberon/releases/download/v{version}/uberon-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/UBERON_0000105'
        ]
    },
    {
        'field': 'cell_type',
        'ontology_name': 'CL',
        'old_version': '2025-07-30',
        'new_version': '2026-03-26',
        'url_template': 'https://github.com/obophenotype/cell-ontology/releases/download/v{version}/cl-base.json'
    },
    {
        'field': 'disease',
        'ontology_name': 'MONDO',
        'old_version': '2025-09-02',
        'new_version': '2026-05-05',
        'url_template': 'https://github.com/monarch-initiative/mondo/releases/download/v{version}/mondo-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/MONDO_0000001',
            'http://purl.obolibrary.org/obo/MONDO_0021178'
        ]
    },
    {
        'field': 'cell_type',
        'ontology_name': 'FBbt',
        'old_version': '2025-08-07',
        'new_version': '2026-04-03',
        'url_template': 'https://github.com/FlyBase/drosophila-anatomy-developmental-ontology/releases/download/v{version}/fbbt-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/FBbt_00007002'
        ]
    },
    {
        'field': 'tissue',
        'ontology_name': 'FBbt',
        'old_version': '2025-08-07',
        'new_version': '2026-04-03',
        'url_template': 'https://github.com/FlyBase/drosophila-anatomy-developmental-ontology/releases/download/v{version}/fbbt-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/FBbt_10000000'
        ],
        'exclude_terms': [
            'http://purl.obolibrary.org/obo/FBbt_00007002'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'FBdv',
        'old_version': '2025-05-29',
        'new_version': '2026-04-02',
        'url_template': 'https://github.com/FlyBase/drosophila-developmental-ontology/releases/download/v{version}/fbdv-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/FBdv_00007014',
            'http://purl.obolibrary.org/obo/FBdv_00005259'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'HsapDv',
        'old_version': '2025-01-23',
        'new_version': '2025-01-23',
        'url_template': 'https://github.com/obophenotype/developmental-stage-ontologies/releases/download/v{version}/hsapdv.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/HsapDv_0000001'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'MmusDv',
        'old_version': '2025-01-23',
        'new_version': '2025-01-23',
        'url_template': 'https://github.com/obophenotype/developmental-stage-ontologies/releases/download/v{version}/mmusdv.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/MmusDv_0000001'
        ]
    },
    {
        'field': 'self_reported_ethnicity',
        'ontology_name': 'HANCESTRO',
        'old_version': '2025-04-01',
        'new_version': '2025-10-14',
        'url_template': 'https://github.com/EBISPOT/hancestro/archive/refs/tags/v{version}.zip',
        'json_path': 'hancestro-{version}/hancestro-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/HANCESTRO_0601',
            'http://purl.obolibrary.org/obo/HANCESTRO_0602'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'WBls',
        'old_version': 'WS298',
        'new_version': '2026-04-16',
        'url_template': 'https://github.com/obophenotype/c-elegans-development-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'c-elegans-development-ontology-{version}/wbls-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/WBls_0000803',
            'http://purl.obolibrary.org/obo/WBls_0000804'
        ]
    },
    {
        'field': 'tissue',
        'ontology_name': 'ZFA',
        'old_version': '2025-09-05',
        'new_version': '2026-03-31',
        'url_template': 'https://github.com/ZFIN/zebrafish-anatomical-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'zebrafish-anatomical-ontology-{version}/zfa-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/ZFA_0100000'
        ],
        'exclude_terms': [
            'http://purl.obolibrary.org/obo/ZFA_0009000'
        ]
    },
    {
        'field': 'cell_type',
        'ontology_name': 'ZFA',
        'old_version': '2025-09-05',
        'new_version': '2026-03-31',
        'url_template': 'https://github.com/ZFIN/zebrafish-anatomical-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'zebrafish-anatomical-ontology-{version}/zfa-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/ZFA_0009000'
        ]
    },
    {
        'field': 'development_stage',
        'ontology_name': 'ZFS',
        'old_version': '2025-09-05',
        'new_version': '2026-03-31',
        'url_template': 'https://github.com/ZFIN/zebrafish-anatomical-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'zebrafish-anatomical-ontology-{version}/zfa.json', #can use zfs-base.json in future
        'target_terms': [
            'http://purl.obolibrary.org/obo/ZFS_0100000'
        ]
    },
    {
        'field': 'tissue',
        'ontology_name': 'WBbt',
        'old_version': '2025-08-18',
        'new_version': '2025-08-18',
        'url_template': 'https://github.com/obophenotype/c-elegans-gross-anatomy-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'c-elegans-gross-anatomy-ontology-{version}/wbbt-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/WBbt_0005766'
        ],
        'exclude_terms': [
            'http://purl.obolibrary.org/obo/WBbt_0004017',
            'http://purl.obolibrary.org/obo/WBbt_00006803'
        ]
    },
    {
        'field': 'cell_type',
        'ontology_name': 'WBbt',
        'old_version': '2025-08-18',
        'new_version': '2025-08-18',
        'url_template': 'https://github.com/obophenotype/c-elegans-gross-anatomy-ontology/archive/refs/tags/v{version}.zip',
        'json_path': 'c-elegans-gross-anatomy-ontology-{version}/wbbt-base.json',
        'target_terms': [
            'http://purl.obolibrary.org/obo/WBbt_0004017'
        ],
        'exclude_terms': [
            'http://purl.obolibrary.org/obo/WBbt_0006803'
        ]
    }
]

def format_term_id(term_url):
    term_id = term_url.split('/')[-1]
    return term_id.replace('_', ':')
    

def load_ontology_from_archive(ontology_name, version, archive_url, json_path):
    # Format the internal path with version
    internal_path = json_path.format(version=version)
    
    response = requests.get(archive_url)
    if response.status_code != 200:
        raise Exception(f"Failed to download archive: {response.status_code}")
    
    # Detect archive type from URL
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        with zf.open(internal_path) as f:
            return json.load(f)


def load_ontology_data(ontology_name, version, url_template, json_path):
    url = url_template.format(version=version)
    print(f'Loading {ontology_name} {version} from {url}...')
    
    if json_path:
        return load_ontology_from_archive(ontology_name, version, url, json_path)
    
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f'Failed to load {ontology_name} version {version}: {response.status_code}')


def get_all_children(edges, term, pred_type, ontology_name, exclude_terms=None):
    """Get all descendants of a term."""
    # Build edges dictionary for specified predicate
    pred_edges = {
        e['sub']: e['obj'] for e in edges 
        if e['pred'] == pred_type
        and e['sub'].split('/')[-1].split('_')[0] == ontology_name
    }
    
    # Reverse to get parent -> children mapping
    children = defaultdict(list)
    for child, parent in pred_edges.items():
        children[parent].append(child)
    
    # Recursive function to get all descendants
    def get_descendants(t):
        result = set()
        for child in children[t]:
            if exclude_terms and child in exclude_terms:
                continue  # Don't traverse this branch
            result.add(child)
            result.update(get_descendants(child))
        return result
    
    return get_descendants(term)


def get_combined_descendants(edges, target_terms, pred_type, ontology_name, exclude_terms=None):
    """Get combined descendants for multiple target terms."""
    all_descendants = set()
    for term in target_terms:
        all_descendants.update(get_all_children(edges, term, pred_type, ontology_name, exclude_terms))
    return all_descendants


def get_term_labels(nodes):
    """Extract term labels from nodes."""
    labels = {}
    for node in nodes:
        node_id = node.get('id')
        label = node.get('lbl') or node.get('label') or 'No label'
        if node_id:
            labels[node_id] = label
    return labels


def compare_ontology_versions(
    field,
    ontology_name,
    old_version,
    new_version,
    url_template,
    json_path=None,
    target_terms=None,
    exclude_terms=None,
    pred_type='is_a'
):
    if old_version == new_version:
        print(f'{ontology_name} no version update: {old_version}')
        return

    # Load both versions
    old_data = load_ontology_data(ontology_name, old_version, url_template, json_path)
    new_data = load_ontology_data(ontology_name, new_version, url_template, json_path)
    
    old_edges = old_data['graphs'][0]['edges']
    new_edges = new_data['graphs'][0]['edges']

    exclusion_set = set()
    if exclude_terms:
        for term in exclude_terms:
            exclusion_set.add(term)
            exclusion_set.update(get_all_children(old_edges, term, pred_type, ontology_name))
            exclusion_set.update(get_all_children(new_edges, term, pred_type, ontology_name))
        print(f'Excluding {len(exclusion_set)} terms from exclusion list')

    if target_terms:
        # Get descendants
        print(f'Getting descendants from {old_version}...')
        old_terms = get_combined_descendants(old_edges, target_terms, pred_type, ontology_name, exclusion_set)
        
        print(f'Getting descendants from {new_version}...\n')
        new_terms = get_combined_descendants(new_edges, target_terms, pred_type, ontology_name, exclusion_set)
    else:
        # Compare all terms instead of descendants
        old_terms = {
            node['id'] for node in old_data['graphs'][0]['nodes']
            if node['id'].split('/')[-1].split('_')[0] == ontology_name
            and node['id'] not in exclusion_set
        }
        new_terms = {
            node['id'] for node in new_data['graphs'][0]['nodes']
            if node['id'].split('/')[-1].split('_')[0] == ontology_name
            and node['id'] not in exclusion_set
        }
    print(f'\nTotal in {old_version}: {len(old_terms)}')
    print(f'Total in {new_version}: {len(new_terms)}')
    
    # Compare versions
    added = new_terms - old_terms
    removed = old_terms - new_terms
    unchanged = old_terms & new_terms
    
    # Get labels
    old_labels = get_term_labels(old_data['graphs'][0]['nodes'])
    new_labels = get_term_labels(new_data['graphs'][0]['nodes'])
    
    # Print results
    print(f'{ontology_name}-{field} COMPARISON: {old_version} vs {new_version}')
    print(f'{"="*80}')
    print(f'Unchanged: {len(unchanged)}')

    print(f'\nADDED in {new_version} (not in {old_version}): {len(added)}')
    if added:
        print(f'{"="*80}')
    
        lim = 10
        for term in sorted(added)[:lim]:
            label = new_labels.get(term, 'No label')
            print(f'  + {term}')
            print(f'    {label}')
        if len(added) > lim:
            print(f'  ... and {len(added) - lim} more')

    print(f'\nREMOVED in {new_version} (was in {old_version}): {len(removed)}')
    if removed:
        print(f'{"="*80}')
        for term in sorted(removed)[:lim]:
            label = old_labels.get(term, 'No label')
            print(f'  - {term}')
            print(f'    {label}')
        if len(removed) > lim:
            print(f'  ... and {len(removed) - lim} more')

    sample_added = random.sample(sorted(added), min(2, len(added))) if added else []
    sample_removed = random.sample(sorted(removed), min(2, len(removed))) if removed else []

    print(f'{"="*120}')
    print(f'{"="*120}')
    return {
        'sample_added': [format_term_id(t) for t in sample_added],
        'sample_removed': [format_term_id(t) for t in sample_removed]
    }


def compare_multiple_ontologies(ontology_configs, output_json):
    summary = {}

    for config in ontology_configs:
        try:
            results = compare_ontology_versions(**config)
            if results:
                f = config['field']
                if f not in summary:
                    summary[f] = {
                        'deprecated': [],
                        'new': []
                    }
                summary[f]['deprecated'].extend(results['sample_removed'])
                summary[f]['new'].extend(results['sample_added'])
                print('')
        except Exception as e:
            print(f'\nError processing {config["ontology_name"]}: {e}')

    with open(output_json, 'w') as f:
        json.dump(summary, f, indent=2)

compare_multiple_ontologies(configs, 'ontology_diff_summary.json')
