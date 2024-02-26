import re


req_present = {
    'uns': [
        'study_PI',
        'unpublished' #suggesting this is removed
    ],
    'obs': [
        'sample_ID',
        'library_ID',
        'institute',
        'library_preparation_batch',
        'library_sequencing_run',
        'alignment_software'
    ]
}

req_enum = {
    'obs': {
        'manner_of_death': [
            '0','1','2','3','4',
            'unknown',
            'not applicable'
        ],
        'sample_source': [
            'surgical donor',
            'postmortem donor',
            'living organ donor'
        ],
        'sample_collection_method': [
            'brush',
            'scraping',
            'biopsy',
            'surgical resection',
            'blood draw',
            'body fluid',
            'other'
        ],
        'sample_site_condition': [
            'healthy',
            'diseased',
            'adjacent'
        ],
        'sample_preservation_method': [
            'fresh', 'ambient temperature',
            'cut slide', 'paraffin block',
            'frozen at -70C', 'frozen at -80C', 'frozen at -150C',
            'frozen in liquid nitrogen', 'frozen in vapor phase',
            'RNAlater at 4C', 'RNAlater at 25C', 'RNAlater at -20C',
            'other'
        ],
        'sequenced_fragment': [
            '3 prime tag',
            '5 prime tag',
            'full length'
        ],
        'reference_genome': [
            'GRCh38', 'GRCh37',
            'GRCm39', 'GRCm38', 'GRCm37',
            'not applicable'
        ]
    }
}

req_pattern ={
    'obs': {
        'cell_enrichment': '^CL:(([\d]{7})(\+|\-)|([0]{7}))$',
        'gene_annotation_version': '^v(7[5-9]|[8-9][0-9]|10[0-9]|11[01])|GCF_000001405.(2[5-9]|3[0-9]|40)$'
    }
}

def validate_hca(adata):
    errors = []

    for c in req_present['uns']:
        if c not in adata.uns:
            errors.append(f'uns.{c} missing')

    for c in req_present['obs']:
        if c not in adata.obs.columns:
            errors.append(f'obs.{c} missing')

    for c,e in req_enum['obs'].items():
        if c not in adata.obs.columns:
            errors.append(f'obs.{c} missing')
        else:
            invalid = [v for v in adata.obs[c].unique() if v not in e]
            if invalid:
                errors.append(f'{",".join(invalid)} not valid for obs.{c}')

    for c,p in req_pattern['obs'].items():
        if c not in adata.obs.columns:
            errors.append(f'obs.{c} missing')
        else:
            invalid = [v for v in adata.obs[c].unique() if not re.match(p, v)]
            if invalid:
                errors.append(f'{",".join(invalid)} not valid for obs.{c}')

    return errors
