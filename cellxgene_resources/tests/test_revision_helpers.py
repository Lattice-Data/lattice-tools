"""
Unit tests for the revision helpers in cellxgene_mods (experimental_condition revisions).
Run from cellxgene_resources/ in an env with cellxgene-schema installed:  pytest tests/
"""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import cellxgene_mods as cxgm

EC = cxgm.EC_FIELD


def make_adata():
    rng = np.random.default_rng(0)
    X = sparse.csr_matrix(rng.poisson(1.0, size=(20, 5)).astype(np.float32))
    obs = pd.DataFrame(
        {
            'sample_strain': pd.Categorical(['ctrl'] * 10 + ['treated'] * 8 + ['other'] * 2),
            'donor_id': pd.Categorical(['d1'] * 20),
            'cell_type': pd.Categorical(['T cell'] * 20),  # portal-added label
            'observation_joinid': [f'j{i}' for i in range(20)],  # portal-added
            'is_primary_data': True,
        },
        index=[f'cell{i}' for i in range(20)],
    )
    var = pd.DataFrame(
        {'feature_name': [f'g{i}' for i in range(5)], 'feature_is_filtered': False},
        index=[f'ENSG{i:011d}' for i in range(5)],
    )
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.raw = a.copy()
    a.obsm['X_umap'] = rng.normal(size=(20, 2))
    a.uns['title'] = 'synthetic'
    a.uns['citation'] = 'portal writes this'
    a.uns['organism'] = 'Mus musculus'
    return a


# ---------------------------------------------------------------- add column

def test_add_experimental_condition_basic():
    a = make_adata()
    ct = cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})
    assert a.obs[EC].dtype.name == 'category'
    assert (a.obs.loc[a.obs['sample_strain'] == 'treated', EC] == 'CHEBI:53076').all()
    assert (a.obs.loc[a.obs['sample_strain'] != 'treated', EC] == 'na').all()
    assert ct.loc['treated', 'CHEBI:53076'] == 8
    assert ct.loc['ctrl', 'na'] == 10


def test_add_experimental_condition_multi_term_sorted_and_deduplicated():
    a = make_adata()
    cxgm.add_experimental_condition(
        a, 'sample_strain',
        {'treated': ['uniprot:P05112', 'CHEBI:16412', 'EFO:0001702', 'CHEBI:16412']},
    )
    treated = a.obs.loc[a.obs['sample_strain'] == 'treated', EC].unique()
    assert list(treated) == ['CHEBI:16412 || EFO:0001702 || uniprot:P05112']


def test_add_experimental_condition_two_source_values():
    a = make_adata()
    cxgm.add_experimental_condition(
        a, 'sample_strain', {'treated': ['CHEBI:53076'], 'other': ['EFO:0001702']}
    )
    assert (a.obs.loc[a.obs['sample_strain'] == 'other', EC] == 'EFO:0001702').all()
    assert (a.obs[EC] == 'na').sum() == 10


def test_add_experimental_condition_missing_source_column():
    with pytest.raises(KeyError):
        cxgm.add_experimental_condition(make_adata(), 'nope', {'treated': ['CHEBI:53076']})


def test_add_experimental_condition_unmapped_value():
    with pytest.raises(ValueError, match='not found'):
        cxgm.add_experimental_condition(make_adata(), 'sample_strain', {'typo': ['CHEBI:53076']})


def test_add_experimental_condition_all_na_refused():
    a = make_adata()
    a.obs['sample_strain'] = a.obs['sample_strain'].cat.add_categories(['unused'])
    # mapping a value that no cell carries is caught earlier; force all-na another way
    a.obs.loc[:, 'sample_strain'] = 'ctrl'
    a.obs['sample_strain'] = a.obs['sample_strain'].cat.add_categories(['treated']) \
        if 'treated' not in a.obs['sample_strain'].cat.categories else a.obs['sample_strain']
    with pytest.raises(ValueError):
        cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})


def test_add_experimental_condition_refuses_to_overwrite():
    a = make_adata()
    cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})
    with pytest.raises(ValueError, match='already has'):
        cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})
    cxgm.add_experimental_condition(a, 'sample_strain', {'other': ['CHEBI:53076']}, overwrite=True)
    assert (a.obs.loc[a.obs['sample_strain'] == 'treated', EC] == 'na').all()


def test_format_experimental_condition():
    assert cxgm.format_experimental_condition(['b', 'a', 'b']) == 'a || b'
    assert cxgm.format_experimental_condition(['CHEBI:53076']) == 'CHEBI:53076'


# ---------------------------------------------------------------- term checks

@pytest.mark.parametrize('term', [
    'CHEBI:53076',   # oxazolone
    'CHEBI:16412',   # lipopolysaccharide, from the schema's own example
    'EFO:0001702',   # temperature
    'EFO:0002755',   # diet
    'EFO:0002756',   # fasting
    'EFO:0002757',   # high fat diet, a diet descendant
    'uniprot:P05112',
    'anti-uniprot:Q99467',
])
def test_usable_terms(term):
    assert cxgm.check_experimental_condition_term(term) == []


@pytest.mark.parametrize('term, fragment', [
    ('CHEBI:24431', 'forbidden list'),           # chemical entity itself
    ('CHEBI:25212', 'forbidden list'),           # metabolite
    ('CHEBI:23888', 'descendant of forbidden CHEBI:50906'),  # drug, a role
    ('CHEBI:999999999', 'pinned ChEBI release'),
    ('EFO:0009899', 'limited to'),               # an assay term
    ('EFO:9999999', 'pinned EFO release'),
    ('uniprot:notanaccession', 'UniProt accession'),
    ('MONDO:0002052', 'prefix'),
    ('oxazolone', 'prefix'),
])
def test_unusable_terms(term, fragment):
    problems = cxgm.check_experimental_condition_term(term)
    assert problems, term
    assert any(fragment in p for p in problems), problems


def test_check_terms_returns_labels_and_raises_on_bad():
    labels = cxgm.check_experimental_condition_terms(['CHEBI:53076', 'uniprot:P05112'])
    assert labels['CHEBI:53076'].startswith('4-(ethoxymethylene)')
    assert labels['uniprot:P05112'] == 'uniprot:P05112'
    with pytest.raises(ValueError, match='unusable'):
        cxgm.check_experimental_condition_terms(['CHEBI:53076', 'CHEBI:24431'])


# ---------------------------------------------------------------- spec

def write_spec(tmp_path, spec):
    p = tmp_path / 'spec.json'
    p.write_text(json.dumps(spec))
    return p


def good_spec():
    return {
        'collection_id': 'c',
        'ticket': 'CXG-950',
        'datasets': {'d1': {'source_column': 'sample_strain', 'terms': {'treated': ['CHEBI:53076']}}},
    }


def test_load_revision_spec_ok(tmp_path):
    spec = cxgm.load_revision_spec(write_spec(tmp_path, good_spec()))
    assert cxgm.spec_terms(spec) == ['CHEBI:53076']


@pytest.mark.parametrize('mutate', [
    lambda s: s.pop('collection_id'),
    lambda s: s.__setitem__('datasets', {}),
    lambda s: s['datasets']['d1'].pop('source_column'),
    lambda s: s['datasets']['d1'].__setitem__('terms', {}),
    lambda s: s['datasets']['d1'].__setitem__('terms', {'treated': 'CHEBI:53076'}),
    lambda s: s['datasets']['d1'].__setitem__('terms', {'treated': []}),
])
def test_load_revision_spec_rejects_bad(tmp_path, mutate):
    spec = good_spec()
    mutate(spec)
    with pytest.raises(ValueError):
        cxgm.load_revision_spec(write_spec(tmp_path, spec))


# ---------------------------------------------------------------- fingerprints

def test_fingerprint_unchanged_passes():
    a = make_adata()
    cxgm.assert_only_changes(cxgm.fingerprint_adata(a), cxgm.fingerprint_adata(a))


def test_fingerprint_full_revision_passes():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})
    a.obs.drop(columns=['cell_type', 'observation_joinid'], inplace=True)
    a.var.drop(columns=['feature_name'], inplace=True)
    a.raw.var.drop(columns=['feature_name'], inplace=True)
    del a.uns['citation']
    del a.uns['organism']
    cxgm.assert_only_changes(
        before, cxgm.fingerprint_adata(a),
        added_obs=[EC], removed_obs=cxgm.OBS_PORTAL_ALL,
        removed_var=cxgm.VAR_PORTAL_REQUIRED, removed_uns=cxgm.UNS_PORTAL_REQUIRED,
    )


def test_fingerprint_catches_matrix_change():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.X[0, 0] = a.X[0, 0] + 1
    with pytest.raises(AssertionError, match='X changed'):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))


def test_fingerprint_catches_raw_change():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    raw = a.raw.to_adata()
    raw.X[1, 1] = raw.X[1, 1] + 1
    a.raw = raw
    with pytest.raises(AssertionError, match='raw.X changed'):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))


def test_fingerprint_catches_unexpected_obs_column_and_value_change():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.obs['surprise'] = 1
    with pytest.raises(AssertionError, match="obs\\['surprise'\\] added"):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.obs['donor_id'] = a.obs['donor_id'].cat.add_categories(['d2'])
    a.obs.iloc[0, a.obs.columns.get_loc('donor_id')] = 'd2'
    with pytest.raises(AssertionError, match="obs\\['donor_id'\\] changed"):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))


def test_fingerprint_catches_embedding_and_index_change():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.obsm['X_umap'][0, 0] += 1
    with pytest.raises(AssertionError, match="obsm\\['X_umap'\\] changed"):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.obs.index = [f'x{i}' for i in range(20)]
    with pytest.raises(AssertionError, match='__index__'):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))


def test_fingerprint_catches_disallowed_removal_and_missing_addition():
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    a.obs.drop(columns=['donor_id'], inplace=True)
    with pytest.raises(AssertionError, match="obs\\['donor_id'\\] removed"):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a))
    a = make_adata()
    before = cxgm.fingerprint_adata(a)
    with pytest.raises(AssertionError, match='not added'):
        cxgm.assert_only_changes(before, cxgm.fingerprint_adata(a), added_obs=[EC])


def test_fingerprint_survives_write_read_round_trip(tmp_path):
    a = make_adata()
    cxgm.add_experimental_condition(a, 'sample_strain', {'treated': ['CHEBI:53076']})
    before = cxgm.fingerprint_adata(a)
    p = tmp_path / 'a.h5ad'
    a.write(p, compression='gzip')
    cxgm.assert_only_changes(before, cxgm.fingerprint_adata(ad.read_h5ad(p)), ignore_uns=True)


# ---------------------------------------------------------------- misc

def test_build_upload_manifest():
    assert cxgm.build_upload_manifest({'anndata': 's3://old'}, 's3://new') == {'anndata': 's3://new'}
    m = cxgm.build_upload_manifest({'anndata': 's3://old', 'atac_fragment': 's3://frag'}, 's3://new')
    assert m == {'anndata': 's3://new', 'atac_fragment': 's3://frag'}


def test_revised_path():
    assert cxgm.revised_path('/x/abc.h5ad').name == 'abc_revised.h5ad'


def test_download_skips_when_size_matches(tmp_path, monkeypatch):
    class FakeResponse:
        def __init__(self, payload):
            self._payload = payload
        def raise_for_status(self):
            pass
        def json(self):
            return self._payload
    existing = tmp_path / 'd1.h5ad'
    existing.write_bytes(b'x' * 10)
    monkeypatch.setattr(
        cxgm.requests, 'get',
        lambda url, **kw: FakeResponse({'assets': [{'filetype': 'H5AD', 'url': 'u', 'filesize': 10}]}),
    )
    assert cxgm.download_dataset_h5ad('c', 'd1', tmp_path) == existing
