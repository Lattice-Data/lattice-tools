"""
Controlled term resolution, gatherer through flattener.

Controlled term references arrive in two shapes: some fields embed a dict
carrying '@id' and 'term_name', others hold only a path string like
'/controlled_terms/MONDO:0005148/'. The gatherer fetches the ControlledTerm
objects so both shapes resolve to the same dict, and the flattener keeps that
dict intact so split_controlled_term_columns() can emit _term_id / _term_name.
"""

import pandas as pd
import pytest

import DB2lattice
from constants import Configs, DEFAULT_DISEASE_TERM
from DB2Flattener import DB2Flattener
from DB2Gatherer import DB2Gatherer
from DB2_utils import split_controlled_term_columns


LUNG = {'@id': '/controlled_terms/UBERON:0002048/', 'term_name': 'lung'}
DIABETES = {'@id': '/controlled_terms/MONDO:0005148/', 'term_name': 'type 2 diabetes mellitus'}
CANCER = {'@id': '/controlled_terms/MONDO:0004992/', 'term_name': 'cancer'}
EUROPEAN = {'@id': '/controlled_terms/HANCESTRO:0005/', 'term_name': 'European'}
FLEX = {'@id': '/controlled_terms/EFO:0009922/', 'term_name': "10x 3' v3"}

RESOLVED_TERMS = {t['@id']: t for t in (LUNG, DIABETES, CANCER, EUROPEAN)}

UNRESOLVABLE = '/controlled_terms/MONDO:9999999/'


CONFIGS = Configs(
    FIELD_TYPES={
        'diseases': {'type': 'array', 'elements': 'string'},
        'sample_terms': {'type': 'array', 'elements': 'string'},
        'donors': {'type': 'array', 'elements': 'string'},
        'ethnicity': {'type': 'string'},
    },
    OBJECT_CONFIG={
        'droplet_based_libraries': {
            'api_type': 'DropletBasedLibrary',
            'fields': ['@id', 'feature_types', 'aliases'],
            'references': {},
        },
        'tissues': {
            'api_type': 'Tissue',
            'fields': ['@id', 'aliases', 'diseases', 'sample_terms', 'preservation_method'],
            'references': {
                'diseases': 'controlled_terms',
                'sample_terms': 'controlled_terms',
                'donors': ['human_donors'],
            },
        },
        'human_donors': {
            'api_type': 'HumanDonor',
            'fields': ['@id', 'ethnicity'],
            'references': {'ethnicity': 'controlled_terms'},
        },
        'controlled_terms': {
            'api_type': 'ControlledTerm',
            'fields': ['@id', 'term_id', 'term_name'],
            'references': {},
        },
    },
)


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = CONFIGS
    f.gatherer = None
    return f


# ---------------------------------------------------------------------------
# Gatherer: fetching ControlledTerm objects
# ---------------------------------------------------------------------------


@pytest.fixture
def recorded_requests(monkeypatch):
    """Stub DB2lattice.get_report and record every request it receives."""
    calls = []

    def fake_get_report(obj_type, filter_url, field_lst, connection):
        calls.append({'type': obj_type, 'filter': filter_url, 'fields': field_lst})
        if obj_type == 'HumanDonor':
            return [{'@id': '/human_donors/d1/', 'ethnicity': EUROPEAN['@id']}]
        if obj_type == 'ControlledTerm':
            return [t for t in RESOLVED_TERMS.values() if f"@id={t['@id']}" in filter_url]
        return []

    monkeypatch.setattr(DB2lattice, 'get_report', fake_get_report)
    return calls


def make_gatherer():
    return DB2Gatherer(connection=None, configs=CONFIGS)


def _sample_with_donor():
    return {
        '/tissues/t1/': {
            '@id': '/tissues/t1/',
            'diseases': [DIABETES['@id']],   # controlled term on the sample itself
            'donors': ['/human_donors/d1/'],  # donor carries a controlled term of its own
        }
    }


def test_controlled_terms_fetched_by_id_not_uuid(recorded_requests):
    """A controlled term reference is a term id, not a uuid, so @id is the filter."""
    g = make_gatherer()
    g.resolve_references_for_samples(_sample_with_donor())

    ct_call = next(c for c in recorded_requests if c['type'] == 'ControlledTerm')
    assert '@id=/controlled_terms/' in ct_call['filter']
    assert 'uuid=' not in ct_call['filter']


def test_other_types_still_fetched_by_uuid(recorded_requests):
    g = make_gatherer()
    g.resolve_references_for_samples(_sample_with_donor())

    donor_call = next(c for c in recorded_requests if c['type'] == 'HumanDonor')
    assert donor_call['filter'].startswith('&uuid=')


def test_controlled_terms_request_only_id_and_term_name(recorded_requests):
    """split_term_cell derives the term id from @id, so the rest is dead weight."""
    g = make_gatherer()
    g.resolve_references_for_samples(_sample_with_donor())

    ct_call = next(c for c in recorded_requests if c['type'] == 'ControlledTerm')
    assert ct_call['fields'] == ['@id', 'term_name']


def test_terms_reachable_only_via_a_fetched_object_are_included(recorded_requests):
    """
    The donor's ethnicity term cannot be seen until the donor itself has been
    fetched, so the controlled term fetch has to run after that second pass.
    """
    g = make_gatherer()
    g.resolve_references_for_samples(_sample_with_donor())

    resolved = g.resolved_objects['ControlledTerm']
    assert DIABETES['@id'] in resolved      # found on the sample, first pass
    assert EUROPEAN['@id'] in resolved      # found on the donor, second pass


def test_resolved_controlled_terms_keyed_by_id(recorded_requests):
    g = make_gatherer()
    g.resolve_references_for_samples(_sample_with_donor())

    assert g.resolved_objects['ControlledTerm'][DIABETES['@id']] == DIABETES


def test_no_controlled_term_refs_makes_no_request(recorded_requests):
    g = make_gatherer()
    assert g.fetch_controlled_terms(set()) == {}
    assert not [c for c in recorded_requests if c['type'] == 'ControlledTerm']


def test_oversized_response_raises(monkeypatch):
    """
    An ignored '@id' filter would return every controlled term in the instance.
    That is silently wrong rather than an error, so refuse it.
    """
    def dump_everything(obj_type, filter_url, field_lst, connection):
        return [
            {'@id': f'/controlled_terms/X:{i}/', 'term_name': f't{i}'}
            for i in range(5000)
        ]

    monkeypatch.setattr(DB2lattice, 'get_report', dump_everything)
    g = make_gatherer()

    with pytest.raises(RuntimeError, match="does not appear to have been applied"):
        g.fetch_controlled_terms({DIABETES['@id']})


def test_unresolved_terms_warn_but_do_not_raise(monkeypatch, capsys):
    monkeypatch.setattr(DB2lattice, 'get_report', lambda **kwargs: [])
    g = make_gatherer()

    assert g.fetch_controlled_terms({UNRESOLVABLE}) == {}
    assert 'did not resolve' in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Flattener: turning resolved terms into _term_id / _term_name columns
# ---------------------------------------------------------------------------


def _tissue(sample_id='/tissues/t1/', alias='lab:sample1', **overrides):
    tissue = {
        '@id': sample_id,
        '@type': ['Tissue'],
        'aliases': [alias],
        'preservation_method': 'frozen',
    }
    tissue.update(overrides)
    return tissue


def _lib():
    return {
        '@id': '/droplet_based_libraries/gex/',
        '@type': ['DropletBasedLibrary'],
        'feature_types': ['Gene Expression'],
        'aliases': ['lab:gex'],
    }


def _rmf(sample_id='/tissues/t1/', rmf_id='/raw_matrix_files/r1/', alias='lab:rmf1'):
    return {
        '@id': rmf_id,
        'aliases': [alias],
        'samples': [sample_id],
    }


def _main_df(tissues, extra_lib_data=None):
    """
    create_dataframe for one library with one raw matrix file per sample,
    with the controlled term columns split as flatten_matrix_file_set does.
    """
    rmfs = [
        _rmf(sample_id=t['@id'], rmf_id=f'/raw_matrix_files/r{i}/', alias=f'lab:rmf{i}')
        for i, t in enumerate(tissues)
    ]
    lib_data = {
        'library': _lib(),
        'samples': tissues,
        'raw_matrix_files': rmfs,
    }
    lib_data.update(extra_lib_data or {})

    main_df, _ = make_flattener().create_dataframe({
        'libraries': {'uuid-0': lib_data},
        'resolved_objects': {'ControlledTerm': RESOLVED_TERMS},
    })
    return split_controlled_term_columns(main_df)


def _main_row(tissue, extra_lib_data=None):
    """_main_df for a single sample, returning the one row."""
    return _main_df([tissue], extra_lib_data).iloc[0]


def test_array_field_resolves_to_parallel_id_and_name_columns():
    row = _main_row(_tissue(diseases=[DIABETES['@id'], CANCER['@id']]))

    assert list(row['tissues_diseases_term_id']) == ['MONDO:0005148', 'MONDO:0004992']
    assert list(row['tissues_diseases_term_name']) == [
        'type 2 diabetes mellitus',
        'cancer',
    ]


def test_single_element_array_collapses_to_scalar():
    row = _main_row(_tissue(sample_terms=[LUNG['@id']]))

    assert row['tissues_sample_terms_term_id'] == 'UBERON:0002048'
    assert row['tissues_sample_terms_term_name'] == 'lung'


def test_original_dict_column_is_dropped():
    row = _main_row(_tissue(sample_terms=[LUNG['@id']]))

    assert 'tissues_sample_terms' not in row.index


def test_non_controlled_term_field_untouched():
    row = _main_row(_tissue(sample_terms=[LUNG['@id']]))

    assert row['tissues_preservation_method'] == 'frozen'


def test_embedded_dict_passes_through_without_lookup():
    """Fields that already arrive embedded need no resolution."""
    f = make_flattener()
    assert f._resolve_controlled_term(FLEX, {}) is FLEX


def test_reference_on_a_related_object_resolves():
    """
    Path through _flatten_resolved_references: the term lives on the donor, not
    on the sample, and must not be flattened to a string by _join_unique.
    """
    donor = {'@id': '/human_donors/d1/', 'ethnicity': EUROPEAN['@id']}
    row = _main_row(
        _tissue(donors=['/human_donors/d1/']),
        extra_lib_data={'human_donors': [donor]},
    )

    assert row['human_donors_ethnicity_term_id'] == 'HANCESTRO:0005'
    assert row['human_donors_ethnicity_term_name'] == 'European'


# --- disease defaulting ---


@pytest.mark.parametrize(
    'tissue_kwargs',
    [
        pytest.param({}, id='diseases_key_absent'),
        pytest.param({'diseases': []}, id='diseases_empty_list'),
    ],
)
def test_absent_disease_defaults_to_normal(tissue_kwargs):
    row = _main_row(_tissue(**tissue_kwargs))

    assert row['tissues_diseases_term_id'] == 'PATO:0000461'
    assert row['tissues_diseases_term_name'] == 'normal'


def test_recorded_disease_is_not_overwritten_by_the_default():
    row = _main_row(_tissue(diseases=[DIABETES['@id']]))

    assert row['tissues_diseases_term_name'] == 'type 2 diabetes mellitus'


def test_unresolvable_disease_stays_empty_rather_than_normal():
    """
    A disease that failed to resolve is not the same as no disease. Reporting it
    as normal would hide the gatherer's resolution warning behind healthy-looking
    data.

    Uses two samples so the column is mixed: a column that is empty for every
    row is never detected as a controlled term column and so is never split.
    """
    df = _main_df([
        _tissue(sample_id='/tissues/t1/', alias='lab:s1', diseases=[DIABETES['@id']]),
        _tissue(sample_id='/tissues/t2/', alias='lab:s2', diseases=[UNRESOLVABLE]),
    ]).set_index('raw_file_samples')

    assert df.loc['s1', 'tissues_diseases_term_name'] == 'type 2 diabetes mellitus'
    assert pd.isna(df.loc['s2', 'tissues_diseases_term_name'])


def test_default_disease_term_is_not_shared_between_rows():
    """Each row gets its own copy, so a downstream mutation cannot leak."""
    f = make_flattener()
    first = dict(DEFAULT_DISEASE_TERM)
    second = dict(DEFAULT_DISEASE_TERM)

    assert first == second
    assert first is not second
    assert f._resolve_controlled_term_value(None, RESOLVED_TERMS) is None


# --- helpers ---


@pytest.mark.parametrize(
    'value,expected',
    [
        pytest.param(None, None, id='none'),
        pytest.param([], None, id='empty_list'),
        pytest.param([UNRESOLVABLE], None, id='all_unresolvable'),
    ],
)
def test_resolve_controlled_term_value_empty_cases(value, expected):
    f = make_flattener()
    assert f._resolve_controlled_term_value(value, RESOLVED_TERMS) is expected


def test_resolve_controlled_term_value_keeps_list_shape_for_arrays():
    f = make_flattener()
    result = f._resolve_controlled_term_value(
        [DIABETES['@id'], CANCER['@id']], RESOLVED_TERMS
    )
    assert result == [DIABETES, CANCER]


def test_resolve_controlled_term_value_scalar_for_non_arrays():
    f = make_flattener()
    assert f._resolve_controlled_term_value(LUNG['@id'], RESOLVED_TERMS) == LUNG


def test_dedupe_terms_drops_duplicates_and_sorts_by_id():
    f = make_flattener()
    result = f._dedupe_terms([CANCER, DIABETES, CANCER])

    assert [t['@id'] for t in result] == [CANCER['@id'], DIABETES['@id']]


def test_dedupe_terms_single_term_returns_scalar():
    f = make_flattener()
    assert f._dedupe_terms([LUNG, LUNG]) == LUNG


def test_dedupe_terms_empty_returns_none():
    f = make_flattener()
    assert f._dedupe_terms([]) is None
