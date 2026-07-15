import pandas as pd
import pytest

from constants import Configs
from DB2Flattener import DB2Flattener


MIN_CONFIGS = Configs(
    FIELD_TYPES={},
    OBJECT_CONFIG={
        'droplet_based_libraries': {
            'api_type': 'DropletBasedLibrary',
            'fields': ['@id', 'feature_types', 'CRO_group_identifier', 'aliases'],
            'references': {},
        },
        'tissues': {
            'api_type': 'Tissue',
            'fields': ['@id', 'aliases', 'author_metadata'],
            'references': {},
        },
    },
)


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = MIN_CONFIGS
    f.gatherer = None
    return f


def _tissue(sample_id='/tissues/s1/', alias='lab:sample1', author_metadata=None):
    return {
        '@id': sample_id,
        '@type': ['Tissue'],
        'aliases': [alias],
        'author_metadata': author_metadata,
    }


def _lib(lib_id, feature_types, cro='RSJS_fast_1', alias=None):
    return {
        '@id': lib_id,
        '@type': ['DropletBasedLibrary'],
        'feature_types': feature_types,
        'CRO_group_identifier': cro,
        'aliases': [alias or f'lab:{feature_types[0].replace(" ", "_")}'],
    }


def _rmf(rmf_id='/raw_matrix_files/r1/', alias='lab:rmf1', sample_id='/tissues/s1/'):
    return {
        '@id': rmf_id,
        'aliases': [alias],
        'samples': [sample_id],
    }


def _complete_data(libs_with_rmfs):
    """libs_with_rmfs: list of (library, [raw_matrix_files], [samples])"""
    libraries = {}
    for i, (lib, rmfs, samples) in enumerate(libs_with_rmfs):
        libraries[f'uuid-{i}'] = {
            'library': lib,
            'raw_matrix_files': rmfs,
            'samples': samples,
        }
    return {'libraries': libraries, 'resolved_objects': {'ControlledTerm': {}}}


# --- _row_is_gex (new method) ---

@pytest.mark.parametrize(
    'row,expected',
    [
        ({'droplet_based_libraries_feature_types': ['Gene Expression']}, True),
        ({'droplet_based_libraries_feature_types': ['CRISPR Guide Capture']}, False),
        ({'droplet_based_libraries_feature_types': 'Gene Expression'}, True),
        ({'plate_based_libraries_feature_types': ['Gene Expression']}, True),
        # missing FT: plate assumed GEX
        ({'plate_based_libraries_@id': '/plate_based_libraries/p1/'}, True),
        # missing FT: droplet assumed non-GEX
        ({'droplet_based_libraries_@id': '/droplet_based_libraries/d1/'}, False),
    ],
)
def test_row_is_gex(row, expected):
    f = make_flattener()
    assert f._row_is_gex(pd.Series(row)) is expected


# --- create_dataframe: one row per (RMF, library) ---

def test_create_dataframe_one_row_per_rmf_library():
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib('/droplet_based_libraries/gex/', ['Gene Expression'], alias='lab:gex')
    cri = _lib('/droplet_based_libraries/cri/', ['CRISPR Guide Capture'], alias='lab:cri')

    main_df, sample_df = f.create_dataframe(_complete_data([
        (gex, [rmf], [sample]),
        (cri, [rmf], [sample]),
    ]))

    assert len(main_df) == 2
    assert set(main_df['droplet_based_libraries_@id']) == {
        '/droplet_based_libraries/gex/',
        '/droplet_based_libraries/cri/',
    }
    assert set(main_df['raw_matrix_file_alias']) == {'rmf1'}
    ftypes = set()
    for ft in main_df['droplet_based_libraries_feature_types']:
        ftypes.update(ft if isinstance(ft, list) else [ft])
    assert 'Gene Expression' in ftypes
    assert 'CRISPR Guide Capture' in ftypes


def test_create_dataframe_dedupes_same_library_twice_on_rmf():
    """Same lib attached via two library buckets should not double MAIN rows."""
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib('/droplet_based_libraries/gex/', ['Gene Expression'])

    main_df, _ = f.create_dataframe(_complete_data([
        (gex, [rmf], [sample]),
        (gex, [rmf], [sample]),
    ]))

    assert len(main_df) == 1
    assert main_df.iloc[0]['droplet_based_libraries_@id'] == '/droplet_based_libraries/gex/'


def test_sample_df_not_inflated_when_rmf_has_two_libraries():
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib('/droplet_based_libraries/gex/', ['Gene Expression'], alias='lab:gex')
    cri = _lib('/droplet_based_libraries/cri/', ['CRISPR Guide Capture'], alias='lab:cri')

    main_df, sample_df = f.create_dataframe(_complete_data([
        (gex, [rmf], [sample]),
        (cri, [rmf], [sample]),
    ]))

    assert len(main_df) == 2
    assert sample_df is not None
    assert sample_df.reset_index()['raw_matrix_file_alias'].nunique() == 1
    assert len(sample_df) == 1


# --- author_metadata explode ---

def test_author_metadata_dict_exploded_to_columns():
    f = make_flattener()
    sample = _tissue(author_metadata={'mouse litter batch': 'A1', 'diet': 'fasted'})
    rmf = _rmf()
    gex = _lib('/droplet_based_libraries/gex/', ['Gene Expression'])

    main_df, _ = f.create_dataframe(_complete_data([
        (gex, [rmf], [sample]),
    ]))

    assert 'tissues_mouse_litter_batch' in main_df.columns
    assert 'tissues_diet' in main_df.columns
    assert main_df.iloc[0]['tissues_mouse_litter_batch'] == 'A1'
    assert main_df.iloc[0]['tissues_diet'] == 'fasted'
    assert 'tissues_author_metadata' not in main_df.columns


def test_author_metadata_non_dict_kept_as_normal_field():
    f = make_flattener()
    sample = _tissue(author_metadata=None)
    rmf = _rmf()
    gex = _lib('/droplet_based_libraries/gex/', ['Gene Expression'])

    main_df, _ = f.create_dataframe(_complete_data([
        (gex, [rmf], [sample]),
    ]))

    assert 'tissues_author_metadata' in main_df.columns
    assert 'tissues_mouse_litter_batch' not in main_df.columns
