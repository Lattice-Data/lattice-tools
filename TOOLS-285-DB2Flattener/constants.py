from generate_constants import load_and_return_constant_dicts


# Contains information about objects
# Their types, and what fields can be found in each

# URL length limit for chunking (includes base URL overhead)
MAX_URL_LENGTH = 3800
# Base URL overhead for chunking calculations (base URL + field params + safety margin)
BASE_URL_OVERHEAD = 700

# Global field type definitions
# What type to expect when getting each value
FIELD_TYPES, OBJECT_CONFIG = load_and_return_constant_dicts()

# Keys use _term_name suffix for columns produced by df_utils.split_controlled_term_columns
PROP_MAP_GEO = {
    'droplet_based_libraries_CRO_group_identifier': '*library name',
    'droplet_based_libraries_library_construction_technology_term_name': '*library strategy',
    'non_human_donors_taxa': '*organism',
    'tissues_sample_terms_term_name': '**tissue',
    'organoids_sample_terms_term_name': '**tissue',
    'cell_lines_sample_terms_term_name': '**cell_line',
    'tissues_enriched_cell_types_term_name': '**cell_type',
    'primary_cell_cultures_enriched_cell_types_term_name': '**cell_type',
    'raw_matrix_file_alias': 'raw_file',
    'droplet_based_libraries_library_cardinality': 'single or paired-end'
}

# API types to resolve when walking sample references.
# Controlled terms are always parsed from paths and never fetched here.
FETCHED_SAMPLE_REFERENCE_API_TYPES = frozenset({
    'HumanDonor',
    'NonHumanDonor',
})


