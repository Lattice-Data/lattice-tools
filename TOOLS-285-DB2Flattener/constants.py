from dataclasses import dataclass
from typing import Any, TypeAlias


# Contains information about objects
# Their types, and what fields can be found in each

# typing for various objects used with JSON profile parsing
Hierarchy: TypeAlias = dict[str, dict[str, dict]]
JSONProfile: TypeAlias = dict[str, Any]
FieldTypes: TypeAlias = dict[str, dict[str, str]]
ObjectConfig: TypeAlias = dict[str, dict[str, Any]]


@dataclass
class Configs:
    """
    Use as a container to hold the parsed configs from profile schemas
    Structure:
    FIELD_TYPES: {
        {field}: {
            "type": {datatype value},
            "elements {optional}": {datatype of collection items},
        }
    }

    OBJECT_CONFIG: {
        {object url_prefix}: {
            "api_type": {object API Name},
            "fields": list[fields],
            "references": {
                "{field}": {object url_prefix}
            }
        }
    }
        
    """
    FIELD_TYPES: FieldTypes
    OBJECT_CONFIG: ObjectConfig


# URL length limit for chunking (includes base URL overhead)
MAX_URL_LENGTH = 3800
# Base URL overhead for chunking calculations (base URL + field params + safety margin)
BASE_URL_OVERHEAD = 700

# Keys use _term_name suffix for columns produced by DB2_utils.split_controlled_term_columns
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