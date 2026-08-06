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


# Audit/provenance fields present on nearly every Lattice schema profile.
# Excluded from OBJECT_CONFIG so they don't get flattened into a column
# (and, in submitted_by's case, resolved as a reference) for every object type.
EXCLUDED_FIELDS = {'creation_timestamp', 'submitted_by'}

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

PROP_MAP_BIOHUB = {
    "raw_file_samples": "sample_name",
    "non_human_donors_cxg_donor_id":"donor_id",
    "human_donors_cxg_donor_id":"donor_id",
    "non_human_donors_taxa":"organism",
    "human_donors_taxa":"organism",
    "non_human_donors_sex":"sex",
    "human_donors_sex":"sex",
    "human_donors_ethnicity":"self_reported_ethnicity",
    "droplet_based_libraries_library_construction_technology_term_name":"assay",
    "tissues_upper_bound_age": "tissues_upper_bound_age",
    "tissues_lower_bound_age": "tissues_lower_bound_age",
    "tissues_age_units": "tissues_age_units",
    "tissues_diseases": "disease",
    "tissues_sample_terms_term_name":"tissue",
    "tissues_developmental_stages_term_name":"development_stage",
    "tissues_multiplexing_barcodes": "sample_probe_barcode",
    "tissues_@type":"tissue_type",
    "tissues_selection_markers":"suspension_enrichment_factors",
    "tissues_suspension_type":"suspension_type",
    "tissues_preservation_method":"preservation_method",
    "experimental_conditions_condition": "experimental_condition",
    "experimental_conditions_text_value": "experimental_perturbation",
    "experimental_conditions_upper_bound_duration": "experimental_conditions_upper_bound_duration",
    "experimental_conditions_lower_bound_duration": "experimental_conditions_lower_bound_duration",
    "experimental_conditions_duration_units": "experimental_conditions_duration_units",
    "genetic_modifications_strategy":"genetic_perturbation_strategy"
}

TISSUE_TYPE_MAP = {
    "CellLine": "cell line",
    "Organoid": "organoid",
    "PrimaryCellCulture": "primary cell culture",
    "Tissue": "tissue"
}

GENETIC_PERTURBATION_MAP = {
    "activation screen":"CRISPR activation screen",
    "interference screen":"CRISPR interference screen",
    "knockout mutation":"CRISPR knockout mutant",
    "knockout screen":"CRISPR knockout screen",
}

REFORMAT_LIST = [
    "sample_probe_barcode",
    "suspension_enrichment_factors",
]
