# Contains information about objects
# Their types, and what fields can be found in each

# URL length limit for chunking (includes base URL overhead)
MAX_URL_LENGTH = 3800
# Base URL overhead for chunking calculations (base URL + field params + safety margin)
BASE_URL_OVERHEAD = 700

# Global field type definitions
# What type to expect when getting each value
FIELD_TYPES = {
    'uuid': {'type': 'string'},
    '@id': {'type': 'string'},
    'date_obtained': {'type': 'string'},
    'age_units': {'type': 'string'},
    'lower_bound_age': {'type': 'number'},
    'upper_bound_age': {'type': 'number'},
    'genetic_modification': {'type': 'string'},
    'sex': {'type': 'string'},
    'ethnicity': {'type': 'string'},
    'term_name': {'type': 'string'}, # Don't think this is implemented yet in the API right?
    'strategy': {'type': 'string'},
    'ontological_term': {'type': 'string'},
    'controlled_terms': {'type': 'string'},
    'host': {'type': 'string'},
    'host_tissue': {'type': 'string'},
    'taxa': {'type': 'string'},
    'library_cardinality': {'type': 'string'},
    'chemistry_version': {'type': 'string'},
    'feature_types': {'type': 'string'},
    'multiplexing_method': {'type': 'string'},
    'kit_version': {'type': 'string'},
    's3_uri': {'type': 'string'},
    'crc64nvme_base64': {'type': 'string'},
    'suspension_type': {'type': 'string'},
    'library_construction_technology' : {'type': 'string'},
    'condition' : {'type':'string'},
    'file_format' : {'type':'string'},
    'read_count' : {'type': 'number'},
    'CRO_order' : {'type': 'string'},
    'library' : {'type': 'string'},
    'read1' : {'type': 'string'},
    'read2' : {'type': 'string'},
    'read3' : {'type': 'string'},
    'index1' : {'type': 'string'},
    'index2' : {'type': 'string'},
    'sequencing_platform' : {'type': 'string'},
    'trimmed_cram' : {'type': 'string'},
    'untrimmed_cram' : {'type': 'string'},
    'run_cardinality' : {'type': 'string'},
    'genome_annotation' : {'type': 'string'},
    'genome_assembly' : {'type': 'string'},
    'software' : {'type': 'string'},
    'software_version' : {'type': 'string'},
    'CRO_group_identifier' : {'type': 'string'},
    
    # Array fields with element types
    'sample_terms': {'type': 'array', 'elements': 'string'},
    'donors': {'type': 'array', 'elements': 'string'},
    'aliases': {'type': 'array', 'elements': 'string'},
    'treatments': {'type': 'array', 'elements': 'string'},
    'diseases': {'type': 'array', 'elements': 'string'},
    'enriched_cell_types': {'type': 'array', 'elements': 'string'},
    'depleted_cell_types': {'type': 'array', 'elements': 'string'},
    'intended_cell_types': {'type': 'array', 'elements': 'string'},
    'samples': {'type': 'array', 'elements': 'string'},
    'derived_from': {'type': 'array', 'elements': 'string'},
    'feature_keys': {'type': 'array', 'elements': 'string'},
    'experimental_conditions' : {'type':'array', 'elements': 'string'},
    'source_sequence_file_sets' : {'type':'array', 'elements': 'string'},
    'processed_matrix_files' : {'type':'array', 'elements': 'string'},
    'raw_matrix_files' : {'type':'array', 'elements': 'string'},
    'sequence_file_sets' : {'type':'array', 'elements': 'string'},
    
    # Special handling cause object
    'feature_counts': {'type': 'array', 'elements': 'object'}
}

OBJECT_CONFIG = {
    # Biosamples
    'tissues': {
        'api_type': 'Tissue',
        'fields': [
            'uuid', '@id', 'aliases','date_obtained', 'sample_terms', 'donors', 
            'age_units', 'lower_bound_age', 'upper_bound_age', 
            'genetic_modification', 'treatments', 'suspension_type', 
            'experimental_conditions', 'enriched_cell_types',
            'depleted_cell_types', 'diseases', 'selection_markers'
            
        ],
        'references': {
            'sample_terms': 'controlled_terms',
            'donors': ['human_donors', 'non_human_donors'],
            'genetic_modification': 'genetic_modifications',
            'treatments': 'treatments',
            'experimental_conditions' : 'experimental_conditions',
            'enriched_cell_types' : 'controlled_terms',
            'depleted_cell_types' : 'controlled_terms',
            'diseases' : 'controlled_terms'
            
        }
    },
    
    'cell_lines': {
        'api_type': 'CellLine',
        'fields': [
            'uuid', '@id', 'aliases', 'date_obtained', 'sample_terms', 'donors', 
            'genetic_modification', 'treatments', 'suspension_type', 
            'experimental_conditions', 'enriched_cell_types',
            'depleted_cell_types', 'diseases', 'host', 'host_tissue', 'selection_markers'
            
        ],
        'references': {
            'sample_terms': 'controlled_terms',
            'donors': ['human_donors', 'non_human_donors'],
            'genetic_modification': 'genetic_modifications',
            'treatments': 'treatments',
            'experimental_conditions' : 'experimental_conditions',
            'enriched_cell_types' : 'controlled_terms',
            'depleted_cell_types' : 'controlled_terms',
            'diseases' : 'controlled_terms',
            'host' : ['human_donors', 'non_human_donors'],
            'host_tissue' : 'controlled_terms'
            
        }
    },
    
    'organoids': {
        'api_type': 'Organoid',
        'fields': [
            'uuid', '@id', 'aliases', 'date_obtained', 'sample_terms', 'donors', 
            'genetic_modification', 'treatments', 'suspension_type', 
            'experimental_conditions', 'enriched_cell_types',
            'depleted_cell_types', 'diseases', 'intended_cell_types', 'selection_markers'
            
        ],
        'references': {
            'sample_terms': 'controlled_terms',
            'donors': ['human_donors', 'non_human_donors'],
            'genetic_modification': 'genetic_modifications',
            'treatments': 'treatments',
            'experimental_conditions' : 'experimental_conditions',
            'enriched_cell_types' : 'controlled_terms',
            'depleted_cell_types' : 'controlled_terms',
            'diseases' : 'controlled_terms',
            'intended_cell_types' : 'controlled_terms',
            
        }
    },
    
    'primary_cell_cultures': {
        'api_type': 'PrimaryCellCulture',
        'fields': [
            'uuid', '@id', 'aliases', 'date_obtained', 'sample_terms', 'donors', 
            'age_units', 'lower_bound_age', 'upper_bound_age', 
            'genetic_modification', 'treatments', 'suspension_type', 
            'experimental_conditions', 'enriched_cell_types',
            'depleted_cell_types', 'diseases', 'selection_markers'
            
        ],
        'references': {
            'sample_terms': 'controlled_terms',
            'donors': ['human_donors', 'non_human_donors'],
            'genetic_modification': 'genetic_modifications',
            'treatments': 'treatments',
            'experimental_conditions' : 'experimental_conditions',
            'enriched_cell_types' : 'controlled_terms',
            'depleted_cell_types' : 'controlled_terms',
            'diseases' : 'controlled_terms'
            
        }
    },
    
    # Donors
    'human_donors': {
        'api_type': 'HumanDonor',
        'fields': ['uuid', '@id', 'sex', 'ethnicity', 'taxa'],
        'references': {
            'ethnicity': 'controlled_terms'
        }
    },
    'non_human_donors': {
        'api_type': 'NonHumanDonor',
        'fields': ['uuid', '@id', 'sex', 'taxa'],
        'references': {
        }
    },
    
    # Libraries
    'droplet_based_libraries': {
        'api_type': 'DropletBasedLibrary',
        'fields': ['samples', 'uuid', 'aliases', 'library_construction_technology', 'library_cardinality',
                  'chemistry_version', 'feature_types', 'multiplexing_method', 'CRO_group_identifier'],
        'references': {
            'samples' : ['primary_cell_cultures',
                         'organoids',
                         'cell_lines',
                         'tissues'],
            'library_construction_technology' : 'controlled_terms'
        }
    },
    'plate_based_libraries': {
        'api_type': 'PlateBasedLibrary',
        'fields': ['samples', 'uuid', 'aliases', 'library_construction_technology',
                  'kit_version', 'multiplexing_method', 'CRO_group_identifier'],
        'references': {
            'samples' : ['primary_cell_cultures',
                         'organoids',
                         'cell_lines',
                         'tissues'],
            'library_construction_technology' : 'controlled_terms'
        }
    },
    
    # Files
    'raw_matrix_files': {
        'api_type': 'RawMatrixFile',
        'fields': ['uuid', 'aliases', 's3_uri', 'crc64nvme_base64', 'feature_keys', 'derived_from', 'feature_counts',
                  'file_format', 'samples'],
        'references': {
            'derived_from' : 'sequence_files',
            'samples' : ['primary_cell_cultures',
                         'organoids',
                         'cell_lines',
                         'tissues'],
        }
    },
    'processed_matrix_files': {
        'api_type': 'ProcessedMatrixFile',
        'fields': ['uuid', 'aliases', 's3_uri', 'crc64nvme_base64', 'feature_keys', 'derived_from', 'feature_counts',
                  'file_format'],
        'references': {
            'derived_from' : 'sequence_files' # Or should this be raw matrix files? or both?
        }
    },
    'sequence_files': {
        'api_type': 'SequenceFile',
        'fields': ['uuid', 'aliases', 's3_uri', 'crc64nvme_base64', 'feature_keys', 'feature_counts',
                  'file_format', 'sequence_file_sets'],
        'references': {
            'sequence_file_sets':'sequence_file_sets'
        }
    },

    # FileSets
    'sequence_file_sets': {
        'api_type': 'SequenceFileSet',
        'fields': ['uuid', 'aliases', 'CRO_order', 'index1', 'index2', 'library', 
                   'read1', 'read2', 'read3', 'run_cardinality', 'sequencing_platform',
                  'trimmed_cram', 'untrimmed_cram'],
        'references': {
            'index1':'sequence_files',
            'index2':'sequence_files',
            'library':['droplet_based_libraries','plate_based_libraries'],
            'read1':'sequence_files',
            'read2':'sequence_files',
            'read3':'sequence_files',
            'trimmed_cram' : 'sequence_files',
            'untrimmed_cram' : 'sequence_files'
        }
    },
    'matrix_file_sets': {
        'api_type': 'MatrixFileSet',
        'fields': ['uuid', 'aliases', 'genome_annotation', 'genome_assembly',
                  'processed_matrix_files', 'raw_matrix_files', 'software', 'software_version',
                  'source_sequence_file_sets'],
        'references': {
            'processed_matrix_files' : 'processed_matrix_files',
            'raw_matrix_files' : 'raw_matrix_files',
            'source_sequence_file_sets' : 'sequence_file_sets',
        }
    },

    # References
    'controlled_terms': { # Special handling, need to get term_id directly from string
        'api_type': 'ControlledTerm',
        'fields': ['@id', 'term_name'],
        'references': {}  # No nested references
    },
    
    'genetic_modifications': {
        'api_type': 'GeneticModification',
        'fields': ['@id', 'strategy'],
        'references': {}
    },
    
    'treatments': {
        'api_type': 'Treatment', 
        'fields': ['@id', 'ontological_term'],
        'references': {}
    },
    'experimental_conditions': {
        'api_type': 'ExperimentalCondition', 
        'fields': ['@id', 'condition', 'text_value'],
        'references': {}
    }
    
}