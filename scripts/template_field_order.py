priorityFields = {
	'publication': {
		'order': [],
		'preferred': []
	},
	'dataset': {
		'order': [],
		'preferred': []
	},
	'human_postnatal_donor': {
		'order': [
			'donor_id',
			'sex',
			'development_ontology.term_name',
			'development_ontology.term_id',
			'age',
			'age_units',
			'description',
			'reported_ethnicity',
			'ethnicity',
			'living_at_sample_collection',
			'causes_of_death',
			'death_type',
			'body_mass_index',
			'diseases',
			'treatments',
			'disease_state',
			'ancestry_method',
			'ancestry',
			'ancestry.fraction',
			'ancestry.ancestry_group.term_id',
			'ancestry.ancestry_group.term_name',
			'genotype',
			'weight',
			'weight_unit'
		],
		'preferred': [
			'age_units',
			'description',
			'reported_ethnicity',
			'ethnicity',
			'living_at_sample_collection',
			'body_mass_index',
			'diseases',
		]
	},
	'human_prenatal_donor': {
		'order': [
			'donor_id',
			'development_ontology.term_name',
			'development_ontology.term_id',
			'sex',
			'conceptional_age',
			'conceptional_age_units',
			'description'
		],
		'preferred': [
			'description'
		]
	},
	'tissue': {
		'order': [
			'derivation_process',
			'derived_from',
			'biosample_ontology.term_name',
			'biosample_ontology.term_id',
			'description',
			'spatial_information',
			'source',
			'preservation_method',
			'preservation_time',
			'preservation_time_units',
			'death_to_preservation_interval',
			'death_to_preservation_interval_units',
			'collection_to_preservation_interval',
			'collection_to_preservation_interval_units',
			'date_obtained',
			'diseases',
			'dbxrefs'
		],
		'preferred': [
			'description',
			'spatial_information',
			'source',
			'preservation_method',
			'preservation_time',
			'preservation_time_units',
			'death_to_preservation_interval',
			'death_to_preservation_interval_units',
			'collection_to_preservation_interval',
			'collection_to_preservation_interval_units',
			'date_obtained',
			'diseases',
			'dbxrefs'
		]
	},
	'cell_culture': {
		'order': [
			'derivation_process',
			'derived_from',
			'biosample_ontology.term_name',
			'biosample_ontology.term_id',
			'description',
			'growth_medium',
			'starting_quantity',
			'starting_quantity_units',
			'treatments',
			'genetic_modifications',
			'source',
			'product_id',
			'dbxrefs'
		],
		'preferred': [
			'description',
			'growth_medium',
			'starting_quantity',
			'starting_quantity_units',
			'treatments',
			'genetic_modifications',
			'source',
			'product_id',
			'dbxrefs'
		]
	},
	'organoid': {
		'order': [
			'derivation_process',
			'derived_from',
			'biosample_ontology.term_name',
			'biosample_ontology.term_id',
			'description',
			'growth_medium',
			'embedded_in_matrigel',
			'treatments',
			'starting_quantity',
			'starting_quantity_units',
			'genetic_modifications',
			'preservation_method',
			'preservation_time',
			'preservation_time_units',
			'dbxrefs'
		],
		'preferred': [
			'description',
			'growth_medium',
			'embedded_in_matrigel',
			'treatments',
			'starting_quantity',
			'starting_quantity_units',
			'genetic_modifications',
			'preservation_method',
			'preservation_time',
			'preservation_time_units',
			'dbxrefs'
		]
	},
	'suspension': {
		'order': [
			'suspension_type',
			'derivation_process',
			'derived_from',
			'description',
			'collection_to_dissociation_interval',
			'collection_to_dissociation_interval_units',
			'death_to_dissociation_interval',
			'death_to_dissociation_interval_units',
			'dissociation_reagent',
			'dissociation_time',
			'dissociation_time_units',
			'red_blood_cell_removal',
			'enrichment_factors',
			'enriched_cell_types',
			'cell_depletion_factors',
			'depleted_cell_types',
			'cell_viability_method',
			'percent_cell_viability',
			'single_cell_isolation_method',
			'single_cell_isolation_method_version',
			'dbxrefs',
			'feature_antibodies',
			'treatments',
			'starting_quantity',
			'starting_quantity_units'
		],
		'preferred': [
			'description',
			'collection_to_dissociation_interval',
			'collection_to_dissociation_interval_units',
			'death_to_dissociation_interval',
			'death_to_dissociation_interval_units',
			'dissociation_reagent',
			'dissociation_time',
			'dissociation_time_units',
			'red_blood_cell_removal',
			'enrichment_factors',
			'enriched_cell_types',
			'cell_depletion_factors',
			'depleted_cell_types',
			'cell_viability_method',
			'percent_cell_viability',
			'single_cell_isolation_method',
			'single_cell_isolation_method_version',
			'dbxrefs',
		]
	},
	'tissue_section': {
		'order': [],
		'preferred': []
	},
	'library': {
		'order': [
			'derived_from',
			'protocol',
			'dataset',
			'lab',
			'description',
			'starting_quantity',
			'starting_quantity_units',
			'date_constructed',
			'dbxrefs'
		],
		'preferred': [
			'description',
			'starting_quantity',
			'starting_quantity_units',
			'date_constructed',
			'dbxrefs'
		]
	},
	'sequencing_run': {
		'order': [
			'derived_from',
			'notes',
			'documents',
			'dbxrefs',
			'demultiplexed_link'
		],
		'preferred': [
			'dbxrefs'
		]
	},
	'raw_sequence_file': {
		'order': [
			'derived_from',
			'dataset',
			's3_uri',
			'read_type',
			'lab'
		],
		'preferred': [
			's3_uri',
			'read_type',
			'lab'
		]
	},
	'raw_matrix_file': {
		'order': [
			'background_barcodes_included',
			'derivation_process',
			'derived_from',
			'output_types',
			'dataset',
			'lab',
			'assembly',
			'genome_annotation',
			'cellranger_assay_chemistry',
			'software',
			's3_uri',
			'intronic_reads_counted',
			'feature_keys',
			'value_units'
		],
		'preferred': [
			'lab',
			'assembly',
			'genome_annotation',
			'cellranger_assay_chemistry',
			'software',
			's3_uri',
			'intronic_reads_counted',
			'feature_keys',
			'value_units'
		]
	},
	'processed_matrix_file': {
		'order': [],
		'preferred': []
	},
	'cell_annotation': {
		'order': [],
		'preferred': []
	},
	'library_protocol': {
		'order': [],
		'preferred': []
	},
	'treatment': {
		'order': [
			'treatment_type',
			'treatment_term_name',
			'treatment_term_id'
		],
		'preferred': [
			'treatment_term_id'
		]
	},
	'antibody': {
		'order': [],
		'preferred': []
	},
	'target': {
		'order': [],
		'preferred': []
	},
	'lab': {
		'order': [],
		'preferred': []
	},
	'user': {
		'order': [
			'first_name',
			'last_name',
			'institute_name',
			'email'
		]
	}
}
