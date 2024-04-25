UNREPORTED_VALUE = 'unknown'
SCHEMA_VERSION = '5_1'
MTX_DIR = 'matrix_files'

# Reference files by which the flattener will filter var features
REF_FILES = {
	'ercc':'genes_ercc.csv',
	'human':'genes_homo_sapiens.csv',
	'mouse':'genes_mus_musculus.csv',
	'sars':'genes_sars_cov_2.csv'
}

# Metadata to be gathered for each object type
CELL_METADATA = {
	'donor': [
		'donor_id',
		'age_display',
		'sex',
		'ethnicity',
		'causes_of_death.term_name',
		'diseases.term_id',
		'diseases.term_name',
		'family_medical_history',
		'living_at_sample_collection',
		'menopausal_status',
		'organism.taxon_id',
		'risk_score_tyrer_cuzick_lifetime',
		'smoker',
		'times_pregnant'
		],
	'sample': [
		'age_development_stage_redundancy',
		'uuid',
		'preservation_method',
		'biosample_ontology.term_id',
		'biosample_ontology.organ_slims',
		'summary_development_ontology_at_collection.development_slims',
		'summary_development_ontology_at_collection.term_id',
		'derivation_process',
		'diseases.term_id',
		'diseases.term_name',
		'disease_state',
		'menstrual_phase_at_collection',
		'source',
		'summary_body_mass_index_at_collection',
		'treatment_summary',
		'growth_medium',
		'genetic_modifications',
		'@type'
		],
	'tissue_section': [
		'uuid',
		'thickness',
		'thickness_units'
	],
	'suspension': [
		'cell_depletion_factors',
		'depleted_cell_types.term_name',
		'derivation_process',
		'dissociation_reagent',
		'dissociation_time',
		'dissociation_time_units',
		'enriched_cell_types.term_name',
		'enrichment_factors',
		'percent_cell_viability',
		'uuid',
		'suspension_type',
		'tissue_handling_interval',
		'@id'
		],
	'library': [
		'uuid',
		'protocol.assay_ontology.term_id',
		'starting_quantity',
		'starting_quantity_units',
		'@id'
	],
	'raw_matrix': [
		'assembly',
		'genome_annotation',
		'software'
	],
	'seq_run': [
		'platform'
	]
}

DATASET_METADATA = {
	'final_matrix': [
		'description',
		'default_embedding',
		'is_primary_data'
		]
	}

ANNOT_FIELDS = [
	'cell_ontology.term_id',
	'author_cell_type',
	'cell_state'
]

ANTIBODY_METADATA = {
	'antibody': [
		'oligo_sequence',
		'host_organism',
		'source',
		'product_ids',
		'clone_id',
		'control',
		'isotype'
	],
	'target': [
		'label',
		'organism.scientific_name'
	]
}


# Mapping of field name (object_type + "_" + property) and what needs to be in the final cxg h5ad
PROP_MAP = {
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'sample_summary_development_ontology_at_collection_term_id': 'development_stage_ontology_term_id',
	'sample_age_development_stage_redundancy': 'donor_age_redundancy',
	'sample_disease_state': 'disease_state',
	'sample_summary_body_mass_index_at_collection': 'donor_BMI_at_collection',
	'sample_growth_medium': 'growth_medium',
	'sample_genetic_modifications': 'genetic_modifications',
	'sample_menstrual_phase_at_collection': 'menstrual_phase_at_collection',
	'library_protocol_assay_ontology_term_id': 'assay_ontology_term_id',
	'donor_sex': 'sex',
	'sample_@type': 'tissue_type',
	'donor_donor_id': 'donor_id',
	'donor_organism_taxon_id': 'organism_ontology_term_id',
	'donor_ethnicity': 'self_reported_ethnicity_ontology_term_id',
	'donor_age_display': 'donor_age',
	'donor_risk_score_tyrer_cuzick_lifetime': 'tyrer_cuzick_lifetime_risk',
	'donor_smoker': 'donor_smoking_status',
	'donor_causes_of_death_term_name': 'donor_cause_of_death',
	'matrix_description': 'title',
	'matrix_default_embedding': 'default_embedding',
	'matrix_is_primary_data': 'is_primary_data',
	'cell_annotation_author_cell_type': 'author_cell_type',
	'cell_annotation_cell_ontology_term_id': 'cell_type_ontology_term_id',
	'cell_annotation_cell_state': 'cell_state',
	'suspension_suspension_type': 'suspension_type',
	'suspension_enriched_cell_types_term_name': 'suspension_enriched_cell_types',
	'suspension_depleted_cell_types_term_name': 'suspension_depleted_cell_types',
	'suspension_cell_depletion_factors': 'suspension_depletion_factors',
	'suspension_tissue_handling_interval': 'tissue_handling_interval',
	'antibody_oligo_sequence': 'barcode',
	'antibody_source': 'vendor',
	'antibody_product_ids': 'vender_product_ids',
	'antibody_clone_id': 'clone_id',
	'antibody_isotype': 'isotype',
	'antibody_host_organism': 'host_organism',
	'target_organism_scientific_name': 'target_organism',
	'raw_matrix_software': 'alignment_software',
	'raw_matrix_genome_annotation': 'mapped_reference_annotation',
	'raw_matrix_assembly': 'mapped_reference_assembly',
	'seq_run_platform': 'sequencing_platform'
}