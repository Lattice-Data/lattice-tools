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
		'@type',
		'date_obtained'
		],
	'tissue_section': [
		'uuid',
		'thickness',
		'thickness_units'
	],
	'suspension': [
		'cell_depletion_factors',
		'depleted_cell_types.term_name',
		'depleted_cell_types.term_id',
		'derivation_process',
		'dissociation_reagent',
		'dissociation_time',
		'dissociation_time_units',
		'enriched_cell_types.term_name',
		'enriched_cell_types.term_id',
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
		'protocol.end_bias',
		'starting_quantity',
		'starting_quantity_units',
		'@id',
		'lab.institute_name',
		'dbxrefs'
	],
	'raw_matrix': [
		'assembly',
		'genome_annotation',
		'software',
		'intronic_reads_counted'
	],
	'seq_run': [
		'platform'
	],
	'raw_seq': [
		'flowcell_details'
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
	'sample_derivation_process': 'sample_collection_method',
	'sample_disease_state': 'disease_state',
	'sample_summary_body_mass_index_at_collection': 'donor_BMI_at_collection',
	'sample_growth_medium': 'growth_medium',
	'sample_genetic_modifications': 'genetic_modifications',
	'sample_menstrual_phase_at_collection': 'menstrual_phase_at_collection',
	'sample_source': 'tissue_source',
	'sample_date_obtained' : 'sample_collection_year',
	'library_protocol_assay_ontology_term_id': 'assay_ontology_term_id',
	'library_lab_institute_name': 'institute',
	'library_protocol_end_bias': 'sequenced_fragment',
	'library_dbxrefs' : 'library_id_repository',
	'library_starting_quantity':'cell_number_loaded',
	'library_starting_quantity_units':'cell_number_loaded_units',
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
	'suspension_enriched_cell_types_term_id': 'suspension_enriched_cell_terms',
	'suspension_depleted_cell_types_term_id': 'suspension_depleted_cell_terms',
	'suspension_cell_depletion_factors': 'suspension_depletion_factors',
	'suspension_tissue_handling_interval': 'tissue_handling_interval',
	'suspension_percent_cell_viability':'cell_viability_percentage',
	'antibody_oligo_sequence': 'barcode',
	'antibody_source': 'vendor',
	'antibody_product_ids': 'vender_product_ids',
	'antibody_clone_id': 'clone_id',
	'antibody_isotype': 'isotype',
	'antibody_host_organism': 'host_organism',
	'target_organism_scientific_name': 'target_organism',
	'raw_matrix_software': 'alignment_software',
	'raw_matrix_genome_annotation': 'gene_annotation_version',
	'raw_matrix_assembly': 'reference_genome',
	'raw_matrix_intronic_reads_counted':'intronic_reads_counted',
	'seq_run_platform': 'sequencing_platform'
	'raw_seq_flowcell_details': 'library_sequencing_run'
}

GENCODE_MAP = {
	'GENCODE 44': 'v110',
	'GENCODE 43': 'v109',
	'GENCODE 42': 'v108',
	'GENCODE 41': 'v107',
	'GENCODE 40': 'v106',
	'GENCODE 39': 'v105',
	'GENCODE 38': 'v104',
	'GENCODE 37': 'v103',
	'GENCODE 36': 'v102',
	'GENCODE 35': 'v101',
	'GENCODE 34': 'v100',
	'GENCODE 33': 'v99',
	'GENCODE 32': 'v98',
	'GENCODE 31': 'v97',
	'GENCODE 30': 'v96',
	'GENCODE 29': 'v94',
	'GENCODE 28': 'v92',
	'GENCODE 27': 'v90',
	'GENCODE 26': 'v88',
	'GENCODE 25': 'v85',
	'GENCODE 24': 'v83',
	'GENCODE 23': 'v81',
	'GENCODE 22': 'v79',
	'GENCODE 21': 'v77',
	'GENCODE 20': 'v76',
	'GENCODE 19': 'v74',
}

SAMPLE_COLLECTION_MAP = {
	'percutaneous biopsy': 'biopsy',
	'open biopsy': 'biopsy',
	'resection': 'surgical resection',
	'dissection': 'surgical resection',
	'swab': 'brush',
	'bronchoalveolar lavage': 'bodily fluid',
	'aspiration': 'biopsy',
	'density centrifugation': 'other',
	'enzymatic digestion': 'other',
	'cryosection': 'other'
}

SAMPLE_PRESERVATION_MAP = {
	'cryopreservation': 'frozen at -80C',
	'flash-freezing': 'frozen in liquid nitrogen',
	'n/a (fresh)': 'fresh',
	'paraffin embedding': 'paraffin block',
	'OCT embedding': 'frozen at -80C'
}

OPTIONAL_COLUMNS = [
	'alignment_software',
	'cell_state',
	'cell_viability_percentage',
	'cell_number_loaded',
	'cell_number_loaded_units',
	'disease_state',
	'donor_BMI_at_collection',
	'donor_cause_of_death',
	'donor_family_medical_history',
	'donor_living_at_sample_collection',
	'donor_menopausal_status',
	'donor_smoking_status',
	'donor_times_pregnant',
	'gene_annotation_version',
	'genetic_modifications',
	'growth_medium',
	'menstrual_phase_at_collection',
	'reference_genome',
	'reported_diseases',
	'sample_treatment_summary',
	'sample_collection_year',
	'sequencing_platform',
	'suspension_dissociation_reagent',
	'suspension_dissociation_time',
	'suspension_dissociation_time_units',
	'suspension_depleted_cell_types',
	'suspension_derivation_process',
	'suspension_enriched_cell_types',
	'suspension_enrichment_factors',
	'suspension_depletion_factors', 
	'suspension_uuid',
	'tissue_section_thickness',
	'tissue_section_thickness_units',
	'tissue_handling_interval',
	'tyrer_cuzick_lifetime_risk',
	'library_id_repository',
	'intronic_reads_counted'
]

COLUMNS_TO_DROP = [
	'author_donor_@id',
	'author_donor_x',
	'author_donor_y',
	'batch',
	'donor_age_redundancy',
	'donor_diseases_term_id',
	'donor_diseases_term_name',
	'library_@id_x',
	'library_@id_y',
	'library_authordonor',
	'library_donor_@id',
	'library_@id',
	'raw_matrix_accession',
	'sample_biosample_ontology_cell_slims',
	'sample_summary_development_ontology_at_collection_development_slims',
	'sample_diseases_term_id',
	'sample_diseases_term_name',
	'sample_biosample_ontology_organ_slims',
	'sex',
	'suspension_@id'
]


# Accepted accessions for library dbxrefs

ACCEPTED_ACCESSIONS = {
	'EGA:EGAX',
	'SRA:SRX',
	'ENA:ERX'
}

