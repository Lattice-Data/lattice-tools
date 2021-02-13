lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'uuid': 'provenance.document_id'
		},
	'HumanPostnatalDonor': {
		'class': 'donor_organism',
		'age': 'organism_age',
		'age_units': 'organism_age_unit.text',
		'alcohol_history': 'medical_history.alcohol_history',
		'body_mass_index': 'human_specific.body_mass_index',
		'cause_of_death': 'death.cause_of_death',
		'ethnicity.term_id': 'human_specific.ethnicity.ontology',
		'ethnicity.term_name': 'human_specific.ethnicity.ontology_label',
		'genotype': 'biomaterial_core.genotype',
		'height': 'height',
		'height_unit': 'height_unit.text',
		'life_stage': 'development_stage.ontology_label',
		'life_stage_term_id': 'development_stage.ontology',
		'living_at_sample_collection': 'is_living',
		'organism.ncbi_taxon': 'genus_species.ontology',
		'organism.scientific_name': 'genus_species.ontology_label',
		'sex': 'sex',
		'smoking_history': 'medical_history.smoking_history',
		'test_results': 'medical_history.test_results',
		'uuid': 'provenance.document_id',
		'weight': 'weight',
		'weight_unit': 'weight_unit.text'
		},
	'HumanPrenatalDonor': {
		'class': 'donor_organism',
		'ethnicity.term_id': 'human_specific.ethnicity.ontology',
		'ethnicity.term_name': 'human_specific.ethnicity.ontology_label',
		'gestational_age': 'organism_age',
		'gestational_age_units': 'organism_age_unit.text',
		'life_stage': 'development_stage.ontology_label',
		'life_stage_term_id': 'development_stage.ontology',
		'organism.ncbi_taxon': 'genus_species.ontology',
		'organism.scientific_name': 'genus_species.ontology_label',
		'uuid': 'provenance.document_id',
		'weight': 'weight',
		'weight_unit': 'weight_unit.text'
		},
	'MousePrenatalDonor': {
		'class': 'donor_organism',
		'gestational_age': 'organism_age',
		'gestational_age_units': 'organism_age_unit.text',
		'life_stage': 'development_stage.ontology_label',
		'life_stage_term_id': 'development_stage.ontology',
		'organism.ncbi_taxon': 'genus_species.ontology',
		'organism.scientific_name': 'genus_species.ontology_label',
		'strain_term_id': 'mouse_specific.strain.ontology',
		'strain_term_name': 'mouse_specific.strain.ontology_label',
		'uuid': 'provenance.document_id',
		'weight': 'weight',
		'weight_unit': 'weight_unit.text'
		},
	'MousePrenatalDonor': {
		'class': 'donor_organism',
		'age': 'organism_age',
		'age_units': 'organism_age_unit.text',
		'life_stage': 'development_stage.ontology_label',
		'life_stage_term_id': 'development_stage.ontology',
		'organism.ncbi_taxon': 'genus_species.ontology',
		'organism.scientific_name': 'genus_species.ontology_label',
		'strain_term_id': 'mouse_specific.strain.ontology',
		'strain_term_name': 'mouse_specific.strain.ontology_label',
		'uuid': 'provenance.document_id',
		'weight': 'weight',
		'weight_unit': 'weight_unit.text'
		},
	'Tissue': {
		'class': 'specimen_from_organism',
		'uuid': 'provenance.document_id'
		},
	'CellCulture': {
		'class': 'cell_line',
		'uuid': 'provenance.document_id'
		},
	'Organoid': {
		'class': 'organoid',
		'uuid': 'provenance.document_id'
		},
	'Suspension': {
		'class': 'cell_suspension',
		'uuid': 'provenance.document_id'
		},
	'Library': {
		'class': 'library_preparation_protocol',
		'uuid': 'provenance.document_id'
		},
	'RawSequenceFile': {
		'class': 'sequence_file',
		'file_format': 'file_core.format',
		'md5sum': 'file_core.checksum',
		'read_index': 'read_type',
		'read_length': 'read_length',
		'submitted_file_name': 'file_core.file_name',
		'uuid': 'provenance.document_id'
		},
	'SequencingRun': {
		'class': 'sequencing_protocol',
		'uuid': 'provenance.document_id'
		}
}
