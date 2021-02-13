lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'publications': {
			'lattice': 'references'
		}
	},
	'HumanPostnatalDonor': {
		'class': 'donor_organism',
		'biomaterial_core.biomaterial_id': {
			'lattice': 'uuid'
		},
		'biomaterial_core.genotype': {
			'lattice': 'genotype'
		},
		'biomaterial_core.ncbi_taxon_id': {
			'lattice': 'organism.taxon_id',
			'value_map': {
				'NCBI:9606': '9606',
				'NCBI:10090': '10090'
			}
		},
		'death.cause_of_death': {
			'lattice': 'cause_of_death'
		},
		'development_stage.ontology': {
			'lattice': 'life_stage_term_id'
		},
		'development_stage.ontology_label': {
			'lattice': 'life_stage'
		},
		'genus_species.ontology': {
			'lattice': 'organism.taxon_id'
		},
		'genus_species.ontology_label': {
			'lattice': 'organism.scientific_name'
		},
		'height': {
			'lattice': 'height'
		},
		'height_unit.text': {
			'lattice': 'height_unit'
		},
		'human_specific.body_mass_index': {
			'lattice': 'body_mass_index'
		},
		'human_specific.ethnicity.ontology': {
			'lattice': 'ethnicity.term_id'
		},
		'human_specific.ethnicity.ontology_label': {
			'lattice': 'ethnicity.term_name'
		},
		'is_living': {
			'lattice': 'living_at_sample_collection',
			'value_map': {
				True: 'yes',
				False: 'no'
			}
		},
		'medical_history.alcohol_history': {
			'lattice': 'alcohol_history'
		},
		'medical_history.smoking_history': {
			'lattice': 'smoking_history'
		},
		'medical_history.test_results': {
			'lattice': 'test_results'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'organism_age': {
			'lattice': 'age'
		},
		'organism_age_unit.text': {
			'lattice': 'age_units'
		},
		'sex': {
			'lattice': 'sex'
		},
		'weight': {
			'lattice': 'weight'
		},
		'weight_unit.text': {
			'lattice': 'weight_unit'
		}
	},
	'HumanPrenatalDonor': {
		'class': 'donor_organism',
		'biomaterial_core.biomaterial_id': {
			'lattice': 'uuid'
		},
		'biomaterial_core.ncbi_taxon_id': {
			'lattice': 'organism.taxon_id',
			'value_map': {
				'NCBI:9606': '9606',
				'NCBI:10090': '10090'
			}
		},
		'development_stage.ontology': {
			'lattice': 'life_stage_term_id'
		},
		'development_stage.ontology_label': {
			'lattice': 'life_stage'
		},
		'genus_species.ontology': {
			'lattice': 'organism.taxon_id'
		},
		'genus_species.ontology_label': {
			'lattice': 'organism.scientific_name'
		},
		'human_specific.ethnicity.ontology': {
			'lattice': 'ethnicity.term_id'
		},
		'human_specific.ethnicity.ontology_label': {
			'lattice': 'ethnicity.term_name'
		},
		'organism_age': {
			'lattice': 'gestational_age'
		},
		'organism_age_unit.text': {
			'lattice': 'gestational_age_units'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'weight': {
			'lattice': 'weight'
		},
		'weight_unit.text': {
			'lattice': 'weight_unit'
		}
	},
	'MousePrenatalDonor': {
		'class': 'donor_organism',
		'development_stage.ontology': {
			'lattice': 'life_stage_term_id'
		},
		'development_stage.ontology_label': {
			'lattice': 'life_stage'
		},
		'genus_species.ontology': {
			'lattice': 'organism.taxon_id'
		},
		'genus_species.ontology_label': {
			'lattice': 'organism.scientific_name'
		},
		'organism_age': {
			'lattice': 'gestational_age'
		},
		'organism_age_unit.text': {
			'lattice': 'gestational_age_units'
		},
		'mouse_specific.strain.ontology': {
			'lattice': 'strain_term_id'
		},
		'mouse_specific.strain.ontology_label': {
			'lattice': 'strain_term_name'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'weight': {
			'lattice': 'weight'
		},
		'weight_unit.text': {
			'lattice': 'weight_unit'
		}
	},
	'MousePostnatalDonor': {
		'class': 'donor_organism',
		'development_stage.ontology': {
			'lattice': 'life_stage_term_id'
		},
		'development_stage.ontology_label': {
			'lattice': 'life_stage'
		},
		'genus_species.ontology': {
			'lattice': 'organism.taxon_id'
		},
		'genus_species.ontology_label': {
			'lattice': 'organism.scientific_name'
		},
		'organism_age': {
			'lattice': 'age'
		},
		'organism_age_unit.text': {
			'lattice': 'age_units'
		},
		'mouse_specific.strain.ontology': {
			'lattice': 'strain_term_id'
		},
		'mouse_specific.strain.ontology_label': {
			'lattice': 'strain_term_name'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'weight': {
			'lattice': 'weight'
		},
		'weight_unit.text': {
			'lattice': 'weight_unit'
		}
	},
	'Tissue': {
		'class': 'specimen_from_organism',
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	},
	'CellCulture': {
		'class': 'cell_line',
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	},
	'Organoid': {
		'class': 'organoid',
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	},
	'Suspension': {
		'class': 'cell_suspension',
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	},
	'Library': {
		'class': 'library_preparation_protocol',
		'cdna_library_amplification_method.text': {
			'lattice': 'protocol.library_amplification_method'
		},
		'end_bias': {
			'lattice': 'protocol.end_bias'
		},
		'input_nucleic_acid_molecule.text': {
			'lattice': 'protocol.biological_macromolecule'
		},
		'library_construction_method.text': {
			'lattice': 'protocol.title'
		},
		'library_preamplification_method.text': {
			'lattice': 'protocol.library_preamplification_method'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'strand': {
			'lattice': 'protocol.strand_specificity'
		}
	},
	'RawSequenceFile': {
		'class': 'sequence_file',
		'file_core.checksum': {
			'lattice': 'md5sum'
		},
		'file_core.file_name': {
			'lattice': 'submitted_file_name'
		},
		'file_core.format': {
			'lattice': 'file_format'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'read_index': {
			'lattice': 'read_type',
			'value_map': {
				'Read 1': 'read 1',
				'Read 2': 'read 2',
				'Read 1N': 'read 1',
				'Read 2N': 'read 2'
			}
		},
		'read_length': {
			'lattice': 'read_length'
		}
	},
	'SequencingRun': {
		'class': 'sequencing_protocol',
		'instrument_manufacturer.text': {
			'lattice': 'platform'
		},
		'local_machine_name': {
			'lattice': 'flowcell_details.machine'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	}
}
