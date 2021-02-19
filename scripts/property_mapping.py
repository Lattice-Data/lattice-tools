donor = {
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
	'development_stage.ontology': {
		'lattice': 'life_stage_term_id'
	},
	'development_stage.ontology_label': {
		'lattice': 'life_stage'
	},
	'diseases': {
		'lattice': 'diseases',
		'future_subprop_map': { # need disease embedded in Donor
			'ontology': {
				'lattice': 'term_id',
			},
			'ontology_label': {
				'lattice': 'term_name'
			}
		}
	},
	'genus_species.ontology': {
		'lattice': 'organism.taxon_id'
	},
	'genus_species.ontology_label': {
		'lattice': 'organism.scientific_name'
	},
	'provenance.document_id': {
		'lattice': 'uuid'
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
}

human_donor = {
	'human_specific.ethnicity.ontology': {
		'lattice': 'ethnicity.term_id'
	},
	'human_specific.ethnicity.ontology_label': {
		'lattice': 'ethnicity.term_name'
	}
}

mouse_donor = {
	'mouse_specific.strain.ontology': {
		'lattice': 'strain_term_id'
	},
	'mouse_specific.strain.ontology_label': {
		'lattice': 'strain_term_name'
	}
}

prenatal_donor = {
	'organism_age': {
		'lattice': 'gestational_age'
	},
	'organism_age_unit.text': {
		'lattice': 'gestational_age_units'
	}
}

postnatal_donor = {
	'organism_age': {
		'lattice': 'age'
	},
	'organism_age_unit.text': {
		'lattice': 'age_units'
	}
}

biosample = {
	'biomaterial_core.biomaterial_id': {
		'lattice': 'uuid'
	},
	'dbxrefs': {
		'lattice': 'dbxrefs',
		'BioSample': 'biomaterial_core.biosamples_accession',
		'SRA': 'biomaterial_core.insdc_sample_accession'	
	},
	'provenance.document_id': {
		'lattice': 'uuid'
	}
}

lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'dbxrefs': {
			'lattice': 'dbxrefs',
			'SRA': 'insdc_project_accessions',
			'GEO': 'geo_series_accessions',
			'ArrayExpress': 'array_express_accessions',
			'BioProject': 'insdc_study_accessions',
			'BioStudies': 'biostudies_accessions'
		},
		'project_core.project_description': {
			'lattice': 'description'
		},
		'project_core.project_title': {
			'lattice': 'title'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'publications': {
			'lattice': 'references',
			'subprop_map': {
				'authors': {
					'lattice': 'authors',
				},
				'title': {
					'lattice': 'title'
				},
				'dbxrefs': {
					'lattice': 'identifiers',
					'doi': 'doi',
					'PMID': 'pmid'
				}
			}
		}
	},
	'HumanPostnatalDonor': {
		**donor,
		**human_donor,
		**postnatal_donor,
		'death.cause_of_death': {
			'lattice': 'cause_of_death'
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
		}
	},
	'HumanPrenatalDonor': {
		**donor,
		**human_donor,
		**prenatal_donor
	},
	'MousePostnatalDonor': {
		**donor,
		**mouse_donor,
		**postnatal_donor
	},
	'MousePrenatalDonor': {
		**donor,
		**mouse_donor,
		**prenatal_donor
	},
	'Tissue': {
		**biosample,
		'class': 'specimen_from_organism',
		'diseases': {
			'lattice': 'diseases',
			'future_subprop_map': {
				'ontology': {
					'lattice': 'term_id',
				},
				'ontology_label': {
					'lattice': 'term_name'
				}
			}
		},
		'organ_parts.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'organ_parts.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		},
		'preservation_storage.preservation_method': {
			'lattice': 'preservation_method',
			'value_map': {
				'n/a (fresh)': 'fresh',
				'paraffin embedding': 'paraffin block'
			}
		},
		'preservation_storage.storage_method': {
			'lattice': 'preservation_method',
			'value_map': {
				'n/a (fresh)': 'fresh',
				'cryopreservation': 'cryopreservation, other', # confirm with DCP
				'paraffin embedding': 'formalin fixed and paraffin embedded'
			}
		},
		'preservation_storage.storage_time': {
			'lattice': 'preservation_time'
		},
		'preservation_storage.storage_time_unit.text': {
			'lattice': 'preservation_time_units'
		},
		'purchased_specimen.catalog_number': {
			'lattice': 'product_id'
		},
		'purchased_specimen.lot_number': {
			'lattice': 'lot_id'
		},
		'purchased_specimen.manufacturer': {
			'lattice': 'source'
		},
		'state_of_specimen.ischemic_temperature': {
			'lattice': 'ischemic_temperature',
			'value_map': {
				'warm': 'warm',
				'cold': 'cold'
			}
		},
		'state_of_specimen.ischemic_time': {
			'lattice': 'ischemic_time'
		},
		'state_of_specimen.postmortem_interval': {
			'lattice': 'death_to_preservation_interval'
		}
	},
	'CellCulture': {
		**biosample,
		'class': 'cell_line',
		'catalog_number': {
			'lattice': 'product_id'
		},
		'cell_type.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'cell_type.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		},
		'diseases': {
			'lattice': 'diseases',
			'future_subprop_map': {
				'ontology': {
					'lattice': 'term_id',
				},
				'ontology_label': {
					'lattice': 'term_name'
				}
			}
		},
		'growth_conditions.growth_medium': {
			'lattice': 'growth_medium'
		},
		'growth_conditions.mycoplasma_testing_method': {
			'lattice': 'mycoplasma_testing_method'
		},
		'growth_conditions.mycoplasma_testing_results': {
			'lattice': 'mycoplasma_testing_results'
		},
		'lot_number': {
			'lattice': 'lot_id'
		},
		'supplier': {
			'lattice': 'source'
		}
	},
	'Organoid': {
		**biosample,
		'class': 'organoid',
		'age': {
			'lattice': 'post_differentiation_time'
		},
		'age_unit.text': {
			'lattice': 'post_differentiation_time_units'
		},
		'embedded_in_matrigel': {
			'lattice': 'embedded_in_matrigel'
		},
		'growth_environment': {
			'lattice': 'growth_medium'
		},
		'model_organ_part.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'model_organ_part.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		}
	},
	'Suspension': {
		**biosample,
		'class': 'cell_suspension',
		'cell_morphology.cell_size': {
			'lattice': 'cell_size'
		},
		'cell_morphology.cell_size_unit.text': {
			'lattice': 'cell_size_unit'
		},
		'cell_morphology.cell_viability_method': {
			'lattice': 'cell_viability_method'
		},
		'cell_morphology.percent_cell_viability': {
			'lattice': 'percent_cell_viability'
		},
		'growth_conditions.culture_environment': {
			'lattice': 'medium'
		},
		'selected_cell_types': {
			'lattice': 'enriched_cell_types',
			'subprop_map': {
				'ontology': {
					'lattice': 'term_id'
				},
				'ontology_label': {
					'lattice': 'term_name'
				}
			}
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
		'nucleic_acid_source': {
			'lattice': 'assay_type',
			'value_map': {
				'scATAC-seq': 'single cell',
				'snATAC-seq': 'single nucleus',
				'scRNA-seq': 'single cell',
				'snRNA-seq': 'single nucleus',
				'CITE-seq': 'single cell',
				'bulk ATAC-seq': 'bulk',
				'bulk RNA-seq': 'bulk'
			}
		},
		'primer': {
			'lattice': 'protocol.first_strand_primer',
			'value_map': {
				'poly(dT)': 'poly-dT'
			}
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
				'Read 1N': 'read 1', # need to confirm
				'Read 2N': 'read 2' # need to confirm
			}
		},
		'read_length': {
			'lattice': 'read_length'
		}
	},
	'SequencingRun': {
		'class': 'sequencing_protocol',
		'instrument_manufacturer_model.text': {
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
