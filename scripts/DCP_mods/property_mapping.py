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
			'NCBITaxon:9606': '9606',
			'NCBITaxon:10090': '10090'
		}
	},
	'development_stage.ontology': {
		'lattice': 'development_ontology.term_id'
	},
	'development_stage.ontology_label': {
		'lattice': 'development_ontology.term_name'
	},
	'development_stage.text': {
		'lattice': 'development_ontology.term_name'
	},
	'diseases': {
		'lattice': 'diseases',
		'subprop_map': {
			'ontology': {
				'lattice': 'term_id',
			},
			'ontology_label': {
				'lattice': 'term_name'
			},
			'text': {
				'lattice': 'term_name'
			}
		}
	},
	'genus_species.ontology': {
		'lattice': 'organism.taxon_id'
	},
	'genus_species.ontology_label': {
		'lattice': 'organism.taxon_id',
		'value_map': {
			'NCBITaxon:9606': 'Homo sapiens',
			'NCBITaxon:10090': 'Mus musculus'
		}
	},
	'genus_species.text': {
		'lattice': 'organism.taxon_id',
		'value_map': {
			'NCBITaxon:9606': 'Homo sapiens',
			'NCBITaxon:10090': 'Mus musculus'
		}
	},
	'medical_history.treatment': {
		'lattice': 'treatment_summary'
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
	},
	'human_specific.ethnicity.text': {
		'lattice': 'ethnicity.term_name'
	}
}

mouse_donor = {
	'mouse_specific.strain.ontology': {
		'lattice': 'strain_term_id'
	},
	'mouse_specific.strain.ontology_label': {
		'lattice': 'strain_term_name'
	},
	'mouse_specific.strain.text': {
		'lattice': 'strain_term_name'
	}
}

prenatal_donor = {
	'gestational_age': {
		'lattice': 'conceptional_age' # 2 weeks, or 14 days, are added later
	},
	'gestational_age_unit.text': {
		'lattice': 'conceptional_age_units'
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
	'biomaterial_core.ncbi_taxon_id': {
		'lattice': 'donors' # later pull organism.taxon_id
	},
	'dbxrefs': {
		'lattice': 'dbxrefs',
		'BioSample': 'biomaterial_core.biosamples_accession',
		'SRA': 'biomaterial_core.insdc_sample_accession'	
	},
	'provenance.document_id': {
		'lattice': 'uuid'
	},
	'treatment_summary': {
		'lattice': 'treatment_summary'
	}
}

lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'contributors': {
			'lattice': 'contributors',
			'subprop_map': {
				'contact_name': {
					'lattice': 'title'
				},
				'institution': {
					'lattice': 'institute_name'
				},
				'project_role': {
					'lattice': 'job_title'
				}
			}
		},
		'corresponding_contributors': { # DCP_mapper merges these with contributors
			'lattice': 'corresponding_contributors',
			'subprop_map': {
				'contact_name': {
					'lattice': 'title'
				},
				'institution': {
					'lattice': 'institute_name'
				},
				'email': {
					'lattice': 'email'
				},
				'project_role': {
					'lattice': 'job_title'
				}
			}
		},
		'dbxrefs': {
			'lattice': 'dbxrefs',
			'SRA': 'insdc_project',
			'GEO': 'geo_series',
			'ArrayExpress': 'array_express_investigation',
			'BioProject': 'insdc_study'
		},
		'funders.funder_name': {
			'lattice': 'funding_organizations'
		},
		'project_core.project_description': {
			'lattice': 'description'
		},
		'project_core.project_title': {
			'lattice': 'dataset_title'
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
				'doi': {
					'lattice': 'doi'
				},
				'pmid': {
					'lattice': 'pmid'
				},
				'publication_title': {
					'lattice': 'title'
				}
			}
		},
		'supplementary_links': {
			'lattice': 'urls'
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
		'biomaterial_core.biomaterial_description': {
			'lattice': 'summary'
		},
		'diseases': {
			'lattice': 'diseases',
			'subprop_map': {
				'ontology': {
					'lattice': 'term_id',
				},
				'ontology_label': {
					'lattice': 'term_name'
				},
				'text': {
					'lattice': 'term_name'
				}
			}
		},
		'organ.text': {
			'lattice': 'biosample_ontology.organ_slims'
		},
		'organ_parts.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'organ_parts.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		},
		'organ_parts.text': {
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
				'cryopreservation': 'frozen at -80C',
				'flash-freezing': 'frozen in liquid nitrogen',
				'paraffin embedding': 'paraffin block'
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
		'spatial_information': {
			'lattice': 'spatial_information'
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
		'state_of_specimen.ischemic_time_units': { # deleted before report
			'lattice': 'ischemic_time_units'
		},
		'state_of_specimen.postmortem_interval': {
			'lattice': 'death_to_preservation_interval'
		},
		'state_of_specimen.postmortem_interval_units': { # deleted before report
			'lattice': 'death_to_preservation_interval_units'
		}
	},
	'CellCulture': {
		**biosample,
		'class': 'cell_line',
		'biomaterial_core.biomaterial_description': {
			'lattice': 'summary'
		},
		'cell_type.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'cell_type.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		},
		'cell_type.text': {
			'lattice': 'biosample_ontology.term_name'
		},
		'growth_conditions.growth_medium': {
			'lattice': 'growth_medium'
		},
		'growth_conditions.mycoplasma_testing_method': {
			'lattice': 'mycoplasma_testing_method'
		},
		'growth_conditions.mycoplasma_testing_results': {
			'lattice': 'mycoplasma_testing_results',
			'value_map': {
				'negative': 'pass',
				'positive': 'fail'
			}
		},
		'lot_number': {
			'lattice': 'lot_id'
		},
		'supplier': {
			'lattice': 'source'
		},
		'tissue.text': {
			'lattice': 'biosample_ontology.organ_slims'
		},
		'cell_line_type': {
			'lattice': 'biosample_ontology.cell_slims'
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
		'biomaterial_core.biomaterial_description': {
			'lattice': 'summary'
		},
		'embedded_in_matrigel': {
			'lattice': 'embedded_in_matrigel'
		},
		'growth_environment': {
			'lattice': 'growth_medium'
		},
		'model_organ.text': {
			'lattice': 'biosample_ontology.organ_slims'
		},
		'model_organ_part.ontology': {
			'lattice': 'biosample_ontology.term_id'
		},
		'model_organ_part.ontology_label': {
			'lattice': 'biosample_ontology.term_name'
		},
		'model_organ_part.text': {
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
		'dissociation_reagent': {
			'lattice': 'dissociation_reagent'
		},
		'dissociation_time': {
			'lattice': 'dissociation_time'
		},
		'dissociation_time_units': {
			'lattice': 'dissociation_time_units'
		},
		'enrichment_factors': {
			'lattice': 'enrichment_factors'
		},
		'estimated_cell_count': {
			'lattice': 'starting_quantity'
		},
		'estimated_count_units': { # deleted before report
			'lattice': 'starting_quantity_units'
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
				},
				'text': {
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
			'lattice': 'protocol.end_bias',
			'value_map': {
				"3'": "3 prime end bias",
				"5'": "5 prime end bias",
				"full length": "full length"
			}
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
			'lattice': 'assay',
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
				'poly(dT)': 'poly-dT',
				'adapter': 'random'
			}
		},
		'protocol_core.protocol_id': {
			'lattice': 'protocol.name'
		},
		'provenance.document_id': {
			'lattice': 'protocol.uuid'
		},
		'strand': {
			'lattice': 'protocol.strand_specificity'
		}
	},
	'RawSequenceFile': {
		'class': 'sequence_file',
		'crc32c': {
			'lattice': 'crc32c'
		},
		'file_core.checksum': {
			'lattice': 'md5sum'
		},
		'file_core.file_name': {
			'lattice': 'submitted_file_name'
		},
		'file_core.format': {
			'lattice': 'file_format'
		},
		'file_size': {
			'lattice': 'file_size' # used for file_descriptor, not other metadata
		},
		'insdc_run_accessions': {
			'lattice': 'derived_from' # DCP_mapper pulls out just dbxrefs
		},
		'library_prep_id': {
			'lattice': 'libraries'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'read_index': {
			'lattice': 'read_type',
			'value_map': {
				'i5 index': 'index2',
				'i7 index': 'index1',
				'Read 1': 'read1',
				'Read 2': 'read2',
				'Read 1N': 'read1',
				'Read 2N': 'read2'
			}
		},
		'read_length': {
			'lattice': 'read_length'
		},
		'sha256': {
			'lattice': 'sha256'
		},
		's3_uri': {
			'lattice': 's3_uri' # used to transfer the file but deleted from the metadata
		},
		'external_uri': {
			'lattice': 'external_uri' # used to transfer the file but deleted from the metadata
		}
	},
	'Document': {
		'class': 'supplementary_file',
		'content_type': {
			'lattice': 'attachment.type' # used in file_descriptor
		},
		'file_core.file_name': {
			'lattice': 'attachment.download'
		},
		'file_description': {
			'lattice': 'description'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	}
}

donor_stages = {
	'MmusDv:0000002': 'embryonic mouse stage',
	'MmusDv:0000031': 'fetal stage',
	'MmusDv:0000036': 'Theiler stage 27',
	'MmusDv:0000112': 'premature stage',
	'MmusDv:0000110': 'mature stage',
	'HsapDv:0000002': 'embryonic human stage',
	'HsapDv:0000037': 'fetal stage',
	'HsapDv:0000082': 'newborn human stage',
	'HsapDv:0000083': 'infant stage',
	'HsapDv:0000081': 'child stage',
	'HsapDv:0000086': 'adolescent stage',
	'HsapDv:0000087': 'human adult stage'
}

dcp_versions = {
  'file_descriptor': '2.0.0',
  'links': '2.1.1',
  'cell_line': '9.0.5',
  'cell_suspension': '13.3.0',
  'donor_organism': '15.5.0',
  'imaged_specimen': '3.3.0',
  'organoid': '11.3.0',
  'specimen_from_organism': '10.4.0',
  'analysis_file': '6.2.0',
  'image_file': '2.2.0',
  'reference_file': '3.2.0',
  'sequence_file': '9.2.0',
  'supplementary_file': '2.2.0',
  'analysis_process': '12.0.0',
  'process': '9.2.0',
  'project': '9.0.4',
  'analysis_protocol': '9.1.0',
  'aggregate_generation_protocol': '2.1.0',
  'collection_protocol': '9.2.0',
  'differentiation_protocol': '2.2.0',
  'dissociation_protocol': '6.2.0',
  'enrichment_protocol': '3.1.0',
  'ipsc_induction_protocol': '3.2.0',
  'imaging_preparation_protocol': '2.2.0',
  'imaging_protocol': '11.2.0',
  'protocol': '7.1.0',
  'library_preparation_protocol': '6.2.0',
  'sequencing_protocol': '10.1.0'
}
