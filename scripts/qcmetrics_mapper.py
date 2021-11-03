cellranger = {
	'schema_mapping': {
		"valid_barcodes": "percent_valid_barcodes",
		"sequencing_saturation": "percent_sequencing_saturation",
		"q30_bases_in_barcode": "percent_q30_bases_in_barcode",
		"q30_bases_in_rna_read": "percent_q30_bases_in_rna_read",
		"q30_bases_in_rna_read_2": "percent_q30_bases_in_rna_read_2",
		"q30_bases_in_sample_index": "percent_q30_bases_in_sample_index",
		"q30_bases_in_umi": "percent_q30_bases_in_umi",
		"reads_mapped_to_genome": "percent_reads_mapped_to_genome",
		"reads_mapped_confidently_to_genome": "percent_reads_mapped_confidently_to_genome",
		"reads_mapped_confidently_to_intergenic_regions": "percent_reads_mapped_confidently_to_intergenic_regions",
		"reads_mapped_confidently_to_intronic_regions": "percent_reads_mapped_confidently_to_intronic_regions",
		"reads_mapped_confidently_to_exonic_regions": "percent_reads_mapped_confidently_to_exonic_regions",
		"reads_mapped_confidently_to_transcriptome": "percent_reads_mapped_confidently_to_transcriptome",
		"reads_mapped_antisense_to_gene": "percent_reads_mapped_antisense_to_gene",
		"fraction_reads_in_cells": "percent_reads_in_cells",
		"frac_valid_barcode": "frac_valid_barcodes",
		"frac_waste_duplicate": "frac_waste_dup",
		"estimated_fraction_cells_annotated": "estimated_frac_cells_annotated",
		"fraction_cell_calling_noise": "frac_cell_calling_noise",
		"fraction_gelbead_doublets_cells": "frac_gelbead_doublets_cells",
		"fraction_of_genome_within_2000bp_of_peaks": "frac_of_genome_within_2000bp_of_peaks",
		"bc_q30_bases_fract": "frac_q30_bases_in_barcode",
		"r1_q30_bases_fract": "frac_q30_bases_in_read1",
		"r2_q30_bases_fract": "frac_q30_bases_in_read2",
		"si_q30_bases_fract": "frac_q30_bases_in_sample_index",
		"median_frags_overlapping_peaks_per_cell": "median_fragments_overlapping_peaks_per_cell",
		"frac_waste_mitochondrial": "frac_waste_mito",
		"waste_non_cell_barcode_fragments": "waste_noncell_barcode_fragments",
		"frac_waste_non_cell_barcode": "frac_waste_noncell_barcode",
		"frac_valid_noncell": "frac_valid_noncells",
		"num_fragments": "number_fragments",
		"num_reads": "number_reads",
		"r1_tot_bases": "total_bases_in_read1",
		"r2_tot_bases": "total_bases_in_read2",
		"si_tot_bases": "total_bases_in_sample_index",
		"bc_tot_bases": "total_bases_in_barcode",
		"si_q30_bases": "q30_bases_in_sample_index",
		"r2_q30_bases": "q30_bases_in_read2",
		"r1_q30_bases": "q30_bases_in_read1",
		"bc_q30_bases": "q30_bases_in_barcode",
		"waste_duplicate_fragments": "waste_dup_fragments",
		"waste_mitochondrial_fragments": "waste_mito_fragments"
	},
	'value_mapping': {
		"Single Cell 3' v1": "SC3Pv1",
		"Single Cell 3' v2": "SC3Pv2",
		"Single Cell 3' v3": "SC3Pv3",
		"Single Cell 5' PE": "SC5P-PE",
		"Single Cell 5' R2-only": "SC5P-R2",
		"Single Cell 5' R1-only": "SC5P-R1"
	},
	'should_match': {
		"annotated_cells": "cells_detected",
		"waste_cell_dup_fragments": "waste_dup_fragments",
		"waste_cell_lowmapq_fragments": "waste_lowmapq_fragments",
		"waste_cell_mito_fragments": "waste_mito_fragments",
		"waste_cell_unmapped_fragments": "waste_unmapped_fragments"
	}
}

dragen = {
	'schema_mapping': {
		'Passing cells': 'estimated_number_of_cells',
		'Mean reads per cell': 'mean_reads_per_cell',
		'Median genes per cell': 'median_genes_per_cell',
		'Median UMIs per cell': 'median_umi_counts_per_cell',
		'Total input reads': 'number_of_reads',
		'Q30 bases R1 %': 'percent_q30_bases_in_rna_read',
		'Q30 bases R2 %': 'percent_q30_bases_in_rna_read_2',
		'Mapped reads %': 'percent_reads_mapped_to_genome',
		'Sequencing saturation': 'percent_sequencing_saturation',
		'Total genes detected': 'total_genes_detected'
	},
	'value_mapping': {
		'percent_sequencing_saturation': {
			'factor': 100,
			'action': 'multiply'
		}
	}
}


