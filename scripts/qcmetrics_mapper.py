mappings = {
	"cellranger": {
		"value_mapping": {
			"Single Cell 3' v1": "SC3Pv1",
			"Single Cell 3' v2": "SC3Pv2",
			"Single Cell 3' v3": "SC3Pv3",
			"Single Cell 5' PE": "SC5P-PE",
			"Single Cell 5' R2-only": "SC5P-R2",
			"Single Cell 5' R1-only": "SC5P-R1",
			"Single Cell Multiome ATAC + Gene Expression v1": "SCMuv1",
			"fillInSpatial": "Sp3Pv1"
		},
		"rna" : {
			"schema_mapping": {
				"q30_bases_in_barcode": "frac_q30_bases_in_barcode",
				"q30_bases_in_rna_read": "frac_q30_bases_in_rna_read",
				"q30_bases_in_read_2": "frac_q30_bases_in_rna_read2",
				"q30_bases_in_rna_read_2": "frac_q30_bases_in_rna_read2",
				"q30_bases_in_sample_index": "frac_q30_bases_in_sample_index",
				"q30_bases_in_sample_index_i1": "frac_q30_bases_in_sample_index",
				"q30_bases_in_sample_index_i2": "frac_q30_bases_in_sample_index2",
				"q30_bases_in_umi": "frac_q30_bases_in_umi",
				"fraction_reads_in_cells": "frac_reads_in_cells",
				"reads_mapped_antisense_to_gene": "frac_reads_mapped_antisense_to_gene",
				"reads_mapped_confidently_to_exonic_regions": "frac_reads_mapped_confidently_to_exonic_regions",
				"reads_mapped_confidently_to_genome": "frac_reads_mapped_confidently_to_genome",
				"reads_mapped_confidently_to_intergenic_regions": "frac_reads_mapped_confidently_to_intergenic_regions",
				"reads_mapped_confidently_to_intronic_regions": "frac_reads_mapped_confidently_to_intronic_regions",
				"reads_mapped_confidently_to_transcriptome": "frac_reads_mapped_confidently_to_transcriptome",
				"reads_mapped_to_genome": "frac_reads_mapped_to_genome",
				"reads_with_tso": "frac_reads_with_tso",
				"sequencing_saturation": "frac_sequencing_saturation",
				"fraction_of_transcriptomic_reads_in_cells": "frac_transcriptomic_reads_in_cells",
				"valid_barcodes": "frac_valid_barcodes",
				"valid_umis": "frac_valid_umis",
				"percent_duplicates": "frac_waste_dup",
				"mean_raw_reads_per_cell": "mean_reads_per_cell",
				"estimated_number_of_cells": "total_cells_detected",
				"number_of_reads": "total_reads",
				"sequenced_read_pairs": "total_reads"
			},
			"perc_to_frac": [
				"frac_q30_bases_in_barcode",
				"frac_q30_bases_in_rna_read",
				"frac_q30_bases_in_rna_read2",
				"frac_q30_bases_in_sample_index",
				"frac_q30_bases_in_umi",
				"frac_reads_in_cells",
				"frac_reads_mapped_antisense_to_gene",
				"frac_reads_mapped_confidently_to_exonic_regions",
				"frac_reads_mapped_confidently_to_genome",
				"frac_reads_mapped_confidently_to_intergenic_regions",
				"frac_reads_mapped_confidently_to_intronic_regions",
				"frac_reads_mapped_confidently_to_transcriptome",
				"frac_reads_mapped_to_genome",
				"frac_sequencing_saturation",
				"frac_valid_barcodes",
				"frac_q30_bases_in_antibody_read",
				"frac_reads_in_aggregate_barcodes",
				"frac_reads_recognized_antibody",
				"frac_reads_unrecognized_antibody",
				"frac_reads_usable"
			]
		},
		"atac": {
			"schema_mapping": {
				"bulk_unique_fragments_at_30000000_reads": "bulk_total_unique_fragments_at_30000000_reads",
				"estimated_fraction_cells_annotated": "estimated_frac_cells_annotated",
				"estimated_cells_present": "estimated_total_cells_present",
				"fraction_cell_calling_noise": "frac_cell_calling_noise",
				"fraction_of_transposition_events_in_peaks_in_cells": "frac_cut_fragments_in_peaks_in_cells",
				"fraction_gelbead_doublets_cells": "frac_gelbead_doublets_cells",
				"fraction_of_genome_in_peaks": "frac_genome_in_peaks",
				"fraction_of_genome_within_2000bp_of_peaks": "frac_genome_within_2000bp_of_peaks",
				"fraction_of_high-quality_fragments_in_cells": "frac_high_quality_fragments_in_cells",
				"fraction_of_high-quality_fragments_overlapping_peaks": "frac_high_quality_fragments_overlapping_peaks",
				"fraction_of_high-quality_fragments_overlapping_tss": "frac_high_quality_fragments_overlapping_tss",
				"bc_q30_bases_fract": "frac_q30_bases_in_barcode",
				"r1_q30_bases_fract": "frac_q30_bases_in_read1",
				"q30_bases_in_read_1": "frac_q30_bases_in_read1",
				"q30_bases_in_read_2": "frac_q30_bases_in_read2",
				"r2_q30_bases_fract": "frac_q30_bases_in_read2",
				"q30_bases_in_sample_index_i1": "frac_q30_bases_in_sample_index",
				"si_q30_bases_fract": "frac_q30_bases_in_sample_index",
				"confidently_mapped_read_pairs": "frac_reads_mapped_confidently_to_genome",
				"frac_mapped_confidently": "frac_reads_mapped_confidently_to_genome",
				"frac_valid_barcode": "frac_valid_barcodes",
				"valid_barcodes": "frac_valid_barcodes",
				"frac_valid_noncell": "frac_valid_noncells",
				"frac_waste_duplicate": "frac_waste_dup",
				"percent_duplicates": "frac_waste_dup",
				"frac_waste_mitochondrial": "frac_waste_mito",
				"frac_waste_non_cell_barcode": "frac_waste_noncell_barcode",
				"non-nuclear_read_pairs": "frac_waste_nonnuclear",
				"unmapped_read_pairs": "frac_waste_unmapped",
				"mean_raw_read_pairs_per_cell": "mean_fragments_per_cell",
				"median_frags_overlapping_peaks_per_cell": "median_fragments_overlapping_peaks_per_cell",
				"median_high-quality_fragments_per_cell": "median_high_quality_fragments_per_cell",
				"median_per_cell_unique_fragments_at_10000_hq_rpc": "median_per_cell_unique_fragments_at_10000_HQ_RPC",
				"median_per_cell_unique_fragments_at_10000_rrpc": "median_per_cell_unique_fragments_at_10000_RRPC",
				"median_per_cell_unique_fragments_at_30000_hq_rpc": "median_per_cell_unique_fragments_at_30000_HQ_RPC",
				"median_per_cell_unique_fragments_at_30000_rrpc": "median_per_cell_unique_fragments_at_30000_RRPC",
				"median_per_cell_unique_fragments_at_50000_rrpc": "median_per_cell_unique_fragments_at_50000_RRPC",
				"bc_q30_bases": "q30_bases_in_barcode",
				"r1_q30_bases": "q30_bases_in_read1",
				"r2_q30_bases": "q30_bases_in_read2",
				"si_q30_bases": "q30_bases_in_sample_index",
				"bc_tot_bases": "total_bases_in_barcode",
				"bases_in_peaks": "total_bases_in_peaks",
				"r1_tot_bases": "total_bases_in_read1",
				"r2_tot_bases": "total_bases_in_read2",
				"si_tot_bases": "total_bases_in_sample_index",
				"cells_detected": "total_cells_detected",
				"num_fragments": "total_fragments",
				"sequenced_read_pairs": "total_fragments",
				"number_of_low_targeting_barcodes": "total_low_targeting_barcodes",
				"number_of_peaks": "total_peaks_detected",
				"num_reads": "total_reads",
				"waste_duplicate_fragments": "waste_dup_fragments",
				"waste_mitochondrial_fragments": "waste_mito_fragments",
				"waste_non_cell_barcode_fragments": "waste_noncell_barcode_fragments"
			},
			"should_match": {
				"annotated_cells": "total_cells_detected",
				"waste_cell_dup_fragments": "waste_dup_fragments",
				"waste_cell_lowmapq_fragments": "waste_lowmapq_fragments",
				"waste_cell_mito_fragments": "waste_mito_fragments",
				"waste_cell_unmapped_fragments": "waste_unmapped_fragments"
				}
		},
		"antibody_capture": {
			"schema_mapping": {
				"q30_bases_in_antibody_read": "frac_q30_bases_in_antibody_read",
				"q30_bases_in_antibody_read_2": "frac_q30_bases_in_antibody_read2",
				"q30_bases_in_barcode": "frac_q30_bases_in_barcode",
				"q30_bases_in_umi": "frac_q30_bases_in_umi",
				"fraction_antibody_reads_in_aggregate_barcodes": "frac_reads_in_aggregate_barcodes",
				"fraction_reads_in_barcodes_with_high_umi_counts": "frac_reads_in_barcodes_with_high_umi_counts",
				"antibody_reads_in_cells": "frac_reads_in_cells",
				"fraction_antibody_reads": "frac_reads_recognized_antibody",
				"fraction_unrecognized_antibody": "frac_reads_unrecognized_antibody",
				"fraction_antibody_reads_usable": "frac_reads_usable",
				"sequencing_saturation": "frac_sequencing_saturation",
				"valid_barcodes": "frac_valid_barcodes",
				"antibody_reads_usable_per_cell": "mean_antibody_reads_usable_per_cell",
				"median_umis_per_cell_(summed_over_all_recognized_antibody_barcodes)": "median_umi_counts_per_cell",
				"number_of_reads": "total_reads"
			},
			"should_match": {},
			"perc_to_frac": [
				"frac_q30_bases_in_antibody_read",
				"frac_q30_bases_in_barcode",
				"frac_q30_bases_in_umi",
				"frac_reads_in_aggregate_barcodes",
				"frac_reads_in_cells",
				"frac_reads_recognized_antibody",
				"frac_reads_unrecognized_antibody",
				"frac_reads_usable",
				"frac_sequencing_saturation",
				"frac_valid_barcodes"
			]
		},
		"multiome": {
			"schema_mapping": {
				"feature_linkages_detected": "total_feature_linkages_detected",
				"linked_genes": "total_genes_linked",
				"linked_peaks": "total_peaks_linked",
				"estimated_number_of_cells": "total_cells_detected"
			}
		},
		"spatial": {
			"schema_mapping": {
				"fraction_reads_in_spots_under_tissue": "frac_reads_in_spots_under_tissue",
				"fraction_of_spots_under_tissue": "frac_spots_under_tissue",
				"number_of_spots_under_tissue": "total_spots_under_tissue"
			}
		}
	},
	"dragen": {
		"rna": {
			"schema_mapping": {
				"Passing cells": "total_cells_detected",
				"Mean reads per cell": "mean_reads_per_cell",
				"Median genes per cell": "median_genes_per_cell",
				"Median UMIs per cell": "median_umi_counts_per_cell",
				"Total input reads": "total_reads",
				"Q30 bases R1 %": "frac_q30_bases_in_rna_read",
				"Q30 bases R2 %": "frac_q30_bases_in_rna_read2",
				"Mapped reads %": "frac_reads_mapped_to_genome",
				"Sequencing saturation": "frac_sequencing_saturation",
				"Total genes detected": "total_genes_detected"
			},
			"perc_to_frac": [
				"frac_q30_bases_in_rna_read",
				"frac_q30_bases_in_rna_read2",
				"frac_reads_mapped_to_genome"
			]
		}
	},
	"star": {
		"rna": {
			"schema_mapping": {
				"q30_bases_in_cb+umi": "frac_q30_bases_in_barcode_or_umi",
				"q30_bases_in_rna_read": "frac_q30_bases_in_rna_read",
				"fraction_of_unique_reads_in_cells": "frac_reads_in_cells",
				"reads_mapped_to_genome:_unique": "frac_reads_mapped_confidently_to_genome",
				"reads_mapped_to_gene:_unique_gene": "frac_reads_mapped_confidently_to_transcriptome",
				"reads_mapped_to_genome:_unique+multiple": "frac_reads_mapped_to_genome",
				"reads_mapped_to_gene:_unique+multipe_gene": "frac_reads_mapped_to_transcriptome",
				"sequencing_saturation": "frac_sequencing_saturation",
				"reads_with_valid_barcodes": "frac_valid_barcodes",
				"mean_gene_per_cell": "mean_genes_per_cell",
				"mean_umi_per_cell": "mean_umi_counts_per_cell",
				"median_gene_per_cell": "median_genes_per_cell",
				"median_umi_per_cell": "median_umi_counts_per_cell",
				"estimated_number_of_cells": "total_cells_detected",
				"total_gene_detected": "total_genes_detected",
				"number_of_reads": "total_reads",
				"unique_reads_in_cells_mapped_to_gene": "total_reads_mapped_confidently_to_transcriptome",
				"umis_in_cells": "total_umis_in_cells"
			},
			"perc_to_frac": []
		}
	}
}
