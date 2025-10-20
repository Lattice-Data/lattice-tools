#!/usr/bin/env python3
"""
Annotate CRISPR guides with gene identifiers from GTF file based on coordinate overlap.
Hybrid implementation: one row per overlap, includes gene names, handles edge cases.
"""

import csv
import argparse
from collections import defaultdict


def parse_gtf(gtf_file):
    """
    Parse GTF file and extract gene coordinates.
    Returns a dictionary: {chromosome: [(start, end, gene_id, gene_name, strand), ...]}
    """
    genes_by_chr = defaultdict(list)
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Only process gene features
            if feature == 'gene':
                # Extract gene_id and gene_name from attributes
                gene_id = None
                gene_name = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                
                if gene_id:
                    genes_by_chr[chrom].append((start, end, gene_id, gene_name or 'NA', strand))
    
    # Sort genes by start position for efficient searching
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort()
    
    return genes_by_chr


def normalize_chromosome(chrom):
    """
    Normalize chromosome format by adding 'chr' prefix if not present.
    """
    if chrom and chrom != 'NA' and not chrom.startswith('chr'):
        return f'chr{chrom}'
    return chrom


def find_overlapping_genes(chrom, start, end, genes_by_chr):
    """
    Find all genes that overlap with the given coordinates.
    Returns a list of (gene_id, gene_name, strand) tuples.
    """
    # Normalize chromosome format
    chrom = normalize_chromosome(chrom)
    
    if chrom not in genes_by_chr:
        return []
    
    overlapping = []
    for gene_start, gene_end, gene_id, gene_name, gene_strand in genes_by_chr[chrom]:
        # Check for overlap: guide overlaps gene if guide_end >= gene_start AND guide_start <= gene_end
        if end >= gene_start and start <= gene_end:
            overlapping.append((gene_id, gene_name, gene_strand))
    
    return overlapping


def annotate_guides(guide_file, gtf_file, output_file):
    """
    Annotate guides with gene identifiers from GTF file.
    One row per overlap - if multiple genes overlap, creates multiple rows.
    """
    print(f"Parsing GTF file: {gtf_file}")
    genes_by_chr = parse_gtf(gtf_file)
    print(f"Found genes on {len(genes_by_chr)} chromosomes")
    
    print(f"\nAnnotating guides from: {guide_file}")
    
    with open(guide_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        
        # Add new columns for gene annotation
        fieldnames = list(reader.fieldnames) + ['gene_id', 'gene_name']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        
        guides_processed = 0
        guides_with_genes = 0
        guides_intergenic = 0
        guides_unmapped = 0
        total_overlaps = 0
        
        for row in reader:
            guides_processed += 1
            guide_id = row['id']
            chrom = row['chromosome']
            
            # Handle NA or missing coordinates
            if chrom == 'NA' or row.get('start') == 'NA' or row.get('end') == 'NA':
                output_row = row.copy()
                output_row['gene_id'] = 'NA'
                output_row['gene_name'] = 'NA'
                writer.writerow(output_row)
                guides_unmapped += 1
            else:
                try:
                    start = int(row['start'])
                    end = int(row['end'])
                    
                    # Find overlapping genes
                    overlapping_genes = find_overlapping_genes(chrom, start, end, genes_by_chr)
                    
                    if overlapping_genes:
                        # Create one row per overlapping gene
                        guides_with_genes += 1
                        total_overlaps += len(overlapping_genes)
                        
                        for gene_id, gene_name, gene_strand in overlapping_genes:
                            output_row = row.copy()
                            output_row['gene_id'] = gene_id
                            output_row['gene_name'] = gene_name
                            writer.writerow(output_row)
                    else:
                        # No overlapping genes - intergenic region
                        output_row = row.copy()
                        output_row['gene_id'] = 'intergenic'
                        output_row['gene_name'] = 'intergenic'
                        writer.writerow(output_row)
                        guides_intergenic += 1
                        
                except ValueError:
                    # Handle invalid coordinate values
                    output_row = row.copy()
                    output_row['gene_id'] = 'INVALID_COORDS'
                    output_row['gene_name'] = 'INVALID_COORDS'
                    writer.writerow(output_row)
                    guides_unmapped += 1
    
    print(f"\nAnnotation complete!")
    print(f"  Total guides processed: {guides_processed}")
    print(f"  Guides with gene overlaps: {guides_with_genes}")
    print(f"  Guides in intergenic regions: {guides_intergenic}")
    print(f"  Guides with unmapped/invalid coordinates: {guides_unmapped}")
    print(f"  Total gene overlaps found: {total_overlaps}")
    if guides_with_genes > 0:
        print(f"  Average overlaps per guide (for guides with overlaps): {total_overlaps/guides_with_genes:.2f}")
    print(f"\nOutput written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Annotate CRISPR guides with gene identifiers from GTF file. '
                    'Creates one row per gene overlap.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py guides.csv gencode.gtf results.csv
  
Output format:
  - One row per guide-gene overlap
  - If a guide overlaps multiple genes, multiple rows are created
  - Includes gene_id (ENSG ID), gene_name, and strand information
        """
    )
    parser.add_argument('guide_file', help='Input CSV file with guide mappings')
    parser.add_argument('gtf_file', help='GTF file with gene annotations')
    parser.add_argument('output_file', help='Output CSV file with gene annotations')
    
    args = parser.parse_args()
    
    annotate_guides(args.guide_file, args.gtf_file, args.output_file)


if __name__ == '__main__':
    main()