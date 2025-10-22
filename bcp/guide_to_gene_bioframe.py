#!/usr/bin/env python3
"""
Annotate CRISPR guides with gene identifiers from GTF file using bioframe for overlap detection.
Bioframe implementation: one row per overlap, includes gene names, handles edge cases.
"""

import csv
import argparse
import pandas as pd
import bioframe as bf
from pathlib import Path


def csv_to_bed_dataframe(csv_file):
    """
    Convert exact_matches_formatted.csv to BED format for bioframe.
    
    Args:
        csv_file: Path to CSV file with guide coordinates
        
    Returns:
        pandas.DataFrame: BED format dataframe with columns ['chrom', 'start', 'end', 'name', 'score', 'strand']
    """
    bed_data = []
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            guide_id = row['id']
            chrom = row['chromosome']
            start = row['start']
            end = row['end']
            strand = row['sense']
            
            # Handle NA or missing coordinates
            if chrom == 'NA' or start == 'NA' or end == 'NA':
                continue  # Skip guides with NA coordinates
            
            try:
                start_int = int(start)
                end_int = int(end)
                
                # Normalize chromosome format (add 'chr' prefix if not present)
                if chrom and chrom != 'NA' and not chrom.startswith('chr'):
                    chrom = f'chr{chrom}'
                
                # Convert to BED format (0-based start, 1-based end)
                bed_data.append({
                    'chrom': chrom,
                    'start': start_int - 1,  # Convert to 0-based
                    'end': end_int,          # Keep 1-based end
                    'name': guide_id,
                    'score': 0,
                    'strand': strand
                })
                
            except ValueError:
                # Skip guides with invalid coordinates
                continue
    
    return pd.DataFrame(bed_data)


def load_genes_gtf(gtf_file):
    """
    Load and filter GTF to genes only using bioframe.
    
    Args:
        gtf_file: Path to GTF file
        
    Returns:
        pandas.DataFrame: Filtered genes dataframe with gene_id and gene_name extracted
    """
    print(f"Parsing GTF file: {gtf_file}")
    
    # Load GTF file with bioframe
    genes_df = bf.read_table(gtf_file, schema='gtf')
    
    # Filter to gene features only
    genes_only = genes_df[genes_df['feature'] == 'gene'].copy()
    
    # Extract gene_id and gene_name from attributes column
    def extract_attribute(attr_string, key):
        """Extract a specific attribute from GTF attributes string."""
        import re
        pattern = f'{key} "([^"]+)"'
        match = re.search(pattern, attr_string)
        return match.group(1) if match else 'NA'
    
    genes_only['gene_id'] = genes_only['attributes'].apply(lambda x: extract_attribute(x, 'gene_id'))
    genes_only['gene_name'] = genes_only['attributes'].apply(lambda x: extract_attribute(x, 'gene_name'))
    
    # Convert start and end to int64 for bioframe compatibility
    genes_only['start'] = genes_only['start'].astype('int64')
    genes_only['end'] = genes_only['end'].astype('int64')
    
    print(f"Found genes on {len(genes_only['chrom'].unique())} chromosomes")
    
    return genes_only


def find_overlaps_bioframe(genes_df, guides_df):
    """
    Use bioframe.overlap to find gene-guide overlaps.
    
    Args:
        genes_df: DataFrame with gene annotations
        guides_df: DataFrame with guide coordinates in BED format
        
    Returns:
        pandas.DataFrame: Overlap results from bioframe with cleaned column names
    """
    if guides_df.empty:
        return pd.DataFrame()
    
    # Ensure proper data types for bioframe
    genes_df = genes_df.copy()
    guides_df = guides_df.copy()
    
    genes_df['start'] = genes_df['start'].astype('int64')
    genes_df['end'] = genes_df['end'].astype('int64')
    guides_df['start'] = guides_df['start'].astype('int64')
    guides_df['end'] = guides_df['end'].astype('int64')
    
    # Find overlaps using bioframe
    # bioframe.overlap adds suffixes to columns from the second dataframe
    overlaps = bf.overlap(genes_df, guides_df, how='inner', 
                         cols1=['chrom', 'start', 'end'],
                         cols2=['chrom', 'start', 'end'])
    
    # Rename columns with '_' suffix from the second dataframe (guides)
    # bioframe adds '_' suffix to overlapping column names
    if not overlaps.empty and 'name_' in overlaps.columns:
        overlaps = overlaps.rename(columns={'name_': 'name'})
    
    return overlaps


def format_output(overlaps_df, original_guides_df, output_file):
    """
    Format bioframe results to match current output format.
    
    Args:
        overlaps_df: DataFrame with overlap results from bioframe
        original_guides_df: DataFrame with original guide information
        output_file: Path to output CSV file
    """
    # Create a mapping from guide_id to original guide data
    guide_mapping = {}
    for _, row in original_guides_df.iterrows():
        guide_id = row['id']
        if guide_id not in guide_mapping:
            guide_mapping[guide_id] = row.to_dict()
    
    # Prepare output data
    output_data = []
    
    if not overlaps_df.empty:
        for _, overlap_row in overlaps_df.iterrows():
            guide_id = overlap_row['name']
            gene_id = overlap_row['gene_id']
            gene_name = overlap_row['gene_name']
            
            # Get original guide data
            if guide_id in guide_mapping:
                original_guide = guide_mapping[guide_id]
                
                output_row = {
                    'id': guide_id,
                    'sequence': original_guide['sequence'],
                    'pam': original_guide['pam'],
                    'chromosome': original_guide['chromosome'],
                    'start': original_guide['start'],
                    'end': original_guide['end'],
                    'sense': original_guide['sense'],
                    'gene_id': gene_id,
                    'gene_name': gene_name
                }
                output_data.append(output_row)
    
    # Handle guides with no overlaps (intergenic)
    processed_guide_ids = set(overlaps_df['name']) if not overlaps_df.empty else set()
    
    for guide_id, original_guide in guide_mapping.items():
        if guide_id not in processed_guide_ids:
            # Check if guide has valid coordinates
            chrom = str(original_guide['chromosome'])
            start = str(original_guide['start'])
            end = str(original_guide['end'])
            
            # Check for NA values (as string or pandas NA)
            is_na = (chrom == 'NA' or chrom == 'nan' or pd.isna(original_guide['chromosome']) or
                     start == 'NA' or start == 'nan' or pd.isna(original_guide['start']) or
                     end == 'NA' or end == 'nan' or pd.isna(original_guide['end']))
            
            if is_na:
                # Unmapped guide
                output_row = {
                    'id': guide_id,
                    'sequence': original_guide['sequence'],
                    'pam': original_guide['pam'],
                    'chromosome': 'NA',
                    'start': 'NA',
                    'end': 'NA',
                    'sense': original_guide['sense'],
                    'gene_id': 'NA',
                    'gene_name': 'NA'
                }
            else:
                try:
                    int(start)
                    int(end)
                    # Valid coordinates but no gene overlap - intergenic
                    output_row = {
                        'id': guide_id,
                        'sequence': original_guide['sequence'],
                        'pam': original_guide['pam'],
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                        'sense': original_guide['sense'],
                        'gene_id': 'intergenic',
                        'gene_name': 'intergenic'
                    }
                except ValueError:
                    # Invalid coordinates
                    output_row = {
                        'id': guide_id,
                        'sequence': original_guide['sequence'],
                        'pam': original_guide['pam'],
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                        'sense': original_guide['sense'],
                        'gene_id': 'INVALID_COORDS',
                        'gene_name': 'INVALID_COORDS'
                    }
            
            output_data.append(output_row)
    
    # Write output CSV
    if output_data:
        output_df = pd.DataFrame(output_data)
        # Ensure all columns are strings to match current implementation output
        for col in output_df.columns:
            output_df[col] = output_df[col].astype(str)
        output_df.to_csv(output_file, index=False)
    else:
        # Create empty CSV with headers
        fieldnames = ['id', 'sequence', 'pam', 'chromosome', 'start', 'end', 'sense', 'gene_id', 'gene_name']
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()


def annotate_guides_bioframe(guide_file, gtf_file, output_file):
    """
    Annotate guides with gene identifiers from GTF file using bioframe.
    One row per overlap - if multiple genes overlap, creates multiple rows.
    
    Args:
        guide_file: Path to CSV file with guide mappings
        gtf_file: Path to GTF file with gene annotations
        output_file: Path to output CSV file
    """
    print(f"\nAnnotating guides from: {guide_file}")
    
    # Load original guide data for reference
    original_guides_df = pd.read_csv(guide_file)
    
    # Convert guides to BED format
    guides_bed_df = csv_to_bed_dataframe(guide_file)
    
    # Load and filter genes
    genes_df = load_genes_gtf(gtf_file)
    
    # Find overlaps
    overlaps_df = find_overlaps_bioframe(genes_df, guides_bed_df)
    
    # Format and write output
    format_output(overlaps_df, original_guides_df, output_file)
    
    # Calculate and print statistics
    guides_processed = len(original_guides_df)
    guides_with_genes = len(overlaps_df['name'].unique()) if not overlaps_df.empty else 0
    total_overlaps = len(overlaps_df) if not overlaps_df.empty else 0
    
    # Count intergenic guides
    processed_guide_ids = set(overlaps_df['name']) if not overlaps_df.empty else set()
    guides_intergenic = 0
    guides_unmapped = 0
    
    for _, row in original_guides_df.iterrows():
        guide_id = row['id']
        if guide_id not in processed_guide_ids:
            chrom = str(row['chromosome'])
            start = str(row['start'])
            end = str(row['end'])
            
            # Check for NA values (as string or pandas NA)
            is_na = (chrom == 'NA' or chrom == 'nan' or pd.isna(row['chromosome']) or
                     start == 'NA' or start == 'nan' or pd.isna(row['start']) or
                     end == 'NA' or end == 'nan' or pd.isna(row['end']))
            
            if is_na:
                guides_unmapped += 1
            else:
                try:
                    int(start)
                    int(end)
                    guides_intergenic += 1
                except ValueError:
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
        description='Annotate CRISPR guides with gene identifiers from GTF file using bioframe. '
                    'Creates one row per gene overlap.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python guide_to_gene_bioframe.py guides.csv gencode.gtf results.csv
  
Output format:
  - One row per guide-gene overlap
  - If a guide overlaps multiple genes, multiple rows are created
  - Includes gene_id (ENSG ID), gene_name, and strand information
  - Uses bioframe for efficient overlap detection
        """
    )
    parser.add_argument('guide_file', help='Input CSV file with guide mappings')
    parser.add_argument('gtf_file', help='GTF file with gene annotations')
    parser.add_argument('output_file', help='Output CSV file with gene annotations')
    
    args = parser.parse_args()
    
    annotate_guides_bioframe(args.guide_file, args.gtf_file, args.output_file)


if __name__ == '__main__':
    main()
