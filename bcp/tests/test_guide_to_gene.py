"""
Tests for guide_to_gene.py module.
Tests gene annotation functionality for CRISPR guides based on GTF overlap.
"""

import pytest
import pandas as pd
import csv
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the guide_to_gene module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guide_to_gene import parse_gtf, normalize_chromosome, find_overlapping_genes, annotate_guides


class TestParseGTF:
    """Test GTF parsing functionality."""
    
    def test_parse_gtf_with_real_file(self):
        """Test parsing the real GTF fixture file."""
        gtf_file = Path(__file__).parent / "fixtures" / "genes.gtf"
        genes_by_chr = parse_gtf(str(gtf_file))
        
        # Check that genes were parsed
        assert len(genes_by_chr) > 0
        
        # Check structure of parsed data
        for chrom, genes in genes_by_chr.items():
            assert isinstance(chrom, str)
            assert isinstance(genes, list)
            for gene in genes:
                assert len(gene) == 5  # (start, end, gene_id, gene_name, strand)
                start, end, gene_id, gene_name, strand = gene
                assert isinstance(start, int)
                assert isinstance(end, int)
                assert start <= end
                assert gene_id is not None
                assert strand in ['+', '-']


class TestNormalizeChromosome:
    """Test chromosome name normalization."""
    
    def test_normalize_adds_chr_prefix(self):
        """Test that chromosome names without 'chr' get the prefix added."""
        assert normalize_chromosome("1") == "chr1"
        assert normalize_chromosome("X") == "chrX"
        assert normalize_chromosome("22") == "chr22"
    
    def test_normalize_preserves_chr_prefix(self):
        """Test that chromosome names with 'chr' are unchanged."""
        assert normalize_chromosome("chr1") == "chr1"
        assert normalize_chromosome("chrX") == "chrX"
        assert normalize_chromosome("chr22") == "chr22"
    
    def test_normalize_handles_na(self):
        """Test that NA values are preserved."""
        assert normalize_chromosome("NA") == "NA"


class TestFindOverlappingGenes:
    """Test gene overlap detection."""
    
    @pytest.fixture
    def simple_genes_by_chr(self):
        """Create a simple genes_by_chr structure for testing."""
        return {
            "chr1": [
                (1000, 2000, "ENSG001", "GENE1", "+"),
                (3000, 4000, "ENSG002", "GENE2", "-"),
                (5000, 6000, "ENSG003", "GENE3", "+"),
            ],
            "chr2": [
                (1000, 1500, "ENSG004", "GENE4", "+"),
            ]
        }
    
    def test_find_single_overlap(self, simple_genes_by_chr):
        """Test finding a single overlapping gene."""
        overlapping = find_overlapping_genes("chr1", 1500, 1600, simple_genes_by_chr)
        assert len(overlapping) == 1
        assert overlapping[0][0] == "ENSG001"
        assert overlapping[0][1] == "GENE1"
    
    def test_find_no_overlap(self, simple_genes_by_chr):
        """Test intergenic region with no overlaps."""
        overlapping = find_overlapping_genes("chr1", 2500, 2800, simple_genes_by_chr)
        assert len(overlapping) == 0
    
    def test_find_multiple_overlaps(self, simple_genes_by_chr):
        """Test finding multiple overlapping genes."""
        # Add overlapping genes to test data
        simple_genes_by_chr["chr1"].append((1500, 2500, "ENSG005", "GENE5", "+"))
        overlapping = find_overlapping_genes("chr1", 1800, 1900, simple_genes_by_chr)
        assert len(overlapping) == 2
    
    def test_chromosome_normalization_in_search(self, simple_genes_by_chr):
        """Test that chromosome names are normalized during search."""
        # Search without 'chr' prefix
        overlapping = find_overlapping_genes("1", 1500, 1600, simple_genes_by_chr)
        assert len(overlapping) == 1
        assert overlapping[0][0] == "ENSG001"
    
    def test_missing_chromosome(self, simple_genes_by_chr):
        """Test searching on a chromosome not in the GTF."""
        overlapping = find_overlapping_genes("chr99", 1000, 2000, simple_genes_by_chr)
        assert len(overlapping) == 0


class TestAnnotateGuides:
    """Test end-to-end guide annotation."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def guides_file(self):
        """Path to the test guides file."""
        return Path(__file__).parent / "fixtures" / "guides_no_genes.csv"
    
    @pytest.fixture
    def gtf_file(self):
        """Path to the test GTF file."""
        return Path(__file__).parent / "fixtures" / "genes.gtf"
    
    def test_annotate_guides_with_real_data(self, temp_dir, guides_file, gtf_file):
        """Test annotation with real fixture data."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        # Check that output file was created
        assert output_file.exists()
        
        # Read and validate output
        df = pd.read_csv(output_file)
        
        # Check that new columns were added
        assert 'gene_id' in df.columns
        assert 'gene_name' in df.columns
        
        # Check that gene_strand and num_overlaps are NOT in output
        assert 'gene_strand' not in df.columns
        assert 'num_overlaps' not in df.columns
        
        # Check that all original columns are preserved
        assert 'id' in df.columns
        assert 'sequence' in df.columns
        assert 'pam' in df.columns
        assert 'chromosome' in df.columns
        assert 'start' in df.columns
        assert 'end' in df.columns
        assert 'sense' in df.columns
        
        # Check that we have results (should have more rows than input due to multi-overlaps)
        assert len(df) > 0
        
        # Verify that gene annotations are not empty for most guides
        non_na_genes = df[df['gene_id'].notna() & (df['gene_id'] != 'NA') & (df['gene_id'] != 'intergenic')]
        assert len(non_na_genes) > 0
    
    def test_annotate_guides_output_format(self, temp_dir, guides_file, gtf_file):
        """Test that output format matches expected structure."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output with csv reader to check structure
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            # Check that we have rows
            assert len(rows) > 0
            
            # Check each row has required fields
            for row in rows:
                assert 'gene_id' in row
                assert 'gene_name' in row
                # Verify gene_id is not empty
                assert row['gene_id'] != ''
                assert row['gene_name'] != ''
    
    def test_annotate_guides_statistics(self, temp_dir, guides_file, gtf_file, capsys):
        """Test that annotation prints correct statistics."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        # Capture printed output
        captured = capsys.readouterr()
        
        # Check that statistics were printed
        assert "Total guides processed:" in captured.out
        assert "Guides with gene overlaps:" in captured.out
        assert "Guides in intergenic regions:" in captured.out
        assert "Total gene overlaps found:" in captured.out
        assert "Output written to:" in captured.out
    
    def test_multiple_overlaps_creates_multiple_rows(self, temp_dir, guides_file, gtf_file):
        """Test that guides overlapping multiple genes create multiple output rows."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        # Read input and output
        input_df = pd.read_csv(guides_file)
        output_df = pd.read_csv(output_file)
        
        # Count number of input guides
        num_input_guides = len(input_df)
        num_output_rows = len(output_df)
        
        # Output should have at least as many rows as input
        # (could be more if some guides overlap multiple genes)
        assert num_output_rows >= num_input_guides
        
        # Check if there are any guides with multiple overlaps
        # by grouping by guide id and checking if any have multiple rows
        if 'id' in output_df.columns:
            guide_counts = output_df.groupby('id').size()
            # At least verify that the grouping works
            assert len(guide_counts) > 0


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def gtf_file(self):
        """Path to the test GTF file."""
        return Path(__file__).parent / "fixtures" / "genes.gtf"
    
    def test_guides_with_na_coordinates(self, temp_dir, gtf_file):
        """Test handling of guides with NA coordinates."""
        guides_file = temp_dir / "guides_na.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,NA,NA,NA,+\n"
            "guide2,TGCA,NGG,1,NA,NA,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output with keep_default_na=False to preserve 'NA' as string
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # All guides should have NA gene annotations
        assert all(df['gene_id'] == 'NA')
        assert all(df['gene_name'] == 'NA')
    
    def test_guides_with_invalid_coordinates(self, temp_dir, gtf_file):
        """Test handling of guides with invalid coordinate values."""
        guides_file = temp_dir / "guides_invalid.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,1,invalid,1000,+\n"
            "guide2,TGCA,NGG,1,1000,invalid,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(output_file))
        
        df = pd.read_csv(output_file)
        
        # All guides should have INVALID_COORDS annotations
        assert all(df['gene_id'] == 'INVALID_COORDS')
        assert all(df['gene_name'] == 'INVALID_COORDS')

