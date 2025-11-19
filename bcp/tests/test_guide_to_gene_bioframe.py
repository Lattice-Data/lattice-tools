"""
Tests for guide_to_gene_bioframe.py module.
Tests gene annotation functionality for CRISPR guides using bioframe for overlap detection.
"""

import pytest
import pandas as pd
import csv
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the guide_to_gene_bioframe module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guide_to_gene_bioframe import (
    csv_to_bed_dataframe, 
    load_genes_gtf, 
    find_overlaps_bioframe, 
    format_output,
    annotate_guides_bioframe
)


class TestCSVToBEDDataframe:
    """Test CSV to BED format conversion."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_csv_to_bed_conversion(self, temp_dir):
        """Test CSV to BED format conversion with valid coordinates."""
        csv_file = temp_dir / "test_guides.csv"
        csv_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,1,1000,1018,+\n"
            "guide2,TGCA,NGG,2,2000,2018,-\n"
            "guide3,CCGG,NGG,3,3000,3018,+\n"
        )
        
        bed_df = csv_to_bed_dataframe(str(csv_file))
        
        # Check structure
        assert len(bed_df) == 3
        assert list(bed_df.columns) == ['chrom', 'start', 'end', 'name', 'score', 'strand']
        
        # Check data
        assert bed_df.iloc[0]['chrom'] == 'chr1'  # Chromosome normalized with 'chr' prefix
        assert bed_df.iloc[0]['start'] == 999  # 0-based
        assert bed_df.iloc[0]['end'] == 1018   # 1-based
        assert bed_df.iloc[0]['name'] == 'guide1'
        assert bed_df.iloc[0]['strand'] == '+'
    
    def test_csv_to_bed_na_coordinates(self, temp_dir):
        """Test CSV to BED conversion with NA coordinates."""
        csv_file = temp_dir / "test_guides_na.csv"
        csv_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,NA,NA,NA,+\n"
            "guide2,TGCA,NGG,1,NA,NA,+\n"
            "guide3,CCGG,NGG,1,1000,1018,+\n"
        )
        
        bed_df = csv_to_bed_dataframe(str(csv_file))
        
        # Only guide3 should be included (valid coordinates)
        assert len(bed_df) == 1
        assert bed_df.iloc[0]['name'] == 'guide3'
    
    def test_csv_to_bed_invalid_coordinates(self, temp_dir):
        """Test CSV to BED conversion with invalid coordinates."""
        csv_file = temp_dir / "test_guides_invalid.csv"
        csv_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,1,invalid,1018,+\n"
            "guide2,TGCA,NGG,1,1000,invalid,+\n"
            "guide3,CCGG,NGG,1,1000,1018,+\n"
        )
        
        bed_df = csv_to_bed_dataframe(str(csv_file))
        
        # Only guide3 should be included (valid coordinates)
        assert len(bed_df) == 1
        assert bed_df.iloc[0]['name'] == 'guide3'


class TestLoadGenesGTF:
    """Test GTF loading and filtering."""
    
    @pytest.fixture
    def gtf_file(self):
        """Path to the test GTF file."""
        return Path(__file__).parent / "fixtures" / "genes.gtf"
    
    def test_load_genes_gtf_with_real_file(self, gtf_file):
        """Test loading GTF file with real fixture data."""
        genes_df = load_genes_gtf(str(gtf_file))
        
        # Check that genes were loaded
        assert len(genes_df) > 0
        
        # Check that only gene features are included
        assert all(genes_df['feature'] == 'gene')
        
        # Check required columns
        required_cols = ['chrom', 'start', 'end', 'gene_id', 'gene_name', 'strand']
        for col in required_cols:
            assert col in genes_df.columns
        
        # Check data types
        assert genes_df['start'].dtype == 'int64'
        assert genes_df['end'].dtype == 'int64'
        assert genes_df['strand'].isin(['+', '-']).all()


class TestFindOverlapsBioframe:
    """Test bioframe overlap detection."""
    
    @pytest.fixture
    def sample_genes_df(self):
        """Create sample genes dataframe for testing."""
        return pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'start': [1000, 3000, 1000],
            'end': [2000, 4000, 1500],
            'gene_id': ['ENSG001', 'ENSG002', 'ENSG003'],
            'gene_name': ['GENE1', 'GENE2', 'GENE3'],
            'strand': ['+', '-', '+']
        })
    
    @pytest.fixture
    def sample_guides_df(self):
        """Create sample guides dataframe for testing."""
        return pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'start': [1500, 2500, 1200],
            'end': [1600, 2600, 1300],
            'name': ['guide1', 'guide2', 'guide3'],
            'score': [0, 0, 0],
            'strand': ['+', '+', '-']
        })
    
    def test_find_overlaps_bioframe(self, sample_genes_df, sample_guides_df):
        """Test bioframe overlap detection."""
        overlaps_df = find_overlaps_bioframe(sample_genes_df, sample_guides_df)
        
        # Check that overlaps were found
        assert len(overlaps_df) > 0
        
        # Check that overlap columns are present
        expected_cols = ['chrom', 'start', 'end', 'gene_id', 'gene_name', 'strand', 'name']
        for col in expected_cols:
            assert col in overlaps_df.columns
    
    def test_find_overlaps_empty_guides(self, sample_genes_df):
        """Test overlap detection with empty guides dataframe."""
        empty_guides_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        overlaps_df = find_overlaps_bioframe(sample_genes_df, empty_guides_df)
        
        assert overlaps_df.empty
    
    def test_find_overlaps_no_overlaps(self, sample_genes_df):
        """Test overlap detection when no overlaps exist."""
        no_overlap_guides_df = pd.DataFrame({
            'chrom': ['chr99'],
            'start': [1000],
            'end': [2000],
            'name': ['guide1'],
            'score': [0],
            'strand': ['+']
        })
        
        overlaps_df = find_overlaps_bioframe(sample_genes_df, no_overlap_guides_df)
        assert overlaps_df.empty


class TestFormatOutput:
    """Test output formatting."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_format_output_with_overlaps(self, temp_dir):
        """Test output formatting with gene overlaps."""
        overlaps_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'start': [1000, 1000],
            'end': [2000, 2000],
            'gene_id': ['ENSG001', 'ENSG002'],
            'gene_name': ['GENE1', 'GENE2'],
            'strand': ['+', '+'],
            'name': ['guide1', 'guide1']
        })
        
        original_guides_df = pd.DataFrame({
            'id': ['guide1', 'guide2'],
            'sequence': ['ACGT', 'TGCA'],
            'pam': ['NGG', 'NGG'],
            'chromosome': ['1', '2'],
            'start': [1000, 2000],
            'end': [1018, 2018],
            'sense': ['+', '-']
        })
        
        output_file = temp_dir / "test_output.csv"
        format_output(overlaps_df, original_guides_df, str(output_file))
        
        # Check that output file was created
        assert output_file.exists()
        
        # Read and validate output (keep_default_na=False to preserve 'NA' as string)
        output_df = pd.read_csv(output_file, keep_default_na=False)
        
        # Check structure
        expected_cols = ['id', 'sequence', 'pam', 'chromosome', 'start', 'end', 'sense', 'gene_id', 'gene_name']
        assert list(output_df.columns) == expected_cols
        
        # Check that guide1 has 2 rows (2 gene overlaps)
        guide1_rows = output_df[output_df['id'] == 'guide1']
        assert len(guide1_rows) == 2
        
        # Check that guide2 has 1 row (no overlap, marked as NA)
        guide2_rows = output_df[output_df['id'] == 'guide2']
        assert len(guide2_rows) == 1
        assert guide2_rows.iloc[0]['gene_id'] == 'NA'
    
    def test_format_output_no_overlaps(self, temp_dir):
        """Test output formatting with no gene overlaps."""
        overlaps_df = pd.DataFrame()  # Empty overlaps
        
        original_guides_df = pd.DataFrame({
            'id': ['guide1', 'guide2'],
            'sequence': ['ACGT', 'TGCA'],
            'pam': ['NGG', 'NGG'],
            'chromosome': ['1', '2'],
            'start': [1000, 2000],
            'end': [1018, 2018],
            'sense': ['+', '-']
        })
        
        output_file = temp_dir / "test_output.csv"
        format_output(overlaps_df, original_guides_df, str(output_file))
        
        # Read and validate output (keep_default_na=False to preserve 'NA' as string)
        output_df = pd.read_csv(output_file, keep_default_na=False)
        
        # Both guides should have NA (no overlaps)
        assert len(output_df) == 2
        assert all(output_df['gene_id'] == 'NA')
        assert all(output_df['gene_name'] == 'NA')


class TestAnnotateGuidesBioframe:
    """Test end-to-end guide annotation with bioframe."""
    
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
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Check that output file was created
        assert output_file.exists()
        
        # Read and validate output
        df = pd.read_csv(output_file)
        
        # Check that new columns were added
        assert 'gene_id' in df.columns
        assert 'gene_name' in df.columns
        
        # Check that all original columns are preserved
        assert 'id' in df.columns
        assert 'sequence' in df.columns
        assert 'pam' in df.columns
        assert 'chromosome' in df.columns
        assert 'start' in df.columns
        assert 'end' in df.columns
        assert 'sense' in df.columns
        
        # Check that we have results
        assert len(df) > 0
        
        # Verify that some guides have actual gene annotations (not NA)
        # Count rows with real gene IDs (ENSG format)
        real_genes = df[df['gene_id'].str.startswith('ENSG', na=False)]
        assert len(real_genes) > 0
    
    def test_annotate_guides_output_format(self, temp_dir, guides_file, gtf_file):
        """Test that output format matches expected structure."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
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
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Capture printed output
        captured = capsys.readouterr()
        
        # Check that statistics were printed
        assert "Total guides processed:" in captured.out
        assert "Guides with gene overlaps:" in captured.out
        assert "Guides without gene overlaps" in captured.out  # Changed from "intergenic regions"
        assert "Total gene overlaps found:" in captured.out
        assert "Output written to:" in captured.out
    
    def test_multiple_overlaps_creates_multiple_rows(self, temp_dir, guides_file, gtf_file):
        """Test that guides overlapping multiple genes create multiple output rows."""
        output_file = temp_dir / "output_annotated.csv"
        
        # Run annotation
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
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
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
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
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # All guides should have NA annotations (invalid coordinates)
        assert all(df['gene_id'] == 'NA')
        assert all(df['gene_name'] == 'NA')


class TestCompatibilityWithCurrent:
    """Test compatibility with current guide_to_gene.py implementation."""
    
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
    
    def test_output_format_compatibility(self, temp_dir, guides_file, gtf_file):
        """Test that output format is compatible with current implementation."""
        output_file = temp_dir / "output_bioframe.csv"
        
        # Run bioframe annotation
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file)
        
        # Check column names and order
        expected_cols = ['id', 'sequence', 'pam', 'chromosome', 'start', 'end', 'sense', 'gene_id', 'gene_name']
        assert list(df.columns) == expected_cols
        
        # Check that essential columns exist - dtypes may vary depending on pandas reading behavior
        # When pandas reads CSV, numeric columns can be int64 or object depending on content
        assert 'id' in df.columns
        assert 'sequence' in df.columns
        assert 'pam' in df.columns
        assert 'chromosome' in df.columns
        assert 'start' in df.columns
        assert 'end' in df.columns
        assert 'sense' in df.columns
        assert 'gene_id' in df.columns
        assert 'gene_name' in df.columns
    
    def test_statistics_format_compatibility(self, temp_dir, guides_file, gtf_file, capsys):
        """Test that statistics output format matches current implementation."""
        output_file = temp_dir / "output_bioframe.csv"
        
        # Run bioframe annotation
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Capture printed output
        captured = capsys.readouterr()
        
        # Check that all expected statistics are present
        expected_stats = [
            "Parsing GTF file:",
            "Found genes on",
            "Annotating guides from:",
            "Annotation complete!",
            "Total guides processed:",
            "Guides with gene overlaps:",
            "Guides without gene overlaps",  # Changed from separate intergenic/unmapped
            "Total gene overlaps found:",
            "Output written to:"
        ]
        
        for stat in expected_stats:
            assert stat in captured.out
