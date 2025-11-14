"""
Compatibility tests comparing bioframe implementation vs current guide_to_gene.py.
Tests that both implementations produce identical results.
"""

import pytest
import pandas as pd
import tempfile
from pathlib import Path
import sys

# Add parent directory to path to import both modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from guide_to_gene import annotate_guides
from guide_to_gene_bioframe import annotate_guides_bioframe


class TestCompatibilityBioframeVsCurrent:
    """Test compatibility between bioframe and current implementations."""
    
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
    
    def test_identical_output_format(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations produce identical output format."""
        # Run current implementation
        current_output = temp_dir / "current_output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Check that column names are identical
        assert list(current_df.columns) == list(bioframe_df.columns)
        
        # Check that data types are compatible
        for col in current_df.columns:
            assert current_df[col].dtype == bioframe_df[col].dtype
    
    def test_same_number_of_output_rows(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations produce the same number of output rows."""
        # Run current implementation
        current_output = temp_dir / "current_output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Check that both produce the same number of rows
        assert len(current_df) == len(bioframe_df)
    
    def test_same_guide_coverage(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations process the same guides."""
        # Run current implementation
        current_output = temp_dir / "current_output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Check that both process the same unique guides
        current_guides = set(current_df['id'].unique())
        bioframe_guides = set(bioframe_df['id'].unique())
        
        assert current_guides == bioframe_guides
    
    def test_same_gene_annotations(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations produce the same gene annotations."""
        # Run current implementation
        current_output = temp_dir / "current_output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Sort both dataframes by id and gene_id for comparison
        current_sorted = current_df.sort_values(['id', 'gene_id']).reset_index(drop=True)
        bioframe_sorted = bioframe_df.sort_values(['id', 'gene_id']).reset_index(drop=True)
        
        # Check that gene annotations are identical
        assert current_sorted['gene_id'].equals(bioframe_sorted['gene_id'])
        assert current_sorted['gene_name'].equals(bioframe_sorted['gene_name'])
    
    def test_same_statistics(self, temp_dir, guides_file, gtf_file, capsys):
        """Test that both implementations produce similar statistics."""
        # Run current implementation
        current_output = temp_dir / "current_output.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        current_captured = capsys.readouterr()
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        bioframe_captured = capsys.readouterr()
        
        # Check that both produce statistics
        assert "Total guides processed:" in current_captured.out
        assert "Total guides processed:" in bioframe_captured.out
        
        assert "Guides with gene overlaps:" in current_captured.out
        assert "Guides with gene overlaps:" in bioframe_captured.out
        
        assert "Guides in intergenic regions:" in current_captured.out
        assert "Guides in intergenic regions:" in bioframe_captured.out
    
    def test_edge_cases_compatibility(self, temp_dir, gtf_file):
        """Test that both implementations handle edge cases identically."""
        # Create test file with edge cases
        guides_file = temp_dir / "edge_cases.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGT,NGG,NA,NA,NA,+\n"
            "guide2,TGCA,NGG,1,NA,NA,+\n"
            "guide3,CCGG,NGG,1,invalid,1000,+\n"
            "guide4,TTAA,NGG,1,1000,invalid,+\n"
            "guide5,GGCC,NGG,1,1000,1018,+\n"
        )
        
        # Run current implementation
        current_output = temp_dir / "current_edge_cases.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_edge_cases.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Check that both handle edge cases identically
        assert len(current_df) == len(bioframe_df)
        
        # Check that both produce the same gene annotations for edge cases
        current_sorted = current_df.sort_values('id').reset_index(drop=True)
        bioframe_sorted = bioframe_df.sort_values('id').reset_index(drop=True)
        
        assert current_sorted['gene_id'].equals(bioframe_sorted['gene_id'])
        assert current_sorted['gene_name'].equals(bioframe_sorted['gene_name'])
    
    def test_multiple_overlaps_compatibility(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations handle multiple overlaps identically."""
        # Run current implementation
        current_output = temp_dir / "current_multiple.csv"
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        
        # Run bioframe implementation
        bioframe_output = temp_dir / "bioframe_multiple.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        
        # Read both outputs
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        # Check that both handle multiple overlaps the same way
        # Group by guide ID and count overlaps
        current_overlaps = current_df.groupby('id').size()
        bioframe_overlaps = bioframe_df.groupby('id').size()
        
        # Check that each guide has the same number of overlaps in both implementations
        for guide_id in current_overlaps.index:
            assert current_overlaps[guide_id] == bioframe_overlaps[guide_id]
    
    def test_performance_comparison(self, temp_dir, guides_file, gtf_file):
        """Test that both implementations complete successfully (basic performance check)."""
        import time
        
        # Time current implementation
        current_output = temp_dir / "current_perf.csv"
        start_time = time.time()
        annotate_guides(str(guides_file), str(gtf_file), str(current_output))
        current_time = time.time() - start_time
        
        # Time bioframe implementation
        bioframe_output = temp_dir / "bioframe_perf.csv"
        start_time = time.time()
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(bioframe_output))
        bioframe_time = time.time() - start_time
        
        # Both should complete successfully
        assert current_output.exists()
        assert bioframe_output.exists()
        
        # Both should produce results
        current_df = pd.read_csv(current_output)
        bioframe_df = pd.read_csv(bioframe_output)
        
        assert len(current_df) > 0
        assert len(bioframe_df) > 0
        
        # Log performance comparison (for debugging)
        print(f"Current implementation time: {current_time:.2f}s")
        print(f"Bioframe implementation time: {bioframe_time:.2f}s")
        print(f"Speedup: {current_time/bioframe_time:.2f}x")
