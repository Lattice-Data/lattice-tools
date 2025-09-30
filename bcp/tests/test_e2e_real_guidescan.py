"""
End-to-end tests using actual guidescan binary execution.

These tests are marked with @pytest.mark.e2e and are skipped by default.
Run with: pytest -m e2e
"""

import pytest
import pandas as pd
import subprocess
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


def is_guidescan_available():
    """Check if guidescan is available in PATH."""
    try:
        result = subprocess.run(
            ["guidescan", "--help"],
            capture_output=True,
            check=True
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


# Skip all E2E tests if guidescan not available
pytestmark = [
    pytest.mark.e2e,
    pytest.mark.skipif(
        not is_guidescan_available(),
        reason="guidescan not available in PATH"
    )
]


@pytest.fixture
def toy_genome_index():
    """Get path to toy genome index for E2E tests."""
    fixtures_dir = Path(__file__).parent / "fixtures"
    index_path = fixtures_dir / "toy.fa.index"
    
    # Verify all index files exist
    assert index_path.with_suffix(".index.forward").exists()
    assert index_path.with_suffix(".index.reverse").exists()
    assert index_path.with_suffix(".index.gs").exists()
    
    return index_path


@pytest.fixture
def test_guides_file():
    """Get path to test guides file."""
    return Path(__file__).parent / "fixtures" / "test_guides.txt"


@pytest.fixture
def expected_output():
    """Get path to expected guidescan output."""
    return Path(__file__).parent / "fixtures" / "expected_all_matches.csv"


class TestGuidescanAvailability:
    """Test that guidescan is properly installed and accessible."""
    
    def test_guidescan_executable_exists(self):
        """Verify guidescan command is available."""
        result = subprocess.run(
            ["guidescan", "--help"],
            capture_output=True
        )
        assert result.returncode == 0
        assert b"guidescan" in result.stdout or b"guidescan" in result.stderr
    
    def test_guidescan_enumerate_command_exists(self):
        """Verify guidescan enumerate subcommand is available."""
        result = subprocess.run(
            ["guidescan", "enumerate", "--help"],
            capture_output=True
        )
        assert result.returncode == 0


class TestDirectGuidescanExecution:
    """Test running guidescan directly without the pipeline."""
    
    def test_direct_guidescan_command(self, toy_genome_index, test_guides_file, temp_dir):
        """Test executing guidescan command directly."""
        output_file = temp_dir / "direct_output.csv"
        
        cmd = [
            "guidescan", "enumerate",
            str(toy_genome_index),
            "-f", str(test_guides_file),
            "-o", str(output_file),
            "--format", "csv",
            "--mode", "complete",
            "-m", "1",
            "--max-off-targets", "1"
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        # Should complete successfully
        assert result.returncode == 0
        
        # Output file should exist
        assert output_file.exists()
        
        # Should be readable as CSV
        df = pd.read_csv(output_file)
        assert not df.empty
        assert 'id' in df.columns
        assert 'sequence' in df.columns
    
    def test_guidescan_output_format(self, toy_genome_index, test_guides_file, temp_dir):
        """Verify guidescan output has expected column structure."""
        output_file = temp_dir / "format_test.csv"
        
        cmd = [
            "guidescan", "enumerate",
            str(toy_genome_index),
            "-f", str(test_guides_file),
            "-o", str(output_file),
            "--format", "csv",
            "--mode", "complete",
            "-m", "1",
            "--max-off-targets", "1"
        ]
        
        subprocess.run(cmd, capture_output=True, check=True)
        
        df = pd.read_csv(output_file)
        
        # Check expected columns
        expected_columns = [
            'id', 'sequence', 'match_chrm', 'match_position',
            'match_strand', 'match_distance', 'match_sequence',
            'rna_bulges', 'dna_bulges', 'specificity'
        ]
        
        for col in expected_columns:
            assert col in df.columns, f"Missing column: {col}"


class TestFullPipelineWithRealGuidescan:
    """Test the complete pipeline using real guidescan execution."""
    
    def test_full_pipeline_execution(self, toy_genome_index, test_guides_file, temp_dir):
        """Test complete pipeline with actual guidescan on toy genome."""
        output_dir = temp_dir / "pipeline_output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(toy_genome_index),
            guides_file=str(test_guides_file),
            output_dir=str(output_dir),
            keep_intermediate=True,
            verbose=True
        )
        
        # Run the pipeline
        result_df = pipeline.run()
        
        # Verify pipeline completed
        assert result_df is not None
        assert isinstance(result_df, pd.DataFrame)
        
        # Check all expected output files exist
        assert pipeline.all_matches_file.exists()
        assert pipeline.best_matches_file.exists()
        assert pipeline.final_output_file.exists()
        
        # Verify output structure
        assert 'id' in result_df.columns
        assert 'sequence' in result_df.columns
        assert 'pam' in result_df.columns
        assert 'chromosome' in result_df.columns
        assert 'start' in result_df.columns
        assert 'end' in result_df.columns
        assert 'sense' in result_df.columns
        
        # Should have results for all input guides
        assert len(result_df) == 5, "Should have 5 guides"
        
        # At least one guide should have a match (guide5)
        matches = result_df[result_df['chromosome'] != 'NA']
        assert len(matches) >= 1, "At least one guide should have genomic coordinates"
    
    def test_pipeline_output_properties(self, toy_genome_index, test_guides_file, temp_dir):
        """Test that pipeline output has valid properties."""
        output_dir = temp_dir / "properties_test"
        
        pipeline = GuidescanPipeline(
            genome_index=str(toy_genome_index),
            guides_file=str(test_guides_file),
            output_dir=str(output_dir),
            keep_intermediate=False
        )
        
        result_df = pipeline.run()
        
        # Property checks
        for _, row in result_df.iterrows():
            # If has genomic match, coordinates should be valid
            if row['chromosome'] != 'NA':
                assert pd.notna(row['start'])
                assert pd.notna(row['end'])
                assert row['start'] <= row['end'], "Start should be <= end"
                assert row['sense'] in ['+', '-'], "Sense should be + or -"
            else:
                # If no match, coordinates should be NA
                assert row['start'] == 'NA'
                assert row['end'] == 'NA'
                assert row['sense'] == 'NA'
            
            # Sequence and PAM should always be present
            assert pd.notna(row['sequence'])
            assert pd.notna(row['pam'])
            assert len(row['pam']) > 0
    
    def test_pipeline_with_cleanup(self, toy_genome_index, test_guides_file, temp_dir):
        """Test that cleanup removes intermediate files correctly."""
        output_dir = temp_dir / "cleanup_test"
        
        pipeline = GuidescanPipeline(
            genome_index=str(toy_genome_index),
            guides_file=str(test_guides_file),
            output_dir=str(output_dir),
            keep_intermediate=False  # Should cleanup
        )
        
        result_df = pipeline.run()
        
        # Final output should exist
        assert pipeline.final_output_file.exists()
        
        # Intermediate files should be removed
        assert not pipeline.all_matches_file.exists()
        assert not pipeline.best_matches_file.exists()


class TestE2EEdgeCases:
    """Test edge cases with real guidescan execution."""
    
    def test_guides_with_matches_and_no_matches(self, toy_genome_index, test_guides_file, temp_dir):
        """Test handling of both guides that match and don't match the toy genome."""
        output_dir = temp_dir / "mixed_matches_test"
        
        pipeline = GuidescanPipeline(
            genome_index=str(toy_genome_index),
            guides_file=str(test_guides_file),
            output_dir=str(output_dir)
        )
        
        result_df = pipeline.run()
        
        # Check for guides with no matches
        no_match_guides = result_df[result_df['chromosome'] == 'NA']
        
        # Most guides in test set should have no matches in toy genome
        assert len(no_match_guides) >= 3, "Most guides should have no matches"
        
        # Verify NA guides have proper structure
        for _, row in no_match_guides.iterrows():
            assert row['chromosome'] == 'NA'
            assert row['start'] == 'NA'
            assert row['end'] == 'NA'
            assert row['sense'] == 'NA'
            assert pd.notna(row['sequence'])
            assert pd.notna(row['pam'])
        
        # Check for guides with matches
        match_guides = result_df[result_df['chromosome'] != 'NA']
        
        # At least one guide should have a match (guide5)
        assert len(match_guides) >= 1, "At least one guide should match"
        
        # Verify matched guides have valid coordinates
        for _, row in match_guides.iterrows():
            assert row['chromosome'] != 'NA'
            assert pd.notna(row['start'])
            assert pd.notna(row['end'])
            assert row['start'] <= row['end']
            assert row['sense'] in ['+', '-']

