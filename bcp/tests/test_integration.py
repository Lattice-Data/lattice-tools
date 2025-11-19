"""
Integration tests for the guidescan pipeline.
Tests pipeline orchestration, file management, and error handling.
"""

import pytest
import pandas as pd
from pathlib import Path
import tempfile
import sys
from unittest.mock import patch

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


class TestPipelineInitialization:
    """Test pipeline initialization and configuration."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_pipeline_initialization(self, temp_dir):
        """Test that pipeline initializes with correct parameters."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
            keep_intermediate=True,
            verbose=False
        )
        
        # Check attributes
        assert pipeline.genome_index == genome_index
        assert pipeline.guides_file == guides_file
        assert pipeline.output_dir == output_dir
        assert pipeline.keep_intermediate is True
        assert pipeline.max_mismatches == 1
        assert pipeline.max_off_targets == 1
        assert pipeline.pam_length is None  # Not set until validation
    
    def test_output_directory_created(self, temp_dir):
        """Test that output directory is created automatically."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output" / "nested" / "path"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        # Output directory should be created
        assert pipeline.output_dir.exists()
        assert pipeline.output_dir.is_dir()
    
    def test_output_file_paths_set(self, temp_dir):
        """Test that output file paths are correctly set."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        # Check file paths
        assert pipeline.all_matches_file == output_dir / "all_matches.csv"
        assert pipeline.best_matches_file == output_dir / "best_matches.csv"
        assert pipeline.final_output_file == output_dir / "exact_matches_formatted.csv"


class TestInputValidation:
    """Test input validation and error handling."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_missing_genome_index_raises_error(self, temp_dir):
        """Test that missing genome index file raises FileNotFoundError."""
        genome_index = temp_dir / "nonexistent.fa.index"
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        with pytest.raises(FileNotFoundError) as excinfo:
            pipeline.validate_inputs()
        
        assert "Genome index not found" in str(excinfo.value)
    
    def test_missing_guides_file_raises_error(self, temp_dir):
        """Test that missing guides file raises FileNotFoundError."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "nonexistent.txt"
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        with pytest.raises(FileNotFoundError) as excinfo:
            pipeline.validate_inputs()
        
        assert "Guides file not found" in str(excinfo.value)
    
    def test_validation_sets_pam_length(self, temp_dir):
        """Test that validation sets PAM length from input file."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        # Mock guidescan availability check
        with patch('subprocess.run'):
            pipeline.validate_inputs()
        
        assert pipeline.pam_length == 3


class TestCleanupFunctionality:
    """Test intermediate file cleanup functionality."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def pipeline_with_files(self, temp_dir):
        """Create pipeline with intermediate files."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
            keep_intermediate=False  # Should cleanup
        )
        
        # Create intermediate files
        pipeline.all_matches_file.parent.mkdir(parents=True, exist_ok=True)
        pipeline.all_matches_file.write_text("dummy,data\n")
        pipeline.best_matches_file.write_text("dummy,data\n")
        pipeline.final_output_file.write_text("dummy,data\n")
        
        return pipeline
    
    def test_cleanup_removes_intermediate_files(self, pipeline_with_files):
        """Test that cleanup removes intermediate files when keep_intermediate=False."""
        pipeline = pipeline_with_files
        
        # Verify files exist before cleanup
        assert pipeline.all_matches_file.exists()
        assert pipeline.best_matches_file.exists()
        assert pipeline.final_output_file.exists()
        
        # Run cleanup
        pipeline.cleanup_intermediate_files()
        
        # Intermediate files should be removed
        assert not pipeline.all_matches_file.exists()
        assert not pipeline.best_matches_file.exists()
        
        # Final output should remain
        assert pipeline.final_output_file.exists()
    
    def test_keep_intermediate_preserves_files(self, temp_dir):
        """Test that intermediate files are kept when keep_intermediate=True."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
            keep_intermediate=True  # Should keep files
        )
        
        # Create intermediate files
        pipeline.all_matches_file.parent.mkdir(parents=True, exist_ok=True)
        pipeline.all_matches_file.write_text("dummy,data\n")
        pipeline.best_matches_file.write_text("dummy,data\n")
        
        # Run cleanup
        pipeline.cleanup_intermediate_files()
        
        # All files should still exist
        assert pipeline.all_matches_file.exists()
        assert pipeline.best_matches_file.exists()


class TestPipelineWorkflow:
    """Test the complete pipeline workflow using fixture data."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_workflow_from_all_matches_to_final_output(self, temp_dir):
        """Test pipeline workflow from all_matches.csv to final output."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
            keep_intermediate=True
        )
        
        # Set PAM length manually (normally done in validate_inputs)
        pipeline.pam_length = 3
        
        # Create fake all_matches.csv (simulating guidescan output)
        all_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1,TGCCTCGCGCAGCTCGCGGTGG,0,0,0.400000
guide2,GAGTTCGCTGCGCGCTGTTNGG,3,5000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000"""
        
        pipeline.all_matches_file.write_text(all_matches_content)
        
        # Run post-guidescan workflow
        pipeline.select_best_matches()
        result_df = pipeline.format_exact_matches()
        
        # Verify outputs
        assert pipeline.best_matches_file.exists()
        assert pipeline.final_output_file.exists()
        
        # Verify results
        assert len(result_df) == 2
        
        # guide1 should have exact match from chr1 (not chr2)
        guide1 = result_df[result_df['id'] == 'guide1'].iloc[0]
        assert str(guide1['chromosome']) == '1' or guide1['chromosome'] == 1
        assert guide1['start'] == 1000
        assert guide1['end'] == 1018
        
        # guide2 should have exact match
        guide2 = result_df[result_df['id'] == 'guide2'].iloc[0]
        assert str(guide2['chromosome']) == '3' or guide2['chromosome'] == 3
    
    def test_pipeline_with_mocked_guidescan(self, temp_dir):
        """Test full pipeline with mocked guidescan execution."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
            keep_intermediate=False
        )
        
        # Mock guidescan subprocess call
        def mock_run_guidescan(self):
            """Mock guidescan by creating fake output."""
            all_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000"""
            self.all_matches_file.write_text(all_matches_content)
        
        # Replace run_guidescan with mock
        with patch.object(GuidescanPipeline, 'run_guidescan', mock_run_guidescan):
            # Mock subprocess.run for guidescan availability check
            with patch('subprocess.run'):
                result_df = pipeline.run()
        
        # Verify pipeline completed
        assert result_df is not None
        assert len(result_df) == 1
        assert result_df.iloc[0]['id'] == 'guide1'
        
        # Final output should exist
        assert pipeline.final_output_file.exists()
        
        # Intermediate files should be cleaned up (keep_intermediate=False)
        assert not pipeline.all_matches_file.exists()
        assert not pipeline.best_matches_file.exists()
    
    def test_error_in_workflow_propagates(self, temp_dir):
        """Test that errors in workflow steps propagate correctly."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")
        
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        # Set PAM length
        pipeline.pam_length = 3
        
        # Don't create all_matches.csv - should cause error in select_best_matches
        with pytest.raises(FileNotFoundError):
            pipeline.select_best_matches()

