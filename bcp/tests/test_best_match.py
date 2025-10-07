"""
Tests for best match selection logic in the guidescan pipeline.
"""

import pytest
import pandas as pd
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


class TestBestMatchSelection:
    """Test the pick_best_match method that selects optimal genomic matches."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def pipeline(self, temp_dir):
        """Create a pipeline instance for testing."""
        # Create dummy files
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
        
        return pipeline
    
    def test_single_exact_match(self, pipeline):
        """Test with a single exact match (distance=0)."""
        data = {
            'id': ['guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['1'],
            'match_position': [1000],
            'match_strand': ['+'],
            'match_distance': [0]
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        assert result['id'] == 'guide1'
        assert result['match_distance'] == 0
        assert result['match_chrm'] == '1'
    
    def test_multiple_matches_pick_smallest_distance(self, pipeline):
        """Test that smallest distance match is selected from multiple matches."""
        data = {
            'id': ['guide1', 'guide1', 'guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['1', '2', '3'],
            'match_position': [1000, 2000, 3000],
            'match_strand': ['+', '+', '-'],
            'match_distance': [2, 0, 1]
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        # Should pick the match with distance=0 (chr2)
        assert result['match_distance'] == 0
        assert result['match_chrm'] == '2'
        assert result['match_position'] == 2000
    
    def test_prefer_actual_match_over_na(self, pipeline):
        """Test that actual genomic match is preferred over NA (no match) result."""
        data = {
            'id': ['guide1', 'guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['NA', '1'],
            'match_position': ['NA', 1000],
            'match_strand': ['+', '+'],
            'match_distance': ['NA', 0]
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        # Should prefer the actual match (chr1) over NA
        assert result['match_chrm'] == '1'
        assert result['match_distance'] == 0
        assert result['match_position'] == 1000
    
    def test_all_na_matches(self, pipeline):
        """Test handling when all matches are NA (no genomic matches found)."""
        data = {
            'id': ['guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['NA'],
            'match_position': ['NA'],
            'match_strand': ['+'],
            'match_distance': ['NA']
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        # Should return the NA row
        assert result['id'] == 'guide1'
        assert result['match_chrm'] == 'NA'
    
    def test_distance_comparison_with_invalid_values(self, pipeline):
        """Test that invalid distance values are handled (treated as infinity)."""
        data = {
            'id': ['guide1', 'guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['1', '2'],
            'match_position': [1000, 2000],
            'match_strand': ['+', '+'],
            'match_distance': ['invalid', 1]
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        # Should pick the match with valid distance=1 (chr2)
        assert result['match_distance'] == 1
        assert result['match_chrm'] == '2'
    
    def test_exact_match_preferred_over_near_match(self, pipeline):
        """Test that exact match (distance=0) is preferred over near matches."""
        data = {
            'id': ['guide1', 'guide1', 'guide1'],
            'sequence': ['TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG', 'TGCCTCGCGCAGCTCGCGGNGG'],
            'match_chrm': ['1', '2', '3'],
            'match_position': [1000, 2000, 3000],
            'match_strand': ['+', '-', '+'],
            'match_distance': [1, 0, 2]
        }
        df = pd.DataFrame(data)
        
        result = pipeline.pick_best_match(df)
        
        # Should pick exact match (distance=0, chr2)
        assert result['match_distance'] == 0
        assert result['match_chrm'] == '2'
        assert result['match_strand'] == '-'
    
    def test_select_best_matches_integration(self, pipeline, temp_dir):
        """Test the full select_best_matches workflow with fixture data."""
        # Create all_matches.csv with multiple matches per guide
        all_matches_file = temp_dir / "output" / "all_matches.csv"
        all_matches_file.parent.mkdir(exist_ok=True)
        
        all_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1,TGCCTCGCGCAGCTCGCGGTGG,0,0,0.400000
guide2,GAGTTCGCTGCGCGCTGTTNGG,1,5000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide3,GCCCGCTCCCCGCGATCCCNGG,NA,NA,+,NA,NA,NA,NA,NA"""
        
        all_matches_file.write_text(all_matches_content)
        
        # Update pipeline to use this file
        pipeline.all_matches_file = all_matches_file
        pipeline.best_matches_file = temp_dir / "output" / "best_matches.csv"
        
        # Run select_best_matches
        pipeline.select_best_matches()
        
        # Read the result
        result_df = pd.read_csv(pipeline.best_matches_file)
        
        # Verify results
        assert len(result_df) == 3, "Should have 3 guides"
        
        # guide1 should have exact match from chr1 (distance=0, not chr2 with distance=1)
        guide1_row = result_df[result_df['id'] == 'guide1'].iloc[0]
        assert guide1_row['match_distance'] == 0
        # Chromosome can be read as int or string, both are valid
        assert str(guide1_row['match_chrm']) == '1' or guide1_row['match_chrm'] == 1
        
        # guide2 should have its only match
        guide2_row = result_df[result_df['id'] == 'guide2'].iloc[0]
        assert guide2_row['match_distance'] == 0
        assert str(guide2_row['match_chrm']) == '1' or guide2_row['match_chrm'] == 1
        
        # guide3 should have NA (no match)
        guide3_row = result_df[result_df['id'] == 'guide3'].iloc[0]
        # NA can be read as string 'NA' or NaN
        assert guide3_row['match_chrm'] == 'NA' or pd.isna(guide3_row['match_chrm'])

