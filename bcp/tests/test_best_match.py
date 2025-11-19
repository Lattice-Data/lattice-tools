"""
Tests for best match selection logic in the guidescan pipeline.
Tests the new behavior: retains ALL exact matches, or single best non-exact match.
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
    """Test the select_best_matches method that retains all exact matches."""
    
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
        output_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        return pipeline
    
    def create_all_matches_file(self, pipeline, content):
        """Helper to create an all_matches.csv file with given content."""
        pipeline.all_matches_file.write_text(content)
    
    def run_and_read_results(self, pipeline):
        """Helper to run select_best_matches and read results."""
        pipeline.select_best_matches()
        return pd.read_csv(pipeline.best_matches_file)
    
    def test_single_exact_match(self, pipeline):
        """Test with a single exact match (distance=0)."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        assert len(result_df) == 1
        assert result_df.iloc[0]['id'] == 'guide1'
        assert result_df.iloc[0]['match_distance'] == 0
        assert str(result_df.iloc[0]['match_chrm']) == '1'
    
    def test_multiple_exact_matches_keep_all(self, pipeline):
        """Test that ALL exact matches are retained when multiple exist."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,0
guide1,TGCCTCGCGCAGCTCGCGGNGG,3,3000,+,0"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should keep ALL 3 exact matches
        assert len(result_df) == 3
        guide1_matches = result_df[result_df['id'] == 'guide1']
        assert len(guide1_matches) == 3
        assert all(guide1_matches['match_distance'] == 0)
        chromosomes = set(str(c) for c in guide1_matches['match_chrm'])
        assert chromosomes == {'1', '2', '3'}
    
    def test_mixed_exact_and_nonexact_keep_only_exact(self, pipeline):
        """Test that only exact matches are kept when both exact and non-exact exist."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1
guide1,TGCCTCGCGCAGCTCGCGGNGG,3,3000,+,0"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should keep only the 2 exact matches (chr1 and chr3), not chr2 with distance=1
        assert len(result_df) == 2
        guide1_matches = result_df[result_df['id'] == 'guide1']
        assert len(guide1_matches) == 2
        assert all(guide1_matches['match_distance'] == 0)
        chromosomes = set(str(c) for c in guide1_matches['match_chrm'])
        assert chromosomes == {'1', '3'}
    
    def test_only_nonexact_matches_result_in_na(self, pipeline):
        """Test that guides with only non-exact matches result in NA (treated as no match)."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,2
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1
guide1,TGCCTCGCGCAGCTCGCGGNGG,3,3000,+,3"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should have NA result since no exact matches exist
        assert len(result_df) == 1
        assert result_df.iloc[0]['id'] == 'guide1'
        # Check for NA - can be string 'NA' or NaN
        assert result_df.iloc[0]['match_chrm'] == 'NA' or pd.isna(result_df.iloc[0]['match_chrm'])
        assert result_df.iloc[0]['match_distance'] == 'NA' or pd.isna(result_df.iloc[0]['match_distance'])
    
    def test_prefer_actual_match_over_na(self, pipeline):
        """Test that actual genomic match is preferred over NA (no match) result."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,NA,NA,+,NA
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should prefer the actual exact match (chr1) over NA
        assert len(result_df) == 1
        # Chromosome can be '1', 1, or 1.0 depending on pandas parsing
        chrm = str(result_df.iloc[0]['match_chrm'])
        assert chrm == '1' or chrm == '1.0'
        assert result_df.iloc[0]['match_distance'] == 0
        assert result_df.iloc[0]['match_position'] == 1000
    
    def test_all_na_matches(self, pipeline):
        """Test handling when all matches are NA (no genomic matches found)."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,NA,NA,+,NA"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should return one NA row
        assert len(result_df) == 1
        assert result_df.iloc[0]['id'] == 'guide1'
        assert result_df.iloc[0]['match_chrm'] == 'NA' or pd.isna(result_df.iloc[0]['match_chrm'])
    
    def test_distance_comparison_with_invalid_values(self, pipeline):
        """Test that invalid distance values are handled correctly - treated as non-exact matches."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,invalid
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,+,1"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # Should result in NA since neither is an exact match (distance=0)
        assert len(result_df) == 1
        assert result_df.iloc[0]['id'] == 'guide1'
        # Check for NA - can be string 'NA' or NaN
        assert result_df.iloc[0]['match_chrm'] == 'NA' or pd.isna(result_df.iloc[0]['match_chrm'])
        assert result_df.iloc[0]['match_distance'] == 'NA' or pd.isna(result_df.iloc[0]['match_distance'])
    
    def test_multiple_guides_integration(self, pipeline):
        """Test the full select_best_matches workflow with multiple guides."""
        content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1,TGCCTCGCGCAGCTCGCGGTGG,0,0,0.400000
guide2,GAGTTCGCTGCGCGCTGTTNGG,1,5000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide2,GAGTTCGCTGCGCGCTGTTNGG,5,6000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide3,GCCCGCTCCCCGCGATCCCNGG,NA,NA,+,NA,NA,NA,NA,NA"""
        
        self.create_all_matches_file(pipeline, content)
        result_df = self.run_and_read_results(pipeline)
        
        # guide1: has 1 exact and 1 non-exact match -> keep only exact (chr1)
        guide1_matches = result_df[result_df['id'] == 'guide1']
        assert len(guide1_matches) == 1
        assert guide1_matches.iloc[0]['match_distance'] == 0
        # Chromosome can be '1', 1, or 1.0 depending on pandas parsing
        chrm = str(guide1_matches.iloc[0]['match_chrm'])
        assert chrm == '1' or chrm == '1.0'
        
        # guide2: has 2 exact matches -> keep BOTH
        guide2_matches = result_df[result_df['id'] == 'guide2']
        assert len(guide2_matches) == 2
        assert all(guide2_matches['match_distance'] == 0)
        # Chromosome can be parsed as int or float, handle both
        chromosomes = set(str(c).rstrip('0').rstrip('.') for c in guide2_matches['match_chrm'])
        assert chromosomes == {'1', '5'}
        
        # guide3: has only NA -> keep 1 NA result
        guide3_matches = result_df[result_df['id'] == 'guide3']
        assert len(guide3_matches) == 1
        assert guide3_matches.iloc[0]['match_chrm'] == 'NA' or pd.isna(guide3_matches.iloc[0]['match_chrm'])

