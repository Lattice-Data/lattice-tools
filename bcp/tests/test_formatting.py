"""
Tests for exact match formatting and coordinate calculation in the guidescan pipeline.
"""

import pytest
import pandas as pd
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


class TestExactMatchFormatting:
    """Test the format_exact_matches method that processes and formats genomic coordinates."""
    
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
        
        # Set PAM length (normally done in validate_inputs)
        pipeline.pam_length = 3
        
        return pipeline
    
    def test_sequence_splitting_3bp_pam(self, pipeline):
        """Test splitting full sequence into protospacer and 3bp PAM."""
        full_sequence = "TGCCTCGCGCAGCTCGCGGNGG"  # 19bp + 3bp
        
        protospacer_length = len(full_sequence) - pipeline.pam_length
        protospacer = full_sequence[:protospacer_length]
        pam = full_sequence[protospacer_length:]
        
        assert protospacer == "TGCCTCGCGCAGCTCGCGG"
        assert pam == "NGG"
        assert len(protospacer) == 19
        assert len(pam) == 3
    
    def test_sequence_splitting_4bp_pam(self, pipeline):
        """Test splitting full sequence into protospacer and 4bp PAM."""
        pipeline.pam_length = 4
        full_sequence = "TGCCTCGCGCAGCTCGCGGNNGG"  # 19bp + 4bp
        
        protospacer_length = len(full_sequence) - pipeline.pam_length
        protospacer = full_sequence[:protospacer_length]
        pam = full_sequence[protospacer_length:]
        
        assert protospacer == "TGCCTCGCGCAGCTCGCGG"
        assert pam == "NNGG"
        assert len(protospacer) == 19
        assert len(pam) == 4
    
    def test_coordinate_calculation_plus_strand(self, pipeline):
        """Test genomic coordinate calculation for + strand match."""
        position = 1000  # Match position from guidescan
        protospacer_length = 19
        strand = "+"
        
        # For + strand: start = position, end = position + length - 1
        start_pos = position
        end_pos = position + protospacer_length - 1
        
        assert start_pos == 1000
        assert end_pos == 1018
        assert end_pos - start_pos + 1 == protospacer_length
    
    def test_coordinate_calculation_minus_strand(self, pipeline):
        """Test genomic coordinate calculation for - strand match."""
        position = 1000  # Match position from guidescan (PAM_start + 1)
        protospacer_length = 19
        strand = "-"

        # For - strand: guidescan reports position as (PAM_start + 1)
        # Protospacer starts 2 bases to the RIGHT of this position
        start_pos = position + 2
        end_pos = position + 2 + protospacer_length - 1

        assert start_pos == 1002
        assert end_pos == 1020
        assert end_pos - start_pos + 1 == protospacer_length
    
    def test_exact_match_only_distance_zero(self, pipeline, temp_dir):
        """Test that only exact matches (distance=0) are included in output."""
        # Create best_matches.csv with exact and near matches
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide2,GAGTTCGCTGCGCGCTGTTNGG,2,2000,+,1,GAGTTCGCTGCGCGCTGTTAGG,0,0,0.400000
guide3,GCCCGCTCCCCGCGATCCCNGG,3,3000,-,0,GCCCGCTCCCCGCGATCCCCGG,0,0,0.500000"""
        
        pipeline.best_matches_file.write_text(best_matches_content)
        
        # Run format_exact_matches
        result_df = pipeline.format_exact_matches()
        
        # Only guides with distance=0 should have genomic coordinates
        exact_matches = result_df[result_df['chromosome'] != 'NA']
        
        assert len(exact_matches) == 2, "Should have 2 exact matches (guide1 and guide3)"
        assert 'guide1' in exact_matches['id'].values
        assert 'guide3' in exact_matches['id'].values
        
        # guide2 should have NA coordinates (distance=1, not exact)
        guide2_row = result_df[result_df['id'] == 'guide2'].iloc[0]
        assert guide2_row['chromosome'] == 'NA'
        assert guide2_row['start'] == 'NA'
        assert guide2_row['end'] == 'NA'
    
    def test_no_match_results_in_na(self, pipeline, temp_dir):
        """Test that guides with no genomic matches get NA coordinates."""
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,NA,NA,+,NA,NA,NA,NA,NA"""
        
        pipeline.best_matches_file.write_text(best_matches_content)
        
        result_df = pipeline.format_exact_matches()
        
        assert len(result_df) == 1
        row = result_df.iloc[0]
        
        assert row['id'] == 'guide1'
        assert row['sequence'] == 'TGCCTCGCGCAGCTCGCGG'  # Protospacer only
        assert row['pam'] == 'NGG'
        assert row['chromosome'] == 'NA'
        assert row['start'] == 'NA'
        assert row['end'] == 'NA'
        assert row['sense'] == 'NA'
    
    def test_format_exact_matches_integration_plus_strand(self, pipeline, temp_dir):
        """Test full formatting workflow with + strand exact match."""
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000"""
        
        pipeline.best_matches_file.write_text(best_matches_content)
        
        result_df = pipeline.format_exact_matches()
        
        assert len(result_df) == 1
        row = result_df.iloc[0]
        
        # Check all output columns
        assert row['id'] == 'guide1'
        assert row['sequence'] == 'TGCCTCGCGCAGCTCGCGG'  # 19bp protospacer
        assert row['pam'] == 'NGG'  # 3bp PAM
        assert str(row['chromosome']) == '1' or row['chromosome'] == 1
        assert row['start'] == 1000
        assert row['end'] == 1018  # 1000 + 19 - 1
        assert row['sense'] == '+'
    
    def test_format_exact_matches_integration_minus_strand(self, pipeline, temp_dir):
        """Test full formatting workflow with - strand exact match."""
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,-,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000"""

        pipeline.best_matches_file.write_text(best_matches_content)

        result_df = pipeline.format_exact_matches()

        assert len(result_df) == 1
        row = result_df.iloc[0]

        # Check coordinate calculation for - strand
        assert row['id'] == 'guide1'
        assert row['sequence'] == 'TGCCTCGCGCAGCTCGCGG'
        assert row['pam'] == 'NGG'
        assert str(row['chromosome']) == '1' or row['chromosome'] == 1
        assert row['start'] == 1002  # position + 2
        assert row['end'] == 1020  # position + 2 + 19 - 1
        assert row['sense'] == '-'
    
    def test_format_exact_matches_mixed_results(self, pipeline, temp_dir):
        """Test formatting with mix of exact matches, near matches, and no matches."""
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide2,GAGTTCGCTGCGCGCTGTTNGG,2,2000,-,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide3,GCCCGCTCCCCGCGATCCCNGG,3,3000,+,1,GCCCGCTCCCCGCGATCCCAGG,0,0,0.400000
guide4,CTGTGGGGCCCTGTCCATGNGG,NA,NA,+,NA,NA,NA,NA,NA"""
        
        pipeline.best_matches_file.write_text(best_matches_content)
        
        result_df = pipeline.format_exact_matches()
        
        assert len(result_df) == 4
        
        # guide1: exact match, + strand
        g1 = result_df[result_df['id'] == 'guide1'].iloc[0]
        assert g1['chromosome'] != 'NA'
        assert g1['start'] == 1000
        assert g1['end'] == 1018
        assert g1['sense'] == '+'
        
        # guide2: exact match, - strand
        g2 = result_df[result_df['id'] == 'guide2'].iloc[0]
        assert g2['chromosome'] != 'NA'
        assert g2['start'] == 2002  # position + 2
        assert g2['end'] == 2020  # position + 2 + 19 - 1
        assert g2['sense'] == '-'
        
        # guide3: near match (distance=1), should have NA
        g3 = result_df[result_df['id'] == 'guide3'].iloc[0]
        assert g3['chromosome'] == 'NA'
        assert g3['sequence'] == 'GCCCGCTCCCCGCGATCCC'
        assert g3['pam'] == 'NGG'
        
        # guide4: no match, should have NA
        g4 = result_df[result_df['id'] == 'guide4'].iloc[0]
        assert g4['chromosome'] == 'NA'
        assert g4['sequence'] == 'CTGTGGGGCCCTGTCCATG'
        assert g4['pam'] == 'NGG'
    
    def test_output_file_created(self, pipeline, temp_dir):
        """Test that exact_matches_formatted.csv is created with correct structure."""
        best_matches_content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000"""
        
        pipeline.best_matches_file.write_text(best_matches_content)
        
        result_df = pipeline.format_exact_matches()
        
        # Check file was created
        assert pipeline.final_output_file.exists()
        
        # Read the file and verify structure
        saved_df = pd.read_csv(pipeline.final_output_file)
        
        expected_columns = ['id', 'sequence', 'pam', 'chromosome', 'start', 'end', 'sense']
        assert list(saved_df.columns) == expected_columns
        assert len(saved_df) == 1

