"""
Tests for PAM validation and length calculation in the guidescan pipeline.
"""

import pytest
from pathlib import Path
import tempfile
import sys

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


class TestPAMValidation:
    """Test PAM validation and length calculation from CSV column."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def dummy_genome_index(self, temp_dir):
        """Create a dummy genome index file."""
        index_file = temp_dir / "genome.fa.index"
        index_file.write_text("# Dummy genome index\n")
        return index_file
    
    def test_uniform_pam_ngg(self, temp_dir, dummy_genome_index):
        """Test validation with uniform NGG PAM (3bp)."""
        guides_file = temp_dir / "guides_ngg.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NGG,,,\n"
            "guide3,GCCCGCTCCCCGCGATCCC,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        # PAM length should be determined during validate_inputs
        pipeline.validate_inputs()
        
        assert pipeline.pam_length == 3, f"Expected PAM length 3, got {pipeline.pam_length}"
    
    def test_uniform_pam_nngg(self, temp_dir, dummy_genome_index):
        """Test validation with uniform NNGG PAM (4bp)."""
        guides_file = temp_dir / "guides_nngg.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NNGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NNGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        pipeline.validate_inputs()
        
        assert pipeline.pam_length == 4, f"Expected PAM length 4, got {pipeline.pam_length}"
    
    def test_uniform_pam_nnngg(self, temp_dir, dummy_genome_index):
        """Test validation with uniform NNNGG PAM (5bp)."""
        guides_file = temp_dir / "guides_nnngg.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NNNGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NNNGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        pipeline.validate_inputs()
        
        assert pipeline.pam_length == 5, f"Expected PAM length 5, got {pipeline.pam_length}"
    
    def test_mixed_pams_raises_error(self, temp_dir, dummy_genome_index):
        """Test that mixed PAM sequences raise ValueError."""
        guides_file = temp_dir / "guides_mixed.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
            "guide2,GAGTTCGCTGCGCGCTGTT,NNGG,,,\n"
            "guide3,GCCCGCTCCCCGCGATCCC,NNNGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        with pytest.raises(ValueError) as excinfo:
            pipeline.validate_inputs()
        
        assert "All PAM sequences must be identical" in str(excinfo.value)
        assert "NGG" in str(excinfo.value) or "NNGG" in str(excinfo.value)
    
    def test_empty_guides_file_raises_error(self, temp_dir, dummy_genome_index):
        """Test that empty guides file (just header) raises ValueError."""
        guides_file = temp_dir / "guides_empty.txt"
        guides_file.write_text("id,sequence,pam,chromosome,position,sense\n")
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        with pytest.raises(ValueError) as excinfo:
            pipeline.validate_inputs()
        
        assert "No PAM sequences found" in str(excinfo.value)
    
    def test_missing_pam_column_raises_error(self, temp_dir, dummy_genome_index):
        """Test that missing PAM column raises ValueError."""
        guides_file = temp_dir / "guides_no_pam.txt"
        guides_file.write_text(
            "id,sequence,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        with pytest.raises(ValueError) as excinfo:
            pipeline.validate_inputs()
        
        assert "must have a 'pam' column" in str(excinfo.value)
    
    def test_pam_length_used_in_sequence_splitting(self, temp_dir, dummy_genome_index):
        """Test that calculated PAM length correctly splits protospacer from PAM."""
        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        )
        
        output_dir = temp_dir / "output"
        
        pipeline = GuidescanPipeline(
            genome_index=str(dummy_genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )
        
        pipeline.validate_inputs()
        
        # Full sequence from guidescan would be protospacer + PAM
        # Example: TGCCTCGCGCAGCTCGCGGNGG (19bp + 3bp = 22bp)
        full_sequence = "TGCCTCGCGCAGCTCGCGGNGG"
        protospacer_length = len(full_sequence) - pipeline.pam_length
        
        protospacer = full_sequence[:protospacer_length]
        pam = full_sequence[protospacer_length:]
        
        assert len(protospacer) == 19, f"Expected protospacer length 19, got {len(protospacer)}"
        assert len(pam) == 3, f"Expected PAM length 3, got {len(pam)}"
        assert pam == "NGG", f"Expected PAM 'NGG', got '{pam}'"

