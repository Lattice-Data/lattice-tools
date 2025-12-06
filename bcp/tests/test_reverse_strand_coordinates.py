"""
Unit tests for reverse strand coordinate calculation.

These tests verify the critical fix for negative strand coordinate handling
where guidescan reports position as (PAM_start + 1), requiring a +2 adjustment
to get the actual protospacer start position.
"""

import pytest
import pandas as pd
from pathlib import Path
import tempfile
import csv
import sys

# Add parent directory to path to import the pipeline module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guidescan_pipeline import GuidescanPipeline


class TestReverseStrandCoordinates:
    """Test reverse strand coordinate calculation."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def pipeline(self, temp_dir):
        """Create a pipeline instance for testing."""
        genome_index = temp_dir / "genome.fa.index"
        genome_index.write_text("# Dummy genome index\n")

        guides_file = temp_dir / "guides.txt"
        guides_file.write_text(
            "id,sequence,pam,chromosome,position,sense\n"
            "test_guide,ACGCGCCGCGCACCGACGT,NGG,,,\n"
        )

        output_dir = temp_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        pipeline = GuidescanPipeline(
            genome_index=str(genome_index),
            guides_file=str(guides_file),
            output_dir=str(output_dir)
        )

        pipeline.pam_length = 3
        return pipeline

    def create_best_matches_file(self, pipeline, guide_id, sequence, chromosome,
                                 position, strand, distance):
        """Helper to create a best_matches.csv file for testing."""
        with open(pipeline.best_matches_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['id', 'sequence', 'match_chrm', 'match_position',
                           'match_strand', 'match_distance', 'match_sequence',
                           'rna_bulges', 'dna_bulges', 'specificity'])
            writer.writerow([guide_id, sequence, chromosome, position, strand,
                           distance, sequence.replace('N', 'A'), 0, 0, 0.5])

    def test_negative_strand_coordinate_calculation(self, pipeline):
        """
        Test that negative strand coordinates are calculated correctly.

        Guidescan reports position as (PAM_start + 1) for negative strand.
        For a 19bp protospacer, we need to add 2 to get the actual start.

        Example: CD24-1
        - Guidescan position: 18992714
        - PAM location: 18992713-18992715 (3bp)
        - Protospacer start: 18992716 (position + 2)
        - Protospacer end: 18992734 (start + 19 - 1)
        """
        # Create test data - simulating CD24-1
        guide_id = "CD24-1"
        sequence = "ACGCGCCGCGCACCGACGTNGG"  # 19bp + NGG
        chromosome = "Y"
        guidescan_position = 18992714  # This is PAM_start + 1
        strand = "-"

        self.create_best_matches_file(
            pipeline, guide_id, sequence, chromosome,
            guidescan_position, strand, distance=0
        )

        # Run format_exact_matches
        result_df = pipeline.format_exact_matches()

        # Verify results
        assert len(result_df) == 1
        row = result_df.iloc[0]

        # Expected: start = position + 2, end = position + 2 + 19 - 1
        expected_start = 18992716
        expected_end = 18992734

        assert int(row['start']) == expected_start, \
            f"Expected start={expected_start}, got {row['start']}"
        assert int(row['end']) == expected_end, \
            f"Expected end={expected_end}, got {row['end']}"
        assert row['sense'] == "-"
        assert row['chromosome'] == chromosome

        # Verify span is correct (19bp protospacer)
        span = int(row['end']) - int(row['start']) + 1
        assert span == 19, f"Expected span of 19bp, got {span}bp"

    def test_positive_strand_coordinate_calculation(self, pipeline):
        """
        Test that positive strand coordinates are calculated correctly.

        For positive strand, guidescan reports the actual protospacer start.

        Example: ARMC5-1
        - Guidescan position: 31459567
        - Protospacer start: 31459567 (same as position)
        - Protospacer end: 31459585 (start + 19 - 1)
        """
        guide_id = "ARMC5-1"
        sequence = "TGCCTCGCGCAGCTCGCGGNGG"  # 19bp + NGG
        chromosome = "16"
        guidescan_position = 31459567
        strand = "+"

        self.create_best_matches_file(
            pipeline, guide_id, sequence, chromosome,
            guidescan_position, strand, distance=0
        )

        result_df = pipeline.format_exact_matches()

        assert len(result_df) == 1
        row = result_df.iloc[0]

        # Expected: start = position, end = position + 19 - 1
        expected_start = 31459567
        expected_end = 31459585

        assert int(row['start']) == expected_start
        assert int(row['end']) == expected_end
        assert row['sense'] == "+"
        assert row['chromosome'] == chromosome

        span = int(row['end']) - int(row['start']) + 1
        assert span == 19

    def test_multiple_negative_strand_guides(self, pipeline):
        """Test multiple negative strand guides to ensure consistent calculation."""
        test_cases = [
            # (guide_id, guidescan_pos, expected_start, expected_end)
            ("DDB1-1", 61333045, 61333047, 61333065),
            ("AP3S1-1", 115842122, 115842124, 115842142),
            ("CD24-1", 18992714, 18992716, 18992734),
        ]

        for guide_id, gs_pos, exp_start, exp_end in test_cases:
            sequence = "ACGCGCCGCGCACCGACGTNGG"
            chromosome = "test"

            self.create_best_matches_file(
                pipeline, guide_id, sequence, chromosome,
                gs_pos, "-", distance=0
            )

            result_df = pipeline.format_exact_matches()
            row = result_df.iloc[0]

            assert int(row['start']) == exp_start, \
                f"{guide_id}: Expected start={exp_start}, got {row['start']}"
            assert int(row['end']) == exp_end, \
                f"{guide_id}: Expected end={exp_end}, got {row['end']}"
            assert row['sense'] == "-"

    def test_negative_strand_offset_is_exactly_20bp(self, pipeline):
        """
        Verify that the coordinate shift for negative strand is exactly 20bp.

        Old calculation: start = position - 19 + 1 = position - 18
        New calculation: start = position + 2
        Difference: (position + 2) - (position - 18) = 20
        """
        guide_id = "test_guide"
        sequence = "ACGCGCCGCGCACCGACGTNGG"
        guidescan_position = 1000000

        # Old (wrong) calculation
        old_start = guidescan_position - 19 + 1  # = 999982
        old_end = guidescan_position  # = 1000000

        # New (correct) calculation via pipeline
        self.create_best_matches_file(
            pipeline, guide_id, sequence, "chr1",
            guidescan_position, "-", distance=0
        )

        result_df = pipeline.format_exact_matches()
        row = result_df.iloc[0]

        new_start = int(row['start'])
        new_end = int(row['end'])

        # Verify the offset is exactly 20bp
        assert new_start - old_start == 20, \
            f"Expected 20bp offset, got {new_start - old_start}bp"
        assert new_end - old_end == 20, \
            f"Expected 20bp offset for end, got {new_end - old_end}bp"

    def test_protospacer_length_variations(self, pipeline):
        """Test that coordinate calculation works for different protospacer lengths."""
        test_cases = [
            (17, "CGCGCCGCGCACCGACG", 17),  # 17bp protospacer
            (18, "ACGCGCCGCGCACCGACG", 18),  # 18bp protospacer
            (19, "ACGCGCCGCGCACCGACGT", 19),  # 19bp protospacer
            (20, "ACGCGCCGCGCACCGACGTA", 20),  # 20bp protospacer
        ]

        for length, protospacer, expected_span in test_cases:
            sequence = protospacer + "NGG"
            guidescan_position = 1000000

            self.create_best_matches_file(
                pipeline, f"guide_{length}bp", sequence, "chr1",
                guidescan_position, "-", distance=0
            )

            result_df = pipeline.format_exact_matches()
            row = result_df.iloc[0]

            # For negative strand: start = position + 2, end = position + 2 + length - 1
            expected_start = guidescan_position + 2
            expected_end = guidescan_position + 2 + length - 1

            assert int(row['start']) == expected_start
            assert int(row['end']) == expected_end

            span = int(row['end']) - int(row['start']) + 1
            assert span == expected_span, \
                f"Expected {expected_span}bp span, got {span}bp"

    def test_na_handling_preserved(self, pipeline):
        """Verify that NA values are handled correctly in coordinate calculation."""
        guide_id = "unmapped_guide"
        sequence = "ACGCGCCGCGCACCGACGTNGG"

        # Create best_matches with NA values
        with open(pipeline.best_matches_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['id', 'sequence', 'match_chrm', 'match_position',
                           'match_strand', 'match_distance', 'match_sequence',
                           'rna_bulges', 'dna_bulges', 'specificity'])
            writer.writerow([guide_id, sequence, "NA", "NA", "NA", "NA",
                           "NA", "NA", "NA", "1.0"])

        result_df = pipeline.format_exact_matches()

        assert len(result_df) == 1
        row = result_df.iloc[0]

        assert row['chromosome'] == "NA"
        assert row['start'] == "NA"
        assert row['end'] == "NA"
        assert row['sense'] == "NA"


class TestCoordinateSystemDocumentation:
    """Tests that serve as documentation for the coordinate system."""

    def test_guidescan_position_meaning_negative_strand(self):
        """
        Document what guidescan's 'position' field means for negative strand.

        For negative strand Cas9 guides, the genomic structure is:
        5'---[PAM 3bp][Protospacer Nbp]---3'

        Guidescan reports: position = PAM_start + 1 (the 2nd base of PAM)

        To get protospacer coordinates:
        - Protospacer starts 2 bases to the RIGHT of guidescan position
        - start = position + 2
        - end = position + 2 + protospacer_length - 1
        """
        # This test documents the expected behavior
        guidescan_position = 100
        pam_length = 3
        protospacer_length = 19

        # PAM location
        pam_start = guidescan_position - 1
        pam_end = guidescan_position + 1
        assert pam_end - pam_start + 1 == pam_length

        # Protospacer location
        protospacer_start = guidescan_position + 2
        protospacer_end = protospacer_start + protospacer_length - 1
        assert protospacer_end - protospacer_start + 1 == protospacer_length

        # Verify PAM is immediately to the LEFT of protospacer
        assert pam_end + 1 == protospacer_start

    def test_guidescan_position_meaning_positive_strand(self):
        """
        Document what guidescan's 'position' field means for positive strand.

        For positive strand Cas9 guides:
        Position = protospacer start (straightforward, 1-based)
        """
        guidescan_position = 100
        protospacer_length = 19

        # For positive strand, position IS the protospacer start
        protospacer_start = guidescan_position
        protospacer_end = guidescan_position + protospacer_length - 1

        assert protospacer_end - protospacer_start + 1 == protospacer_length
