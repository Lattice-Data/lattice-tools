"""
Tests for multiple gene overlap scenarios with bioframe implementation.
Verifies that guides overlapping 2+ genes are correctly identified and reported.
"""

import pytest
import pandas as pd
import tempfile
from pathlib import Path
import sys

# Add parent directory to path to import the guide_to_gene_bioframe module
sys.path.insert(0, str(Path(__file__).parent.parent))

from guide_to_gene_bioframe import annotate_guides_bioframe


class TestMultipleGeneOverlaps:
    """Test scenarios where a single guide overlaps multiple genes."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_guide_overlapping_two_genes(self, temp_dir):
        """Test guide that overlaps exactly 2 genes."""
        # Create a GTF with 2 overlapping genes
        gtf_file = temp_dir / "overlapping_genes.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
            "chr1\ttest\tgene\t1500\t2500\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n"
        )
        
        # Create a guide that overlaps both genes (position 1600-1618)
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1600,1618,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # Should have 2 rows (one for each gene)
        assert len(df) == 2
        
        # Both rows should be for guide1
        assert all(df['id'] == 'guide1')
        
        # Should have both gene IDs
        gene_ids = set(df['gene_id'])
        assert gene_ids == {'GENE1', 'GENE2'}
        
        # Should have both gene names
        gene_names = set(df['gene_name'])
        assert gene_names == {'GENE1', 'GENE2'}
    
    def test_guide_overlapping_three_genes(self, temp_dir):
        """Test guide that overlaps 3 genes."""
        # Create a GTF with 3 overlapping genes
        gtf_file = temp_dir / "three_overlapping_genes.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
            "chr1\ttest\tgene\t1400\t2200\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n"
            "chr1\ttest\tgene\t1600\t2500\t.\t-\t.\tgene_id \"GENE3\"; gene_name \"GENE3\";\n"
        )
        
        # Create a guide that overlaps all 3 genes (position 1700-1718)
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1700,1718,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # Should have 3 rows (one for each gene)
        assert len(df) == 3
        
        # All rows should be for guide1
        assert all(df['id'] == 'guide1')
        
        # Should have all 3 gene IDs
        gene_ids = set(df['gene_id'])
        assert gene_ids == {'GENE1', 'GENE2', 'GENE3'}
    
    def test_minimum_1bp_overlap(self, temp_dir):
        """Test that even 1bp overlap is detected."""
        # Create a GTF with a gene (GTF uses 1-based coordinates)
        # Gene spans positions 1000-2000 inclusive
        gtf_file = temp_dir / "gene.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
        )
        
        # Create guides with different overlap amounts
        # Input CSV uses 1-based inclusive coordinates
        # After conversion: CSV start-1 to CSV end (0-based start, 1-based end)
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            # guide1: 1000-1019 overlaps gene 1000-2000 (overlaps at position 1000)
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1000,1019,+\n"
            # guide2: 2000-2019 overlaps gene 1000-2000 (overlaps at position 2000)
            "guide2,ACGTACGTACGTACGTACG,NGG,1,2000,2019,+\n"
            # guide3: 980-999 NO overlap with gene 1000-2000
            "guide3,ACGTACGTACGTACGTACG,NGG,1,980,999,+\n"
            # guide4: 2001-2020 NO overlap with gene 1000-2000
            "guide4,ACGTACGTACGTACGTACG,NGG,1,2001,2020,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # guide1 and guide2 should overlap (at boundary positions)
        overlapping = df[df['gene_id'] == 'GENE1']
        assert len(overlapping) == 2, f"Expected 2 overlapping guides, got {len(overlapping)}: {overlapping['id'].tolist()}"
        assert set(overlapping['id']) == {'guide1', 'guide2'}
        
        # guide3 and guide4 should not overlap (should be NA)
        non_overlapping = df[df['gene_id'] == 'NA']
        assert len(non_overlapping) == 2, f"Expected 2 non-overlapping guides, got {len(non_overlapping)}: {non_overlapping['id'].tolist()}"
        assert set(non_overlapping['id']) == {'guide3', 'guide4'}
    
    def test_multiple_guides_multiple_genes(self, temp_dir):
        """Test multiple guides with varying numbers of gene overlaps."""
        # Create a GTF with overlapping genes
        gtf_file = temp_dir / "genes.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
            "chr1\ttest\tgene\t1500\t2500\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n"
            "chr1\ttest\tgene\t3000\t4000\t.\t+\t.\tgene_id \"GENE3\"; gene_name \"GENE3\";\n"
        )
        
        # Create guides with different overlap patterns
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            # Overlaps only GENE1
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1100,1118,+\n"
            # Overlaps GENE1 and GENE2
            "guide2,ACGTACGTACGTACGTACG,NGG,1,1600,1618,+\n"
            # Overlaps only GENE2
            "guide3,ACGTACGTACGTACGTACG,NGG,1,2100,2118,+\n"
            # Overlaps only GENE3
            "guide4,ACGTACGTACGTACGTACG,NGG,1,3500,3518,+\n"
            # No overlaps (intergenic)
            "guide5,ACGTACGTACGTACGTACG,NGG,1,2600,2618,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # Total should be 6 rows: 1+2+1+1+1 = 6
        assert len(df) == 6
        
        # Check guide1: 1 row, GENE1
        guide1_rows = df[df['id'] == 'guide1']
        assert len(guide1_rows) == 1
        assert guide1_rows.iloc[0]['gene_id'] == 'GENE1'
        
        # Check guide2: 2 rows, GENE1 and GENE2
        guide2_rows = df[df['id'] == 'guide2']
        assert len(guide2_rows) == 2
        assert set(guide2_rows['gene_id']) == {'GENE1', 'GENE2'}
        
        # Check guide3: 1 row, GENE2
        guide3_rows = df[df['id'] == 'guide3']
        assert len(guide3_rows) == 1
        assert guide3_rows.iloc[0]['gene_id'] == 'GENE2'
        
        # Check guide4: 1 row, GENE3
        guide4_rows = df[df['id'] == 'guide4']
        assert len(guide4_rows) == 1
        assert guide4_rows.iloc[0]['gene_id'] == 'GENE3'
        
        # Check guide5: 1 row, NA (intergenic)
        guide5_rows = df[df['id'] == 'guide5']
        assert len(guide5_rows) == 1
        assert guide5_rows.iloc[0]['gene_id'] == 'NA'
    
    def test_overlaps_preserve_original_guide_info(self, temp_dir):
        """Test that all overlaps preserve original guide information."""
        # Create a GTF with 2 overlapping genes
        gtf_file = temp_dir / "genes.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
            "chr1\ttest\tgene\t1500\t2500\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n"
        )
        
        # Create a guide that overlaps both genes
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1600,1618,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # Both rows should have identical guide information
        assert all(df['id'] == 'guide1')
        assert all(df['sequence'] == 'ACGTACGTACGTACGTACG')
        assert all(df['pam'] == 'NGG')
        # Chromosome should be string (pandas may read as int, so convert)
        assert all(df['chromosome'].astype(str) == '1')
        assert all(df['start'].astype(str) == '1600')
        assert all(df['end'].astype(str) == '1618')
        assert all(df['sense'] == '+')
        
        # Only gene_id and gene_name should differ
        assert df.iloc[0]['gene_id'] != df.iloc[1]['gene_id']
        assert df.iloc[0]['gene_name'] != df.iloc[1]['gene_name']
    
    def test_real_data_multiple_overlaps(self, temp_dir):
        """Test with real fixture data that has known multiple overlaps."""
        guides_file = Path(__file__).parent / "fixtures" / "guides_no_genes.csv"
        gtf_file = Path(__file__).parent / "fixtures" / "genes.gtf"
        output_file = temp_dir / "output.csv"
        
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Read output
        df = pd.read_csv(output_file, keep_default_na=False)
        
        # Count guides with multiple overlaps
        guide_counts = df.groupby('id').size()
        guides_with_multiple_overlaps = guide_counts[guide_counts > 1]
        
        # Verify that some guides have multiple overlaps
        # (based on the real data, we know there are overlaps)
        assert len(guides_with_multiple_overlaps) > 0
        
        # Verify that for each guide with multiple overlaps,
        # all rows have the same guide info but different gene info
        for guide_id, count in guides_with_multiple_overlaps.items():
            guide_rows = df[df['id'] == guide_id]
            
            # Same guide info
            assert len(guide_rows['sequence'].unique()) == 1
            assert len(guide_rows['pam'].unique()) == 1
            assert len(guide_rows['chromosome'].unique()) == 1
            assert len(guide_rows['start'].unique()) == 1
            assert len(guide_rows['end'].unique()) == 1
            assert len(guide_rows['sense'].unique()) == 1
            
            # Different gene info (at least one should be different)
            assert len(guide_rows['gene_id'].unique()) == count
            
    def test_statistics_with_multiple_overlaps(self, temp_dir, capsys):
        """Test that statistics correctly report multiple overlaps."""
        # Create a GTF with overlapping genes
        gtf_file = temp_dir / "genes.gtf"
        gtf_file.write_text(
            "chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n"
            "chr1\ttest\tgene\t1500\t2500\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n"
        )
        
        # Create guides
        guides_file = temp_dir / "guides.csv"
        guides_file.write_text(
            "id,sequence,pam,chromosome,start,end,sense\n"
            # Overlaps both genes
            "guide1,ACGTACGTACGTACGTACG,NGG,1,1600,1618,+\n"
            # Overlaps only GENE1
            "guide2,ACGTACGTACGTACGTACG,NGG,1,1100,1118,+\n"
        )
        
        output_file = temp_dir / "output.csv"
        annotate_guides_bioframe(str(guides_file), str(gtf_file), str(output_file))
        
        # Capture output
        captured = capsys.readouterr()
        
        # Verify statistics
        assert "Total guides processed: 2" in captured.out
        assert "Guides with gene overlaps: 2" in captured.out
        assert "Total gene overlaps found: 3" in captured.out  # 2 + 1
        assert "Average overlaps per guide (for guides with overlaps): 1.50" in captured.out

