"""
Tests for loading file manifests from CSV/TSV files.
"""
import os
import pytest
from qa_mods import load_files_from_manifest


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')


class TestLoadFilesFromManifest:
    """Tests for load_files_from_manifest function."""

    def test_load_tsv_no_header(self):
        """Test loading TSV manifest without header, S3 in column 0."""
        manifest_path = os.path.join(FIXTURES_DIR, 'test_manifest.tsv')
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=manifest_path,
            delimiter='\t',
            s3_column=0,
            has_header=False
        )
        
        # Should have 6 raw files (4 from GROUP1, 2 from GROUP2)
        assert len(all_raw_files) == 6
        
        # Should have 2 groups with processed files
        assert len(all_proc_files) == 2
        assert 'GROUP1' in all_proc_files
        assert 'GROUP2' in all_proc_files
        
        # GROUP1 should have 2 processed files
        assert len(all_proc_files['GROUP1']) == 2
        # GROUP2 should have 2 processed files
        assert len(all_proc_files['GROUP2']) == 2

    def test_load_csv_with_header(self):
        """Test loading CSV manifest with header, S3 in column 1."""
        manifest_path = os.path.join(FIXTURES_DIR, 'test_manifest_with_header.csv')
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=manifest_path,
            delimiter=',',
            s3_column=1,
            has_header=True
        )
        
        # Should have 3 raw files
        assert len(all_raw_files) == 3
        
        # Should have 2 groups with processed files
        assert len(all_proc_files) == 2
        assert 'GROUP1' in all_proc_files
        assert 'GROUP2' in all_proc_files

    def test_s3_prefix_stripping(self):
        """Verify s3://bucket/ prefix is correctly removed."""
        manifest_path = os.path.join(FIXTURES_DIR, 'test_manifest.tsv')
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=manifest_path,
            delimiter='\t',
            s3_column=0,
            has_header=False
        )
        
        # None of the paths should start with 's3://'
        for f in all_raw_files:
            assert not f.startswith('s3://'), f"Path still has S3 prefix: {f}"
        
        for group_files in all_proc_files.values():
            for f in group_files:
                assert not f.startswith('s3://'), f"Path still has S3 prefix: {f}"
        
        # Paths should start with the key after bucket
        assert all_raw_files[0] == 'proj/order/GROUP1/raw/RUN01-GROUP1_GEX-UG01-BC01.csv'

    def test_raw_processed_separation(self):
        """Verify raw files go to list, processed files grouped by group name."""
        manifest_path = os.path.join(FIXTURES_DIR, 'test_manifest.tsv')
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=manifest_path,
            delimiter='\t',
            s3_column=0,
            has_header=False
        )
        
        # All raw files should contain '/raw/'
        for f in all_raw_files:
            assert '/raw/' in f, f"Raw file missing '/raw/' in path: {f}"
        
        # All processed files should contain '/processed/'
        for group_files in all_proc_files.values():
            for f in group_files:
                assert '/processed/' in f, f"Processed file missing '/processed/' in path: {f}"

    def test_group_extraction(self):
        """Verify group name is correctly extracted from processed file paths."""
        manifest_path = os.path.join(FIXTURES_DIR, 'test_manifest.tsv')
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=manifest_path,
            delimiter='\t',
            s3_column=0,
            has_header=False
        )
        
        # GROUP1 processed files should have GROUP1 in path
        for f in all_proc_files['GROUP1']:
            assert 'GROUP1/processed/' in f
        
        # GROUP2 processed files should have GROUP2 in path
        for f in all_proc_files['GROUP2']:
            assert 'GROUP2/processed/' in f

    def test_empty_results_for_no_matching_paths(self, tmp_path):
        """Test that empty manifest returns empty results."""
        # Create an empty manifest
        empty_manifest = tmp_path / "empty.tsv"
        empty_manifest.write_text("")
        
        # This should handle empty file gracefully
        # Note: pandas will raise an error for truly empty file, 
        # so we test with a file that has no raw/processed paths
        other_manifest = tmp_path / "other.tsv"
        other_manifest.write_text("s3://bucket/some/other/path.txt\t/orig\n")
        
        all_raw_files, all_proc_files = load_files_from_manifest(
            manifest_path=str(other_manifest),
            delimiter='\t',
            s3_column=0,
            has_header=False
        )
        
        assert len(all_raw_files) == 0
        assert len(all_proc_files) == 0
