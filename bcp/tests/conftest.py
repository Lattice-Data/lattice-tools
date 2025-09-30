"""
Shared pytest fixtures and configuration for guidescan pipeline tests.
"""

import pytest
import tempfile
from pathlib import Path


@pytest.fixture
def temp_dir():
    """
    Create a temporary directory for test files.
    
    Yields:
        Path: Temporary directory path that will be cleaned up after test.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def dummy_genome_index(temp_dir):
    """
    Create a dummy genome index file for testing.
    
    Args:
        temp_dir: Temporary directory fixture
        
    Returns:
        Path: Path to dummy genome index file
    """
    index_file = temp_dir / "genome.fa.index"
    index_file.write_text("# Dummy genome index for testing\n")
    return index_file


@pytest.fixture
def sample_guides_file_ngg(temp_dir):
    """
    Create a sample guides file with uniform NGG PAM.
    
    Args:
        temp_dir: Temporary directory fixture
        
    Returns:
        Path: Path to guides file
    """
    guides_file = temp_dir / "guides_ngg.txt"
    guides_file.write_text(
        "id,sequence,pam,chromosome,position,sense\n"
        "guide1,TGCCTCGCGCAGCTCGCGG,NGG,,,\n"
        "guide2,GAGTTCGCTGCGCGCTGTT,NGG,,,\n"
        "guide3,GCCCGCTCCCCGCGATCCC,NGG,,,\n"
    )
    return guides_file


@pytest.fixture
def sample_all_matches_csv(temp_dir):
    """
    Create a sample all_matches.csv file (simulating guidescan output).
    
    Args:
        temp_dir: Temporary directory fixture
        
    Returns:
        Path: Path to all_matches.csv file
    """
    all_matches_file = temp_dir / "all_matches.csv"
    content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide1,TGCCTCGCGCAGCTCGCGGNGG,2,2000,-,1,TGCCTCGCGCAGCTCGCGGTGG,0,0,0.400000
guide2,GAGTTCGCTGCGCGCTGTTNGG,1,5000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide3,GCCCGCTCCCCGCGATCCCNGG,NA,NA,+,NA,NA,NA,NA,NA"""
    all_matches_file.write_text(content)
    return all_matches_file


@pytest.fixture
def sample_best_matches_csv(temp_dir):
    """
    Create a sample best_matches.csv file.
    
    Args:
        temp_dir: Temporary directory fixture
        
    Returns:
        Path: Path to best_matches.csv file
    """
    best_matches_file = temp_dir / "best_matches.csv"
    content = """id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity
guide1,TGCCTCGCGCAGCTCGCGGNGG,1,1000,+,0,TGCCTCGCGCAGCTCGCGGCGG,0,0,0.500000
guide2,GAGTTCGCTGCGCGCTGTTNGG,1,5000,+,0,GAGTTCGCTGCGCGCTGTTGGG,0,0,0.500000
guide3,GCCCGCTCCCCGCGATCCCNGG,NA,NA,+,NA,NA,NA,NA,NA"""
    best_matches_file.write_text(content)
    return best_matches_file


# Configure pytest markers
def pytest_configure(config):
    """
    Register custom pytest markers.
    """
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )

