"""
CRISPR Guide Processing Pipeline
Processes guide sequences using guidescan and formats the results.
"""

import argparse
import subprocess
import csv
import pandas as pd
import sys
import os
import logging
from pathlib import Path

# Configure logging - simple format without colors
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger = logging.getLogger(__name__)


class GuidescanPipeline:
    """Pipeline for processing CRISPR guide sequences."""
    
    def __init__(self, 
                 genome_index: str,
                 guides_file: str,
                 output_dir: str,
                 keep_intermediate: bool = False,
                 verbose: bool = False):
        """
        Initialize the pipeline.
        
        Args:
            genome_index: Path to genome index file
            guides_file: Path to input guides file
            output_dir: Directory for output files
            keep_intermediate: Keep intermediate files
            verbose: Enable verbose logging
        """
        self.genome_index = Path(genome_index)
        self.guides_file = Path(guides_file)
        self.output_dir = Path(output_dir)
        self.keep_intermediate = keep_intermediate
        
        # Fixed parameters - not configurable
        self.max_mismatches = 1
        self.max_off_targets = 1
        
        # PAM length will be determined from input file
        self.pam_length = None
        
        # Set logging level
        if verbose:
            logger.setLevel(logging.DEBUG)
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define output file paths
        self.all_matches_file = self.output_dir / "all_matches.csv"
        self.best_matches_file = self.output_dir / "best_matches.csv"
        self.final_output_file = self.output_dir / "exact_matches_formatted.csv"
    
    def validate_and_get_pam_length(self) -> int:
        """
        Validate that all PAMs in the guides file are identical and return PAM length.
        Reads the PAM column from the CSV file and checks for uniformity.
        
        Returns:
            Length of the PAM sequence
            
        Raises:
            ValueError: If PAMs are not uniform across all guides
        """
        logger.info("Validating PAM sequences from input guides...")
        
        pams = set()
        
        try:
            import csv
            with open(self.guides_file, 'r') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    if 'pam' not in row:
                        raise ValueError("Input file must have a 'pam' column")
                    
                    pam = row['pam'].strip()
                    if pam:  # Skip empty PAM values
                        pams.add(pam)
            
            if len(pams) == 0:
                raise ValueError("No PAM sequences found in input file")
            
            if len(pams) > 1:
                raise ValueError(
                    f"All PAM sequences must be identical. Found multiple PAMs: {sorted(pams)}. "
                    f"Please ensure all guides use the same PAM sequence."
                )
            
            # All PAMs are identical, get the single PAM
            pam_sequence = pams.pop()
            pam_length = len(pam_sequence)
            
            logger.info(f"PAM sequence validated: {pam_sequence} (length: {pam_length})")
            
            return pam_length
            
        except Exception as e:
            logger.error(f"Failed to validate PAM sequences: {e}")
            raise
    
    def validate_inputs(self) -> None:
        """Validate input files and parameters."""
        # Check if genome index files exist
        # For real guidescan indices, check for .forward, .reverse, .gs files
        # For test/dummy indices, allow simple file existence
        forward_file = Path(str(self.genome_index) + ".forward")
        reverse_file = Path(str(self.genome_index) + ".reverse")
        gs_file = Path(str(self.genome_index) + ".gs")
        
        has_real_index = (forward_file.exists() and reverse_file.exists() and gs_file.exists())
        has_dummy_index = self.genome_index.exists()
        
        if not (has_real_index or has_dummy_index):
            raise FileNotFoundError(
                f"Genome index not found: {self.genome_index}\n"
                f"For real indices, expected: {forward_file}, {reverse_file}, {gs_file}"
            )
        
        if not self.guides_file.exists():
            raise FileNotFoundError(f"Guides file not found: {self.guides_file}")
        
        # Check if guidescan is available
        try:
            subprocess.run(["guidescan", "--help"], 
                         capture_output=True, 
                         check=False)
        except FileNotFoundError:
            raise RuntimeError("guidescan not found in PATH. Please install guidescan2.")
        
        # Validate PAM uniformity and get PAM length
        self.pam_length = self.validate_and_get_pam_length()
    
    def run_guidescan(self) -> None:
        """Run guidescan enumerate command."""
        logger.info("Running guidescan enumerate...")
        
        cmd = [
            "guidescan", "enumerate",
            str(self.genome_index),
            "-f", str(self.guides_file),
            "-o", str(self.all_matches_file),
            "--format", "csv",
            "--mode", "complete",
            "-m", str(self.max_mismatches),
            "--max-off-targets", str(self.max_off_targets)
        ]
        
        logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                logger.debug(f"guidescan output:\n{result.stdout}")
            
            logger.info(f"guidescan completed. Results saved to {self.all_matches_file}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"guidescan failed with exit code {e.returncode}")
            logger.error(f"Error output: {e.stderr}")
            raise RuntimeError(f"guidescan enumeration failed: {e.stderr}")
    
    def select_best_matches(self) -> None:
        """
        Select matches for each guide.
        
        For each guide:
        - If exact matches exist (distance = 0): Keep ALL exact matches
        - If no exact matches exist: Keep one NA result (treat as no match)
        """
        logger.info("Selecting matches for each guide (retaining all exact matches only)...")
        
        try:
            df = pd.read_csv(self.all_matches_file)
            
            if df.empty:
                logger.warning("No matches found in guidescan output")
                # Create empty best matches file
                df.to_csv(self.best_matches_file, index=False)
                return
            
            # Distance is column 6 (index 5)
            distance_col = df.columns[5]
            guide_id_col = df.columns[0]
            chromosome_col = df.columns[2]
            
            # Convert distance to numeric for comparison
            df['distance_numeric'] = pd.to_numeric(df[distance_col], errors='coerce')
            
            results = []
            for guide_id, group in df.groupby(guide_id_col):
                # Check for NA matches (no genomic match found)
                has_na_only = (group[chromosome_col] == "NA").all()
                
                if has_na_only:
                    # Keep one NA result
                    results.append(group.head(1))
                else:
                    # Filter out NA results
                    real_matches = group[group[chromosome_col] != "NA"]
                    
                    # Filter for exact matches (distance = 0)
                    exact_matches = real_matches[real_matches['distance_numeric'] == 0]
                    
                    if len(exact_matches) > 0:
                        # Keep ALL exact matches
                        results.append(exact_matches)
                    else:
                        # No exact matches - keep one NA result (treat as no match)
                        # Use the first NA row if available, otherwise create NA result
                        na_rows = group[group[chromosome_col] == "NA"]
                        if len(na_rows) > 0:
                            results.append(na_rows.head(1))
                        else:
                            # Create a NA result from the first real match by replacing match data with NA
                            na_result = real_matches.head(1).copy()
                            na_result[chromosome_col] = "NA"
                            # Set position, distance to NA (columns 3, 5)
                            na_result[df.columns[3]] = "NA"
                            na_result[distance_col] = "NA"
                            na_result['distance_numeric'] = pd.NA
                            results.append(na_result)
            
            result_df = pd.concat(results, ignore_index=True)
            
            # Drop the temporary numeric distance column before saving
            result_df = result_df.drop(columns=['distance_numeric'])
            
            result_df.to_csv(self.best_matches_file, index=False)
            
            # Log statistics
            guide_count = df[guide_id_col].nunique()
            result_count = len(result_df)
            exact_match_count = len(result_df[result_df[chromosome_col] != "NA"])
            logger.info(f"Matches saved to {self.best_matches_file}")
            logger.info(f"Processed {guide_count} guides: {exact_match_count} exact matches, {guide_count - exact_match_count} no matches")
            
        except Exception as e:
            logger.error(f"Failed to select matches: {e}")
            raise
    
    def format_exact_matches(self) -> pd.DataFrame:
        """
        Process exact matches from best matches and format output.
        
        Returns:
            DataFrame with formatted results
        """
        logger.info("Formatting exact matches...")
        
        results = []
        
        try:
            with open(self.best_matches_file, 'r') as infile:
                reader = csv.reader(infile)
                
                # Skip header
                header = next(reader, None)
                
                for row in reader:
                    if not row or all(cell.strip() == '' for cell in row):
                        continue
                    
                    if len(row) < 6:
                        continue
                    
                    id_val = row[0]
                    full_sequence = row[1]
                    chromosome = row[2] if len(row) > 2 else "NA"
                    position_str = row[3] if len(row) > 3 else "NA"
                    strand = row[4] if len(row) > 4 else "+"
                    match_distance_str = row[5] if len(row) > 5 else "NA"
                    
                    # Split sequence into protospacer and PAM
                    sequence_length = len(full_sequence)
                    protospacer_length = sequence_length - self.pam_length
                    
                    protospacer = full_sequence[:protospacer_length]
                    pam = full_sequence[protospacer_length:]
                    
                    # Process exact matches only (distance = 0)
                    if chromosome != "NA" and position_str != "NA" and match_distance_str != "NA":
                        try:
                            match_distance = int(float(match_distance_str))
                            position = int(float(position_str))
                            
                            if match_distance == 0:
                                # Calculate genomic coordinates based on strand
                                if strand == "+":
                                    start_pos = position
                                    end_pos = position + protospacer_length - 1
                                    sense = "+"
                                else:
                                    # For negative strand, guidescan reports (PAM_start + 1)
                                    # which is 2 bases to the LEFT of the protospacer start
                                    # So we add 2 to get the actual protospacer start position
                                    start_pos = position + 2
                                    end_pos = position + 2 + protospacer_length - 1
                                    sense = "-"
                                
                                results.append([
                                    id_val, protospacer, pam, 
                                    chromosome, start_pos, end_pos, sense
                                ])
                            else:
                                # Non-exact match
                                results.append([
                                    id_val, protospacer, pam, 
                                    "NA", "NA", "NA", "NA"
                                ])
                        except (ValueError, TypeError) as e:
                            logger.warning(f"Error processing row for {id_val}: {e}")
                            results.append([
                                id_val, protospacer, pam, 
                                "NA", "NA", "NA", "NA"
                            ])
                    else:
                        # No match found
                        results.append([
                            id_val, protospacer, pam, 
                            "NA", "NA", "NA", "NA"
                        ])
            
            # Convert to DataFrame
            df = pd.DataFrame(
                results, 
                columns=['id', 'sequence', 'pam', 'chromosome', 'start', 'end', 'sense']
            )
            
            # Save to CSV
            df.to_csv(self.final_output_file, index=False)
            logger.info(f"Processing complete. Output saved to: {self.final_output_file}")
            
            return df
            
        except Exception as e:
            logger.error(f"Failed to format exact matches: {e}")
            raise
    
    def cleanup_intermediate_files(self) -> None:
        """Remove intermediate files if requested."""
        if not self.keep_intermediate:
            logger.info("Cleaning up intermediate files...")
            
            for file in [self.all_matches_file, self.best_matches_file]:
                if file.exists():
                    file.unlink()
                    logger.debug(f"Removed {file}")
    
    def run(self) -> pd.DataFrame:
        """
        Run the complete pipeline.
        
        Returns:
            DataFrame with final results
        """
        logger.info("Starting CRISPR guide processing pipeline...")
        
        try:
            # Validate inputs and detect PAM length
            self.validate_inputs()
            
            # Step 1: Run guidescan
            self.run_guidescan()
            
            # Step 2: Select best matches
            self.select_best_matches()
            
            # Step 3: Format exact matches
            df = self.format_exact_matches()
            
            # Cleanup if requested
            self.cleanup_intermediate_files()
            
            logger.info("Pipeline completed successfully!")
            
            # Print summary
            total_guides = len(df)
            exact_matches = len(df[df['chromosome'] != 'NA'])
            logger.info(f"Summary: {total_guides} guides processed, {exact_matches} exact matches found")
            
            return df
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="Process CRISPR guide sequences using guidescan",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i genome.fa.index -g guides.txt -o results/
  %(prog)s -i genome.fa.index -g guides.txt -o results/ --keep-intermediate
  %(prog)s -i genome.fa.index -g guides.txt -o results/ -v
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--genome-index',
        required=True,
        help='Path to genome index file (.fa.index)'
    )
    
    parser.add_argument(
        '-g', '--guides',
        required=True,
        help='Path to input guides file (text file with one guide per line)'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for results'
    )
    
    # Optional arguments
    parser.add_argument(
        '--keep-intermediate',
        action='store_true',
        help='Keep intermediate files (all_matches.csv, best_matches.csv)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize and run pipeline
        pipeline = GuidescanPipeline(
            genome_index=args.genome_index,
            guides_file=args.guides,
            output_dir=args.output_dir,
            keep_intermediate=args.keep_intermediate,
            verbose=args.verbose
        )
        
        results = pipeline.run()
        
        # Print first few results if verbose
        if args.verbose and not results.empty:
            print("\nFirst 5 results:")
            print(results.head())
        
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("Pipeline interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()