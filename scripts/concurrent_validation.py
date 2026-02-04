import argparse
import multiprocessing
import os
import subprocess
from pathlib import Path


CPU_COUNT = os.cpu_count()
FILE = Path(__file__).resolve()
DIR = FILE.parent
SCRIPT_NAME = FILE.name 
PRINT_WIDTH = 113

EPILOG = f"""
Script to run CXG validation in parallel.
By default will collect h5ad files from lattice-tools/scripts directory
Use argument -d --directory to specify other location
Use argument -r --revised to only collect h5ad files with the '_revised.h5ad' suffix
Use argument -t --testfile to select h5ads from a test_flattener input txt file

Examples:
    python {SCRIPT_NAME} -d /mnt/test_files -r
    python {SCRIPT_NAME} --revised --directory /Users/me/Documents/curation/files
    python {SCRIPT_NAME} --testfile test_processed_matrix_files.txt
"""

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--directory", 
        "-d",
        type=Path,
        help="Directory to collect h5ad files, default is lattice-tools/scripts or location of this script",
        default=DIR,
    )
    parser.add_argument(
        "--testfile", 
        "-t",
        help="Filter h5ads based on flattener testing input file",
        default="",
    )
    parser.add_argument(
        "--revised", 
        "-r",
        help="Add -r/--revised flag to filter for only '_revised.h5ad' files in directory",
        action="store_true",
    )
    args = parser.parse_args()
    return args


def make_file_list(args):
    file_suffix = "_revised.h5ad" if args.revised else ".h5ad"

    if args.testfile:
        testfile_path = Path(DIR) / args.testfile
        assert testfile_path.suffix == ".txt", "Test file needs to be txt file"
        with open(testfile_path, "r") as f:
            test_files = [line.partition("#")[0].strip() for line in f if not line.startswith("#")]

        tested_h5ads = [
            f 
            for f in args.directory.iterdir() 
            if any(
                test_file in f 
                for test_file in test_files
            ) 
            and f.name.endswith(file_suffix) 
        ]

        return tested_h5ads
    else:
        files = [f for f in args.directory.iterdir() if f.name.endswith(file_suffix)]
        return files


ARGS = getArgs()
files = make_file_list(ARGS)
workers = min(len(files), CPU_COUNT)

def validate(file_name):
    full_path = ARGS.directory / file_name
    validate = subprocess.run(
        ["cellxgene-schema", "validate", full_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    print(f"Validation for {file_name}...")
    for line in validate.stdout.decode('utf-8').split('\n'):
        print(line)
    for line in validate.stderr.decode('utf-8').split('\n'):
        print(line)
    print("=" * PRINT_WIDTH + "\n")


def validate_all_files(files):
    with multiprocessing.Pool(processes=workers) as pool:
        pool.map(validate, files)


if __name__ == "__main__":
    revised_str = " REVISED" if ARGS.revised else ""
    print(f"\nFound {len(files)}{revised_str} h5ad(s) in {ARGS.directory} to validate")
    print("=" * PRINT_WIDTH + "\n")
    validate_all_files(files)
