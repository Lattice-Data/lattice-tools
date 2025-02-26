import argparse
import lattice
import multiprocessing
import traceback
from flattener import main


DEFAULT_MATRIX_TXT = "test_processed_matrix_files.txt"

EPILOG = f"""
Script to test flattener in parallel.
Will create new process per each ProcessedMatrixFile given as input
Can provide txt file with name of one ProcessedMatrixFile per line
Txt file supports commented lines staring with '#' or inline comments with '#'
Use -p/--processed-matrices argument to list ProcMatrixFiles on command line

Final terminal print out will show either:
    ProcessedMatrixFile: SUCCESS
    ProcessedMatrixFile: FAILURE
        Traceback of exception/error

Examples:
    python test_flattener.py --mode prod -p LATDF190KNY LATDF477OUM
    python test_flattener.py --mode prod -f my_matrix_files.txt

Default 'python test_flattener.py' == 'python test_flattener.py -m demo -f {DEFAULT_MATRIX_TXT}'

Currently no guarantee for data integrity when using this for actual flattening in parallel.
Could be race conditions when different processes try to simultaneously
    read the same RawMatrixFile
    download the same RawMatrixFile 
    concurrently test prod and demo with the same set of ProcessedMatrixFiles
"""

connection = None

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--mode", 
        "-m",
        help="Lattice env to test flattening, prod or demo, default demo",
        default="demo"
    )
    parser.add_argument(
        "--processed-matrices", 
        "-p",
        dest="matrices",
        help="Command line list of ProcessedMatrixFiles to test, overrides txt file list if present",
        nargs="+",
        default=list()
    )
    parser.add_argument(
        "--file", 
        "-f",
        help=f"Txt file with list of ProcessedMatrixFiles to test. One ProcMatrix per line. Default is {DEFAULT_MATRIX_TXT}",
        default=DEFAULT_MATRIX_TXT
    )
    parser.add_argument(
        "--hcatier1",
        "-t",
        help="Add flag to include HCA Tier 1 requirements",
        action="store_true",
    )
    args = parser.parse_args()
    return args


args = getArgs()

def set_global_connection(mode):
    global connection
    if not connection:
        connection = lattice.Connection(mode)


def flatten(file):
    try:
        main(file, connection, args.hcatier1)
        return (file, "SUCCESS")
    except SystemExit as e:
        return (file, e)
    except Exception as e:
        return (file, e)


def make_file_list(args):
    if args.matrices:
        return args.matrices
    else:
        with open(args.file, "r") as f:
            return [line.partition("#")[0].strip() for line in f if not line.startswith("#")]


if __name__ == "__main__":
    results = list()
    files = make_file_list(args)
    workers = len(files)
    with multiprocessing.Pool(initializer=set_global_connection, initargs=(args.mode,), processes=workers) as pool:
        iterator = pool.imap(flatten, files)
        for matrix in iterator:
            try:
                results.append(matrix)
            except StopIteration:
                break
            except Exception as e:
                print(f"ERROR: {e}")
                results.append((matrix, e))

    print("FINAL RESULTS:")
    print("=" * 80)
    for file, returned_obj in results:
        if isinstance(returned_obj, str):
            print(f"{file}: SUCCESS")
        else:
            print(f"{file}: FAILURE")
            traceback.print_exception(None, returned_obj, returned_obj.__traceback__)
            print("=" * 80)
