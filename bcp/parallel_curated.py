import argparse
import multiprocessing
import traceback
import sys
import subprocess
import os


EPILOG = f"""
Script run curate_matrices.py in a parallel fashion

Example:
    python curate_matrices.py --bucket czi-psomagen --sheet your_sheet_id --project marson-mapping-grns-perturb-seq --groupfile groups.txt
    python curate_matrices.py --bucket czi-psomagen --sheet your_sheet_id --project marson-mapping-grns-perturb-seq --grouplist CD4i_R1L01 CD4i_R1L02

For more details:
    python %(prog)s --help

"""


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bucket", 
        "-b",
        help="bucket where h5 files live",
        default="czi-psomagen"
    )
    parser.add_argument(
        "--sheet", 
        "-s",
        help="google sheet id of metadata"
    )
    parser.add_argument(
        "--project", 
        "-p",
        help="Name of the project directory on S3 bucket"
    )
    parser.add_argument(
        "--groupfile",
        "-f",
        help="the file containing GroupID of the sample that for curated matrices generation"
    )
    parser.add_argument(
        "--grouplist",
        "-l",
        help="a commandline list containing GroupID of the sample that for curated matrices generation",
        nargs="+",
        default=list()
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit()
    return args


def create_arguments(args, groups):
    '''
    Creates list of tuples that are the commandlines for curate_matrices.py

    return: list of tuple commandlines
    '''
    commands_list = []
    for i in range(len(groups)):
        new_tuple = ("--bucket", args.bucket, "--sheet", args.sheet, "--project", args.project, "--group", groups[i])
        commands_list.append(new_tuple)
    return commands_list


def make_groups_list(args):
    '''
    Allows comments in the groupfile, which will be skipped

    return: list of GroupIDs
    '''
    if args.grouplist:
        return args.grouplist
    else:
        with open(args.groupfile, "r") as f:
            return [line for line in f if not line.startswith("#")]


# Define the script you want to run in parallel
CURATE_SCRIPT = "curate_matrices.py"

def run_curate(commands_tuple):
    command = [sys.executable, CURATE_SCRIPT] + list(commands_tuple)
    print(f"Starting command: {' '.join(command)}")

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, 'PYTHONUNBUFFERED': '1'}
        )
        return f"SUCCESS: Command `{' '.join(command)}` finished.\n{result.stdout.strip()}"
    except subprocess.CalledProcessError as e:
        return f"ERROR: Command `{' '.join(command)}` failed with exit code {e.returncode}.\n{e.stderr.strip()}"
    except Exception as e:
        return f"ERROR: An unexpected error occurred for command `{' '.join(command)}`: {e}"



args = getArgs()


if __name__ == "__main__":
    results = list()
    groups = []
    groups = make_groups_list(args)
    commands_list = create_arguments(args, groups)
    workers = min(len(groups), 50)
    with multiprocessing.Pool(processes=workers) as pool:
        iterator = pool.imap(run_curate, commands_list)
        for matrix in iterator:
            try:
                results.append(matrix)
            except StopIteration:
                break
            except Exception as e:
                print(f"ERROR: {e}")
                results.append((matrix, e))

    print("FINAL RESULTS:")
    success_counts = 0
    fail_counts = 0
    for result in results:
        print(result)
            


    
