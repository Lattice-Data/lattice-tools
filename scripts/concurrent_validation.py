import argparse
import logging
import logging.config
import logging.handlers
import multiprocessing
import os
import subprocess
import sys
from datetime import datetime
from multiprocessing import Queue
from pathlib import Path

import yaml

CPU_COUNT = os.cpu_count()
FILE = Path(__file__).resolve()
DIR = FILE.parent
LOG_QUEUE = None
SCRIPT_NAME = FILE.name
PRINT_WIDTH = 113
TERMINAL_ERROR_LIMIT = 100

EPILOG = f"""
Script to run CXG validation in parallel.
By default will collect h5ad files from lattice-tools/scripts directory
Use argument -d --directory to specify other location
Use argument -r --revised to only collect h5ad files with the '_revised.h5ad' suffix
Use argument -t --testfile to select h5ads from a test_flattener input txt file
Use argument -pa --pre-analysis to run pre-analysis validation
Use arugment -i --ignore-labels to ignore ontological labels when validating
Use arugment -lo --log-output to log validation output to file

Examples:
    python {SCRIPT_NAME} -d /mnt/test_files -r
    python {SCRIPT_NAME} --revised --directory /Users/me/Documents/curation/files
    python {SCRIPT_NAME} --testfile test_processed_matrix_files.txt
    python {SCRIPT_NAME} --directory /home/shared/home/curated_matrices --pre-analysis --ignore-labels
    python {SCRIPT_NAME} --directory /home/shared/home/curated_matrices -pa -i --log-output
"""

def getArgs() -> argparse.Namespace:
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
        help="Filter h5ads based on flattener testing input file. Test file should be in same directory as h5ads",
        default="",
    )
    parser.add_argument(
        "--revised", 
        "-r",
        help="Add -r/--revised flag to filter for only '_revised.h5ad' files in directory",
        action="store_true",
    )
    parser.add_argument(
        "--ignore-labels",
        "-i",
        help="Ignore ontology labels when validating",
        action="store_const",
        const="-i",
        default="",
    )
    parser.add_argument(
        "--pre-analysis",
        "-pa",
        help="Include pre-analysis validation requirements in validating the data",
        action="store_const",
        const="-pa",
        default="",
    )
    parser.add_argument(
        "--log-output", 
        "-lo",
        help="Log validation output to file",
        action="store_true",
    )
    args = parser.parse_args()
    return args


def create_logger(queue: Queue) -> None:
    """
    Function to create logger in process workers
    For stream and file handling, best to set root logger in process
    and then let the listener handle logger from queue

    Now with yaml logger config, using fragment curator yaml with some 
    modifications in the final __main__ block
    """
    handler = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG)
    root.propagate = False


def worker_init(queue: Queue) -> None:
    """
    Init function to get logging queue to individual process workers.
    Different from fragment curator where manager is used
    """
    global LOG_QUEUE
    LOG_QUEUE = queue
    create_logger(LOG_QUEUE)


def make_file_list(args: argparse.Namespace) -> list[Path]:
    """
    Create list of paths for h5ad files to be validated
    """
    file_suffix = "_revised.h5ad" if args.revised else ".h5ad"

    if args.testfile:
        testfile_path = Path(args.directory) / args.testfile
        assert testfile_path.suffix == ".txt", "Test file needs to be txt file"
        with open(testfile_path, "r") as f:
            test_files = [line.partition("#")[0].strip() for line in f if not line.startswith("#")]

        tested_h5ads = [
            f 
            for f in args.directory.iterdir()
            if any(
                test_file in f.name 
                for test_file in test_files
            ) 
            and f.name.endswith(file_suffix)
        ]

        return tested_h5ads
    else:
        files = [f for f in args.directory.iterdir() if f.name.endswith(file_suffix)]
        return files


def validate(file_name: Path) -> subprocess.CompletedProcess | None:
    """
    Worker function to validate with CXG validator
    Supports pass-through pre-analysis and ignore labels flags
    """
    full_path = ARGS.directory / file_name

    # only add pre-analysis and ignore-labels if they exist, otherwise validtor errors
    command = ["cellxgene-schema", "validate"]
    for arg in [ARGS.pre_analysis, ARGS.ignore_labels]:
        if arg:
            command.append(arg)
    command.append(full_path)

    try:
        validate = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logging.exception(f"Error with subprocess: {e}")
    except Exception as e:
        logging.exception(f"Exception occurred: {e}")

    if validate:
        logging.info(f"Validation for {file_name}...")
        validation_status = validate.stderr.splitlines()[-1]
        logging.info(validation_status)

    logging.info("=" * PRINT_WIDTH + "\n")

    return validate


def validate_all_files(files: list[Path]) -> list[subprocess.CompletedProcess]:
    """
    Function to create worker pool.
    Uses init function to allow workers access to logging queue
    """
    with multiprocessing.Pool(
        processes=workers,
        initializer=worker_init,
        initargs=(logging_queue,),
    ) as pool:
        results = pool.map(validate, files)

    return results


# need in global namespace for process workers to access values
ARGS = getArgs()


if __name__ == "__main__":
    files = make_file_list(ARGS)
    workers = min(len(files), CPU_COUNT)
    logging_queue = Queue(maxsize=-1)

    # resuing fragment curator logger with some changes for here
    with open(DIR / "fragment_curator_mods" / "log_config_fragment_curator.yaml", "rt") as f:
        logging_config = yaml.safe_load(f.read())

    time_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_name = f"{time_date}_outfile_concurrent_validation.log"
    logging_config["handlers"]["file"]["filename"] = log_file_name
    # no message format to just pass CXG validation logging to console and file
    logging_config["formatters"]["detailed"]["format"] = "%(message)s"
    logging_config["formatters"]["simple"]["format"] = "%(message)s"
    logging_config["loggers"] = {
        "concurrent_validation": {
            "level": "DEBUG",
            "handlers": ["console", "file"],
            "propagate": False
        }
    }
    logging.config.dictConfig(logging_config)

    # set up concurrent validation logger and first log messages
    logger = logging.getLogger("concurrent_validation")
    listener = logging.handlers.QueueListener(logging_queue, *logger.handlers, respect_handler_level=True)
    listener.start()
    logger.debug("Logger thread started")
    logger.debug(f"Command line args: {' '.join(sys.argv)}")

    revised_str = " REVISED" if ARGS.revised else ""
    logger.info(f"\nFound {len(files)}{revised_str} h5ad(s) in {ARGS.directory} to validate")
    logger.info("=" * PRINT_WIDTH + "\n")
    logger.info("===== Summary Validation Results =====\n")

    results = validate_all_files(files)

    logger.info("===== Final Validation Results =====\n")
    over_error_limit = False
    for completed in results:
        # file will be last arg in args list, preserved as Path, so get name attribute
        file_name = completed.args[-1].name
        errors = completed.stderr.splitlines()
        num_errors = len(errors)
        logger.info(f"Validation for {file_name}...")
        for error in errors[:TERMINAL_ERROR_LIMIT]:
            print(error)
        if num_errors > TERMINAL_ERROR_LIMIT:
            over_error_limit = True
            print("...")
            print(f"{num_errors} total validation errors")
            print(f"Full error log available in log file '{log_file_name}'")
        for error in errors:
            logger.debug(error)

        logger.info("=" * PRINT_WIDTH + "\n")

    logger.debug("Logger thread shutdown")
    listener.stop()

    # just remove log file at end if unwanted, instead of complicated logging config modification
    delete_log = (not ARGS.log_output) and (not over_error_limit)
    if delete_log:
        Path(log_file_name).unlink()
