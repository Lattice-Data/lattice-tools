import argparse
import boto3.session
import fsspec
import h5py
import logging
import logging.config
import logging.handlers
import matplotlib.pyplot as plt
import multiprocessing
import os
import pandas as pd
import re
import subprocess
import sys
import threading
import traceback
import yaml

from anndata._io.specs import read_elem
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from multiprocessing import (
    Queue,
    Manager
)
from pathlib import Path
from typing import Callable

from lattice import (
    Connection,
    get_report,
    parse_ids,
)


EPILOG = """
This script will filter and concatenate ATAC fragment files from Lattice.
This is done in parallel as much as possible.

Examples:
    Coming soon

    python %(prog)s --help
"""
BARCODE_PATTERN = r"[ACGT]{16}([-_]\d{1,2})?"
REPLACE_WITH = "B@RCODE"
FRAGMENT_DIR = Path("atac_fragments/")
FS = fsspec.filesystem("s3")
DOWNLOAD_THREADS = 8
NUM_PROCESS_WORKERS = os.cpu_count() // 2
PRINT_WIDTH = 57
PROCESSED_MATRIX_FIELD_LIST = [
    "accession",
    "s3_uri",
    "cell_label_location",
    "cell_label_mappings",
]
RAW_MATRIX_FIELD_LIST = [
    "accession",
    "s3_uri",
    "fragment_file_s3_uri"
]
FRAGMENT_COL_NAMES = [
    "chrom",
    "start",
    "end",
    "barcode",
    "readSupport"
]

@dataclass
class URIPath:
    """
    Dataclass to hold various parsed information for each S3 URI
    Only needs full uri to create the various attributes
    """

    full_uri: str

    def __post_init__(self):
        self.bucket_name: str = self.full_uri.split("/")[2]
        self.file_path = self.full_uri.replace("s3://{}/".format(self.bucket_name), "")
        self.file_name = self.full_uri.split("/")[-1]
        self.parent_path = self.file_path.replace(self.file_name, "")
        self.full_parent_path = self.full_uri.replace(self.file_name, "")


@dataclass
class WorkerQueues:
    """
    Queues for workers to use
    Organization for FragmentFileMeta.queues
    """
    logging_queue: Queue
    compression_queue: Queue


@dataclass
class FragmentFileMeta:
    """
    Dataclass to pass fragment file meta to workers for download and processing
    """

    cell_label_location: str
    label: str
    accession: str
    uri: URIPath
    barcodes: pd.Series
    queues: WorkerQueues

    def __post_init__(self):
        self.download_file_name = self.accession + "_" + self.uri.file_name
        # probably better to add this logic somewhere else
        self.filtered_fragment_path_name = FRAGMENT_DIR / self.download_file_name.replace("fragments.tsv.gz", "filtered_fragments.tsv")
        self.is_file_local = (FRAGMENT_DIR / self.download_file_name).is_file()
        # saved file will likely be compressed
        self.is_filtered_file_local = (FRAGMENT_DIR / self.filtered_fragment_path_name.name.replace("tsv", "tsv.gz")).is_file()


@dataclass
class FragmentWorkerResult:
    """
    Dataclass for returning worker results, plot currently unused
    """
    accession: str
    file_saved: bool
    stats: dict
    plot: plt.figure
    output_path: Path
    checked_duplicates: bool = False
    has_duplicates: bool = False
    duplicates: pd.DataFrame | None = None


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--file",
        "-f",
        help="Any identifier for the matrix of interest."
    )
    parser.add_argument(
        "--mode",
        "-m",
        help="The machine to run on."
    )
    parser.add_argument(
        "--deduplicate",
        "-d",
        help="Remove duplicates from filtered fragment files",
        action="store_true",
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def logger_thread(queue: Queue):
    """
    Logger queue thread to listen and take in logging records
    To get queue to process workers, it's easiest to pass queue 
    as FragmentFileMeta attribute
    Process workers their own logger created with create_logger()
    Thread workers and other functions can look to global scope
    """
    while True:
        try:
            record = queue.get()
            if record is None:
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)
        except Exception:
            import sys
            import traceback
            print("Whoops! Problem:", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            # need to clean up listener to exit infinite exception loop
            logger_thread_shutdown_and_exit()


def logger_thread_shutdown_and_exit():
    """
    Shutdown and cleanup to prevent broken process pipe/queue during sys.exit()
    Call this whenever sys.exit() is needed
    Prevents infinite loop of exception catching in logger process
    Can also use atexit callback, but this function seems cleaner in terms of 
    limiting other exceptions and tracebacks
    """
    logger.critical("sys.exit() call, shutting down")
    logger.critical("Logger thread executing error shutdown")
    queues.logging_queue.put_nowait(None)
    sys.exit()


def create_logger(queue: Queue):
    """
    Function to create logger in process workers
    For stream and file handling, best to set root logger in process
    and then let the listener handle logger from queue

    Now with yaml logger config, and concatenator logger variable for
    global scope. This allows for main process/main thread level to INFO
    and concatenator logging to DEBUG. This prevents tens of thousands of 
    boto3 and s3 debug messages during multithreaded download, and other possible
    module-level debug messages from flooding the log. Worker module debug messages
    could also flood logs in this setup, but that seems acceptable to troublshoot
    any weird parallelism issues within the workers
    """
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    # with worker pool, recycled worker already has handler, don't need to add another
    # without check, messages will double and triple
    if not root.handlers:
        root.addHandler(h)
    root.setLevel(logging.DEBUG)
    root.propagate = False


def download_object(s3_client, fragment_meta: FragmentFileMeta):
    """
    Thread worker to download files from S3
    """
    download_path = FRAGMENT_DIR / fragment_meta.download_file_name
    logger.info(f"Downloading {fragment_meta.download_file_name} to {download_path}")
    s3_client.download_file(
        fragment_meta.uri.bucket_name,
        fragment_meta.uri.file_path,
        str(download_path)
    )
    return "Success"


def download_parallel_multithreading(files_to_download: list[FragmentFileMeta]):
    """
    Function to orchestrate thread workers for parallel download
    Logging looks to global scope
    """
    # Create a session and use it to make our client
    session = boto3.session.Session()
    s3_client = session.client("s3")

    if not FRAGMENT_DIR.is_dir():
        FRAGMENT_DIR.mkdir()

    # Dispatch work tasks with our s3_client
    with ThreadPoolExecutor(max_workers=DOWNLOAD_THREADS) as executor:
        future_to_key = {
            executor.submit(download_object, s3_client, key): key.download_file_name
            for key in files_to_download
            if not key.is_file_local
        }

        num_all_files = len(files_to_download)
        num_locally = num_all_files - len(future_to_key)
        logger.info(f"{num_all_files} raw fragment files in Lattice")
        logger.info(f"Found {num_locally} raw files locally, downloading {len(future_to_key)} files")

        if not future_to_key:
            logger.info("All raw files local, no downloading needed")

        for future in futures.as_completed(future_to_key):
            key = future_to_key[future]
            exception = future.exception()

            if not exception:
                yield key, future.result()
            else:
                yield key, exception


def query_lattice(processed_matrix_accession: str, connection: Connection, queues: WorkerQueues) -> list[FragmentFileMeta]:
    """
    Query Lattice for metadata on processed matrix, raw matrices, fragment file S3 URIs

    Return FragmentFileMeta dataclass for use with downstream functions
    """
    processed_matrix_report = get_report(
        "ProcessedMatrixFile",
        f"&@id=/processed-matrix-files/{processed_matrix_accession}/",
        PROCESSED_MATRIX_FIELD_LIST,
        connection
    )[0]

    h5ad_uri = URIPath(processed_matrix_report["s3_uri"])
    with h5py.File(FS.open(h5ad_uri.full_uri)) as f:
        barcodes = read_elem(f["obs"]).index.to_series()

    logger.debug(f"h5ad for processed matrix file: {h5ad_uri.full_uri}")

    index_patterns = {re.sub(BARCODE_PATTERN, REPLACE_WITH, b) for b in barcodes}

    logger.debug("set of barcode patterns:")
    for pattern in index_patterns:
        logger.debug(pattern)

    master_fragment_file_meta = []
    files_missing = False
    for raw_matrix_meta in processed_matrix_report["cell_label_mappings"]:
        obj_type, filter_url = parse_ids([raw_matrix_meta["raw_matrix"]])
        raw_matrix_report = get_report(
            obj_type,
            filter_url,
            field_lst=RAW_MATRIX_FIELD_LIST,
            connection=connection
        )[0]

        accession = raw_matrix_report["accession"]
        fragment_uri = raw_matrix_report.get("fragment_file_s3_uri", None)

        if fragment_uri is None:
            logger.error(f"raw matrix {accession} does not have fragment file S3 URI")
            files_missing = True

        if not FS.isfile(fragment_uri):
            logger.error(f"raw matrix {accession} fragment file does not exist on S3")
            files_missing = True
            continue

        master_fragment_file_meta.append(
            FragmentFileMeta(
                cell_label_location=processed_matrix_report["cell_label_location"],
                label=raw_matrix_meta["label"],
                accession=accession,
                uri=URIPath(fragment_uri),
                barcodes=barcodes,
                queues=queues
            )
        )

    if files_missing:
        logger.critical("Missing S3 URIs and/or files, exiting...")
        logger_thread_shutdown_and_exit()

    return master_fragment_file_meta


def download_fragment_files(files_to_download: list[FragmentFileMeta]) -> None:
    """
    Call multi-threaded download function and check all files are now local
    """
    for key, result in download_parallel_multithreading(files_to_download):
        logger.info(f"{key} download result: {result}")

    # maybe this check should be out of the function?
    if all([(FRAGMENT_DIR / file.download_file_name).exists() for file in files_to_download]):
        logger.info("All raw files local, processing fragments now")
    else:
        logger.error("Some raw files not found locally, please rerun")
        logger_thread_shutdown_and_exit()


def filter_worker(fragment_meta: FragmentFileMeta) -> FragmentWorkerResult:
    """
    Worker function to filter raw fragment files.
    Each process needs to initialize logging to report to queue and subsequent
    handling by the logger listener thread
    """
    create_logger(fragment_meta.queues.logging_queue)
    #read in the fragments
    logging.info(f"Starting filtering of {fragment_meta.download_file_name}...")
    if fragment_meta.cell_label_location == "suffix":
        a = f"{REPLACE_WITH}{fragment_meta.label}"
    else:
        a = f"{fragment_meta.label}{REPLACE_WITH}"

    logging.debug(f"{fragment_meta.accession} barcode replace with: {a}")
    logging.debug(f"{fragment_meta.accession} label: {fragment_meta.label}")
    logging.debug(f"{fragment_meta.accession} cell_label_location: {fragment_meta.cell_label_location}")

    file_path = FRAGMENT_DIR / fragment_meta.download_file_name
    frags_df = pd.read_csv(
        file_path,
        comment="#",
        sep="\t",
        names=FRAGMENT_COL_NAMES
    )

    # TODO: figure out how to report QA stuff, initial attempt crashed comp
    #plot for QA
    counts = frags_df["barcode"].value_counts()
    # fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    # axes[0].hist(counts, range=(0,1000), bins=200)
    # axes[0].set_ylim(ymin=0)
    # axes[0].set_title("raw")

    #store stats for QA
    raw_min = counts.min()
    raw_mean = round(counts.mean())

    #update the barcode to match the CxG matrix obx index
    frags_df["barcode"] = frags_df["barcode"].apply(lambda x: re.sub(REPLACE_WITH, re.search(BARCODE_PATTERN, x).group(), a))

    #filter down to only barcodes in the CxG matrix
    frags_df = frags_df[frags_df["barcode"].isin(fragment_meta.barcodes)]

    #plot for QA
    counts = frags_df["barcode"].value_counts()
    # axes[1].hist(counts, range=(0,1000), bins=200)
    # axes[1].set_ylim(ymin=0)
    # axes[1].set_title("filtered")

    #store stats for QA
    stats = {
        "raw_matrix": fragment_meta.accession,
        "raw min": raw_min,
        "filt min": counts.min(),
        "raw mean": raw_mean,
        "filt mean": round(counts.mean()),
        "unique barcodes": len(counts)
    }
    figure = "test_string"

    #write the filtered fragments file
    output = fragment_meta.filtered_fragment_path_name
    frags_df.to_csv(
        output,
        sep="\t",
        index=False,
        header=False
    )
    file_saved = True if output.exists() else False

    logging.info(f"Finished filtering of {fragment_meta.download_file_name}. SUCCESS: {file_saved}")

    if file_saved:
        fragment_meta.queues.compression_queue.put(str(output))
        logging.debug(f"Put {output} on compression queue")

    return FragmentWorkerResult(
        accession=fragment_meta.accession,
        file_saved=file_saved,
        stats=stats,
        plot=figure,
        output_path=output,
    )


def duplicate_worker(fragment_meta: FragmentFileMeta) -> FragmentWorkerResult:
    """
    Worker to remove duplicates. Call after concatenator run and CXG validation finds duplicates
    in the final concatenated fragment file.
    Like filter worker, need to initialize logger to report to queue and listener thread
    """
    #read in the fragments
    create_logger(fragment_meta.queues.logging_queue)
    file_saved = False
    compressed_file_name = fragment_meta.filtered_fragment_path_name.name.replace("tsv", "tsv.gz")
    logging.info(f"Starting de-duplication of {compressed_file_name}...")

    file_path = FRAGMENT_DIR / compressed_file_name
    filtered_df = pd.read_csv(
        file_path,
        comment="#",
        sep="\t",
        names=FRAGMENT_COL_NAMES
    )

    duplicates = filtered_df[filtered_df.duplicated(keep=False)]
    filtered_df = filtered_df.drop_duplicates(keep="first")

    has_duplicates = True if not duplicates.empty else False
    logging.info(f"For {fragment_meta.accession}, found duplicates: {has_duplicates}")

    output = fragment_meta.filtered_fragment_path_name

    # only need to save file again if duplicates found
    if has_duplicates:
        filtered_df.to_csv(
            output,
            sep="\t",
            index=False,
            header=False
        )
        file_saved = output.exists()
        logging.info(f"Finished saving deduplicated {fragment_meta.filtered_fragment_path_name}. SUCCESS: {file_saved}")
        if file_saved:
            fragment_meta.queues.compression_queue.put(str(output))
            logging.debug(f"Put {output} on compression queue")


    return FragmentWorkerResult(
        accession=fragment_meta.accession,
        file_saved=file_saved,
        stats={},
        plot=None,
        output_path=output,
        checked_duplicates=True,
        has_duplicates=has_duplicates,
        duplicates=duplicates if not duplicates.empty else None
    )


def compress_worker(file_path: str | os.PathLike) -> None:
    logger.debug(f"Filtered file to compress: {file_path}")
    logger.info(f"Started gzipping {file_path}...")
    s = subprocess.run(["gzip", "-f", file_path])
    file_name = f"{s.args[-1].split('/')[-1]}.gz"
    logger.info(f"Finished gzipping {file_name}")


def compress_thread_target(queues: WorkerQueues, worker: Callable):
    logger.debug("Started compress thread target")
    while True:
        file_path = queues.compression_queue.get()
        if file_path is None:
            logger.debug("Put None on compression queue from thread target")
            queues.compression_queue.put(None)
            logger.debug("Stopped compress thread target")
            break
        worker(file_path)


def gzip_thread(queues: WorkerQueues, num_threads: int):
    logger.debug("Starting compression thread")

    # initialize threads
    threads = []
    for _ in range(num_threads):
        thread = threading.Thread(
            target=compress_thread_target,
            args=(queues, compress_worker),
        )
        thread.start()
        threads.append(thread)

    # wait for threads to finish
    for thread in threads:
        thread.join()

    logger.info("All compression completed")


def concat_files(filtered_files: list[str | os.PathLike]) -> None:
    """
    Concat all filtered, compressed files into final file. Will always run
    and compress all files associated with Processed Matrix accession
    """
    ind_frag_files_gz = [str(f) + '.gz' for f in filtered_files]
    for f in ind_frag_files_gz:
        logger.debug(f"File to add to final concatenated file: {f}")

    concat_frags = FRAGMENT_DIR / f"{args.file}_concatenated_filtered_fragments.tsv.gz"
    logger.info(f"Concatenating {len(ind_frag_files_gz)} compressed filtered files into final file...")
    subprocess.run(["cat " + " ".join(ind_frag_files_gz) + " > " + str(concat_frags)], shell=True)
    logger.info(f"Final file saved as {concat_frags.parent / concat_frags.name}")


def run_processing_pool(worker_function: Callable, fragment_meta: list[FragmentFileMeta]) -> list[FragmentWorkerResult]:
    """
    Run process pool for fragment processing.
    Currently can run pool to filter fragments or remove duplicates.
    Can add lift over worker function at some point if that is needed
    """
    results = list()
    with multiprocessing.Pool(processes=NUM_PROCESS_WORKERS) as pool:
        iterator = pool.imap(worker_function, fragment_meta)
        for meta in iterator:
            try:
                results.append(meta)
            except StopIteration:
                break
            except Exception as e:
                logger.error(f"Worker/ProcessPool error: {e}")
                results.append((meta, e))

    return results


def print_results(results: list[FragmentWorkerResult]) -> None:
    """
    Function to display worker results
    Could probably be cleaner with logic, works for both filtering
    and de-duplication results
    """
    logger.info("Filter results")
    logger.info("=" * PRINT_WIDTH)
    for meta in results:
        if isinstance(meta, FragmentWorkerResult):
            file = meta.output_path.name + ".gz" if meta.checked_duplicates else meta.output_path.name
            logger.info(f"Fragment file: {file}")
            logger.info(f"Filtered file saved: {meta.file_saved}")
            if meta.checked_duplicates:
                logger.info(f"Duplicates results for {meta.accession}: {meta.has_duplicates}")
                if meta.duplicates is not None:
                    logger.info(meta.duplicates)
        else:
            logger.error("FAILURE")
            traceback.print_exception(None, meta, meta.__traceback__)
        logger.info("=" * PRINT_WIDTH)

    filtered = all(not meta.checked_duplicates for meta in results)
    if filtered:
        stats_df = pd.DataFrame([item.stats for item in results])
        logger.info(stats_df)
        logger.info(f"Fragment unique barcodes: {stats_df['unique barcodes'].sum()}")
        logger.info(f"AnnData unique barcodes: {fragment_meta_list[0].barcodes.shape[0]}")


if __name__ == "__main__":
    # might be able to use just queue and not multiprocesser manager
    # even with thread listener, manager queue best to handle process inheritance
    with Manager() as manager:
        args = getArgs()
        connection = Connection(args.mode)

        # set up queues and logging listener thread
        queues = WorkerQueues(
            logging_queue = manager.Queue(),
            compression_queue = manager.Queue(),
        )

        with open("log_config_concatenator.yaml", "rt") as f:
            logging_config = yaml.safe_load(f.read())

        log_file_name = f"{args.file}_outfile_concatenator.log"
        logging_config["handlers"]["file"]["filename"] = log_file_name
        logging.config.dictConfig(logging_config)

        listener = threading.Thread(target=logger_thread, args=(queues.logging_queue,))
        listener.start()

        # set up concatenator logger and first log messages
        logger = logging.getLogger("concatenator")
        logger.debug("Logger thread started")
        logger.debug(f"Command line args: {' '.join(sys.argv)}")
        logger.debug(f"Running concatenator on env: {args.mode}")

        fragment_meta_list = query_lattice(args.file, connection, queues)
        download_fragment_files(fragment_meta_list)

        worker_function = filter_worker
        if args.deduplicate:
            if all([meta.is_filtered_file_local for meta in fragment_meta_list]):
                logger.info("Found filtered fragment files for all raw matrices, starting duplicate check...")
                worker_function = duplicate_worker
            else:
                logger.error("Rerun concatenator to generate all filtered fragment files")
                missing_files = [
                    meta.filtered_fragment_path_name.name + ".gz" 
                    for meta in fragment_meta_list 
                    if not meta.is_filtered_file_local
                ]
                logger.error(f"Missing following filtered files, {len(missing_files)} total:")
                for file in missing_files:
                    logger.error(file)
                logger_thread_shutdown_and_exit()

        total_fragment_files = len(fragment_meta_list)
        compression_thread = threading.Thread(
            target=gzip_thread, 
            args=(queues, total_fragment_files,)
        )
        compression_thread.start()

        logger.debug(f"Worker pool function: {worker_function.__name__}()")
        results = run_processing_pool(worker_function, fragment_meta_list)
        print_results(results)

        queues.compression_queue.put(None)
        logger.debug("Put None on compression queue from main thread")
        compression_thread.join()
        logger.debug("Compression thread shutdown")

        # always concat all compressed filtered files into final file
        compressed_files = [file.output_path.absolute() for file in results]
        concat_files(compressed_files)

        # shut down logger thread
        logger.debug("Logger thread shutdown")
        queues.logging_queue.put(None)
        listener.join()
