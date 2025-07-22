import argparse
import traceback
from typing import Callable
import boto3.session
import fsspec
import h5py
import matplotlib.pyplot as plt
import multiprocessing
import os
import pandas as pd
import re
import subprocess
import sys

from anndata._io.specs import read_elem
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

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
NUM_FILTER_WORKERS = os.cpu_count() // 2
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
class FragmentFileMeta:
    """
    Dataclass to pass fragment file meta to workers for download and processing
    """

    cell_label_location: str
    label: str
    accession: str
    uri: URIPath
    barcodes: pd.Series

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
    Dataclass for returning filter results, plot currently unused
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


def download_object(s3_client, fragment_meta: FragmentFileMeta):
    download_path = FRAGMENT_DIR / fragment_meta.download_file_name
    print(f"Downloading {fragment_meta.download_file_name} to {download_path}")
    s3_client.download_file(
        fragment_meta.uri.bucket_name,
        fragment_meta.uri.file_path,
        str(download_path)
    )
    return "Success"


def download_parallel_multithreading(files_to_download: list[FragmentFileMeta]):
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
        print(f"{num_all_files} raw fragment files in Lattice")
        print(f"Found raw {num_locally} files locally, downloading {len(future_to_key)} files")

        if not future_to_key:
            print("All raw files local, no downloading needed")

        for future in futures.as_completed(future_to_key):
            key = future_to_key[future]
            exception = future.exception()

            if not exception:
                yield key, future.result()
            else:
                yield key, exception


def query_lattice(processed_matrix_accession: str, connection: Connection) -> list[FragmentFileMeta]:
    processed_matrix_report = get_report(
        "ProcessedMatrixFile",
        f"&@id=/processed-matrix-files/{processed_matrix_accession}/",
        PROCESSED_MATRIX_FIELD_LIST,
        connection
    )[0]

    h5ad_uri = URIPath(processed_matrix_report["s3_uri"])
    with h5py.File(FS.open(h5ad_uri.full_uri)) as f:
        barcodes = read_elem(f["obs"]).index.to_series()

    master_fragment_file_meta = []
    for raw_matrix_meta in processed_matrix_report["cell_label_mappings"]:
        obj_type, filter_url = parse_ids([raw_matrix_meta["raw_matrix"]])
        raw_matrix_report = get_report(
            obj_type,
            filter_url,
            field_lst=RAW_MATRIX_FIELD_LIST,
            connection=connection
        )[0]

        # TODO: current workaround until fragment uris placed in Lattice
        fragment_uri = raw_matrix_report["s3_uri"].replace("filtered_feature_bc_matrix.h5", "atac_fragments.tsv.gz")
        accession = raw_matrix_report["accession"]
        raw_matrix_meta["fragment_file_s3_uri"] = fragment_uri
        raw_matrix_meta["accession"] = raw_matrix_report["accession"]

        assert FS.isfile(raw_matrix_meta["fragment_file_s3_uri"]), f"raw matrix{accession} does not have fragment file s3"
        master_fragment_file_meta.append(
            FragmentFileMeta(
                cell_label_location=processed_matrix_report["cell_label_location"],
                label=raw_matrix_meta["label"],
                accession=accession,
                uri=URIPath(fragment_uri),
                barcodes=barcodes
            )
        )

    return master_fragment_file_meta


def download_fragment_files(files_to_download: list[FragmentFileMeta]) -> None:
    for key, result in download_parallel_multithreading(files_to_download):
        print(f"{key} download result: {result}")

    # maybe this check should be out of the function?
    if all([(FRAGMENT_DIR / file.download_file_name).exists() for file in files_to_download]):
        print("All raw files local, processing fragments now")
    else:
        print("Some raw files not found locally, please rerun")
        sys.exit()


def filter_worker(fragment_meta: FragmentFileMeta) -> FragmentWorkerResult:
    #read in the fragments
    print(f"Starting filtering of {fragment_meta.download_file_name}...")
    if fragment_meta.cell_label_location == "suffix":
        a = f"{REPLACE_WITH}{fragment_meta.label}"
    else:
        a = f"{fragment_meta.label}{REPLACE_WITH}"

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
    print(f"Finished filtering of {fragment_meta.download_file_name}. SUCCESS: {file_saved}")

    return FragmentWorkerResult(
        accession=fragment_meta.accession,
        file_saved=file_saved,
        stats=stats,
        plot=figure,
        output_path=output,
    )


def duplicate_worker(fragment_meta: FragmentFileMeta) -> FragmentWorkerResult:
    #read in the fragments
    file_saved = False
    compressed_file_name = fragment_meta.filtered_fragment_path_name.name.replace("tsv", "tsv.gz")
    print(f"Starting de-duplication of {compressed_file_name}...")

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
    print(f"For {fragment_meta.accession}, found duplicates: {has_duplicates}")

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
        print(f"Finished saving deduplicated {fragment_meta.filtered_fragment_path_name}. SUCCESS: {file_saved}")

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


def compress_files(filtered_files: list[str | os.PathLike]) -> None:
    if not filtered_files:
        print("All filtered files already compressed")
        return

    processes = []
    print(f"Starting gzip compression of {len(filtered_files)} filtered files...")
    for f in filtered_files:
        p = subprocess.Popen(["gzip","-f",f])
        processes.append(p)

    for p in processes:
        p.wait()


def concat_files(filtered_files: list[str | os.PathLike]) -> None:
    ind_frag_files_gz = [str(f) + '.gz' for f in filtered_files]
    concat_frags = FRAGMENT_DIR / f"{args.file}_concatenated_filtered_fragments.tsv.gz"
    print(f"Concatenating {len(ind_frag_files_gz)} compressed filtered files into final file...")
    subprocess.run(["cat " + " ".join(ind_frag_files_gz) + " > " + str(concat_frags)], shell=True)


def run_processing_pool(worker_function: Callable, fragment_meta: list[FragmentFileMeta]) -> list[FragmentWorkerResult]:
    results = list()
    with multiprocessing.Pool(processes=NUM_FILTER_WORKERS) as pool:
        iterator = pool.imap(worker_function, fragment_meta)
        for meta in iterator:
            try:
                results.append(meta)
            except StopIteration:
                break
            except Exception as e:
                print(f"ERROR: {e}")
                results.append((meta, e))

    return results


def print_results(results: list[FragmentWorkerResult]) -> None:
    print("Filter results")
    print("=" * PRINT_WIDTH)
    for meta in results:
        if isinstance(meta, FragmentWorkerResult):
            file = meta.output_path.name + ".gz" if meta.checked_duplicates else meta.output_path.name
            print(f"Fragment file: {file}")
            print(f"Filtered file saved: {meta.file_saved}")
            if meta.checked_duplicates:
                print(f"Duplicates results for {meta.accession}: {meta.has_duplicates}")
                if meta.duplicates is not None:
                    print(meta.duplicates)
        else:
            print("FAILURE")
            traceback.print_exception(None, meta, meta.__traceback__)
        print("=" * PRINT_WIDTH)

    filtered = all(not meta.checked_duplicates for meta in results)
    if filtered:
        stats_df = pd.DataFrame([item.stats for item in results])
        print(stats_df)
        print(f"Fragment unique barcodes: {stats_df['unique barcodes'].sum()}")
        print(f"AnnData unique barcodes: {fragment_meta[0].barcodes.shape[0]}")


if __name__ == "__main__":
    args = getArgs()
    connection = Connection(args.mode)

    fragment_meta = query_lattice(args.file, connection)
    download_fragment_files(fragment_meta)

    worker_function = filter_worker
    if args.deduplicate:
        if all([meta.is_filtered_file_local for meta in fragment_meta]):
            print("Found filtered fragment files for all raw matrices, starting duplicate check...")
            worker_function = duplicate_worker
        else:
            print("Rerun concatenator to generate all filtered fragment files")
            missing_files = [
                meta.filtered_fragment_path_name.name + ".gz" 
                for meta in fragment_meta if not meta.is_filtered_file_local
            ]
            print(f"Missing following filtered files, {len(missing_files)} total:")
            for file in missing_files:
                print(file)
            sys.exit()

    results = run_processing_pool(worker_function, fragment_meta)
    print_results(results)

    # only need to gzip newly saved files
    non_compressed_files = [file.output_path.absolute() for file in results if file.file_saved]
    compress_files(non_compressed_files)

    compressed_files = [file.output_path.absolute() for file in results]
    concat_files(compressed_files)
