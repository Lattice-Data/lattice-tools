import argparse
import traceback
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
BARCODE_PATTERN = r"[ACGT]{16}"
REPLACE_WITH = "B@RCODE"
FRAGMENT_DIR = Path("atac_fragments/")
FS = fsspec.filesystem("s3")
DOWNLOAD_THREADS = 8
NUM_FILTER_WORKERS = os.cpu_count() // 2
PROCESSED_MATRIX_ACCESSION = "LATDF393MGJ"
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
        self.is_file_local = (FRAGMENT_DIR / self.download_file_name).is_file()


@dataclass
class FragmentFilterResult:
    """
    Dataclass for returning filter results, plot currently unused
    """
    accession: str
    success: bool
    stats: dict
    plot: plt.figure
    output_path: Path


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
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def download_object(s3_client, fragment_meta: FragmentFileMeta):
    download_path = Path(FRAGMENT_DIR) / fragment_meta.download_file_name
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
        print(f"{num_all_files} fragment files in Lattice")
        print(f"Found {num_locally} files locally, downloading {len(future_to_key)} files")

        if not future_to_key:
            print("All files local, no downloading needed")

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
        print("All files local, filtereing fragments now")
    else:
        print("Some files not found locally, please rerun")
        sys.exit()


def filter_worker(fragment_meta: FragmentFileMeta) -> FragmentFilterResult:
    #read in the fragments
    print(f"Starting filtering of {fragment_meta.download_file_name}...")
    if fragment_meta.cell_label_location == "suffix":
        a = f"{REPLACE_WITH}-1{fragment_meta.label}"
    else:
        a = f"{fragment_meta.label}{REPLACE_WITH}-1"

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
    output = file_path.parent / file_path.name.replace("fragments.tsv.gz", "filtered_fragments.tsv")
    frags_df.to_csv(
        output,
        sep="\t",
        index=False,
        header=False
    )
    success = True if output.exists() else False
    print(f"Finished filtering of {fragment_meta.download_file_name}. SUCCESS: {success}")

    return FragmentFilterResult(
        accession=fragment_meta.accession,
        success=success,
        stats=stats,
        plot=figure,
        output_path=output,
    )


def compress_and_concat(filtered_files: list[str | os.PathLike]) -> None:
    processes = []
    print("Starting gzip compression of filtered files...")
    for f in filtered_files:
        p = subprocess.Popen(["gzip",f])
        processes.append(p)

    for p in processes:
        p.wait()

    ind_frag_files_gz = [str(f) + '.gz' for f in filtered_files]
    concat_frags = FRAGMENT_DIR / f"{args.file}_concatenated_filtered_fragments.tsv.gz"
    print("Concatenating final file...")
    subprocess.run(["cat " + " ".join(ind_frag_files_gz) + " > " + str(concat_frags)], shell=True)


if __name__ == "__main__":
    args = getArgs()
    connection = Connection(args.mode)

    fragment_meta = query_lattice(args.file, connection)
    download_fragment_files(fragment_meta)

    results = list()
    with multiprocessing.Pool(processes=NUM_FILTER_WORKERS) as pool:
        iterator = pool.imap(filter_worker, fragment_meta)
        for meta in iterator:
            try:
                results.append(meta)
            except StopIteration:
                break
            except Exception as e:
                print(f"ERROR: {e}")
                results.append((meta, e))

    print("Filter results")
    for meta in results:
        if isinstance(meta, FragmentFilterResult):
            print(f"Fragment file: {meta.output_path.name}")
            print(f"Filtered file saved: {meta.success}")
        else:
            print("FAILURE")
            traceback.print_exception(None, meta, meta.__traceback__)
        print("=" * 40)

    filtered_files = [file.output_path.absolute() for file in results if file.success]
    compress_and_concat(filtered_files)
