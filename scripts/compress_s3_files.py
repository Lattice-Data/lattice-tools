import argparse
import boto3
import logging
import os
import scanpy as sc
import subprocess
import sys

from dataclasses import dataclass
from datetime import datetime


EPILOG = """
This script will take a list of S3 URIs of uncompressed h5ad files in a txt file as input, it will create a compressed h5ad
in a temp directory, and then upload and replace the original S3 h5ad with the compressed version. If the 
compressed file is the same size or larger, the orginal file will be left in place.

Examples:

    python %(prog)s -f list_of_uri.txt

For more details:

    python %(prog)s --help
"""

S3_CLIENT = boto3.client("s3")
TEMP_DIR = "temp_dir"

# initial logging setup to then use within functions
# console print and logging all handled through logging setup,
# no repetition of print and logging statements
# debug level will print to console for updates while running
# info level to console and log to record most important info
logging.captureWarnings(capture=True)
logger = logging.getLogger('compress_s3_files')

@dataclass
class URIMetaInfo:
    """
    Dataclass to hold various parsed information for each S3 URI
    Only needs full uri to create the various attributes
    Can add __repr__ and __str__ methods for print debugging or other
    usefulness when working with these objects
    """

    full_uri: str
    metadata = None

    def __post_init__(self):
        self.bucket_name: str = self.full_uri.split("/")[2]
        self.file_path = self.full_uri.replace("s3://{}/".format(self.bucket_name), "")
        self.file_name = self.full_uri.split("/")[-1]
        self.new_file_name = self.file_name.replace(".h5ad", "_compressed.h5ad")
        self.new_file_path = self.file_path.replace(".h5ad", "_compressed.h5ad")


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--file", "-f", help="A text list of S3 URIs of files to be compressed"
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def download_file(uri_info: URIMetaInfo):
    """
    Given an S3 URI, the file will be downloaded locally

    :param uri_info: URIMetaInfo object, uses proper attributes below
    :return: None, downloads file
    """
    logger.debug(uri_info.file_name + " downloading ...")

    S3_CLIENT.download_file(
        uri_info.bucket_name, 
        uri_info.file_path, 
        TEMP_DIR + "/" + uri_info.file_name
    )
    logger.debug(uri_info.file_name + " downloaded")


def get_file_metadata(uri_info: URIMetaInfo):
    """
    Gets S3 object metadata and adds to instance URIMetaInfo.metadata attribute

    :param uri_info: URIMetaInfo object, uses proper attributes below
    :return: None, object attribute updated
    """
    uri_info.metadata = S3_CLIENT.get_object_attributes(
        Bucket=uri_info.bucket_name, 
        Key=uri_info.file_path,
        ObjectAttributes=["ObjectSize"]
    )
    logger.info(f"Processing file {uri_info.file_name} ...")
    logger.info(f"Last Modified: {uri_info.metadata['LastModified']}")
    logger.info(f"Object Size: {uri_info.metadata['ObjectSize']:,} bytes")


def compress_h5ad(h5ad: URIMetaInfo):
    """
    Takes downloaded h5ad and saves gzipped h5ad

    :param h5ad: URIMetaInfo, uses attributes to look for local h5ad
    :return: None, creates compressed h5ad in TEMP_DIR
    """
    adata = sc.read_h5ad(os.path.join(TEMP_DIR, h5ad.file_name))
    sc.write(os.path.join(TEMP_DIR, h5ad.new_file_name), adata, compression="gzip")
    logger.debug(f"Compressed {h5ad.new_file_name} to {TEMP_DIR}")


def log_files_lists(s3_uris, index, files_changed, files_not_changed, log_file):
    """
    Appends to log file the S3 URIs not changed in place and the URIs not processed
    The URIs are appended one per line so it is easy to copy into another file to rerun this script
    or use elsewhere 

    :param s3_uris: list of URIMetaInfo objects generated from input txt
    :param index: int, use index from enumerate in main loop
    :param files_changed: list of strs of the full URI for files uploaded to S3
    :param files_not_changed: list of strs of the full URI for files not changed on S3
    :param log_file: str, name of log file for this run
    :return: None, appends log file created for current run of script
    """
    remaining_uris = [uri.full_uri for uri in s3_uris[index:]]
    logging_sections = {
        "Files uploaded to S3:\n": files_changed,
        "Files NOT uploaded due to original size <= compressed size:\n": files_not_changed,
        "Files NOT processed due to error and/or script exit:\n" : remaining_uris
    }
    with open(log_file, 'a') as f:
        for title, uris in logging_sections.items():
            f.write(title)
            for uri in uris:
                f.write(f"{uri}\n")
            f.write("\n")


def main(s3_uri_file):
    # set up logging
    log_file = f"{s3_uri_file}_outfile_s3_compress.log"
    time_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    file_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)-7s - %(message)s")
    console_formatter = logging.Formatter("%(name)s - %(message)s")
    ch.setFormatter(console_formatter)
    fh.setFormatter(file_formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.info("Date and time of s3 compression run: " + time_date)

    if not os.path.exists(TEMP_DIR):
        os.mkdir(TEMP_DIR)

    with open(s3_uri_file, "r") as f:
        s3_uris = [URIMetaInfo(line.strip()) for line in f]

    number_of_files = len(s3_uris)
    files_not_changed = []
    files_changed = []

    for index, uri in enumerate(s3_uris):
        logger.debug(f"Remaining files to process: {number_of_files}")

        # get object metadata
        try:
            get_file_metadata(uri)
        except S3_CLIENT.exceptions.NoSuchKey as e:
            logger.error(f"{e}: Object key {uri.file_path} does not exist")
            log_files_lists(s3_uris, index, files_changed, files_not_changed, log_file)
            sys.exit(f"{e}: Exiting script")

        # download
        try:
            download_file(uri)
        except subprocess.CalledProcessError as e:
            logger.error("ERROR: {} Failed to find file {} on S3".format(e, uri.full_uri))
            log_files_lists(s3_uris, index, files_changed, files_not_changed, log_file)
            sys.exit("Exiting script")

        # compress
        compress_h5ad(uri)

        compressed_size = os.path.getsize(os.path.join(TEMP_DIR, uri.new_file_name))
        original_size = uri.metadata["ObjectSize"]

        # only upload if compressed smaller than original
        if original_size > compressed_size:
            logger.info(f"Original size {original_size:,} > {compressed_size:,} compressed, uploading to S3")
            logger.debug(f"Uploading compressed h5ad to this object key: {uri.file_path}")
            with open(os.path.join(TEMP_DIR, uri.new_file_name), "rb") as f:
                S3_CLIENT.upload_fileobj(f, uri.bucket_name, uri.file_path)
            logger.debug(f"File {uri.file_name} uploaded to S3")
            files_changed.append(uri.full_uri)
        else:
            logger.info(f"Original size {original_size:,} <= {compressed_size:,} compressed, NOT uploading to S3")
            files_not_changed.append(uri.full_uri)

        # remove h5ads
        os.remove(os.path.join(TEMP_DIR, uri.file_name))
        os.remove(os.path.join(TEMP_DIR, uri.new_file_name))
        logger.debug(f"Removed files {uri.file_name} and {uri.new_file_name} from {TEMP_DIR}")

        number_of_files -= 1
        logger.info("===============================================")

    # final logging list after loop completion
    log_files_lists(s3_uris, len(s3_uris), files_changed, files_not_changed, log_file)


args = getArgs()
if __name__ == "__main__":
    main(args.file)
