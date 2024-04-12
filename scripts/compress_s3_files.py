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
This script will take a list of s3 uris of uncompressed h5ad files as input, it will create a compressed h5ad
in a temp directory, and then upload and replace the original s3 h5ad with the compressed version. If the 
compressed file is the same size or smaller, the orginal file will be left in place.

Examples:

    python %(prog)s -f list_of_uri.txt

For more details:

    python %(prog)s --help
"""

S3_CLIENT = boto3.client("s3")
TEMP_DIR = "temp_dir"


@dataclass
class URIMetaInfo:
    """
    Dataclass to hold various parsed information for each S3 uri
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
        self.new_file_name = self.file_name.replace(".h5ad", "_curated.h5ad")
        self.new_file_path = self.file_path.replace(".h5ad", "_curated.h5ad")


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--file", "-f", help="A text list of s3 uri of files to be compressed"
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def download_file(uri_info: URIMetaInfo):
    """
    Given an S3 uri, the file will be downloaded locally

    :param uri_info: URIMetaInfo object, uses proper attributes below
    :return: None, downloads file
    """
    print(uri_info.file_name + " downloading ...")

    try:
        S3_CLIENT.download_file(
            uri_info.bucket_name, uri_info.file_path, TEMP_DIR + "/" + uri_info.file_name
        )
    except subprocess.CalledProcessError as e:
        logging.error("ERROR: {} Failed to find file {} on s3".format(e, uri_info.full_uri))
        sys.exit("ERROR: {} Failed to find file {} on s3".format(e, uri_info.full_uri))
    else:
        print(uri_info.file_name + " downloaded")


def get_file_metadata(uri_info: URIMetaInfo):
    """
    Gets s3 object metadata and adds to URIMetaInfo.metadata attribute

    :param uri_info: URIMetaInfo object, uses proper attributes below
    :return: None, object attribute updated
    """
    try:
        uri_info.metadata = S3_CLIENT.get_object_attributes(
            Bucket=uri_info.bucket_name, 
            Key=uri_info.file_path,
            ObjectAttributes=["ObjectSize"]
        )
    except S3_CLIENT.exceptions.NoSuchKey as e:
        logging.error(f"{e}: Object key {uri_info.file_path} does not exist")
        sys.exit(f"{e}: Object key {uri_info.file_path} does not exist")
    else:
        print(f"Processing file {uri_info.file_name} ...")
        print(f"Last Modified: {uri_info.metadata['LastModified']}")
        print(f"Object Size: {uri_info.metadata['ObjectSize']}")


def compress_h5ad(h5ad: URIMetaInfo):
    """
    Takes downloaded h5ad and saves gzipped h5ad

    :param h5ad: URIMetaInfo, uses attributes to look for local h5ad
    :return: None, creates compressed h5ad in TEMP_DIR
    """
    adata = sc.read_h5ad(os.path.join(TEMP_DIR, h5ad.file_name))
    sc.write(os.path.join(TEMP_DIR, h5ad.new_file_name), adata, compression="gzip")
    print(f"Compressed {h5ad.new_file_name} to {TEMP_DIR}")


def local_exists(*files: list[str]) -> bool:
    """
    Helper function to try to clarify logic in main() for existence of local
    files/directories

    :param *files: list of strings of local directories and file
    :return: True if directory/file exists, else False

    Make sure to use list with single item to check directory:
    local_exists([TEMP_DIR]) -> argument will be TEMP_DIR
    local_exists([TEMP_DIR, "file.h5ad"]) -> argument will be TEMP_DIR/file.h5ad
    local_exists(TEMP_DIR) -> argument will be T/E/M/P/_D/I/R
    """
    return os.path.exists("/".join(*files))


def s3_exists(uri_info: URIMetaInfo) -> bool:
    """
    Checks if file exists on S3

    :param uri_info: URIMetaInfo object, uses proper attributes below
    :return: True if object key exists on S3 bucket, else False
    """
    try:
        S3_CLIENT.get_object(
            Bucket=uri_info.bucket_name,
            Key=uri_info.new_file_path,
        )
        return True
    except S3_CLIENT.exceptions.NoSuchKey:
        return False


def main(s3_uri_file):
    log_file = f"{s3_uri_file}_outfile_s3_compress.log"
    logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)
    time_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    logging.info("Date and time of s3 compression run: " + time_date)
    logging.captureWarnings(True)

    files_not_changed = []

    if not local_exists([TEMP_DIR]):
        os.mkdir(TEMP_DIR)

    with open(s3_uri_file, "r") as f:
        s3_uris = [URIMetaInfo(line.strip()) for line in f]

    total_files = len(s3_uris)

    for uri in s3_uris:
        print(f"Remaining files to process: {total_files}")
        # get object metadata
        get_file_metadata(uri)

        # download
        download_file(uri)

        # compress
        compress_h5ad(uri)

        # upload
        compressed_size = os.path.getsize(os.path.join(TEMP_DIR, uri.new_file_name))
        original_size = uri.metadata["ObjectSize"]

        if original_size > compressed_size:
            print(f"Uploading compressed h5ad to this object key: {uri.file_path}")
            with open(os.path.join(TEMP_DIR, uri.new_file_name), "rb") as f:
                S3_CLIENT.upload_fileobj(f, uri.bucket_name, uri.file_path)
        else:
            print(f"INFO: Original file size {original_size} <= {compressed_size}, not uploading to S3")
            files_not_changed.append(uri.full_uri)

        # remove h5ads
        if local_exists([TEMP_DIR, uri.file_name]):
            os.remove(os.path.join(TEMP_DIR, uri.file_name))
            os.remove(os.path.join(TEMP_DIR, uri.new_file_name))
            print(f"Removed files {uri.file_name} and {uri.new_file_name} from {TEMP_DIR}")
        else:
            print(f"File {uri.file_name} not found in {TEMP_DIR}")

        total_files -= 1
        print("=====================================")

    with open(log_file, 'a') as f:
        f.write("Files not uploaded due to original size <= compressed size :\n")
        for uri in files_not_changed:
            f.write(f"{uri}\n")


args = getArgs()
if __name__ == "__main__":
    main(args.file)
