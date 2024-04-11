import argparse
import boto3
import logging
import os
import scanpy as sc
import subprocess
import sys

from botocore.exceptions import ClientError
from dataclasses import dataclass


EPILOG = """
This script will take a list of s3 uri of uncompressed h5ad files as input, and it will create a compressed h5ad
in the same S3 directory, with '_curated' added to file name.

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
    print(uri_info.file_name + " downloading")

    try:
        uri_info.metadata = S3_CLIENT.get_object_attributes(
            Bucket=uri_info.bucket_name, 
            Key=uri_info.file_path,
            ObjectAttributes=["ObjectSize"]
        )
        print(f"Last Modified: {uri_info.metadata['LastModified']}")
        print(f"Object Size: {uri_info.metadata['ObjectSize']}")
        S3_CLIENT.download_file(
            uri_info.bucket_name, uri_info.file_path, TEMP_DIR + "/" + uri_info.file_name
        )
    except subprocess.CalledProcessError as e:
        logging.error("ERROR: {} Failed to find file {} on s3".format(e, uri_info.full_uri))
        sys.exit("ERROR: {} Failed to find file {} on s3".format(e, uri_info.full_uri))
    else:
        print(uri_info.file_name + " downloaded")


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
    local_exists([TEMP_DIR]) -> argument = TEMP_DIR
    local_exists(TEMP_DIR) -> argument = T/E/M/P/_D/I/R
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
    if not local_exists([TEMP_DIR]):
        os.mkdir(TEMP_DIR)

    with open(s3_uri_file, "r") as f:
        s3_uris = [URIMetaInfo(line.strip()) for line in f]

    for uri in s3_uris:
        # download
        download_file(uri)

        # compress
        compress_h5ad(uri)

        # upload
        print(f"Uploading compressed h5ad to this object key: {uri.file_path}")
        with open(os.path.join(TEMP_DIR, uri.new_file_name), "rb") as f:
            S3_CLIENT.upload_fileobj(f, uri.bucket_name, uri.file_path)

        # remove h5ads
        if local_exists([TEMP_DIR, uri.file_name]):
            os.remove(os.path.join(TEMP_DIR, uri.file_name))
            os.remove(os.path.join(TEMP_DIR, uri.new_file_name))
            print(f"Removed files {uri.file_name} and {uri.new_file_name} from {TEMP_DIR}")
        else:
            print(f"File {uri.file_name} not found in {TEMP_DIR}")


args = getArgs()
if __name__ == "__main__":
    main(args.file)
