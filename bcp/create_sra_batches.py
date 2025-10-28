import argparse
import csv
import fsspec
import pandas as pd
import requests
import sys

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote


EPILOG = """
Script to break SRA uploads into smaller batches below the 5 TB limit

Inputs:
    txt/csv of Fastq S3 URIs for upload
    Wrangling Sheet ID

Outputs:
    txts/csvs of batches for upload
    updated csv of SRA_metadata sheet for SRA submission
"""

FS = fsspec.filesystem('s3')
BATCH_LIMIT = 5_000_000_000_000
FILE_LIMIT = 100_000_000_000


@dataclass
class FastqMeta:
    """
    Parsed metadata from S3 directory/filename structure until Lattice v2 exists

    Works for first upload project, skip S3 info if using for predicted split file

    TODO:
        Validation
        Logic for dealing with key or full S3 URI; currently only key works
        Maybe base class and childs for other CZI bucket schema/SOP
        Mechanism for dealing/tracking paired files?
        Better operability with upload script
    """
    s3_key: str
    predicted_file: bool = False
    size: int = 0

    def __post_init__(self):
        # base meta for building other attributes
        self._path = Path(self.s3_key)

        # skip if file predicted and S3 meta does not exist
        if not self.predicted_file:
            self._s3_meta = FS.info(self.s3_key)
        
            # get S3 meta and set as attributes
            for key, value in self._s3_meta.items():
                setattr(self, key, value)

        # CZI raw file attributes
        self.file_name = self._path.name
        self.provider = self._path.parts[0].split("-")[-1]
        self.lastname = self._path.parts[1].split("-")[0]
        self.projectname = "-".join(self._path.parts[1].split("-")[1:])
        self.order_number = self._path.parts[2]
        self.GroupID = self._path.parts[3]
        self.sample_name = self._path.parts[5].split("-")[1]
        self.library_ID = self.sample_name + "_lib"
        self.RunID = self._path.parts[5].split("-")[0]
        self.Assay = self.sample_name.split("_")[-1]
        self.UG_BC = "-".join(self._path.parts[5].split("_")[2].split("-")[1:])
        self.SampleNumber = self._path.parts[5].split("_")[3]
        self.LaneNumber = self._path.parts[5].split("_")[4]
        self.ReadNumber = self._path.parts[5].split("_")[5]
        self.OptionalPartNumber = self._path.parts[5].split("_")[6].split(".")[0]

    def display_meta(self):
        """
        Print out display of S3 and CZI attributes
        """
        for attr_name, attr_value in self.__dict__.items():
            if not attr_name.startswith("_"):
                # better handling of large nums to see byte size
                is_int = isinstance(attr_value, int) and not isinstance(attr_value, bool)
                attr_formatted = f"{attr_value:_}" if is_int else f"{attr_value}"
                print(f"{attr_name}: {attr_formatted}")
        print()


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--s3uris',
        '-s',
        dest='s3uris',
        help='File with list of S3 URIs'
    )
    parser.add_argument(
        '--googlesheet',
        '-g',
        dest='googlesheet',
        help='google sheet ID of wrangling sheet'
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def get_wrangling_sheet(sheet_id: str, tab_name: str = "SRA_metadata") -> pd.DataFrame:
    """Get SRA_metadat wrangling sheet"""
    url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={quote(tab_name)}'
    sra_meta = pd.read_csv(url)
    return sra_meta


def create_split_file_meta(fastq_keys: list[str]) -> list[FastqMeta]:
    """
    Create list of FastqMeta that takes into account proper splitting.
    This should mirror the upload script's determination of splitting and
    file naming.

    """
    fastq_metas = []

    for key in fastq_keys:
        meta = FastqMeta(key)
        if meta.size > FILE_LIMIT:
            if meta.ReadNumber == "R1":
                # will add R1 during R2 add
                continue
            read1_key = key.replace("R2_001.fastq.gz", "R1_001.fastq.gz")
            read2_key = key
            meta1 = FastqMeta(read1_key)
            meta2 = meta
            # looks like trick to cause floor or ceiling rounding
            split_number = -(-meta2.size // FILE_LIMIT)
            file1_size = meta1.size // split_number
            file2_size = meta2.size // split_number
            stop = split_number + 1
            for num in range(1, stop):
                fastq_metas.append(
                    FastqMeta(
                        read1_key.replace(".fastq.gz", f"_split{num}.fastq.gz"),
                        predicted_file=True,
                        size=file1_size,
                    )
                )
                fastq_metas.append(
                    FastqMeta(
                        read2_key.replace(".fastq.gz", f"_split{num}.fastq.gz"),
                        predicted_file=True,
                        size=file2_size,
                    )
                )
            continue

        fastq_metas.append(FastqMeta(key))

    return fastq_metas


def create_fastq_batches(fastq_metas: list[FastqMeta]) -> dict[int, list[FastqMeta]]:
    """
    After generation of the FastqMeta for all files and/or their splits,
    this function calculates the batches.

    Will likely need to call twice:
        With generation of predicted split files for wrangling
        Without predicted splits for feeding into uploading script

    Using GroupID for the moment to keep GEX and CRI fastqs together
    """
    batch_num = 1
    current_size_sum = 0
    sample_dict = defaultdict(list)
    batches = defaultdict(list)

    for meta in fastq_metas:
        sample_dict[meta.GroupID].append(meta)

    for group_id, metas in sample_dict.items():
        group_size = sum([meta.size for meta in metas])
        current_size_sum += group_size
        if current_size_sum > BATCH_LIMIT:
            # reset size and increase to next batch
            current_size_sum = group_size
            batch_num += 1
        for meta in metas:
            batches[batch_num].append(meta)

    return batches


def create_final_sra_meta_df(fastq_split_meta: list[FastqMeta], sra_meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create the final df for the wrangling sheet that takes into account the predicted split files
    """
    new_file_dict = defaultdict(list)

    for meta in fastq_split_meta:
        new_file_dict[meta.sample_name].append(meta.file_name)

    df_list = []
    for sample_name, files in new_file_dict.items():
        temp_dict = {"sample_name": sample_name}
        for file_num, file in enumerate(files, start=1):
            if file_num == 1:
                temp_dict["filename"] = file
                continue
            temp_dict[f"filename{file_num}"] = file
        df_list.append(temp_dict)

    file_df = pd.DataFrame(df_list)
    stripped_sra_meta = sra_meta_df.loc[ : , "sample_name": "filetype"]
    final_sra_meta_df = stripped_sra_meta.merge(file_df, on="sample_name", how="left")
    return final_sra_meta_df

 
if __name__ == "__main__":
    args = getArgs()
    sra_meta_df = get_wrangling_sheet(args.googlesheet)
    s3_uri_path = Path(args.s3uris)

    with open(s3_uri_path, "r") as f:
        s3_uris = [line.lstrip("s3://").strip() for line in f]
    
    fastq_metas = [FastqMeta(uri) for uri in s3_uris]
    print("Found base file fastq metadata")

    split_fastq_metas = create_split_file_meta(s3_uris)
    print("Created split file fastq metadata")

    uploader_batches = create_fastq_batches(fastq_metas)
    sra_meta_batches = create_fastq_batches(split_fastq_metas)
    print("Created batches")

    final_sra_meta_df = create_final_sra_meta_df(split_fastq_metas, sra_meta_df)
    final_sra_meta_df.to_csv("final_sra_metadata.csv")
    print("Saved updated SRA_metadata as 'final_sra_metadata.csv'")

    last_name = fastq_metas[0].lastname
    project_name = fastq_metas[0].projectname

    print("Writing to batch files")
    for batch, metas in uploader_batches.items():
        # does this need to be csv? could just be txt with one item per row
        with open(f"{last_name}_{project_name}_batch{batch}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            for meta in metas:
                writer.writerow([meta.name])
