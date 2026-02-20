import json
import os
import re
import sys
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from urllib.parse import ParseResult, urlparse
from urllib.request import Request, urlopen

import fsspec
import pandas as pd
import requests
from bs4 import BeautifulSoup

sys.path.append(os.path.dirname(os.path.abspath('../cellxgene_resources')))


@dataclass
class BCPBasePathMetadata:
    provider: str
    project_name: str
    order_number: str

    @staticmethod
    def pick_metadata_class(uri_path: URIPath) -> BCPBasePathMetadata:
        if "cellranger" in uri_path.parts:
            return TenXProcessedPathMeta

        if "processed" in uri_path.parts:
            return ScaleProcessedPathMeta
        
        if "raw" in uri_path.parts:
            if len(uri_path.parts) > 6:
                return ScaleRawPathMeta
            return TenXRawPathMeta

        return BCPBasePathMetadata

        
@dataclass
class ScalePathMeta(BCPBasePathMetadata):
    experiment_id: str
    is_raw: bool

    
@dataclass
class ScaleRawPathMeta(ScalePathMeta):
    run_id: str


@dataclass
class ScaleProcessedPathMeta(ScalePathMeta):
    run_date: str


@dataclass
class TenXPathMeta(BCPBasePathMetadata):
    group_id: str


@dataclass
class TenXRawPathMeta(TenXPathMeta):
    is_raw: bool


@dataclass
class TenXProcessedPathMeta(TenXPathMeta):
    is_raw: bool
    cellranger: str
    run_date: str


TEMP_DIR = Path("temp_cellranger/")
FS = fsspec.filesystem("s3")


@dataclass
class URIPath:
    """
    Dataclass to deal with S3 URIs and common metadata associated with them.
    full_uri should be in the format {s3://bucket/key/...}

    Since fsspec is used to get info, type, and size, if S3 URI is not valid,
    this will raise a FileNotFoundError exception
    """
    full_uri: str
    local_dir: Path = TEMP_DIR

    def __post_init__(self):
        self._parsed: ParseResult = urlparse(self.full_uri)
        self.info: dict = FS.info(self.full_uri)
        self.type: str = self.info["type"]
        self.size: int = self.info["size"]

    @property
    def bucket(self) -> str:
        return self._parsed.netloc
    
    @property
    def key(self) -> str:
        return self._parsed.path.strip("/")
    
    @property
    def parts(self) -> list[str]:
        """All S3 URI parts, including bucket and key"""
        return [self.bucket, *[part.strip("/") for part in self.key.split("/")]]   

    @property
    def file_name(self) -> str | None:
        if self.type == "file":
            return self.full_uri.split("/")[-1]
        return None

    @property
    def local_path(self) -> Path:
        if self.file_name:
            return self.local_dir / self.file_name
        return self.local_dir

    @property
    def bcp_path_metadata(self) -> BCPBasePathMetadata | None:
        """
        Trying dataclass inheritance to parse BCP metadata from paths.
        Uses trick to get list of input args as order and depth of S3 URI to match
        part of URI to corresponding metadata
        
        So far seems to work for 10x and Scale. Will return None if not a general BCP
        S3 URI
        """
        if (
            len(self.parts) > 2 and
            self.parts[0].startswith("czi-") and
            self.parts[2].startswith(("AN", "NV"))
        ):
            # pick class to generate standards input list
            unbound_class = BCPBasePathMetadata.pick_metadata_class(self)

            # pull args from correct class to get ordered standards input list
            standards_input = unbound_class.__dict__["__match_args__"]

            # traverse minimal depth between parts of URI and standards
            traversal_depth = min(
                len(self.parts),
                len(standards_input),
            )
            
            kwargs = {}
            for standards_level, uri_part in zip(standards_input[:traversal_depth], self.parts):
                match standards_level:
                    case "is_raw":
                        kwargs["is_raw"] = True if uri_part == "raw" else False
                    case str(standards_level):
                        kwargs[standards_level] = uri_part
            
            # construct instance with correct kwargs
            return unbound_class(**kwargs)

        return None


@dataclass
class LatticeMetadata:
    '''
    Dataclass to hold and map ontologies to a Lattice metadata spreadsheet.
    The functions allow for this object to read in a specific tab for a google sheet into a dataframe
    Subset should be a dict with key as column to subset, and value to filter column on
    '''
    
    sheet_id: str
    tab_name: str
    subset: dict[str, str] | None = None

    def get_gid(self):
        '''
        Given sheet id and tab name, return gid
        '''
        sheet_url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}'
        req = Request(sheet_url, headers={'User-Agent' : "Magic Browser"})
        s = urlopen(req)
        soup = BeautifulSoup(s, 'html.parser')
        tab_ids = {}
        pattern = re.compile('var bootstrapData = (.*?)};')
        for s in soup.find_all('script'):
            if pattern.search(str(s)):
                d = pattern.search(str(s)).group()[20:-1]
                data = json.loads(d)
                for t in data['changes']['topsnapshot']:
                    u = t[1].split('"')
                    if len(u) > 5:
                        tab_ids[u[5]] = u[1]
        return tab_ids[self.tab_name]

    def get_metadata_df(self):
        '''
        Given sheet id and gid, return lattice metadata in a dataframe
        Subset if subset dictionary is present
        '''
        url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}/export?format=csv&gid={self.gid}'
        response = requests.get(url)
        sample_df = pd.read_csv(BytesIO(response.content), comment="#", dtype=str).dropna(axis=1,how='all')

        if self.subset:
            sample_df = sample_df[sample_df[self.subset["column"]] == self.subset["filter_value"]]

        return sample_df

    def __post_init__(self):
        self.gid = self.get_gid()
        self.metadata_df = self.get_metadata_df()
