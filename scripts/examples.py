import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Protocol, Self, TypeAlias

import pandas as pd


# value object, use some validation to check it is correct, retains properties of initial type
# since it inherited from it, can use to type arguments and return values
class S3Uri(str):
    def __new__(cls, value: Any) -> Self:
        val = str(value)
        if not val.startswith("s3://"):
            raise ValueError(f"String input '{val}' for '{cls.__name__}' does not start with '{val}'")
        return super().__new__(cls, val)


# simple dataclass used as a structure to logically group together data/attributes
# can have nested classes/objects for deeper types
@dataclass
class FileMeta:
    file_name: str
    path: Path | str
    s3_uri: S3Uri | str


# type alias to keep more complicated typing cleaner in function signatures
# dictionary with string keys and a list of values with either ints or floats
LatticeJSONReport: TypeAlias = dict[str, str | list[Any]]
TissueJSON: TypeAlias = dict[str, list[str]]

# now we know type and structure of some variable
tissue_results: TissueJSON = {}


# string literal typing for str that should be just a few values
Mode: TypeAlias = Literal["prod", "demo"]


# combining some constants with value object validation and typing
DATASET_RE = re.compile(r"LATDS\d{3}[A-Z]{3}")
LIBRARY_RE = re.compile(r"LATLB\d{3}[A-Z]{3}")
OTHER_RE = re.compile(r"LATDF\d{3}[A-Z]{3}")

STRING_RE = {
    "DataSetString": DATASET_RE,
    "RawMatrixString": OTHER_RE,
    "RawSequenceString": OTHER_RE,
    "ProcessedMatrixString": OTHER_RE,
    "SequenceAlignmentString": OTHER_RE,
    "LibraryString": LIBRARY_RE,
}


# value object combined with inheritance. dict lookup takes name of class to make sure regex is valid
# child classes used to have correct name of string, doesn't need further implementation due to inheritance
class FileObjectString(str):
    def __new__(cls, value: Any) -> Self:
        val = str(value)
        string_re = STRING_RE[cls.__name__]
        if not string_re.match(val):
            raise ValueError(f"String input '{val}' for '{cls.__name__}' does not match pattern '{string_re}'")
        return super().__new__(cls, val)


class DataSetString(FileObjectString):
    """
    String that should be "LATDS{3 digits}{3 uppercase letters}"
    """
    

class RawMatrixString(FileObjectString):
    """
    String that should be "LATDF{3 digits}{3 uppercase letters}"
    """


class RawSequenceString(FileObjectString):
    """
    String that should be "LATDF{3 digits}{3 uppercase letters}"
    """


class ProcessedMatrixString(FileObjectString):
    """
    String that should be "LATDF{3 digits}{3 uppercase letters}"
    """


class SequenceAlignmentString(FileObjectString):
    """
    String that should be "LATDF{3 digits}{3 uppercase letters}"
    """


class LibraryString(FileObjectString):
    """
    String that should be "LATLB{3 digits}{3 uppercase letters}"
    """


# protocols can be a nice way to provide containers/expectations of data and methods without inheritance
# anything typed for LoadResource should have a data attribute and the 2 methods available
# can add additional classes to handle different data like JSONs, csvs, google sheets, etc
class LoadResource(Protocol):
    data: LatticeJSONReport | pd.DataFrame

    def load_data(self, input: Path | str):
        ...

    def format_data(self) -> pd.DataFrame:
        ...


def export_to_geo(resource: LoadResource, path: Path | str):
    resource.load_data(path)
    resource.format_data()


# function signature typed out for arguments and returns, mix of custom types and 
# objectes defined by other packages
def main(dataset: DataSetString, mode: Mode) -> pd.DataFrame:
    print(dataset)
    print(mode)

    processed_matrix: ProcessedMatrixString = ProcessedMatrixString("LATDF384UIO")

    if processed_matrix.startswith("LATDS"):
        raise ValueError("Check processed matrix input")

    return pd.DataFrame(dataset)


if __name__ == "__main__":
    main(dataset=DataSetString("LATDS977NKW"), mode="prod")
