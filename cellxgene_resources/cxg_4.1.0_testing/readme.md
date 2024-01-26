Running tests
----------------
QA validation tests for moving towards schema 4.1.0 migration.


General process
---------------- 
To get proper testing environment: 
    make sure CZI single-cell-curation repo is up to date 
    cellxgene-schema pip package is installed/updated to at least 4.0.1 from this repo.
    pytest is part of test env

Clone the lattice conda environment or create new environment as follows:

$ conda create --name cxg4testing python=3.9
$ conda activate cxg4testing
$ conda install -c conda-forge pint pandas jsonschema boto3 jupyter bs4
$ pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic crcmod cellxgene-schema lxml pytest
$ pip install git+https://github.com/chanzuckerberg/single-cell-curation/@main#subdirectory=cellxgene_schema_cli


Running tests
---------------- 
Navigate to this directory and run pytest from command line with desired flags
$ pytest -vv
