QA Testing
----------------
QA validation tests for moving towards schema 4.1.0 migration. Almost all tests will be through pytest; `validate_notebook_workflow.ipynb` checks that validation through the curator notebook workflow does not break.


General process
---------------- 
Setup test env, run pytest through the command line, run the notebook.


Installation
---------------- 
To get proper testing environment: 
- Local, up-to-date [CZI single-cell-curation repo](https://github.com/chanzuckerberg/single-cell-curation)
- `pytest` is part of test env
- `cellxgene-schema` pip package is installed/updated to at least 4.0.1 from this repo.

Create a seperate test environment since schema will be > 4.0.0. This can be done by cloning the lattice env or creating a new one.

### 1. Clone conda lattice env:
```
conda create --name cxg4testing --clone lattice
```

### OR Create a new test env
```
conda create --name cxg4testing python=3.9
```
Activate this env
```
conda activate cxg4testing
```
Install further dependencies
```
conda install -c conda-forge pint pandas jsonschema boto3 jupyter bs4
```

```
pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic crcmod lxml pytest
```
### 2. Check pytest installed with `pip list | grep pytest`

### 3. Latest version of `cellxgene-schema`
```
pip install git+https://github.com/chanzuckerberg/single-cell-curation/@main#subdirectory=cellxgene_schema_cli
```
If cellxgene-schema does not seem to update, uninstall the package via:
```
pip uninstall cellxgene-schema
```
and reinstall like above

The tests assume the standard location of cloned repos that Lattice uses:
```
~/GitClones/CZI/single-cell-curation/
```

Running tests
---------------- 
Make sure the test env is activated.
Navigate to this directory and run pytest from command line with desired flags
```
pytest -vv
```