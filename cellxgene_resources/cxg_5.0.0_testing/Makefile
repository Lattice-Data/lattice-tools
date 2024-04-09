GIT_REF = main
CONDA_NAME = cxg51testing
CLONE_ENV = lattice
PYTHON_VERSION = 3.11

venv:
	python$(PYTHON_VERSION) -m venv venv;

venv-install: requirements.txt 
	pip install -r requirements.txt;
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli
	
conda-clone:
	conda create --name $(CONDA_NAME) --clone $(CLONE_ENV)

conda:
	conda create --name $(CONDA_NAME) python=$(PYTHON_VERSION)

conda-install: 
	conda install -c conda-forge pint pandas jsonschema boto3 jupyter bs4
	pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic crcmod lxml pytest
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli

conda-vars:
	cp /Users/$(USER)/anaconda3/envs/$(CLONE_ENV)/conda-meta/state /Users/$(USER)/anaconda3/envs/$(CONDA_NAME)/conda-meta/state

cxg-cli:
	pip uninstall cellxgene-schema
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli

tests: cxg-cli
	pytest -vvv

clean:
	rm -rf __pycache__ 
	rm -rf venv
