CONDA_NAME = cxg51testing
CLONE_ENV = lattice
GIT_REF = main
PYTHON_VERSION = 3.11

.PHONY: help
help:
	@echo "Makefile for 5.1.0 QA testing"
	@echo "============ Default values ============"
	@echo "PYTHON_VERSION="$(PYTHON_VERSION)
	@echo "To create: CONDA_NAME="$(CONDA_NAME)
	@echo "To clone: CLONE_ENV="$(CLONE_ENV)
	@echo "Make cellxgene-schema from: GIT_REF="$(GIT_REF)
	@echo "========================================"
	@echo "Change any of the above as needed:"
	@echo "========================================"
	@echo "make conda PYTHON_VERSION=3.10"
	@echo "make conda-install GIT_REF=v5.0.0"
	@echo "make conda-vars CLONE_ENV=myconda"
	@echo "make cxg-schema GIT_REF=v5.0.1"
	@echo "make tests GIT_REF=v5.0.2"
	@echo "========================================"
	@echo "Or run defaults:"
	@echo "========================================"
	@echo "make conda"
	@echo "conda activate cxg51testing"
	@echo "make conda-install"
	@echo "make tests"

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

cxg-schema:
	pip uninstall cellxgene-schema
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli

tests: cxg-schema
	pytest -vvv

clean:
	rm -rf __pycache__ 
	rm -rf venv
