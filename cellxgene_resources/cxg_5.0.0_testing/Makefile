GIT_REF = main
CLONE_ENV = lattice

venv:
	python3.10 -m venv venv;

venv-install: requirements.txt 
	pip install -r requirements.txt;
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli
	
conda-clone:
	conda create --name cxg51testing --clone $(CLONE_ENV)

conda:
	conda create --name cxg51testing python=3.10

conda-install: 
	conda install -c conda-forge pint pandas jsonschema boto3 jupyter bs4
	pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic crcmod lxml pytest
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli

cxg-cli:
	pip uninstall cellxgene-schema
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@$(GIT_REF)#subdirectory=cellxgene_schema_cli

tests: cxg-cli
	pytest -vvv

clean:
	rm -rf __pycache__ 
	rm -rf venv
