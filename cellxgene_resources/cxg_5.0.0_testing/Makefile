venv:
	python3.9 -m venv venv;

tests: latest-cli
	pytest -vvv

venv-install: requirements.txt
	pip install -r requirements.txt;
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@main#subdirectory=cellxgene_schema_cli

latest-cli:
	pip uninstall cellxgene-schema
	pip install git+https://github.com/chanzuckerberg/single-cell-curation/@main#subdirectory=cellxgene_schema_cli

clean:
	rm -rf __pycache__ 
	rm -rf venv
