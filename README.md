# lattice-tools
External scripts used to interact with the Lattice Database

Environment configuration
---------------- 
1. Create a virtual environment. This example uses anaconda. Other options would also work, like venv or pyenv
    ```
    conda create --name lattice python=3.9
    ```
    You will need to be in this environment for the following instructions
    ```
    conda activate lattice
    ```

2. Install the following packages
    ```
    conda install -c conda-forge pint pandas jsonschema boto3 jupyter bs4
    ```
    ```
    pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic-bin==0.4.14 crcmod cellxgene-schema
    ```
3. Define variables in your environment based on the various servers you might submit to based on an alias for each server (`ALIAS_KEY`, `ALIAS_SECRET`, `ALIAS_SERVER`). For example, when submitting to the production instance of Lattice, you might call this `prod`.
So you'd define the following three variables.

	`$ conda env config vars set PROD_KEY=<key>`

	`$ conda env config vars set PROD_SECRET=<secret>`

	`$ conda env config vars set PROD_SERVER=https://www.lattice-data.org/`

Your demo access will be the same, but the demo server will change with each new demo.

	`$ conda env config vars set DEMO_KEY=<key>`

	`$ conda env config vars set DEMO_SECRET=<secret>`

4. After defining those, you'll need to reactivate your environment
    ```
    conda activate lattice
    ```
	You can then confirm that they are defined
    ```
    conda env config vars list
    ```

Available tools
---------------- 
* **qcmetrics_reader.py**
Transforms quality metrics and other processing information from various files of a standard CellRanger outs/ directory into the Lattice schema

* **checkfiles.py**
Gathers data file content information and compares with submitted metadata [run instructions](docs/checkfiles.md)

* **DCP_mapper.py**
Transforms a Lattice Dataset into HCA DCP-approved schema and stages at the DCP for submission to the HCA Portal [run instructions](docs/DCP_mapper.md)

* **flattener.py**
Transforms a contributor matrix, raw count data, and Lattice metadata into a cellxgene-approved matrix file [run instructions](docs/flattener.md)

* **make_template.py**
Produces a tabular representation of Lattice schema submittable properties, for ease of wrangling

* **submit_metadata.py**
Transforms tabulated metadata into json objects and posts/patches to the Lattice DB [use instructions](docs/submit_metadata.md)
