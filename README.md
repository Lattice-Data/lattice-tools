# lattice-tools
Scripts used by the Lattice data coordination team for single cell data wrangling

## Environment configuration

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
    pip install requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy python-magic-bin crcmod cellxgene-schema
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

## Available tools

### cellxgene_resources/<br>*for curating towards [CZ CELLxGENE Discover](cellxgene.cziscience.com)*
* **curation_qa.ipynb**
Quality assurance checks on an AnnData object

* **curation_sample_code.ipynb**
Various samples of how to manipulate an AnnData object during curation

* **upload_local.ipynb**
Submitting local files to CELLxGENE<br>
Please note: <br>
    This script utilizes the [single-cell-curation](https://github.com/chanzuckerberg/single-cell-curation/tree/main) repo which should be cloned to the following directory `~/GitClones/CZI/` and CXG API keys should be stored in `~/Documents/keys/cxg-api-key.txt`


### scripts/<br>*for curating towards or out of [Lattice DB](lattice-data.org)*
* **checkfiles.py**
Gathers data file content information and compares with submitted metadata [run instructions](docs/checkfiles.md)
If running locally, may need to install [Homebrew](https://brew.sh/) and `brew install md5sha1sum` so `md5sum` can run from checkfiles

* **DCP_mapper.py**
Transforms a Lattice Dataset into HCA DCP-approved schema and stages at the DCP for submission to the HCA Portal [run instructions](docs/DCP_mapper.md)<br>
Requires additional steps:
    ```
    pip install google-api-python-client google-cloud-storage
    ```
    `$ conda env config vars set GOOGLE_APPLICATION_CREDENTIALS=<creds.json>`

* **DCP_project_ready.ipynb**
Validates a project staged for submission to the HCA Data Portal.
Requires additional step:
    ```
    $ conda install -c anaconda more-itertools
    ```

* **flattener.py**
Transforms a contributor matrix, raw count data, and Lattice metadata into a cellxgene-approved matrix file [run instructions](docs/flattener.md)

* **geo_metadata.py**
Transforms a Lattice Dataset into GEO submission format

* **make_template.py**
Produces a tabular representation of Lattice schema submittable properties, for ease of wrangling<br>
Requires additional steps:<br>
    Follow instructions [here](https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html) to enable API & generate credentials<br>
    `$ conda env config vars set CLIENT_SECRET_FILE=<creds.json>`

* **qcmetrics_reader.py**
Transforms quality metrics and other processing information from various files of a standard CellRanger outs/ directory into the Lattice schema

* **query_by_dataset_lab.ipynb**
Return Donor, Sample, or Suspension objects from the Lattice DB for a given Dataset or Lab

* **s3_recent_uploads.ipynb**
Return files recently uploaded to the submitter S3 buckets

* **submit_metadata.py**
Transforms tabulated metadata into json objects and posts/patches to the Lattice DB [use instructions](docs/submit_metadata.md)

* **validate_demo.ipynb**
Compares various aspects of the production DB and a specified demo DB to identify potential bugs.

* **validate_checksums.py**
Identifies any duplicated files in the Lattice DB. To be executed after each checkfiles run.