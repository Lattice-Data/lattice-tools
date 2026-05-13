## fragment_curator.py
This script takes in a Lattice identififer of a final processed matrix and creates the filtered, concatenated fragment file that can be submitted to CELLxGENE.

This curation is done in parallel as much as possible due to the massive fragment files that correspond to raw matrix files.

Installation requirements
----------------
Follow instructions in the `lattice-tools` README for creating a local conda environment. 

If running on JupyterHub, all required packages are already present. Make sure to export env variables for key, secret, and server to connect to Lattice DB:
```
export PROD_KEY={key string}
export PROD_SECRET={secret string}
export PROD_SERVER={server string}
```

Lattice Requirements Before Running
----------------
Make sure the following metadata exists on Lattice:
 - ProcessedMatrixFile
    - s3_uri
    - cell_label_location
    - cell_label_mappings
 - For each RawMatrixFile:
    - s3_uri
    - fragment_file_s3_uri

The script should provide helpful errors if this metadata cannot be found when querying Lattice

Running fragment_curator.py
----------------
```
$ python fragment_curator.py --mode prod --file LATDF366MLP [OPTIONAL] --deduplicate [OPTIONAL] --golang
```
--mode: Use 'local' or 'prod' to use the local or production database instance, respectively

--file: Any identifier for the matrix of interest

--deduplicate [OPTIONAL]: use to deduplicate fragment files if CELLxGENE validation finds duplicate entries exist

--golang [OPTIONAL]: by default, script uses python/pandas to filter fragment files. This option uses a binary written in Go for filtering. This is faster and way more memory efficent. Currently this binary cannot use regex for barcode replacement and the script will error if this is necessary.

The script will produce a filtered, concatenated, gzipped tsv file in this format:
```
{ProcessedMatrixFile Accession}_concatenated_filtered_fragments.tsv.gz
``` 
in the `scripts/atac_fragments` directory. 

There will be summary logging to the terminal and detailed logging info saved to the following file in the `scripts/` directory:
```
{ProcessedMatrixFile Accession}_outfile_fragment_curator.log
```
