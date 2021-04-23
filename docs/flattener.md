## flattener.py
This script takes in a Lattice identififer of a final matrix and creates the corresponding h5ad that conforms to cellxgene requirements (https://github.com/chanzuckerberg/single-cell-curation/tree/main/docs).

Installation requirements
----------------
Create and activate lattice\_submit environment as documented on https://github.com/Lattice-Data/lattice-tools. Additional python library to install is:
```
$ pip install rpy2
```

For converting a Seurat object to h5ad format, R is required to be installed on the machine (https://www.r-project.org/). The required libraries are:
```

```

Running flattener.py
----------------
```
$ python flattener.py --mode local --file LATDF119AAA
```
--mode: Use 'local' or 'prod' to use the local or production database instance, respectively
--file: Any identifier for the matrix of interest

The script will produce a h5ad file in the current directory where the script is being run from. The file name corresponds to the accession of the final matrix, appended with the version of the flattener.py used to create the file. A temporary directory 'matrixi\_files/' will be created to hold downloaded and intermediate files, and, therefore, make sure there is no such directory present.

Version update logging
----------------
**Version 2**: 
-Add ability to demultiplex metadata from experiments pooled at the library entity. The requirement is that the the 'author\_donor\_column' metadata field is filled out in the final matrix object. 
-Convert development_stage to term name so that corresponds with term id.
-Ethnicity is empty string when ethnicity is unknown.


**Version 1**: Initial version, which can take a h5ad or Seurat object as input from RNA-seq and ATAC-seq assays. For RNA-seq assays, the raw matrix is subsetted from the Cell Ranger filtered raw counts. For ATAC-seq assays, the corresponding raw matrix from the activity gene matrix is used.
