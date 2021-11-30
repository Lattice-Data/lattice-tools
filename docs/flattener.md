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
Seurat
Signac
SeuratDisk
reticulate
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
**Version 4**:
- Corresponds with https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/schema.md
- Add tyrer_cuzick_lifetime_risk, enriched_cell_types, mapped_reference_annotation, and enrichment-factors as optional metadata fields for obs
- Add is_primary_data, organism_ontology_term_id, sex_ontology_term_id as required metadata fields for obs.
- Removed *_onotology_term_name fields as they will be populated by cxg portal
- Add feature_is_filtered to var and feature_biotype to both var and raw.var
- Filter both var and raw.var to pinned gene annotation (GENCODE v38)
- var index must be Ensembl IDs
- Pad matrix with implied zeros to make X and raw.X the same shapes
- For datasets with raw matrices mapped to multiple annotations, will do out join for raw.X and inner join with padded implied zeros for X  
- Remove reported_disease and donor_age if disease_ontology_term_name and development_stage_ontology_name are redundant
- For uns, add schema_version and removed organism, organism_ontology_term_id, deafult_field, version.corpora_encoding_version, and version.corpora_schema_version.
- Add enrichment_factors and cell_state to optional fields
- Use development_ontology_at_collection and age_development_stage_redundancy at Tissue rather than development ontology at Donor
- Remove cell_type_category
- Handles multiplexed donor datasets correctly
- Looks at feature_keys of ProcessedMatrix to determine whether or not the there needs to be mapping of Ensembl IDs
- Distinguish between serially linked vs pooled suspensions, and ignores first suspension if there are serially linked suspensions

**Version 3**:
- Allow reading from h5ad file format for raw count matrices
- Raw matrix will be an outer join to allow for merging of matrices with varying feature counts
- Add ability to read data from spatial transcriptomic assays
- Transfer additional layers from the 'layers' attribute of the contributor final h5ad matrix to cxg h5ad
- Permit final matrices that do not have prefix/suffix added to cell barcodes in cell_label_mappings
- Transfer X_spatial embedding
- Looks for TissueSection objects in experimental graph


**Version 2**: 
- Add ability to demultiplex metadata from experiments pooled at the library entity. The requirement is that the the 'author\_donor\_column' metadata field is filled out in the final matrix object. 
- Convert development_stage to term name so that corresponds with term id.
- Ethnicity is empty string when ethnicity is unknown.
- Add optional columns, which are removed if they are all empty or unreported values.
- Update assay ontology logic to obtain information from linked OntologyTerm
- Add logic for cell culture tissue ontology to be from organ slims
- If cell culture ontology is not UBERON, will go get most specific tissue slim


**Version 1**: Initial version, which can take a h5ad or Seurat object as input from RNA-seq and ATAC-seq assays. For RNA-seq assays, the raw matrix is subsetted from the Cell Ranger filtered raw counts. For ATAC-seq assays, the corresponding raw matrix from the activity gene matrix is used.
