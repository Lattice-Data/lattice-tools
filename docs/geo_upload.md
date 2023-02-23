# Submitting data to GEO
This document describes the process of taking raw matrix and fastq files for a Dataset in Lattice and submitting them to GEO along with the associated metadata.

Create metadata submission spreadsheet
----------------
- Run geo_metadata.py for the desired Dataset with the following commandline:
```
python geo_metadata.py -m local -d <Dataset_accession>
```
	This will create three files
	- \<Dataset_accession>\_metadata.csv: This is a csv of all the study, sample, protocol, and fasta file metadata.
	- \<Dataset_accession>\_md5sum.csv: This is a csv of all the fasta and raw matrix file names and checksums.
	- \<Dataset_accession>\_s3_uri.csv: This is a list of the s3 URI for all fastq and raw matrix files.
