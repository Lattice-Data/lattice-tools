# **Billion Cell Data Generation SOP**

The CZI Cell Science team has partnered with industry leaders in the single-cell genomics space to drive down costs and increase data generation efficiencies. This document defines the process and details for how to utilize the pipelines in place.

### **Sequencing**

### *Researchers will send 10X or Scale libraries directly to Psomogen, where they will sequence libraries using Ultima sequencers. A project tracker sheet has been created to track the progress of the steps below. Please work to keep the sheet as updated as possible (sheet to be linked once contracts are finalized).* 

* Researchers connect with CZI and Psomagen representatives to confirm what libraries will be prepared and when they can be shipped  
* Psomagen to send the necessary intake form to the researcher to collect the necessary experimental metadata  
* Psomagen to ensure the intake form is completed and provide the shipping address to the researcher  
* The researcher will send libraries to the specified address  
* Psomagen to confirm receipt and estimated sequencing start/completion date with CZI and researcher  
* Psomogen to sequence library and notify researcher of completion will sequence data and upload the raw data to 1)CZI AWS and 2\) 10X CLI for processing.  
  * CZI AWS info:  
    * Name: czi-psomagen  
    * ARN: aws:kms:us-east-1:577638397886:key/37a976fb-65d5-404e-a7bf-defb7f34206d

## **Data upload**

All projects will have a subdirectory within the S3 bucket and additional subdirectories for raw data files and cellranger outputs.

* Directory structure \- s3://czi-{provider}/{lastname}-{projectname}/  
  * should be all lowercase & use “-” (hyphen) delimiters  
* Each Psomagen order spreadsheet submitted by the lab should be uploaded to the project directory   
  * one spreadsheet per order (aka: PSOM sample manifest)  
  * order\_numer format: AN00012345  
  * s3://czi-{provider}/{lastname}-{projectname}/{order\_numer}/{}.xlsx

## **Raw Sequence Data Sharing**

1. For each Billion Cell project,   
   1. The raw sequence data will be grouped by Sample and put into a raw/ directory within the project directory \- s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{GroupID}/raw/  
      1. File naming convention – {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S{SampleNumber}\_L{LaneNumber}\_R{ReadNumber}\_{OptionalPartNumber}.fastq.gz  
         1. Assay is **GEX**, **CRI**,  **ATAC, or viral\_ORF** based on the “Analysis Target” column of the order template  
      2. Example filepath: s3://czi-psomagen/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/AN00012345/CD4i\_R1L01/raw/416640-CD4i\_R1L01\_GEX-Z0238-CTGCACATTGTAGAT\_S1\_L001\_R1\_001.fastq.gz  
      3. UG: format is “Z0123” (UG 100 output)  
   2. The file names noted below should be the same as when they are input into any pipeline. Any file renaming should be done prior to both processing and S3 upload.  
   3. In addition to the raw data, please ensure that the run details, including QC metrics,  are included and the original metadata sheet is uploaded to the raw data directory.  Confirmed files expected:

| File name | Desc |
| :---- | :---- |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.csv | General run stats \- CSV format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.json | General run stats \- JSON format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.scRNA.applicationQC.h5 | scRNA specific QC \- h5 format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.scRNA.applicationQC.html | scRNA specific QC \- html format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_Log.final.out | STAR Solo log from the QC sample |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_Log.out | STAR Solo log from the QC sample |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_Log.progress.out | STAR Solo log from the QC sample |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_ReadsPerGene.out.tab | STAR Solo gene count from the QC sample |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R1\_001.csv | General stats for R1 FQ \- CSV format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R1\_001.fastq.gz | R1 FQ file for use with CellRanger |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R1\_001.json | General stats for R1 FQ \- JSON format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R1\_001\_sample.fastq.gz | R1 Sample FQ file used for QC  |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R2\_001.csv | General stats for R1 FQ \- CSV format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R2\_001.fastq.gz | R2 FQ file for use with CellRanger |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R2\_001.json | General stats for R2 FQ \- JSON format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S1\_L001\_R2\_001\_sample.fastq.gz | R2 Sample FQ file used for QC  |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_SJ.out.tab | STAR solo read counts |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-failure\_codes.csv | Sample Trimmer stats failure code |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-stats.csv | Sample Trimmer stats code |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_unmatched.cram | CRAM containing unmatched reads |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_unmatched.csv | CSV stats for unmatched reads |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_unmatched.json | JSON stats for unmatched reads |

### **Data processing (Additional details TBD by 10X)**

### *To increase interoperability across Billion Cell Projects, sequencing providers will run CellRanger on all raw data to centralize data processing steps using the parameters defined below.*

1. The 10X CLI provides CellRanger specifications and parameters based on the assay selected on the CLI. Additional details can be found below

### [**Cell Ranger**](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments) **Estimation of Expression**

* CellRanger version: 9.0.1

### **Genomes**

* Reference genomes ([aligned to 10X Genomics CellRanger pipeline](https://www.10xgenomics.com/support/software/cell-ranger/downloads)):

| organism | 10x flex | non-flex Cell Rangerand other aligner not counting multi-mapped reads |
| :---- | :---- | :---- |
| human | Probe set v1.1.0: 2024-A | Biotype-filtered Ensembl 114 (GENCODE v48) |
| mouse | Probe set v1.1.1: 2024-A | Biotype-filtered Ensembl 114 (GENCODE M37) |

	  
**Run Parameters**

* \--id \= from sequencing run metadata  
* \--sample \= from sequencing run metadata  
* \--create-bam=true *\#non-10x flex only*  
* \--include-introns=true (default)  
* Min-crispr-umi   
  * Default \= 3

## **Processed Data Sharing**

*Below are the steps to access and share the processed data. Please note that the process to share processed data directly from the 10X CLI is still TBD.*

* See separate tab “[AWS Permission Steps](?tab=t.4musrumcw9a0)” for additional instructions  
* The alignment output data will be grouped by Sample and put into a processed/ directory and a subdirectory of the pipeline run (cellranger or scalerna). Outputs will be grouped by Run Date as follows \- s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{GroupID}/processed/cellranger/{Run\_Date}/  
  * Run\_Date is the date that the Cell Ranger run was started \- Run\_*YYYY-MM-DD*  
  * The unique GroupID should match that of the raw fastq files from the raw/ subdirectory,  input into cellranger.  
  * File naming convention \- the Cell Ranger outs/ directory should be uploaded to the above directory without changing file or folder names  
    * Example filepath:  
      s3://czi-psomagen/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/AN00012345/CD4i\_R1L01/processed/cellranger/Run\_2025-02-01/outs/multi/count/raw\_feature\_bc\_matrix.h5

## **File Manifest**

### **When to Submit a Manifest?**

After both the raw and processed data for a given Sample have been uploaded to s3.

### **Where to Submit a Manifest?**

Directly to s3 in a specific path that will trigger the workflow. The data providers already have permission to upload to this bucket which simplifies implementation.	

Format: s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{GroupID}/file-manifest.json

### **Manifest Schema**

```json
{
  "type": "array",
  "items": {
    "type": "string",
    "description": "List of S3_URI to check in the manifest."
  }
}

```

Example Fille Manifest

```json
[
  "s3://czi-psomagen/project/AN1234/group-ID/raw/file-1.fastq.gz",
  "s3://czi-psomagen/project/AN1234/group-ID/raw/file-2.csv",
  ...
]

```

