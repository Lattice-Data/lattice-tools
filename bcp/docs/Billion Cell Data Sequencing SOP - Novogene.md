# **Billion Cell Data Generation SOP**

The CZI Cell Science team has partnered with industry leaders in the single-cell genomics space to drive down costs and increase data generation efficiencies. This document defines the process and details for how to utilize the pipelines in place.

### **Sequencing**

### *Researchers will send 10X or Scale libraries directly to Psomogen, where they will sequence libraries using Ultima sequencers. A project tracker sheet has been created to track the progress of the steps below. Please work to keep the sheet as updated as possible (sheet to be linked once contracts are finalized).* 

* Researchers connect with CZI and Novogene representatives to confirm what libraries will be prepared and when they can be shipped  
* Novogene to send the necessary intake form to the researcher to collect the necessary experimental metadata  
* Novogene to ensure the intake form is completed and provide the shipping address to the researcher  
* The researcher will send libraries to the specified address  
* Novogene to confirm receipt and estimated sequencing start/completion date with CZI and researcher  
* Psomogen to sequence library and notify researcher of completion will sequence data and upload the raw data to 1)CZI AWS and 2\) 10X CLI for processing.  
  * CZI AWS info:  
    * Name: czi-Novogene  
    * ARN: aws:kms:us-east-1:577638397886:key/37a976fb-65d5-404e-a7bf-defb7f34206d

## **Data upload**

All projects will have a subdirectory within the S3 bucket and additional subdirectories for raw data files and cellranger outputs.

* Directory structure \- s3://czi-{provider}/{lastname}-{projectname}/  
  * should be all lowercase & use “-” (hyphen) delimiters  
* Each Novogene order spreadsheet submitted by the lab should be uploaded to the project directory   
  * one spreadsheet per order (aka: Novogene sample manifest)  
  * order number format: NVUS123456789-01  
  * s3://czi-{provider}/{lastname}-{projectname}/{order\_numer}/{}.xlsx

## 

## **Raw Sequence Data Sharing for 10x**

1. For each Billion Cell project,   
   1. The raw sequence data will be grouped by Sample and put into a raw/ directory within the project directory \- s3://czi-{provider}/{lastname}-{projectname}/{order\_numer}/{GroupID}/raw/  
      1. File naming convention – {RunID}-{GroupID}\_{Assay}-{UG-BC}\_S{SampleNumber}\_L{LaneNumber}\_R{ReadNumber}\_{OptionalPartNumber}.fastq.gz  
         1. Assay is **GEX**, **CRI**, or **ATAC** based on the “Analysis Target” column of the order template  
* Clarification about SampleID, LibraryID, GroupID 

  ![][image1]

  2. Example filepath: s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/AN00012345/CD4i\_R1L01/raw/416640-CD4i\_R1L01\_GEX-Z0238-CTGCACATTGTAGAT\_S1\_L001\_R1\_001.fastq.gz  
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

## 

## **Raw Sequence Data Sharing for sci-RNA-seq and Hash sci-RNA-seq**

2. For each Billion Cell project,   
   1. The raw sequence data will be grouped by RunID and put into a raw/ directory within the project directory- s3://czi-{provider}/{lastname}-{projectname}/{order\_numer}/{ExperimentID}/raw/{RunID}/  
      1. File naming convention – {RunID}-{GroupID}\_{Assay}-{UG-BC}.cram  
         1. Assay is one of the following:  GEX, hash\_oligo, GEX\_hash\_oligo  
      2. Example filepath:  s3://czi-novogene/shendure-seahub-bcp/NVUS123456789-01/2025-01-17-BIRTH/raw/427082/427082-PD-011\_GEX-Z1471-CAGTCGCCTGACGAT.cram  
         s3://czi-novogene/shendure-seahub-bcp/NVUS123456789-01/2025-01-17-BIRTH/raw/427082/427082-PD-011\_hash\_oligo-Z1471-CAGTCGCCTGACGAT.cram  
         s3://czi-novogene/shendure-seahub-bcp/NVUS123456789-01/2025-01-17-BIRTH/raw/427082/427082-PD-011\_GEX\_hash\_oligo-Z1471-CAGTCGCCTGACGAT.cram  
      3. UG-BC: format is “Z0123-CAGTCGCCTGACGAT” (UG 100 output)  
* Comments specific to order NVUS2024101701-09 Library name from SIF could be R100A (that corresponds to GroupID), the assay for this order is GEX\_hash\_oligo.  
  2. The file names noted below should be the same as when they are input into any pipeline. Any file renaming should be done prior to both processing and S3 upload.  
  3. In addition to the raw data, please ensure that the run details, including QC metrics,  are included and the original metadata sheet is uploaded to the raw data directory.  Confirmed files expected:

| File name | Desc |
| :---- | :---- |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.cram | Cram files containing reads for the library |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.json | General run stats \- JSON format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}.csv | General run stats \- CSV format |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-failure\_codes.csv | Sample Trimmer stats failure code |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-stats.csv | Sample Trimmer stats code |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-stats\_FlowQ.metric |  |
| {RunID}-{GroupID}\_{Assay}-{UG-BC}\_trimmer-stats\_SNVQ.metrix |  |

## 

## **Raw Sequence Data Sharing for ScaleBio Quantum and Scaleplex**

3. For each Billion Cell project,   
   1. The raw sequence data will be grouped by Ultima barcode and put into a raw/ directory within the project directory- s3://czi-{provider}/{lastname}-{projectname}/{order\_numer}/{ExperimentID}/raw/{RunID}/  
      1. File naming convention – {RunID}-{GroupID}\_{Assay}\_{UG\_RT}.cram  
         1. Assay is one of the following:  GEX, hash\_oligo, GEX\_hash\_oligo  
      2. Example filepath: s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-04/RNA3\_096/raw/430351/430351-RNA3-096H\_hash\_oligo\_QSR-8-SCALEPLEX\_10A.cram  
         s3://czi-novogene/trapnell-seahub-bcp/NVUS2024101701-04/RNA3\_096/raw/430351/430351-RNA3-096H\_GEX\_QSR-8\_10A.cram  
      3. **UG\_RT**: format is “QSR-8-{SCALEPLEX}\_10A” (UG 100 output), where “SCALEPLEX” is present in the hash\_oligo associated files.  
           
   2. The file names noted below should be the same as when they are input into any pipeline. Any file renaming should be done prior to both processing and S3 upload.  
   3. In addition to the raw data, please ensure that the run details, including QC metrics,  are included and the original metadata sheet is uploaded to the raw data directory.  Confirmed files expected:

| File name | Desc |
| :---- | :---- |
| {RunID}-{GroupID}\_{Assay}\_{UG\_RT}.cram | Cram files containing reads for the library |
| {RunID}-{GroupID}\_{Assay}\_{UG\_RT}.json | General run stats \- JSON format |
| {RunID}-{GroupID}\_{Assay}\_{UG\_RT}.csv | General run stats \- CSV format |
| {RunID}-{GroupID}\_{Assay}\_{UG}\_trimmer-failure\_codes.csv | Sample Trimmer stats failure code |
| {RunID}-{GroupID}\_{Assay}\_{UG}\_trimmer-stats.csv | Sample Trimmer stats code |
| {RunID}-{GroupID}\_{Assay}\_{UG}\_unmatched.cram | CRAM containing unmatched reads |
| {RunID}-{GroupID}\_{Assay}\_{UG}\_unmatched.csv | CSV stats for unmatched reads |
| {RunID}-{GroupID}\_{Assay}\_{UG}\_unmatched.json | JSON stats for unmatched reads |

### **Data processing (Additional details TBD by 10X)**

### *To increase interoperability across Billion Cell Projects, sequencing providers will run CellRanger on all raw data to centralize data processing steps using the parameters defined below.*

1. The 10X CLI provides CellRanger specifications and parameters based on the assay selected on the CLI. Additional details can be found below

### [**Cell Ranger**](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments) **Estimation of Expression**

* CellRanger version: 9.0.1

### **Genomes**

* Reference genomes ([aligned to 10X Genomics CellRanger pipeline](https://www.10xgenomics.com/support/software/cell-ranger/downloads)):

| organism | 10x flex | non-flex Cell Rangerand other aligner not counting multi-mapped reads |
| :---- | :---- | :---- |
| human | Probe set v1.1.0: 2024-A | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/GRCh38-v48.tar.gz |
| mouse | Probe set v1.1.1: 2024-A | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/GRCm39-vM37.tar.gz |
| zebrafish |  | reference package yet to be determined, but based on [this version](http://Danio_rerio.GRCz11.114.gtf.gz) |
| quail (Coturnix japonica) |  | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/Coturnix\_japonica\_2.1.tar.gz |
| anole (Anolis carolinensis) |  | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/Anolis.GCF\_035594765.1.tar.gz |
| frog (Ranitomeya imitator) |  | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/Ranitomeya.GCF\_032444005.1.tar.gz |
| turtle (Chrysemys picta bellii) |  | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/Chrysemys.GCF\_011386835.1.tar.gz |
| alligator (Alligator mississippiensis) |  | https://submissions-lattice.s3.us-west-1.amazonaws.com/BCP/references/Alligator.GCF\_030867095.1.tar.gz |

	

**Run Parameters**

* \--id \= from sequencing run metadata  
* \--sample \= from sequencing run metadata  
* \--create-bam=true ***\#non-10x flex only***  
* \--include-introns=true (default)  
* Min-crispr-umi   
  * Default \= 3

## **Processed Data Sharing**

*Below are the steps to access and share the processed data. Please note that the process to share processed data directly from the 10X CLI is still TBD.*

* See separate tab “[AWS Permission Steps](https://docs.google.com/document/d/1muWCbEj3XYD5B-zOC0-gPEjqOrYnJn3HUVAaa5f7kn4/edit?tab=t.4musrumcw9a0)” for additional instructions  
* The alignment output data will be grouped by GroupID and put into a processed/ directory and a subdirectory of the pipeline run (cellranger or scalerna). Outputs will be grouped by Run Date as follows \- s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{GroupID}/processed/cellranger/{Run\_Date}/  
  * Run\_Date is the date that the Cell Ranger run was started \- Run\_*YYYY-MM-DD*  
  * The unique GroupID should match that of the raw fastq files from the raw/ subdirectory,  input into cellranger.  
  * File naming convention \- the Cell Ranger outs/ directory should be uploaded to the above directory without changing file or folder names  
    * Example filepath:  
      s3://czi-novogene/weissman-scaling-in-vivo-perturb-seq-in-the-liver-and-beyond/AN00012345/CD4i\_R1L01/processed/cellranger/Run\_2025-02-01/outs/multi/count/raw\_feature\_bc\_matrix.h5

## **File Manifest**

### **When to Submit a Manifest?**

After both the raw and processed data for a given Sample have been uploaded to s3.

### **Where to Submit a Manifest?**

Directly to s3 in a specific path that will trigger the workflow. The data providers already have permission to upload to this bucket which simplifies implementation.	

Format:

\- sci-RNA-seq: s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{RunID}/file-manifest.json

\- 10x: s3://czi-{provider}/{lastname}-{projectname}/{order\_number}/{GroupID}/file-manifest.json

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

[image1]: <data:image/png;base64,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>