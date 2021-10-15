General process
---------------- 
cxg_audit.py was run using --refresh option to download the h5ad version of each published dataset in the cellxgene portal, and upload the file to the Lattice S3 original/ directory using the file naming collectionID_datasetID.h5ad

The audit report was imported into a Google Sheet for collaborative examination and curation. The actions to take were noted in the sheet based on information found in associated publications, code repositories, data archive records, etc.

make_guides.py was run to extract the actions noted for each dataset and write into a json-formatted guide, stored in guides/.

obs_uns_migrator.py was run to manipulate the h5ad in original/ as noted in the corresponding json guide, and uploaded the resulting file in the Lattice S3 working/ directory.

Using information in associated publications, code repositories, data archive records, etc., a gtf or a data source file (an alternate matrix file) was identified that contained a mapping between gene symbols and Ensembl IDs for that particular dataset (in some cases, a gtf mapped symbols to RefSeq IDs, which were then mapped to Ensembl IDs using the HGNC BioMart). The gtf or source file was parsed into a 2-column mapping file with gene_symbols and gene_ids.

var_migrator.py was run to manipulate the h5ad from the Lattice S3 working/ directory by matching values in the index of var and raw.var, if present, to gene_symbols in the specified mapping file, and replacing with the corresponding gene_id. This tool also added other schema 2.0.0 requirements for var and raw.var, if present, before running cellxgene-schema validate and uploading the resulting file in the Lattice S3 final/ directory.


Running checkfiles on an EC2 instance
---------------- 
Create a new spot instance from the checkfiles template (r5ad.8xlarge)
```
aws ec2 run-instances --launch-template LaunchTemplateName=cxg
```
Identify the Public DNS for the instance and ssh to instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-04b6c97b14c54de18 --query Reservations[*].Instances[*].[PublicDnsName] --output text
```
```
ssh -i lattice_ec2.pem ec2-user@<Public DNS>
```
Mount the storage device to the instance & create directory to hold AWS credentials & set your variable to the location of the credentials & install necessary modules
```
sudo mkfs.ext4 /dev/nvme1n1 ; sudo mount /dev/nvme1n1 /mnt ; sudo chown ec2-user /mnt ; cd /mnt ; mkdir .aws ; export AWS_SHARED_CREDENTIALS_FILE=credentials ; sudo pip3 install boto3 requests pandas scanpy cellxgene-schema; pip3 install numpy==1.20 anndata==0.7.5
```
Copy relevant script & aws credentials from your local machine to the instance. If running var_migrator.py, then gene_symbol_custom.py is already required.
```
scp -i lattice_ec2.pem cxg_audit.py ~/.aws/credentials ec2-user@<Public DNS>:/mnt
```
If running obs_uns_migrator.py, then the guides must also be copied.
```
scp -i lattice_ec2.pem -r guides/ ec2-user@<Public DNS>:/mnt
```
Run the thing
```
python3 cxg_audit.py
```
Upon completion, download the various report files
```
scp -i lattice_ec2.pem ec2-user@<Public DNS>:/mnt/report.txt ./
```
Get the instance ID and terminate the instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-04b6c97b14c54de18 --query Reservations[*].Instances[*].[InstanceId] --output text
```
```
aws ec2 terminate-instances --instance-ids <instance ID>
```
