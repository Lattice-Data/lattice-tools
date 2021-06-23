Running DCP_mapper on an EC2 instance
---------------- 
This should only be needed if there are files being transferred via an external ftp site. Files hosted at the Lattice S3 storage will not need to be downloaded between S3 and GCP, so can be done on your local machine.
 
Create a new spot instance from the checkfiles template
```
aws ec2 run-instances --launch-template LaunchTemplateName=checkfiles
```
Identify the Public DNS for the instance and ssh to instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-0e4035ae3f70c400f --query Reservations[*].Instances[*].[PublicDnsName] --output text
```
```
ssh -i lattice_ec2.pem ec2-user@<Public DNS>
```
Mount the storage device to the instance
```
sudo mkfs.ext4 /dev/nvme1n1 ; sudo mount /dev/nvme1n1 /mnt ; sudo chown ec2-user /mnt ; cd /mnt
```
Install python and necessary modules
```
sudo yum install python3
```
```
sudo pip3 install requests Pint google-api-python-client google-cloud-storage crcmod boto3
```
Define 3 variables for Lattice db permissions
```
export PROD_KEY=<> ; export PROD_SECRET=<> ; export PROD_SERVER=https://www.lattice-data.org
```
Copy 4 files from your local lattice-tools clone to the instance
```
scp -i lattice_ec2.pem lattice.py DCP_mapper.py ec2-user@<Public DNS>:/mnt
```
Copy DCP_mods directory your local lattice-tools clone to the instance
```
scp -i lattice_ec2.pem -r DCP_mods/ ec2-user@<Public DNS>:/mnt
```

Copy the GCP credentials json file from your local computer to the instance
```
scp -i lattice_ec2.pem gcp_creds.json ec2-user@<Public DNS>:/mnt
```
If you will be checking files hosted in the Lattice S3 storage, create directory to hold AWS credentials & set your variables to the location of the AWS & GCP credentials
```
mkdir .aws
```
```
export AWS_SHARED_CREDENTIALS_FILE=.aws/credentials ; export GOOGLE_APPLICATION_CREDENTIALS=/mnt/gcp_creds.json
```
... and copy credentials from your local machine to the instance
```
scp -i lattice_ec2.pem ~/.aws/credentials ec2-user@<Public DNS>:/mnt/.aws
```
Run the script on a given Dataset
```
python3 DCP_mapper.py -m prod -d <dataset_identifier>
```
