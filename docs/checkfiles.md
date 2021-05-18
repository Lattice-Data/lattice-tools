Running checkfiles on an EC2 instance
---------------- 
Create a new spot instance from the checkfiles template
```
aws ec2 run-instances --launch-template LaunchTemplateName=checkfiles
```
If only matrix files are being checked, you can limit the duration of the instance to 1 hour
```
aws ec2 run-instances --launch-template LaunchTemplateName=checkfiles --instance-market-options MarketType=spot,SpotOptions={BlockDurationMinutes=60}
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
sudo pip3 install boto3 tables requests h5py crcmod
```
Define 3 variables for Lattice db permissions
```
export PROD_KEY=<> ; export PROD_SECRET=<> ; export PROD_SERVER=https://www.lattice-data.org
```
Copy 2 files from your local machine to the instance
```
scp -i lattice_ec2.pem lattice.py checkfiles.py ec2-user@<Public DNS>:/mnt
```
If you will be checking files hosted in the Lattice S3 storage, create directory to hold AWS credentials & set your variable to the location of the credentials
```
mkdir .aws ; export AWS_SHARED_CREDENTIALS_FILE=.aws/credentials
```
... and copy credentials from your local machine to the instance
```
scp -i lattice_ec2.pem ~/.aws/credentials ec2-user@<Public DNS>:/mnt/.aws
```
Run checkfiles.py on the files
```
python3 checkfiles.py -m prod --update -q "report/?type=RawSequenceFile&audit.ERROR.category=file+not+validated" &
```
If you wish to exit your ssh session and leave checkfiles running, first disown the job so it won't receive SIGHUP when you logout
```
disown -h %1
```
```
exit
```
Upon completion, download the report
```
scp -i lattice_ec2.pem ec2-user@<Public DNS>:/mnt/<report txt> ./
```
Put the report in the [Checkfiles Drive](https://drive.google.com/drive/u/2/folders/1iomrTnd11hAH6S2iOMKciU6Llg1BZorP)
Get the instance ID and terminate the instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-0e4035ae3f70c400f --query Reservations[*].Instances[*].[InstanceId] --output text
```
```
aws ec2 terminate-instances --instance-ids <instance ID>
```