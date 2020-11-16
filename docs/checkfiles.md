Running checkfiles on an EC2 instance
---------------- 
Create a new spot instance from the checkfiles template
```
$ aws ec2 run-instances --launch-template LaunchTemplateName=checkfiles
```
If only matrix files are being checked, you can limit the duration of the instance to 1 hour
```
$ aws ec2 run-instances --launch-template LaunchTemplateName=checkfiles --instance-market-options MarketType=spot,SpotOptions={BlockDurationMinutes=60}
```
Identify the Public DNS for the instance and ssh to instance
```
$ aws ec2 describe-instances --filters Name=image-id,Values=ami-0e4035ae3f70c400f --query Reservations[*].Instances[*].[PublicDnsName] --output text
$ ssh -i "lattice_ec2.pem" ec2-user@<Public DNS>
```
Create directory to hold AWS credentials
```
$ mkdir .aws
```
Mount the storage device to the instance
```
$ sudo mkfs.ext4 /dev/nvme1n1
$ sudo mount /dev/nvme1n1 /mnt
$ sudo chown ec2-user /mnt
$ cd /mnt
```
Install python and necessary modules
```
$ sudo yum install python3
$ sudo pip3 install boto3 tables requests h5py
```
Define 3 variables for Lattice db permissions
```
$ export PROD_KEY=<> ; export PROD_SECRET=<> ; export PROD_SERVER=https://www.lattice-data.org
```
Copy 3 files from your local machine to the instance
```
$ scp -i lattice_ec2.pem ~/.aws/credentials ec2-user@<Public DNS>:.aws
$ scp -i lattice_ec2.pem checkfiles.py ec2-user@<Public DNS>:/mnt
$ scp -i lattice_ec2.pem lattice.py ec2-user@<Public DNS>:/mnt
```
Run checkfiles.py on the files
Upon completion, download the report
```
$ scp -i lattice_ec2.pem ec2-user@<Public DNS>:/mnt/<report txt> ./
```
Put the report in the [Checkfiles Drive](https://drive.google.com/drive/u/2/folders/1iomrTnd11hAH6S2iOMKciU6Llg1BZorP)
Get the instance ID and terminate the instance
```
$ aws ec2 describe-instances --filters Name=image-id,Values=ami-0e4035ae3f70c400f --query Reservations[*].Instances[*].[InstanceId] --output text
$ aws ec2 terminate-instances --instance-ids <instance ID>
```