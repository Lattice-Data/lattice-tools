Running checkfiles on an EC2 instance
---------------- 
Create a new spot instance from the checkfiles template
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
sudo mkfs.ext4 /dev/nvme1n1 ; sudo mount /dev/nvme1n1 /mnt ; sudo chown ec2-user /mnt ; cd /mnt ; mkdir .aws ; export AWS_SHARED_CREDENTIALS_FILE=credentials ; sudo pip3 install boto3 requests pandas scanpy ; pip3 install numpy==1.20
```
Copy audit or migration script & aws credentials from your local machine to the instance
```
scp -i lattice_ec2.pem cxg_audit.py ~/.aws/credentials ec2-user@<Public DNS>:/mnt
```
```
scp -i lattice_ec2.pem -r guides/ ec2-user@<Public DNS>:/mnt
```
Run the thing
```
python3 cxg_audit.py &
```
OR
```
python3 migrator.py &
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
Get the instance ID and terminate the instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-04b6c97b14c54de18 --query Reservations[*].Instances[*].[InstanceId] --output text
```
```
aws ec2 terminate-instances --instance-ids <instance ID>
```
