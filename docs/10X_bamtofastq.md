Running 10X bamtofastq on an EC2 instance
---------------- 
Create a new spot instance from the bamtofastq template
```
aws ec2 run-instances --launch-template LaunchTemplateName=bamtofastq
```
Identify the Public DNS for the instance and ssh to instance
```
aws ec2 describe-instances --filters Name=image-id,Values=ami-070c8ca4ac77fae0b --query 'Reservations[*].Instances[*].[PublicDnsName]' --output text
```
```
ssh -i lattice_ec2.pem ec2-user@<Public DNS>
```
Add on the extra storage and set AWS credentials for copying to/from Lattice S3
```
sudo mkfs.ext4 /dev/nvme1n1 ; sudo mount /dev/nvme1n1 /mnt ; sudo chown ec2-user /mnt ; cd /mnt ; export AWS_SHARED_CREDENTIALS_FILE=/mnt/credentials
```
Install rust (#1 to proceed)
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh 
```
...and other libraries
```
source $HOME/.cargo/env ; sudo yum group install "Development Tools" ; wget https://cmake.org/files/v3.12/cmake-3.12.3.tar.gz ; tar zxvf cmake-3.* ; cd cmake-3.* ; ./bootstrap --prefix=/usr/local ; make -j$(nproc) ; sudo make install
```
Install the bamtofastq tool
```
cd .. ; git clone https://github.com/10XGenomics/bamtofastq.git ; cd bamtofastq/ ; cargo build --release
```
From your local machine, copy up your credentials
```
scp -i lattice_ec2.pem ~/.aws/credentials ec2-user@<Public DNS>:/mnt
```
Copy the first bam
```
aws s3 cp s3://submissions-lattice-sra/SRR10080331/immvarYE_0907_4.bam ./
```
Run the conversion
```
target/release/bamtofastq --reads-per-fastq=500000000 immvarYE_0907_4.bam ./fastqs
```
Likely need to change the name of the output fastqs so they'll be unique file names
```
cd fastqs/* ; for a in `ls *` ; do mv $a immvarYE_0907_4_$a ; done ; cd .. ; cd ..
```
Before copying them to s3
```
aws s3 cp --recursive fastqs/ s3://submissions-lattice-sra/SRR16227577/
```
Remove the bam & the fastq directory before doing the next one
```
rm C70_TST_possorted_genome_bam.bam ; rm -r fastqs
```
Command templates for doing many conversions
https://docs.google.com/spreadsheets/d/1M0bNaRoiYwWBZFkaqZABY6V65UHetOlbZ0g7KUzmg-c/edit#gid=1217473508