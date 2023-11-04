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
source $HOME/.cargo/env ; sudo yum group install "Development Tools" ; wget https://cmake.org/files/v3.12/cmake-3.12.3.tar.gz ; tar zxvf cmake-3.12.3.tar.gz ; cd cmake-3.12.3 ; ./bootstrap --prefix=/usr/local ; make -j$(nproc) ; sudo make install
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
aws s3 cp s3://submissions-lattice-sra/PRJNA1020837/SRR26173958/C70_SC_5pr_possorted_genome_bam.bam ./
```
Run the conversion
```
target/release/bamtofastq --reads-per-fastq=500000000 C70_SC_5pr_possorted_genome_bam.bam ./fastqs
```
Each set of bam files requires thorough review as they vary may in terms of how many fastqs are produced, how they are named, and what the output file structure is. What follows outlines the steps to take, but the specific code may need to be adjusted.

Likely need to change the name of the output fastqs so they'll be unique file names
```
cd fastqs/* ; for a in `ls *` ; do mv $a SRR26173958_$a ; done ; cd .. ; cd ..
```
...before copying them to s3
```
aws s3 cp --recursive fastqs/ s3://submissions-lattice-sra/PRJNA1020837/SRR26173958/
```
Remove the bam & the fastq directory before doing the next one
```
rm C70_SC_5pr_possorted_genome_bam.bam ; rm -r fastqs
```
[Command templates](https://docs.google.com/spreadsheets/d/1M0bNaRoiYwWBZFkaqZABY6V65UHetOlbZ0g7KUzmg-c/edit#gid=1217473508) for doing many conversions