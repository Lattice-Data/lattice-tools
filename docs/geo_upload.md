# Submitting data to GEO
This document describes the process of taking raw matrix and fastq files for a Dataset in Lattice and submitting them to GEO along with the associated metadata.

Create metadata submission spreadsheet
----------------
- Run geo_metadata.py for the desired Dataset with the following commandline:
```
python geo_metadata.py -m local -d <Dataset_accession>
```
	This will create three files
	- <Dataset_accession>_metadata.csv: This is a csv of all the study, sample, protocol, and fasta file metadata.
	- <Dataset_accession>_md5sum.csv: This is a csv of all the fasta and raw matrix file names and checksums.
	- <Dataset_accession>_s3_uri.csv: This is a list of the s3 URI for all fastq and raw matrix files.
- Upload the metadata and md5sum csv files in a google spreadsheet. Have wrangler add any additional contributor or protocol information.
- Make sure all raw matrix and fastq files listed have unique files names.

Launch EC2
----------------
- Launch an EC2 instance with enough disk space for 2 upload processes running in parallel.
- Install c compiler:
```
sudo yum install gcc
```
- Install ncftp on machine:
```
wget ftp://ftp.ncftp.com/ncftp/ncftp-3.2.6-src.tar.gz
tar -zxvf ncftp-3.2.6-src.tar.gz
cd ncftp-3.2.6
./configure
make
make install
make clean
```

Upload fastq files
----------------
- Make sure all raw matrix and fastq files listed have unique files names.
- Have the wrangler obtain user, password, and directory from GEO. Create subdirectories for raw matrix files and raw sequence files.
- Create a shell script upload_fastq.sh to download file from s3, upload to GEO, and remove the file locally with following example:
```
#!/bin/bash

# Open file of s3 uri

while read URI
do
        IFS='/' read -ra uri_array <<< "$URI"
        aws s3 cp $URI .
        array_length=${#uri_array[@]}
        filename_loc=$((array_length - 1))
        filename="${uri_array[$filename_loc]}"
        ncftpput -B 33554432 -u <user> -p <password> ftp-private.ncbi.nlm.nih.gov <directory> $filename
        rm $filename
done < $1
```
- To be make the script executable, run:
```
chmod a+x upload_fastq.sh
```
- Subset the s3_uri csv files into two files, and run the shell script in the background in multiple screen windows on the EC2:
```
screen -S geo
./upload_fastq.sh part1_s3_uri.csv &
```
- Check the files on GEO using ncftpls:
```
ncftpls -u <user> -p <password> ftp://ftp-private.ncbi.nlm.nih.gov/uploads/<remainder_of_path>/
```


