# Submitting data to GEO
This document describes the process of taking raw matrix and fastq files for a Dataset in Lattice and submitting them to GEO along with the associated metadata.

Create metadata submission spreadsheet
----------------
- Run geo_metadata.py for the desired Dataset with the following commandline:
```
python geo_metadata.py -m local -d <Dataset_accession>
```
- This will create three files:
	- <Dataset_accession>_metadata.csv: This is a csv of all the study, sample, protocol, and fasta file metadata.
	- <Dataset_accession>_md5sum.csv: This is a csv of all the fasta and raw matrix file names and checksums.
	- <Dataset_accession>_s3_uri.csv: This is a list of the s3 URI for all fastq and raw matrix files.
- Upload the metadata and md5sum csv files in a google spreadsheet as two separate tabs. Have wrangler add any additional contributor or protocol information.

Launch EC2
----------------
- Launch an EC2 instance with enough disk space for at least 2 upload processes running in parallel.
- Install c compiler:
```
sudo yum install gcc
```
- Install ncftp on machine:
```
wget https://www.ncftp.com/downloads/ncftp/binaries/ncftp-3.2.6-linux-x86_64-glibc2.17-export.tar.gz; tar -zxvf ncftp-3.2.6-linux-x86_64-glibc2.17-export.tar.gz; cd ncftp-3.2.6/; sudo make install
```

Upload processed files and metadata spreadsheet
----------------
- Make sure all raw matrix files listed have unique files names.
- Have the wrangler create a GEO profile, if needed, and start a new submission. Obtain the user, password, and directory from the GEO account information. Create subdirectories for raw matrix files and raw sequence files using the ncftp function:
```
ncftp
```
```
open ftp://geoftp:<password>@ftp-private.ncbi.nlm.nih.gov
```
```
cd uploads/<directory>
```
```
mkdir raw_data
```
```
mkdir processed_data
```
- Create a shell script upload_fastq.sh to download file from s3, upload to GEO, and remove the file locally with following example, updating the "user", "password", and "directory" variables:
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
        ncftpput -B 33554432 -u <user> -p <password> ftp-private.ncbi.nlm.nih.gov /uploads/<directory>/processed_data/ $filename
        rm $filename
done < $1
```
- To be make the script executable, run:
```
chmod a+x upload_fastq.sh
```
- Subset the s3_uri.csv file for matrix files, and run the bash script to upload to the GEO server. Replace "h5$" with "h5ad$" if working with h5ad files:
```
grep 'h5$' <Dataset_accession>_s3_uri.csv > h5_s3_uri.csv
```
```
./upload_fastq.sh h5_s3_uri.csv &
```
- Download metadata google spreadsheet as excel spreadsheet. Upload to GEO server using ncftpput, replacing the "user", "password", "directory", and "excel_filename" variables:
```
ncftpput -B 33554432 -u <user> -p <password> ftp-private.ncbi.nlm.nih.gov /uploads/<directory>/processed_data/ <excel_filename>
```

Upload fastq files
----------------
- Make sure all fastq files listed have unique files names.
- To run the long processes of uploading to the fastq file to the GEO server, create a 'screen' session, and reattach when needed:
```
screen -S geo
```
```
screen -r geo
```
- Update the shell script using "raw_data" as part of the directory path, replacing "user", "password", and "directory" variables:
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
        ncftpput -B 33554432 -u <user> -p <password> ftp-private.ncbi.nlm.nih.gov /uploads/<directory>/raw_data/ $filename
        rm $filename
done < $1
```
- Subset the s3_uri csv to only fastq files and split into at least two files. Run the shell script in the background in multiple screen windows on the EC2:
```
./upload_fastq.sh part1_s3_uri.csv &
```
- Check the files on GEO using ncftpls:
```
ncftpls -u <user> -p <password> ftp://ftp-private.ncbi.nlm.nih.gov/uploads/<directory>/raw_data/
```


