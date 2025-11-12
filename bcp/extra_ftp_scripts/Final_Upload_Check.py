import sys
import argparse
import subprocess
from dataclasses import dataclass, astuple
import csv
import boto3
from botocore.exceptions import ClientError, BotoCoreError
import logging
import os.path
import ftplib

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--address',
        '-a',
        dest='address',
        help='FTP address'
    )
    parser.add_argument(
        '--username',
        '-u',
        dest='username',
        help='Username for FTP'
    )
    parser.add_argument(
        '--password',
        '-p',
        dest='password',
        help='Password for FTP'
    )
    parser.add_argument(
        '--folder',
        '-f',
        dest='folder',
        help='New folder to put files in'
    )
    parser.add_argument(
        '--list',
        '-l',
        dest='input_file_list',
        help='CSV file list'
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


@dataclass
class FTPUploadInfo:
    """
    Dataclass for storing server connection information for connecting to FTP server via ncftpput
    """
    address: str
    username: str
    password: str
    folder: str

@dataclass
class SingleFastQFile:
    """
    Dataclass for storing information about a single fastq file
    """
    
    S3_Path: str
    file_name: str
    size: str

def get_file_info(S3_Link):
    """
    Function for gathering file information, most importantly size, from S3
    """
    try:
        completed_process_s3 = subprocess.run(['aws','s3','ls', S3_Link],capture_output=True,text=True, check=True)
        result_s3 = completed_process_s3.stdout.split()
        return SingleFastQFile(
            S3_Path = S3_Link,
            file_name = result_s3[3],
            size = result_s3[2]
        )
    except Exception as e:
        logging.error(f'An error occured while trying to get file info for {S3_Link}, stopping run: {e}')
        print(f'An error occured while trying to get file info for {S3_Link}, stopping run: {e}')
        sys.exit()


def get_S3_sizes(S3_list):
    """
    Function that takes in list of S3 file links and gets sizes

    Outputs list of SingleFastQFile
    """
    file_list = []
    with open(S3_list, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            file_list.extend(row)
        
    single_files = []
    sorted_file_list = sorted(file_list)

    for file_index_1 in range(len(sorted_file_list)):
        single_files.append(get_file_info(sorted_file_list[file_index_1]))
    return single_files


def main(ftp_server_info, original_list):
    '''
    Takes in full list of S3 URIs for files wanted to be uploaded, and server info. Checks what files have been uploaded fully.
    '''
    print('Starting Final Check of uploaded files:')
    files_with_sizes = get_S3_sizes(original_list) # need to parse out file names?
    confirmed_uploaded = []
    ftp_uploaded = []
    try:
        with ftplib.FTP(ftp_server_info.address) as ftp:
            ftp.login(user=ftp_server_info.username, passwd=ftp_server_info.password)
            ftp.cwd(ftp_server_info.folder)
            dir_list = []
            ftp.dir(dir_list.append)
            for line in dir_list:
                ftp_uploaded.append(line[29:].strip().split(' '))
    except ftplib.all_errors as e:
        print(f'Error occured obtaining size of file {file.file_name} on FTP: \n {e}')
        logging.info(f'Error occured obtaining size of file {file.file_name} on FTP: \n {e}')
    for file in files_with_sizes:
        for line in ftp_uploaded:
            if file.file_name == line[-1]:
                if file.size == line[3]:
                    print('Confirmed Uploaded: ' + file.S3_Path)
                    confirmed_uploaded.append(file)
                else:
                    print('File sizes do not match')
                    print('On S3 ' + file.file_name + 'is ' + file.size + 'bytes')
                    print('On FTP ' + line[-1] + 'is ' + line[3] + 'bytes')
            else:
                print('File not present on FTP: ' + file.file_name)
    with open('unuploaded.csv', 'w', newline='') as csvfile:
        # Create a CSV writer object
        writer = csv.writer(csvfile)
        for file in files_with_sizes:
            if file not in confirmed_uploaded:
                writer.writerows([[file.file_name]])
    with open('unuploaded_full_uris.csv', 'w', newline='') as csvfile:
        # Create a CSV writer object
        writer = csv.writer(csvfile)
        for file in files_with_sizes:
            if file not in confirmed_uploaded:
                writer.writerows([[file.S3_Path]])
    with open('confirmed_uploaded.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for file in files_with_sizes:
            if file in confirmed_uploaded:
                writer.writerows([[file.file_name]])
    with open('confirmed_uploaded_full_uris.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for file in files_with_sizes:
            if file in confirmed_uploaded:
                writer.writerows([[file.S3_Path]])


if __name__ == '__main__':
    args = getArgs()
    ftp_info = FTPUploadInfo(
        address = args.address,
        username = args.username,
        password = args.password,
        folder = args.folder)
    main(ftp_info,args.input_file_list)


