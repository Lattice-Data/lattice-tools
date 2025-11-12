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
        '--current',
        '-c',
        dest='current',
        help='Folder files are currently in'
    )
    parser.add_argument(
        '--desired',
        '-d',
        dest='desired',
        help='Folder files are desired to be in'
    )
    parser.add_argument(
        '--list',
        '-l',
        dest='input_file_list',
        help='list of files to move'
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
    current: str
    desired: str



def main(ftp_server_info, move_list):
    '''
    Takes in full list of files desired to be moved to new directory in FTP and current and desired directory names
    Moves files to desired directory
    '''
    print(f'Moving files from {ftp_server_info.current} to {ftp_server_info.desired}')
    file_list = []
    with open(move_list, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            file_list.extend(row)
    try:
        with ftplib.FTP(ftp_server_info.address) as ftp:
            ftp.login(user=ftp_server_info.username, passwd=ftp_server_info.password)
            for file in file_list:
                    location = ftp_server_info.current + file
                    destination = ftp_server_info.desired + file
                    ftp.rename(location,destination)
    except ftplib.all_errors as e:
        print(f'Error occured: {e}')
        logging.info(f'Error occured: {e}')


if __name__ == '__main__':
    args = getArgs()
    ftp_info = FTPUploadInfo(
        address = args.address,
        username = args.username,
        password = args.password,
        current = args.current,
        desired = args.desired)
    main(ftp_info,args.input_file_list)


