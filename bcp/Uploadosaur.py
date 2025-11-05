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
try:
    import asyncio
except ImportError:
    print('asyncio not installed')
try:
    import aioboto3
except ImportError:
    print('aioboto3 not installed')

"""
Uploadosaur expects the following packages to be installed prior to running:
 - fastqsplitter: https://fastqsplitter.readthedocs.io/en/stable/#fastqsplitter
 - ncftp: https://www.ncftp.com/ncftp/

Logging will be sent to Uploadosaur.log

Files will be downloaded to the same directory that Uploadosaur is run from

To Run:
python Uploadosaur.py -a [Address of FTP Server] -u [Username] -p [Password] -f [FTP Upload Folder] -l [List of S3_Links]

Please direct complaints about Uploadosaur to anyone except Jim Chaffer and Brian Mott
"""

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
    size: int
    split_amount: int | None = None

@dataclass
class PairedFastQFiles:
    """
    Dataclass for keeping paired read fastq files together for split amount setting
    """
    file_1: SingleFastQFile
    file_2: SingleFastQFile
    

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
            size = int(result_s3[2])
        )
    except Exception as e:
        logging.error(f'An error occured while trying to get file info for {S3_Link}, stopping run: {e}')
        print(f'An error occured while trying to get file info for {S3_Link}, stopping run: {e}')
        sys.exit()
        


def parse_file_list(S3_list):
    """
    Function that takes in list of S3 file links and determines which files are paired and which are not.

    To determine which files are paired, checks if file names are the same except for one character,
    and that character before the different character is 'R' signifying 'Read'.

    Outputs list of SingleFastQFile and list of PairedFastQFile dataclasses
    """
    file_list = []
    with open(S3_list, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            file_list.extend(row)
        
    single_files = []
    paired_files = []
    sorted_file_list = sorted(file_list)

    for file_index_1 in range(len(sorted_file_list)):
        for file_index_2 in range(file_index_1 + 1,len(sorted_file_list)):
            file_1 = sorted_file_list[file_index_1]
            file_2 = sorted_file_list[file_index_2]
            if len(file_1) == len(file_2):
                diff_count = 0
                for i in range(len(file_1)):
                    if file_1[i] != file_2[i]:
                        char_before = file_1[i-1]
                        diff_count += 1
                        if diff_count > 1 and file_index_2 == (len(sorted_file_list) - 1):
                            if sorted_file_list[file_index_1] not in paired_files:
                                single_files.append(get_file_info(sorted_file_list[file_index_1]))
                                break
                        elif diff_count > 1:
                            break
                if diff_count == 1 and char_before == 'R':
                    paired_files.append(PairedFastQFiles(
                        file_1 = get_file_info(sorted_file_list[file_index_1]),
                        file_2 = get_file_info(sorted_file_list[file_index_2])
                    ))
    return single_files, paired_files


def split_setter(file_list, paired = False):
    """
    Function that determines the split amount for either a SingleFastQFile, or a pair of SingleFastQFiles
    If a pair, sets split amount to be the same between the two SingleFastQFiles so split the same
    """
    split_amount = 100000000000 # 100gb in bytes
    if paired:
        for pair in file_list:
            if pair.file_1.size >= split_amount or pair.file_2.size >= split_amount:
                if pair.file_1.size >= pair.file_2.size:
                    pair.file_1.split_amount = -(-pair.file_1.size // split_amount)
                    pair.file_2.split_amount = -(-pair.file_1.size // split_amount)
                else:
                    pair.file_1.split_amount = -(-pair.file_2.size // split_amount)
                    pair.file_2.split_amount = -(-pair.file_2.size // split_amount)
    else:
        for file in file_list:
            if file.size >= split_amount:
                print(file.file_name)
                file.split_amount = -(-file.size // split_amount)
    return file_list
        


async def splitter(original_file, sem):
    """
    Asynchronous function that splits fastq file according to split amount set in split_setter()
    Returns list of new split file names, and removes original file to save space
    """
    print(f'Starting split of {original_file.file_name}')
    logging.info(f'Starting split of {original_file.file_name}')
    output_files = [(original_file.file_name.split('.')[0] + '_split' + str(i+1) + '.fastq.gz') for i in range(0,original_file.split_amount)]
    split_command = ['fastqsplitter','-c', '7', '-t', '3', '-i', original_file.file_name]
    for o in output_files:
        split_command.append('-o')
        split_command.append(o)
    logging.info(f'Full split command: {split_command}')
    async with sem:
        split = await asyncio.create_subprocess_exec(
            *split_command,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
    
        stdout, stderr = await split.communicate()
        return_code = split.returncode
    
        if return_code == 0:
            print(f'{original_file.file_name} split successfully')
            logging.info(f'{original_file.file_name} split successfully')
            await asyncio.create_subprocess_exec('rm','./'+original_file.file_name)
            logging.info(f'{original_file.file_name} removed from local directory')
            print(f'{original_file.file_name} removed from local directory')
            return output_files
        else:
            print(f'Failed to split {original_file.file_name}')
            print(f'Error output:\n{stderr.decode().strip()}')
            logging.error(f'Failed to split {original_file.file_name}')
            logging.error(f'Error output:\n{stderr.decode().strip()}')
            return None


async def downloader(file, sem):
    """
    Asynchronous function that downloads file from S3, checking that downloaded file size matches expected
    Returns same SingleFastQFile object that was input if download was successful, None if not.

    Prior to downloading, checks for presence of file locally, if present sends directly to upload.
    """
    
    async with sem:
        if ((os.path.exists(file.file_name)) and (os.path.getsize(file.file_name)) == file.file_size):
            logging.info(f'{file.file_name} already present locally, skipping download')
            print(f'{file.file_name} already present locally, skipping download')
            return file
        else:
            logging.info(f'Starting download of {file.file_name}')
            print(f'Starting download of {file.file_name}')
            session = aioboto3.Session()
            bucket_name = file.S3_Path.split('/',1)[0]
            key = file.S3_Path.split('/',1)[1]
            async with session.client("s3") as s3_client:
                try:
                    await s3_client.download_file(bucket_name, key, ('./' + file.file_name))
                    size_result = subprocess.run(['ls', '-l', file.file_name],capture_output=True,text=True, check=True) # Might need to improve this
                    if int(size_result.stdout.split()[4]) == int(file.size):
                        logging.info(f"Download Completed successfully: {'s3://' + file.S3_Path} -> {'./' + file.file_name}")
                        print(f"Download Completed successfully: {'s3://' + file.S3_Path} -> {'./' + file.file_name}")
                        return file
                    else:
                        logging.error(f'File sizes between S3 {file.size} and local {size_result.stdout.split()[4]} do not match for file {file.S3_Path}')
                        print(f'File sizes between S3 {file.size} and local {size_result.stdout.split()[4]} do not match for file {file.S3_Path}')
                        return None
                except ClientError as e:
                    logging.error(f'An AWS error occured for file {file.S3_Path} during download: {e.response.get("Error", {}).get("Code")}')
                    print(f'An AWS error occured for file {file.S3_Path} during download: {e.response.get("Error", {}).get("Code")}')
                    return None
                except BotoCoreError as e:
                    logging.error(f'BotoCoreError for file {file.S3_Path} during download: {e}')
                    print(f'BotoCoreError for file {file.S3_Path} during download: {e}')
                    return None
                except Exception as e:
                    logging.error(f'Some sort of AWS Error occured for file {file.S3_Path} during download: {e}')
                    print(f'Some sort of AWS Error occured for file {file.S3_Path} during download: {e}')
                    return None


async def uploader(ftp_server_info, file_name, sem):
    """
    Asynchronous function for uploading FastQ file to FTP server using ncftpput
    Removes local copy of file upon successful upload
    """
    async with sem:
        logging.info(f'{file_name} sent to upload')
        print(f'{file_name} sent to upload')
        try:
            upload = await asyncio.create_subprocess_exec(
                'ncftpput','-z', '-V', '-B', '33554432', '-u', ftp_server_info.username, '-p', ftp_server_info.password, 
                ftp_server_info.address, ftp_server_info.folder, file_name)
        
            await upload.wait()
            await asyncio.create_subprocess_exec('rm','./'+file_name)
            logging.info(f'{file_name} appears to have successfully uploaded, file removed from local')
            print(f'{file_name} appears to have successfully uploaded, file removed from local')
        except Exception as e:
            logging.info(f'An error occured uploading {file_name} to ftp. Further info: {e}')
            print(f'An error occured uploading {file_name} to ftp. Further info: {e}')
            


async def Path_Setter(ftp_server_info, Full_File_List):
    """
    Main asynchronous handling function. Awaits completion of various coprocesses, upon completion, sends
    results to next step in workflow. Files are downloaded, then either split and uploaded or just directly uploaded.
    Semaphores are used to limit number of coprocesses occuring for any step.
    """
    download_sem = asyncio.Semaphore(10) # How many download processes can occur at once
    upload_sem = asyncio.Semaphore(5) # How many upload processes can occur at once
    split_sem = asyncio.Semaphore(3) # How many splitting processes can occur at once
    split_tasks = []
    upload_tasks = []
    download_tasks = [asyncio.create_task(downloader(file, download_sem)) for file in Full_File_List]
    for completed_download in asyncio.as_completed(download_tasks):
        downloaded_file = await completed_download
        if downloaded_file == None:
            logging.error('No downloaded file returned to path setter. Nothing to pass to upload or splitter')
            print(f'Some sort of issue with file, not passing to upload or splitter')
            pass
        elif downloaded_file.split_amount == None:
            upload_tasks.append(asyncio.create_task(uploader(ftp_server_info, downloaded_file.file_name, upload_sem)))
        else:
            split_tasks.append(asyncio.create_task(splitter(downloaded_file, split_sem)))

    for completed_split in asyncio.as_completed(split_tasks):
        split_files = await completed_split
        if split_files == None:
            logging.error(f'An error occured while splitting {file_to_delete.file_name}')
            print(f'An error occured while splitting {file_to_delete.file_name}')
            pass
        for split_file in split_files:
            upload_tasks.append(asyncio.create_task(uploader(ftp_server_info, split_file, upload_sem)))

    for completed_upload in asyncio.as_completed(upload_tasks):
        await completed_upload


def main(ftp_server_info, csv_file):
    logging.basicConfig(
        filename='Uploadosaur.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    Full_File_List = []
    Single_File_List, Paired_File_List = parse_file_list(csv_file)
    if Single_File_List:
        Single_File_List = split_setter(Single_File_List)
        Full_File_List.extend(Single_File_List)
    if Paired_File_List:
        Paired_File_List = split_setter(Paired_File_List, paired=True)
        for pair in Paired_File_List:
            Full_File_List.append(pair.file_1)
            Full_File_List.append(pair.file_2)# Breaking up pairs now that have same split amounts
    asyncio.run(Path_Setter(ftp_server_info, Full_File_List))


if __name__ == '__main__':
    args = getArgs()
    ftp_info = FTPUploadInfo(
        address = args.address,
        username = args.username,
        password = args.password,
        folder = args.folder)
    main(ftp_info,args.input_file_list)