import sys
import argparse
import subprocess
from dataclasses import dataclass, astuple, is_dataclass
import csv
import boto3
from botocore.exceptions import ClientError, BotoCoreError
import logging
import os.path
import shutil
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
 - Aspera: https://www.ibm.com/support/fixcentral/swg/selectFixes?parent=ibm%7EOther%20software&product=ibm/Other+software/IBM+Aspera+High-Speed+Transfer+Server&release=4.4.7&platform=All&function=all

Logging will be sent to Uploadosaur.log

Files will be downloaded to the same directory that Uploadosaur is run from

To Run:
python Uploadosaur.py -k [Full path to Aspera Key file] -a [account folder ex. uploads/jlz_stanford.edu_xxxxxx/]  -f [name of final destination folder within the account folder ex. upload_folder_1] -l [List of S3_Links]

Please direct complaints about Uploadosaur to anyone except Jim Chaffer and Brian Mott
"""

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--key_path',
        '-k',
        dest='key_path',
        help='FULL path to Aspera openssh key'
    )
    parser.add_argument(
        '--account',
        '-a',
        dest='account',
        help='Account folder for aspera ex. uploads/jlz_stanford.edu_xxxxxxx/'
    )
    parser.add_argument(
        '--folder',
        '-f',
        dest='folder',
        help='Name of upload destination folder ex. upload_folder_1'
    )
    parser.add_argument(
        '--list',
        '-l',
        dest='input_file_list',
        help='CSV file list containing full URI S3 paths of files'
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args

@dataclass
class AsperaUploadInfo:
    """
    Dataclass for storing server connection information for connecting to FTP server via Aspera
    """
    key_path: str
    account: str
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
    with open(S3_list, encoding='utf-8', newline='') as f:
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
                            if sorted_file_list[file_index_1] not in [pair.file_2.S3_Path for pair in paired_files]:
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
                file.split_amount = -(-file.size // split_amount)
    return file_list
        


async def splitter(original_file, download_folder, split_folder, sem):
    """
    Asynchronous function that splits fastq file according to split amount set in split_setter()
    Returns list of new split file names, and removes original file to save space
    """
    print(f'Starting split of {original_file.file_name}')
    logging.info(f'Starting split of {original_file.file_name}')
    try:
        if not os.path.exists(split_folder):
            os.makedirs(split_folder)
        shutil.move(download_folder + original_file.file_name, split_folder)
        logging.info(f"File '{original_file.file_name}' moved successfully from '{download_folder}' to '{split_folder}'")
    except FileNotFoundError:
        print(f"Error: Original file '{original_file.file_name}' not found in '{download_folder}'")
        logging.info(f"Error: Original file '{original_file.file_name}' not found in '{download_folder}'")
    except Exception as e:
        print(f"An error occurred in moving '{original_file.file_name}' from download folder ",
        f"'{download_folder}' to split folder '{split_folder}': {e}")
        logging.info(f"An error occurred in moving '{original_file.file_name}' from download folder ",
        f"'{download_folder}' to split folder '{split_folder}': {e}")
    output_files = [(download_folder + original_file.file_name.split('.')[0] + '_split' + str(i+1) + '.fastq.gz') for i in range(0,original_file.split_amount)]
    split_command = ['fastqsplitter','-c', '7', '-t', '3', '-i', split_folder + original_file.file_name]
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
            return output_files
        else:
            print(f'Failed to split {original_file.file_name}')
            print(f'Error output:\n{stderr.decode().strip()}')
            logging.error(f'Failed to split {original_file.file_name}')
            logging.error(f'Error output:\n{stderr.decode().strip()}')
            return original_file.file_name


async def downloader(file, download_folder, sem):
    """
    Asynchronous function that downloads file from S3, checking that downloaded file size matches expected
    Returns same SingleFastQFile object that was input if download was successful, None if not.

    Prior to downloading, checks for presence of file locally, if present skips.
    """
    
    async with sem:
        if ((os.path.exists(download_folder + file.file_name)) and (os.path.getsize(download_folder + file.file_name)) == file.size):
            logging.info(f'{file.file_name} already present in {download_folder}, skipping download')
            print(f'{file.file_name} already present in {download_folder}, skipping download')
            return file
        else:
            logging.info(f'Starting download of {file.file_name} to {download_folder}')
            print(f'Starting download of {file.file_name} to {download_folder}')
            session = aioboto3.Session()
            bucket_name = file.S3_Path.split('/',1)[0]
            key = file.S3_Path.split('/',1)[1]
            async with session.client("s3") as s3_client:
                try:
                    await s3_client.download_file(bucket_name, key, (download_folder + file.file_name))
                    size_result = subprocess.run(['ls', '-l', download_folder + file.file_name],capture_output=True,text=True, check=True)
                    if int(size_result.stdout.split()[4]) == int(file.size):
                        logging.info(f"Download Completed successfully: {'s3://' + file.S3_Path} -> {download_folder + file.file_name}")
                        print(f"Download Completed successfully: {'s3://' + file.S3_Path} -> {download_folder + file.file_name}")
                        return file
                    else:
                        logging.error(f'File sizes between S3 {file.size} and ',
                        f'local {size_result.stdout.split()[4]} do not match for file {file.S3_Path}')
                        print(f'File sizes between S3 {file.size} and local ',
                        f'{size_result.stdout.split()[4]} do not match for file {file.S3_Path}')
                        return file.file_name
                except ClientError as e:
                    logging.error(f'An AWS error occured for file {file.S3_Path} during download: {e.response.get("Error", {}).get("Code")}')
                    print(f'An AWS error occured for file {file.S3_Path} during download: {e.response.get("Error", {}).get("Code")}')
                    return file.file_name
                except BotoCoreError as e:
                    logging.error(f'BotoCoreError for file {file.S3_Path} during download: {e}')
                    print(f'BotoCoreError for file {file.S3_Path} during download: {e}')
                    return file.file_name
                except Exception as e:
                    logging.error(f'Some sort of AWS Error occured for file {file.S3_Path} during download: {e}')
                    print(f'Some sort of AWS Error occured for file {file.S3_Path} during download: {e}')
                    return file.file_name


def uploader(ftp_server_info, download_folder):
    """
    Synchronous function for uploading FastQ file directory using Aspera
    Uploads all files in download folder directory to same name directory made on Aspera
    Does not delete local directory after upload
    """
    logging.info(f'Starting upload of all files in {download_folder}')
    print(f'Starting upload of all files in {download_folder}')
    full_aspera_dest = 'subasp@upload.ncbi.nlm.nih.gov:' + ftp_server_info.account
    upload_complete = False
    while upload_complete != True:
        try:
            upload = subprocess.run(['/home/jovyan/.aspera/connect/bin/ascp', '-i', ftp_server_info.key_path, '-QT', 
                                '-l', '600m', '-k1', '-d', download_folder, full_aspera_dest]) 
            if upload.returncode == 0:
                logging.info(f'The upload of all files has completed')
                upload_complete = True
                print(f'The upload of all files has completed')
            else:
                logging.info(f'An error occured during the upload, restarting')
                print(f'An error occured during the upload, restarting')
                
        except Exception as e:
            logging.info(f'An error occured when trying to upload \n Further info: {e}')
            print(f'An error occured when trying to upload \n Further info: {e}')
            


async def Path_Setter(Full_File_List, download_folder):
    """
    Main asynchronous handling function. 
    Awaits completion of various coprocesses, upon completion, sends results to next step in workflow. 
    Files are downloaded, then either split and uploaded or just directly uploaded.
    Semaphores are used to limit number of coprocesses occuring for any step.
    """
    download_sem = asyncio.Semaphore(5) # How many download processes can occur at once
    download_folder = './' + download_folder + '/'
    split_folder = './Need_Splitting/'
    os.makedirs(download_folder, exist_ok=True)
    split_sem = asyncio.Semaphore(3) # How many splitting processes can occur at once
    split_tasks = []
    failed_downloads = []
    failed_splits = []
    download_tasks = [asyncio.create_task(downloader(file, download_folder, download_sem)) for file in Full_File_List]
    for completed_download in asyncio.as_completed(download_tasks):
        downloaded_file = await completed_download
        if is_dataclass(downloaded_file) != True:
            logging.error('No downloaded file returned to path setter. Nothing to pass to upload or splitter')
            print(f'Some sort of issue with file, not passing to upload or splitter')
            failed_downloads.append(downloaded_file)
            pass
        elif downloaded_file.split_amount == None:
            pass
        else:
            split_tasks.append(asyncio.create_task(splitter(downloaded_file, download_folder, split_folder, split_sem)))

    for completed_split in asyncio.as_completed(split_tasks):
        split_files = await completed_split
        if not is_instance(split_files, list):
            logging.error(f'An error occured while splitting {split_files}')
            print(f'An error occured while splitting {split_files}')
            failed_splits.append(split_files)
            pass
    return failed_downloads, failed_splits, download_folder


def main(aspera_info, csv_file):
    logging.basicConfig(
        filename='Uploadosaur.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info('STARTING NEW RUN OF UPLOADOSAUR')
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
    failed_downloads, failed_splits, download_folder = asyncio.run(Path_Setter(Full_File_List, aspera_info.folder))
    if failed_downloads or failed_splits:
        if failed_downloads:
            print('The following files had download failures: \n')
            logging.info('The following files had download failures: \n')
            for f in failed_downloads:
                print(f'{f} \n')
                logging.info(f'{f} \n')
        if failed_splits:
            print('The following files failed during splitting: \n')
            logging.info('The following files failed during splitting: \n')
            for f in failed_splits:
                print(f'{f} \n')
                logging.info(f'{f} \n')
        while True:
            user_choice = input('With this in mind, do you still want to upload? (y/n): ').lower()
            logging.info('With this in mind, do you still want to upload? (y/n): ')
            logging.info(user_choice)
            if user_choice == 'y':
                print('Proceeding to upload')
                logging.info('Proceeding to upload')
                uploader(aspera_info, download_folder)
                break 
            elif user_choice == 'n':
                print('Not proceeding with upload, terminating.')
                logging.info('Not proceeding with upload, terminating')
                break
            else:
                print_help("Invalid input. Please enter 'y' or 'n'.")
                logging.info("Invalid input. Please enter 'y' or 'n'.")
    else:
        while True:
            user_choice = input('No download or split failures, proceed with upload? (y/n): ').lower()
            logging.info('No download or split failures, proceed with upload? (y/n): ')
            if user_choice == 'y':
                print("Proceeding to upload")
                logging.info('Proceeding to upload')
                uploader(aspera_info, download_folder)
                break 
            elif user_choice == 'n':
                print("Not proceeding with upload, terminating.")
                logging.info('Not proceeding with upload, terminating')
                break
            else:
                print_help("Invalid input. Please enter 'y' or 'n'.")
                logging.info("Invalid input. Please enter 'y' or 'n'.")
        
    


if __name__ == '__main__':
    args = getArgs()
    aspera_info = AsperaUploadInfo(
        key_path = args.key_path,
        account = args.account,
        folder = args.folder
        )
    main(aspera_info,args.input_file_list)