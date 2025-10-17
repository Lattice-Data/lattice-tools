import sys
import argparse
import subprocess
from dataclasses import dataclass, astuple
import csv
import asyncio

#subprocess.run(['ls','-l'])

#subprocess.run(['aws','s3','cp', filepath])

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--address",
        "-a",
        dest='address',
        help="FTP address"
    )
    parser.add_argument(
        "--username",
        "-u",
        dest='username',
        help="Username for FTP"
    )
    parser.add_argument(
        "--password",
        "-p",
        dest='password',
        help="Password for FTP"
    )
    parser.add_argument(
        "--folder",
        "-f",
        dest='folder',
        help="New folder to put files in"
    )
    parser.add_argument(
        "--list",
        "-l",
        dest='input_file_list',
        help="CSV? file list"
    )
    #? new folder to make and put files in?
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args

@dataclass
class FTPUploadInfo:
    address: str
    username: str
    password: str
    folder: str
    
@dataclass
class SingleFastQFile:
    S3_Path: str
    file_name: str
    size: int
    split_amount: int | None = None

@dataclass
class PairedFastQFiles:
    file_1: SingleFastQFile
    file_2: SingleFastQFile


async def run(cmd: str):
    # From async documentation
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stderr=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE
    )

    stdout, stderr = await proc.communicate()
    if stdout:
        return stdout

    #print(f'[{cmd!r} exited with {proc.returncode}]')
    #if stdout:
    #    print(f'[stdout]\n{stdout.decode()}')
    #if stderr:
    #    print(f'[stderr]\n{stderr.decode()}')

def get_file_info(S3_Link):
    completed_process_s3 = subprocess.run(['aws','s3','ls', S3_Link],capture_output=True,text=True, check=True)
    result_s3 = completed_process_s3.stdout.split()
    # Should add error catch here for if file isn't present?
    return SingleFastQFile(
        S3_Path = S3_Link,
        file_name = result_s3[3],
        size = int(result_s3[2])
    )

def parse_file_list(S3_list):
    # Take in csv of S3 links to files
    # Output list of SingleFastQFile and PairedFastQFiles
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



# Sets split amount for files
def split_setter(file_list, paired = False):
    split_amount = 100000000000
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
            



def sort_files(project, bucket, group, initial_file_list):
    filepath = bucket + '/' + project + '/' + group + '/'
    completed_process_s3 = subprocess.run(['aws','s3','ls', filepath],capture_output=True,text=True, check=True)
    result_s3 = completed_process_s3.stdout
    lines_s3 = result_s3.strip().split('\n')
    samples = []
    for line in lines_s3:
        parts = line.split()
        if parts:
            if parts[0] == 'PRE':
                samples.append(parts[1])

    file_list = []
    for sample in samples:
        full_filepath = filepath + sample + 'raw/'
        completed_process_samp = subprocess.run(['aws','s3','ls', full_filepath],capture_output=True,text=True, check=True)
        result_samp = completed_process_samp.stdout
        lines_samp = result_samp.strip().split('\n')
        for line in lines_samp:
            parts = line.split()
            if parts:
                if parts[3].endswith('.fastq.gz'):
                    file_list.append({'filename':parts[3],'size':parts[2]})

async def splitter(file):
    while True:
        original_file = await split_queue.get()
        output_files = [original_file.file_name.split('.')[0] + '_split' + str(i+1) + '.fastq.gz') for i in range(0,original_file.split_amount)]
        output_file_command = [('-o ' + o_file) for o_file in output_files]
        output_file_command = ' '.join(output_files)
        await asyncio.run(run('fastqsplitter ' + '-i ' + original_file.file_name + ' ' + output_file_command))
        # This is where you could print out message saying file is getting split
        return file.file_name, output_files # Need to return the original filename and list of split files to then send to uploader


async def downloader(file):
    # Probably want to multithread for this, ask BrianM
    return file


async def uploader(ftp_server_info, file_name):
    # Probably want to multithread for this, ask BrianM
    # file_name is just local name of file upload directly to FTP
    # Return file name so can delete local copy
        


async def Path_Setter(ftp_server_info, Full_File_List):
    split_tasks = []
    upload_tasks = []
    download_tasks = [downloader(file)) for file in Full_File_List]

    # Need to make it so checking when download tasks are complete
    # When download task complete, either add to split_tasks or upload_tasks and then be checking if those are complete
    # When split_task complete add to upload task
    # When upload task complete delete local copy of file
    for completed_download in asyncio.as_completed(download_tasks):
        downloaded_file = await completed_download
        if downloaded_file.split_amount == None:
            upload_tasks.append(uploader(ftp_server_info, downloaded_file))
        else:
            split_tasks.append(splitter(downloaded_file))

    for completed_split in asyncio.as_completed(split_tasks):
        file_to_delete, split_files = await completed_split
        subprocess.run('rm',file_to_delete) # Remove original file to save space? Check if works
        for split_file in split_files:
            upload_tasks.append(uploader(ftp_server_info, split_file))

    for completed_upload in asyncio.as_completed(upload_tasks):
        file_to_delete = await completed_upload
        subprocess.run('rm',file_to_delete) # Remove original file to save space? Check if works
        # Probably print out confirmation file was uploaded?


def main(ftp_server_info, csv_file):
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
    #original_file = Full_File_List[2]
    #output_files = [('-o ' + original_file.file_name.split('.')[0] + '_split' + str(i+1) + '.fastq.gz') for i in range(0,original_file.split_amount)]
    #output_files = ' '.join(output_files)
    #print(output_files)
    #Path_Setter(ftp_server_info, Full_File_List)
    

if __name__ == '__main__':
    args = getArgs()
    # Probably throw this in a function or something
    ftp_info = FTPUploadInfo(
        address = args.address,
        username = args.username,
        password = args.password,
        folder = args.folder)
    main(ftp_info,args.input_file_list)