import subprocess

#subprocess.run(['ls','-l'])

#subprocess.run(['aws','s3','cp', filepath])

# S3 info
bucket='czi-psomagen'
project='marson-mapping-grns-perturb-seq'
group='AN00024637'

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

sorted_file_list = sorted(file_list, key=lambda x: x['filename'])
non_paired_files_to_split = []
files_to_upload = []
paired_files_to_split = []
for file_index_1 in range(len(sorted_file_list)):
    for file_index_2 in range(file_index_1 + 1,len(sorted_file_list)):
        file_1 = sorted_file_list[file_index_1]['filename']
        file_2 = sorted_file_list[file_index_2]['filename']
        if len(file_1) == len(file_2):
            diff_count = 0
            for i in range(len(file_1)):
                if file_1[i] != file_2[i]:
                    char_before = file_1[i-1]
                    # Maybe check here that character before  (i-1) is 'R'?
                    diff_count += 1
                    if diff_count > 1 and file_index_2 == (len(sorted_file_list) - 1):
                        if int(sorted_file_list[file_index_1]['size']) >= 100000000000:
                            if sorted_file_list[file_index_1] not in paired_files_to_split:
                                non_paired_files_to_split.append(sorted_file_list[file_index_1])
                            break
                        else:
                            files_to_upload.append(sorted_file_list[file_index_1])
                            break
                    elif diff_count > 1:
                        break
            if diff_count == 1 and char_before == 'R':
                if int(sorted_file_list[file_index_1]['size']) >= 100000000000 or \
                    int(sorted_file_list[file_index_2]['size']) >= 100000000000:
                        paired_files_to_split.append([sorted_file_list[file_index_1],sorted_file_list[file_index_2]])
                else:
                    files_to_upload.append(sorted_file_list[file_index_1])
                break

print("NON-PAIRED FILES TO SPLIT: \n")
for file in non_paired_files_to_split:
    print(file)
    print('\n')

print("PAIRED FILES TO SPLIT: \n")
for file in paired_files_to_split:
    print(file)
    print('\n')

print("FILES TO UPLOAD DIRECT: \n")
for file in files_to_upload:
    print(file)
    print('\n')
    
                
            
#    if int(file['size']) >= 107374182400:
#        final_files.append(file)
#print(final_files)
    


# ncftp info
#username='testuser1'
#password='Testpasscode7'
#server='ftp://13.56.229.27'
#dest_path='/files'
#file='output.txt'



#subprocess.run(['ncftpput','-u',username,'-p',password,server,dest_path,file])

#completed_process = subprocess.run(['ls','-l'],capture_output=True,text=True, check=True)
#result = completed_process.stdout
#lines = result.strip().split('\n')
#data = []
#for line in lines[1:]:
#    parts = line.split()
#    if parts:
#        data.append({'size':parts[4],'filename':parts[8]})
#print(data)