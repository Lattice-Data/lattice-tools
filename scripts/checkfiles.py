import argparse
import boto3
import botocore
import crcmod
import hashlib
import h5py
import json
import lattice
import logging
import os
import re
import requests
import scanpy as sc
import shutil
import socket
import subprocess
import sys
from datetime import datetime
from ftplib import error_perm, FTP
from urllib.parse import urljoin


EPILOG = '''
Run sanity checks on files.

Examples:

    python %(prog)s --mode production --accessions accessions.txt
    python %(prog)s --mode production --accesssions LATDF101HHH,LATDF102HHH,LATDF100HHH
    python %(prog)s --mode production --query "report/?type=RawSequenceFile&derived_from=/sequencing-runs/2a12eb7b-ed78-466a-9552-7512bdd7f45f/"
    python %(prog)s --s3-file s3://submissions-czi012eye/chen_2020/19D014_NeuNT_2_outs/raw_feature_bc_matrix.h5 --file-format hdf5

This relies on local variables to be defined based on the --mode you provide
to direct the updates to a server and to provide permissions
For example, if specifying --mode production, to make the changes on a local instance,
the following variables need to be defined...
PRODUCTION_KEY, PRODUCTION_SECRET, PRODUCTION_SERVER

For more details:

        python %(prog)s --help
'''


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--mode', '-m',
                        help='The machine to run on.')
    parser.add_argument('--update',
                        default=False,
                        action='store_true',
                        help='Let the script proceed with the changes.  Default is False'),
    parser.add_argument('--query', '-q',
                        help="override the file search query, e.g. 'accession=ENCFF000ABC'")
    parser.add_argument('--accessions', '-a',
                        help='one or more file accessions to check, comma separated or a file containing a list of file accessions to check')
    parser.add_argument('--include-validated',
                        default=False,
                        action='store_true',
                        help='Check all files even if they are validated=True in the Lattice database')
    parser.add_argument('--s3-file',
                        help="path to a file at s3 to check, comma separated or a file containing a list of file accessions to check")
    parser.add_argument('--ext-file',
                        help="path to a file elsewhere to check")
    parser.add_argument('--file-format',
                        help='the specified file format if an s3-file or local-file is being checked')
    args = parser.parse_args()
    return args


checkfiles_version = '2.0'

abstract_file_types = [
    'Item',
    'File',
    'DataFile',
    'AnalysisFile'
]

GZIP_TYPES = [
    "fastq",
    "mex"
]

machine_pattern = '^(@[a-zA-Z\d]+[a-zA-Z\d_-]*'
run_id_pattern = '[a-zA-Z\d-]+'
flowcell_pattern = '[a-zA-Z\d_-]+'

output_dir = 'output'
chkfls_dir = 'checkfiles'

read_name_prefix = re.compile(
    machine_pattern + ':' + run_id_pattern + ':' + flowcell_pattern +
    ':\d+:\d+:\d+:\d+)$'
)

illumina_read_name_pattern = re.compile(
    machine_pattern + ':' + run_id_pattern + ':' + flowcell_pattern +
    ':\d+:\d+:\d+:\d+)'
)

illumina_read_name_pattern_no_flowcell = re.compile(
    machine_pattern +
    ':\d+:\d+:\d+:\d+)'
)

srr_read_name_pattern = re.compile(
    '^(@[S|E]RR[\d.]+)$'
)


def is_path_gzipped(path):
    try:
        f = open(path, 'rb')
    except IsADirectoryError:
        return 'IsADirectoryError'
    else:
        magic_number = f.read(2)
        return magic_number == b'\x1f\x8b'


def process_illumina_read_name_pattern(read_name, read_numbers_set, signatures_set, signatures_no_smp_in_set, srr_flag):
    read_name_array = re.split(r'[:\s]', read_name)
    machine = read_name_array[0].replace('@','')
    flowcell = read_name_array[2]
    lane_number = read_name_array[3]
    if srr_flag:
        member_pair = list(read_numbers_set)[0]
    else:
        member_pair = read_name_array[-4]
        read_numbers_set.add(member_pair)
    smp_index = read_name_array[-1]
    signatures_set.add(
        machine + ':' + flowcell + ':' +
        lane_number + ':' + member_pair + ':' + smp_index)
    signatures_no_smp_in_set.add(
        machine + ':' + flowcell + ':' +
        lane_number + ':' + member_pair)


def process_no_flowcell_read_name_pattern(read_name, read_numbers_set, signatures_set, signatures_no_smp_in_set, srr_flag):
    read_name_array = re.split(r'[:\s]', read_name)
    machine = read_name_array[0].replace('@','')
    lane_number = read_name_array[1]
    if srr_flag:
        member_pair = list(read_numbers_set)[0]
    else:
        member_pair = read_name_array[-4]
        read_numbers_set.add(member_pair)
    smp_index = read_name_array[-1]
    signatures_set.add(
        machine + ':' + lane_number + ':' + member_pair + ':' + smp_index)
    signatures_no_smp_in_set.add(
        machine + ':' + lane_number + ':' + member_pair)


def process_illumina_prefix(read_name, signatures_set, old_illumina_current_prefix, read_numbers_set, srr_flag):
    if srr_flag:
        member_pair = list(read_numbers_set)[0]
    else:
        read_number = '1'
        read_numbers_set.add(member_pair)
    read_name_array = re.split(r':', read_name)

    if len(read_name_array) > 3:
        flowcell = read_name_array[2]
        lane_number = read_name_array[3]

        prefix = flowcell + ':' + lane_number
        if prefix != old_illumina_current_prefix:
            old_illumina_current_prefix = prefix

            signatures_set.add(
                flowcell + ':' + lane_number + ':' +
                member_pair + '::' + read_name)

    return old_illumina_current_prefix


def process_read_name_line(read_name_line, old_illumina_current_prefix, read_numbers_set, signatures_no_smp_in_set, signatures_set, read_lengths_dictionary, errors, srr_flag):
    read_name = read_name_line.strip()
    words_array = re.split(r'\s', read_name)
    if illumina_read_name_pattern.match(read_name) is not None:
        process_illumina_read_name_pattern(
            read_name,
            read_numbers_set,
            signatures_set,
            signatures_no_smp_in_set,
            srr_flag)
    elif illumina_read_name_pattern_no_flowcell.match(read_name) is not None:
        process_no_flowcell_read_name_pattern(
            read_name,
            read_numbers_set,
            signatures_set,
            signatures_no_smp_in_set,
            srr_flag)
    elif srr_read_name_pattern.match(read_name.split(' ')[0]) is not None:
        # in case the readname is following SRR format, read number will be
        # defined using SRR format specifications, and not by the illumina portion of the read name
        # srr_flag is used to distinguish between srr and "regular" readname formats
        srr_portion = read_name.split(' ')[0]
        if srr_portion.count('.') == 2:
            read_numbers_set.add(srr_portion[-1])
        else:
            read_numbers_set.add('1')
        illumina_portion = read_name.split(' ')[1]
        old_illumina_current_prefix = process_read_name_line('@'+illumina_portion,
                                                            old_illumina_current_prefix,
                                                            read_numbers_set,
                                                            signatures_no_smp_in_set,
                                                            signatures_set,
                                                            read_lengths_dictionary,
                                                            errors, True)
    else:
        # unrecognized read_name_format
        # current convention is to include WHOLE
        # readname at the end of the signature
        if len(words_array) == 1:
            if read_name_prefix.match(read_name) is not None:
                # new illumina without second part
                old_illumina_current_prefix = process_illumina_prefix(
                    read_name,
                    signatures_set,
                    old_illumina_current_prefix,
                    read_numbers_set,
                    srr_flag)

    return old_illumina_current_prefix


def process_fastq_file(job):
    logging.info('Getting fastq metadata')
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri', item.get('external_uri'))
    local_path = download_url.split('/')[-1]

    read_numbers_set = set()
    signatures_set = set()
    signatures_no_smp_in_set = set()
    read_lengths_dictionary = {}
    read_count = 0
    old_illumina_current_prefix = 'empty'
    try:
        fastq_data_stream = subprocess.Popen(['gunzip --stdout {}'.format(
                                                        local_path)],
                                                        shell=True,
                                                        executable='/bin/bash',
                                                        stdout=subprocess.PIPE)
        line_index = 0
        for encoded_line in fastq_data_stream.stdout:
            try:
                line = encoded_line.decode('utf-8')
            except UnicodeDecodeError:
                errors['readname_encoding'] = 'Error occured, while decoding the readname string.'
            else:
                line_index += 1
                if line_index == 1:
                    
                    # may be from here deliver a flag about the presence/absence of the readnamedetails
                    old_illumina_current_prefix = \
                        process_read_name_line(
                            line,
                            old_illumina_current_prefix,
                            read_numbers_set,
                            signatures_no_smp_in_set,
                            signatures_set,
                            read_lengths_dictionary,
                            errors, False)
                elif line_index == 2:
                    read_count += 1
                    length = len(line.strip())
                    if length not in read_lengths_dictionary:
                        read_lengths_dictionary[length] = 0
                    read_lengths_dictionary[length] += 1

            line_index = line_index % 4
    except IOError:
        errors['unzipped_fastq_streaming'] = 'Error occured, while streaming unzipped fastq.'
    else:
        if line_index != 0:
            errors['incomplete fastq file'] = 'Last line was not a factor of 4, but has a remainder of {}'.format(line_index)
        # read_count update
        logging.info('Determining read count')
        results['read_count'] = read_count

        # read1/read2
        logging.info('Determining read number')
        if len(read_numbers_set) == 1:
            read_number = next(iter(read_numbers_set))
            results['member_pair'] = read_number

        # read_length
        logging.info('Determining read lengths')
        results['read_length'] = max(read_lengths_dictionary, key=read_lengths_dictionary.get)

        read_lengths_list = []
        for k in sorted(read_lengths_dictionary.keys()):
            read_lengths_list.append({'read_length': k, 'read_count': read_lengths_dictionary[k]})

        results['read_lengths_list'] = ', '.join(map(str, read_lengths_list))

        # signatures
        logging.info('Determining read signatures')
        all_flowcells = set()
        all_smp_ins = set()
        for sign in signatures_set:
            sign_array = sign.split(':')
            if len(sign_array) > 4:
                flowcell = (sign_array[0], sign_array[1], sign_array[2])
                all_flowcells.add(flowcell)
                all_smp_ins.add(sign_array[4])

                #results['sample indices'] = all_smp_ins

                flowcell_details = []
                for flowcell in all_flowcells:
                    flowcell_details.append({'machine': flowcell[0], 'flowcell': flowcell[1], 'lane': flowcell[2]})
                results['flowcell_details'] = flowcell_details
            else:
                flowcell = (sign_array[0], sign_array[1])
                all_flowcells.add(flowcell)
                all_smp_ins.add(sign_array[3])

                #results['sample indices'] = all_smp_ins

                flowcell_details = []
                for flowcell in all_flowcells:
                    flowcell_details.append({'machine': flowcell[0], 'lane': flowcell[1]})
                results['flowcell_details'] = flowcell_details
        results['signature'] = list(signatures_set)


def process_h5matrix_file(job):
    logging.info('Getting h5 matrix metadata')
    feature_type_mapping = {
        'Gene Expression': 'gene',
        'Peaks': 'peak',
        'Antibody Capture': 'antibody capture',
    }

    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri', item.get('external_uri'))
    local_path = download_url.split('/')[-1]

    hdf5_validate = h5py.is_hdf5(local_path)
    results['is_hdf5'] = hdf5_validate
    if hdf5_validate != True:
        errors['is_hdf5'] = hdf5_validate

    if local_path.endswith('.h5'):
        adata = sc.read_10x_h5(local_path, gex_only=False)
    else:
        adata = sc.read_h5ad(local_path, backed='r')

    results['observation_count'] = adata.obs.shape[0]

    if 'genome' in adata.var.columns:
        results['genomes'] = adata.var['genome'].unique().tolist()

    if 'feature_types' not in adata.var.columns:
        errors['var.feature_types'] = 'is missing'
    else:
        feature_counts = []
        for k,v in adata.var['feature_types'].value_counts().to_dict().items():
            key = feature_type_mapping.get(k,k)
            feature_counts.append({'feature_type': key, 'feature_count': v})
            if k not in feature_type_mapping:
                errors[f'var.feature_types[{k}]'] = 'not a valid feature type'
        results['feature_counts'] = feature_counts

    if 'gene_ids' not in adata.var.columns:
        errors['var.gene_ids'] = 'is missing'
    else:
        if 'feature_types' in adata.var.columns:
            not_ensg = [g for g in adata.var[adata.var['feature_types'] == 'Gene Expression']['gene_ids'] if g.startswith('ENS') == False]
            if len(not_ensg) > 0:
                errors['var.gene_ids'] = str(len(not_ensg)) + ' features in var.gene_ids not ENS-formatted'
            with_version = [g for g in adata.var[adata.var['feature_types'].isin(['Peaks', 'Antibody Capture']) == False]['gene_ids'] if '.' in g]
        else:
            with_version = [g for g in adata.var['gene_ids'] if '.' in g]
        if len(with_version) > 0:
            errors['ENSG format'] = str(len(with_version)) + ' versions in var.gene_ids'

        pary_genes = [g for g in adata.var['gene_ids'] if 'PAR_Y' in g]
        for g in pary_genes:
            symb = adata.var.loc[adata.var['gene_ids'] == g].index[0]
            if adata.var[adata.var.index == symb].shape[0] < 2:
                errors['PAR_Y symbols'] = 'expecting PAR_Y symbols to be duplicated'
            if g.endswith('PAR_Y') != True:
                errors['PAR_Y ID'] = 'expecting PAR_Y genes to end with PAR_Y'
            if 'gene_versions' in adata.var.columns:
                ver = adata.var.loc[adata.var['gene_ids'] == g]['gene_versions'][0]
                if g != '_'.join([ver.split('.')[0]] + ver.split('_')[1:]):
                    errors['PAR_Y version'] = 'PAR_Y gene version does not match ID'

    if 'gene_versions' in adata.var.columns:
        without_version = [g for g in adata.var['gene_versions'] if '.' not in g]
        if len(with_version) > 0:
            errors['ENSG.N format'] = str(len(without_version)) + ' IDs without version in var.gene_versions'

    if len([g for g in adata.var.index if g.startswith('ENS')]) == adata.var.shape[0]:
            results['feature_keys'] = ['Ensembl gene ID']

    check_dates = ['Mar-1','1-Mar','Dec-1','1-Dec']
    date_symbols = [s for s in check_dates if s in adata.var.index]
    if len(date_symbols) > 0:
        errors['date-formatted symbols'] = 'symbols present that have been reformatted: ' + ','.join(date_symbols)

    incorrect_cell_ids = [c for c in adata.obs.index if c.endswith('.1')]
    if len(incorrect_cell_ids) > 0:
        errors['cell_id suffix'] = 'obs.index values should not end in .1'


def download_s3_directory(job):
    item = job['item']
    errors = job['errors']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    dir_path = download_url.replace('s3://{}/'.format(bucket_name), '')

    s3client = boto3.client("s3")
    tmp_dir = 'raw_feature_bc_matrix'
    os.mkdir(tmp_dir)
    logging.info(file_name + ' downloading')
    for file_name in ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']:
        try:
            s3client.download_file(bucket_name, dir_path + '/' + file_name, '{}/{}'.format(tmp_dir, file_name))
        except botocore.exceptions.ClientError as e:
            errors['s3 uri error'] = e.response['Error']['Message']
        else:
            logging.info(file_name + ' downloaded')

    return tmp_dir, job


def download_s3_file(job):
    item = job['item']
    errors = job['errors']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
    file_name = download_url.split('/')[-1]

    job['download_start'] = datetime.now()

    s3client = boto3.client("s3")
    logging.info(file_name + ' downloading')
    try:
        s3client.download_file(bucket_name, file_path, file_name)
    except botocore.exceptions.ClientError as e:
        errors['s3 uri error'] = e.response['Error']['Message']
    else:
        logging.info(file_name + ' downloaded')
        job['download_stop'] = datetime.now()
        job['download_time'] = job['download_stop'] - job['download_start']

    return file_name, job


def download_external(job):
    item = job['item']
    errors = job['errors']

    download_url = item.get('external_uri')
    ftp_server = download_url.split('/')[2]
    ftp = FTP(ftp_server)
    ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 75)
    #in some cases the next line needs to be commented out
    ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
    ftp.login(user='anonymous', passwd = 'password')

    file_path = download_url.replace('ftp://{}/'.format(ftp_server), '')
    file_name = download_url.split('/')[-1]

    job['download_start'] = datetime.now()

    logging.info(file_name + ' downloading')
    try:
        ftp.retrbinary('RETR ' + file_path, open(file_name, 'wb').write)
    except error_perm as e:
        errors['file download error'] = e
        os.remove(file_name)
    else:
        ftp.quit()
        logging.info(file_name + ' downloaded')
        job['download_stop'] = datetime.now()
        job['download_time'] = job['download_stop'] - job['download_start']

    return file_name, job


def set_s3_tags(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    file_path = download_url.replace('s3://{}/'.format(bucket_name), '')

    tagging = {'TagSet': [
        {
            'Key': 'validated',
            'Value': checkfiles_version
        }
    ]}

    s3client = boto3.client("s3")
    response = s3client.put_object_tagging(Bucket=bucket_name, Key=file_path, Tagging=tagging)
    if response['ResponseMetadata']['HTTPStatusCode'] == 200:
        job['patch_s3tag'] = 'success'
    else:
        job['patch_s3tag'] = 'status code:{}'.format(response['HTTPStatusCode'])


def compare_with_db(job, connection):
    server = connection.server

    file = job['item']
    errors = job['errors']
    results = job['results']
    post_json = {}

    logging.info('Comparing results with metadata in Lattice DB')

    metadata_consistency = []
    metadata_inconsistency = []

    obj_types = file['@type']
    schema = next(iter(set(obj_types) - set(abstract_file_types)))
    schema_url = urljoin(server, 'profiles/' + schema + '/?format=json')
    schema_properties = requests.get(schema_url).json()['properties']

    for key in results.keys():
        file_value = file.get(key)
        results_value = results.get(key)
        # if it's a schema property currently absent, prepare to patch it
        if not file_value and schema_properties.get(key):
            post_json[key] = results_value
        # if the file information matches the current database metadata, log it
        elif results_value == file_value:
            outcome = '{} consistent ({})'.format(key, results_value)
            metadata_consistency.append(outcome)
        elif file_value != None:
            #first check for embedded properties, specifically for flowcell_details
            if isinstance(results_value, list) and isinstance(file_value, list) \
                and len(results_value) == 1 and len(file_value) == 1:
                if schema_properties.get(key) and schema_properties[key].get('items'):
                    post_flag = True
                    results_obj = results[key][0]
                    file_obj = file[key][0]
                    for subkey in results_obj:
                        if schema_properties[key]['items']['properties'].get(subkey):
                            if results_obj.get(subkey) == file_obj.get(subkey):
                                outcome = '{}.{} consistent ({})'.format(key, subkey, results_obj.get(subkey))
                                metadata_consistency.append(outcome)
                            elif file_obj.get(subkey) and results_obj.get(subkey) != file_obj.get(subkey):
                                post_flag = False
                                outcome = '{}.{} inconsistent ({}-s3file, {}-submitted)'.format(key, subkey, results_obj.get(subkey), file_obj.get(subkey))
                                metadata_inconsistency.append(outcome)
                        else:
                            del results_obj[subkey]
                    if post_flag == True:
                        post_json[key] = [results_obj]
            elif isinstance(results_value, list) and isinstance(file_value, list):
                if len(results_value) == len(file_value) and len([x for x in results_value if x in file_value]) == len(results_value):
                    outcome = '{} consistent ({})'.format(key, results_value)
                    metadata_consistency.append(outcome)
                else:
                    outcome = '{} inconsistent ({}-s3file, {}-submitted)'.format(key, results_value, file_value)
                    metadata_inconsistency.append(outcome)
            # we have an inconsistency to log
            else:
                outcome = '{} inconsistent ({}-s3file, {}-submitted)'.format(key, results_value, file_value)
                metadata_inconsistency.append(outcome)
    if len(metadata_inconsistency) == 0 and schema_properties.get('validated') and file.get('validated') != True:
        post_json['validated'] = True

    job['post_json'] = post_json
    if metadata_consistency:
        results['metadata_consistency'] = metadata_consistency
    if metadata_inconsistency:
        errors['metadata_inconsistency'] = metadata_inconsistency


def check_file(job):
    file = job['item']
    errors = job['errors']
    results = job['results']

    download_url = file.get('s3_uri', file.get('external_uri'))
    local_path = download_url.split('/')[-1]

    job['check_start'] = datetime.now()

    # check file size
    logging.info('Getting file size')
    file_stat = os.stat(local_path)
    results['file_size'] = file_stat.st_size
    if local_path.endswith('.h5ad'):
        adata = sc.read_h5ad(local_path, backed='r')
        adata.write(filename='temp.h5ad', compression='gzip')
        compressed_file_stat = os.stat('temp.h5ad')
        if (file_stat.st_size / 1.5) > compressed_file_stat.st_size:
            errors['compression'] = f'file size is {file_stat.st_size} but compressed file size is {compressed_file_stat.st_size}'
        os.remove('temp.h5ad')

    # get the md5sum
    # faster than doing it in Python
    try:
        logging.info('Getting file md5')
        output = subprocess.check_output('md5sum {}'.format(local_path),
            shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        errors['md5sum'] = e.output.decode(errors='replace').rstrip('\n')
    else:
        results['md5sum'] = output[:32].decode(errors='replace')
        try:
            int(results['md5sum'], 16)
        except ValueError:
            errors['md5sum'] = output.decode(errors='replace').rstrip('\n')


    # get the sha256
    sha256_hash = hashlib.sha256()
    with open(local_path,"rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            sha256_hash.update(byte_block)
        results['sha256'] = sha256_hash.hexdigest()


    # get the crc32c
    with open(local_path, 'rb') as f:
        crc32c_func = crcmod.predefined.Crc('crc-32c')
        while True:
            s = f.read(65536)
            if not s:
                break
            crc32c_func.update(s)
        results['crc32c'] = crc32c_func.hexdigest().lower()


    # check for correct gzip status
    try:
        is_gzipped = is_path_gzipped(local_path)
    except Exception as e:
        return job
    else:
        if file.get('file_format') not in GZIP_TYPES:
            if is_gzipped:
                errors['gzip'] = 'Expected un-gzipped file'
        elif not is_gzipped:
            errors['gzip'] = 'Expected gzipped file'
        else:
            logging.info('Getting file content md5')
            # May want to replace this with something like:
            # $ cat $local_path | tee >(md5sum >&2) | gunzip | md5sum
            # or http://stackoverflow.com/a/15343686/199100    
            # if gzipped, unzip, get md5sum
            job['content_md5sum_start'] = datetime.now()
            try:
                output = subprocess.check_output(
                    'set -o pipefail; gunzip --stdout {} | md5sum'.format(local_path),
                    shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                errors['content_md5sum'] = e.output.decode(errors='replace').rstrip('\n')
            else:
                results['content_md5sum'] = output[:32].decode(errors='replace')
                job['content_md5sum_stop'] = datetime.now()
                job['content_md5sum_time'] = job['content_md5sum_stop'] - job['content_md5sum_start']

            # do format-specific validation
            if file.get('file_format') == 'fastq':
                process_fastq_file(job)
                os.remove(local_path)
                logging.info(local_path + ' removed')
    if file.get('file_format') == 'hdf5':
        process_h5matrix_file(job)
        os.remove(local_path)
        logging.info(local_path + ' removed')

    job['check_stop'] = datetime.now()
    job['check_time'] = job['check_stop'] - job['check_start']

    return job


def fetch_files(report_out, connection=None, query=None, accessions=None, s3_file=None, ext_file=None, file_format=None, include_validate=False):
    logging.info('Fetching files that need validation')
    if accessions or query:
        server = connection.server
        if accessions:
            # checkfiles using a file with a list of file accessions to be checked
            if '.' in accessions:
                r = None
                if os.path.isfile(accessions):
                    ACCESSIONS = [line.rstrip('\n') for line in open(accessions)]
            # checkfiles using a list of accessions in the command, comma separated
            else:
                ACCESSIONS = accessions.split(',')

        # checkfiles using a query
        else:
            query_url = urljoin(server, query.replace('report', 'search') + '&format=json&limit=all&field=accession')
            r = requests.get(query_url, auth=connection.auth)
            try:
                r.raise_for_status()
            except requests.HTTPError:
                return
            else:
                ACCESSIONS = [x['accession'] for x in r.json()['@graph']]

        jobs = []
        for acc in ACCESSIONS:
            item_url = urljoin(server, acc + '/?frame=object')
            fileObject = requests.get(item_url, auth=connection.auth)
            if not fileObject.ok:
                errors['file_HTTPError'] = ('HTTP error: unable to get file object')
            else:
                check_me_flag = True
                blockers = []
                file_json = fileObject.json()
                if file_json.get('file_format') == 'fastq' and len(file_json.get('derived_from')) == 1:
                    seq_run = file_json.get('derived_from')[0]
                    item_url = urljoin(server, seq_run + '/?frame=object')
                    seqrunObject = requests.get(item_url, auth=connection.auth).json()
                    file_json['derived_from'] = seqrunObject
                if file_json.get('no_file_available') == True:
                    blockers.append('marked as no_file_available')
                    check_me_flag = False
                if not file_json.get('s3_uri') and not file_json.get('external_uri'):
                    blockers.append('uri not submitted')
                    check_me_flag = False
                if file_json.get('validated') == True and not include_validate:
                    blockers.append('already validated')
                    check_me_flag = False
                if check_me_flag == False:
                    uri = file_json.get('s3_uri',file_json.get('external_uri'))
                    if uri == None:
                        uri = ''
                    out = open(output_dir + '/' + chkfls_dir + '/' + report_out, 'a')
                    out.write(acc + '\t' + uri + '\t' + ','.join(blockers) + '\n')
                    out.close()
                elif check_me_flag == True:
                    job = {
                        'item': file_json,
                        'results': {},
                        'errors': {}
                    }
                    jobs.append(job)
        return jobs

    # checkfiles on a file that is not in the Lattice database but is at s3
    elif s3_file:
        if os.path.isfile(s3_file):
            s3FILES = [line.rstrip('\n') for line in open(s3_file)]
        else:
            s3FILES = s3_file.split(',')
        jobs = []
        for file in s3FILES:
            job = {
                'item': {
                    'accession': 'not yet submitted',
                    's3_uri': file
                },
                'results': {},
                'errors': {}
            }
            if file_format:
                job['item']['file_format'] = file_format
            jobs.append(job)
        return jobs

    # checkfiles on a file that is not in the Lattice database but is available elsewhere
    elif ext_file:
        jobs = [{
            'item': {
                'accession': 'not yet submitted',
                'external_uri': ext_file
            },
            'results': {},
            'errors': {}
        }]
        if file_format:
            jobs[0]['item']['file_format'] = file_format
        return jobs


def report(job):
    file_obj = job.get('item')
    tab_report = '\t'.join([
        file_obj.get('accession', file_obj.get('uuid')),
        str(file_obj.get('s3_uri', file_obj.get('external_uri'))),
        str(job['errors']),
        str(job['results']),
        str(job.get('post_json')),
        job.get('patch_result', 'n/a'),
        job.get('patch_s3tag', 'n/a'),
        str(job.get('download_time')),
        str(job.get('check_time')),
        str(job.get('content_md5sum_time'))
    ])
    return tab_report + '\n'


def main():
    logging.basicConfig(filename='checkfiles.log', level=logging.INFO)
    logging.info('Started')

    args = getArgs()
    if (args.query or args.accessions) and not args.mode:
        sys.exit('ERROR: --mode is required with --query/--accessions')

    arg_count = 0
    for arg in [args.query, args.accessions, args.s3_file, args.ext_file]:
        if arg:
            arg_count += 1
    if arg_count != 1:
        sys.exit('ERROR: exactly one of --query, --accessions, --s3-file, --ext-file is required, {} given'.format(arg_count))


    if args.mode:
        connection = lattice.Connection(args.mode)
    else:
        connection = ''

    initiating_run = 'STARTING Checkfiles version {}'.format(checkfiles_version)
    logging.info(initiating_run)

    # Checking for presence of / creating output folder and associated subfolder
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
    if os.path.exists(output_dir + '/' + chkfls_dir) == False:
        os.mkdir(output_dir + '/' + chkfls_dir)

    timestr = datetime.now().strftime('%Y_%m_%d-%H_%M_%S')
    report_out = 'report_{}.tsv'.format(timestr)
    logging.info('Writing results to {}'.format(report_out))
    report_headers = '\t'.join([
        'identifier',
        'uri',
        'errors',
        'results',
        'json_patch',
        'Lattice patched?',
        'S3 tag patched?',
        'download_time',
        'check_time',
        'content_md5sum_time'
    ])
    with open(output_dir + '/' + chkfls_dir + '/' + report_out, 'w') as out:
        out.write(report_headers + '\n')

    jobs = fetch_files(report_out, connection, args.query, args.accessions, args.s3_file, args.ext_file, args.file_format, args.include_validated)

    if jobs:
        logging.info('CHECKING {} files'.format(len(jobs)))
        for job in jobs:
            file_obj = job.get('item')
            logging.info('Starting {}'.format(file_obj.get('@id', 'File not in DB')))
            if file_obj.get('external_uri'):
                local_file, job = download_external(job)
            elif file_obj.get('file_format') == 'mex':
                local_file, job = download_s3_directory(job)
            else:
                local_file, job = download_s3_file(job)
            if os.path.exists(local_file):
                check_file(job)
                if not args.s3_file and not args.ext_file:
                    compare_with_db(job, connection)
                    if job['results'].get('flowcell_details') and file_obj.get('derived_from'):
                        dets = job['results']['flowcell_details']
                        sorted_dets = sorted(dets, key=lambda k: (k.get('machine'), k.get('flowcell'), k.get('lane')))
                    if job['post_json'] and not job['errors'] and args.update:
                        logging.info('PATCHING {}'.format(file_obj.get('accession')))
                        patch = lattice.patch_object(file_obj.get('accession'), connection, job['post_json'])
                        job['patch_result'] = patch['status']
                        if file_obj.get('s3_uri'):
                            set_s3_tags(job)
            out = open(output_dir + '/' + chkfls_dir + '/' + report_out, 'a')
            out.write(report(job))
            out.close()

        finishing_run = 'FINISHED Checkfiles at {}'.format(datetime.now())
        logging.info(finishing_run)
    else:
        logging.info('FINISHED No files to check, see report*.tsv for details')

    logging.info('Results written to {}'.format(report_out))
    logging.info('Finished')

if __name__ == '__main__':
    main()
