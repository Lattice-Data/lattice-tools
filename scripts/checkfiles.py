import argparse
import boto3
from datetime import datetime
import h5py #conda install h5py
import json
import lattice
import os
import re
import requests
import shutil
import subprocess
import sys
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
                        help="one or more file accessions to check, comma separated or a file containing a list of file accessions to check")
    parser.add_argument('--s3-file',
                        help="path to a file at s3 to check")
    parser.add_argument('--file-format',
                        help="the specified file format if an s3-file or local-file is being checked")
    args = parser.parse_args()
    return args


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

read_name_prefix = re.compile(
    machine_pattern + ':' + run_id_pattern + ':' + flowcell_pattern +
    ':\d+:\d+:\d+:\d+)$')

illumina_read_name_pattern = re.compile(
    machine_pattern + ':' + run_id_pattern + ':' + flowcell_pattern +
    ':\d+:\d+:\d+:\d+[\s_][123]:[YXN]:[0-9]+:([ACNTG\+]*|[0-9]*))$'
)

illumina_read_name_pattern_no_flowcell = re.compile(
    machine_pattern +
    ':\d+:\d+:\d+:\d+[\s_][123]:[YXN]:[0-9]+:([ACNTG\+]*|[0-9]*))$'
)

srr_read_name_pattern = re.compile(
    '^(@SRR[\d.]+)$'
)


def is_path_gzipped(path):
    try:
        f = open(path, 'rb')
    except IsADirectoryError:
        return 'IsADirectoryError'
    else:
        magic_number = f.read(2)
        return magic_number == b'\x1f\x8b'


def process_illumina_read_name_pattern(read_name, read_numbers_set, signatures_set, signatures_no_barcode_set, srr_flag):
    read_name_array = re.split(r'[:\s]', read_name)
    machine = read_name_array[0].replace('@','')
    flowcell = read_name_array[2]
    lane_number = read_name_array[3]
    if srr_flag:
        member_pair = list(read_numbers_set)[0]
    else:
        member_pair = read_name_array[-4]
        read_numbers_set.add(member_pair)
    barcode_index = read_name_array[-1]
    signatures_set.add(
        machine + ':' + flowcell + ':' +
        lane_number + ':' + member_pair + ':' + barcode_index)
    signatures_no_barcode_set.add(
        machine + ':' + flowcell + ':' +
        lane_number + ':' + member_pair)


def process_no_flowcell_read_name_pattern(read_name, read_numbers_set, signatures_set, signatures_no_barcode_set, srr_flag):
    read_name_array = re.split(r'[:\s]', read_name)
    machine = read_name_array[0].replace('@','')
    lane_number = read_name_array[1]
    if srr_flag:
        member_pair = list(read_numbers_set)[0]
    else:
        member_pair = read_name_array[-4]
        read_numbers_set.add(member_pair)
    barcode_index = read_name_array[-1]
    signatures_set.add(
        machine + ':' + lane_number + ':' + member_pair + ':' + barcode_index)
    signatures_no_barcode_set.add(
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


def process_read_name_line(read_name_line, old_illumina_current_prefix, read_numbers_set, signatures_no_barcode_set, signatures_set, read_lengths_dictionary, errors, srr_flag):
    read_name = read_name_line.strip()
    words_array = re.split(r'\s', read_name)
    if illumina_read_name_pattern.match(read_name) is not None:
        process_illumina_read_name_pattern(
            read_name,
            read_numbers_set,
            signatures_set,
            signatures_no_barcode_set,
            srr_flag)
    elif illumina_read_name_pattern_no_flowcell.match(read_name) is not None:
        process_no_flowcell_read_name_pattern(
            read_name,
            read_numbers_set,
            signatures_set,
            signatures_no_barcode_set,
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
                                                            signatures_no_barcode_set,
                                                            signatures_set,
                                                            read_lengths_dictionary,
                                                            errors, True, read_name_details)
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
            else:
                errors['fastq_format_readname'] = read_name
                # the only case to skip update content error - due to the changing
                # nature of read names
        else:
            errors['fastq_format_readname'] = read_name
    # found a match to the regex of "almost" illumina read_name

    return old_illumina_current_prefix


def process_fastq_file(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri')
    local_path = download_url.split('/')[-1]

    read_numbers_set = set()
    signatures_set = set()
    signatures_no_barcode_set = set()
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
                            signatures_no_barcode_set,
                            signatures_set,
                            read_lengths_dictionary,
                            errors, False)
            if line_index == 2:
                read_count += 1
                length = len(line.strip())
                if length not in read_lengths_dictionary:
                    read_lengths_dictionary[length] = 0
                read_lengths_dictionary[length] += 1

            line_index = line_index % 4
    except IOError:
        errors['unzipped_fastq_streaming'] = 'Error occured, while streaming unzipped fastq.'
    else:
        # read_count update
        results['read_count'] = read_count

        # read1/read2
        if len(read_numbers_set) == 1:
            read_number = next(iter(read_numbers_set))
            results['member_pair'] = read_number
        elif len(read_numbers_set) > 1:
            errors['inconsistent_read_numbers'] = \
                'fastq file contains mixed read numbers ' + \
                '{}.'.format(', '.join(sorted(list(read_numbers_set))))

        # read_length
        results['read_length'] = max(read_lengths_dictionary, key=read_lengths_dictionary.get)

        read_lengths_list = []
        for k in sorted(read_lengths_dictionary.keys()):
            read_lengths_list.append({'read_length': k, 'read_count': read_lengths_dictionary[k]})

        results['read_lengths_list'] = ', '.join(map(str, read_lengths_list))

        # signatures
        all_flowcells = set()
        all_barcodes = set()
        for sign in signatures_set:
            sign_array = sign.split(':')
            if len(sign_array) > 4:
                flowcell = (sign_array[0], sign_array[1], sign_array[2])
                all_flowcells.add(flowcell)
                all_barcodes.add(sign_array[4])

                results['barcodes'] = all_barcodes

                flowcell_details = []
                for flowcell in all_flowcells:
                    flowcell_details.append({'machine': flowcell[0], 'flowcell': flowcell[1], 'lane': flowcell[2]})
                results['flowcell_details'] = flowcell_details
            else:
                flowcell = (sign_array[0], sign_array[1])
                all_flowcells.add(flowcell)
                all_barcodes.add(sign_array[3])

                results['barcodes'] = all_barcodes

                flowcell_details = []
                for flowcell in all_flowcells:
                    flowcell_details.append({'machine': flowcell[0], 'lane': flowcell[1]})
                results['flowcell_details'] = flowcell_details


def process_h5matrix_file(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri')
    local_path = download_url.split('/')[-1]

    # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
    # https://docs.h5py.org/en/stable/
    f = h5py.File(local_path, 'r')
    for k in list(f.keys()):
        dset = f[k]
        results['barcode_count'] = dset['barcodes'].shape[0]
        results['features'] = dset['features']['id'].shape[0]


def process_mexmatrix_file(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    tmp_dir = 'raw_feature_bc_matrix'
    mtx_path = 'matrix.mtx.gz'
    feature_path = 'features.tsv.gz'
    barcode_path = 'barcodes.tsv.gz'

    # do the matrix file in the directory
    try:
        matrix_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, mtx_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['matrix_information_extraction'] = 'Failed to extract information from ' + mtx_path

    try:
        line_index = 0
        for encoded_line in matrix_stream.stdout:
            line_index += 1
            if line_index == 3:
                count_summary = encoded_line.decode('utf-8').strip().split()
                features_count_frommatrix = int(count_summary[0])
                barcodes_count_frommatrix = int(count_summary[1])
                break
    except IOError:
        errors['unzipped_matrix_streaming'] = 'Error occured, while streaming unzipped matrix file.'

    # do the features file in the directory
    try:
        features_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, feature_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['features_information_extraction'] = 'Failed to extract information from ' + feature_path

    features_row_count = 0
    try:
        for encoded_line in features_stream.stdout:
            features_row_count += 1
    except IOError:
        errors['unzipped_features_streaming'] = 'Error occured, while streaming unzipped features file.'
    else:
        # read_count update
        if features_row_count == features_count_frommatrix:
            results['feature_count'] = features_row_count
        else:
            errors['feature_count_discrepancy'] = 'Feature count from matrix ({}) does not match row count in features.tsv ({})'.format(
                features_count_frommatrix, features_row_count)

    # do the barcodes file in the directory
    try:
        barcodes_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, barcode_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['barcodes_information_extraction'] = 'Failed to extract information from ' + barcode_path

    barcodes_row_count = 0
    try:
        for encoded_line in barcodes_stream.stdout:
            barcodes_row_count += 1
    except IOError:
        errors['unzipped_barcodes_streaming'] = 'Error occured, while streaming unzipped barcodes file.'
    else:
        # read_count update
        if barcodes_row_count == barcodes_count_frommatrix:
            results['barcode_count'] = barcodes_row_count
        else:
            errors['barcode_count_discrepancy'] = 'Barcode count from matrix ({}) does not match row count in barcodes.tsv ({})'.format(
                barcodes_count_frommatrix, barcodes_row_count)

    shutil.rmtree(tmp_dir)


    mtx_path = 'matrix.mtx.gz'
    feature_path = 'features.tsv.gz'
    barcode_path = 'barcodes.tsv.gz'

    # do the matrix file in the directory
    try:
        matrix_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, mtx_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['matrix_information_extraction'] = 'Failed to extract information from ' + mtx_path

    try:
        line_index = 0
        for encoded_line in matrix_stream.stdout:
            line_index += 1
            if line_index == 3:
                count_summary = encoded_line.decode('utf-8').strip().split()
                features_count_frommatrix = int(count_summary[0])
                barcodes_count_frommatrix = int(count_summary[1])
                break
    except IOError:
        errors['unzipped_matrix_streaming'] = 'Error occured, while streaming unzipped matrix file.'

    # do the features file in the directory
    try:
        features_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, feature_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['features_information_extraction'] = 'Failed to extract information from ' + feature_path

    features_row_count = 0
    try:
        for encoded_line in features_stream.stdout:
            features_row_count += 1
    except IOError:
        errors['unzipped_features_streaming'] = 'Error occured, while streaming unzipped features file.'
    else:
        # read_count update
        if features_row_count == features_count_frommatrix:
            results['feature_count'] = features_row_count
        else:
            errors['feature_count_discrepancy'] = 'Feature count from matrix ({}) does not match row count in features.tsv ({})'.format(
                features_count_frommatrix, features_row_count)

    # do the barcodes file in the directory
    try:
        barcodes_stream = subprocess.Popen(['gunzip --stdout {}/{}'.format(tmp_dir, barcode_path)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        errors['barcodes_information_extraction'] = 'Failed to extract information from ' + barcode_path

    barcodes_row_count = 0
    try:
        for encoded_line in barcodes_stream.stdout:
            barcodes_row_count += 1
    except IOError:
        errors['unzipped_barcodes_streaming'] = 'Error occured, while streaming unzipped barcodes file.'
    else:
        # read_count update
        if barcodes_row_count == barcodes_count_frommatrix:
            results['barcode_count'] = barcodes_row_count
        else:
            errors['barcode_count_discrepancy'] = 'Barcode count from matrix ({}) does not match row count in barcodes.tsv ({})'.format(
                barcodes_count_frommatrix, barcodes_row_count)


def download_s3_directory(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    dir_path = download_url.replace('s3://{}/'.format(bucket_name), '')

    s3client = boto3.client("s3")
    tmp_dir = 'raw_feature_bc_matrix'
    os.mkdir(tmp_dir)
    for file_name in ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']:
        try:
            s3client.download_file(bucket_name, dir_path + '/' + file_name, '{}/{}'.format(tmp_dir, file_name))
        except subprocess.CalledProcessError as e:
            errors['file not found'] = 'Failed to find file on s3'
        else:
            print(file_name + ' downloaded')


def download_s3_file(job):
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
    file_name = download_url.split('/')[-1]

    job['download_start'] = datetime.now()

    s3client = boto3.client("s3")
    try:
        s3client.download_file(bucket_name, file_path, file_name)
    except subprocess.CalledProcessError as e:
        errors['file not found'] = 'Failed to find file on s3'
    else:
        print(file_name + ' downloaded')
        job['download_stop'] = datetime.now()
        job['download_time'] = job['download_stop'] - job['download_start']


def compare_with_db(job, connection):
    server = connection.server

    file = job['item']
    errors = job['errors']
    results = job['results']
    post_json = {}

    metadata_consistency = []
    metadata_inconsistency = []

    obj_types = file['@type']
    schema = next(iter(set(obj_types) - set(abstract_file_types)))
    schema_url = urljoin(server, 'profiles/' + schema + '/?format=json')
    schema_properties = requests.get(schema_url).json()['properties']

    for key in results.keys():
        if not file.get(key) and schema_properties.get(key):
            post_json[key] = results.get(key)
        elif results.get(key) == file.get(key):
            outcome = '{} consistent ({})'.format(key, results.get(key))
            metadata_consistency.append(outcome)
        elif file.get(key) != None:
            outcome = '{} inconsistent ({}-s3file, {}-submitted)'.format(key, results.get(key), file.get(key))
            metadata_inconsistency.append(outcome)
    if len(metadata_inconsistency) == 0 and schema_properties.get('validated'):
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

    download_url = file.get('s3_uri')
    local_path = download_url.split('/')[-1]

    job['check_start'] = datetime.now()

    # check file size & md5sum
    file_stat = os.stat(local_path)
    results['file_size'] = file_stat.st_size
    # Faster than doing it in Python.
    try:
        output = subprocess.check_output('md5sum-lite {}'.format(local_path),
            shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        errors['md5sum'] = e.output.decode(errors='replace').rstrip('\n')
    else:
        results['md5sum'] = output[:32].decode(errors='replace')
        try:
            int(results['md5sum'], 16)
        except ValueError:
            errors['md5sum'] = output.decode(errors='replace').rstrip('\n')
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
            # May want to replace this with something like:
            # $ cat $local_path | tee >(md5sum >&2) | gunzip | md5sum
            # or http://stackoverflow.com/a/15343686/199100    
            # if gzipped, unzip, get md5sum
            try:
                output = subprocess.check_output(
                    'set -o pipefail; gunzip --stdout {} | md5sum-lite'.format(local_path),
                    shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                errors['content_md5sum'] = e.output.decode(errors='replace').rstrip('\n')
            else:
                results['content_md5sum'] = output[:32].decode(errors='replace')

            # do format-specific validation
            if file.get('file_format') == 'fastq':
                process_fastq_file(job)
                os.remove(local_path)
                print(local_path + ' removed')
    if file.get('file_format') == 'mex':
        process_mexmatrix_file(job)
        shutil.rmtree(tmp_dir)
        print(tmp_dir + ' removed')
    if file.get('file_format') == 'hdf5':
        process_h5matrix_file(job)
        os.remove(local_path)
        print(local_path + ' removed')

    job['check_stop'] = datetime.now()
    job['check_time'] = job['check_stop'] - job['check_start']

    return job


def fetch_files(out, connection=None, query=None, accessions=None, s3_file=None, file_format=None):
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
                if not file_json.get('s3_uri'):
                    blockers.append('s3_uri not submitted')
                    check_me_flag = False
                if file_json.get('validated') == True:
                    blockers.append('already validated')
                    check_me_flag = False
                if check_me_flag == False:
                    out.write(acc + '\t' + file_json.get('s3_uri','') + '\t' + ','.join(blockers) + 3*('\t') + 'n/a' + '\n')
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
        jobs = [{
            'item': {
                'accession': 'not yet submitted',
                's3_uri': s3_file
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
        file_obj.get('s3_uri', ''),
        str(job['errors']),
        str(job['results']),
        str(job['post_json']),
        job.get('patch_result', 'n/a'),
        str(job.get('download_time')),
        str(job.get('check_time'))
    ])
    return tab_report + '\n'


def main():
    args = getArgs()
    if (args.query or args.accessions) and not args.mode:
        sys.exit('ERROR: --mode is required with --query/--accessions')

    arg_count = 0
    for arg in [args.query, args.accessions, args.s3_file]:
        if arg:
            arg_count += 1
    if arg_count != 1:
        sys.exit('ERROR: exactly one of --query, --accessions, --s3-file is required, {} given'.format(arg_count))


    connection = lattice.Connection(args.mode)
    version = '0.9'

    initiating_run = 'STARTING Checkfiles version {}'.format(version)
    print(initiating_run)

    timestr = datetime.now().strftime('%Y_%m_%d-%H_%M_%S')
    out = open('report_{}.txt'.format(timestr), 'w')
    report_headers = '\t'.join([
        'identifier',
        's3_uri',
        'errors',
        'results',
        'json_patch',
        'Lattice patched?',
        'download_time',
        'check_time'
    ])
    out.write(report_headers + '\n')

    jobs = fetch_files(out, connection, args.query, args.accessions, args.s3_file, args.file_format)

    if jobs:
        seq_run_jobs = []
        print('CHECKING {} files'.format(len(jobs)))
        for job in jobs:
            file_obj = job.get('item')
            if file_obj.get('file_format') == 'mex':
                download_s3_directory(job)
            else:
                download_s3_file(job)
            check_file(job)
            compare_with_db(job, connection)
            if job['results'].get('flowcell_details') and file_obj.get('derived_from'):
                seq_run_jobs.append({
                                    'item': file_obj['derived_from'],
                                    'results': {'flowcell_details': job['results']['flowcell_details']},
                                    'errors': {}
                                    })
            if job['post_json'] and not job['errors'] and args.update:
                print('PATCHING {}'.format(file_obj.get('accession')))
                patch = lattice.patch_object(file_obj.get('accession'), connection, job['post_json'])
                job['patch_result'] = patch['status']
            out.write(report(job))
        if seq_run_jobs:
            print('CHECKING {} sequencing_runs'.format(len(seq_run_jobs)))
        for job in seq_run_jobs:
            compare_with_db(job, connection)
            if job['post_json'] and not job['errors'] and args.update:
                print('PATCHING {}'.format(job['item'].get('uuid')))
                patch = lattice.patch_object(job['item'].get('uuid'), connection, job['post_json'])
                job['patch_result'] = patch['status']
            out.write(report(job))


        finishing_run = 'FINISHED Checkfiles at {}'.format(datetime.now())
        print(finishing_run)
        out.close()
    else:
        print('FINISHED No files to check, see report.txt for details')


if __name__ == '__main__':
    main()
