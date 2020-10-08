import argparse
import boto3
import datetime
import h5py
import json
import os.path
import re
import requests
import shutil
import subprocess
import sys
from urllib.parse import urljoin
from shlex import quote

EPILOG = '''
Run sanity checks on files.

Example:

    %(prog)s --mode production --output check_files.log https://www.encodeproject.org

This relies on local variables to be defined based on the --mode you provide
to direct the updates to a server and to provide permissions
For example, if specifying --mode production, to make the changes on a local instance,
the following variables need to be defined...
PRODUCTION_KEY, PRODUCTION_SECRET, PRODUCTION_SERVER

For more details:

        %(prog)s --help
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

# should turrn this into a shared doc for all scripts to use
class Connection(object):
    def __init__(self, mode):
        if not (os.environ.get(mode.upper() + '_KEY') 
            and os.environ.get(mode.upper() + '_SECRET') 
            and os.environ.get(mode.upper() + '_SERVER')):
            sys.exit('ERROR: ' + mode.upper() + '_KEY, ' + mode.upper() + '_SECRET, ' + mode.upper() + "_SERVER not all defined. Try 'conda env config vars list' to list existing variables")
        self.authid = os.environ.get(mode.upper() + '_KEY')
        self.authpw = os.environ.get(mode.upper() + '_SECRET')
        self.server = os.environ.get(mode.upper() + '_SERVER')
        if not self.server.endswith('/'):
            self.server += '/'
        self.headers = {'content-type': 'application/json',
                        'accept': 'application/json'}
        self.auth = (self.authid, self.authpw)


GZIP_TYPES = [
    "fastq"
]

read_name_prefix = re.compile(
    '^(@[a-zA-Z\d]+[a-zA-Z\d_-]*:[a-zA-Z\d-]+:[a-zA-Z\d_-]' +
    '+:\d+:\d+:\d+:\d+)$')

read_name_pattern = re.compile(
    '^(@[a-zA-Z\d]+[a-zA-Z\d_-]*:[a-zA-Z\d-]+:[a-zA-Z\d_-]' +
    '+:\d+:\d+:\d+:\d+[\s_][123]:[YXN]:[0-9]+:([ACNTG\+]*|[0-9]*))$'
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
    read_name_array = re.split(r'[:\s_]', read_name)
    flowcell = read_name_array[2]
    lane_number = read_name_array[3]
    if srr_flag:
        read_number = list(read_numbers_set)[0]
    else:
        read_number = read_name_array[-4]
        read_numbers_set.add(read_number)
    barcode_index = read_name_array[-1]
    signatures_set.add(
        flowcell + ':' + lane_number + ':' +
        read_number + ':' + barcode_index + ':')
    signatures_no_barcode_set.add(
        flowcell + ':' + lane_number + ':' +
        read_number + ':')


def process_illumina_prefix(read_name, signatures_set, old_illumina_current_prefix, read_numbers_set, srr_flag):
    if srr_flag:
        read_number = list(read_numbers_set)[0]
    else:
        read_number = '1'
        read_numbers_set.add(read_number)
    read_name_array = re.split(r':', read_name)

    if len(read_name_array) > 3:
        flowcell = read_name_array[2]
        lane_number = read_name_array[3]

        prefix = flowcell + ':' + lane_number
        if prefix != old_illumina_current_prefix:
            old_illumina_current_prefix = prefix

            signatures_set.add(
                flowcell + ':' + lane_number + ':' +
                read_number + '::' + read_name)

    return old_illumina_current_prefix


def process_read_name_line(read_name_line, old_illumina_current_prefix, read_numbers_set, signatures_no_barcode_set, signatures_set, read_lengths_dictionary, errors, srr_flag):
    read_name = read_name_line.strip()
    words_array = re.split(r'\s', read_name)
    if read_name_pattern.match(read_name) is None:
        if srr_read_name_pattern.match(read_name.split(' ')[0]) is not None:
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
    else:
        process_illumina_read_name_pattern(
            read_name,
            read_numbers_set,
            signatures_set,
            signatures_no_barcode_set,
            srr_flag)

    return old_illumina_current_prefix


def process_sequence_line(sequence_line, read_lengths_dictionary):
    length = len(sequence_line.strip())
    if length not in read_lengths_dictionary:
        read_lengths_dictionary[length] = 0
    read_lengths_dictionary[length] += 1


def process_barcodes(signatures_set):
    set_to_return = set()
    flowcells_dict = {}
    for entry in signatures_set:
        (f, l, r, b, rest) = entry.split(':')
        if (f, l, r) not in flowcells_dict:
            flowcells_dict[(f, l, r)] = {}
        if b not in flowcells_dict[(f, l, r)]:
            flowcells_dict[(f, l, r)][b] = 0
        flowcells_dict[(f, l, r)][b] += 1
    for key in flowcells_dict.keys():
        barcodes_dict = flowcells_dict[key]
        total = 0
        for b in barcodes_dict.keys():
            total += barcodes_dict[b]
        for b in barcodes_dict.keys():
            if ((float(total)/float(barcodes_dict[b])) < 100):
                set_to_return.add(key[0] + ':' +
                                  key[1] + ':' +
                                  key[2] + ':' +
                                  b + ':')
    return set_to_return


def process_read_lengths(read_lengths_dict, lengths_list, submitted_read_length, read_count, threshold_percentage, errors_to_report, result):
    reads_quantity = sum([count for length, count in read_lengths_dict.items()
                          if (submitted_read_length - 2) <= length <= (submitted_read_length + 2)])
    if ((threshold_percentage * read_count) > reads_quantity):
        errors_to_report['read_length'] = \
            'in file metadata the read_length is {}, '.format(submitted_read_length) + \
            'however the uploaded fastq file contains reads of following length(s) ' + \
            '{}. '.format(', '.join(map(str, lengths_list)))


def process_fastq_file(job):
    item = job['item']
    errors = job['errors']
    result = job['result']

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
                process_sequence_line(line, read_lengths_dictionary)

            line_index = line_index % 4
    except IOError:
        errors['unzipped_fastq_streaming'] = 'Error occured, while streaming unzipped fastq.'
    else:

        # read_count update
        result['read_count'] = read_count

        # read1/read2
        if len(read_numbers_set) > 1:
            errors['inconsistent_read_numbers'] = \
                'fastq file contains mixed read numbers ' + \
                '{}.'.format(', '.join(sorted(list(read_numbers_set))))

        # read_length
        read_lengths_list = []
        for k in sorted(read_lengths_dictionary.keys()):
            read_lengths_list.append((k, read_lengths_dictionary[k]))

        if 'read_length' in item and item['read_length'] > 2:
            process_read_lengths(read_lengths_dictionary,
                                 read_lengths_list,
                                 item['read_length'],
                                 read_count,
                                 0.9,
                                 errors,
                                 result)
        else:
            errors['read_length'] = 'no specified read length in the uploaded fastq file, ' + \
                                    'while read length(s) found in the file were {}. '.format(
                                        ', '.join(map(str, read_lengths_list)))
        # signatures
        signatures_for_comparison = set()
        is_UMI = False
        if 'flowcell_details' in item and len(item['flowcell_details']) > 0:
            for entry in item['flowcell_details']:
                if 'barcode' in entry and entry['barcode'] == 'UMI':
                    is_UMI = True
                    break
        if old_illumina_current_prefix == 'empty' and is_UMI:
            for entry in signatures_no_barcode_set:
                signatures_for_comparison.add(entry + 'UMI:')
        else:
            if old_illumina_current_prefix == 'empty' and len(signatures_set) > 100:
                signatures_for_comparison = process_barcodes(signatures_set)
                if len(signatures_for_comparison) == 0:
                    for entry in signatures_no_barcode_set:
                        signatures_for_comparison.add(entry + 'mixed:')

            else:
                signatures_for_comparison = signatures_set

        result['fastq_signature'] = sorted(list(signatures_for_comparison))

    os.remove(local_path)


def process_h5matrix_file(job):
    item = job['item']
    errors = job['errors']
    result = job['result']

    download_url = item.get('s3_uri')
    local_path = download_url.split('/')[-1]

    # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
    # https://docs.h5py.org/en/stable/
    f = h5py.File(local_path, 'r')
    for k in list(f.keys()):
        dset = f[k]
        result['barcode_count'] = dset['barcodes'].shape[0]
        result['features'] = dset['features']['id'].shape[0]

    os.remove(local_path)


def process_matrix_file(job):
    item = job['item']
    errors = job['errors']
    result = job['result']

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
            result['feature_count'] = features_row_count
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
            result['barcode_count'] = barcodes_row_count
        else:
            errors['barcode_count_discrepancy'] = 'Barcode count from matrix ({}) does not match row count in barcodes.tsv ({})'.format(
                barcodes_count_frommatrix, barcodes_row_count)

    shutil.rmtree(tmp_dir)


def download_s3_directory(job):
    item = job['item']
    errors = job['errors']
    result = job['result']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    dir_path = download_url.replace('s3://{}/'.format(bucket_name), '')

    s3client = boto3.client("s3")
    tmp_dir = 'raw_feature_bc_matrix'
    os.mkdir(tmp_dir)
    for file_name in ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']:
        try:
            s3client.download_file(bucket_name, dir_path + '/' + file_name, '{}/{}'.format(tmp_dir, file_name))
            '''
            subprocess.Popen(['aws s3 cp {} ./'.format(download_url)],
                shell=True, executable='/bin/bash', stdout=subprocess.PIPE)'''
        except subprocess.CalledProcessError as e:
            errors['file not found'] = 'Failed to find file on s3'
        else:
            #once it has successfully downloaded, then go back
            print(file_name + ' downloaded')


def download_s3_file(job):
    item = job['item']
    errors = job['errors']
    result = job['result']

    download_url = item.get('s3_uri')
    bucket_name = download_url.split('/')[2]
    file_path = download_url.replace('s3://{}/'.format(bucket_name), '')
    file_name = download_url.split('/')[-1]

    s3client = boto3.client("s3")
    try:
        s3client.download_file(bucket_name, file_path, file_name)
        '''
        subprocess.Popen(['aws s3 cp {} ./'.format(download_url)],
            shell=True, executable='/bin/bash', stdout=subprocess.PIPE)'''
    except subprocess.CalledProcessError as e:
        errors['file not found'] = 'Failed to find file on s3'
    else:
            #once it has successfully downloaded, then go back
            print(file_name + ' downloaded')


def check_file(job):
    file = job['item']
    errors = job['errors']
    result = job['result']

    download_url = file.get('s3_uri')
    local_path = download_url.split('/')[-1]

    # check file size & md5sum
    file_stat = os.stat(local_path)
    result['file_size'] = file_stat.st_size
    result['last_modified'] = datetime.datetime.utcfromtimestamp(
        file_stat.st_mtime).isoformat() + 'Z'
    # Faster than doing it in Python.
    try:
        output = subprocess.check_output('md5sum-lite {}'.format(local_path),
            shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        errors['md5sum'] = e.output.decode(errors='replace').rstrip('\n')
    else:
        result['md5sum'] = output[:32].decode(errors='replace')
        try:
            int(result['md5sum'], 16)
        except ValueError:
            errors['md5sum'] = output.decode(errors='replace').rstrip('\n')
        if file.get('md5sum') and result['md5sum'] != file.get('md5sum'):
            errors['md5sum'] = \
                'checked %s does not match item %s' % (result['md5sum'], file.get('md5sum'))
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
                result['content_md5sum'] = output[:32].decode(errors='replace')

            # do format-specific validation
            if file.get('file_format') == 'fastq':
                process_fastq_file(job)
    if file.get('file_format') == 'mex':
        process_matrix_file(job)
    if file.get('file_format') == 'h5':
        process_h5matrix_file(job)

    return job


def fetch_files(out, mode=None, query=None, accessions=None, s3_file=None, file_format=None):
    if accessions or query:
        connection = Connection(mode)
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
            query_url = urljoin(connection.server, query.replace('report', 'search') + '&format=json&limit=all&field=accession')
            r = requests.get(query_url, auth=connection.auth)
            try:
                r.raise_for_status()
            except requests.HTTPError:
                return
            else:
                ACCESSIONS = [x['accession'] for x in r.json()['@graph']]

        jobs = []
        for acc in ACCESSIONS:
            item_url = urljoin(connection.server, acc + '/?frame=object')
            fileObject = requests.get(item_url, auth=connection.auth)
            if not fileObject.ok:
                errors['file_HTTPError'] = ('HTTP error: unable to get file object')
            else:
                file_json = fileObject.json()
                if file_json.get('no_file_available') == True:
                    out.write(acc + '\t' + ' marked as no_file_available' + '\n')
                elif not file_json.get('s3_uri'):
                    out.write(acc + '\t' + ' s3_uri not submitted' + '\n')
                else:
                    job = {
                        'item': file_json,
                        'result': {},
                        'errors': {}
                    }
                    jobs.append(job)
        return jobs

    # checkfiles on a file that is not in the Lattice database but is at s3
    elif s3_file:
        # NEED TO CHECK IF FILE IS ACCESSIBLE
        jobs = [{
            'item': {
                'accession': 'not yet submitted',
                's3_uri': s3_file
            },
            'result': {},
            'errors': {}
        }]
        if file_format:
            jobs[0]['item']['file_format'] = file_format
        return jobs


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

    version = '0.9'

    initiating_run = 'STARTING Checkfiles version {} at {}'.format(
        version, datetime.datetime.now())
    print(initiating_run)

    out = open('report.txt', 'w')
    report_headers = '\t'.join([
        'accession',
        'errors',
        'results',
        's3_uri'
    ])
    out.write(report_headers + '\n')

    jobs = fetch_files(out, args.mode, args.query, args.accessions, args.s3_file, args.file_format)

    if jobs:
        print('CHECKING {} files'.format(len(jobs)))

        for job in jobs:
            if job['item'].get('file_format') == 'mex':
                download_s3_directory(job)
            else:
                download_s3_file(job)
            check_file(job)
            file_obj = job.get('item')
            tab_report = '\t'.join([
                file_obj.get('accession', 'UNKNOWN'),
                str(job['errors']),
                str(job['result']),
                file_obj.get('s3_uri', '')
            ])
            out.write(tab_report + '\n')

        finishing_run = 'FINISHED Checkfiles at {}'.format(datetime.datetime.now())
        print(finishing_run)
        out.close()


if __name__ == '__main__':
    main()
