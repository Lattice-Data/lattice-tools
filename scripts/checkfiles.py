import argparse
import datetime
import time
import os.path
import json
import sys
import subprocess
import re
import requests
import copy
import functools
import multiprocessing
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
    parser.add_argument('--local-file',
                        help="path to a single local file to check")
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
    "bam",
    "bed",
    "bedpe",
    "fastq",
    "gff",
    "gtf",
    "tar",
    "txt",
    "wig",
    "vcf",
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
    with open(path, 'rb') as f:
        magic_number = f.read(2)
    return magic_number == b'\x1f\x8b'


def process_illumina_read_name_pattern(read_name,
                                       read_numbers_set,
                                       signatures_set,
                                       signatures_no_barcode_set,
                                       srr_flag):
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


def process_illumina_prefix(read_name,
                                signatures_set,
                                old_illumina_current_prefix,
                                read_numbers_set,
                                srr_flag):
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


def process_read_name_line(read_name_line,
                           old_illumina_current_prefix,
                           read_numbers_set,
                           signatures_no_barcode_set,
                           signatures_set,
                           read_lengths_dictionary,
                           errors, srr_flag, read_name_details):
    read_name = read_name_line.strip()
    if read_name_details:
        #extract fastq signature parts using read_name_detail
        read_name_array = re.split(r'[:\s]', read_name)

        flowcell = read_name_array[read_name_details['flowcell_id_location']]
        lane_number = read_name_array[read_name_details['lane_id_location']]
        if not read_name_details.get('read_number_location'):
            read_number = "1"
        else:
            read_number = read_name_array[read_name_details['read_number_location']]
        read_numbers_set.add(read_number)
        
        if not read_name_details.get('barcode_location'):
            barcode_index = ''
        else:
            barcode_index = read_name_array[read_name_details['barcode_location']]
        
        signatures_set.add(
            flowcell + ':' + lane_number + ':' +
            read_number + ':' + barcode_index + ':')
        signatures_no_barcode_set.add(
            flowcell + ':' + lane_number + ':' +
            read_number + ':')
    else:
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


def process_fastq_file(job, fastq_data_stream, session, url):
    item = job['item']
    errors = job['errors']
    result = job['result']

    read_name_details = get_read_name_details(job.get('@id'), errors, session, url)

    read_numbers_set = set()
    signatures_set = set()
    signatures_no_barcode_set = set()
    read_lengths_dictionary = {}
    read_count = 0
    old_illumina_current_prefix = 'empty'
    try:
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
                            errors, False,
                            read_name_details)
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


def process_read_lengths(read_lengths_dict,
                         lengths_list,
                         submitted_read_length,
                         read_count,
                         threshold_percentage,
                         errors_to_report,
                         result):
    reads_quantity = sum([count for length, count in read_lengths_dict.items()
                          if (submitted_read_length - 2) <= length <= (submitted_read_length + 2)])
    if ((threshold_percentage * read_count) > reads_quantity):
        errors_to_report['read_length'] = \
            'in file metadata the read_length is {}, '.format(submitted_read_length) + \
            'however the uploaded fastq file contains reads of following length(s) ' + \
            '{}. '.format(', '.join(map(str, lengths_list)))


def create_a_list_of_barcodes(details):
    barcodes = set()
    for entry in details:
        barcode = entry.get('barcode')
        lane = entry.get('lane')
        if lane and barcode:
            barcodes.add((lane, barcode))
    return barcodes


def get_read_name_details(job_id, errors, session, url):
    query = job_id +'?datastore=database&frame=object&format=json'
    try:
        r = session.get(urljoin(url, query))
    except requests.exceptions.RequestException as e:
        errors['lookup_for_read_name_detaisl'] = ('Network error occured, while looking for '
                                                  'file read_name details on the portal. {}').format(str(e))
    else:
        details = r.json().get('read_name_details')
        if details:
            return details


def check_file(connection, file):
    download_url = file.get('s3_uri', file.get('file_path'))
    if download_url.startswith('s3:'):
        local_path = os.path.join('/s3', download_url[len('s3://'):])
    else:
        local_path = download_url

    result = {}
    errors = {}

    try:
        file_stat = os.stat(local_path)
    #  When file is not on S3 we are getting FileNotFoundError
    except FileNotFoundError:
        errors['file_not_found'] = (
            'File not found, check the file path, likely in the s3_uri field.'
        )
        print(errors)
    #  Happens when there is S3 connectivity issue: "OSError: [Errno 107] Transport endpoint is not connected"
    except OSError:
        errors['file_check_skipped_due_to_s3_connectivity'] = (
            'File check was skipped due to temporary S3 connectivity issues'
        )
        print(errors)
    else:
        result['file_size'] = file_stat.st_size
        result['last_modified'] = datetime.datetime.utcfromtimestamp(
            file_stat.st_mtime).isoformat() + 'Z'
        print(result)
        quit()
        # Faster than doing it in Python.
        try:
            output = subprocess.check_output(
                ['md5sum', local_path], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            errors['md5sum'] = e.output.decode(errors='replace').rstrip('\n')
        else:
            result['md5sum'] = output[:32].decode(errors='replace')
            try:
                int(result['md5sum'], 16)
            except ValueError:
                errors['md5sum'] = output.decode(errors='replace').rstrip('\n')
            if result['md5sum'] != item['md5sum']:
                errors['md5sum'] = \
                    'checked %s does not match item %s' % (result['md5sum'], item['md5sum'])
        try:
            is_gzipped = is_path_gzipped(local_path)
        except Exception as e:
            return job
        else:
            if item['file_format'] not in GZIP_TYPES:
                if is_gzipped:
                    errors['gzip'] = 'Expected un-gzipped file'
            elif not is_gzipped:
                errors['gzip'] = 'Expected gzipped file'
            else:
                # May want to replace this with something like:
                # $ cat $local_path | tee >(md5sum >&2) | gunzip | md5sum
                # or http://stackoverflow.com/a/15343686/199100
                try:
                    output = subprocess.check_output(
                        'set -o pipefail; gunzip --stdout %s | md5sum' % quote(local_path),
                        shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    errors['content_md5sum'] = e.output.decode(errors='replace').rstrip('\n')

                if item['file_format'] == 'fastq':
                    try:
                        process_fastq_file(job,
                                        subprocess.Popen(['gunzip --stdout {}'.format(
                                                            local_path)],
                                                            shell=True,
                                                            executable='/bin/bash',
                                                            stdout=subprocess.PIPE),
                                        session, url)
                    except subprocess.CalledProcessError as e:
                        errors['fastq_information_extraction'] = 'Failed to extract information from ' + \
                                                                local_path
        if errors:
            errors['gathered information'] = 'Gathered information about the file was: {}.'.format(
                str(result))

        return job


def fetch_files(connection, out, query=None, accessions=None, s3_file=None, local_file=None):
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
    elif query:
        query_url = urljoin(connection.server, query.replace('report', 'search') + '&format=json&limit=all&field=accession')
        r = requests.get(query_url, auth=connection.auth)
        try:
            r.raise_for_status()
        except requests.HTTPError:
            return
        else:
            ACCESSIONS = [x['accession'] for x in r.json()['@graph']]

    # checkfiles on a file that is not in the Lattice database but is at s3
    elif s3_file:
        # NEED TO CHECK IF FILE IS ACCESSIBLE
        file_objects = [
            {
            'accession': 'not yet submitted',
            's3_uri': s3_file
            }
        ]
        return file_objects

    # checkfiles on a file that is not in the Lattice database and is local
    elif local_file:
        # NEED TO CHECK IF FILE IS ACCESSIBLE
        file_objects = [
            {
            'accession': 'not yet submitted',
            'file_path': local_file
            }
        ]
        return file_objects

    file_objects = []
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
                file_objects.append(file_json)
    return file_objects


def main():
    args = getArgs()
    if not args.mode:
        sys.exit('ERROR: --mode is required')
    connection = Connection(args.mode)

    arg_count = 0
    for arg in [args.query, args.accessions, args.s3_file, args.local_file]:
        if arg:
            arg_count += 1
    if arg_count != 1:
        sys.exit('ERROR: exactly one of --query, --accessions, --s3-file --local-file is required, {} given'.format(arg_count))

    version = '0.9'

    initiating_run = 'STARTING Checkfiles version {} on {} at {}'.format(
        version, connection.server, datetime.datetime.now())
    print(initiating_run)

    out = open('report.txt', 'w')
    report_headers = '\t'.join([
        'accession',
        'errors',
        's3_uri'
    ])
    out.write(report_headers + '\n')

    file_objects = fetch_files(connection, out, args.query, args.accessions, args.s3_file, args.local_file)
    print(file_objects)

    print('CHECKING {} files'.format(len(file_objects)))
    for file_obj in file_objects:
        check_file(connection, file_obj)
        tab_report = '\t'.join([
            file_obj.get('accession', 'UNKNOWN'),
            str(job['errors']),
            file_obj.get('s3_uri', '')
        ])
        out.write(tab_report + '\n')

    finishing_run = 'FINISHED Checkfiles at {}'.format(datetime.datetime.now())
    print(finishing_run)
    out.close()


if __name__ == '__main__':
    main()
