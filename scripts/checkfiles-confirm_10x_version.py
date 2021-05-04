import argparse
import boto3
import hashlib
import h5py
import json
import lattice
import logging
import os
import re
import requests
import shutil
import socket
import subprocess
import sys
import tables
import zlib
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
    parser.add_argument('--s3-file',
                        help="path to a file at s3 to check, comma separated or a file containing a list of file accessions to check")
    parser.add_argument('--ext-file',
                        help="path to a file elsewhere to check")
    parser.add_argument('--file-format',
                        help='the specified file format if an s3-file or local-file is being checked')
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


def is_path_gzipped(path):
    try:
        f = open(path, 'rb')
    except IsADirectoryError:
        return 'IsADirectoryError'
    else:
        magic_number = f.read(2)
        return magic_number == b'\x1f\x8b'


def process_fastq_file(job, v2_list, v3_list):
    logging.info('Getting fastq metadata')
    item = job['item']
    errors = job['errors']
    results = job['results']

    download_url = item.get('s3_uri', item.get('external_uri'))
    local_path = download_url.split('/')[-1]

    barcode_list = []
    barcode_counts = {
        'v2': 0,
        'v3': 0,
        'both': 0,
        'neither': 0
    }
    assigned_count = 0
    assigned_limit = 100
    read_count = 0
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
                if line_index == 2:
                    read_count += 1
                    barcode = line.strip()[0:16]
                    if barcode not in barcode_list:
                        v2 = False
                        v3 = False
                        if barcode in v2_list:
                            v2 = True
                        if barcode in v3_list:
                            v3 = True
                        if v2 == True and v3 == True:
                            barcode_counts['both'] += 1
                        elif v2 == True:
                            barcode_counts['v2'] += 1
                            assigned_count += 1
                        elif v3 == True:
                            barcode_counts['v3'] += 1
                            assigned_count += 1
                        else:
                            barcode_counts['neither'] += 1
                        barcode_list.append(barcode)
                        if assigned_count == assigned_limit:
                            break

                line_index = line_index % 4
    except IOError:
        errors['unzipped_fastq_streaming'] = 'Error occured, while streaming unzipped fastq.'
    else:
        if barcode_counts['neither'] != 0:
            errors['barcode not assigned'] = '{} barcodes not found in v2 or v3 barcode lists'.format(str(barcode_counts['neither']))

        if barcode_counts['v2'] > 0 and barcode_counts['v3'] == 0:
            results['10x version'] = 'v2'
        elif barcode_counts['v3'] > 0 and barcode_counts['v2'] == 0:
            results['10x version'] = 'v3'
        else:
            results['10x version'] = 'undetermined'

        barcode_list = []
        for k in sorted(barcode_counts.keys()):
            barcode_list.append({'barcode': k, 'read_count': barcode_counts[k]})

        results['barcode_list'] = ', '.join(map(str, barcode_list))


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
    except subprocess.CalledProcessError as e:
        errors['file not found'] = 'Failed to find file on s3'
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
    #ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
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
        # if it's a schema property currently absent, prepare to patch it
        if not file.get(key) and schema_properties.get(key):
            post_json[key] = results.get(key)
        # if the file information matches the current database metadata, log it
        elif results.get(key) == file.get(key):
            outcome = '{} consistent ({})'.format(key, results.get(key))
            metadata_consistency.append(outcome)
        elif file.get(key) != None:
            #first check for embedded properties, specifically for flowcell_details
            if isinstance(results.get(key), list) and isinstance(file.get(key), list) \
                and len(results.get(key)) == 1 and len(file.get(key)) == 1:
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
            # we have an inconsistency to log
            else:
                outcome = '{} inconsistent ({}-s3file, {}-submitted)'.format(key, results.get(key), file.get(key))
                metadata_inconsistency.append(outcome)
    if len(metadata_inconsistency) == 0 and schema_properties.get('validated') and file.get('validated') != True:
        post_json['validated'] = True

    job['post_json'] = post_json
    if metadata_consistency:
        results['metadata_consistency'] = metadata_consistency
    if metadata_inconsistency:
        errors['metadata_inconsistency'] = metadata_inconsistency


def check_file(job, v2_list, v3_list):
    file = job['item']
    errors = job['errors']
    results = job['results']

    download_url = file.get('s3_uri', file.get('external_uri'))
    local_path = download_url.split('/')[-1]

    job['check_start'] = datetime.now()

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
            # do format-specific validation
            if file.get('file_format') == 'fastq':
                process_fastq_file(job, v2_list, v3_list)
                os.remove(local_path)
                logging.info(local_path + ' removed')

    job['check_stop'] = datetime.now()
    job['check_time'] = job['check_stop'] - job['check_start']

    return job


def fetch_files(report_out, connection=None, query=None, accessions=None, s3_file=None, ext_file=None, file_format=None):
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
                if check_me_flag == False:
                    uri = file_json.get('s3_uri',file_json.get('external_uri'))
                    if uri == None:
                        uri = ''
                    out = open(report_out, 'a')
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


def import_txt_to_list(s3_uri):
    bucket_name = s3_uri.split('/')[2]
    file_path = s3_uri.replace('s3://{}/'.format(bucket_name), '')
    file_name = s3_uri.split('/')[-1]

    s3client = boto3.client("s3")
    logging.info(file_name + ' downloading')
    try:
        s3client.download_file(bucket_name, file_path, file_name)
    except subprocess.CalledProcessError as e:
        sys.exit('ERROR: {} not found, check uri'.format(file_name))
    else:
        logging.info(file_name + ' downloaded')

    barcode_list = [line.strip() for line in open(file_name, 'r')]
    os.remove(file_name)

    return barcode_list


def main():
    logging.basicConfig(filename='checkfiles.log', level=logging.INFO)
    logging.info('Started')

    args = getArgs()
    if (args.query or args.accessions) and not args.mode:
        sys.exit('ERROR: --mode is required with --query/--accessions')

    arg_count = 0
    for arg in [args.query, args.accessions, args.s3_file]:
        if arg:
            arg_count += 1
    if arg_count != 1:
        sys.exit('ERROR: exactly one of --query, --accessions, --s3-file is required, {} given'.format(arg_count))


    if args.mode:
        connection = lattice.Connection(args.mode)
    else:
        connection = ''

    initiating_run = 'STARTING Checkfiles version'
    logging.info(initiating_run)

    timestr = datetime.now().strftime('%Y_%m_%d-%H_%M_%S')
    report_out = 'report_{}.txt'.format(timestr)
    logging.info('Writing results to {}'.format(report_out))
    report_headers = '\t'.join([
        'identifier',
        'uri',
        'errors',
        'results'
    ])
    with open(report_out, 'w') as out:
        out.write(report_headers + '\n')

    v2_barcodes_uri = 's3://submissions-lattice/cellranger-whitelist/737K-august-2016.txt'
    v3_barcodes_uri = 's3://submissions-lattice/cellranger-whitelist/3M-february-2018.txt'

    jobs = fetch_files(report_out, connection, args.query, args.accessions, args.s3_file, args.ext_file, args.file_format)

    if jobs:
        v2_list = import_txt_to_list(v2_barcodes_uri)
        v3_list = import_txt_to_list(v3_barcodes_uri)
        all_seq_runs = []
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
                check_file(job, v2_list, v3_list)
            out = open(report_out, 'a')
            out.write(report(job))
            out.close()

        finishing_run = 'FINISHED Checkfiles at {}'.format(datetime.now())
        logging.info(finishing_run)
    else:
        logging.info('FINISHED No files to check, see report.txt for details')

    logging.info('Results written to {}'.format(report_out))
    logging.info('Finished')

if __name__ == '__main__':
    main()
