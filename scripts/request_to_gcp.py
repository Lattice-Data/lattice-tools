import argparse
import boto3
import glob
import googleapiclient.discovery
import json
import os
import socket
from datetime import datetime, timezone
from ftplib import error_perm, FTP
from google.cloud import storage
from google.oauth2 import service_account


def ftp_file_transfer(dataset_id, file_uris):
    sink_path = 'staging/{}/data/'.format(dataset_id)
    for uri in file_uris:
        ftp_download(uri)
        local_path = uri.split('/')[-1]
        local_file_transfer(local_path, sink_path)
        os.remove(local_path)


def ftp_download(uri):
    ftp_server = uri.split('/')[2]
    ftp = FTP(ftp_server)
    ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 75)
    ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
    ftp.login(user='anonymous', passwd = 'password')

    file_path = uri.replace('ftp://{}/'.format(ftp_server), '')
    file_name = uri.split('/')[-1]

    try:
        ftp.retrbinary('RETR ' + file_path, open(file_name, 'wb').write)
    except error_perm as e:
        errors['file download error'] = e
        os.remove(file_name)
    else:
        ftp.quit()


# adapted from https://github.com/GoogleCloudPlatform/python-docs-samples/blob/master/storage/transfer_service/aws_request.py
def aws_file_transfer(dataset_id, file_uris):
    s3_session = boto3.Session(profile_name='read-only')
    s3_credentials = s3_session.get_credentials()

    d_now = datetime.now(tz=timezone.utc)

    sink_path = 'staging/{}/data/'.format(dataset_id)

    files_per_buck = {}
    for uri in file_uris:
        path = uri.split('/')
        bucket = path[2]
        if files_per_buck.get(bucket):
            files_per_buck[bucket].append('/'.join(path[3:]))
        else:
            files_per_buck[bucket] = ['/'.join(path[3:])]

    """Create a one-time transfer from Amazon S3 to Google Cloud Storage."""
    storagetransfer = googleapiclient.discovery.build('storagetransfer', 'v1')

    for source_bucket,files in files_per_buck.items():
        description = 'Transfer of {} files for dataset {}'.format(len(files),dataset_id)
        transfer_job = {
            'description': description,
            'projectId': 'broad-dsp-monster-hca-dev',
            'status': 'ENABLED',
            'schedule': {
                'scheduleStartDate': {
                    'day': d_now.day,
                    'month': d_now.month,
                    'year': d_now.year
                }
            },
            'transferSpec': {
                'awsS3DataSource': {
                    'bucketName': source_bucket,
                    'awsAccessKey': {
                        'accessKeyId': s3_credentials.access_key,
                        'secretAccessKey': s3_credentials.secret_key
                    }
                },
                'objectConditions': {
                    'includePrefixes': files
                },
                'gcsDataSink': {
                    'bucketName': 'broad-dsp-monster-hca-dev-lattice',
                    'path': sink_path
                }
            }
        }
        result = storagetransfer.transferJobs().create(body=transfer_job).execute()
        print('Returned transferJob: {}'.format(
            json.dumps(result, indent=4)))


def local_dir_transfer(local_path, gcs_path=None):
    if not gcs_path:
        gcs_path = 'staging/' + local_path

    bucket_name = 'broad-dsp-monster-hca-dev-lattice'
    bucket = storage.Client().bucket(bucket_name)

    if os.path.isdir(local_path):
        for local_file in glob.glob(local_path + '/**'):
            if not os.path.isfile(local_file):
                local_dir_transfer(local_file, gcs_path + "/" + os.path.basename(local_file))
            else:
                remote_path = os.path.join(gcs_path, local_file[1 + len(local_path):])
                blob = bucket.blob(remote_path)
                blob.upload_from_filename(local_file)
    elif os.path.isfile(local_path):
        remote_path = os.path.join(gcs_path, local_file[1 + len(local_path):])
        blob = bucket.blob(remote_path)
        blob.upload_from_filename(local_file)


def local_file_transfer(local_file, gcs_path):
    assert os.path.isfile(local_file)

    bucket_name = 'broad-dsp-monster-hca-dev-lattice'
    bucket = storage.Client().bucket(bucket_name)

    remote_path = os.path.join(gcs_path, local_file)
    blob = bucket.blob(remote_path)
    blob.upload_from_filename(local_file)
