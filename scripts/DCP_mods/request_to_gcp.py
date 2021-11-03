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


def list_bucket_contents():
    files = []
    client = storage.Client()
    bucket = storage.Bucket(client, 'broad-dsp-monster-hca-prod-lattice')
    all_blobs = list(client.list_blobs(bucket))
    for blob in all_blobs:
        files.append(blob.name)
    return files


def ftp_file_transfer(dataset_id, file_uris):
    sink_path = 'staging/{}/data/'.format(dataset_id)
    gcp_files = list_bucket_contents()
    for uri in file_uris:
        local_path = uri.split('/')[-1]
        if (sink_path + local_path) not in gcp_files:
            ftp_download(uri)
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

    files_per_buckpath = {}
    for uri in file_uris:
        path = uri.split('/')
        aws_bucket = path[2]
        aws_path = '/'.join(path[3:-1])
        aws_file = path[-1]
        aws_buckpath = (aws_bucket, aws_path)
        if files_per_buckpath.get(aws_buckpath):
            files_per_buckpath[aws_buckpath].append(aws_file)
        else:
            files_per_buckpath[aws_buckpath] = [aws_file]

    """Create a one-time transfer from Amazon S3 to Google Cloud Storage."""
    storagetransfer = googleapiclient.discovery.build('storagetransfer', 'v1')

    for source,files in files_per_buckpath.items():
        source_bucket = source[0]
        source_path = source[1] + '/'
        description = 'Transfer of {} files for dataset {}'.format(len(files),dataset_id)
        transfer_job = {
            'description': description,
            'projectId': 'mystical-slate-284720',
            'status': 'ENABLED',
            'schedule': {
                'scheduleStartDate': {
                    'day': d_now.day,
                    'month': d_now.month,
                    'year': d_now.year
                },
                'scheduleEndDate': {
                    'day': d_now.day,
                    'month': d_now.month,
                    'year': d_now.year
                }
            },
            'transferSpec': {
                'awsS3DataSource': {
                    'bucketName': source_bucket,
                    'path': source_path,
                    'awsAccessKey': {
                        'accessKeyId': s3_credentials.access_key,
                        'secretAccessKey': s3_credentials.secret_key
                    }
                },
                'objectConditions': {
                    'includePrefixes': files
                },
                'gcsDataSink': {
                    'bucketName': 'broad-dsp-monster-hca-prod-lattice',
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

    bucket_name = 'broad-dsp-monster-hca-prod-lattice'
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

    bucket_name = 'broad-dsp-monster-hca-prod-lattice'
    bucket = storage.Client().bucket(bucket_name)

    remote_path = os.path.join(gcs_path, local_file)
    blob = bucket.blob(remote_path)
    blob.upload_from_filename(local_file)
