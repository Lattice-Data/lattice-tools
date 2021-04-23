import argparse
import boto3
import glob
import googleapiclient.discovery
import json
import os
from google.cloud import storage
from google.oauth2 import service_account
from datetime import datetime, timezone


# adapted from https://github.com/GoogleCloudPlatform/python-docs-samples/blob/master/storage/transfer_service/aws_request.py
def aws_file_transfer(dataset_id, file_uris):
    session = boto3.Session(profile_name='read-only')
    credentials = session.get_credentials()

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
            'status': 'ENABLED',
            'projectId': 'project_id', # I AM NEEDED
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
                        'accessKeyId': credentials.access_key,
                        'secretAccessKey': credentials.secret_key
                    }
                },
                'ObjectConditions': {
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


def directory_transfer(local_path, gcs_path=None):
    assert os.path.isdir(local_path)
    if not gcs_path:
        gcs_path = 'staging/' + local_path
    for local_file in glob.glob(local_path + '/**'):
        if not os.path.isfile(local_file):
            directory_transfer(local_file, gcs_path + "/" + os.path.basename(local_file))
        else:
            bucket_name = 'broad-dsp-monster-hca-dev-lattice'
            bucket = storage.Client().bucket(bucket_name)
            remote_path = os.path.join(gcs_path, local_file[1 + len(local_path):])
            blob = bucket.blob(remote_path)
            blob.upload_from_filename(local_file)
