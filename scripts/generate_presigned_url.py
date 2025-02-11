import argparse
import sys
import boto3
import botocore


EPILOG = f"""
Script to create an upload presigned URL CURL command that uploads a file to S3 bucket.

Examples:

	python %(prog)s -f file_to_upload.h5ad -b bucketname -n "dir1/dir2/uploaded_file_name.h5ad" -t timeout_limit
	python %(prog)s -f file_to_upload.h5ad -b bucketname

"""

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, 
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
    	"--file",
    	"-f",
    	help="file that needs to be uploaded to S3"
    )
    parser.add_argument(
    	"--bucket",
    	"-b",
    	help="name of bucket on S3"
    )
    parser.add_argument(
    	"--name",
    	"-n",
    	help="name of destination file on S3",
    	default = None
    )
    parser.add_argument(
    	"--time",
    	"-t",
    	help="length of time in seconds that presigned URL will last",
    	default = 3600
    )
    args = parser.parse_args()
    if len(sys.argv) < 2:
    	parser.print_help()
    	sys.exit()
    return args

ARGS = getArgs()

if __name__ == '__main__':
	s3 = boto3.client('s3')
	
	if not ARGS.name:
		ARGS.name = ARGS.file

	# Check to see if there is an existing bucket or file already in bucket
	try:
		s3.head_bucket(Bucket=ARGS.bucket)
	except botocore.exceptions.ClientError:
		print("There is no such bucket: {}".format(ARGS.bucket))
		sys.exit()
	else:
		try:
			s3.head_object(Bucket = ARGS.bucket, Key=ARGS.name)
		except botocore.exceptions.ClientError:
			url = s3.generate_presigned_url('put_object', Params={'Bucket':ARGS.bucket,'Key':ARGS.name}, ExpiresIn=ARGS.time, HttpMethod='PUT')
			print("curl --request PUT --upload-file {} '{}'".format(ARGS.file, url))
		else:
			url = s3.generate_presigned_url('put_object', Params={'Bucket':ARGS.bucket,'Key':ARGS.name}, ExpiresIn=ARGS.time, HttpMethod='PUT')
			print("WARNING: File already exists on S3: {}/{}".format(ARGS.bucket,ARGS.name))
			print("curl --request PUT --upload-file {} '{}'".format(ARGS.file, url))



