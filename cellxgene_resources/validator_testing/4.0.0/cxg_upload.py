import os
import sys


scc_repo_loc = os.path.expanduser('~/GitClones/CZI/')
sys.path.append(os.path.abspath(scc_repo_loc + 'single-cell-curation/notebooks/curation_api/python/'))

from src.utils.config import set_api_access_config
from src.dataset import create_dataset,upload_local_datafile


api_key_file_path = os.path.expanduser('~/Documents/keys/cxg-api-key-dev.txt')
set_api_access_config(api_key_file_path, env='dev')

def upload(file):
    collection_id = 'aec7f284-3cce-48f2-b256-af685dd40a38'
    d = create_dataset(collection_id)
    upload_local_datafile(file, collection_id, d)
