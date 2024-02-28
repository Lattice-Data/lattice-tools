import cellxgene_census
import concurrent.futures
import h5py
import json
import os
import pandas as pd
import requests
import sys
from anndata._io.specs import read_elem
from urllib.parse import quote
from time import perf_counter, sleep

# set repo and api key paths, might move this later
czi_repo_loc = os.path.expanduser('~/GitClones/CZI/')
api_key_file_path = os.path.expanduser('~/Documents/keys/cxg-api-key.txt')
sys.path.append(os.path.abspath(czi_repo_loc + 'single-cell-curation/notebooks/curation_api/python/'))

from src.utils.config import set_api_access_config
from src.collection import get_collections, get_collection
from src.dataset import get_dataset,get_datasets

# importing rna assays supported in census, might just hardcode instead
# from the census builder directory so requires python 3.11 and other dependencies
sys.path.append(os.path.abspath(czi_repo_loc + 'cellxgene-census/tools/cellxgene_census_builder/src/'))
from cellxgene_census_builder.build_soma.globals import RNA_SEQ

set_api_access_config(api_key_file_path)

# make sure czi single-cell-curation branch set to lattice/dev-ont-migration to get json files from their repo
mig_loc = czi_repo_loc + 'single-cell-curation/cellxgene_schema_cli/cellxgene_schema/'


class ApiData:
    def __init__(self, from_file: bool = False):
        self.from_file = from_file

        self.private_collections = None
        self.private_datasets = None
        self.public_collections = None
        self.public_datasets = None
        self.current_dev_terms = None

        self.private_collection_ids = None

        self.fetch_api_data(self.from_file)

    def fetch_api_data(self, from_file):
        if from_file:
            print('Loading CXG API JSONs from file...')
            self._load_from_file()
            print('CXG API JSONs loaded successfully')
        else:
            print('Loading CXG API JSONs from website...')
            self._load_from_czi()
            print('CXG API JSONs loaded successfully')


        print(f"{len(self.public_collections)} Public Collections")
        print(f"{len(self.private_collection_ids)} Private Collections")

        print('Loading current dev terms from website...')
        self._fetch_dev_terms()
        print('Current dev terms loaded successfully')

    def save_api_jsons(self):
        jsons = {
            'api_private_datasets.json': self.private_datasets,
            'api_private_collections.json': self.private_collections,
            'api_public_datasets.json': self.public_datasets,
            'api_public_collections.json': self.public_collections
        }

        for file_name, json_attr in jsons.items():
            with open(file_name, 'w', encoding='utf-8') as f:
                json.dump(json_attr, f, ensure_ascii=False, indent=2)

    def _fetch_dev_terms(self):
        self.current_dev_terms = ['unknown']
        
        for ont in ['hsapdv', 'mmusdv']:
            url = f'https://www.ebi.ac.uk/ols4/api/ontologies/{ont}/terms?obsoletes=false&size=500'
            r = requests.get(url).json()
            terms = [t['obo_id'] for t in r['_embedded']['terms']]
            self.current_dev_terms.extend(terms)

    def _load_from_file(self):
        try:
            self.public_collections = json.load(open('api_public_collections.json'))
            self.public_datasets = json.load(open('api_public_datasets.json'))
            self.private_collections = json.load(open('api_private_collections.json'))
            self.private_datasets = json.load(open('api_private_datasets.json'))
        except FileNotFoundError as e:
            print(f'File not found: {e}')

        self.private_collection_ids = {c['collection_id'] for c in self.private_collections if not c.get('revision_of')}

    def _load_from_czi(self):
            self.public_collections = get_collections()
            self.public_datasets = get_datasets()
            self.private_collections = get_collections(visibility='PRIVATE')

            self.private_collection_ids = {c['collection_id'] for c in self.private_collections if not c.get('revision_of')}
            self.private_datasets = [get_collection(c) for c in self.private_collection_ids]
        

class CensusData:
    def __init__(self, census_version: str = 'latest'):
        self.census_version = census_version

        self._collection_dict = None
        self.obs_df = None

        self.get_data()

    def get_data(self):
        self._get_collection_dict()
        self._get_obs_df()

    def _get_collection_dict(self):
        '''
        The SOMA 'census_data' object does not contain collection_id, so need to generate dict that
        can map collection_id from dataset_id
        returns dict with {dataset_id: collection_id} from census for further identification
        '''
        print('Generating dict for mapping collection_id from dataset_id...')
        with cellxgene_census.open_soma(census_version=self.census_version) as census:
            dataset_info = census['census_info']['datasets'].read().concat().to_pandas()
            
        mapping = dataset_info.groupby('dataset_id')['collection_id'].agg('unique')
        df = pd.DataFrame(mapping)
        df['collection_id'] = df['collection_id'].apply(lambda x: x[0])
        df_dict = df.to_dict('index')
        final_dict = {k: v['collection_id'] for k, v in df_dict.items()}
        print('Successfully created dictionary')
        self._collection_dict = final_dict

    def _get_obs_df(self):
        print('Generating obs df from census with required columns...')
        def get_specieis_obs_df(species):
            with cellxgene_census.open_soma(census_version=self.census_version) as census:
                dev_obs = census['census_data'][species] \
                .obs.read(column_names=['dataset_id', 'development_stage', 'development_stage_ontology_term_id', 'donor_id']) \
                .concat().to_pandas()

            dev_obs['collection_id'] = dev_obs['dataset_id'].map(self._collection_dict)
            return dev_obs
        
        dfs = [get_specieis_obs_df(x) for x in ['mus_musculus', 'homo_sapiens']]
        final_df = pd.concat(dfs)
        print('Successfully created obs_df')
        self.obs_df = final_df


class DevMigrationTool:
    def __init__(self, sheet_id, repo_path, automigrate_json, donor_updates_json):
        self.sheet_id: str = sheet_id
        self.repo_path = repo_path
        self.automigrate_json = automigrate_json
        self.donor_updates_json = donor_updates_json

        self.donor_updates_json_collection_ids = None
        self.extended_donor_updates: dict = None

        self.api_private_collections: dict = None
        self.api_private_datasets: dict = None
        self.api_public_collections: dict = None
        self.api_public_datasets: dict = None

        self.api_private_names_dict: dict = None
        self.private_name_to_collection_id = None
        self.api_public_dataset_ids_dict: dict = None
        self.api_current_dev_terms: list = None
        self.api_private_collection_ids: list = None
        self.private_dataset_versions_dict = None

        self.load_methods = [m for m in dir(DevMigrationTool) if callable(getattr(DevMigrationTool, m)) and m.startswith('load')]

        self.api_data_loaded = False
        self.private_uncovered_terms = None
        self.public_uncovered_terms = None

        self.google_sheet: pd.DataFrame = None

        self.census_collection_dict: dict = None
        self.census_obs_df: pd.DataFrame = None
        self.census_all_uncovered_donors: pd.DataFrame = None
        self.census_assays = set(RNA_SEQ)

        self.private_all_obs_df = None
        self.public_all_obs_df = None
        self.private_all_uncovered_donors = None
        self.private_urls = None
        self.uncovered_census_datasets = []
        self.uncovered_public_datasets = []
        self.public_all_uncovered_donors = None
        self.api_public_collection_ids = None
    
    def load_jsons(self):
        self.automigrate_terms = json.load(open(self.repo_path + self.automigrate_json))
        self.donor_updates_json = json.load(open(self.repo_path + self.donor_updates_json))

        temp_dict = {}
        for k, v in self.donor_updates_json.items():
            if k in self.api_public_collection_ids:
                temp_dict[k] = v
            else:
                value = self.donor_updates_json[k]
                private_collection_id = self.private_name_to_collection_id[k]
                temp_dict[private_collection_id] = value

        self.donor_updates_json_collection_ids = temp_dict

    def load_api_data(self):
        if self.from_file:
            self.api_public_collections = json.load(open('api_public_collections.json'))
            self.api_public_datasets = json.load(open('api_public_datasets.json'))
            self.api_private_collections = json.load(open('api_private_collections.json'))
            self.api_private_datasets = json.load(open('api_private_datasets.json'))
        else:
            self.api_public_collections = get_collections()
            self.api_public_datasets = get_datasets()
            self.api_private_collections = get_collections(visibility='PRIVATE')

        self.api_public_dataset_ids_dict = {d['dataset_id']: d for d in self.api_public_datasets}
        self.public_dataset_id_to_collection_id = {d['dataset_id']: d['collection_id'] for d in self.api_public_datasets}
        self.api_private_names_dict = {c['collection_id']: c['name'] for c in self.api_private_collections if not c.get('revision_of')}
        self.private_name_to_collection_id = {c['name']: c['collection_id'] for c in self.api_private_collections if not c.get('revision_of')}
        self.api_private_collection_ids = [c for c in self.api_private_names_dict]
        self.api_public_collection_ids = [c['collection_id'] for c in self.api_public_collections]

        if not self.api_private_datasets:
            self.api_private_datasets = [get_collection(c) for c in self.api_private_collection_ids]

        dv_id_to_dsi = {}
        for c in self.api_private_datasets:
            datasets = c['datasets']
            temp_dict = {d['dataset_version_id']: d['dataset_id'] for d in datasets}
            dv_id_to_dsi.update(temp_dict)

        public_dataset_dict = {d['dataset_version_id']: d['dataset_id'] for d in self.api_public_datasets} 
        public_dataset_dict.update(dv_id_to_dsi)
        self.api_dataset_version_id_to_dataset_id = public_dataset_dict
        
        print(f"{len(self.api_public_collections)} Public Collections")
        print(f"{len(self.api_private_collection_ids)} Private Collections")

    def load_dev_terms(self):
        self.api_current_dev_terms = ['unknown']
        
        for ont in ['hsapdv', 'mmusdv']:
            url = f'https://www.ebi.ac.uk/ols4/api/ontologies/{ont}/terms?obsoletes=false&size=500'
            r = requests.get(url).json()
            terms = [t['obo_id'] for t in r['_embedded']['terms']]
            self.api_current_dev_terms.extend(terms)

    def fetch_api_data(self, *, from_file: bool = False):
        self.from_file = from_file

        if self.api_data_loaded:
            print('Data already fetched')
            return

        print("Loading data...")
        for method in self.load_methods:
            print(f"Running {method}()...")
            fxn = getattr(DevMigrationTool, method)
            fxn(self)

        print("Data loaded successfully")
        self.api_data_loaded = True

    def update_donor_json(self, tab_names: list[str] = ['public donors 2024', 'private donors 2024']):
        """
        Function to load and update donor_updates.json file. Reads from google 
        sheet tabs and appends info into a dictionary in the format of 
        donor_updates.json:
            {
                "collection_id" OR "name": {
                    "donor_id": "new_dev_term"
                }
            }
        Will update with further docstring later once rest of class is fleshed out
        """
        if not self.api_data_loaded:
            print("Please fetch api data")
            return

        # load donor_updates.json or create new
        donor_updates = self.donor_updates_json

        for tab in tab_names:

            url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}/gviz/tq?tqx=out:csv&sheet={quote(tab)}'
            donor_meta = pd.read_csv(url)
            col_names = ['donor_id', 'new_dev_stage', 'collection_id']

            for _, row in donor_meta[col_names].iterrows():
                c = row['collection_id']
                d = row['donor_id']
                t = row['new_dev_stage']

                # collection_id won't change, title can, so change at last possible moment
                if c in self.api_private_names_dict.keys():
                    c = self.api_private_names_dict[c]
                else:
                    c = c

                if c in donor_updates:
                    donor_updates[c][d] = t
                else:
                    donor_updates[c] = {d:t}
        
        self.extended_donor_updates = donor_updates

    def save_json(self, input_dict, file_name, indent=4):
        with open(file_name, 'w', encoding='utf-8') as f:
            json.dump(input_dict, f, ensure_ascii=False, indent=indent)

    def save_api_jsons(self):
        jsons = {
            'api_private_datasets.json': self.api_private_datasets,
            'api_private_collections.json': self.api_private_collections,
            'api_public_datasets.json': self.api_public_datasets,
            'api_public_collections.json': self.api_public_collections
        }

        for file_name, json_attr in jsons.items():
            with open(file_name, 'w', encoding='utf-8') as f:
                json.dump(json_attr, f, ensure_ascii=False, indent=2)

    def get_google_sheet(self):
        dfs = []
        for tab in ['private donors 2024', 'public donors 2024']:
            url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}/gviz/tq?tqx=out:csv&sheet={quote(tab)}'
            donor_meta = pd.read_csv(url)[['donor_id', 'development_stage_ontology_term_id', 'development_stage', 'collection_id', 'dataset_ids', 'new_dev_stage']]
            donor_meta['tab_name'] = tab
            dfs.append(donor_meta)

        final_df = pd.concat(dfs)
        self.google_sheet = final_df

    def make_public_donor_dev_terms(self):
        pub_collections = {}
        for d in self.api_public_datasets:
            c_id = d['collection_id']
            dev_stages = {t['ontology_term_id']:t['label'] for t in d['development_stage']}
            if c_id in pub_collections:
                pub_collections[c_id]['donor_id'].extend(d['donor_id'])
                pub_collections[c_id]['development_stage'].update(dev_stages)
            else:
                pub_collections[c_id] = {
                    'donor_id': d['donor_id'],
                    'development_stage': dev_stages
                }
        self.public_donor_ids_dev_terms = pub_collections

    def make_private_donor_dev_terms(self):
        private_collections = {}
        for c in self.api_private_datasets:
            c_id = c['collection_id']
            datasets = c['datasets']
            for d in datasets:
                # skip datasets that have not met validation
                if d['processing_status'] != 'SUCCESS':
                    print(f"Current processing error with dataset {d['dataset_version_id']} in collection {c_id}")
                    continue
                dev_stages = {t['ontology_term_id']:t['label'] for t in d['development_stage']}
                if c_id in private_collections:
                    private_collections[c_id]['donor_id'].extend(d['donor_id'])
                    private_collections[c_id]['development_stage'].update(dev_stages)
                else:
                    private_collections[c_id] = {
                        'donor_id': d['donor_id'],
                        'development_stage': dev_stages
                    }
        self.private_donor_ids_dev_terms = private_collections

    def make_private_dataset_versions_dict(self):
        dataset_version_dict = {}
        for c in self.api_private_datasets:
            c_id = c['collection_id']
            datasets = c['datasets']
            for d in datasets:
                # skip datasets that have not met validation
                if d['processing_status'] != 'SUCCESS':
                    print(f"Current processing error with dataset {d['dataset_version_id']} in collection {c_id}")
                    continue
                dataset_version = d['dataset_version_id']
                dataset_version_dict[dataset_version] = c_id
        self.private_dataset_versions_dict = dataset_version_dict

    def uncovered_terms(self, visibility: str = 'public'):
        ''' 
        collection_dict format:
        {collection_id:
            {
                donor_id: list[donors],
                development_stage: dict{dev_term_id: term}
            }    
        }
        returns dict in this format collection_id: {dev_term_id: dev_term}
        param current_terms: list of current dev terms, generated above
        param automigrate_terms: dict from JSON file on single-cell-curation repo, {current_term: new migrated term}
        '''
        if visibility == 'public':
            collection_dict = self.public_donor_ids_dev_terms
        else:
            collection_dict = self.private_donor_ids_dev_terms

        collections = {}
        for k,v in collection_dict.items():
            dev_stages = v['development_stage']
            deprecated = [t for t in dev_stages.keys() if t not in self.api_current_dev_terms and t.startswith('UBERON:') is False]
            not_migrated = [t for t in deprecated if t not in self.automigrate_terms.keys()]
            if not_migrated and k not in self.donor_updates_json_collection_ids:
                print(k)
                for t in not_migrated:
                    entry = {t: dev_stages[t]}
                    # print('--',t,':',dev_stages[t])
                    if k in collections:
                        collections[k].update(entry)
                    else:
                        collections[k] = entry                
                    print(entry)
                print('')
        
        if visibility == 'public':
            self.public_uncovered_terms = collections
        else:
            self.private_uncovered_terms = collections

    def uncovered_census_terms(self, visibility: str):
        if not self.public_uncovered_terms:
            print("Please generate public uncovered terms first with .uncovered_terms('public')")
            return

        if visibility == 'uncovered public':
            uncovered_terms = self.public_uncovered_terms
        elif visibility == 'google sheet all public':
            uncovered_terms = self.create_google_sheet_aggregated_dict(sheet_tab='public donors 2024', key_col='collection_id', value_col='development_stage_ontology_term_id')
            
        public_covered_donors = self.create_google_sheet_aggregated_dict(
            sheet_tab='public donors 2024',
            key_col='collection_id',
            value_col='donor_id'
        )

        print("Datasets not fully in Census that contain possible donors with uncovered dev terms:")
        for collection, uncovered in uncovered_terms.items():
            uncovered_set = {k for k in uncovered.keys()} if visibility == 'uncovered public' else {k for k in uncovered}
            json_info = get_collection(collection)
            datasets = json_info['datasets']
            covered_donors = set(public_covered_donors[collection])
            for dataset in datasets:
                dataset_id = dataset['dataset_id']
                assays = dataset['assay']
                donors = set(dataset['donor_id'])
                dev_stages = {t['ontology_term_id'] for t in dataset['development_stage']}
                assays_set = {a['ontology_term_id'] for a in assays}
                if (len(assays_set - self.census_assays) > 0 and \
                    len(uncovered_set.intersection(dev_stages)) > 0 and \
                    len(donors - covered_donors) > 0
                        ):
                    print(f"Unapproved Census assay(s) {assays_set - self.census_assays} in dataset {dataset_id} from collection: {collection}")
                    if visibility == 'uncovered public':
                        self.uncovered_census_datasets.append(dataset_id)
                    else:
                        self.uncovered_public_datasets.append(dataset_id)

    def estimate_private_dataset_download(self):
        urls_dict = {}
        
        if not self.private_uncovered_terms:
            print("Please generate private uncovered terms first with uncovered_terms('private')")
            return

        # private_dataset_dict = {c: get_collection(c) for c in self.api_private_collection_ids}
        private_dataset_dict = {c['collection_id']: c for c in self.api_private_datasets}

        for cid, terms in self.private_uncovered_terms.items():
            c_info = private_dataset_dict[cid]
            uncovered_dev_terms = [k for k in terms.keys()]
            for dataset in c_info['datasets']:
                if dataset['processing_status'] != 'SUCCESS':
                    print(f"Dataset {dataset['dataset_id']} from Collection {cid} not fully processed")
                    continue
                dev_stages = dataset['development_stage']
                all_dev_terms = [d['ontology_term_id'] for d in dev_stages]
                for t in uncovered_dev_terms:
                    if t in all_dev_terms:
                        url = dataset['assets'][0]['url']
                        filesize = dataset['assets'][0]['filesize']
                        urls_dict[url] = filesize

        print(f'Total number of h5ad files to download: {len(urls_dict)}')
        print(f'Total filesize of download: {sum(urls_dict.values()) / 2 ** 30:.2f} GB')

        self.private_urls = urls_dict

    def estimate_public_dataset_download(self):
        public_urls = {}
        for d in self.uncovered_census_datasets:
            url = self.api_public_dataset_ids_dict[d]['assets'][0]['url']
            filesize = self.api_public_dataset_ids_dict[d]['assets'][0]['filesize']
            public_urls[url] = filesize

        print(f'Total number of h5ad files to download: {len(public_urls)}')
        print(f'Total filesize of download: {sum(public_urls.values()) / 2 ** 30:.2f} GB')

        self.public_urls = public_urls

    def download_datasets(self, which_urls: str, save_dir: str, mutlithreading: bool = False):

        # multithreading approach
        # much quicker but no exception handling, absolute chaos print stream with download percent
        # now reworked to write chunks instead of holding file in memory before saving
        def multi_get_file(url):
            file_from_url = url.split('/')[-1].replace('.h5ad', '')
            file_name = self.api_dataset_version_id_to_dataset_id[file_from_url] + '.h5ad'
            with requests.get(url, stream=True) as res:
                res.raise_for_status()
                filesize = int(res.headers["Content-Length"])
                with open(save_dir + file_name, 'wb') as file:
                    total_bytes_received = 0
                    for chunk in res.iter_content(chunk_size=1024 * 1024):
                        file.write(chunk)
                        total_bytes_received += len(chunk)
                        percent_of_total_upload = float("{:.1f}".format(total_bytes_received / filesize * 100))
                        print(f"{percent_of_total_upload}% downloaded from {file_name}")
            sleep(res.elapsed.total_seconds())


        # based off of scc repo download_assets()
        def single_get_file(url):
            try:
                file_from_url = url.split('/')[-1].replace('.h5ad', '')
                file_name = self.api_dataset_version_id_to_dataset_id[file_from_url] + '.h5ad'
                print(f"\nDownloading {file_name}... ")
                with requests.get(url, stream=True) as res:
                    res.raise_for_status()
                    filesize = int(res.headers["Content-Length"])
                    with open(save_dir + file_name, "wb") as df:
                        total_bytes_received = 0
                        for chunk in res.iter_content(chunk_size=1024 * 1024):
                            df.write(chunk)
                            total_bytes_received += len(chunk)
                            percent_of_total_upload = float("{:.1f}".format(total_bytes_received / filesize * 100))
                            color = "\033[38;5;10m" if percent_of_total_upload == 100 else ""
                            print(f"\033[1m{color}{percent_of_total_upload}% downloaded\033[0m\r", end="")
            except requests.HTTPError as e:
                print(f'Download failed for url {url}: ', e)
                raise e


        if which_urls == 'private':
            url_list = self.private_urls.keys()
        elif which_urls == 'public':
            url_list = self.public_urls.keys()
        elif which_urls == 'both':
            url_list = [k for k in self.public_urls.keys()]
            url_list.extend([k for k in self.private_urls.keys()])
        else:
            print("Please enter 'private', 'public', or 'both' for the which_urls argument")
            return

        start = perf_counter()

        # 5 concurrent downloads, could probably go higher, not sure when CZI 
        # gets upset with spamming or robo traffic
        if mutlithreading:
            with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                executor.map(multi_get_file, url_list)
        else:
            for url in url_list:
                single_get_file(url)

        stop = perf_counter()

        print(f"Took {stop - start} seconds to download all files")

    def get_census_collection_dict(self, census_version: str = 'latest'):
        '''
        param: census_version: str 'latest' | 'stable'
        returns dict with {dataset_id: collection_id} from census for further identification
        '''
        with cellxgene_census.open_soma(census_version=census_version) as census:
            dataset_info = census['census_info']['datasets'].read().concat().to_pandas()
            
        mapping = dataset_info.groupby('dataset_id')['collection_id'].agg('unique')
        df = pd.DataFrame(mapping)
        df['collection_id'] = df['collection_id'].apply(lambda x: x[0])
        df_dict = df.to_dict('index')
        final_dict = {k: v['collection_id'] for k, v in df_dict.items()}
        self.census_collection_dict = final_dict
        
    def get_census_obs_df(self):
        if not self.census_collection_dict:
            self.get_census_collection_dict()

        def get_obs_df(species):
            with cellxgene_census.open_soma(census_version='latest') as census:
                dev_obs = census['census_data'][species] \
                .obs.read(column_names=['dataset_id', 'development_stage', 'development_stage_ontology_term_id', 'donor_id']) \
                .concat().to_pandas()

            dev_obs['collection_id'] = dev_obs['dataset_id'].map(self.census_collection_dict)
            return dev_obs
        
        dfs = [get_obs_df(x) for x in ['mus_musculus', 'homo_sapiens']]
        final_df = pd.concat(dfs)
        self.census_obs_df = final_df

    def create_obs_df_from_files(self, visibility: str, h5ad_location: str):
        h5ads = [f for f in os.listdir(h5ad_location) if f.endswith('.h5ad')]

        def generate_obs_df(h5ad_file):
            dataset_id = h5ad_file.split('.')[0]

            with h5py.File(h5ad_location + h5ad_file) as f:
                obs = read_elem(f['obs'])
                uns = read_elem(f['uns'])

            collection_id = uns['citation'].split('/')[-1]

            donor_df = obs[['donor_id', 'development_stage_ontology_term_id', 'development_stage']] \
                .value_counts() \
                .to_frame() \
                .rename(columns={0: 'counts'}) \
                .reset_index() 
            donor_df['collection_id'] = collection_id
            donor_df['dataset_id'] = dataset_id

            return donor_df

        obs_dfs = [generate_obs_df(file) for file in h5ads]

        final_df = pd.concat(obs_dfs)
        
        if visibility == 'private':
            self.private_all_obs_df = final_df
        else:
            self.public_all_obs_df = final_df

    def dev_query(self, uncovered_terms):
        stages = [f"development_stage == '{v}'" for v in uncovered_terms.values()]
        return ' | '.join(stages)


    def get_all_uncovered_donors(self, visibility):
        # starting with proper uncovered terms
        if visibility == 'census':
            uncovered_terms = self.public_uncovered_terms
        elif visibility == 'private':
            uncovered_terms = self.private_uncovered_terms
        else:
            collection_set = {self.public_dataset_id_to_collection_id[d] for d in self.uncovered_census_datasets}
            uncovered_terms = {k: v for k, v in self.public_uncovered_terms.items() if k in collection_set}

        def get_obs_dev_df(collection_id, uncovered_terms, visibility):
            if visibility == 'census':
                initial_df = self.census_obs_df
            elif visibility == 'private':
                initial_df = self.private_all_obs_df
            else:
                initial_df = self.public_all_obs_df

            filt = initial_df['collection_id'] == collection_id
            df = initial_df[filt][['collection_id', 'donor_id', 'dataset_id', 'development_stage_ontology_term_id', 'development_stage']] \
                .value_counts() \
                .to_frame() \
                .rename(columns={0: 'counts'}) \
                .reset_index() \
                .set_index('donor_id') \
                .query(self.dev_query(uncovered_terms)) \
                .sort_index() \
                .sort_values(by=['collection_id', 'dataset_id', 'donor_id'])
            return df


        dfs = [get_obs_dev_df(k, v, visibility) for k, v in uncovered_terms.items()]
        all_donor_devs = pd.concat(dfs)
        datasets_df = all_donor_devs.groupby(['donor_id', 'development_stage_ontology_term_id', 'development_stage', 'collection_id'])['dataset_id'] \
            .aggregate(list) \
            .reset_index() \
            .sort_values(by=['collection_id', 'donor_id'])

        if visibility == 'census':
            self.census_all_uncovered_donors = datasets_df
        elif visibility == 'private':
            self.private_all_uncovered_donors = datasets_df
        else:
            self.public_all_uncovered_donors = datasets_df

    def validate_google_sheet(self):
        if self.google_sheet.empty:
            print("Please load google sheet by using .get_google_sheet()")

        for idx, row in self.google_sheet.iterrows():
            donor = row['donor_id']
            collection_id = row['collection_id']
            new_dev_stage = row['new_dev_stage']
            visibilty = 'private' if row['tab_name'].startswith('private') else 'public'

            if visibilty == 'private':
                donor_id_dev_terms = self.private_donor_ids_dev_terms
            else:
                donor_id_dev_terms = self.public_donor_ids_dev_terms

            try:
                assert collection_id in donor_id_dev_terms, f'{collection_id = } not in {visibilty} donor_id_dev_terms'
                assert donor in donor_id_dev_terms[collection_id]['donor_id'], f'{donor = } not in {collection_id = }'
                assert new_dev_stage in self.api_current_dev_terms, f'{new_dev_stage} not in new dev terms for {donor = } in {collection_id = }'
            except AssertionError as e:
                print(f"Issue with row {idx}: {e}")
                
    def create_google_sheet_aggregated_dict(self, sheet_tab: str, key_col: str, value_col: str) -> dict:
        '''
        Creates dict from google sheet tab with key column and a list of values aggregated based on key
        '''
        if self.google_sheet.empty:
            print("Please load google sheet by using .get_google_sheet()")
        
        filt = self.google_sheet['tab_name'] == sheet_tab
        all_public_collections_dict = self.google_sheet[filt][[key_col, value_col]] \
           .value_counts() \
           .to_frame() \
           .reset_index() \
           .drop(columns='count') \
           .groupby(key_col) \
           .aggregate(list) \
           .to_dict()[value_col]

        return all_public_collections_dict
