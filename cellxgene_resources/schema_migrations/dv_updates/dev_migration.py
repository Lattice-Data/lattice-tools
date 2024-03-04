import cellxgene_census
import concurrent.futures
import h5py
import json
import os
import pandas as pd
import requests
import sys
from anndata._io.specs import read_elem
from collections import defaultdict
from urllib.parse import quote
from time import perf_counter, sleep
from typing import Collection

# set repo and api key paths, might move this later
CZI_REPO_LOC = os.path.expanduser('~/GitClones/CZI/')
API_KEY_FILE_PATH = os.path.expanduser('~/Documents/keys/cxg-api-key.txt')
sys.path.append(os.path.abspath(CZI_REPO_LOC + 'single-cell-curation/notebooks/curation_api/python/'))

from src.utils.config import set_api_access_config
from src.collection import get_collections, get_collection
from src.dataset import get_dataset,get_datasets

# importing rna assays supported in census, might just hardcode instead
# from the census builder directory so requires python 3.11 and other dependencies
sys.path.append(os.path.abspath(CZI_REPO_LOC + 'cellxgene-census/tools/cellxgene_census_builder/src/'))
from cellxgene_census_builder.build_soma.globals import RNA_SEQ

set_api_access_config(API_KEY_FILE_PATH)

# make sure czi single-cell-curation branch set to lattice/dev-ont-migration to get json files from their repo
MIGRATION_LOCATION = CZI_REPO_LOC + 'single-cell-curation/cellxgene_schema_cli/cellxgene_schema/'

# other constants defined here
AGG_DATASET_DIR = '/Users/brianmott/Documents/Curation/CellxGene/CXG-346/private_uncovered/MigToolDownload/'
GOOGLE_SHEET_ID = '1bELrjC18WH7wVyxlfKPvWjvUKKqy7y4iFav9ddNooAg'
GOOGLE_SHEET_COLS = [
    'donor_id', 
    'development_stage_ontology_term_id', 
    'development_stage', 
    'collection_id', 
    'dataset_ids', 
    'new_dev_stage'
]


class ApiData:
    def __init__(self, from_file: bool = False):
        self.from_file = from_file
        self.private_collections: list[dict] | None = None
        self.private_datasets: list[dict] | None = None
        self.public_datasets: list[dict] | None  = None
        self.current_dev_terms: list | None = None

        self.private_collection_ids: set | None = None

        self.fetch_api_data(self.from_file)

    def fetch_api_data(self, from_file: bool) -> None:
        cxg_api_fxn = self._load_cxg_api_from_file if from_file else self._load_cxg_api_from_web

        for fxn in [self._fetch_dev_terms, cxg_api_fxn]:
            name = fxn.__name__.replace('_', '', 1)
            print(f'Running {name}()...')
            fxn()
            print(f'{name} loaded successfully')

        # only using public collections object for number, don't need to fetch
        public_collections = {d['collection_id'] for d in self.public_datasets}
        print(f"{len(public_collections)} Public Collections")
        print(f"{len(self.private_collection_ids)} Private Collections")
        print(f"{len(self.public_datasets)} Public Datasets")
        print(f"{len(self.private_datasets)} Private Datasets")

    def save_api_jsons(self) -> None:
        jsons = {
            'api_private_datasets.json': self.private_datasets,
            'api_private_collections.json': self.private_collections,
            'api_public_datasets.json': self.public_datasets,
        }

        for file_name, json_attr in jsons.items():
            with open(file_name, 'w', encoding='utf-8') as f:
                json.dump(json_attr, f, ensure_ascii=False, indent=2)

    def _fetch_dev_terms(self) -> None:
        self.current_dev_terms = ['unknown']
        
        for ont in ['hsapdv', 'mmusdv']:
            url = f'https://www.ebi.ac.uk/ols4/api/ontologies/{ont}/terms?obsoletes=false&size=500'
            r = requests.get(url).json()
            terms = [t['obo_id'] for t in r['_embedded']['terms']]
            self.current_dev_terms.extend(terms)

    def _load_cxg_api_from_file(self) -> None:
        try:
            self.public_datasets = json.load(open('api_public_datasets.json'))
            self.private_collections = json.load(open('api_private_collections.json'))
            self.private_datasets = json.load(open('api_private_datasets.json'))
        except FileNotFoundError as e:
            print(f'File not found: {e}')

        self.private_collection_ids = {c['collection_id'] for c in self.private_collections if not c.get('revision_of')}

    def _load_cxg_api_from_web(self) -> None:
            self.public_datasets = get_datasets()
            self.private_collections = get_collections(visibility='PRIVATE')

            self.private_collection_ids = {c['collection_id'] for c in self.private_collections if not c.get('revision_of')}
            # different format for dataset schema, below to get consistent format for downstream functions
            full_private_info = [get_collection(c) for c in self.private_collection_ids]
            private_datasets = []
            for collection in full_private_info:
                collection_id = collection['collection_id']
                datasets = collection['datasets']
                for dataset in datasets:
                    dataset['collection_id'] = collection_id
                    private_datasets.append(dataset)

            self.private_datasets = private_datasets
        

class CensusData:
    def __init__(self, census_version: str = 'latest'):
        self.census_version = census_version

        self._collection_dict: dict[str, str] | None = None
        self.obs_df: pd.DataFrame | None = None

        self.get_data()

    def get_data(self) -> None:
        self._get_collection_dict()
        self._get_obs_df()

    def _get_collection_dict(self) -> None:
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

    def _get_obs_df(self) -> None:
        print('Generating obs_df from census with required columns...')

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


class AggregatedDatasetsDF:
    def __init__(self, visibility: str, h5ad_location: str):
        self.visibility = visibility
        self.h5ad_location = h5ad_location
        self.all_obs_df: pd.DataFrame | None = None

        self.create_obs_df_from_files(self.h5ad_location)

    def create_obs_df_from_files(self, h5ad_location: str) -> None:
        print(f'Generating aggregated obs_df from {self.visibility} datasets...')
        h5ads = [f for f in os.listdir(h5ad_location) if f.endswith('.h5ad')]

        def generate_obs_df(h5ad_file: str) -> pd.DataFrame:
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
        
        print(f'Successfully created {self.visibility} aggregated obs_df')
        self.all_obs_df = final_df


class RepoJSONs:
    def __init__(self, repo_path: str):
        self.repo_path = repo_path
        self.automigrate_terms: dict | None = None
        self.donor_updates: dict | None = None
        self.updated_donor_updates: dict | None = None
        
        self._load_jsons()

    def _load_jsons(self) -> None:
        self.automigrate_terms = json.load(open(self.repo_path + 'automigrate_terms.json'))
        self.donor_updates = json.load(open(self.repo_path + 'donor_updates.json'))
        print(f"Loaded repo JSONs from {'lattice' if self.repo_path == '' else 'single-cell-curation'}")


class GoogleSheet:
    def __init__(self, sheet_id: str, sheet_cols: list[str]):
        self.sheet_id = sheet_id
        self.sheet_cols = sheet_cols
        self.full_sheet: pd.DataFrame | None = None
        self.public_donors: pd.DataFrame | None = None
        self.private_donors: pd.DataFrame | None = None

        self._get_google_sheet()

    def _get_google_sheet(self) -> None:
        dfs = []
        print('Loading data from HsapDv Updates Google sheet...')
        for tab in ['private donors 2024', 'public donors 2024']:
            url = f'https://docs.google.com/spreadsheets/d/{self.sheet_id}/gviz/tq?tqx=out:csv&sheet={quote(tab)}'
            donor_meta = pd.read_csv(url)[self.sheet_cols]
            donor_meta['tab_name'] = tab
            dfs.append(donor_meta)

        final_df = pd.concat(dfs)
        self.full_sheet = final_df
        self.public_donors = final_df[final_df['tab_name'] == 'public donors 2024']
        self.private_donors = final_df[final_df['tab_name'] == 'private donors 2024']
        print('Google sheet data loaded successfully')


class AllData:
    def __init__(self, cxg_from_file: bool = True, load_agg_dfs: bool = False):
        self.cxg_from_file = cxg_from_file
        self.load_agg_dfs = load_agg_dfs

        self.api: ApiData = ApiData(from_file=self.cxg_from_file)
        self.census: CensusData = CensusData()
        self.jsons: RepoJSONs = RepoJSONs('')
        self.google: GoogleSheet = GoogleSheet(GOOGLE_SHEET_ID, GOOGLE_SHEET_COLS)

        if self.load_agg_dfs:
            self.agg_private: AggregatedDatasetsDF = AggregatedDatasetsDF('private', AGG_DATASET_DIR)
            self.agg_public: AggregatedDatasetsDF = AggregatedDatasetsDF('public', AGG_DATASET_DIR + 'public/')
        else:
            self.agg_private, self.agg_public = None, None


#TODO: need to rework this; better to append or create new? how to handle working off of complete build or update?
def update_donor_json(donor_updates: dict, tab_names: list[str] = ['public donors 2024', 'private donors 2024']):
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
    
    return donor_updates


def _make_donor_dev_terms(api_datasets: list[dict]) -> dict:
    donor_dev_terms = {}
    for d in api_datasets:
        c_id = d['collection_id']
        if d['processing_status'] != 'SUCCESS':     # processing private dataset will not have all keys
            print(f"Current processing error with dataset {d['dataset_version_id']} in collection {c_id}")
            continue
        dev_stages = {t['ontology_term_id']:t['label'] for t in d['development_stage']}
        if c_id in donor_dev_terms:
            donor_dev_terms[c_id]['donor_id'].extend(d['donor_id'])
            donor_dev_terms[c_id]['development_stage'].update(dev_stages)
        else:
            donor_dev_terms[c_id] = {
                'donor_id': d['donor_id'],
                'development_stage': dev_stages
            }

    return donor_dev_terms


def make_all_donor_dev_terms(private_datasets: list[dict], public_datasets: list[dict]) -> tuple[dict, dict, dict]:
    all_donor_dev_terms = {}

    private = _make_donor_dev_terms(private_datasets)
    public = _make_donor_dev_terms(public_datasets)
    
    for dd_terms in [private, public]:
        all_donor_dev_terms.update(dd_terms)

    return all_donor_dev_terms, private, public


# TODO: better way to get collection ids from donor_updates_json
def uncovered_terms(
    donor_dev_terms_dict: dict, 
    current_dev_terms: list,
    automigrate_terms: dict,
    donor_updates_json: dict
) -> dict:
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
    collection_dict = donor_dev_terms_dict

    uncovered_collections = {}
    for k,v in collection_dict.items():
        dev_stages = v['development_stage']
        deprecated = [t for t in dev_stages.keys() if t not in current_dev_terms and t.startswith('UBERON:') is False]
        not_migrated = [t for t in deprecated if t not in automigrate_terms.keys()]
        if not_migrated and k not in donor_updates_json:
            print(k)
            for t in not_migrated:
                entry = {t: dev_stages[t]}
                if k in uncovered_collections:
                    uncovered_collections[k].update(entry)
                else:
                    uncovered_collections[k] = entry                
                print(entry)
            print('')
    
    return uncovered_collections


def uncovered_census_terms(uncovered_terms: dict, donor_source: dict) -> list:
    """
    Provide donor source to check against, can be from google sheet or json file/dict
    donor_source dict format: collection_id: list[donor_ids]
    below code will make dict in proper format from google sheet

    public_covered_donors = create_google_sheet_aggregated_dict(
        google_sheet=donor_source,
        sheet_tab='public donors 2024',
        key_col='collection_id',
        value_col='donor_id'
    )

    collections_donors = {}

    for collection, values in data.jsons.donor_updates.items():
        donors = [k for k in values.keys()]
        collections_donors[collection] = donors

    collections_donors

    """
    census_assays = set(RNA_SEQ)
    uncovered_datasets = []

    print("Datasets not fully in Census that contain possible donors with uncovered dev terms:")
    for collection, uncovered in uncovered_terms.items():
        uncovered_set = {k for k in uncovered.keys()}
        json_info = get_collection(collection)
        datasets = json_info['datasets']
        covered_donors = set(donor_source.get(collection, {}))
        for dataset in datasets:
            dataset_id = dataset['dataset_id']
            assays = dataset['assay']
            donors = set(dataset['donor_id'])
            dev_stages = {t['ontology_term_id'] for t in dataset['development_stage']}
            assays_set = {a['ontology_term_id'] for a in assays}
            if (len(assays_set - census_assays) > 0 and 
                len(uncovered_set.intersection(dev_stages)) > 0 and 
                len(donors - covered_donors) > 0
                    ):
                print(f"Unapproved Census assay(s) {assays_set - census_assays} in dataset {dataset_id} from collection: {collection}")
                uncovered_datasets.append(dataset_id)

    return uncovered_datasets


def estimate_private_dataset_download(private_datasets: list[dict], private_uncovered_terms: dict) -> dict:
    urls_dict = {}
    collection_datasets = defaultdict(list)
    for dataset in private_datasets:
        collection_datasets[dataset['collection_id']].append(dataset)

    for cid, terms in private_uncovered_terms.items():
        datasets = collection_datasets[cid]
        uncovered_dev_terms = [k for k in terms.keys()]
        for dataset in datasets:
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

    print(f'Total number of private h5ad files to download: {len(urls_dict)}')
    print(f'Total private dataset download size: {sum(urls_dict.values()) / 2 ** 30:.2f} GB')

    return urls_dict


def estimate_public_dataset_download(public_datasets: list[dict], uncovered_census_datasets: list) -> dict:
    public_urls = {}
    dataset_ids_dict = {d['dataset_id']: d for d in public_datasets}
    for d in uncovered_census_datasets:
        url = dataset_ids_dict[d]['assets'][0]['url']
        filesize = dataset_ids_dict[d]['assets'][0]['filesize']
        public_urls[url] = filesize

    print(f'Total number of public h5ad files to download: {len(public_urls)}')
    print(f'Total public dataset download size: {sum(public_urls.values()) / 2 ** 30:.2f} GB')

    return public_urls


def estimate_full_download(public_urls_dict: dict, private_urls_dict: dict) -> None:
    public_urls_dict.update(private_urls_dict)
    all_urls = public_urls_dict
    print(f'Total number of all h5ad files to download: {len(all_urls)}')
    print(f'Total filesize of all downloads: {sum(all_urls.values()) / 2 ** 30:.2f} GB')


def download_datasets(
    urls_dict: dict[str, int], 
    save_dir: str, 
    private_datasets: list[dict], 
    public_datasets: list[dict], 
    mutlithreading: bool = False
) -> None:
    urls = [k for k in urls_dict.keys()]

    dataset_version_id_to_dataset_id = {}

    for ds in [private_datasets, public_datasets]:
        temp_dict = {d['dataset_version_id']: d['dataset_id'] for d in ds} 
        dataset_version_id_to_dataset_id.update(temp_dict)

    # multithreading approach
    # much quicker but no exception handling, absolute chaos print stream with download percent
    # now reworked to write chunks instead of holding file in memory before saving
    def multi_get_file(url: str) -> None:
        file_from_url = url.split('/')[-1].replace('.h5ad', '')
        file_name = dataset_version_id_to_dataset_id[file_from_url] + '.h5ad'
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
    def single_get_file(url: str) -> None:
        try:
            file_from_url = url.split('/')[-1].replace('.h5ad', '')
            file_name = dataset_version_id_to_dataset_id[file_from_url] + '.h5ad'
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

    start = perf_counter()

    # 5 concurrent downloads, could probably go higher, not sure when CZI 
    # gets upset with spamming or robo traffic
    if mutlithreading:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            executor.map(multi_get_file, urls)
    else:
        for url in urls:
            single_get_file(url)

    stop = perf_counter()

    print(f"Took {stop - start} seconds to download all files")


def _dev_query(uncovered_terms: dict) -> str:
    stages = [f"development_stage == '{v}'" for v in uncovered_terms.values()]
    return ' | '.join(stages)


def get_all_uncovered_donors(uncovered_terms: dict, all_obs_df: pd.DataFrame) -> pd.DataFrame:
    # starting with proper uncovered terms
    def get_obs_dev_df(collection_id: str, uncovered_terms: dict, all_obs_df: pd.DataFrame) -> pd.DataFrame:
        filt = all_obs_df['collection_id'] == collection_id
        df = all_obs_df[filt][['collection_id', 'donor_id', 'dataset_id', 'development_stage_ontology_term_id', 'development_stage']] \
            .value_counts() \
            .to_frame() \
            .rename(columns={0: 'counts'}) \
            .reset_index() \
            .set_index('donor_id') \
            .query(_dev_query(uncovered_terms)) \
            .sort_index() \
            .sort_values(by=['collection_id', 'dataset_id', 'donor_id'])

        return df

    dfs = [get_obs_dev_df(k, v, all_obs_df) for k, v in uncovered_terms.items()]
    all_donor_devs = pd.concat(dfs)
    datasets_df = all_donor_devs.groupby(['donor_id', 'development_stage_ontology_term_id', 'development_stage', 'collection_id'])['dataset_id'] \
        .aggregate(list) \
        .reset_index() \
        .sort_values(by=['collection_id', 'donor_id'])

    return datasets_df


def create_aggregated_dict_from_df(dataframe: pd.DataFrame, key_col: str, value_col: str, value_collection: Collection[set | list]) -> dict:
    '''
    Creates dict from dataframe with key column and a list of values aggregated based on key
    Returns dict[key_col]: list[value_col]
    '''
    agg_dict = dataframe[[key_col, value_col]] \
       .value_counts() \
       .to_frame() \
       .reset_index() \
       .drop(columns='count') \
       .groupby(key_col) \
       .aggregate(value_collection) \
       .to_dict()[value_col]

    return agg_dict


def validate_google_sheet(google_sheet: pd.DataFrame, current_dev_terms: list, donor_dev_terms: dict) -> None:
    for idx, row in google_sheet.iterrows():
        donor = row['donor_id']
        collection_id = row['collection_id']
        new_dev_stage = row['new_dev_stage']
        visibilty = 'private' if row['tab_name'].startswith('private') else 'public'

        try:
            assert collection_id in donor_dev_terms, f'{collection_id = } not in {visibilty} donor_id_dev_terms'
            assert donor in donor_dev_terms[collection_id]['donor_id'], f'{donor = } not in {collection_id = }'
            assert new_dev_stage in current_dev_terms, f'{new_dev_stage} not in new dev terms for {donor = } in {collection_id = }'
        except AssertionError as e:
            print(f"Issue with row {idx}: {e}")
