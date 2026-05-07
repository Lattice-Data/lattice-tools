import lattice, os, json
import pandas as pd
from dataclasses import dataclass, asdict, field, InitVar
import argparse
import sys
from typing import ClassVar, Optional

"""
Initial Rough Draft Script for GEO Flattening

Takes as input MatrixFileSet UUID
"""

# Constants
os.environ['DEMO_KEY'] = ''
os.environ['DEMO_SECRET'] = ''
os.environ['DEMO_SERVER'] = 'https://lattice-api-dev.demo.lattice-data.org'
conn = lattice.Connection('demo')

def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--UUID',
        '-U',
        dest='uuid',
        help='UUID of Matrix File Set to flatten'
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return args


def determine_sample_obj_type(sample_filter):
    """
    Determine the correct object type based on the sample filter
    """
    if 'tissue' in sample_filter.lower():
        return 'Tissue'
    elif 'cell_line' in sample_filter.lower() or 'cellline' in sample_filter.lower():
        return 'CellLine'
    elif 'organoid' in sample_filter.lower():
        return 'Organoid'
    elif 'primary_cell' in sample_filter.lower() or 'primarycell' in sample_filter.lower():
        return 'PrimaryCellCulture'
    else:
        # Default or try to detect from other patterns
        return 'Tissue'  # or raise an exception

def determine_donor_obj_type(donor_filter):
    """
    Determine the correct donor object type
    """
    if 'human' in donor_filter.lower():
        return 'HumanDonor'
    elif 'non_human' in donor_filter.lower() or 'nonhuman' in donor_filter.lower():
        return 'NonHumanDonor'
    else:
        return 'HumanDonor'


def simple_chunk_filter(filter_url, max_length=4000):
    """
    Simple chunking: split filter_url into chunks that stay under max_length
    Handles both complete URLs and query parameter strings
    """
    if len(filter_url) <= max_length:
        return [filter_url]
    
    # Handle case where filter_url is just query parameters (starts with &)
    if filter_url.startswith('&'):
        # Remove the leading & and split by &
        params = filter_url[1:].split('&')
    elif '?' in filter_url:
        # Normal URL with base and query
        base_url, query = filter_url.split('?', 1)
        params = query.split('&')
    else:
        # Just parameters without leading &
        params = filter_url.split('&')
        base_url = ''
    
    # Separate UUID params from other params
    uuid_params = []
    other_params = []
    
    for p in params:
        if p.startswith('uuid='):
            # Fix the UUID format - extract just the UUID from the path
            if '/sequence_files/' in p and p.endswith('/'):
                # Convert uuid=/sequence_files/XXX/ to uuid=XXX
                uuid_part = p.split('/sequence_files/')[1].rstrip('/')
                uuid_params.append(f'uuid={uuid_part}')
            else:
                # Already in correct format
                uuid_params.append(p)
        else:
            other_params.append(p)
    
    # If no UUID params to chunk, return as-is
    if not uuid_params:
        return [filter_url]
    
    # Calculate base size for other parameters
    other_query = '&'.join(other_params)
    base_length = len(other_query)
    if other_query:
        base_length += 1  # for the '&' before first UUID
    
    chunks = []
    current_chunk = []
    current_length = base_length
    
    for uuid_param in uuid_params:
        param_length = len(uuid_param) + 1  # +1 for '&'
        
        if current_length + param_length > max_length:
            # Save current chunk and start new one
            if current_chunk:
                if other_query:
                    chunk_query = '&' + other_query + '&' + '&'.join(current_chunk)
                else:
                    chunk_query = '&' + '&'.join(current_chunk)
                chunks.append(chunk_query)
            
            # Start new chunk
            current_chunk = [uuid_param]
            current_length = base_length + len(uuid_param) + 1
        else:
            current_chunk.append(uuid_param)
            current_length += param_length
    
    # Add final chunk
    if current_chunk:
        if other_query:
            chunk_query = '&' + other_query + '&' + '&'.join(current_chunk)
        else:
            chunk_query = '&' + '&'.join(current_chunk)
        chunks.append(chunk_query)
    
    return chunks


def safe_get_report(obj_type, filter_url, field_lst, connection):
    """Simple safe wrapper - check length upfront and chunk if needed"""
    
    # Ensure filter_url starts with & if it doesn't already
    if not filter_url.startswith('&') and not filter_url.startswith('?'):
        filter_url = '&' + filter_url
    
    if len(filter_url) <= 4000:
        # Under limit, make single request
        return lattice.get_report(
            obj_type=obj_type,
            filter_url=filter_url,
            field_lst=field_lst,
            connection=connection
        )
    
    # Over limit, chunk it
    chunks = simple_chunk_filter(filter_url, 4000)
    all_results = []
    
    for chunk_url in chunks:
        chunk_results = lattice.get_report(
            obj_type=obj_type,
            filter_url=chunk_url,
            field_lst=field_lst,
            connection=connection
        )
        
        if chunk_results:
            all_results.extend(chunk_results)
    
    return all_results


@dataclass
class SequenceFile:
    """
    Dataclass to hold parsed information for SequenceFile
    """
    uuid: str
    file_type: str = field(init=False)
    library_id: str = field(init=False)
    sequence_file_sets: list[str] = field(init=False)  # List of sequence file set IDs
    
    _instances: ClassVar[dict[str, "SequenceFile"]] = {}

    def __post_init__(self):
        # This will be called when created without bulk data
        seq_file_json = lattice.get_object(self.uuid, conn, frame='object')
        self._populate_from_json(seq_file_json)

    def _populate_from_json(self, seq_file_json: dict):
        """Populate fields from JSON data."""
        self.file_type = seq_file_json['file_format']
        full_library_id = seq_file_json['aliases'][0].split(":")[1]
    
        # Extract base library ID for linking to Library objects
        self.library_id = full_library_id.split('_S1_')[0] if '_S1_' in full_library_id else full_library_id
        self.sequence_file_sets = seq_file_json.get('sequence_file_sets', [])
    
        # Link to library if it exists
        library = Library.get_by_library_id(self.library_id)
        if library is not None:
            if not hasattr(library, 'sequence_files'):
                library.sequence_files = []
            if self not in library.sequence_files:
                library.sequence_files.append(self)

    @classmethod
    def get_or_create(cls, uuid: str) -> "SequenceFile":
        """Return existing SequenceFile if UUID exists, otherwise create and register a new one."""
        if uuid not in cls._instances:
            cls._instances[uuid] = cls(uuid=uuid)
        return cls._instances[uuid]

    @classmethod
    def get_or_create_with_data(cls, seq_file_data: dict) -> "SequenceFile":
        """Create or return existing SequenceFile using prefetched data."""
        uuid = seq_file_data['uuid']
        
        if uuid not in cls._instances:
            # Create instance without calling __post_init__
            instance = cls.__new__(cls)
            instance.uuid = uuid
            instance._populate_from_json(seq_file_data)
            cls._instances[uuid] = instance
        
        return cls._instances[uuid]


@dataclass
class Donor:
    uuid: str
    ethnicity: str = field(init=False)
    sex: str = field(init=False)
    json_object: InitVar[dict]

    _instances: ClassVar[dict[str, "Donor"]] = {}

    def __post_init__(self, json_object: dict):
        self.sex = json_object.get('sex', 'n/a')
        self.ethnicity = json_object.get('ethnicity_term', 'n/a')

    @classmethod
    def get_or_create(cls, json_object: dict) -> "Donor":
        uuid = json_object['uuid']
        if uuid not in cls._instances:
            cls._instances[uuid] = cls(uuid=uuid, json_object=json_object)
        return cls._instances[uuid]

@dataclass
class Biosample:
    uuid: str
    collection_date: str = field(init=False)
    tissue: list[str] = field(init=False)
    donors: list[Donor] = field(init=False)
    age_units: str = field(init=False)
    lower_bound_age: str = field(init=False)
    upper_bound_age: str = field(init=False)
    genetic_modifications: list[str] = field(init=False)
    treatments: list[str] = field(init=False)
    json_object: InitVar[dict]
    donor_map: InitVar[dict]
    genetic_modification_map: InitVar[dict]
    treatment_map: InitVar[dict]

    _instances: ClassVar[dict[str, "Biosample"]] = {}

    def __post_init__(self, json_object: dict, donor_map: dict, genetic_modification_map: dict, treatment_map: dict):
        self.collection_date = json_object.get('date_obtained')
        self.tissue = [t['term_name'] for t in json_object['sample_terms']]
        self.donors = [Donor.get_or_create(donor_map[d['@id']]) for d in json_object['donors']]
        self.age_units = json_object.get('age_units') or 'n/a'
        self.lower_bound_age = json_object.get('lower_bound_age') or 'n/a'
        self.upper_bound_age = json_object.get('upper_bound_age') or 'n/a'
        # Genetic Modification information is stored within modality property of the genetic_modification object
        genetic_mod_id = json_object.get('genetic_modification')
    
        if genetic_mod_id:
            modality = genetic_modification_map.get(genetic_mod_id, {}).get('modality', 'n/a')
            self.genetic_modifications = [modality]
        else:
            self.genetic_modifications = ['n/a']
        
        # Treatment handling, different than Genetic Modification because it's an array field, but similar in that it's objects
        treatment_ids = json_object.get('treatments', [])
        if treatment_ids:
            self.treatments = [
                treatment_map.get(treatment_id, {}).get('ontological_term', 'n/a') 
                for treatment_id in treatment_ids
            ]
        else:
            self.treatments = ['n/a']

    @classmethod
    def get_or_create(cls, json_object: dict, donor_map: dict, 
                      genetic_modification_map: dict = None, treatment_map: dict = None) -> "Biosample":
        uuid = json_object['uuid']
        if uuid not in cls._instances:
            cls._instances[uuid] = cls(
                uuid=uuid, 
                json_object=json_object, 
                donor_map=donor_map, 
                genetic_modification_map=genetic_modification_map or {},
                treatment_map=treatment_map or {}
            )
        return cls._instances[uuid]


@dataclass
class Library:
    uuid: str
    library_alias: str = field(init=False)
    library_id: str = field(init=False)
    samples: list[Biosample] = field(init=False)
    library_strategy: str = field(init=False)
    design_description: str = field(init=False)
    sequence_files: list[SequenceFile] = field(init=False) 

    _instances: ClassVar[dict[str, "Library"]] = {}
    _by_library_id: ClassVar[dict[str, "Library"]] = {}

    def __post_init__(self):
        # Fetch library
        library_json = lattice.get_object(self.uuid, conn, frame='object')
        
        # Use the @id as the library_id for mapping
        self.library_id = library_json['@id']
        
        # Extract alias from aliases field
        aliases = library_json.get('aliases', [])
        if aliases:
            # Get the part after the colon
            self.library_alias = aliases[0].split(':')[-1] if ':' in aliases[0] else aliases[0]
        else:
            self.library_alias = self.library_id  # Fallback
        
        self.library_strategy = library_json.get('chemistry_version')
        self.design_description = library_json.get('description')
    
        # Bulk fetch all samples
        sample_filter = ''.join([f"&@id={s}" for s in library_json['samples']])
        sample_obj_type = determine_sample_obj_type(sample_filter)
        sample_results = safe_get_report(
            obj_type=sample_obj_type,
            filter_url=sample_filter,
            field_lst=['uuid', '@id', 'date_obtained', 'sample_terms', 'donors', 'age_units', 
                       'lower_bound_age', 'upper_bound_age', 'genetic_modification', 'treatments'],
            connection=conn
        )
        # Handle getting genetic modifications as strings
        genetic_mod_ids = list({
            s.get('genetic_modification') 
            for s in sample_results 
            if s.get('genetic_modification')
        })

        # Handle getting Treatments as arrays
        treatment_ids = list({
            treatment_id
            for s in sample_results 
            for treatment_id in s.get('treatments', [])
        })

        # Fetch Genetic Modifications
        genetic_modification_map = {}
        if genetic_mod_ids:
            genetic_mod_filter = ''.join([f"&@id={gm_id}" for gm_id in genetic_mod_ids])
            genetic_mod_results = safe_get_report(
                obj_type='GeneticModification',
                filter_url=genetic_mod_filter,
                field_lst=['@id', 'modality'],
                connection=conn
            )
            genetic_modification_map = {gm['@id']: gm for gm in genetic_mod_results}

        # Fetch treatments
        treatment_map = {}
        if treatment_ids:
            treatment_filter = ''.join([f"&@id={treatment_id}" for treatment_id in treatment_ids])
            treatment_results = safe_get_report(
                obj_type='Treatment',
                filter_url=treatment_filter,
                field_lst=['@id', 'ontological_term'],
                connection=conn
            )
            treatment_map = {t['@id']: t for t in treatment_results}

        # Bulk fetch all donors
        donor_ids = list({d['@id'] for s in sample_results for d in s['donors']})
        donor_filter = ''.join([f"&@id={d}" for d in donor_ids])
        donor_obj_type = determine_donor_obj_type(donor_filter)
        donor_results = safe_get_report(
            obj_type=donor_obj_type,
            filter_url=donor_filter,
            field_lst=['uuid', '@id', 'sex', 'ethnicity'],
            connection=conn
        )
        
        # Handle case where ethnicity might be None/empty
        ethnicity_ids = [d['ethnicity'] for d in donor_results if d.get('ethnicity')]
        ethnicity_map = {}
        
        if ethnicity_ids:
            ethnicity_filter = ''.join([f"&@id={e}" for e in ethnicity_ids])
            ethnicity_results = safe_get_report(
                obj_type='ControlledTerm',
                filter_url=ethnicity_filter,
                field_lst=['@id', 'term_name'],
                connection=conn
            )
            ethnicity_map = {e['@id']: e['term_name'] for e in ethnicity_results}
    
        # Donor map creation
        donor_map = {
            d['@id']: {
                **d, 
                'ethnicity_term': ethnicity_map.get(d.get('ethnicity'), 'n/a')
            }
            for d in donor_results
        }

        # Construct nested dataclass objects from prefetched data
        self.samples = [Biosample.get_or_create(s, donor_map, genetic_modification_map, treatment_map) for s in sample_results]

        Library._by_library_id[self.library_id] = self

    @classmethod
    def get_or_create(cls, uuid: str) -> "Library":
        if uuid not in cls._instances:
            cls._instances[uuid] = cls(uuid=uuid)
        return cls._instances[uuid]

    @classmethod
    def get_by_library_id(cls, library_id: str) -> Optional["Library"]:
        """Get library instance by library ID."""
        return cls._by_library_id.get(library_id)

    @property
    def tissues(self) -> list[str]:
        return sorted(set(t for b in self.samples for t in b.tissue))

    @property
    def collection_dates(self) -> list[str]:
        dates = sorted(set(b.collection_date for b in self.samples if b.collection_date is not None))
        return dates if dates else None

    @property
    def donors(self) -> list[Donor]:
        seen = set()
        return [d for b in self.samples for d in b.donors if not (d.uuid in seen or seen.add(d.uuid))]

    @property
    def donor_sexes(self) -> list[str]:
        return sorted(set(d.sex for b in self.samples for d in b.donors))

    @property
    def donor_ethnicities(self) -> list[str]:
        return sorted(set(d.ethnicity for b in self.samples for d in b.donors))

    @property
    def age_units_list(self) -> list[str]:
        """Get all unique age units from samples."""
        units = [b.age_units for b in self.samples if b.age_units != 'n/a']
        return sorted(set(units)) if units else ['n/a']
    
    @property
    def lower_bound_ages(self) -> list[str]:
        """Get all unique lower bound ages from samples."""
        ages = [b.lower_bound_age for b in self.samples if b.lower_bound_age != 'n/a']
        return sorted(set(ages)) if ages else ['n/a']
    
    @property
    def upper_bound_ages(self) -> list[str]:
        """Get all unique upper bound ages from samples."""
        ages = [b.upper_bound_age for b in self.samples if b.upper_bound_age != 'n/a']
        return sorted(set(ages)) if ages else ['n/a']

    @property
    def genetic_modifications(self) -> list[str]:
        """Get all unique genetic modification modalities from samples."""
        mods = [mod for b in self.samples for mod in b.genetic_modifications if mod != 'n/a']
        return sorted(set(mods)) if mods else ['n/a']

    @property
    def treatments(self) -> list[str]:
        """Get all unique treatment ontological terms from samples."""
        treatments = [treatment for b in self.samples for treatment in b.treatments if treatment != 'n/a']
        return sorted(set(treatments)) if treatments else ['n/a']


@dataclass
class RawMatrixFile:
    """
    Dataclass to hold parsed information for RawMatrixFile
    """
    uuid: str
    s3_uri: str = field(init=False)
    sequence_file_sets: dict[str, list[SequenceFile]] = field(init=False)
    sequence_fileset_libraries: dict[str, Library] = field(init=False)
    sequence_fileset_metadata: dict[str, dict] = field(init=False)
    
    _instances: ClassVar[dict[str, "RawMatrixFile"]] = {}

    def __post_init__(self):
        # Fetch raw matrix file data
        matrix_json = lattice.get_object(self.uuid, conn, frame='object')
        self.s3_uri = matrix_json['s3_uri']
        
        # Bulk fetch all sequence files
        sequence_file_uuids = matrix_json['derived_from']
        sequence_file_filter = ''.join([f"&uuid={file_uuid}" for file_uuid in sequence_file_uuids])
        
        sequence_file_results = safe_get_report(
            obj_type='SequenceFile',
            filter_url=sequence_file_filter,
            field_lst=['uuid', 'file_format', 'aliases', 'sequence_file_sets'],
            connection=conn
        )
        
        # Create SequenceFile objects with prefetched data
        sequence_files = [
            SequenceFile.get_or_create_with_data(sf_data) 
            for sf_data in sequence_file_results
        ]
        
        # Group sequence files by their sequence file sets
        self.sequence_file_sets = self._group_by_sequence_file_sets(sequence_files)
        
        # Create mapping of sequence file sets to libraries
        self.sequence_fileset_libraries = self._map_filesets_to_libraries()

    def _group_by_sequence_file_sets(self, sequence_files: list[SequenceFile]) -> dict[str, list[SequenceFile]]:
        """Group sequence files by their sequence file sets."""
        grouped = {}
        
        for seq_file in sequence_files:
            for fileset_id in seq_file.sequence_file_sets:
                if fileset_id not in grouped:
                    grouped[fileset_id] = []
                grouped[fileset_id].append(seq_file)
        
        return grouped

    def _map_filesets_to_libraries(self) -> dict[str, Library]:
        """Create mapping of sequence file sets to their libraries."""
        fileset_to_library = {}
        
        if not self.sequence_file_sets:
            return fileset_to_library
        
        print(f"Fetching {len(self.sequence_file_sets)} sequence file sets to get library IDs...")
        
        # Bulk fetch sequence file set data
        fileset_filter = ''.join([f"&@id={fileset_id}" for fileset_id in self.sequence_file_sets.keys()])
        
        try:
            fileset_results = safe_get_report(
                obj_type='SequenceFileSet',
                filter_url=fileset_filter,
                field_lst=['@id', 'library', 'run_cardinality', 'sequencing_platform'],
                connection=conn
            )
            
            print(f"Got {len(fileset_results)} sequence file set results")

            self.sequence_fileset_metadata = {
            fs['@id']: {
                'run_cardinality': fs.get('run_cardinality', 'n/a'),
                'sequencing_platform': fs.get('sequencing_platform', 'n/a')
            }
            for fs in fileset_results
        }
            
            # Extract library IDs and group by type
            droplet_library_ids = set()
            plate_library_ids = set()
            fileset_to_library_id = {}
            
            for fileset_data in fileset_results:
                fileset_id = fileset_data['@id']
                library_info = fileset_data.get('library')
                if library_info and '@id' in library_info:
                    library_id = library_info['@id']
                    fileset_to_library_id[fileset_id] = library_id
                    
                    # Determine library type from the path
                    if '/droplet_based_libraries/' in library_id:
                        droplet_library_ids.add(library_id)
                    elif '/plate_based_libraries/' in library_id:
                        plate_library_ids.add(library_id)
            
            total_libraries = len(droplet_library_ids) + len(plate_library_ids)
            print(f"Found {total_libraries} unique library IDs: {len(droplet_library_ids)} droplet-based, {len(plate_library_ids)} plate-based")
            
            # Fetch libraries by type
            all_libraries = []
            
            # Fetch droplet-based libraries
            if droplet_library_ids:
                droplet_filter = ''.join([f"&@id={lib_id}" for lib_id in droplet_library_ids])
                try:
                    droplet_results = safe_get_report(
                        obj_type='DropletBasedLibrary',
                        filter_url=droplet_filter,
                        field_lst=['uuid', '@id', 'aliases', 'sample'],
                        connection=conn
                    )
                    if droplet_results:
                        print(f"Found {len(droplet_results)} droplet_based_library results")
                        all_libraries.extend(droplet_results)
                except Exception as e:
                    print(f"Error fetching droplet_based_library: {e}")
            
            # Fetch plate-based libraries
            if plate_library_ids:
                plate_filter = ''.join([f"&@id={lib_id}" for lib_id in plate_library_ids])
                try:
                    plate_results = safe_get_report(
                        obj_type='PlateBasedLibrary',
                        filter_url=plate_filter,
                        field_lst=['uuid', '@id', 'aliases', 'sample'],
                        connection=conn
                    )
                    if plate_results:
                        print(f"Found {len(plate_results)} plate_based_library results")
                        all_libraries.extend(plate_results)
                except Exception as e:
                    print(f"Error fetching plate_based_library: {e}")
            
            # Create Library objects and map by @id
            library_by_id = {}
            for lib_data in all_libraries:
                library = Library.get_or_create(lib_data['uuid'])
                library_by_id[lib_data['@id']] = library
            
            # Map filesets to libraries
            for fileset_id, library_id in fileset_to_library_id.items():
                if library_id in library_by_id:
                    fileset_to_library[fileset_id] = library_by_id[library_id]
                else:
                    print(f"Warning: Could not find library object for library_id {library_id}")
            
        except Exception as e:
            print(f"Error fetching sequence file sets: {e}")
        
        return fileset_to_library


    @property
    def file_types(self) -> list[str]:
        """Get all unique file types across all sequence files."""
        return sorted(set(
            sf.file_type 
            for fileset in self.sequence_file_sets.values() 
            for sf in fileset
        ))

    @property
    def library_ids(self) -> list[str]:
        """Get all unique library IDs across all sequence files."""
        return sorted(set(
            sf.library_id 
            for fileset in self.sequence_file_sets.values() 
            for sf in fileset
        ))

    @property
    def libraries(self) -> list[Library]:
        """Get all unique libraries across all sequence file sets."""
        return list(self.sequence_fileset_libraries.values())
    
    @classmethod
    def get_or_create(cls, uuid: str) -> "RawMatrixFile":
        """Return existing RawMatrixFile if UUID exists, otherwise create and register a new one."""
        if uuid not in cls._instances:
            cls._instances[uuid] = cls(uuid=uuid)
        return cls._instances[uuid]

def create_dataframe(raw_matrix_files: list[RawMatrixFile], genome_annotation: str = 'n/a', genome_assembly: str = 'n/a') -> pd.DataFrame:
    """
    Create a detailed DataFrame indexed by library alias containing aggregated library information
    from RawMatrixFile objects with comma-separated arrays for multiple values.
    """
    library_data = []
    
    for rmf in raw_matrix_files:
        for fileset_id, library in rmf.sequence_fileset_libraries.items():
            
            # Get sequence files for this fileset
            seq_files = rmf.sequence_file_sets.get(fileset_id, [])

            # Get SequenceFileSet metadata
            fileset_metadata = rmf.sequence_fileset_metadata.get(fileset_id, {})
            
            # Basic library information
            library_info = {
                'library_alias': library.library_alias,
                'library_strategy': library.library_strategy,
                'design_description': library.design_description,
                'genetic_modifications': ','.join(library.genetic_modifications),
                'treatments': ','.join(library.treatments),

                # MatrixFileSet-level information
                'genome_annotation': genome_annotation,
                'genome_assembly': genome_assembly,
                
                # SequenceFileSet-level information
                'run_cardinality': fileset_metadata.get('run_cardinality', 'n/a'),
                'sequencing_platform': fileset_metadata.get('sequencing_platform', 'n/a'),
                
                # Raw matrix file information
                'raw_matrix_s3_uri': rmf.s3_uri,
                
                # Sequence file information - comma-separated arrays
                'sequence_file_types': ','.join(sorted(set(sf.file_type for sf in seq_files))),
                
                # Sample information - comma-separated arrays
                'tissues': ','.join(library.tissues) if library.tissues else 'n/a',
                'collection_dates': ','.join(library.collection_dates) if library.collection_dates else 'n/a',
                'age_units': ','.join(library.age_units_list) if library.age_units_list else 'n/a',
                'lower_bound_age': ','.join(library.lower_bound_ages) if library.lower_bound_ages else 'n/a',
                'upper_bound_age': ','.join(library.upper_bound_ages) if library.upper_bound_ages else 'n/a',
                
                # Donor information - comma-separated arrays
                'donor_sexes': ','.join(library.donor_sexes) if library.donor_sexes else 'n/a',
                'donor_ethnicities': ','.join(library.donor_ethnicities) if library.donor_ethnicities else 'n/a',
            }
            
            library_data.append(library_info)
    
    # Create DataFrame
    if not library_data:
        return pd.DataFrame()
    
    df = pd.DataFrame(library_data)
    
    # Remove duplicate rows if any library appears multiple times
    df = df.drop_duplicates(subset=['library_alias'])
    
    # Set library_alias as index
    df.set_index('library_alias', inplace=True)
    
    return df



if __name__ == "__main__":
    args = getArgs()
    try:
        print(f"Processing MatrixFileSet UUID: {args.uuid}")
        
        # Get MatrixFileSet data
        MatrixFileSet_json = lattice.get_object(args.uuid, conn, frame='object')
        if not MatrixFileSet_json:
            print(f"Error: Could not find MatrixFileSet with UUID {args.uuid}")
            sys.exit(1)

        # Extract genome information from MatrixFileSet
        genome_annotation = MatrixFileSet_json.get('genome_annotation', 'n/a')
        genome_assembly = MatrixFileSet_json.get('genome_assembly', 'n/a')
            
        if not MatrixFileSet_json.get('raw_matrix_files'):
            print("Warning: No raw matrix files found in this MatrixFileSet")
            sys.exit(0)
        
        print(f"Found {len(MatrixFileSet_json['raw_matrix_files'])} raw matrix files")
        
        # Get RawMatrixFile data
        RawMatrix_filter = ''.join([f"&@id={s}" for s in MatrixFileSet_json['raw_matrix_files']])
        raw_matrix_results = safe_get_report(
            obj_type='RawMatrixFile',
            filter_url=RawMatrix_filter,
            field_lst=['uuid', '@id', 's3_uri', 'derived_from'],
            connection=conn
        )
        
        print(f"Processing {len(raw_matrix_results)} raw matrix files...")
        
        # Create RawMatrixFile objects (this will fetch all nested data)
        raw_matrix_files = [RawMatrixFile.get_or_create(rmf['uuid']) for rmf in raw_matrix_results]
        
        print("Creating DataFrame...")
        
        # Create DataFrame
        df = create_dataframe(raw_matrix_files, genome_annotation, genome_assembly)
        
        # Display results
        print(f"\nDataFrame created!")
        print(f"Shape: {df.shape}")
        print(f"Libraries: {len(df)}")
        print(f"Columns: {len(df.columns)}")
        
        # Show sample of the data
        print("\nFirst few rows:")
        print(df.head())
        
        print("\nColumn names:")
        for i, col in enumerate(df.columns):
            print(f"  {i+1:2d}. {col}")
        
        # Save to CSV
        output_filename = f'matrix_fileset_{args.uuid}_detailed.csv'
        df.to_csv(output_filename)
        print(f"\nDataFrame saved to: {output_filename}")
        
        # Optional: Save to Excel for easier viewing
        excel_filename = f'matrix_fileset_{args.uuid}_detailed.xlsx'
        df.to_excel(excel_filename)
        print(f"Also saved as Excel file: {excel_filename}")
        
    except Exception as e:
        import traceback
        print(f"Error: {e}")
        print(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)
    
    
    # Make report request for MatrixFileSet info
    # Get list of raw matrix files from json
    # Create raw matrix file objects, calculate all sequence files needed, make report request for all sequence files
    # Get sequence file sets, group sequence file objects by sequence file set, get report for sequence file sets
    # get libraries from sequence file set report, make library objects

    # RawMatrixFile dataclass will contain list of sequence file sets containing sequence file dataclass objects, list of libraries, map for library->sequencefileset
    # Library dataclass will then contain donor, sample info etc.
    
    