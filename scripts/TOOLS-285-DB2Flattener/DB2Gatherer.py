from constants import OBJECT_CONFIG, FIELD_TYPES, MAX_URL_LENGTH, BASE_URL_OVERHEAD
import sys
import os
# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

import lattice

class DB2Gatherer:
    def __init__(self, connection):
        self.connection = connection
        self.resolved_objects = {}  # {object_type: {id: object}}

    def extract_uuid_from_id(self, object_id):
        """Extract UUID from @id path like '/tissues/uuid/' -> 'uuid'"""
        if '/' in object_id:
            return object_id.split('/')[-2] if object_id.endswith('/') else object_id.split('/')[-1]
        return object_id
    
    def get_api_type_from_id(self, object_id):
        """Determine API type from @id path"""
        for path_key, config in OBJECT_CONFIG.items():
            if f'/{path_key}/' in object_id:
                return config['api_type']
        return None
    
    def chunk_and_fetch(self, obj_type, object_ids):
        """Fetch objects efficiently with URL chunking if needed"""
        if not object_ids:
            return []
        
        config = None
        for cfg in OBJECT_CONFIG.values():
            if cfg['api_type'] == obj_type:
                config = cfg
                break
        
        if not config:
            print(f"Warning: No config found for {obj_type}")
            return []
        
        print(f"Fetching {len(object_ids)} {obj_type} objects...")
        
        # Remove duplicates
        unique_ids = list(set(object_ids))
        
        # Build filter URL
        filter_url = '&' + '&'.join([f"uuid={oid}" for oid in unique_ids])
        
        # Single request if under limit
        if len(filter_url) <= MAX_URL_LENGTH:
            try:
                results = lattice.get_report(
                    obj_type=obj_type,
                    filter_url=filter_url,
                    field_lst=config['fields'],
                    connection=self.connection
                )
                return results or []
            except Exception as e:
                print(f"Error fetching {obj_type}: {e}")
                return []
        
        # Chunked requests
        print(f"URL too long, chunking...")
        all_results = []
        chunk_size = (MAX_URL_LENGTH - BASE_URL_OVERHEAD) // 50  # Rough estimate
        
        for i in range(0, len(unique_ids), chunk_size):
            chunk_ids = unique_ids[i:i + chunk_size]
            chunk_filter = '&' + '&'.join([f"uuid={oid}" for oid in chunk_ids])
            
            try:
                chunk_results = lattice.get_report(
                    obj_type=obj_type,
                    filter_url=chunk_filter,
                    field_lst=config['fields'],
                    connection=self.connection
                )
                
                if chunk_results:
                    all_results.extend(chunk_results)
                    
            except Exception as e:
                print(f"Chunk error: {e}")
                continue
        
        return all_results
    
    def fetch_controlled_terms(self, reference_ids):
        """Fetch controlled terms by both @id and aliases"""
        if not reference_ids:
            return []
        
        print(f"Fetching controlled terms for {len(reference_ids)} references...")
        
        # Separate @id refs from name refs
        id_refs = [ref for ref in reference_ids if ref.startswith('/controlled_terms/')]
        name_refs = [ref for ref in reference_ids if not ref.startswith('/controlled_terms/')]
        
        all_results = []
        
        # Fetch by @id
        if id_refs:
            id_filter = '&' + '&'.join([f"@id={ref}" for ref in id_refs])
            try:
                id_results = lattice.get_report(
                    obj_type='ControlledTerm',
                    filter_url=id_filter,
                    field_lst=['@id', 'term_name'],
                    connection=self.connection
                )
                if id_results:
                    all_results.extend(id_results)
            except Exception as e:
                print(f"Error fetching controlled terms by @id: {e}")
        
        # Fetch by aliases (name refs)
        if name_refs:
            alias_filter = '&' + '&'.join([f"aliases={ref}" for ref in name_refs])
            try:
                alias_results = lattice.get_report(
                    obj_type='ControlledTerm',
                    filter_url=alias_filter,
                    field_lst=['@id', 'term_name'],
                    connection=self.connection
                )
                if alias_results:
                    all_results.extend(alias_results)
            except Exception as e:
                print(f"Error fetching controlled terms by aliases: {e}")
        
        return all_results
    
    def extract_references_from_field(self, field_value, field_name):
        """Extract reference IDs from a field using FIELD_TYPES for proper handling"""
        if not field_value:
            return []
        
        field_spec = FIELD_TYPES.get(field_name, {'type': 'string'})
        refs = []
        
        if field_spec['type'] == 'array':
            # Handle array fields
            if isinstance(field_value, list):
                for item in field_value:
                    if isinstance(item, dict):
                        # Pre-expanded object like {"@id": "/controlled_terms/xyz/", "term_name": "..."}
                        ref_id = item.get('@id')
                        if ref_id:
                            refs.append(ref_id)
                    elif isinstance(item, str):
                        # String reference like "/controlled_terms/xyz/"
                        refs.append(item)
        else:
            # Handle single value fields (type: 'string')
            if isinstance(field_value, dict):
                # Pre-expanded object
                ref_id = field_value.get('@id')
                if ref_id:
                    refs.append(ref_id)
            elif isinstance(field_value, str):
                # String value or reference
                refs.append(field_value)
        
        return refs
    
    def resolve_references_for_samples(self, all_samples):
        """Collect all references first, then batch fetch by type"""
        all_reference_ids = {}  # {api_type: set(ids)}
        
        # First pass: collect ALL non-controlled-term references from samples
        for sample in all_samples.values():
            sample_api_type = self.get_api_type_from_id(sample['@id'])
            config = None
            for cfg in OBJECT_CONFIG.values():
                if cfg['api_type'] == sample_api_type:
                    config = cfg
                    break
            
            if not config:
                continue
            
            for field_name, ref_types in config.get('references', {}).items():
                field_value = sample.get(field_name)
                refs = self.extract_references_from_field(field_value, field_name)
                
                if isinstance(ref_types, str):
                    ref_types = [ref_types]
                
                for ref in refs:
                    if 'controlled_terms' not in ref_types and ref.startswith('/'):
                        api_type = self.get_api_type_from_id(ref)
                        if api_type:
                            if api_type not in all_reference_ids:
                                all_reference_ids[api_type] = set()
                            all_reference_ids[api_type].add(ref)
        
        # Batch fetch all non-controlled-term references by type
        for api_type, ref_ids in all_reference_ids.items():
            uuids = [self.extract_uuid_from_id(ref_id) for ref_id in ref_ids]
            ref_objects = self.chunk_and_fetch(api_type, uuids)
            
            if api_type not in self.resolved_objects:
                self.resolved_objects[api_type] = {}
            
            for obj in ref_objects:
                self.resolved_objects[api_type][obj['@id']] = obj
        
        # Second pass: collect controlled term references from ALL objects
        controlled_term_refs = set()
        all_objects_to_scan = list(all_samples.values())
        for ref_dict in self.resolved_objects.values():
            all_objects_to_scan.extend(ref_dict.values())
        
        for obj in all_objects_to_scan:
            obj_api_type = self.get_api_type_from_id(obj.get('@id', ''))
            config = None
            for cfg in OBJECT_CONFIG.values():
                if cfg['api_type'] == obj_api_type:
                    config = cfg
                    break
            
            if not config:
                continue
            
            for field_name, ref_types in config.get('references', {}).items():
                if 'controlled_terms' in (ref_types if isinstance(ref_types, list) else [ref_types]):
                    field_value = obj.get(field_name)
                    refs = self.extract_references_from_field(field_value, field_name)
                    
                    for ref in refs:
                        if ref.startswith('/controlled_terms/'):
                            controlled_term_refs.add(ref)
        
        # Batch fetch all controlled terms
        if controlled_term_refs:
            controlled_terms = self.fetch_controlled_terms(list(controlled_term_refs))
            self.resolved_objects['ControlledTerm'] = {}
            for ct in controlled_terms:
                self.resolved_objects['ControlledTerm'][ct['@id']] = ct
    
    def add_references_to_library(self, library_data, samples):
        """Add resolved references to library data based on its samples"""
        added_refs = {
            'donors': set(),
            'treatments': set(),
            'controlled_terms': set(), 
            'genetic_modifications': set(),
            'experimental_conditions': set()
        }
        
        for sample in samples:
            sample_api_type = self.get_api_type_from_id(sample['@id'])
            config = None
            for cfg in OBJECT_CONFIG.values():
                if cfg['api_type'] == sample_api_type:
                    config = cfg
                    break
            
            if not config:
                continue
            
            for field_name, ref_types in config.get('references', {}).items():
                field_value = sample.get(field_name)
                refs = self.extract_references_from_field(field_value, field_name)
                
                for ref in refs:
                    # Find the resolved object
                    resolved_obj = None
                    for api_type, objects in self.resolved_objects.items():
                        if ref in objects:
                            resolved_obj = objects[ref]
                            break
                    
                    if resolved_obj:
                        # Determine which collection to add to
                        collection = None
                        if resolved_obj.get('@id', '').startswith('/human_donors/') or resolved_obj.get('@id', '').startswith('/non_human_donors/'):
                            collection = 'donors'
                        elif resolved_obj.get('@id', '').startswith('/treatments/'):
                            collection = 'treatments'
                        elif resolved_obj.get('@id', '').startswith('/controlled_terms/'):
                            collection = 'controlled_terms'
                        elif resolved_obj.get('@id', '').startswith('/genetic_modifications/'):
                            collection = 'genetic_modifications'
                        elif resolved_obj.get('@id', '').startswith('/experimental_conditions/'):
                            collection = 'experimental_conditions'
                        
                        if collection and ref not in added_refs[collection]:
                            library_data[collection].append(resolved_obj)
                            added_refs[collection].add(ref)
    
    def gather_complete_library_data(self, matrix_file_set_uuid):
        """Main method: gather all data grouped by library"""
        print(f"Gathering library data for MatrixFileSet: {matrix_file_set_uuid}")
        
        # Step 1: Get MatrixFileSet and its raw matrix files
        matrix_file_sets = self.chunk_and_fetch('MatrixFileSet', [matrix_file_set_uuid])
        if not matrix_file_sets:
            print("MatrixFileSet not found")
            return None
        
        matrix_file_set = matrix_file_sets[0]
        
        # Extract raw matrix file UUIDs
        raw_matrix_refs = matrix_file_set.get('raw_matrix_files', [])
        raw_matrix_uuids = []
        for ref in raw_matrix_refs:
            if isinstance(ref, dict):
                ref_id = ref.get('@id', '')
            else:
                ref_id = ref
            if ref_id:
                raw_matrix_uuids.append(self.extract_uuid_from_id(ref_id))
        
        print(f"Found {len(raw_matrix_uuids)} raw matrix files")
        
        # Step 2: Get raw matrix files and their sequence files
        raw_matrix_files = self.chunk_and_fetch('RawMatrixFile', raw_matrix_uuids)
        
        sequence_file_uuids = set()
        for rmf in raw_matrix_files:
            derived_from = rmf.get('derived_from', [])
            for ref in derived_from:
                if isinstance(ref, dict):
                    ref_id = ref.get('@id', '')
                else:
                    ref_id = ref
                if ref_id:
                    sequence_file_uuids.add(self.extract_uuid_from_id(ref_id))
        
        print(f"Found {len(sequence_file_uuids)} sequence files")
        
        # Step 3: Get sequence files and their file sets
        sequence_files = self.chunk_and_fetch('SequenceFile', list(sequence_file_uuids))
        
        file_set_uuids = set()
        for sf in sequence_files:
            file_sets = sf.get('sequence_file_sets', [])
            for ref in file_sets:
                if isinstance(ref, dict):
                    ref_id = ref.get('@id', '')
                else:
                    ref_id = ref
                if ref_id:
                    file_set_uuids.add(self.extract_uuid_from_id(ref_id))
        
        print(f"Found {len(file_set_uuids)} file sets")
        
        # Step 4: Get file sets and their libraries
        file_sets = self.chunk_and_fetch('SequenceFileSet', list(file_set_uuids))
        
        library_uuids = set()
        for fs in file_sets:
            library_ref = fs.get('library', '')
            if library_ref:
                if isinstance(library_ref, dict):
                    lib_id = library_ref.get('@id', '')
                else:
                    lib_id = library_ref
                if lib_id:
                    library_uuids.add(self.extract_uuid_from_id(lib_id))
        
        print(f"Found {len(library_uuids)} libraries")
        
        # Step 5: Get libraries (try both types)
        droplet_libraries = self.chunk_and_fetch('DropletBasedLibrary', list(library_uuids))
        plate_libraries = self.chunk_and_fetch('PlateBasedLibrary', list(library_uuids))
        all_libraries = droplet_libraries + plate_libraries
        
        print(f"Successfully fetched {len(all_libraries)} libraries total")
        
        # Step 6: Get all samples referenced by libraries
        sample_refs = set()
        sample_uuids_by_type = {}  # {api_type: [uuids]}
        
        for library in all_libraries:
            samples = library.get('samples', [])
            for ref in samples:
                if isinstance(ref, dict):
                    ref_id = ref.get('@id', '')
                else:
                    ref_id = ref
                
                if ref_id:
                    sample_refs.add(ref_id)
                    api_type = self.get_api_type_from_id(ref_id)
                    if api_type:
                        if api_type not in sample_uuids_by_type:
                            sample_uuids_by_type[api_type] = []
                        sample_uuids_by_type[api_type].append(self.extract_uuid_from_id(ref_id))
        
        # Fetch all sample types
        all_samples = {}
        for api_type, uuids in sample_uuids_by_type.items():
            samples = self.chunk_and_fetch(api_type, uuids)
            for sample in samples:
                all_samples[sample['@id']] = sample
        
        print(f"Found {len(all_samples)} samples total")
        
        # Step 7: Resolve all references from samples
        self.resolve_references_for_samples(all_samples)
        
        # Step 8: Structure data by library
        libraries_data = {}

        for library in all_libraries:
            lib_uuid = library['uuid']
            libraries_data[lib_uuid] = {
                'library': library,
                'samples': [],
                'donors': [],
                'treatments': [],
                'controlled_terms': [],
                'genetic_modifications': [],
                'experimental_conditions': []
            }
            
            # Add samples for this library
            library_samples = library.get('samples', [])
            library_sample_objects = []
            for ref in library_samples:
                if isinstance(ref, dict):
                    ref_id = ref.get('@id', '')
                else:
                    ref_id = ref
                
                if ref_id in all_samples:
                    sample_obj = all_samples[ref_id]
                    libraries_data[lib_uuid]['samples'].append(sample_obj)
                    library_sample_objects.append(sample_obj)
            
            # Add all resolved references for this library's samples
            self.add_references_to_library(libraries_data[lib_uuid], library_sample_objects)
        
        print(f"Structured data for {len(libraries_data)} libraries")
        
        return {
            'matrix_file_set': matrix_file_set,
            'libraries': libraries_data,
            'resolved_objects': self.resolved_objects
        }
    