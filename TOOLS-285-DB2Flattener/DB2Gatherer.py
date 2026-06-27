from constants import OBJECT_CONFIG, FIELD_TYPES, MAX_URL_LENGTH, BASE_URL_OVERHEAD
import sys
import os
import DB2lattice

# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)


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
                results = DB2lattice.get_report(
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
                chunk_results = DB2lattice.get_report(
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
        
    def extract_controlled_term_id(self, controlled_term_ref):
        """Extract the semantic term ID from controlled term @id path"""
        if '/controlled_terms/' in controlled_term_ref:
            # Extract everything after '/controlled_terms/' and before trailing '/'
            term_id = controlled_term_ref.split('/controlled_terms/')[-1]
            if term_id.endswith('/'):
                term_id = term_id[:-1]  # Remove trailing slash
            return term_id
        return controlled_term_ref
    
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
        """Collect all references first, then batch fetch by type (excluding controlled terms)"""
        all_reference_ids = {}  # {api_type: set(ids)}
        controlled_term_values = {}  # {ref_path: extracted_term_id}
        
        # First pass: collect ALL references from samples
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
                    if ref.startswith('/'):
                        if 'controlled_terms' in ref_types and '/controlled_terms/' in ref:
                            # Extract term ID directly for controlled terms
                            term_id = self.extract_controlled_term_id(ref)
                            controlled_term_values[ref] = term_id
                        else:
                            # Regular UUID-based reference
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
        
        # Store controlled term values (no API calls needed)
        self.resolved_objects['ControlledTerm'] = controlled_term_values
        
        # Second pass: collect controlled term references from non-sample objects
        for ref_dict in self.resolved_objects.values():
            if ref_dict == controlled_term_values:  # Skip the controlled terms dict
                continue
                
            for obj in ref_dict.values():
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
                                term_id = self.extract_controlled_term_id(ref)
                                controlled_term_values[ref] = term_id
    
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
                    # Handle controlled terms specially
                    if '/controlled_terms/' in ref and 'controlled_terms' in \
                    (ref_types if isinstance(ref_types, list) else [ref_types]):
                        if ref not in added_refs['controlled_terms']:
                            # Add the extracted term ID instead of making API calls
                            term_id = self.resolved_objects['ControlledTerm'].get(ref)
                            if term_id:
                                library_data['controlled_terms'].append({
                                    '@id': ref,
                                    'term_id': term_id
                                })
                                added_refs['controlled_terms'].add(ref)
                    else:
                        # Find the resolved object for other reference types
                        resolved_obj = None
                        for api_type, objects in self.resolved_objects.items():
                            if api_type != 'ControlledTerm' and ref in objects:
                                resolved_obj = objects[ref]
                                break
                        
                        if resolved_obj:
                            # Determine which collection to add to
                            collection = None
                            if resolved_obj.get('@id', '').startswith('/human_donors/') \
                            or resolved_obj.get('@id', '').startswith('/non_human_donors/'):
                                collection = 'donors'
                            elif resolved_obj.get('@id', '').startswith('/treatments/'):
                                collection = 'treatments'
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
        
        # Store raw matrix files in resolved_objects for later use
        if 'RawMatrixFile' not in self.resolved_objects:
            self.resolved_objects['RawMatrixFile'] = {}
        for rmf in raw_matrix_files:
            self.resolved_objects['RawMatrixFile'][rmf['@id']] = rmf
        
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
        
        print(f"Found {len(sequence_file_uuids)} sequence file UUIDs referenced by raw matrix files")
        
        # Step 3: Get sequence files and their file sets
        sequence_files = self.chunk_and_fetch('SequenceFile', list(sequence_file_uuids))
        
        # Store sequence files in resolved_objects
        if 'SequenceFile' not in self.resolved_objects:
            self.resolved_objects['SequenceFile'] = {}
        for sf in sequence_files:
            self.resolved_objects['SequenceFile'][sf['@id']] = sf
        
        print(f"Successfully fetched {len(sequence_files)} sequence files")
        
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
        
        print(f"Found {len(file_set_uuids)} file set UUIDs referenced by sequence files")
        
        # Step 4: Get file sets and their libraries
        file_sets = self.chunk_and_fetch('SequenceFileSet', list(file_set_uuids))
        
        # Store file sets in resolved_objects
        if 'SequenceFileSet' not in self.resolved_objects:
            self.resolved_objects['SequenceFileSet'] = {}
        for fs in file_sets:
            self.resolved_objects['SequenceFileSet'][fs['@id']] = fs
        
        print(f"Successfully fetched {len(file_sets)} file sets")
        
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
        
        print(f"Found {len(library_uuids)} library UUIDs referenced by file sets")
        
        # Step 5: Get libraries (determine types first, then fetch only what exists)
        droplet_uuids = set()  # Use sets to avoid duplicates
        plate_uuids = set()
        
        # Check which type each library is by looking at the file sets that reference them
        for file_set in self.resolved_objects.get('SequenceFileSet', {}).values():
            library_ref = file_set.get('library', '')
            if library_ref:
                if isinstance(library_ref, dict):
                    lib_id = library_ref.get('@id', '')
                else:
                    lib_id = library_ref
                
                lib_uuid = self.extract_uuid_from_id(lib_id)
                if lib_uuid in library_uuids:
                    # Determine type from the @id path
                    if '/droplet_based_libraries/' in lib_id:
                        droplet_uuids.add(lib_uuid)  # Use add() for sets
                    elif '/plate_based_libraries/' in lib_id:
                        plate_uuids.add(lib_uuid)
        
        # Only fetch the types that actually exist
        droplet_libraries = []
        plate_libraries = []
        
        if droplet_uuids:
            droplet_libraries = self.chunk_and_fetch('DropletBasedLibrary', list(droplet_uuids))
        
        if plate_uuids:
            plate_libraries = self.chunk_and_fetch('PlateBasedLibrary', list(plate_uuids))
        
        all_libraries = droplet_libraries + plate_libraries
        print(f"Successfully fetched {len(all_libraries)} libraries total " + 
              f"({len(droplet_libraries)} droplet, {len(plate_libraries)} plate)")
                
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
                'raw_matrix_files': [],
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
        
        # Step 9: Map raw matrix files to libraries
        self._map_raw_matrix_files_to_libraries(raw_matrix_files, libraries_data, all_samples)
        
        print(f"Structured data for {len(libraries_data)} libraries")
        
        return {
            'matrix_file_set': matrix_file_set,
            'libraries': libraries_data,
            'resolved_objects': self.resolved_objects
        }

    def _map_raw_matrix_files_to_libraries(self, raw_matrix_files, libraries_data, all_samples):
        """Map raw matrix files to their corresponding libraries using samples field or data flow"""
        
        for raw_file in raw_matrix_files:
            matched_libraries = []
            
            # Raw matrix file -> sequence files (via derived_from)
            derived_from = raw_file.get('derived_from', [])
            for seq_ref in derived_from:
                if isinstance(seq_ref, dict):
                    seq_id = seq_ref.get('@id', '')
                else:
                    seq_id = seq_ref
                
                # Find the sequence file object
                seq_file = self.resolved_objects.get('SequenceFile', {}).get(seq_id)
                if seq_file:
                    # Sequence file -> file sets (via sequence_file_sets)
                    file_set_refs = seq_file.get('sequence_file_sets', [])
                    for fs_ref in file_set_refs:
                        if isinstance(fs_ref, dict):
                            fs_id = fs_ref.get('@id', '')
                        else:
                            fs_id = fs_ref
                        
                        # Find the file set object
                        file_set = self.resolved_objects.get('SequenceFileSet', {}).get(fs_id)
                        if file_set:
                            # File set -> library
                            library_ref = file_set.get('library', '')
                            if library_ref:
                                if isinstance(library_ref, dict):
                                    lib_id = library_ref.get('@id', '')
                                else:
                                    lib_id = library_ref
                                
                                lib_uuid = self.extract_uuid_from_id(lib_id)
                                if lib_uuid in libraries_data and lib_uuid not in matched_libraries:
                                    matched_libraries.append(lib_uuid)
            
            # Add the raw matrix file to ALL matched libraries
            for lib_uuid in matched_libraries:
                libraries_data[lib_uuid]['raw_matrix_files'].append(raw_file)