from constants import (
    Configs,
    MAX_URL_LENGTH,
    BASE_URL_OVERHEAD,
    CONTROLLED_TERM_FIELDS,
)
from DB2_utils import (
    extract_uuid_from_id,
    extract_references_from_field,
    get_api_type_from_id,
    get_url_prefix_from_id,
)
import DB2lattice


class DB2Gatherer:
    def __init__(self, connection, configs: Configs):
        self.connection = connection
        self.configs = configs
        self.resolved_objects = {}  # {object_type: {id: object}}
    
    def chunk_and_fetch(self, obj_type, object_ids, filter_field='uuid', fields=None):
        """
        Fetch objects efficiently with URL chunking if needed

        filter_field: report query parameter used to select objects. Defaults to
            'uuid'. Controlled terms are selected by '@id' instead, because a
            controlled term reference carries a semantic term id rather than a
            uuid - see fetch_controlled_terms().
        fields: overrides the field list from OBJECT_CONFIG, for when only a
            subset of the profile is needed.
        """
        if not object_ids:
            return []

        config = None
        for cfg in self.configs.OBJECT_CONFIG.values():
            if cfg['api_type'] == obj_type:
                config = cfg
                break

        if not config:
            print(f"Warning: No config found for {obj_type}")
            return []

        field_lst = fields or config['fields']

        print(f"Fetching {len(object_ids)} {obj_type} objects...")

        # Remove duplicates
        unique_ids = list(set(object_ids))

        def build_filter(ids):
            return '&' + '&'.join([f"{filter_field}={oid}" for oid in ids])

        # Build filter URL
        filter_url = build_filter(unique_ids)

        # Single request if under limit
        if len(filter_url) <= MAX_URL_LENGTH:
            try:
                results = DB2lattice.get_report(
                    obj_type=obj_type,
                    filter_url=filter_url,
                    field_lst=field_lst,
                    connection=self.connection
                )
                return results or []
            except Exception as e:
                print(f"Error fetching {obj_type}: {e}")
                return []

        # Chunked requests
        print("URL too long, chunking...")
        all_results = []
        chunk_size = (MAX_URL_LENGTH - BASE_URL_OVERHEAD) // 50  # Rough estimate

        for i in range(0, len(unique_ids), chunk_size):
            chunk_ids = unique_ids[i:i + chunk_size]
            chunk_filter = build_filter(chunk_ids)

            try:
                chunk_results = DB2lattice.get_report(
                    obj_type=obj_type,
                    filter_url=chunk_filter,
                    field_lst=field_lst,
                    connection=self.connection
                )

                if chunk_results:
                    all_results.extend(chunk_results)

            except Exception as e:
                print(f"Chunk error: {e}")
                continue

        return all_results
        
    
    def resolve_references_for_samples(self, all_samples):
        """
        Collect all references first, then batch fetch by type

        Controlled terms are gathered last: the second pass below keeps finding
        new controlled term references on objects that only exist once the
        non-controlled-term fetch has completed, so fetching them any earlier
        would miss those.
        """
        all_reference_ids = {}  # {api_type: set(ids)}
        controlled_term_refs = set()  # {ref_path}

        # First pass: collect ALL references from samples
        for sample in all_samples.values():
            sample_api_type = get_api_type_from_id(sample['@id'], self.configs)
            config = None
            for cfg in self.configs.OBJECT_CONFIG.values():
                if cfg['api_type'] == sample_api_type:
                    config = cfg
                    break
            
            if not config:
                continue
            
            for field_name, ref_types in config.get('references', {}).items():
                field_value = sample.get(field_name)
                refs = extract_references_from_field(field_value, field_name, self.configs)
                
                if isinstance(ref_types, str):
                    ref_types = [ref_types]
                
                for ref in refs:
                    if ref.startswith('/'):
                        if 'controlled_terms' in ref_types and '/controlled_terms/' in ref:
                            # Keep the whole path - it is what the ControlledTerm
                            # fetch filters on, and what the flattener looks up by
                            controlled_term_refs.add(ref)
                        else:
                            # Regular UUID-based reference
                            api_type = get_api_type_from_id(ref, self.configs)
                            if api_type:
                                all_reference_ids.setdefault(api_type, set()).add(ref)
        
        # Batch fetch all non-controlled-term references by type
        for api_type, ref_ids in all_reference_ids.items():
            uuids = [extract_uuid_from_id(ref_id) for ref_id in ref_ids]
            ref_objects = self.chunk_and_fetch(api_type, uuids)
            
            if api_type not in self.resolved_objects:
                self.resolved_objects[api_type] = {}
            
            for obj in ref_objects:
                self.resolved_objects[api_type][obj['@id']] = obj
        
        # Second pass: collect controlled term references from non-sample objects
        for api_type, ref_dict in self.resolved_objects.items():
            if api_type == 'ControlledTerm':  # Skip the controlled terms dict
                continue

            for obj in ref_dict.values():
                obj_api_type = get_api_type_from_id(obj.get('@id', ''), self.configs)
                config = None
                for cfg in self.configs.OBJECT_CONFIG.values():
                    if cfg['api_type'] == obj_api_type:
                        config = cfg
                        break
                
                if not config:
                    continue
                
                for field_name, ref_types in config.get('references', {}).items():
                    if 'controlled_terms' in (ref_types if isinstance(ref_types, list) else [ref_types]):
                        field_value = obj.get(field_name)
                        refs = extract_references_from_field(field_value, field_name, self.configs)
                        
                        for ref in refs:
                            if ref.startswith('/controlled_terms/'):
                                controlled_term_refs.add(ref)

        # Now that no further controlled term references can turn up, fetch them
        self.resolved_objects['ControlledTerm'] = self.fetch_controlled_terms(
            controlled_term_refs
        )

    def fetch_controlled_terms(self, controlled_term_refs):
        """
        Fetch ControlledTerm objects, keyed by @id, so term_name is available

        Controlled terms are the one object type filtered by '@id' rather than
        'uuid': the reference carries a semantic term id, not a uuid, so there
        is no uuid available to filter on until the object has been fetched.
        """
        if not controlled_term_refs:
            return {}

        requested = sorted(controlled_term_refs)
        results = self.chunk_and_fetch(
            'ControlledTerm',
            requested,
            filter_field='@id',
            fields=CONTROLLED_TERM_FIELDS,
        )

        # If the '@id' filter were ever ignored, the response would be an
        # unfiltered dump of every controlled term in the instance - silently
        # wrong data rather than an error. Refuse it instead of flattening it.
        if len(results) > len(requested):
            raise RuntimeError(
                f"Requested {len(requested)} controlled terms but the API returned "
                f"{len(results)}. The '@id' filter does not appear to have been applied."
            )

        resolved = {obj['@id']: obj for obj in results if obj.get('@id')}

        missing = [ref for ref in requested if ref not in resolved]
        if missing:
            print(
                f"Warning: {len(missing)} of {len(requested)} controlled terms did not "
                f"resolve, e.g. {missing[:5]}"
            )

        print(f"Resolved {len(resolved)} controlled terms")
        return resolved

    def add_references_to_library(self, library_data, samples):
        """Add resolved non-controlled-term references to library data based on its samples"""
        added_refs= set()
        for sample in samples:
            sample_api_type = get_api_type_from_id(sample['@id'], self.configs)
            config = next(
                (cfg for cfg in self.configs.OBJECT_CONFIG.values() if cfg['api_type'] == sample_api_type),
                None,
            )
            if not config:
                continue
            for field_name, ref_types in config.get('references', {}).items():
                ref_type_list = [ref_types] if isinstance(ref_types, str) else ref_types
                if 'controlled_terms' in ref_type_list:
                    continue
                for ref in extract_references_from_field(sample.get(field_name), field_name, self.configs):
                    resolved_obj = None
                    for api_type, objects in self.resolved_objects.items():
                        if api_type != 'ControlledTerm' and ref in objects:
                            resolved_obj = objects[ref]
                            break
                    if resolved_obj and ref not in added_refs:
                        bucket = get_url_prefix_from_id(ref, self.configs)
                        if bucket:
                            library_data.setdefault(bucket, []).append(resolved_obj)
                            added_refs.add(ref)


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
                raw_matrix_uuids.append(extract_uuid_from_id(ref_id))
        
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
                    sequence_file_uuids.add(extract_uuid_from_id(ref_id))
        
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
                    file_set_uuids.add(extract_uuid_from_id(ref_id))
        
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
                    library_uuids.add(extract_uuid_from_id(lib_id))
        
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
                
                lib_uuid = extract_uuid_from_id(lib_id)
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
        sample_uuids_by_type = {}  # {api_type: [uuids]}
        
        for library in all_libraries:
            samples = library.get('samples', [])
            for ref in samples:
                if isinstance(ref, dict):
                    ref_id = ref.get('@id', '')
                else:
                    ref_id = ref
                
                if ref_id:
                    api_type = get_api_type_from_id(ref_id, self.configs)
                    if api_type:
                        if api_type not in sample_uuids_by_type:
                            sample_uuids_by_type[api_type] = []
                        sample_uuids_by_type[api_type].append(extract_uuid_from_id(ref_id))
        
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
                'raw_matrix_files': []
                # The rest of the library data information is created on demand by setdefault in add_references_to_library()
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
        self._map_raw_matrix_files_to_libraries(raw_matrix_files, libraries_data)
        
        print(f"Structured data for {len(libraries_data)} libraries")
        
        return {
            'matrix_file_set': matrix_file_set,
            'libraries': libraries_data,
            'resolved_objects': self.resolved_objects
        }

    def _map_raw_matrix_files_to_libraries(self, raw_matrix_files, libraries_data):
        """Map raw matrix files to their corresponding libraries using samples field or data flow"""
        
        for raw_file in raw_matrix_files:
            matched_data = {}
            
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
                                
                                lib_uuid = extract_uuid_from_id(lib_id)
                                # Record lib_uuid of libraries found, along with the sequence file and file sets used to get to it
                                if lib_uuid in libraries_data:
                                    bucket = matched_data.setdefault(lib_uuid, {'sequence_files': {}, 
                                                                                'sequence_file_sets': {}})
                                    bucket['sequence_files'][seq_file.get('@id', '')] = seq_file
                                    bucket['sequence_file_sets'][file_set.get('@id', '')] = file_set
            
            # Add a sequence_file and sequence_file_set path specific copy of raw matrix file to each library
            for lib_uuid, bucket in matched_data.items():
                raw_file_copy = dict(raw_file)
                raw_file_copy['sequence_files'] = list(bucket['sequence_files'].values())
                raw_file_copy['sequence_file_sets'] = list(bucket['sequence_file_sets'].values())
                libraries_data[lib_uuid]['raw_matrix_files'].append(raw_file_copy)
