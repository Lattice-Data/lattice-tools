import os
import argparse
import sys
import pandas as pd
from datetime import datetime
from DB2Gatherer import DB2Gatherer

# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

import DB2lattice

class DB2Flattener:
    def __init__(self):
        # Setup connection
        os.environ['DEMO_KEY'] = 'HKA345NO'
        os.environ['DEMO_SECRET'] = 'ar6stvgd7epcxirx'
        os.environ['DEMO_SERVER'] = 'https://lattice-api-dev.demo.lattice-data.org'
        self.connection = DB2lattice.Connection('demo')
        
        # Initialize gatherer
        self.gatherer = DB2Gatherer(self.connection)
    
    def flatten_matrix_file_set(self, matrix_file_set_uuid, output_file=None):
        """
        Flatten a MatrixFileSet into library-indexed data and save as CSV file
        """
        print(f"Processing MatrixFileSet {matrix_file_set_uuid}")
        
        # Gather all data
        complete_data = self.gatherer.gather_complete_library_data(matrix_file_set_uuid)
        
        if not complete_data:
            print("Error: No data gathered")
            return None
        
        print("Creating DataFrame...")
        
        # Generate output filename if not provided
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}.csv"
        
        # Create DataFrame and drop columns with all "n/a"
        df = self.create_dataframe(complete_data)
        cols_to_drop = [col for col in df.columns if df[col].eq('n/a').all()]
        df = df.drop(columns=cols_to_drop)
        
        # Save to CSV
        print(f"Saving to {output_file}...")
        df.to_csv(output_file, index=False)
        
        print(f"✅ CSV file created: {output_file}")
        print(f"   Rows: {len(df)}")
        print(f"   Columns: {len(df.columns)}")
        
        return output_file
        
    def _filter_to_gex_libraries(self, libraries):
        """Filter libraries to only include Gene Expression when there are pairs"""
        gex_libraries = []
        
        for lib in libraries:
            feature_types = lib.get('feature_types', [])
            lib_type = lib.get('@type', [''])[0] if lib.get('@type') else ''
            
            # Check if this library has Gene Expression
            if 'Gene Expression' in feature_types:
                gex_libraries.append(lib)
            elif not feature_types:  # Empty list or missing field
                # For droplet-based libraries without feature_types, assume non-GEX
                if 'DropletBasedLibrary' in lib_type:
                    continue  # Skip droplet libraries without feature_types
                # For plate-based libraries without feature_types, assume GEX until field is added
                elif 'PlateBasedLibrary' in lib_type:
                    gex_libraries.append(lib)
        
        # Always return only GEX libraries
        if gex_libraries:
            print(f"DEBUG: Filtered to {len(gex_libraries)} Gene Expression libraries out of {len(libraries)} total")
            return gex_libraries
        
        # Fallback: if no GEX libraries found, return all (shouldn't happen normally)
        print(f"WARNING: No Gene Expression libraries found, keeping all libraries")
        return libraries
    
    def create_dataframe(self, complete_data):
        """Create a flattened DataFrame with library or raw matrix file as rows"""
        libraries_data = complete_data['libraries']
        resolved_controlled_terms = complete_data['resolved_objects'].get('ControlledTerm', {})

        # Filter libraries to GEX once upfront
        all_libraries = [lib_data['library'] for lib_data in libraries_data.values()]
        filtered_libraries = self._filter_to_gex_libraries(all_libraries)
        filtered_library_ids = {lib.get('@id') for lib in filtered_libraries}
        
        # Filter libraries_data to only include GEX libraries
        filtered_libraries_data = {
            lib_uuid: lib_data for lib_uuid, lib_data in libraries_data.items()
            if lib_data['library'].get('@id') in filtered_library_ids
        }
        
        # Check if any raw matrix files have samples field
        has_raw_file_samples = False
        for lib_data in filtered_libraries_data.values():
            for raw_file in lib_data['raw_matrix_files']:
                if raw_file.get('samples'):
                    has_raw_file_samples = True
                    break
            if has_raw_file_samples:
                break
        
        print(f"Creating DataFrame by {'raw matrix file' if has_raw_file_samples else 'library'}...")
        
        rows = []
        
        if has_raw_file_samples:
            # Raw matrix file-based rows when samples field is present
            # First, collect all raw matrix files and their associated libraries
            raw_file_to_libraries = {}
            
            for lib_data in filtered_libraries_data.values():
                library = lib_data['library']
                samples = lib_data['samples']
                raw_matrix_files = lib_data['raw_matrix_files']
                
                for raw_file in raw_matrix_files:
                    raw_file_id = raw_file.get('@id')
                    if raw_file_id not in raw_file_to_libraries:
                        raw_file_to_libraries[raw_file_id] = {
                            'raw_file': raw_file,
                            'libraries': [],
                            'all_samples': [],
                            'all_donors': [],
                            'all_treatments': [],
                            'all_genetic_modifications': [],
                            'all_experimental_conditions': []
                        }
                    
                    # Add this library's data to the raw matrix file
                    raw_file_to_libraries[raw_file_id]['libraries'].append(library)
                    raw_file_to_libraries[raw_file_id]['all_samples'].extend(samples)
                    raw_file_to_libraries[raw_file_id]['all_donors'].extend(lib_data['donors'])
                    raw_file_to_libraries[raw_file_id]['all_treatments'].extend(lib_data['treatments'])
                    raw_file_to_libraries[raw_file_id]['all_genetic_modifications'].extend(lib_data['genetic_modifications'])
                    raw_file_to_libraries[raw_file_id]['all_experimental_conditions'].extend(lib_data['experimental_conditions'])
          
            # Create one row per raw matrix file
            for file_data in raw_file_to_libraries.values():
                raw_file = file_data['raw_file']
                libraries = file_data['libraries']
                samples = file_data['all_samples']
                donors = file_data['all_donors']
                treatments = file_data['all_treatments']
                genetic_modifications = file_data['all_genetic_modifications']
                experimental_conditions = file_data['all_experimental_conditions']
                
                # Create sample aliases list
                sample_aliases = []
                for sample_ref in raw_file.get('samples', []):
                    sample_obj = next((s for s in samples if s.get('@id') == sample_ref), None)
                    if sample_obj:
                        alias = self._get_clean_alias(sample_obj)
                        sample_aliases.append(alias)
                
                row = {
                    'raw_matrix_file_alias': self._get_clean_alias(raw_file),
                    'library_aliases': self._join_unique([self._get_clean_alias(lib) for lib in libraries]),
                    
                    # Library information (from all associated libraries)
                    'library_protocols': self._join_unique([lib.get('protocol', '') for lib in libraries]),
                    'library_types': self._join_unique([lib.get('@type', [''])[0] if 
                                                        lib.get('@type') else '' for lib in libraries]),
                    
                    # Sample information
                    'raw_file_samples': self._join_unique(sample_aliases),
                    'library_sample_types': self._join_unique([s.get('@type', [''])[0] for s in samples]),
                    'selection_markers': self._join_unique([s.get('selection_markders', [''])[0] for s in samples]),
                    
                    # Donor information
                    'donor_sexes': self._join_unique([d.get('sex', '') for d in donors]),
                    'donor_ethnicities': self._join_unique([
                        self._resolve_controlled_term(d.get('ethnicity'), resolved_controlled_terms)
                        for d in donors
                    ]),
                    'donor_taxa': self._join_unique([d.get('taxa', '') for d in donors]),
                    
                    # Treatment information
                    'treatment_terms': self._join_unique([
                        self._resolve_controlled_term(t.get('ontological_term'), resolved_controlled_terms)
                        for t in treatments
                    ]),
                    
                    # Experimental condition information
                    'experimental_condition': self._join_unique([e.get('condition') for e in experimental_conditions]),
                    'experimental_condition_details': self._join_unique([e.get('text_value') for e in experimental_conditions]),

                    # Genetic modification information
                    'genetic_modifications' : self._join_unique([g.get('strategy') for g in genetic_modifications]),


                    # Cell type information from samples
                    'enriched_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'enriched_cell_types', 
                                                                             resolved_controlled_terms),
                    'depleted_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'depleted_cell_types', 
                                                                             resolved_controlled_terms),
                    'intended_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'intended_cell_types', 
                                                                             resolved_controlled_terms),
                    
                    # Disease information from samples
                    'diseases': self._get_diseases_from_samples(samples, resolved_controlled_terms),
                }
                
                rows.append(row)
        
        else:
            # Library-based rows when no samples field in raw matrix files
            for lib_data in filtered_libraries_data.values():
                library = lib_data['library']
                samples = lib_data['samples']
                donors = lib_data['donors']
                treatments = lib_data['treatments']
                genetic_modifications = lib_data['genetic_modifications']
                raw_matrix_files = lib_data['raw_matrix_files']
                
                # Build the row with library alias first
                row = {
                    'library_alias': self._get_clean_alias(library),
                    'raw_matrix_file_alias': self._get_clean_alias(raw_matrix_files),
                    
                    # Library information
                    'library_protocol': library.get('protocol', 'n/a'),
                    'library_type': library.get('@type', [''])[0] if library.get('@type') else 'n/a',
                    
                    # Sample information
                    'sample_types': self._join_unique([s.get('@type', ['n/a'])[0] for s in samples]),
                    
                    # Donor information
                    'donor_sexes': self._join_unique([d.get('sex', 'n/a') for d in donors]),
                    'donor_ethnicities': self._join_unique([
                        self._resolve_controlled_term(d.get('ethnicity'), resolved_controlled_terms)
                        for d in donors
                    ]),
                    
                    # Treatment information
                    'treatment_terms': self._join_unique([
                        self._resolve_controlled_term(t.get('ontological_term'), resolved_controlled_terms)
                        for t in treatments
                    ]),
                    
                    # Cell type information from samples
                    'enriched_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'enriched_cell_types', 
                                                                             resolved_controlled_terms),
                    'depleted_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'depleted_cell_types', 
                                                                             resolved_controlled_terms),
                    'intended_cell_types': self._get_cell_types_from_samples(samples, 
                                                                             'intended_cell_types', 
                                                                             resolved_controlled_terms),
                    
                    # Disease information from samples
                    'diseases': self._get_diseases_from_samples(samples, resolved_controlled_terms),
                    
                }
                
                rows.append(row)

        return pd.DataFrame(rows)
    
    def _resolve_controlled_term(self, term_ref, resolved_controlled_terms):
        """Resolve a controlled term reference to its term_id (semantic identifier)"""
        if not term_ref:
            return ''
        
        # Handle if it's already a dict with @id
        if isinstance(term_ref, dict):
            term_id = term_ref.get('@id', '')
        else:
            term_id = term_ref
        
        # Look up in the resolved controlled terms dictionary
        return resolved_controlled_terms.get(term_id, '')
    
    def _get_cell_types_from_samples(self, samples, field_name, controlled_terms):
        """Extract cell type terms from samples"""
        cell_types = []
        
        for sample in samples:
            field_value = sample.get(field_name, [])
            for ref in field_value:
                term_name = self._resolve_controlled_term(ref, controlled_terms)
                if term_name:
                    cell_types.append(term_name)
        
        return self._join_unique(cell_types)
    
    def _get_diseases_from_samples(self, samples, controlled_terms):
        """Extract disease terms from samples"""
        diseases = []
        
        for sample in samples:
            disease_refs = sample.get('diseases', [])
            for ref in disease_refs:
                term_name = self._resolve_controlled_term(ref, controlled_terms)
                if term_name:
                    diseases.append(term_name)
        
        return self._join_unique(diseases)
    
    def _join_unique(self, items):
        """Join unique non-empty items with semicolon, return 'n/a' if empty"""
        # Filter out empty/None values
        filtered_items = [str(item).strip() for item in items if item and str(item).strip()]
        
        if not filtered_items:
            return 'n/a'
        
        # Remove duplicates and sort
        unique_items = sorted(set(filtered_items))
        return '; '.join(unique_items)
    
    def _join_list(self, items):
        """Join list items directly, return 'n/a' if empty"""
        if not items:
            return 'n/a'
        return '; '.join(items)
    
    def _get_clean_alias(self, obj_or_list):
        """Extract alias from an object or list of objects, cleaning the result"""
        if not obj_or_list:
            return 'n/a'
        
        aliases = []
        
        # Handle both single objects and lists
        if isinstance(obj_or_list, list):
            # List of objects
            for obj in obj_or_list:
                if isinstance(obj, dict):
                    obj_aliases = obj.get('aliases', [])
                    if obj_aliases:
                        aliases.extend(obj_aliases)
        else:
            # Single object
            if isinstance(obj_or_list, dict):
                obj_aliases = obj_or_list.get('aliases', [])
                if obj_aliases:
                    aliases.extend(obj_aliases)
        
        # If no aliases found, return 'n/a'
        if not aliases:
            return 'n/a'
        
        # Join unique aliases and clean
        result = self._join_unique(aliases)
        
        # Clean the result (remove prefix ending with ':')
        if ':' in result:
            return result.split(':', 1)[1]  # Take everything after the first ':'
        return result

def main():
    parser = argparse.ArgumentParser(
        description="Flatten MatrixFileSet data for DB2 processing and export to CSV"
    )
    
    parser.add_argument(
        '--uuid', '-u',
        required=True,
        help='UUID of the MatrixFileSet to process'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Output CSV filename (optional, auto-generates if not provided)'
    )
    
    args = parser.parse_args()
    
    try:
        flattener = DB2Flattener()
        output_file = flattener.flatten_matrix_file_set(args.uuid, args.output)
        
        if output_file:
            print(f"\n✅ Success! CSV file created: {output_file}")
        else:
            print(f"❌ Failed to process MatrixFileSet {args.uuid}")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()