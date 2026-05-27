import os
import argparse
import sys
import pandas as pd
from datetime import datetime
from DB2Gatherer import DB2Gatherer

# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

import lattice

class DB2Flattener:
    def __init__(self):
        # Setup connection
        os.environ['DEMO_KEY'] = 'HKA345NO'
        os.environ['DEMO_SECRET'] = 'ar6stvgd7epcxirx'
        os.environ['DEMO_SERVER'] = 'https://lattice-api-dev.demo.lattice-data.org'
        self.connection = lattice.Connection('demo')
        
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
        
        # Create DataFrame
        df = self.create_dataframe(complete_data)
        
        # Save to CSV
        print(f"Saving to {output_file}...")
        df.to_csv(output_file, index=False)
        
        print(f"✅ CSV file created: {output_file}")
        print(f"   Rows: {len(df)}")
        print(f"   Columns: {len(df.columns)}")
        
        return output_file
    
    def create_dataframe(self, complete_data):
        """Create a flattened DataFrame with library alias as first column"""
        matrix_file_set = complete_data['matrix_file_set']
        libraries_data = complete_data['libraries']
        resolved_controlled_terms = complete_data['resolved_objects'].get('ControlledTerm', {})
        
        rows = []
        
        for lib_uuid, lib_data in libraries_data.items():
            library = lib_data['library']
            samples = lib_data['samples']
            donors = lib_data['donors']
            treatments = lib_data['treatments']
            genetic_modifications = lib_data['genetic_modifications']
            
            # Build the row with library alias first
            row = {
                # Library alias FIRST
                'library_alias': self._get_clean_alias(library),
                
                # Library info
                'library_chemistry_version': library.get('chemistry_version', 'n/a'),
                'library_cardinality': library.get('library_cardinality', 'n/a'),
                'library_feature_types': self._join_list(library.get('feature_types', [])),
                'library_multiplexing_method': library.get('multiplexing_method', 'n/a'),
                'library_kit_version': library.get('kit_version', 'n/a'),
                
                # Matrix File Set info
                'matrix_file_set_alias': self._get_clean_alias(matrix_file_set),
                'matrix_file_set_genome_assembly': matrix_file_set.get('genome_assembly', 'n/a'),
                'matrix_file_set_genome_annotation': matrix_file_set.get('genome_annotation', 'n/a'),
                'matrix_file_set_software': matrix_file_set.get('software', 'n/a'),
                'matrix_file_set_software_version': matrix_file_set.get('software_version', 'n/a'),
                
                # Donor information
                'donor_sexes': self._join_unique([d.get('sex', '') for d in donors]),
                'donor_ethnicities': self._join_unique([
                    self._resolve_controlled_term(d.get('ethnicity'), resolved_controlled_terms)
                    for d in donors
                ]),
                
                # Treatment information
                'treatment_terms': self._join_unique([
                    self._resolve_controlled_term(t.get('ontological_term'), resolved_controlled_terms)
                    for t in treatments
                ]),
                
                # Sample-level information
                'sample_dates': self._join_unique([
                    s.get('date_obtained', '') for s in samples
                ]),
                'age_lower_bounds': self._join_unique([
                    str(s.get('lower_bound_age', '')) for s in samples 
                    if s.get('lower_bound_age') is not None
                ]),
                'age_upper_bounds': self._join_unique([
                    str(s.get('upper_bound_age', '')) for s in samples 
                    if s.get('upper_bound_age') is not None
                ]),
                'age_units': self._join_unique([s.get('age_units', '') for s in samples]),
                
                 # Cell type information from samples
                'enriched_cell_types': self._get_cell_types_from_samples(samples, 'enriched_cell_types', resolved_controlled_terms),
                'depleted_cell_types': self._get_cell_types_from_samples(samples, 'depleted_cell_types', resolved_controlled_terms),
                'intended_cell_types': self._get_cell_types_from_samples(samples, 'intended_cell_types', resolved_controlled_terms),
                
                # Disease information from samples
                'diseases': self._get_diseases_from_samples(samples, resolved_controlled_terms),
                
                # Genetic modifications
                'genetic_modification_modalities': self._join_unique([
                    gm.get('modality', '') for gm in genetic_modifications
                ]),
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
    
    def _get_clean_alias(self, obj):
        """Extract clean alias from object, stripping lab name prefix"""
        aliases = obj.get('aliases', [])
        if not aliases:
            return obj.get('uuid', 'n/a')
        
        # Take first alias
        full_alias = aliases[0]
        if not full_alias:
            return 'n/a'
        
        # Strip lab name (everything before and including ':')
        if ':' in full_alias:
            clean_alias = full_alias.split(':', 1)[1]
        else:
            clean_alias = full_alias
        
        return clean_alias.strip() or 'n/a'

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