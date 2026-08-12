import argparse
import sys
import re
import pandas as pd
import numpy as np
from datetime import datetime
from DB2Gatherer import DB2Gatherer
from constants import (
    Configs,
    PROP_MAP_GEO,
    PROP_MAP_BIOHUB,
    TISSUE_TYPE_MAP,
    GENETIC_PERTURBATION_MAP,
    REFORMAT_LIST,
)
from DB2_utils import (
    split_controlled_term_columns, 
    collapse_dataframe,
    get_config_obj_type,
    get_url_prefix_from_id,
    extract_references_from_field,
    combine_bound_columns,
    collapse_duplicate_columns,
    strip_author_metadata_column_prefix,
    is_empty,
)
from generate_constants import load_and_return_constant_dicts
import DB2lattice



class DB2Flattener:
    def __init__(self, connection: DB2lattice.Connection, configs: Configs):
        # Setup connection
        self.connection: DB2lattice.Connection = connection
        self.configs = configs
        
        # Initialize gatherer
        self.gatherer = DB2Gatherer(self.connection, configs)

    
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
        
        print("Creating DataFrames...")
        
        # Generate output filename if not provided
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}_MAIN.csv"
        
        # Create main DataFrame and sample DataFrame
        main_df, sample_df = self.create_dataframe(complete_data)
        main_df = main_df.dropna(axis=1, how='all')
        sample_df = sample_df.dropna(axis=1, how='all')

        # Split dict columns into _term_id / _term_name before writing CSV
        main_df = split_controlled_term_columns(main_df)
        if sample_df is not None and not sample_df.empty:
            sample_df = split_controlled_term_columns(sample_df)

        
        # Save main DataFrame to CSV
        print(f"Saving main DataFrame to {output_file}...")
        main_df.to_csv(output_file, index=False)
        
        print(f"✅ Main CSV file created: {output_file}")
        print(f"   Rows: {len(main_df)}")
        print(f"   Columns: {len(main_df.columns)}")

        # Create Biohub DataFrame from main and sample df
        biohub_output = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}_BIOHUB.csv"
        biohub_df = self.create_biohub_dataframe(main_df)
        print(f"Saving biohub DataFrame to {biohub_output}...")
        biohub_df.to_csv(biohub_output, index=False)
        print(f"✅ Biohub CSV file created: {biohub_output}")

        # Create GEO DataFrame from main DataFrame
        geo_output = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}_GEO.csv"
        geo_df = self.create_geo_dataframe(main_df)
        print(f"Saving geo DataFrame to {geo_output}...")
        geo_df.to_csv(geo_output, index=False)
        print(f"✅ GEO CSV file created: {geo_output}")
        
        # Save sample DataFrame if it exists
        if sample_df is not None and not sample_df.empty:
            sample_output = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}_SAMPLES.csv"
            
            print(f"Saving sample DataFrame to {sample_output}...")
            sample_df.to_csv(sample_output, index=True)
            
            print(f"✅ Sample CSV file created: {sample_output}")
            print(f"   Rows: {len(sample_df)}")
            print(f"   Columns: {len(sample_df.columns)}")
        
        return output_file
    
    
    def create_dataframe(self, complete_data):
        """Create main library/raw-file DataFrame and sample DataFrame keyed by raw matrix file"""
        libraries_data = complete_data['libraries']
        resolved_controlled_terms = complete_data['resolved_objects'].get('ControlledTerm', {})
        
        print("Creating DataFrame by raw matrix file...")
        
        rows = []
        sample_df = None
        new_sample_df = None
        
        # Raw matrix file-based rows - samples field is always present
        # First, collect all raw matrix files and their associated libraries
        raw_file_to_libraries = {}
        sample_metadata = {}  # One row per unique sample, keyed by sample_alias
        
        for lib_data in libraries_data.values():
            library = lib_data['library']
            samples = lib_data['samples']
            raw_matrix_files = lib_data['raw_matrix_files']
            
            for raw_file in raw_matrix_files:
                raw_file_id = raw_file.get('@id')
                if raw_file_id not in raw_file_to_libraries:
                    raw_file_to_libraries[raw_file_id] = {
                        'library_entries': [],
                        'all_samples': []
                    }
                
                # Add this library's data to the raw matrix file
                existing_lib_ids = {entry['library'].get('@id') for entry in raw_file_to_libraries[raw_file_id]['library_entries']}
                if library.get('@id') not in existing_lib_ids:
                    raw_file_to_libraries[raw_file_id]['library_entries'].append({'library': library, 'raw_file': raw_file})
                raw_file_to_libraries[raw_file_id]['all_samples'].extend(samples)
                
                # Collect per-sample metadata for later merge with raw matrix files
                for sample_ref in raw_file.get('samples', []):
                    sample_obj = next((s for s in samples if s.get('@id') == sample_ref), None)
                    if sample_obj:
                        sample_alias = self._get_clean_alias(sample_obj)
                        if sample_alias in sample_metadata:
                            continue

                        sample_metadata[sample_alias] = {}

                        sample_type = get_config_obj_type(sample_obj, self.configs)
                        sample_config = self.configs.OBJECT_CONFIG[sample_type]
                        sample_refs = sample_config.get('references', {})
                        for field in sample_config.get('fields', []):
                            value = sample_obj.get(field)
                            # Special handling for author_metadata dictionary
                            if field == 'author_metadata' and isinstance(value, dict):
                                for key, val in value.items():
                                    field_name = f"{sample_type}_{field}_{'_'.join(key.split(' '))}"
                                    sample_metadata[sample_alias][field_name] = val
                            else:
                                if self._is_controlled_term_field(field, sample_refs):
                                    value = self._resolve_controlled_term_value(
                                        value, resolved_controlled_terms
                                    )
                                field_name = f'{sample_type}_{field}'
                                sample_metadata[sample_alias][field_name] = value

                        self._flatten_resolved_references(
                            sample_obj, lib_data, sample_metadata, sample_alias, resolved_controlled_terms
                        )

        # Create one row per (raw matrix file, library)
        for file_data in raw_file_to_libraries.values():
            library_entries = file_data['library_entries']
            samples = file_data['all_samples']

            # For raw_matrix_file_alias and raw_file_samples, it doesn't matter which raw file copy we use
            # As they will all have the same value
            representative_raw_file = library_entries[0]['raw_file']
            sample_aliases = []
            for sample_ref in representative_raw_file.get('samples', []):
                sample_obj = next((s for s in samples if s.get('@id') == sample_ref), None)
                if sample_obj:
                    sample_aliases.append(self._get_clean_alias(sample_obj))

            shared = {
                'raw_matrix_file_alias': self._get_clean_alias(representative_raw_file),
                'raw_file_samples': self._join_unique(sample_aliases),
            }

            for entry in library_entries:
                lib = entry['library']
                # This libraries scoped copy of the raw matrix file
                # Which has the library specific sequence_file and sequence_file_set metadata
                raw_file = entry['raw_file']
                row = dict(shared)
                lib_type = get_config_obj_type(lib, self.configs)
                for field in self.configs.OBJECT_CONFIG[lib_type].get('fields', []):
                    row[f'{lib_type}_{field}'] = lib.get(field)

                for obj_type in ('sequence_files', 'sequence_file_sets'):
                    obj_config = self.configs.OBJECT_CONFIG.get(obj_type, {})
                    for obj_field in obj_config.get('fields', []):
                        values = [
                            obj.get(obj_field) for obj in raw_file.get(obj_type, [])
                            if obj.get(obj_field) not in (None, '')
                        ]
                        row[f'{obj_type}_{obj_field}'] = self._join_unique(values) if values else None
                rows.append(row)
        
        main_df = pd.DataFrame(rows)
        
        # Merge sample metadata with main DataFrame columns
        # Samples is one row per unique (raw matrix file + sample)
        sample_df = None
        if sample_metadata and not main_df.empty:
            per_sample_df = pd.DataFrame.from_dict(sample_metadata, orient='index')
            per_sample_df.index.name = 'sample_alias'
            per_sample_df = per_sample_df.reset_index()

            rmf_sample_keys = main_df[['raw_matrix_file_alias','raw_file_samples']].drop_duplicates()
            
            sample_df = per_sample_df.merge(
                rmf_sample_keys,
                left_on = 'sample_alias',
                right_on ='raw_file_samples', 
                how = 'right'
            )
            
            new_sample_df = sample_df.set_index('raw_matrix_file_alias')

            main_df = main_df.merge(
                per_sample_df,
                left_on = 'raw_file_samples',
                right_on = 'sample_alias',
                how = 'left'
            )
            
            print(f"Creating sample DataFrame with {len(sample_df)} rows "
                  f"({sample_df['raw_matrix_file_alias'].nunique()} raw matrix files)...")
        
        return main_df, new_sample_df

    def _row_is_gex(self, row) -> bool:
        """ Filter df to only GEX libraries"""
        droplet_ft = row.get('droplet_based_libraries_feature_types')
        plate_ft = row.get('plate_based_libraries_feature_types')

        def has_gex(ft):
            if ft is None or (isinstance(ft, float) and pd.isna(ft)):
                return None  # missing
            if isinstance(ft, str):
                return 'Gene Expression' in ft
            if isinstance(ft, list):
                return 'Gene Expression' in ft
            return 'Gene Expression' in str(ft)

        droplet = has_gex(droplet_ft)
        plate = has_gex(plate_ft)

        if droplet is True or plate is True:
            return True
        if droplet is False or plate is False:
            return False

        # Missing feature_types: plate assumed GEX, droplet assumed non-GEX
        if pd.notna(row.get('plate_based_libraries_@id')) or pd.notna(
            row.get('plate_based_libraries_CRO_group_identifier')
        ):
            return True
        if pd.notna(row.get('droplet_based_libraries_@id')) or pd.notna(
            row.get('droplet_based_libraries_CRO_group_identifier')
        ):
            return False
        return True  # if no GEX found, keep all

    def create_geo_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build GEO submission dataframe from already-split main_df, taking GEX libraries only

        Expects _term_name columns (not raw dict columns).
        """
        gex_mask = main_df.apply(self._row_is_gex, axis=1)
        geo_source = main_df[gex_mask].copy()
        print(f"GEO: filtered to {len(geo_source)} GEX rows out of {len(main_df)} MAIN rows")

        subset_keys = [k for k in PROP_MAP_GEO if k in geo_source.columns]
        geo_df = geo_source[subset_keys].copy()
        geo_df.rename(columns=PROP_MAP_GEO, inplace=True)

        group_col = "*library name"
        return collapse_dataframe(geo_df, group_col=group_col)


    def create_biohub_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build Biohub samples dataframe from main_df
        """
        columns_to_keep = [k for k in PROP_MAP_BIOHUB if k in main_df.columns]
        columns_to_keep.extend([k for k in main_df.columns if re.search('_author_metadata_', k)])
        biohub_df = main_df[columns_to_keep].copy()
        biohub_df.rename(columns=PROP_MAP_BIOHUB, inplace=True)
        biohub_df = strip_author_metadata_column_prefix(biohub_df)
        biohub_df = collapse_duplicate_columns(biohub_df)

        # Deduplicate rows, must join lists as they are not hashable
        list_cols = [
            col for col in biohub_df.columns 
            if biohub_df[col].dropna().map(lambda x: isinstance(x, list)).any()
        ]
        for col in list_cols:
            biohub_df[col] = biohub_df[col].apply(lambda x: tuple(x) if isinstance(x, list) else x)
        biohub_df.drop_duplicates(inplace=True)

        # Add columns default values if not present
        if "disease" not in biohub_df.columns:
            biohub_df["disease"] = "normal"
        else:
            biohub_df["disease"] = biohub_df["disease"].apply(
                lambda v: pd.NA if (v is None or v == "" or v == [] or v == ()) else v
            )
            biohub_df["disease"] = biohub_df["disease"].fillna("normal")


        if "self_reported_ethnicity" not in biohub_df.columns:
            biohub_df["self_reported_ethnicity"] = np.where(biohub_df["organism"] == "Homo sapiens", "unknown", "na")
        else:
            biohub_df.loc[biohub_df["self_reported_ethnicity"].isna() & (biohub_df["organism"] == "Homo sapiens"), "self_reported_ethnicity"] = "unknown"
            biohub_df.loc[biohub_df["self_reported_ethnicity"].isna() & (biohub_df["organism"] != "Homo sapiens"), "self_reported_ethnicity"] = "na"

        # Combine multiple columns into one
        biohub_df = combine_bound_columns(
            biohub_df,
            lower_col="experimental_conditions_lower_bound_duration",
            upper_col="experimental_conditions_upper_bound_duration",
            units_col="experimental_conditions_duration_units",
            out_col="experimental_perturbation_time_point",
        )
        biohub_df = combine_bound_columns(
            biohub_df,
            lower_col="tissues_lower_bound_age",
            upper_col="tissues_upper_bound_age",
            units_col="tissues_age_units",
            out_col="age",
        )

        # Update values to match schema
        biohub_df["tissue_type"] = biohub_df["tissue_type"].apply(lambda x: TISSUE_TYPE_MAP.get(x[0], x))
        if "genetic_perturbation_strategy" in biohub_df.columns:
            biohub_df["genetic_perturbation_strategy"] = biohub_df["genetic_perturbation_strategy"].replace(GENETIC_PERTURBATION_MAP)
        for field in REFORMAT_LIST:
            if field in biohub_df.columns:
                biohub_df[field] = biohub_df[field].apply(
                    lambda x: "|".join(map(str, x)) if isinstance(x, (list, tuple)) else x
                )

        return biohub_df


    def _flatten_resolved_references(
        self, sample_obj, lib_data, sample_metadata, sample_alias, resolved_controlled_terms
    ):
        """Flatten resolved reference objects into sample_metadata columns."""
        sample_type = get_config_obj_type(sample_obj, self.configs)
        config = self.configs.OBJECT_CONFIG.get(sample_type, {})

        refs_by_prefix = {}
        for field_name, ref_types in config.get('references', {}).items():
            ref_type_list = [ref_types] if isinstance(ref_types, str) else ref_types
            if 'controlled_terms' in ref_type_list:
                continue

            refs = extract_references_from_field(
                sample_obj.get(field_name), field_name, self.configs
            )
            for ref in refs:
                prefix = get_url_prefix_from_id(ref, self.configs)
                if prefix:
                    refs_by_prefix.setdefault(prefix, []).append(ref)

        for url_prefix, refs in refs_by_prefix.items():
            seen = set()
            unique_refs = [r for r in refs if not (r in seen or seen.add(r))]

            resolved_objs = [
                obj for ref in unique_refs
                for obj in lib_data.get(url_prefix, [])
                if obj.get('@id') == ref
            ]

            obj_config = self.configs.OBJECT_CONFIG.get(url_prefix, {})
            obj_refs = obj_config.get('references', {})

            for obj_field in obj_config.get('fields', []):
                col = f'{url_prefix}_{obj_field}'
                field_ref_types = obj_refs.get(obj_field, [])
                if isinstance(field_ref_types, str):
                    field_ref_types = [field_ref_types]
                is_controlled_term = 'controlled_terms' in field_ref_types

                values = []
                for obj in resolved_objs:
                    value = obj.get(obj_field)
                    if value is None:
                        continue
                    if is_controlled_term:
                        value = self._resolve_controlled_term_value(
                            value, resolved_controlled_terms
                        )
                    if value not in (None, ''):
                        values.append(value)

                # Specific handling for author metadata fields
                if obj_field == 'author_metadata':
                    values_by_key = {}
                    for value in values:
                        if isinstance(value, dict):
                            for key, val in value.items():
                                sanitized_key = '_'.join(key.split(' '))
                                values_by_key.setdefault(sanitized_key, []).append(val)
                    for key, vals in values_by_key.items():
                        col = f"{url_prefix}_{obj_field}_{key}"
                        sample_metadata[sample_alias][col] = self._join_unique(vals)

                elif is_controlled_term:
                    # _join_unique() would str() these dicts into the cell
                    sample_metadata[sample_alias][col] = self._dedupe_terms(values)

                else:
                    sample_metadata[sample_alias][col] = (
                        self._join_unique(values) if values else None
                    )

    def _resolve_controlled_term(self, term_ref, resolved_controlled_terms):
        """
        Resolve a controlled term reference to the ControlledTerm object

        Returns the object ({'@id': ..., 'term_name': ...}) rather than the bare
        term id, so split_controlled_term_columns() can emit both a _term_id and
        a _term_name column. Some reference fields already arrive as embedded
        dicts carrying term_name; those are passed through untouched, so 
        both reference shapes end up identical downstream.
        """
        if not term_ref:
            return None

        if isinstance(term_ref, dict):
            # Already embedded with a term_name - nothing to look up
            if term_ref.get('term_name'):
                return term_ref
            ref_id = term_ref.get('@id', '')
        else:
            ref_id = term_ref

        return resolved_controlled_terms.get(ref_id)

    @staticmethod
    def _is_controlled_term_field(field_name, references):
        """True if OBJECT_CONFIG marks this field as a controlled_terms reference"""
        ref_types = references.get(field_name, [])
        if isinstance(ref_types, str):
            ref_types = [ref_types]
        return 'controlled_terms' in ref_types

    def _resolve_controlled_term_value(self, value, resolved_controlled_terms):
        """
        Resolve a controlled term field value to a dict, or list of dicts

        Deliberately not passed through _join_unique(): keeping dict shape is
        what lets split_controlled_term_columns() split the column into
        _term_id / _term_name. Array-typed fields keep list shape;
        split_term_cell() collapses a single-element list to scalars itself.
        """
        if value is None or value == []:
            return None

        refs = value if isinstance(value, list) else [value]
        resolved = [
            term for term in (
                self._resolve_controlled_term(ref, resolved_controlled_terms)
                for ref in refs
            )
            if term
        ]

        if not resolved:
            return None
        return resolved if isinstance(value, list) else resolved[0]

    @staticmethod
    def _dedupe_terms(values):
        """
        Flatten and de-duplicate resolved controlled terms by @id,
        while keeping them as dicts
        """
        flat = []
        for value in values:
            flat.extend(value if isinstance(value, list) else [value])

        by_id = {
            term['@id']: term
            for term in flat
            if isinstance(term, dict) and term.get('@id')
        }
        terms = sorted(by_id.values(), key=lambda term: term['@id'])

        if not terms:
            return None
        return terms[0] if len(terms) == 1 else terms


    def _join_unique(self, items):
        """Join unique non-empty items with semicolon"""
        # Flatten any list-valued items (e.g. array-typed fields) before stringifying
        flattened_items = []
        for item in items:
            if isinstance(item, list):
                flattened_items.extend(item)
            else:
                flattened_items.append(item)

        # Filter out empty/None values
        filtered_items = [str(item).strip() for item in flattened_items if item and str(item).strip()]
        
        if not filtered_items:
            return
        
        # Remove duplicates and sort
        unique_items = sorted(set(filtered_items))
        return '; '.join(unique_items)
    
    
    def _get_clean_alias(self, obj_or_list):
        """Extract alias from an object or list of objects, cleaning the result"""
        if not obj_or_list:
            return
        
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
        
        # If no aliases found, return
        if not aliases:
            return
        
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
        '--mode', '-m',
        help='mode/instance to run on',
        default='db2_demo'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output CSV filename (optional, auto-generates if not provided)'
    )
    
    args = parser.parse_args()
    
    try:
        connection = DB2lattice.Connection(args.mode)
        field_types, object_config = load_and_return_constant_dicts(args.mode)
        configs = Configs(
            FIELD_TYPES=field_types, 
            OBJECT_CONFIG=object_config
        )
        flattener = DB2Flattener(connection, configs)
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
