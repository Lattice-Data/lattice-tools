#!/usr/bin/env python3
"""
Analyze S3-to-local file mapping:
1. Validate uniqueness (each local file maps to unique S3 path)
2. Discover and validate mapping logic
3. Identify outliers

Usage:
    python analyze_mapping.py <mapping_file> <identifiers> [--mode MODE]
    
Modes:
    experiment (default): For files with experiment names like GENE7, CHEM16, REF3
    groupid: For files with group IDs like CUIMC_003_003
    
Examples:
    python analyze_mapping.py mapping.txt CHEM16
    python analyze_mapping.py mapping.txt GENE7,REF3
    python analyze_mapping.py mapping.txt CUIMC --mode groupid
"""

import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_mapping_file(filepath):
    """Parse the mapping file into list of (s3_path, local_path) tuples."""
    mappings = []
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = None
            # Try comma first (CSV format)
            if ',' in line:
                parts = line.split(',', 1)
            # Try tab
            elif '\t' in line:
                parts = line.split('\t')
            # Then try multiple spaces
            elif '  ' in line:
                parts = re.split(r'\s{2,}', line)
            # Finally try to split S3 path from local path (s3://... /path)
            elif line.startswith('s3://') and ' /' in line:
                # Find the transition from S3 to local path
                match = re.match(r'(s3://\S+)\s+(/\S+)', line)
                if match:
                    parts = [match.group(1), match.group(2)]
            
            if parts is None or len(parts) != 2:
                print(f"WARNING: Line {line_num} has {len(parts) if parts else 1} columns (expected 2): {line[:100]}...")
                continue
            mappings.append((parts[0].strip(), parts[1].strip(), line_num))
    return mappings


def extract_group_ids_from_mappings(mappings):
    """
    Extract all unique group IDs from S3 paths.
    
    Expected S3 path format:
    s3://czi-{provider}/{project}/{order}/{GroupID}/...
    """
    group_ids = set()
    s3_pattern = re.compile(r's3://[^/]+/[^/]+/[^/]+/([^/]+)/')
    
    for s3_path, _, _ in mappings:
        match = s3_pattern.match(s3_path)
        if match:
            group_id = match.group(1)
            # Filter out order-level files (e.g., "436665_LibraryInfo.xml")
            # Group IDs typically don't have underscores and aren't file-like
            if not ('_' in group_id and (group_id.endswith('.xml') or group_id.endswith('.json') or 
                                          group_id.endswith('.csv') or group_id.endswith('.txt'))):
                # Also filter out batch numbers like "436665_merged_trimmer-stats.csv"
                if not group_id[0].isdigit():
                    group_ids.add(group_id)
    
    return sorted(group_ids)


def check_uniqueness(mappings, outs_only=False):
    """Check that each local file maps to a unique S3 path.
    
    If outs_only=True, only check S3 paths containing '/outs/' for uniqueness.
    """
    print("=" * 80)
    print("1. UNIQUENESS VALIDATION")
    if outs_only:
        print("   (checking only S3 paths with '/outs/' per SOP)")
    print("=" * 80)
    
    # Filter mappings if outs_only
    if outs_only:
        filtered_mappings = [(s3, local, ln) for s3, local, ln in mappings if '/outs/' in s3]
        print(f"\nFiltered to {len(filtered_mappings)} mappings with '/outs/' in S3 path (out of {len(mappings)} total)")
    else:
        filtered_mappings = mappings
    
    # Check for duplicate local paths
    local_to_s3 = defaultdict(list)
    for s3_path, local_path, line_num in filtered_mappings:
        local_to_s3[local_path].append((s3_path, line_num))
    
    duplicate_locals = {k: v for k, v in local_to_s3.items() if len(v) > 1}
    
    if duplicate_locals:
        print(f"\nERROR: Found {len(duplicate_locals)} local files with multiple S3 mappings:")
        for local, s3_list in duplicate_locals.items():
            print(f"\n  Local: {local}")
            for s3, ln in s3_list:
                print(f"    Line {ln}: {s3}")
    else:
        print("\n✓ Each local file maps to exactly one S3 path")
    
    # Check for duplicate S3 paths
    s3_to_local = defaultdict(list)
    for s3_path, local_path, line_num in filtered_mappings:
        s3_to_local[s3_path].append((local_path, line_num))
    
    duplicate_s3 = {k: v for k, v in s3_to_local.items() if len(v) > 1}
    
    if duplicate_s3:
        print(f"\nERROR: Found {len(duplicate_s3)} S3 paths with multiple local mappings:")
        for s3, local_list in duplicate_s3.items():
            print(f"\n  S3: {s3}")
            for local, ln in local_list:
                print(f"    Line {ln}: {local}")
    else:
        print("✓ Each S3 path maps to exactly one local file")
    
    print(f"\nSummary:")
    print(f"  Total mappings checked: {len(filtered_mappings)}")
    print(f"  Unique local paths: {len(local_to_s3)}")
    print(f"  Unique S3 paths: {len(s3_to_local)}")
    
    return len(duplicate_locals) == 0 and len(duplicate_s3) == 0


def analyze_mapping_logic(mappings, experiments):
    """Analyze the transformation logic from local to S3 paths."""
    print("\n" + "=" * 80)
    print("2. MAPPING LOGIC ANALYSIS")
    print("=" * 80)
    
    # Build patterns dynamically based on experiment names (allow any of them)
    exp_pattern = '|'.join(re.escape(e) for e in experiments)
    
    # Flexible local pattern - capture everything before batch-date as base_path
    # Pattern for sample-level files
    local_sample_pattern = re.compile(
        r'(?P<base_path>.+?)/'
        r'(?P<batch>\d+)-(?P<date>\d+_\d+)/'
        rf'(?P<sample_id>\d+-(?P<exp>{exp_pattern})_(?P<plate>P\d+(?:_\d+)?)_(?P<well>[A-H]\d+)-(?P<zid>Z\d+)-(?P<barcode>[ACGT]+))/'
        r'(?P<filename>.+)'
    )
    
    # Pattern for batch-level files (no sample subdirectory)
    local_batch_pattern = re.compile(
        r'(?P<base_path>.+?)/'
        r'(?P<batch>\d+)-(?P<date>\d+_\d+)/'
        r'(?P<filename>[^/]+)$'
    )
    
    s3_sample_pattern = re.compile(
        rf's3://czi-novogene/trapnell-seahub-bcp/(?P<order>NVUS\d+(?:-\d+)+)/(?P<exp>{exp_pattern})/raw/(?P<batch>\d+)/'
        rf'(?P<sample_id>\d+-(?P<exp2>{exp_pattern})_(?P<plate>P\d+(?:_\d+)?)_(?P<well>[A-H]\d+)_GEX_hash_oligo-(?P<zid>Z\d+)-(?P<barcode>[ACGT]+))'
        r'(?P<suffix>.+)'
    )
    
    # Pattern for batch-level S3 files
    s3_batch_pattern = re.compile(
        rf's3://czi-novogene/trapnell-seahub-bcp/(?P<order>NVUS\d+(?:-\d+)+)/(?P<exp>{exp_pattern})/raw/(?P<batch>\d+)/'
        r'(?P<filename>[^/]+)$'
    )
    
    transformations = []
    batch_files = []
    outliers = []
    
    for s3_path, local_path, line_num in mappings:
        # Try sample-level pattern first (both local and S3 must match)
        local_sample_match = local_sample_pattern.match(local_path)
        s3_sample_match = s3_sample_pattern.match(s3_path)
        
        # Try batch-level pattern
        local_batch_match = local_batch_pattern.match(local_path)
        s3_batch_match = s3_batch_pattern.match(s3_path)
        
        if local_sample_match and s3_sample_match:
            # Sample-level file
            local_groups = local_sample_match.groupdict()
            s3_groups = s3_sample_match.groupdict()
            
            # Check field mappings
            issues = []
            
            if local_groups['batch'] != s3_groups['batch']:
                issues.append(f"batch mismatch: local={local_groups['batch']}, s3={s3_groups['batch']}")
            if local_groups['exp'] != s3_groups['exp']:
                issues.append(f"experiment mismatch: local={local_groups['exp']}, s3={s3_groups['exp']}")
            if local_groups['plate'] != s3_groups['plate']:
                issues.append(f"plate mismatch: local={local_groups['plate']}, s3={s3_groups['plate']}")
            if local_groups['well'] != s3_groups['well']:
                issues.append(f"well mismatch: local={local_groups['well']}, s3={s3_groups['well']}")
            if local_groups['zid'] != s3_groups['zid']:
                issues.append(f"zid mismatch: local={local_groups['zid']}, s3={s3_groups['zid']}")
            if local_groups['barcode'] != s3_groups['barcode']:
                issues.append(f"barcode mismatch: local={local_groups['barcode']}, s3={s3_groups['barcode']}")
            
            local_base = local_groups['filename']
            s3_suffix = s3_groups['suffix']
            expected_local_base = f"{local_groups['batch']}-{local_groups['exp']}_{local_groups['plate']}_{local_groups['well']}-{local_groups['zid']}-{local_groups['barcode']}{s3_suffix}"
            
            if local_base != expected_local_base:
                issues.append(f"filename transform unexpected: got '{local_base}', expected '{expected_local_base}'")
            
            if issues:
                outliers.append((line_num, 'field_mismatch', local_path, s3_path, issues))
            else:
                transformations.append({
                    'line': line_num,
                    'type': 'sample',
                    'base_path': local_groups['base_path'],
                    'batch': local_groups['batch'],
                    'exp': local_groups['exp'],
                    'plate': local_groups['plate'],
                    'well': local_groups['well'],
                    'zid': local_groups['zid'],
                    'barcode': local_groups['barcode'],
                    'date': local_groups['date'],
                    'order': s3_groups['order'],
                    'extension': Path(local_base).suffix if '.' in local_base else local_base.split('_')[-1]
                })
            continue
        else:
            # Try batch-level pattern
            local_batch_match = local_batch_pattern.match(local_path)
            s3_batch_match = s3_batch_pattern.match(s3_path)
            
            if local_batch_match and s3_batch_match:
                local_groups = local_batch_match.groupdict()
                s3_groups = s3_batch_match.groupdict()
                
                issues = []
                if local_groups['batch'] != s3_groups['batch']:
                    issues.append(f"batch mismatch: local={local_groups['batch']}, s3={s3_groups['batch']}")
                
                if issues:
                    outliers.append((line_num, 'batch_field_mismatch', local_path, s3_path, issues))
                else:
                    batch_files.append({
                        'line': line_num,
                        'type': 'batch',
                        'base_path': local_groups['base_path'],
                        'batch': local_groups['batch'],
                        'date': local_groups['date'],
                        'local_filename': local_groups['filename'],
                        's3_filename': s3_groups['filename'],
                        'order': s3_groups['order'],
                        'exp': s3_groups['exp']
                    })
            else:
                # Neither pattern matched
                outliers.append((line_num, 'parse_fail', local_path, s3_path))
    
    # Report findings
    print(f"\nSuccessfully parsed:")
    print(f"  Sample-level files: {len(transformations)}")
    print(f"  Batch-level files: {len(batch_files)}")
    print(f"  Total: {len(transformations) + len(batch_files)} / {len(mappings)} mappings")
    
    if transformations:
        # Analyze consistent patterns for sample files
        base_paths = set(t['base_path'] for t in transformations)
        orders = set(t['order'] for t in transformations)
        exps = set(t['exp'] for t in transformations)
        plates = set(t['plate'] for t in transformations)
        batches = set(t['batch'] for t in transformations)
        dates = set(t['date'] for t in transformations)
        
        print(f"\nSample-level file fields:")
        print(f"  Local base paths: {sorted(base_paths)}")
        print(f"  Order IDs: {sorted(orders)}")
        print(f"  Experiments: {sorted(exps)}")
        print(f"  Plates: {sorted(plates)}")
        print(f"  Batch IDs: {sorted(batches)}")
        print(f"  Dates: {sorted(dates)}")
        
        # File types
        extensions = defaultdict(int)
        for t in transformations:
            ext = t['extension']
            extensions[ext] += 1
        
        print(f"\n  File types:")
        for ext, count in sorted(extensions.items(), key=lambda x: -x[1]):
            print(f"    {ext}: {count}")
    
    if batch_files:
        # Analyze batch-level files
        batch_filenames = defaultdict(int)
        for bf in batch_files:
            batch_filenames[bf['local_filename']] += 1
        
        print(f"\nBatch-level file types:")
        for fname, count in sorted(batch_filenames.items(), key=lambda x: -x[1]):
            print(f"    {fname}: {count}")
    
    print(f"\n" + "-" * 40)
    print("MAPPING LOGIC DISCOVERED:")
    print("-" * 40)
    exp_str = '|'.join(experiments)
    print(f"""
The transformation from local to S3 follows this pattern:

LOCAL:  {{base_path}}/{{batch}}-{{date}}/{{batch}}-{{exp}}_{{plate}}_{{well}}-{{zid}}-{{barcode}}/{{filename}}
S3:     s3://czi-novogene/trapnell-seahub-bcp/{{order}}/{{exp}}/raw/{{batch}}/{{batch}}-{{exp}}_{{plate}}_{{well}}_GEX_hash_oligo-{{zid}}-{{barcode}}{{suffix}}

Where {{exp}} is one of: {exp_str}

Key transformations:
1. Static S3 prefix: s3://czi-novogene/trapnell-seahub-bcp/
2. Order ID is added to S3 path (not present in local)
3. Local base path and date are dropped in S3 path
4. Directory structure flattens: local nested dirs → single S3 directory per batch
5. Sample ID transformation: 
   - Local: {{batch}}-{{exp}}_{{plate}}_{{well}}-{{zid}}-{{barcode}}
   - S3:    {{batch}}-{{exp}}_{{plate}}_{{well}}_GEX_hash_oligo-{{zid}}-{{barcode}}
   - The hyphen between well and Z-ID becomes "_GEX_hash_oligo-"
""")
    
    return outliers


def report_outliers(outliers):
    """Report any mapping outliers found."""
    print("\n" + "=" * 80)
    print("3. OUTLIERS")
    print("=" * 80)
    
    if not outliers:
        print("\n✓ No outliers found - all mappings follow the same logic")
        return
    
    print(f"\nFound {len(outliers)} outliers:")
    
    # Group by type
    by_type = defaultdict(list)
    for o in outliers:
        by_type[o[1]].append(o)
    
    for outlier_type, items in by_type.items():
        print(f"\n{outlier_type}: {len(items)} cases")
        for item in items:  # Show all outliers
            print(f"  Line {item[0]}:")
            print(f"    Local: {item[2]}")
            print(f"    S3: {item[3]}")
            if len(item) > 4:
                for issue in item[4]:
                    print(f"    Issue: {issue}")


def analyze_mapping_logic_groupid(mappings, groupids):
    """Analyze the transformation logic for groupid-based mappings with SOP compliance."""
    print("\n" + "=" * 80)
    print("2. MAPPING LOGIC ANALYSIS (GroupID mode with SOP compliance)")
    print("=" * 80)
    
    # Base suffixes common to all assays
    BASE_SUFFIXES = [
        '.csv',                          # General run stats - CSV
        '.json',                         # General run stats - JSON
        '_trimmer-failure_codes.csv',    # Trimmer failure codes
        '_trimmer-stats.csv',            # Trimmer stats
        '_unmatched.cram',               # Unmatched reads CRAM
        '_unmatched.csv',                # Unmatched reads CSV
        '_unmatched.json',               # Unmatched reads JSON
    ]
    
    # FASTQ-related suffixes per read type
    def get_read_suffixes(read_type):
        """Get suffixes for a given read type (R1, R2, I1, I2, R3)"""
        return [
            f'_S1_L001_{read_type}_001.csv',
            f'_S1_L001_{read_type}_001.fastq.gz',
            f'_S1_L001_{read_type}_001.json',
            # _sample.fastq.gz variants are optional - not required
        ]
    
    # scRNA/GEX specific suffixes - none are required
    SCRNA_SUFFIXES = []

    # Optional suffixes - valid if present but not required (presence/absence not flagged)
    OPTIONAL_SUFFIXES = [
        '.scRNA.applicationQC.h5',       # scRNA specific QC - h5
        '.scRNA.applicationQC.html',     # scRNA specific QC - html
        '_Log.final.out',                # STAR Solo log
        '_Log.out',                      # STAR Solo log
        '_Log.progress.out',             # STAR Solo log
        '_ReadsPerGene.out.tab',         # STAR Solo gene count
        '_SJ.out.tab',                   # STAR solo read counts
        '_S1_L001_R1_001_sample.fastq.gz',  # sample fastq R1
        '_S1_L001_R2_001_sample.fastq.gz',  # sample fastq R2
    ]
    
    # Assay-specific expected read types
    ASSAY_REQUIRED_READS = {
        'GEX': ['R1', 'R2'],
        'CRI': ['R1', 'R2'],
        'ATAC': ['R1', 'R2', 'I2'],
    }
    
    # Unexpected read types for any assay (flag these)
    UNEXPECTED_READS = ['I1', 'R3']
    
    # Build expected suffixes per assay
    def get_expected_suffixes_for_assay(assay):
        expected = BASE_SUFFIXES.copy()
        if assay in ['GEX', 'CRI']:
            expected.extend(SCRNA_SUFFIXES)
        required_reads = ASSAY_REQUIRED_READS.get(assay, ['R1', 'R2'])
        for read_type in required_reads:
            expected.extend(get_read_suffixes(read_type))
        return set(expected)
    
    # All potentially valid suffixes (for pattern matching)
    ALL_VALID_SUFFIXES = set(BASE_SUFFIXES + SCRNA_SUFFIXES + OPTIONAL_SUFFIXES)
    for read_type in ['R1', 'R2', 'I1', 'I2', 'R3']:
        ALL_VALID_SUFFIXES.update(get_read_suffixes(read_type))
    
    # Valid assay types per SOP
    VALID_ASSAYS = ['GEX', 'CRI', 'ATAC']
    
    # More permissive pattern for catching mismatches
    # GroupID can be simple (CH05) or compound (CUIMC_003_003)
    groupid_permissive = r'[A-Za-z0-9]+(?:_[A-Za-z0-9]+)*'
    # For local paths, groupid may have assay concatenated without underscore (CH05GEX)
    groupid_with_assay = r'[A-Za-z0-9]+(?:_[A-Za-z0-9]+)*'
    
    # SOP S3 path: s3://czi-{provider}/{project}/{order}/{GroupID}/raw/{RunID}-{GroupID}_{Assay}-{UG-BC}{suffix}
    # Local can have two formats:
    #   1. {base_path}/{batch}-{date}/{batch}-{groupid}-{zid}-{barcode}/{filename}  (groupid includes assay without underscore)
    #   2. {base_path}/{batch}-{date}/{batch}-{groupid}_{assay}-{zid}-{barcode}/{filename}  (with underscore before assay)
    
    # Local pattern - groupid may include assay concatenated (CH05GEX) or with underscore (CUIMC_003_003)
    # Format 1: directory includes group ID: {batch}-{groupid}-{zid}-{barcode}/
    local_sample_pattern = re.compile(
        r'(?P<base_path>.+?)/'
        r'(?P<batch>\d+)-(?P<date>\d+_\d+)/'
        rf'(?P<sample_id>\d+-(?P<groupid>{groupid_with_assay})-(?P<zid>Z\d+)-(?P<barcode>[ACGT]+))/'
        r'(?P<filename>.+)'
    )
    
    # Format 2: directory omits group ID: {batch}-{zid}-{barcode}/
    # Group ID is still present in the filename itself: {batch}-{groupid}-{zid}-{barcode}...
    local_sample_pattern_nogroupid = re.compile(
        r'(?P<base_path>.+?)/'
        r'(?P<batch>\d+)-(?P<date>\d+_\d+)/'
        rf'(?P<batch2>\d+)-(?P<zid>Z\d+)-(?P<barcode>[ACGT]+)/'
        rf'(?P<batch3>\d+)-(?P<groupid_in_file>{groupid_with_assay})-(?P=zid)-(?P=barcode)(?P<filename_suffix>.*)'
    )
    
    local_batch_pattern = re.compile(
        r'(?P<base_path>.+?)/'
        r'(?P<batch>\d+)-(?P<date>\d+_\d+)/'
        r'(?P<filename>[^/]+)$'
    )
    
    # SOP-compliant S3 pattern
    # Path groupid is just the base (CH05), filename has groupid_assay (CH05_ATAC)
    s3_sample_pattern = re.compile(
        rf's3://(?P<bucket>czi-[^/]+)/(?P<project>[^/]+)/(?P<order>NVUS\d+(?:-\d+)+|AN\d+)/(?P<groupid>{groupid_permissive})/raw/'
        rf'(?P<runid>\d+)-(?P<groupid2>{groupid_permissive})_(?P<assay>[A-Z]+)-(?P<zid>Z\d+)-(?P<barcode>[ACGT]+)'
        r'(?P<suffix>.+)'
    )
    
    s3_batch_pattern = re.compile(
        rf's3://(?P<bucket>czi-[^/]+)/(?P<project>[^/]+)/(?P<order>NVUS\d+(?:-\d+)+|AN\d+)/(?P<groupid>{groupid_permissive})/raw/'
        r'(?P<filename>[^/]+)$'
    )
    
    # Also try matching batch-level files at order level (no groupid in path)
    s3_order_batch_pattern = re.compile(
        rf's3://(?P<bucket>czi-[^/]+)/(?P<project>[^/]+)/(?P<order>NVUS\d+(?:-\d+)+|AN\d+)/'
        r'(?P<filename>[^/]+)$'
    )
    
    # Keep track of expected groupids for validation
    expected_groupids = set(groupids)
    
    transformations = []
    batch_files = []
    outliers = []
    sop_warnings = []
    
    # Track files per sample for completeness check
    files_per_sample = defaultdict(set)
    
    for s3_path, local_path, line_num in mappings:
        local_sample_match = local_sample_pattern.match(local_path)
        s3_sample_match = s3_sample_pattern.match(s3_path)
        local_batch_match = local_batch_pattern.match(local_path)
        s3_batch_match = s3_batch_pattern.match(s3_path)
        
        # Fallback: if local dir has no group ID but S3 looks like a sample path,
        # try the pattern that extracts group ID from the filename instead
        local_nogroupid_match = None
        if not local_sample_match and s3_sample_match:
            local_nogroupid_match = local_sample_pattern_nogroupid.match(local_path)
            if local_nogroupid_match:
                # Normalize into the same shape as local_sample_match groupdict
                # so the rest of the validation code can run unchanged
                ng = local_nogroupid_match.groupdict()
                # Reconstruct a fake groupdict compatible with local_sample_pattern
                normalized = {
                    'base_path': ng['base_path'],
                    'batch':     ng['batch'],
                    'date':      ng['date'],
                    'groupid':   ng['groupid_in_file'],
                    'zid':       ng['zid'],
                    'barcode':   ng['barcode'],
                    # Reconstruct filename as the full token starting from {batch}-{groupid}-...
                    'filename':  f"{ng['batch3']}-{ng['groupid_in_file']}-{ng['zid']}-{ng['barcode']}{ng['filename_suffix']}",
                    'sample_id': f"{ng['batch3']}-{ng['groupid_in_file']}-{ng['zid']}-{ng['barcode']}",
                    '_dir_had_groupid': False,
                }
                # Use this as local_sample_match by wrapping in a simple namespace
                class _FakeMatch:
                    def __init__(self, d): self._d = d
                    def groupdict(self): return self._d
                local_sample_match = _FakeMatch(normalized)
        
        if local_sample_match and s3_sample_match:
            local_groups = local_sample_match.groupdict()
            s3_groups = s3_sample_match.groupdict()
            
            issues = []
            warnings = []
            
            # Local groupid might be in different formats:
            # 1. Same as S3 path groupid (e.g., CUIMC_003_003)
            # 2. S3 path groupid + assay concatenated (e.g., CH05GEX from CH05 + GEX)
            # 3. Partial component of S3 path groupid (e.g., MS116D from MS116D_MS116DF)
            local_gid = local_groups['groupid']
            s3_path_gid = s3_groups['groupid']
            s3_file_gid = s3_groups['groupid2']
            assay = s3_groups['assay']
            
            # Check if local groupid matches expected patterns
            local_matches_path = local_gid == s3_path_gid
            local_matches_path_plus_assay = local_gid == s3_path_gid + assay
            local_matches_file = local_gid == s3_file_gid
            local_matches_file_plus_assay = local_gid == s3_file_gid + assay
            
            # Validate S3 path groupid is in expected list
            if s3_path_gid not in expected_groupids:
                issues.append(f"unexpected s3 path groupid: {s3_path_gid} (expected one of: {', '.join(sorted(expected_groupids))})")
            
            # Validate S3 filename groupid matches path groupid
            if s3_path_gid != s3_file_gid:
                issues.append(f"groupid mismatch (s3 internal): s3_path={s3_path_gid}, s3_filename={s3_file_gid}")
            
            # Validate local groupid relationship
            if not (local_matches_path or local_matches_path_plus_assay or local_matches_file or local_matches_file_plus_assay):
                # Check if local_gid is a valid partial component of the s3 groupid
                # For compound IDs like MS116D_MS116DF, local can be MS116D or MS116DF
                s3_gid_parts = s3_path_gid.split('_')
                local_is_valid_partial = local_gid in s3_gid_parts
                
                # Also check if any expected groupid starts with local_gid (original behavior)
                matching_expected = [g for g in expected_groupids if local_gid.startswith(g)]
                
                # Also check if any expected groupid contains local_gid as a component
                matching_component = [g for g in expected_groupids if local_gid in g.split('_')]
                
                if not (local_is_valid_partial or matching_expected or matching_component):
                    issues.append(f"local groupid '{local_gid}' doesn't match s3 groupid '{s3_path_gid}' or '{s3_path_gid}{assay}'")
            
            # Validate batch/RunID
            if local_groups['batch'] != s3_groups['runid']:
                issues.append(f"batch/RunID mismatch: local={local_groups['batch']}, s3={s3_groups['runid']}")
            
            # Validate Z-ID and barcode
            if local_groups['zid'] != s3_groups['zid']:
                issues.append(f"zid mismatch: local={local_groups['zid']}, s3={s3_groups['zid']}")
            if local_groups['barcode'] != s3_groups['barcode']:
                issues.append(f"barcode mismatch: local={local_groups['barcode']}, s3={s3_groups['barcode']}")
            
            # SOP compliance checks
            # Check assay type
            if s3_groups['assay'] not in VALID_ASSAYS:
                warnings.append(f"SOP: invalid assay type '{s3_groups['assay']}' (expected one of: {', '.join(VALID_ASSAYS)})")
            
            # Check bucket naming convention (should be czi-{provider})
            if not s3_groups['bucket'].startswith('czi-'):
                warnings.append(f"SOP: bucket should start with 'czi-', got '{s3_groups['bucket']}'")
            
            # Check Z-ID format (should be Z followed by 4 digits)
            if not re.match(r'Z\d{4}$', s3_groups['zid']):
                warnings.append(f"SOP: Z-ID should be Z followed by 4 digits, got '{s3_groups['zid']}'")
            
            # Track file suffix and assay for completeness check
            sample_key = (s3_groups['runid'], s3_groups['groupid2'], s3_groups['zid'], s3_groups['barcode'], assay)
            files_per_sample[sample_key].add(s3_groups['suffix'])
            
            if issues:
                outliers.append((line_num, 'field_mismatch', local_path, s3_path, issues))
            else:
                if warnings:
                    sop_warnings.append((line_num, local_path, s3_path, warnings))
                
                local_base = local_groups['filename']
                transformations.append({
                    'line': line_num,
                    'type': 'sample',
                    'base_path': local_groups['base_path'],
                    'batch': local_groups['batch'],
                    'groupid': local_groups['groupid'],
                    'zid': local_groups['zid'],
                    'barcode': local_groups['barcode'],
                    'date': local_groups['date'],
                    'order': s3_groups['order'],
                    'bucket': s3_groups['bucket'],
                    'project': s3_groups['project'],
                    'assay': s3_groups['assay'],
                    'suffix': s3_groups['suffix'],
                    'extension': Path(local_base).suffix if '.' in local_base else local_base.split('_')[-1]
                })
            continue
        
        if local_batch_match and s3_batch_match:
            local_groups = local_batch_match.groupdict()
            s3_groups = s3_batch_match.groupdict()
            
            issues = []
            if issues:
                outliers.append((line_num, 'batch_field_mismatch', local_path, s3_path, issues))
            else:
                batch_files.append({
                    'line': line_num,
                    'type': 'batch',
                    'base_path': local_groups['base_path'],
                    'batch': local_groups['batch'],
                    'date': local_groups['date'],
                    'local_filename': local_groups['filename'],
                    's3_filename': s3_groups['filename'],
                    'order': s3_groups['order'],
                    'groupid': s3_groups.get('groupid', 'N/A'),
                    'bucket': s3_groups['bucket']
                })
        elif local_batch_match and s3_order_batch_pattern.match(s3_path):
            # Order-level batch file (no groupid in S3 path)
            local_groups = local_batch_match.groupdict()
            s3_match = s3_order_batch_pattern.match(s3_path)
            s3_groups = s3_match.groupdict()
            
            batch_files.append({
                'line': line_num,
                'type': 'order_batch',
                'base_path': local_groups['base_path'],
                'batch': local_groups['batch'],
                'date': local_groups['date'],
                'local_filename': local_groups['filename'],
                's3_filename': s3_groups['filename'],
                'order': s3_groups['order'],
                'groupid': 'N/A',
                'bucket': s3_groups['bucket']
            })
        else:
            # Provide more detail about why parsing failed
            parse_issues = []
            if not local_sample_match and not local_batch_match:
                parse_issues.append("local path doesn't match expected pattern")
            if not s3_sample_match and not s3_batch_match:
                parse_issues.append("S3 path doesn't match expected pattern")
            if local_sample_match and not s3_sample_match:
                parse_issues.append("local looks like sample but S3 doesn't match sample pattern")
            if s3_sample_match and not local_sample_match:
                parse_issues.append("S3 looks like sample but local doesn't match sample pattern")
            outliers.append((line_num, 'parse_fail', local_path, s3_path, parse_issues if parse_issues else None))
    
    # Report findings
    print(f"\nSuccessfully parsed:")
    print(f"  Sample-level files: {len(transformations)}")
    print(f"  Batch-level files: {len(batch_files)}")
    print(f"  Total: {len(transformations) + len(batch_files)} / {len(mappings)} mappings")
    
    if transformations:
        base_paths = set(t['base_path'] for t in transformations)
        buckets = set(t['bucket'] for t in transformations)
        projects = set(t['project'] for t in transformations)
        orders = set(t['order'] for t in transformations)
        groupids_found = set(t['groupid'] for t in transformations)
        batches = set(t['batch'] for t in transformations)
        dates = set(t['date'] for t in transformations)
        assays = set(t['assay'] for t in transformations)
        
        print(f"\nSample-level file fields:")
        print(f"  Local base paths: {sorted(base_paths)}")
        print(f"  S3 buckets: {sorted(buckets)}")
        print(f"  S3 projects: {sorted(projects)}")
        print(f"  Order IDs: {sorted(orders)}")
        print(f"  Group IDs: {sorted(groupids_found)}")
        print(f"  Batch/Run IDs: {sorted(batches)}")
        print(f"  Dates: {sorted(dates)}")
        print(f"  Assay types: {sorted(assays)}")
        
        # Analyze group ID transformation patterns (covers both dir formats)
        print(f"\n  Group ID Transformations:")
        groupid_to_s3 = {}  # local groupid -> (s3 groupid, assay)
        for s3_path, local_path, line_num in mappings:
            s3_match = s3_sample_pattern.match(s3_path)
            if not s3_match:
                continue
            local_match = local_sample_pattern.match(local_path)
            if local_match:
                local_gid = local_match.group('groupid')
            else:
                local_match_ng = local_sample_pattern_nogroupid.match(local_path)
                if local_match_ng:
                    local_gid = local_match_ng.group('groupid_in_file')
                else:
                    continue
            s3_gid = s3_match.group('groupid')
            assay = s3_match.group('assay')
            if local_gid not in groupid_to_s3:
                groupid_to_s3[local_gid] = (s3_gid, assay)
        
        # Categorize transformation patterns
        pattern_concat = []  # local = s3 + assay (e.g., CH05GEX → CH05 + GEX)
        pattern_exact = []   # local = s3 (e.g., CUIMC_003_003 → CUIMC_003_003)
        pattern_other = []   # Other patterns
        
        for local_gid, (s3_gid, assay) in sorted(groupid_to_s3.items()):
            if local_gid == s3_gid + assay:
                pattern_concat.append(f"{local_gid} → {s3_gid} + {assay}")
            elif local_gid == s3_gid:
                pattern_exact.append(f"{local_gid} → {s3_gid}")
            else:
                pattern_other.append(f"{local_gid} → {s3_gid} (+ {assay})")
        
        if pattern_concat:
            print(f"    Local = S3_GroupID + Assay (concatenated without separator):")
            for p in pattern_concat:
                print(f"      {p}")
        
        if pattern_exact:
            print(f"    Local = S3_GroupID (exact match):")
            for p in pattern_exact:
                print(f"      {p}")
        
        if pattern_other:
            print(f"    Other patterns:")
            for p in pattern_other:
                print(f"      {p}")
        
        # Show unique S3 group IDs extracted
        s3_groupids = set(s3_gid for s3_gid, _ in groupid_to_s3.values())
        print(f"\n  S3 Group IDs (extracted from paths): {sorted(s3_groupids)}")
        print(f"  Local Group IDs (from directories): {sorted(groupid_to_s3.keys())}")
        
        # File suffix analysis
        suffixes = defaultdict(int)
        for t in transformations:
            suffixes[t['suffix']] += 1
        
        print(f"\n  File suffixes:")
        for suffix, count in sorted(suffixes.items(), key=lambda x: -x[1]):
            # Check if suffix is valid for any assay
            if suffix in ALL_VALID_SUFFIXES:
                # Check if it's an unexpected read type
                is_unexpected = any(f'_{rt}_' in suffix for rt in UNEXPECTED_READS)
                if is_unexpected:
                    print(f"    ⚠ {suffix}: {count} (unexpected read type)")
                else:
                    print(f"    ✓ {suffix}: {count}")
            else:
                print(f"    ? {suffix}: {count}")
    
    # SOP completeness check
    print(f"\n" + "-" * 40)
    print("SOP COMPLIANCE CHECK")
    print("-" * 40)
    
    if files_per_sample:
        print(f"\nSamples found: {len(files_per_sample)}")
        
        # Group by assay for reporting
        samples_by_assay = defaultdict(list)
        for sample_key in files_per_sample.keys():
            runid, groupid, zid, barcode, assay = sample_key
            samples_by_assay[assay].append(sample_key)
        
        print(f"\nSamples by assay type:")
        for assay, samples in sorted(samples_by_assay.items()):
            required_reads = ASSAY_REQUIRED_READS.get(assay, ['R1', 'R2'])
            print(f"  {assay}: {len(samples)} samples (required reads: {', '.join(required_reads)})")
        
        incomplete_samples = []
        unexpected_read_samples = []
        extra_files = []
        
        for sample_key, found_suffixes in files_per_sample.items():
            runid, groupid, zid, barcode, assay = sample_key
            sample_name = f"{runid}-{groupid}_{assay}_{zid}-{barcode}"
            
            # Get expected suffixes for this assay
            expected_suffixes = get_expected_suffixes_for_assay(assay)
            
            missing = expected_suffixes - found_suffixes
            extra = found_suffixes - expected_suffixes
            
            # Check for unexpected read types (I1, R3)
            unexpected_reads_found = []
            for suffix in found_suffixes:
                for rt in UNEXPECTED_READS:
                    if f'_{rt}_' in suffix:
                        unexpected_reads_found.append(suffix)
                        break
            
            if missing:
                incomplete_samples.append((sample_name, assay, missing))
            if unexpected_reads_found:
                unexpected_read_samples.append((sample_name, assay, unexpected_reads_found))
            # Extra files exclude unexpected reads (reported separately)
            extra_non_unexpected = extra - set(unexpected_reads_found) - set(OPTIONAL_SUFFIXES)
            if extra_non_unexpected:
                extra_files.append((sample_name, assay, extra_non_unexpected))
        
        if incomplete_samples:
            print(f"\n⚠ Samples missing expected files: {len(incomplete_samples)}")
            for sample_name, assay, missing in incomplete_samples:
                print(f"  {sample_name} ({assay}):")
                for m in sorted(missing):
                    print(f"    missing: {m}")
        else:
            print(f"\n✓ All samples have their expected file types for their assay")
        
        if unexpected_read_samples:
            print(f"\n⚠ Samples with unexpected read types (I1, R3): {len(unexpected_read_samples)}")
            for sample_name, assay, unexpected in unexpected_read_samples:
                print(f"  {sample_name} ({assay}):")
                for u in sorted(unexpected):
                    print(f"    unexpected: {u}")
        
        if extra_files:
            print(f"\nℹ Samples with extra files (not in SOP for their assay): {len(extra_files)}")
            for sample_name, assay, extra in extra_files:
                print(f"  {sample_name} ({assay}):")
                for e in sorted(extra):
                    print(f"    extra: {e}")
    
    # Report SOP warnings
    if sop_warnings:
        print(f"\n⚠ SOP warnings: {len(sop_warnings)}")
        for line_num, local_path, s3_path, warnings in sop_warnings:
            print(f"  Line {line_num}:")
            for w in warnings:
                print(f"    {w}")
    
    if batch_files:
        batch_filenames = defaultdict(int)
        for bf in batch_files:
            batch_filenames[bf['local_filename']] += 1
        
        print(f"\nBatch-level file types:")
        for fname, count in sorted(batch_filenames.items(), key=lambda x: -x[1]):
            print(f"    {fname}: {count}")
    
    print(f"\n" + "-" * 40)
    print("MAPPING LOGIC DISCOVERED:")
    print("-" * 40)
    groupid_str = ', '.join(sorted(groupids))
    print(f"""
The transformation from local to S3 follows this pattern:

LOCAL:  {{base_path}}/{{RunID}}-{{date}}/{{RunID}}-{{GroupID}}-{{UG-BC}}/{{filename}}
S3:     s3://czi-{{provider}}/{{project}}/{{order}}/{{GroupID}}/raw/{{RunID}}-{{GroupID}}_{{Assay}}-{{UG-BC}}{{suffix}}

Expected Group IDs: {groupid_str}
Valid Assay types: {', '.join(VALID_ASSAYS)}

Key transformations:
1. S3 bucket format: czi-{{provider}} (e.g., czi-novogene)
2. Order ID added to S3 path
3. Local base path and date dropped in S3 path
4. Sample ID transformation: 
   - Local: {{RunID}}-{{GroupID}}-{{UG-BC}}
   - S3:    {{RunID}}-{{GroupID}}_{{Assay}}-{{UG-BC}}
   - Hyphen before UG-BC becomes "_{{Assay}}-"
5. Local GroupID can be partial (e.g., MS116D from MS116D_MS116DF)
""")
    
    return outliers


def analyze_mapping_logic_processed(mappings, sample_ids=None):
    """
    Analyze the transformation logic for processed files mappings per SOP.
    
    Args:
        mappings: List of (s3_path, local_path, line_num) tuples
        sample_ids: Optional list of sample IDs to validate. If None, auto-extract from S3 paths.
    """
    print("\n" + "=" * 80)
    print("2. MAPPING LOGIC ANALYSIS (Processed files mode)")
    print("=" * 80)
    
    # Auto-extract group IDs if not provided
    if sample_ids is None:
        sample_ids = extract_group_ids_from_mappings(mappings)
        print(f"\nAuto-detected {len(sample_ids)} Group IDs from S3 paths: {', '.join(sample_ids)}")
    else:
        print(f"\nUsing {len(sample_ids)} specified Group IDs: {', '.join(sample_ids)}")
    
    # Valid pipeline types per SOP
    VALID_PIPELINES = ['cellranger', 'scalerna']
    
    # Key Cell Ranger output files/directories that should be present in outs/
    CELLRANGER_KEY_OUTPUTS = [
        'outs/web_summary.html',
        'outs/metrics_summary.csv',
        'outs/cloupe.cloupe',
        'outs/filtered_feature_bc_matrix.h5',
        'outs/raw_feature_bc_matrix.h5',
        'outs/molecule_info.h5',
        'outs/possorted_genome_bam.bam',
        'outs/possorted_genome_bam.bam.bai',
        'outs/possorted_genome_bam.bam.csi',
    ]
    
    # Multi-sample Cell Ranger outputs (outs/multi/count/)
    CELLRANGER_MULTI_KEY_OUTPUTS = [
        'outs/multi/count/raw_feature_bc_matrix.h5',
        'outs/multi/count/filtered_feature_bc_matrix.h5',
    ]
    
    # S3 pattern for processed files:
    # s3://czi-{provider}/{project}/{order}/{GroupID}/processed/{pipeline}/Run_{date}/{rest_of_path}
    s3_processed_pattern = re.compile(
        rf's3://(?P<bucket>czi-[^/]+)/(?P<project>[^/]+)/(?P<order>NVUS\d+(?:-\d+)+|AN\d+)/'
        rf'(?P<groupid>[^/]+)/processed/(?P<pipeline>[^/]+)/(?P<run_id>Run_\d{{4}}-\d{{2}}-\d{{2}})'
        r'(?P<rest_of_path>/.*)?'
    )
    
    transformations = []
    outliers = []
    sop_warnings = []
    
    # Track files per sample
    files_per_sample = defaultdict(list)
    outs_files_per_sample = defaultdict(set)
    pipelines_found = set()
    run_ids_found = set()
    
    expected_groupids = set(sample_ids)
    
    for s3_path, local_path, line_num in mappings:
        s3_match = s3_processed_pattern.match(s3_path)
        
        if s3_match:
            s3_groups = s3_match.groupdict()
            
            issues = []
            warnings = []
            
            # Validate GroupID is in expected list
            if s3_groups['groupid'] not in expected_groupids:
                issues.append(f"unexpected GroupID: {s3_groups['groupid']} (expected one of: {', '.join(sorted(expected_groupids))})")
            
            # Check pipeline type per SOP
            if s3_groups['pipeline'] not in VALID_PIPELINES:
                warnings.append(f"SOP: non-standard pipeline '{s3_groups['pipeline']}' (expected: {', '.join(VALID_PIPELINES)})")
            
            # Validate Run date format (Run_YYYY-MM-DD)
            run_id = s3_groups['run_id']
            if not re.match(r'Run_\d{4}-\d{2}-\d{2}$', run_id):
                warnings.append(f"SOP: run_id format should be Run_YYYY-MM-DD, got '{run_id}'")
            
            # Track for reporting
            pipelines_found.add(s3_groups['pipeline'])
            run_ids_found.add(s3_groups['run_id'])
            
            sample_key = (s3_groups['groupid'], s3_groups['pipeline'], s3_groups['run_id'])
            rest_of_path = s3_groups.get('rest_of_path', '') or ''
            files_per_sample[sample_key].append(rest_of_path)
            
            # Track outs/ files specifically
            if rest_of_path.startswith('/outs/') or rest_of_path == '/outs':
                # Normalize path for comparison
                outs_path = rest_of_path[1:] if rest_of_path.startswith('/') else rest_of_path
                outs_files_per_sample[sample_key].add(outs_path)
            
            if issues:
                outliers.append((line_num, 'field_mismatch', local_path, s3_path, issues))
            else:
                if warnings:
                    sop_warnings.append((line_num, s3_path, warnings))
                
                transformations.append({
                    'line': line_num,
                    'bucket': s3_groups['bucket'],
                    'project': s3_groups['project'],
                    'order': s3_groups['order'],
                    'groupid': s3_groups['groupid'],
                    'pipeline': s3_groups['pipeline'],
                    'run_id': s3_groups['run_id'],
                    'rest_of_path': rest_of_path,
                    'local_path': local_path,
                })
        else:
            # Provide detailed explanation of what doesn't match
            detailed_issues = []
            
            # Check if it starts with s3://
            if not s3_path.startswith('s3://'):
                detailed_issues.append("Path must start with 's3://'")
            else:
                # Parse what we can find
                parts = s3_path.replace('s3://', '').split('/')
                
                expected_format = "s3://czi-{provider}/{project}/{order}/{GroupID}/processed/{pipeline}/Run_{date}/..."
                detailed_issues.append(f"Expected format: {expected_format}")
                detailed_issues.append(f"Got: {s3_path}")
                
                # Check each component
                if len(parts) < 6:
                    detailed_issues.append(f"Too few path components (found {len(parts)}, need at least 6)")
                    detailed_issues.append(f"  Components found: {' / '.join(parts[:6])}")
                else:
                    bucket, project, order, groupid, subdir, pipeline = parts[:6]
                    
                    # Check bucket format
                    if not bucket.startswith('czi-'):
                        detailed_issues.append(f"Bucket must start with 'czi-', got: '{bucket}'")
                    
                    # Check order format
                    if not re.match(r'(NVUS\d+(?:-\d+)+|AN\d+)$', order):
                        detailed_issues.append(f"Order must be NVUS####-##[-##...] or AN####, got: '{order}'")
                    
                    # Check for 'processed' subdirectory
                    if subdir != 'processed':
                        detailed_issues.append(f"5th component must be 'processed', got: '{subdir}'")
                    
                    # Check pipeline
                    if pipeline not in VALID_PIPELINES:
                        detailed_issues.append(f"Pipeline must be one of {VALID_PIPELINES}, got: '{pipeline}'")
                    
                    # Check Run_date format
                    if len(parts) >= 7:
                        run_date = parts[6]
                        if not re.match(r'Run_\d{4}-\d{2}-\d{2}$', run_date):
                            detailed_issues.append(f"Run date must be 'Run_YYYY-MM-DD', got: '{run_date}'")
            
            outliers.append((line_num, 'parse_fail', local_path, s3_path, detailed_issues))
    
    # Report findings
    print(f"\nSuccessfully parsed:")
    print(f"  Processed files: {len(transformations)}")
    print(f"  Total: {len(transformations)} / {len(mappings)} mappings")
    
    if transformations:
        buckets = set(t['bucket'] for t in transformations)
        projects = set(t['project'] for t in transformations)
        orders = set(t['order'] for t in transformations)
        groupids = set(t['groupid'] for t in transformations)
        
        print(f"\nProcessed file fields:")
        print(f"  S3 buckets: {sorted(buckets)}")
        print(f"  S3 projects: {sorted(projects)}")
        print(f"  Order IDs: {sorted(orders)}")
        print(f"  Group IDs: {sorted(groupids)}")
        print(f"  Pipelines: {sorted(pipelines_found)}")
        print(f"  Run IDs: {sorted(run_ids_found)}")
    
    # SOP Compliance Check
    print(f"\n" + "-" * 40)
    print("SOP COMPLIANCE CHECK")
    print("-" * 40)
    
    if files_per_sample:
        print(f"\nSamples with processed files: {len(files_per_sample)}")
        
        for sample_key in sorted(files_per_sample.keys()):
            groupid, pipeline, run_id = sample_key
            files = files_per_sample[sample_key]
            outs_files = outs_files_per_sample.get(sample_key, set())
            
            print(f"\n  {groupid} ({pipeline}, {run_id}):")
            print(f"    Total files: {len(files)}")
            
            # Count files by category
            outs_count = 0
            internal_pipeline_count = 0
            other_count = 0
            
            for f in files:
                if f.startswith('/outs/') or f.startswith('/outs'):
                    outs_count += 1
                elif '/SC_RNA_COUNTER_CS/' in f or '/SC_' in f or f.startswith('/SC_'):
                    internal_pipeline_count += 1
                else:
                    other_count += 1
            
            print(f"    Files in outs/: {outs_count}")
            print(f"    Internal pipeline files: {internal_pipeline_count}")
            if other_count:
                print(f"    Other files: {other_count}")
            
            # SOP Warning: Internal pipeline files should not be uploaded
            if internal_pipeline_count > 0:
                print(f"    ⚠ SOP WARNING: {internal_pipeline_count} internal pipeline files found")
                print(f"       Per SOP, only the Cell Ranger outs/ directory should be uploaded")
            
            # Check for key Cell Ranger outputs
            if pipeline == 'cellranger':
                if outs_count == 0:
                    print(f"    ⚠ SOP WARNING: No outs/ files found!")
                    print(f"       Per SOP, Cell Ranger outs/ directory should be uploaded")
                else:
                    missing_key_outputs = []
                    found_key_outputs = []
                    
                    # Check standard outputs
                    for key_output in CELLRANGER_KEY_OUTPUTS:
                        if key_output in outs_files:
                            found_key_outputs.append(key_output)
                        else:
                            missing_key_outputs.append(key_output)
                    
                    # Check multi outputs
                    for key_output in CELLRANGER_MULTI_KEY_OUTPUTS:
                        if key_output in outs_files:
                            found_key_outputs.append(key_output)
                    
                    if found_key_outputs:
                        print(f"    ✓ Key outputs found: {len(found_key_outputs)}")
                        for ko in found_key_outputs:
                            print(f"        {ko}")
                    
                    if missing_key_outputs:
                        print(f"    ⚠ Key outputs missing: {len(missing_key_outputs)}")
                        for mo in missing_key_outputs:
                            print(f"        {mo}")
            
            # Show internal file breakdown if present
            if internal_pipeline_count > 0:
                internal_types = defaultdict(int)
                for f in files:
                    if '/SC_RNA_COUNTER_CS/' in f or '/SC_' in f or f.startswith('/SC_'):
                        filename = f.split('/')[-1] if f else 'root'
                        if filename.startswith('_'):
                            internal_types[filename] += 1
                        else:
                            internal_types['other'] += 1
                
                print(f"    Internal file types (should not be uploaded per SOP):")
                for itype, count in sorted(internal_types.items(), key=lambda x: -x[1])[:5]:
                    print(f"        {itype}: {count}")
                if len(internal_types) > 5:
                    print(f"        ... and {len(internal_types) - 5} more types")
    
    # Check for expected GroupIDs not found at all
    found_groupids = set(t['groupid'] for t in transformations)
    missing_groupids = expected_groupids - found_groupids
    if missing_groupids:
        print(f"\n⚠ Expected GroupIDs not found in mapping: {', '.join(sorted(missing_groupids))}")
    
    # Check for GroupIDs missing outs/ files specifically
    groupids_with_outs = set()
    for sample_key in outs_files_per_sample.keys():
        groupid, pipeline, run_id = sample_key
        if outs_files_per_sample[sample_key]:  # has at least one outs file
            groupids_with_outs.add(groupid)
    
    groupids_missing_outs = expected_groupids - groupids_with_outs
    if groupids_missing_outs:
        print(f"\n❌ GroupIDs missing outs/ files: {', '.join(sorted(groupids_missing_outs))}")
        print(f"   These groups have no Cell Ranger output files in the mapping!")
    elif expected_groupids:
        print(f"\n✓ All {len(expected_groupids)} expected GroupIDs have outs/ files")
    
    # Report SOP warnings
    if sop_warnings:
        print(f"\n⚠ SOP warnings: {len(sop_warnings)}")
        for line_num, s3_path, warnings in sop_warnings:
            print(f"  Line {line_num}:")
            for w in warnings:
                print(f"    {w}")
    
    print(f"\n" + "-" * 40)
    print("MAPPING LOGIC DISCOVERED:")
    print("-" * 40)
    print(f"""
The transformation from local to S3 follows this pattern (per SOP):

S3:     s3://czi-{{provider}}/{{project}}/{{order}}/{{GroupID}}/processed/{{pipeline}}/Run_{{date}}/{{rest_of_path}}
LOCAL:  {{local_base_path}}/{{rest_of_path}}

Expected Group IDs: {', '.join(sorted(sample_ids))}
Valid Pipelines (per SOP): {', '.join(VALID_PIPELINES)}

SOP Requirements:
1. S3 bucket format: czi-{{provider}}
2. GroupID should match raw fastq files from raw/ subdirectory
3. Run_Date format: Run_YYYY-MM-DD (date Cell Ranger run started)
4. Cell Ranger outs/ directory uploaded without changing file/folder names
5. Example: {{order}}/{{GroupID}}/processed/cellranger/Run_2025-02-01/outs/multi/count/raw_feature_bc_matrix.h5
""")
    
    return outliers



def load_sif(sif_filepath):
    """
    Load a Sample Intake Form CSV and return two lookup dicts:
      sublibrary_to_library: normalized sublibrary -> library name
                             e.g. 'QSR5' -> 'R096E', 'QSR5SCALEPLEX' -> 'R096E'
      index_to_library:      index sequence -> library name
                             e.g. 'CTACATGCA' -> 'R096E'
    Normalization: strips hyphens from sublibrary names so QSR-5 == QSR5.
    """
    import csv

    sublibrary_to_library = {}
    index_to_library = {}

    encodings = ['utf-8', 'latin-1', 'cp1252']
    for enc in encodings:
        try:
            with open(sif_filepath, encoding=enc) as f:
                reader = csv.DictReader(f)
                cols = reader.fieldnames
                if not cols:
                    continue
                lib_col = cols[0]
                sub_col = cols[1]
                idx_col = cols[2]
                for row in reader:
                    lib = row[lib_col].strip()
                    sub = row[sub_col].strip()
                    idx = row[idx_col].strip()
                    if not lib:
                        continue
                    if sub:
                        # normalize: remove hyphens so QSR-5 == QSR5
                        sub_norm = sub.replace('-', '')
                        sublibrary_to_library[sub_norm] = lib
                    if idx:
                        index_to_library[idx] = lib
            break
        except UnicodeDecodeError:
            continue

    return sublibrary_to_library, index_to_library


def validate_sif_completeness_scale(mappings, sif_filepath):
    """
    Validate that all SIF entries (sublibrary combos) are present in the mapping.
    
    SOP format for Scale/ScaleBio filenames:
      {RunID}-{GroupID}_{Assay}_{UG_RT}
    
    Where:
      - RunID: sequencing run ID (e.g., 430351, 441969)
      - GroupID: sample identifier (e.g., R098A, RNA3-096H)
      - Assay: GEX, hash_oligo, or GEX_hash_oligo
      - UG_RT: Ultima barcode output format: QSR-{N}[-SCALEPLEX]_{well_position}
               Example: QSR-8_10A, QSR-8-SCALEPLEX_10A
    
    Example: 430351-RNA3-096H_hash_oligo_QSR-8-SCALEPLEX_10A.cram
    
    Validation checks that all SIF sublibraries appear in the mapping with correct GroupID.
    """
    import csv
    
    # Load SIF to get expected sublibrary -> GroupID mapping
    expected_sublibraries = {}
    
    encodings = ['utf-8', 'latin-1', 'cp1252']
    for enc in encodings:
        try:
            with open(sif_filepath, encoding=enc) as f:
                reader = csv.DictReader(f)
                cols = reader.fieldnames
                if not cols:
                    continue
                
                lib_col = cols[0]        # Library / GroupID column
                sub_col = cols[1]        # Sublibrary column
                idx_col = cols[2]        # Index column
                
                for row in reader:
                    lib = row[lib_col].strip() if lib_col in row else ""
                    sub = row[sub_col].strip() if sub_col in row else ""
                    idx = row[idx_col].strip() if idx_col in row else ""
                    
                    if not lib or not sub or not idx:
                        continue
                    
                    # Normalize sublibrary: remove hyphens
                    sub_norm = sub.replace('-', '').upper()
                    expected_sublibraries[sub_norm] = lib
            break
        except (UnicodeDecodeError, KeyError):
            continue
    
    # Extract found entries from mapping
    found_sublibraries = defaultdict(set)  # sublibrary_norm -> set of (assay, groupid)
    
    # Also build reverse lookup: index_sequence -> sublibrary_norm (from SIF)
    import csv as csv_mod
    index_to_sublibrary = {}
    for enc in ['utf-8', 'latin-1', 'cp1252']:
        try:
            with open(sif_filepath, encoding=enc) as f:
                reader = csv_mod.DictReader(f)
                cols = reader.fieldnames
                if not cols:
                    continue
                lib_col, sub_col, idx_col = cols[0], cols[1], cols[2]
                for row in reader:
                    sub = row[sub_col].strip() if sub_col in row else ""
                    idx = row[idx_col].strip() if idx_col in row else ""
                    if sub and idx:
                        sub_norm = sub.replace('-', '').upper()
                        index_to_sublibrary[idx] = sub_norm
            break
        except (UnicodeDecodeError, KeyError):
            continue
    
    for s3_path, local_path, line_num in mappings:
        filename = s3_path.split('/')[-1]
        filename_base = filename.rsplit('.', 1)[0]  # Remove extension
        
        # Try SOP form first: {RunID}-{GroupID}_{Assay}_{QSR-N[-SCALEPLEX]}-{well}...
        match = re.match(
            r'(\d+)-(.+?)_(GEX|Hash_oligo|GEX_hash_oligo)_(QSR-\d+(?:-SCALEPLEX)?)-([A-Z0-9]+)',
            filename_base
        )
        
        if match:
            groupid = match.group(2)
            assay = match.group(3)
            sublibrary = match.group(4)
            sub_norm = sublibrary.replace('-', '').upper()
            found_sublibraries[sub_norm].add((assay, groupid))
            continue
        
        # Try index-direct form: {RunID}-{GroupID}_{Assay}_{IndexSequence}...
        match_idx = re.match(
            r'(\d+)-(.+?)_(GEX|Hash_oligo|GEX_hash_oligo)_([ACGT]+)',
            filename_base
        )
        
        if match_idx:
            groupid = match_idx.group(2)
            assay = match_idx.group(3)
            index_seq = match_idx.group(4)
            # Map index sequence back to sublibrary via SIF
            sub_norm = index_to_sublibrary.get(index_seq)
            if sub_norm:
                found_sublibraries[sub_norm].add((assay, groupid))
            # else: index not in SIF — will be caught by main validation
    
    # Validate: check each expected sublibrary is in the mapping
    missing_sublibraries = []
    found_count = 0
    
    for sub_norm, expected_groupid in sorted(expected_sublibraries.items()):
        if sub_norm in found_sublibraries:
            found_count += 1
        else:
            missing_sublibraries.append((sub_norm, expected_groupid))
    
    # Check for GroupID mismatches
    groupid_mismatches = []
    for sub_norm, entries in found_sublibraries.items():
        expected_groupid = expected_sublibraries.get(sub_norm)
        if expected_groupid:
            found_groupids = set(e[1] for e in entries)  # e[1] is groupid
            if expected_groupid not in found_groupids:
                groupid_mismatches.append((sub_norm, expected_groupid, found_groupids))
    
    expected_count = len(expected_sublibraries)
    coverage = found_count / expected_count if expected_count else 1.0
    
    return {
        'expected_sublibraries': expected_sublibraries,
        'found_sublibraries': dict(found_sublibraries),
        'missing_sublibraries': missing_sublibraries,
        'found_count': found_count,
        'expected_count': expected_count,
        'coverage': coverage,
        'groupid_mismatches': groupid_mismatches,
    }


def report_sif_completeness_scale(result):
    """Print a formatted report of SIF completeness validation for Scale mode."""
    expected = result['expected_count']
    found = result['found_count']
    coverage = result['coverage']
    
    print(f"\nSIF Completeness Check:")
    print(f"  Expected sublibraries (from SIF): {expected}")
    print(f"  Found in mapping: {found}")
    print(f"  Coverage: {coverage*100:.1f}%")
    
    if coverage == 1.0:
        print(f"\n  ✓ All SIF entries are present in mapping")
    else:
        print(f"\n  ❌ INCOMPLETE: {len(result['missing_sublibraries'])} SIF entries missing from mapping")
        for sub_norm, expected_groupid in result['missing_sublibraries']:
            print(f"     - {sub_norm} (should map to GroupID {expected_groupid})")
    
    if result['groupid_mismatches']:
        print(f"\n  ⚠️  GroupID MISMATCHES ({len(result['groupid_mismatches'])} found):")
        for sub_norm, expected_gid, found_gids in result['groupid_mismatches']:
            print(f"     - {sub_norm}: expected {expected_gid}, found {found_gids}")
    
    if coverage < 1.0:
        print(f"\n  → Check if sequencing is still in progress or if there's a sample manifest error")
        return False
    else:
        return True



def validate_sop_compliance(mappings):
    """
    Validate that S3 paths comply with the official SOP specification.
    
    SOP S3 Path Format:
      s3://czi-{provider}/{lastname}-{projectname}/{order_number}/{ExperimentID}/raw/{RunID}/
      {RunID}-{GroupID}_{Assay}_{UG_RT}{suffix}
    
    Where:
      - provider: novogene, ultima, etc.
      - lastname-projectname: e.g., trapnell-seahub-bcp
      - order_number: e.g., NVUS2024101701-04
      - ExperimentID: e.g., RNA3_096, CHEM13-R096
      - RunID: 6-digit run ID (430351, 441969)
      - GroupID: sample identifier (RNA3-096H, R098A)
      - Assay: GEX, hash_oligo, or GEX_hash_oligo
      - UG_RT: UG 100 output format: QSR-{N}[-SCALEPLEX]_{well_position}
      - suffix: file extension (.cram, .csv, .json, etc.)
    
    SOP Rules:
    1. Assay must be one of: GEX, hash_oligo, GEX_hash_oligo
    2. UG_RT must contain QSR-{N} pattern
    3. If Assay is hash_oligo, UG_RT MUST contain -SCALEPLEX
    4. If Assay is GEX or GEX_hash_oligo, UG_RT MUST NOT contain SCALEPLEX
    5. UG_RT format: QSR-{N}[-SCALEPLEX]_{well_position}
    
    Returns dict with violations found.
    """
    
    # Parse S3 path structure
    s3_structure_pattern = re.compile(
        r's3://czi-([^/]+)/([^/]+)/([^/]+)/([^/]+)/raw/(\d+)/(.+)$'
    )
    
    # Parse S3 filename: {RunID}-{GroupID}_{Assay}_{UG_RT}
    s3_filename_pattern = re.compile(
        r'(\d+)-(.+?)_(GEX|hash_oligo|GEX_hash_oligo)_(QSR-\d+(?:-SCALEPLEX)?_[A-Z0-9]+)'
    )
    
    violations = []
    total_s3_paths = 0
    
    for s3_path, local_path, line_num in mappings:
        # Parse S3 structure
        struct_match = s3_structure_pattern.match(s3_path)
        if not struct_match:
            # Not a valid S3 path structure
            violations.append((line_num, s3_path, "S3 path doesn't match SOP structure"))
            continue
        
        provider = struct_match.group(1)           # e.g., novogene
        project = struct_match.group(2)            # e.g., trapnell-seahub-bcp
        order = struct_match.group(3)              # e.g., NVUS2024101701-04
        experiment_id = struct_match.group(4)      # e.g., RNA3_096
        runid = struct_match.group(5)              # e.g., 430351
        filename = struct_match.group(6)           # e.g., 430351-RNA3-096H_GEX_QSR-8_10A.cram
        
        # Parse filename
        filename_base = filename.rsplit('.', 1)[0]  # Remove suffix
        file_match = s3_filename_pattern.match(filename_base)
        
        if not file_match:
            violations.append((line_num, s3_path, f"Filename doesn't match SOP pattern: {filename}"))
            continue
        
        total_s3_paths += 1
        
        file_runid = file_match.group(1)
        groupid = file_match.group(2)
        assay = file_match.group(3)
        ug_rt = file_match.group(4)
        
        issues = []
        
        # Rule 1: Assay must be one of GEX, hash_oligo, GEX_hash_oligo
        if assay not in ['GEX', 'hash_oligo', 'GEX_hash_oligo']:
            issues.append(f"Invalid Assay '{assay}' (must be GEX, hash_oligo, or GEX_hash_oligo)")
        
        # Rule 2: UG_RT must contain QSR-{N} pattern
        if not re.match(r'QSR-\d+(?:-SCALEPLEX)?_[A-Z0-9]+', ug_rt):
            issues.append(f"UG_RT format invalid: '{ug_rt}' (must be QSR-{{N}}[-SCALEPLEX]_{{well}})")
        
        # Rule 3 & 4: SCALEPLEX presence must match Assay
        has_scaleplex = 'SCALEPLEX' in ug_rt
        
        if assay == 'hash_oligo' and not has_scaleplex:
            issues.append(f"Assay is hash_oligo but UG_RT doesn't contain SCALEPLEX: '{ug_rt}'")
        elif assay == 'GEX' and has_scaleplex:
            issues.append(f"Assay is GEX but UG_RT contains SCALEPLEX: '{ug_rt}'")
        elif assay == 'GEX_hash_oligo' and has_scaleplex:
            issues.append(f"Assay is GEX_hash_oligo but UG_RT contains SCALEPLEX: '{ug_rt}'")
        
        # Rule 5: Check RunID in filename matches path RunID
        if file_runid != runid:
            issues.append(f"RunID mismatch: path has {runid}, filename has {file_runid}")
        
        if issues:
            violations.append((line_num, s3_path, issues))
    
    return {
        'violations': violations,
        'total_s3_paths': total_s3_paths,
        'violation_count': len(violations),
    }


def report_sop_compliance(result):
    """Print a formatted report of SOP compliance validation."""
    total = result['total_s3_paths']
    violations = result['violation_count']
    
    print(f"\nSOP Compliance Check (S3 paths):")
    print(f"  Total S3 paths checked: {total}")
    print(f"  SOP violations found: {violations}")
    
    if violations == 0:
        print(f"\n  ✓ All S3 paths comply with SOP")
    else:
        print(f"\n  ❌ FOUND {violations} SOP violations:")
        for line_num, s3_path, issues in result['violations'][:10]:
            print(f"\n  Line {line_num}:")
            print(f"    S3: {s3_path}")
            if isinstance(issues, str):
                print(f"    Issue: {issues}")
            else:
                for issue in issues:
                    print(f"    Issue: {issue}")
        
        if len(result['violations']) > 10:
            print(f"\n  ... and {len(result['violations']) - 10} more violations")
    
    return violations == 0

    """Print a formatted report of SIF completeness validation for Scale mode."""
    expected = result['expected_count']
    found = result['found_count']
    coverage = result['coverage']
    
    print(f"\nSIF Completeness Check:")
    print(f"  Expected sublibraries (from SIF): {expected}")
    print(f"  Found in mapping: {found}")
    print(f"  Coverage: {coverage*100:.1f}%")
    
    if coverage == 1.0:
        print(f"\n  ✓ All SIF entries are present in mapping")
    else:
        print(f"\n  ❌ INCOMPLETE: {len(result['missing_sublibraries'])} SIF entries missing from mapping")
        for sub_norm, expected_groupid in result['missing_sublibraries']:
            print(f"     - {sub_norm} (should map to GroupID {expected_groupid})")
    
    if result['groupid_mismatches']:
        print(f"\n  ⚠️  GroupID MISMATCHES ({len(result['groupid_mismatches'])} found):")
        for sub_norm, expected_gid, found_gids in result['groupid_mismatches']:
            print(f"     - {sub_norm}: expected {expected_gid}, found {found_gids}")
    
    if coverage < 1.0:
        print(f"\n  → Check if sequencing is still in progress or if there's a sample manifest error")
        return False
    else:
        return True


def validate_s3_local_consistency(mappings, sif_filepath=None):
    """
    Validate that S3 paths and local paths describe the same samples.
    
    Checks:
    - For SOP form: Sublibrary name consistency (QSR-8 in S3 should match QSR8 in local)
    - For index-direct form: Index sequence consistency + SIF-mediated GroupID check
    - SCALEPLEX presence consistency (Hash_oligo in S3 should have SCALEPLEX in local)
    - GroupID consistency (S3 GroupID should match local sublibrary mapping)
    
    Returns dict with:
    - inconsistent_pairs: list of (line_num, s3_path, local_path, issues)
    - total_checked: number of (S3, local) pairs checked
    - inconsistent_count: number of inconsistent pairs found
    """
    
    # Load SIF index->sublibrary mapping if available (needed for index-direct form)
    index_to_sublibrary_norm = {}
    if sif_filepath:
        import csv as csv_mod
        for enc in ['utf-8', 'latin-1', 'cp1252']:
            try:
                with open(sif_filepath, encoding=enc) as f:
                    reader = csv_mod.DictReader(f)
                    cols = reader.fieldnames
                    if not cols:
                        continue
                    sub_col, idx_col = cols[1], cols[2]
                    for row in reader:
                        sub = row[sub_col].strip() if sub_col in row else ""
                        idx = row[idx_col].strip() if idx_col in row else ""
                        if sub and idx:
                            index_to_sublibrary_norm[idx] = sub.replace('-', '').upper()
                break
            except (UnicodeDecodeError, KeyError):
                continue

    # S3 regex: SOP form
    s3_pattern_sop = re.compile(
        r'(\d+)-(.+?)_(GEX|Hash_oligo|GEX_hash_oligo)_(QSR-\d+(?:-SCALEPLEX)?)'
    )
    # S3 regex: index-direct form
    s3_pattern_idx = re.compile(
        r'(\d+)-(.+?)_(GEX|Hash_oligo|GEX_hash_oligo)_([ACGT]+)'
    )
    
    # Local pattern (index-direct): extract QSR-N and optional SCALEPLEX and index
    local_pattern_idx = re.compile(
        r'(\d+)-QSR-(\d+)(?:-(SCALEPLEX))?-([ACGT]+)'
    )
    # Local pattern (SOP): extract sublibNorm and QSR-N[-SCALEPLEX]
    local_pattern_sop = re.compile(
        r'(\d+)-(QSR\d+(?:SCALEPLEX)?)_QSR-(\d+)(?:-(SCALEPLEX))?/'
    )
    
    inconsistent_pairs = []
    total_checked = 0
    
    for s3_path, local_path, line_num in mappings:
        # Extract from S3
        s3_filename = s3_path.split('/')[-1].rsplit('.', 1)[0]
        
        s3_match = s3_pattern_sop.match(s3_filename)
        s3_form = 'sop'
        if not s3_match:
            s3_match = s3_pattern_idx.match(s3_filename)
            s3_form = 'index_direct'
        
        # Extract from local — try index-direct first, then SOP
        local_match = local_pattern_idx.search(local_path)
        local_form = 'index_direct'
        if not local_match:
            local_match = local_pattern_sop.search(local_path)
            local_form = 'sop'
        
        if not s3_match or not local_match:
            continue
        
        total_checked += 1
        issues = []
        
        # Parse local
        if local_form == 'index_direct':
            local_qsr_n = local_match.group(2)
            local_has_scaleplex = local_match.group(3) is not None
            local_index = local_match.group(4)
            local_sublibrary_norm = f"QSR{local_qsr_n}"
            if local_has_scaleplex:
                local_sublibrary_norm += "SCALEPLEX"
        else:
            local_sublibrary_norm = local_match.group(2).upper()
            local_qsr_n = local_match.group(3)
            local_has_scaleplex = local_match.group(4) is not None
            local_index = None
        
        # Parse S3
        s3_groupid = s3_match.group(2)
        s3_assay = s3_match.group(3)
        
        if s3_form == 'sop':
            s3_sublibrary = s3_match.group(4)
            s3_sublibrary_norm = s3_sublibrary.replace('-', '').upper()
            s3_has_scaleplex = 'SCALEPLEX' in s3_sublibrary
            
            # 1. Sublibrary name must match
            if s3_sublibrary_norm != local_sublibrary_norm:
                issues.append(
                    f"Sublibrary mismatch: S3 has {s3_sublibrary_norm}, local has {local_sublibrary_norm}"
                )
            
            # 2. SCALEPLEX presence must match
            if s3_has_scaleplex != local_has_scaleplex:
                s3_desc = "has SCALEPLEX" if s3_has_scaleplex else "no SCALEPLEX"
                local_desc = "has SCALEPLEX" if local_has_scaleplex else "no SCALEPLEX"
                issues.append(
                    f"SCALEPLEX mismatch: S3 {s3_desc}, local {local_desc}"
                )
        else:
            # Index-direct form
            s3_index = s3_match.group(4)
            
            # 1. Index sequences must match between S3 and local (if local has them)
            if local_index and s3_index != local_index:
                issues.append(
                    f"Index mismatch: S3 has {s3_index}, local has {local_index}"
                )
            
            # 2. If SIF available, verify S3 index -> sublibrary matches local sublibrary
            if index_to_sublibrary_norm:
                s3_sublibrary_via_sif = index_to_sublibrary_norm.get(s3_index)
                if s3_sublibrary_via_sif and s3_sublibrary_via_sif != local_sublibrary_norm:
                    issues.append(
                        f"SIF sublibrary mismatch: S3 index {s3_index} -> {s3_sublibrary_via_sif} (SIF), "
                        f"but local has {local_sublibrary_norm}"
                    )
        
        # 3. Assay should match SCALEPLEX presence
        if local_has_scaleplex and s3_assay not in ("Hash_oligo", "GEX_hash_oligo"):
            issues.append(
                f"Assay mismatch: local has SCALEPLEX but S3 assay is {s3_assay} (expected Hash_oligo)"
            )
        elif not local_has_scaleplex and s3_assay == "Hash_oligo":
            issues.append(
                f"Assay mismatch: local has no SCALEPLEX but S3 assay is Hash_oligo"
            )
        
        if issues:
            inconsistent_pairs.append((line_num, s3_path, local_path, issues))
    
    return {
        'inconsistent_pairs': inconsistent_pairs,
        'total_checked': total_checked,
        'inconsistent_count': len(inconsistent_pairs),
    }


def report_s3_local_consistency(result):
    """Print a formatted report of S3/local consistency validation."""
    total = result['total_checked']
    inconsistent = result['inconsistent_count']
    
    print(f"\nS3/Local Path Consistency Check:")
    print(f"  Pairs checked: {total}")
    print(f"  Inconsistent pairs: {inconsistent}")
    
    if inconsistent == 0:
        print(f"\n  ✓ All S3 and local paths are consistent")
    else:
        print(f"\n  ❌ FOUND {inconsistent} inconsistent pairs:")
        for line_num, s3_path, local_path, issues in result['inconsistent_pairs'][:10]:
            print(f"\n  Line {line_num}:")
            print(f"    S3: {s3_path}")
            print(f"    Local: {local_path}")
            for issue in issues:
                print(f"    Issue: {issue}")
        
        if len(result['inconsistent_pairs']) > 10:
            print(f"\n  ... and {len(result['inconsistent_pairs']) - 10} more")
    
    return inconsistent == 0


def analyze_mapping_logic_scale(mappings, sif_filepath):
    """
    Analyze Scale/scRNA-seq mappings that use the trapnell-seahub-bcp bucket.

    Scale data distinguishes between ExperimentID and GroupID:
      - ExperimentID: the path-level directory (e.g. CHEM13-R096), a chip/experiment
        container that holds multiple sub-libraries
      - GroupID: the filename-level library identifier (e.g. R096G), the actual
        per-sublibrary identifier that maps to SIF entries

    S3 format:
      s3://czi-novogene/trapnell-seahub-bcp/{order}/{ExperimentID}/raw/{wafer}/
          {wafer}-{GroupID}_{assay}_{index}{suffix}
      where assay is GEX or Hash_oligo

    Local format (GEX):
      {base}/{wafer}-{date}/{wafer}-QSR-{N}-{index}/{wafer}-QSR-{N}-{index}{suffix}

    Local format (Hash_oligo / SCALEPLEX):
      {base}/{wafer}-{date}/{wafer}-QSR-{N}-SCALEPLEX-{index}/{wafer}-QSR-{N}-SCALEPLEX-{index}{suffix}

    Validation uses the SIF to map sublibrary (QSR5/QSR-5) -> GroupID (R096E).
    The ExperimentID is validated independently as a consistent container.
    Suffix normalization: S3 uses hyphens (-trimmer-stats.csv),
                          local uses underscores (_trimmer-stats.csv).
    """
    print("\n" + "=" * 80)
    print("2. MAPPING LOGIC ANALYSIS (Scale mode)")
    print("=" * 80)

    # Load SIF
    sublibrary_to_library, index_to_library = load_sif(sif_filepath)
    print(f"\nLoaded SIF: {len(sublibrary_to_library)} sublibrary entries, "
          f"{len(index_to_library)} index entries")
    print("  Sublibrary -> Library:")
    for sub, lib in sorted(sublibrary_to_library.items()):
        print(f"    {sub} -> {lib}")
    
    # STEP 1: SIF COMPLETENESS VALIDATION
    print("\n" + "-" * 80)
    print("STEP 1: SIF COMPLETENESS CHECK")
    print("-" * 80)
    sif_completeness = validate_sif_completeness_scale(mappings, sif_filepath)
    report_sif_completeness_scale(sif_completeness)
    
    # STEP 2: S3/LOCAL PATH CONSISTENCY CHECK
    print("\n" + "-" * 80)
    print("STEP 2: S3/LOCAL PATH CONSISTENCY CHECK")
    print("-" * 80)
    s3_local_consistency = validate_s3_local_consistency(mappings, sif_filepath)
    report_s3_local_consistency(s3_local_consistency)

    # S3 pattern per SOP:
    # s3://czi-{provider}/{lastname-projectname}/{order}/{ExperimentID}/raw/{RunID}/
    # {RunID}-{GroupID}_{Assay}_{UG_RT}{suffix}
    #
    # Two known naming conventions for the component after Assay:
    #   1. SOP/UG_RT form: QSR-{N}[-SCALEPLEX]_{wellPosition}
    #      Example: 430351-RNA3-096H_hash_oligo_QSR-8-SCALEPLEX_10A.cram
    #   2. Index-direct form: {IndexSequence}  (Ultima index from SIF)
    #      Example: 441969-R096G_GEX_CTATGCACA.cram
    #
    # We try the SOP form first, then fall back to index-direct.
    s3_sample_pattern_sop = re.compile(
        r's3://(?P<bucket>[^/]+)/(?P<project>[^/]+)/'
        r'(?P<order>[A-Z0-9]+(?:-\d+)*)/(?P<experiment_id>[^/]+)/raw/'
        r'(?P<runid>\d+)/'
        r'(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>GEX|Hash_oligo|GEX_hash_oligo|CRI|ATAC)_(?P<ug_rt>QSR-\d+(?:-SCALEPLEX)?)'
        r'(?P<suffix>.*)'
    )
    s3_sample_pattern_idx = re.compile(
        r's3://(?P<bucket>[^/]+)/(?P<project>[^/]+)/'
        r'(?P<order>[A-Z0-9]+(?:-\d+)*)/(?P<experiment_id>[^/]+)/raw/'
        r'(?P<runid>\d+)/'
        r'(?P<runid2>\d+)-(?P<group_id>[A-Za-z0-9_-]+)_(?P<assay>GEX|Hash_oligo|GEX_hash_oligo|CRI|ATAC)_(?P<index_seq>[ACGT]+)'
        r'(?P<suffix>.*)'
    )

    # S3 order-level batch files (no library/raw/ in path)
    s3_order_pattern = re.compile(
        r's3://(?P<bucket>[^/]+)/(?P<project>[^/]+)/'
        r'(?P<order>NVUS\d+(?:-\d+)+)/(?P<filename>[^/]+)$'
    )

    # S3 run-level batch files (under /raw/{RunID}/ but not sample-level)
    # e.g., 441969_merged_trimmer-stats.csv, 441969_UploadCompleted.json, etc.
    s3_run_batch_pattern = re.compile(
        r's3://(?P<bucket>[^/]+)/(?P<project>[^/]+)/'
        r'(?P<order>[A-Z0-9]+(?:-\d+)*)/(?P<experiment_id>[^/]+)/raw/'
        r'(?P<runid>\d+)/(?P<filename>\d+_.+)$'
    )

    # Local sample pattern (index-direct form):
    # {base}/{wafer}-{date}/{wafer}-QSR-{N}[-SCALEPLEX]-{index}/{wafer}-QSR-{N}[-SCALEPLEX]-{index}{suffix}
    # Example: 441969-QSR-7-CTATGCACA/441969-QSR-7-CTATGCACA.json
    local_sample_pattern_idx = re.compile(
        r'(?P<base_path>.+)/'
        r'(?P<wafer>\d+)-(?P<date>\d+_\d+)/'
        r'(?P<wafer2>\d+)-QSR-(?P<qsr_n>\d+)(?P<scaleplex>-SCALEPLEX)?-(?P<index>[ACGT]+)/'
        r'(?P<wafer3>\d+)-QSR-(?P<qsr_n2>\d+)(?P<scaleplex2>-SCALEPLEX)?-(?P<index2>[ACGT]+)'
        r'(?P<suffix>.*)'
    )

    # Local sample pattern (SOP UG_RT form):
    # {base}/{wafer}-{date}/{wafer}-{sublibNorm}_{QSR-N[-SCALEPLEX]}/{wafer}-{sublibNorm}_{QSR-N[-SCALEPLEX]}_{well_or_suffix}
    # Example: 441969-QSR8_QSR-8/441969-QSR8_QSR-8_12B.csv
    #          441969-QSR8SCALEPLEX_QSR-8-SCALEPLEX/441969-QSR8SCALEPLEX_QSR-8-SCALEPLEX_12A.csv
    #          441969-QSR8_QSR-8/441969-QSR8_QSR-8_trimmer-failure_codes.csv
    local_sample_pattern_sop = re.compile(
        r'(?P<base_path>.+)/'
        r'(?P<wafer>\d+)-(?P<date>\d+_\d+)/'
        r'(?P<wafer2>\d+)-(?P<sublib_norm>QSR\d+(?:SCALEPLEX)?)_QSR-(?P<qsr_n>\d+)(?P<scaleplex>-SCALEPLEX)?/'
        r'(?P<wafer3>\d+)-(?P<sublib_norm2>QSR\d+(?:SCALEPLEX)?)_QSR-(?P<qsr_n2>\d+)(?P<scaleplex2>-SCALEPLEX)?'
        r'(?P<suffix>.*)'
    )

    # Local batch pattern (no sample subdir): {base}/{wafer}-{date}/{filename}
    local_batch_pattern = re.compile(
        r'(?P<base_path>.+)/'
        r'(?P<wafer>\d+)-(?P<date>\d+_\d+)/'
        r'(?P<filename>[^/]+)$'
    )

    def normalize_s3_suffix(s3_suffix):
        """Normalize S3 suffix to local style for comparison."""
        return (s3_suffix
                .replace('-trimmer-failure-codes.csv', '_trimmer-failure_codes.csv')
                .replace('-trimmer-stats.csv', '_trimmer-stats.csv'))

    transformations = []
    batch_files = []
    outliers = []

    # Track files per (wafer, library, assay) for completeness
    files_per_lib = defaultdict(set)

    # Track which S3 naming convention is in use
    s3_naming_convention = None  # will be set to 'sop' or 'index_direct'

    for s3_path, local_path, line_num in mappings:
        # Try SOP form first, then index-direct form for S3
        s3_match = s3_sample_pattern_sop.match(s3_path)
        s3_form = 'sop'
        if not s3_match:
            s3_match = s3_sample_pattern_idx.match(s3_path)
            s3_form = 'index_direct'

        # Try both local patterns: index-direct first (original), then SOP
        local_match = local_sample_pattern_idx.match(local_path)
        local_form = 'index_direct'
        if not local_match:
            local_match = local_sample_pattern_sop.match(local_path)
            local_form = 'sop'

        if s3_match and local_match:
            sg = s3_match.groupdict()
            lg = local_match.groupdict()
            issues = []

            # Track naming convention (detect mixed usage as a warning)
            if s3_naming_convention is None:
                s3_naming_convention = s3_form
            elif s3_naming_convention != s3_form:
                s3_naming_convention = 'mixed'

            assay = sg['assay']

            # --- Extract local sublibrary info (depends on local_form) ---
            if local_form == 'index_direct':
                # Local: {wafer}-QSR-{N}[-SCALEPLEX]-{ACGT_index}/...
                scaleplex = lg['scaleplex'] is not None
                sub_normalized = f"QSR{lg['qsr_n']}{'SCALEPLEX' if scaleplex else ''}"
                local_index = lg.get('index')
            else:
                # Local SOP: {wafer}-{sublibNorm}_QSR-{N}[-SCALEPLEX]/...
                # sublibNorm is already e.g. "QSR8" or "QSR8SCALEPLEX"
                sub_normalized = lg['sublib_norm'].upper()
                scaleplex = lg['scaleplex'] is not None  # -SCALEPLEX in QSR-N part
                local_index = None  # SOP form uses well positions, not index seqs

            # --- S3 validation (depends on s3_form) ---
            if s3_form == 'sop':
                ug_rt = sg['ug_rt']
                has_scaleplex_in_ug = 'SCALEPLEX' in ug_rt

                # Rule: hash_oligo MUST have SCALEPLEX in UG_RT
                if assay == 'Hash_oligo' and not has_scaleplex_in_ug:
                    issues.append(
                        f"SOP violation: Assay is 'Hash_oligo' but UG_RT '{ug_rt}' doesn't contain SCALEPLEX"
                    )
                # Rule: GEX MUST NOT have SCALEPLEX in UG_RT
                if assay == 'GEX' and has_scaleplex_in_ug:
                    issues.append(
                        f"SOP violation: Assay is 'GEX' but UG_RT '{ug_rt}' contains SCALEPLEX"
                    )

                ug_rt_or_index = ug_rt

            else:
                # Index-direct S3 form
                index_seq = sg['index_seq']

                # SIF lookup via index sequence -> expected GroupID
                expected_library_by_idx = index_to_library.get(index_seq)
                if expected_library_by_idx is None:
                    issues.append(
                        f"index sequence '{index_seq}' not found in SIF "
                        f"(known: {', '.join(sorted(index_to_library))})"
                    )
                elif expected_library_by_idx != sg['group_id']:
                    issues.append(
                        f"GroupID mismatch: index '{index_seq}' maps to "
                        f"'{expected_library_by_idx}' (SIF) but S3 filename has '{sg['group_id']}'"
                    )

                # Cross-check: if local has index sequence, it must match S3 index
                if local_index and local_index != index_seq:
                    issues.append(
                        f"index mismatch: local path has '{local_index}', "
                        f"S3 filename has '{index_seq}'"
                    )

                ug_rt_or_index = index_seq

            # --- Cross-validation: local sublibrary -> SIF -> GroupID ---
            expected_library = sublibrary_to_library.get(sub_normalized)
            if expected_library is None:
                issues.append(
                    f"sublibrary '{sub_normalized}' not found in SIF "
                    f"(known: {', '.join(sorted(sublibrary_to_library))})"
                )
            elif expected_library != sg['group_id']:
                issues.append(
                    f"GroupID mismatch: local sublibrary '{sub_normalized}' maps to "
                    f"'{expected_library}' (SIF) but S3 filename has '{sg['group_id']}'"
                )

            # --- SCALEPLEX / assay consistency ---
            if scaleplex and assay not in ('Hash_oligo', 'GEX_hash_oligo'):
                issues.append(
                    f"local has SCALEPLEX but S3 assay is '{assay}' "
                    f"(expected Hash_oligo or GEX_hash_oligo)"
                )
            if not scaleplex and assay == 'Hash_oligo':
                issues.append(
                    f"local has no SCALEPLEX but S3 assay is 'Hash_oligo'"
                )

            # RunID consistency: path RunID should match filename RunID
            if lg['wafer'] != sg['runid']:
                issues.append(f"RunID mismatch: local={lg['wafer']}, s3={sg['runid']}")

            # Note: ExperimentID (path) and GroupID (filename) are intentionally
            # different in Scale data. ExperimentID is a container (e.g. CHEM13-R096)
            # holding multiple GroupIDs (e.g. R096A, R096G, R098A).
            # Consistency is validated via SIF lookups above, not path comparison.

            # Suffix normalization check — only meaningful when local form matches
            # the S3 suffix structure. For SOP local, the suffix includes the well
            # position which is part of the S3 UG_RT, so direct comparison won't work.
            if local_form == 'index_direct':
                s3_suffix_normalized = normalize_s3_suffix(sg['suffix'])
                if lg['suffix'] != s3_suffix_normalized:
                    issues.append(
                        f"suffix mismatch: local={lg['suffix']!r}, "
                        f"s3 (normalized)={s3_suffix_normalized!r}"
                    )

            if issues:
                outliers.append((line_num, 'field_mismatch', local_path, s3_path, issues))
            else:
                files_per_lib[(sg['runid'], sg['group_id'], sg['assay'])].add(sg['suffix'])
                transformations.append({
                    'line': line_num,
                    'runid': sg['runid'],
                    'experiment_id': sg['experiment_id'],
                    'group_id': sg['group_id'],
                    'assay': sg['assay'],
                    'ug_rt': ug_rt_or_index,
                    's3_form': s3_form,
                    'local_form': local_form,
                    'sublibrary': sub_normalized,
                    'order': sg['order'],
                    'suffix': sg['suffix'],
                    'local_base': lg['base_path'],
                    'date': lg['date'],
                })
            continue

        # Try order-level batch file
        s3_order_match = s3_order_pattern.match(s3_path)
        local_batch_match = local_batch_pattern.match(local_path)
        # Also try run-level batch files (under /raw/{RunID}/ but not sample-level)
        s3_run_batch_match = s3_run_batch_pattern.match(s3_path) if not s3_order_match else None

        if s3_order_match and local_batch_match:
            batch_files.append({
                'line': line_num,
                'order': s3_order_match.group('order'),
                'filename': s3_order_match.group('filename'),
                'local_filename': local_batch_match.group('filename'),
            })
            continue
        if s3_run_batch_match and local_batch_match:
            batch_files.append({
                'line': line_num,
                'order': s3_run_batch_match.group('order'),
                'filename': s3_run_batch_match.group('filename'),
                'local_filename': local_batch_match.group('filename'),
            })
            continue

        # Neither matched
        parse_issues = []
        if not s3_match:
            parse_issues.append("S3 path doesn't match Scale sample pattern (tried SOP and index-direct forms)")
        if not local_match:
            parse_issues.append("local path doesn't match Scale sample pattern")
        if s3_match and not local_match:
            parse_issues.append("S3 looks like Scale sample but local doesn't match")
        if local_match and not s3_match:
            parse_issues.append("local looks like Scale sample but S3 doesn't match")
        outliers.append((line_num, 'parse_fail', local_path, s3_path,
                         parse_issues if parse_issues else None))

    # Report
    print(f"\nSuccessfully parsed:")
    print(f"  Sample-level files: {len(transformations)}")
    print(f"  Batch-level files:  {len(batch_files)}")
    print(f"  Total: {len(transformations) + len(batch_files)} / {len(mappings)} mappings")

    if s3_naming_convention:
        conv_desc = {
            'sop': 'SOP (QSR-N[-SCALEPLEX] in S3 filenames)',
            'index_direct': 'Index-direct (Ultima index sequence in S3 filenames)',
            'mixed': 'MIXED (both SOP and index-direct forms found!)',
        }
        print(f"\n  S3 naming convention: {conv_desc.get(s3_naming_convention, s3_naming_convention)}")
        if s3_naming_convention == 'index_direct':
            print(f"  ℹ️  NOTE: S3 filenames use index sequences instead of SOP UG_RT format.")
            print(f"     SOP expects: {{RunID}}-{{GroupID}}_{{Assay}}_QSR-{{N}}[-SCALEPLEX]_{{well}}")
            print(f"     Actual:      {{RunID}}-{{GroupID}}_{{Assay}}_{{IndexSequence}}")
            print(f"     SIF cross-references validate consistency.")
        elif s3_naming_convention == 'mixed':
            print(f"  ⚠️  WARNING: Mixed naming conventions detected — review for consistency.")

    if transformations:
        orders   = sorted(set(t['order']         for t in transformations))
        runids   = sorted(set(t['runid']         for t in transformations))
        exp_ids  = sorted(set(t['experiment_id'] for t in transformations))
        grp_ids  = sorted(set(t['group_id']      for t in transformations))
        sublibs  = sorted(set(t['sublibrary']    for t in transformations))
        assays   = sorted(set(t['assay']         for t in transformations))
        bases    = sorted(set(t['local_base']    for t in transformations))

        print(f"\nFields found:")
        print(f"  Orders:         {orders}")
        print(f"  RunIDs:         {runids}")
        print(f"  ExperimentIDs:  {exp_ids}")
        print(f"  GroupIDs:       {grp_ids}")
        print(f"  Sublibraries:   {sublibs}")
        print(f"  Assays:         {assays}")
        print(f"  Local bases:    {bases}")

        # Sublibrary -> GroupID mapping observed
        print(f"\n  Sublibrary -> GroupID (via SIF):")
        seen = {}
        for t in transformations:
            if t['sublibrary'] not in seen:
                seen[t['sublibrary']] = t['group_id']
        for sub, lib in sorted(seen.items()):
            print(f"    {sub} -> {lib}")

        # File suffix breakdown
        suffixes = defaultdict(int)
        for t in transformations:
            suffixes[t['suffix']] += 1
        print(f"\n  File suffixes:")
        for suf, cnt in sorted(suffixes.items(), key=lambda x: -x[1]):
            print(f"    {suf}: {cnt}")

    # Per-library file completeness
    print(f"\n" + "-" * 40)
    print("SOP COMPLIANCE CHECK")
    print("-" * 40)

    SCALE_EXPECTED_SUFFIXES = {
        'GEX':        {'.cram', '.csv', '.json', '_trimmer-stats.csv',
                       '-trimmer-stats.csv', '_trimmer-failure_codes.csv',
                       '-trimmer-failure-codes.csv'},
        'Hash_oligo': {'.cram', '.csv', '.json', '-trimmer-stats.csv',
                       '-trimmer-failure-codes.csv'},
    }

    if files_per_lib:
        print(f"\nLibrary/assay combinations found: {len(files_per_lib)}")
        for (runid, lib, assay), found in sorted(files_per_lib.items()):
            print(f"  {runid} / {lib} / {assay}: {len(found)} file type(s)")
            for s in sorted(found):
                print(f"    {s}")
    else:
        print("  No sample files parsed successfully.")

    print(f"\n" + "-" * 40)
    print("MAPPING LOGIC DISCOVERED:")
    print("-" * 40)
    print(f"""
S3 path: s3://czi-novogene/trapnell-seahub-bcp/{{order}}/{{ExperimentID}}/raw/{{RunID}}/
S3 file (SOP form):       {{RunID}}-{{GroupID}}_{{assay}}_QSR-{{N}}[-SCALEPLEX]_{{well}}{{suffix}}
S3 file (index-direct):   {{RunID}}-{{GroupID}}_{{assay}}_{{IndexSequence}}{{suffix}}

LOCAL: {{base}}/{{RunID}}-{{date}}/{{RunID}}-QSR-{{N}}[-SCALEPLEX]-{{index}}/
           {{RunID}}-QSR-{{N}}[-SCALEPLEX]-{{index}}{{suffix}}

Key distinctions:
  ExperimentID (S3 path directory): chip/experiment container (e.g. CHEM13-R096)
  GroupID (S3 filename):            per-sublibrary identifier (e.g. R096G)

Key mappings (via SIF):
  SOP form:          S3 QSR-{{N}} -> SIF sublibrary -> GroupID validation
  Index-direct form: S3 IndexSequence -> SIF index -> GroupID validation
  Local sublibrary QSR-{{N}}           -> normalized QSR{{N}}          -> SIF -> GroupID
  Local sublibrary QSR-{{N}}-SCALEPLEX -> normalized QSR{{N}}SCALEPLEX -> SIF -> GroupID

Suffix normalization (S3 -> local):
  -trimmer-stats.csv         -> _trimmer-stats.csv
  -trimmer-failure-codes.csv -> _trimmer-failure_codes.csv
""")

    return outliers


def main():
    # Parse arguments
    args = sys.argv[1:]
    
    if len(args) < 1:
        print("Usage: python analyze_mapping.py <mapping_file> [identifiers] [--mode MODE]")
        print()
        print("Modes:")
        print("  raw-experiment (default): For raw files with experiment names like GENE7, CHEM16")
        print("  raw-groupid: For raw files with group IDs like CUIMC_003_003")
        print("  processed: For processed files only (cellranger output, etc.)")
        print("  mixed: For files with both raw and processed entries")
        print("  scale: For Scale scRNA-seq orders (trapnell-seahub-bcp bucket)")
        print()
        print("For processed and mixed modes, identifiers are optional - auto-extracted from S3 paths.")
        print()
        print("Examples:")
        print("  python analyze_mapping.py mapping.txt GENE7  # raw-experiment mode")
        print("  python analyze_mapping.py mapping.txt REF3,GENE7")
        print("  python analyze_mapping.py mapping.txt CUIMC_003_003 --mode raw-groupid")
        print("  python analyze_mapping.py mapping.txt --mode processed  # Auto-extract group IDs")
        print("  python analyze_mapping.py mapping.txt --mode mixed  # Analyze both raw and processed")
        print("  python analyze_mapping.py mapping.txt --mode scale --sif SIF.csv  # Scale orders")
        sys.exit(1)
    
    filepath = args[0]
    
    # Check for --mode flag
    mode = 'raw-experiment'
    if '--mode' in args:
        mode_idx = args.index('--mode')
        if mode_idx + 1 < len(args):
            mode = args[mode_idx + 1]
    
    # Legacy mode names for backwards compatibility
    if mode == 'experiment':
        mode = 'raw-experiment'
    elif mode == 'groupid':
        mode = 'raw-groupid'
    
    # Parse --sif flag (required for scale mode)
    sif_filepath = None
    if '--sif' in args:
        sif_idx = args.index('--sif')
        if sif_idx + 1 < len(args):
            sif_filepath = args[sif_idx + 1]

    # Validate SIF requirement for scale mode
    if mode == 'scale' and sif_filepath is None:
        print("ERROR: --sif <sif_file> is required for scale mode")
        print("Usage: python analyze_mapping.py <mapping_file> --mode scale --sif <sif_file>")
        sys.exit(1)

    # Parse identifiers - optional for processed and mixed modes
    identifiers = None
    if len(args) >= 2 and args[1] != '--mode':
        identifiers = [e.strip() for e in args[1].split(',')]
    
    # Validate identifiers requirement based on mode
    if mode not in ['processed', 'mixed', 'scale'] and identifiers is None:
        print(f"ERROR: identifiers required for mode '{mode}'")
        print(f"Usage: python analyze_mapping.py <mapping_file> <identifiers> --mode {mode}")
        sys.exit(1)
    
    print("Analyzing mapping file:", filepath)
    print(f"Mode: {mode}")
    if identifiers:
        print(f"Identifiers: {identifiers}")
    else:
        print(f"Identifiers: Will be auto-extracted from S3 paths")
    print()
    
    mappings = parse_mapping_file(filepath)
    print(f"Loaded {len(mappings)} mappings")
    
    # Handle different modes
    if mode == 'mixed':
        # Separate mappings into raw and processed
        raw_mappings = []
        processed_mappings = []
        
        for s3_path, local_path, line_num in mappings:
            if '/processed/' in s3_path:
                processed_mappings.append((s3_path, local_path, line_num))
            elif '/raw/' in s3_path:
                raw_mappings.append((s3_path, local_path, line_num))
            else:
                # Order-level files - could go either way, let's track separately
                raw_mappings.append((s3_path, local_path, line_num))
        
        print(f"\nSeparated into:")
        print(f"  Raw files: {len(raw_mappings)}")
        print(f"  Processed files: {len(processed_mappings)}")
        
        # Check uniqueness for all mappings
        is_unique = check_uniqueness(mappings, outs_only=False)
        
        # Analyze processed files first (if any)
        processed_outliers = []
        if processed_mappings:
            processed_outliers = analyze_mapping_logic_processed(processed_mappings, identifiers)
        
        # Analyze raw files (if any)
        raw_outliers = []
        if raw_mappings:
            # Try to auto-detect whether raw files use group IDs or experiment names
            # For now, default to raw-groupid since we're auto-extracting
            raw_outliers = analyze_mapping_logic_groupid(raw_mappings, identifiers or extract_group_ids_from_mappings(raw_mappings))
        
        # Combine outliers
        outliers = processed_outliers + raw_outliers
        
    elif mode == 'raw-groupid':
        is_unique = check_uniqueness(mappings, outs_only=False)
        outliers = analyze_mapping_logic_groupid(mappings, identifiers)
    elif mode == 'processed':
        is_unique = check_uniqueness(mappings, outs_only=True)
        outliers = analyze_mapping_logic_processed(mappings, identifiers)
    elif mode == 'scale':
        is_unique = check_uniqueness(mappings, outs_only=False)
        outliers = analyze_mapping_logic_scale(mappings, sif_filepath)
    else:  # raw-experiment
        is_unique = check_uniqueness(mappings, outs_only=False)
        outliers = analyze_mapping_logic(mappings, identifiers)
    
    report_outliers(outliers)
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total mappings: {len(mappings)}")
    print(f"Uniqueness check: {'PASS' if is_unique else 'FAIL'}")
    print(f"Outliers: {len(outliers)}")


if __name__ == "__main__":
    main()