"""
Constants for QA pipeline: chemistries, assays, expected cellranger and raw file patterns.
"""

import re

chemistries = {
    "Single Cell 5' R2-only v3": "5p",
    "Single Cell 5' R2-only": "5p",
    "Single Cell 3' v4 (polyA)": "3p",
    "Single Cell 3' v3": "3p",
    "Flex Gene Expression": "flex",
    "GEM-X Flex v2": "flex",
}

valid_assays = ["CRI", "GEX", "ATAC", "viral_ORF", "GEX_hash_oligo", "hash_oligo"]

# Raw pipeline types supported by qa.ipynb / qa_checks (normalize casing via qa_mods.normalize_raw_assay)
ALLOWED_RAW_ASSAYS = frozenset(
    (
        "10x",
        "10x_cram",
        "10x_viral_ORF",
        "sci_jumbo",
        "sci_plex",
        "scale",
        "seahub_sci",
    )
)

# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-flex-outputs-frp
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview
cellranger_ignore = [
    "analysis",
    "cloupe",
    "cell_types",
]

cellranger_expected = {
    "cellranger-9.0.1": {
        "nonflex": {
            "outs": [
                "config.csv",
                "multi/count/raw_feature_bc_matrix.h5",
                "multi/count/raw_feature_bc_matrix.tar.gz",
                "multi/count/raw_molecule_info.h5",
                "multi/count/unassigned_alignments.bam",
                "multi/count/unassigned_alignments.bam.bai",
            ],
            "per_sample": [
                "count/sample_alignments.bam",
                "count/sample_alignments.bam.bai",
                "count/sample_filtered_barcodes.csv",
                "count/sample_filtered_feature_bc_matrix.h5",
                "count/sample_filtered_feature_bc_matrix.tar.gz",
                "count/sample_molecule_info.h5",
                "metrics_summary.csv",
                "web_summary.html",
            ],
        },
        "flex": {
            "outs": [
                "config.csv",
                "multi/count/raw_feature_bc_matrix.h5",
                "multi/count/raw_feature_bc_matrix.tar.gz",
                "multi/count/raw_molecule_info.h5",
                "multi/count/raw_probe_bc_matrix.h5",
            ],
            "per_sample": [
                "count/probe_set.csv",
                "count/sample_filtered_barcodes.csv",
                "count/sample_filtered_feature_bc_matrix.h5",
                "count/sample_filtered_feature_bc_matrix.tar.gz",
                "count/sample_molecule_info.h5",
                "count/sample_raw_feature_bc_matrix.h5",
                "count/sample_raw_feature_bc_matrix.tar.gz",
                "count/sample_raw_probe_bc_matrix.h5",
                "metrics_summary.csv",
                "web_summary.html",
            ],
        },
    },
    "cellranger-10.0.0": {
        "nonflex": {
            "outs": [
                "config.csv",
                "filtered_feature_bc_matrix/barcodes.tsv.gz",
                "filtered_feature_bc_matrix/features.tsv.gz",
                "filtered_feature_bc_matrix/matrix.mtx.gz",
                "filtered_feature_bc_matrix.h5",
                "qc_library_metrics.csv",
                "qc_report.html",
                "qc_sample_metrics.csv",
                "raw_feature_bc_matrix/barcodes.tsv.gz",
                "raw_feature_bc_matrix/features.tsv.gz",
                "raw_feature_bc_matrix/matrix.mtx.gz",
                "raw_feature_bc_matrix.h5",
                "raw_molecule_info.h5",
            ],
            "per_sample": [
                "sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                "sample_filtered_feature_bc_matrix/features.tsv.gz",
                "sample_filtered_feature_bc_matrix/matrix.mtx.gz",
                "sample_raw_feature_bc_matrix/barcodes.tsv.gz",
                "sample_raw_feature_bc_matrix/features.tsv.gz",
                "sample_raw_feature_bc_matrix/matrix.mtx.gz",
                "metrics_summary.csv",
                "sample_filtered_barcodes.csv",
                "sample_filtered_feature_bc_matrix.h5",
                "sample_molecule_info.h5",
                "sample_raw_feature_bc_matrix.h5",
                "web_summary.html",
            ],
        },
        "flex": {
            "outs": [
                "config.csv",
                "filtered_feature_bc_matrix/barcodes.tsv.gz",
                "filtered_feature_bc_matrix/features.tsv.gz",
                "filtered_feature_bc_matrix/matrix.mtx.gz",
                "filtered_feature_bc_matrix.h5",
                "probe_set.csv",
                "qc_library_metrics.csv",
                "qc_report.html",
                "qc_sample_metrics.csv",
                "raw_feature_bc_matrix/barcodes.tsv.gz",
                "raw_feature_bc_matrix/features.tsv.gz",
                "raw_feature_bc_matrix/matrix.mtx.gz",
                "raw_feature_bc_matrix.h5",
                "raw_molecule_info.h5",
                "raw_probe_bc_matrix.h5",
            ],
            "per_sample": [
                "metrics_summary.csv",
                "sample_filtered_barcodes.csv",
                "sample_filtered_feature_bc_matrix.h5",
                "sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                "sample_filtered_feature_bc_matrix/features.tsv.gz",
                "sample_filtered_feature_bc_matrix/matrix.mtx.gz",
                "sample_molecule_info.h5",
                "sample_raw_feature_bc_matrix.h5",
                "sample_raw_feature_bc_matrix/barcodes.tsv.gz",
                "sample_raw_feature_bc_matrix/features.tsv.gz",
                "sample_raw_feature_bc_matrix/matrix.mtx.gz",
                "sample_raw_probe_bc_matrix.h5",
                "web_summary.html",
            ],
        },
    },
    "count": {
        "outs": [
            "filtered_feature_bc_matrix/barcodes.tsv.gz",
            "filtered_feature_bc_matrix/features.tsv.gz",
            "filtered_feature_bc_matrix/matrix.mtx.gz",
            "filtered_feature_bc_matrix.h5",
            "metrics_summary.csv",
            "molecule_info.h5",
            "possorted_genome_bam.bam",
            "possorted_genome_bam.bam.bai",
            "raw_feature_bc_matrix/barcodes.tsv.gz",
            "raw_feature_bc_matrix/features.tsv.gz",
            "raw_feature_bc_matrix/matrix.mtx.gz",
            "raw_feature_bc_matrix.h5",
            "web_summary.html",
        ]
    },
}

raw_expected = {
    "sci_jumbo": [
        ".cram",
        ".cram-metadata.json",
        ".csv",
        ".json",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
        "_FlowQ.metric",
        "_SNVQ.metric",
    ],
    "sci_plex": [
        ".cram",
        ".cram-metadata.json",
        ".csv",
        ".json",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
    ],
    "10x": [
        ".csv",
        ".json",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
        "_unmatched.cram",
        "_unmatched.cram-metadata.json",
        "_unmatched.csv",
        "_unmatched.json",
        "_S1_L001_R1_001.csv",
        "_S1_L001_R1_001.fastq.gz",
        "_S1_L001_R1_001.fastq.gz-metadata.json",
        "_S1_L001_R1_001.json",
        "_S1_L001_R1_001_sample.fastq.gz",
        "_S1_L001_R1_001_sample.fastq.gz-metadata.json",
        "_S1_L001_R2_001.csv",
        "_S1_L001_R2_001.fastq.gz",
        "_S1_L001_R2_001.fastq.gz-metadata.json",
        "_S1_L001_R2_001.json",
        "_S1_L001_R2_001_sample.fastq.gz",
        "_S1_L001_R2_001_sample.fastq.gz-metadata.json",
    ],
    "10x_viral_ORF": [
        ".csv",
        ".json",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
        ".cram",
        ".cram-metadata.json",
        "_FlowQ.metric",
        "_SNVQ.metric",
    ],
    "10x_cram": [
        ".cram",
        ".cram-metadata.json",
        ".csv",
        ".json",
        "_FlowQ.metric",
        "_SNVQ.metric",
        "_trimmer-failure_codes.csv",
        "_trimmer-stats.csv",
    ],
    "scale": [],
    "seahub_sci": [
        ".trim.cram",
        ".trim.csv",
        ".trim.stderr",
        ".trim.stdout",
        ".trim_fail.csv",
    ],
}

raw_optional = {
    "10x": [
        ".scRNA.applicationQC.h5",
        ".scRNA.applicationQC.html",
        "_Log.final.out",
        "_Log.out",
        "_Log.progress.out",
        "_ReadsPerGene.out.tab",
        "_SJ.out.tab",
    ],
    "10x_cram": [
        "_extract_stats.h5",
    ],
    "scale": [],
    "seahub_sci": [
        ".trim.cram-metadata.json",
    ],
}

# ---------------------------------------------------------------------------
# SeaHub lab raw upload patterns (trapnell / hamazaki *-seahub-bcp buckets)
#
# SOP layout:
#   s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/{sublibrary}/{wafer}/
# SOP filename:
#   {wafer}-{sublibrary}[_{well}]_{sublibrary type}-{UG}-{barcode}.trim.*
#
# Known-good examples:
#   .../hamazaki-seahub-bcp/CHEM3-R100/raw/R100E/441389/
#       441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram
#   .../trapnell-seahub-bcp/REF3/raw/REF3_P05_2/436830/
#       436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram
#
# The ExperimentID is not a filename field: it appears in the trapnell example
# only because it is part of that project's sublibrary name, and is absent from
# the hamazaki example.  ExperimentID may contain hyphens (``CHEM3-R100``), so
# path segments must never be hyphen-split.
# ---------------------------------------------------------------------------

# Sublibrary types valid for SeaHub, narrower than ``valid_assays`` (which also
# carries ATAC and viral_ORF).  Ordered longest-first so alternation in
# SEAHUB_STEM_RE prefers ``GEX_hash_oligo`` over ``GEX`` / ``hash_oligo``.
SEAHUB_SUBLIBRARY_TYPES: tuple[str, ...] = (
    "GEX_hash_oligo",
    "hash_oligo",
    "GEX",
    "CRI",
)

# Both suffix tuples are ordered longest-first so ``.cram-metadata.json`` is
# never truncated by ``.cram``.
SEAHUB_TRIM_SUFFIXES: tuple[str, ...] = (
    ".trim.cram-metadata.json",
    ".trim_fail.csv",
    ".trim.cram",
    ".trim.csv",
    ".trim.stderr",
    ".trim.stdout",
)

# Non-compliant family seen in real uploads: the same six trim artifacts with
# the ``.trim`` infix dropped.  Recognised so completeness still runs against
# what was actually delivered (which is how a genuinely absent CRAM surfaces),
# then reported as a ``missing_trim_infix`` SOP violation.
SEAHUB_BARE_SUFFIXES: tuple[str, ...] = (
    ".cram-metadata.json",
    "_fail.csv",
    ".cram",
    ".csv",
    ".stderr",
    ".stdout",
)

# Maps a bare-family suffix to the SOP suffix it should have carried.
SEAHUB_BARE_TO_TRIM_SUFFIX: dict[str, str] = {
    ".cram-metadata.json": ".trim.cram-metadata.json",
    "_fail.csv": ".trim_fail.csv",
    ".cram": ".trim.cram",
    ".csv": ".trim.csv",
    ".stderr": ".trim.stderr",
    ".stdout": ".trim.stdout",
}

SEAHUB_TRIM_TO_BARE_SUFFIX: dict[str, str] = {
    trim: bare for bare, trim in SEAHUB_BARE_TO_TRIM_SUFFIX.items()
}

# Bare-family counterparts of raw_expected / raw_optional["seahub_sci"], derived
# so the required-artifact set has a single source of truth.
SEAHUB_BARE_EXPECTED: tuple[str, ...] = tuple(
    SEAHUB_TRIM_TO_BARE_SUFFIX[ending] for ending in raw_expected["seahub_sci"]
)
SEAHUB_BARE_OPTIONAL: tuple[str, ...] = tuple(
    SEAHUB_TRIM_TO_BARE_SUFFIX[ending] for ending in raw_optional["seahub_sci"]
)

# Suffixes that carry per-well trimmer failure counts, in either family.
SEAHUB_FAIL_SUFFIXES: tuple[str, ...] = (".trim_fail.csv", "_fail.csv")

# fastq_log bucket for wells whose filename omits the sublibrary type, so they
# stay visible in per-sublibrary counts instead of dropping out.  The missing
# type itself is reported as an invalid_sublibrary_type SOP violation.
SEAHUB_UNTYPED_LABEL = "UNTYPED"

_SEAHUB_TYPE_ALT = "|".join(SEAHUB_SUBLIBRARY_TYPES)

# Group is non-greedy so ``R100E_GEX_hash_oligo`` parses as group ``R100E``
# plus type ``GEX_hash_oligo``; a greedy group would instead yield group
# ``R100E_GEX`` plus type ``hash_oligo``.
SEAHUB_STEM_RE = re.compile(
    rf"^(?P<wafer>\d{{6,8}})-(?P<group>.+?)"
    rf"_(?P<assay>{_SEAHUB_TYPE_ALT})"
    r"-(?P<ug>Z\d{4})-(?P<barcode>[ACGT]+)$"
)

# Relaxed variant for stems that omit the sublibrary type entirely, so wafer /
# UG / barcode stay recoverable for reporting and cross-bucket matching.
SEAHUB_STEM_NO_TYPE_RE = re.compile(
    r"^(?P<wafer>\d{6,8})-(?P<group>.+?)"
    r"-(?P<ug>Z\d{4})-(?P<barcode>[ACGT]+)$"
)

SEAHUB_UG_RE = re.compile(r"^Z\d{4}$")
SEAHUB_BARCODE_RE = re.compile(r"^[ACGT]+$")
SEAHUB_WELL_RE = re.compile(r"^[A-H]\d{1,2}$")

SEAHUB_PLATE_SIZES = frozenset({48, 96})

# ---------------------------------------------------------------------------
# Scale raw file patterns
#
# Scale filenames do not follow the "{beginning}{suffix}" convention used by
# 10x / sci assays, so we validate them via regex in check_extra_raw_files
# instead of the raw_expected suffix mechanism.
# ---------------------------------------------------------------------------

# Per-RT file: well position (row 1-12, column A-H) before the extension.
# Production uses underscore (e.g. _5B.cram); some paths use hyphen (-12H.csv).
SCALE_RT_FILE_RE = re.compile(r"[_-]\d{1,2}[A-H]\.(cram|csv|json)$")

# Aggregate file produced per sublibrary: trimmer stats, trimmer failure codes,
# and unmatched reads. Prefix may be _ or -; failure codes file may use
# trimmer-failure-codes or trimmer-failure_codes (underscore before "codes").
SCALE_AGGREGATE_FILE_RE = re.compile(
    r"[_-](trimmer-failure-codes\.csv|trimmer-failure_codes\.csv|trimmer-stats\.csv"
    r"|unmatched\.(cram|csv|json))$"
)

# Relaxed variant: additionally accepts the truncated alias ``..._stats.csv``
# / ``...-stats.csv`` (collaborator-generated files missing the ``trimmer-``
# token). Used only when QARunContext.allow_truncated_stats_name=True.
SCALE_AGGREGATE_FILE_RELAXED_RE = re.compile(
    r"[_-](trimmer-failure-codes\.csv|trimmer-failure_codes\.csv|trimmer-stats\.csv"
    r"|stats\.csv|unmatched\.(cram|csv|json))$"
)

# Wafer-level housekeeping files (no sublibrary / assay prefix).
SCALE_WAFER_MISC_RE = re.compile(
    r"^(\d{6,8}_(SequencingInfo\.json|LibraryInfo\.xml)"
    r"|merged_trimmer-(failure_codes|stats)\.csv)$"
)

SCALE_WORKFLOW_REQUIRED_PARAMS = {
    "bamOut": "true",
    "scalePlex": "true",
    "scalePlexAssignmentMethod": "fc",
}

SCALE_SAMPLES_FORBIDDEN_COLUMNS = frozenset({"scalePlexBarcodes"})
