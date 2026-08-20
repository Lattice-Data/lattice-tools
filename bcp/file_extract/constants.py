from __future__ import annotations

DEFAULT_H5_TARGET_FILENAME = "sample_filtered_feature_bc_matrix.h5"

FASTQ_COLUMNS = [
    "filename",
    "s3_uri",
    "read",
    "lane",
    "size_bytes",
    "crc64nvme_base64",
    "read_count",
    "crc_error",
    "metadata_error",
]

CRAM_COLUMNS = [
    "filename",
    "s3_uri",
    "size_bytes",
    "crc64nvme_base64",
    "read_count",
    "crc_error",
    "metadata_error",
]

# Appended to the diagnostic TSVs so the blank SequenceFileSet.library column can
# be filled by lookup: sample_dir keys the CRO group, set_alias keys the sheet row.
SHEET_HELPER_COLUMNS = [
    "sample_dir",
    "set_alias",
]

# Column order of the Lattice SequenceFile submission sheet, verbatim. Rows are
# projected onto this list by name, so a drifting sheet is a one-line change here.
SEQUENCE_FILE_COLUMNS = [
    "#response",
    "#response_time",
    "uuid",
    "aliases",
    "lab",
    "file_format",
    "s3_uri",
    "crc64nvme_base64",
    "no_file_available",
    "derived_from",
    "status",
    "submitter_comment",
    "description",
    "read_count",
    "file_size",
]

# Column order of the Lattice SequenceFileSet submission sheet, verbatim.
SEQUENCE_FILE_SET_COLUMNS = [
    "#response",
    "#response_time",
    "uuid",
    "aliases",
    "lab",
    "library",
    "run_cardinality",
    "status",
    "submitter_comment",
    "description",
    "CRO_order",
    "is_pilot_order",
    "read1",
    "read2",
    "read3",
    "index1",
    "index2",
    "untrimmed_cram",
    "trimmed_cram",
    "sequencing_platform",
]

# Read designator -> SequenceFileSet slot column.
READ_SLOT_COLUMNS = {
    "R1": "read1",
    "R2": "read2",
    "R3": "read3",
    "I1": "index1",
    "I2": "index2",
}

# Slot fill order within a set; also the sort order of sheet rows within a set.
READ_SLOT_ORDER = ["R1", "R2", "R3", "I1", "I2"]

CRAM_SLOT_COLUMNS = {
    "untrimmed": "untrimmed_cram",
    "trimmed": "trimmed_cram",
}

FILE_FORMAT_FASTQ = "fastq"
FILE_FORMAT_CRAM = "cram"

RUN_CARDINALITY_PAIRED = "paired-end"
RUN_CARDINALITY_SINGLE = "single-end"

STATUS_CURRENT = "current"
NO_FILE_AVAILABLE_FALSE = "FALSE"

H5_BASE_COLUMNS = [
    "library",
    "sample",
    "s3_uri",
    "size_bytes",
    "crc64nvme_base64",
]

SCALE_H5AD_COLUMNS = [
    "filename",
    "s3_uri",
    "crc64nvme_base64",
    "sample",
    "samples",
    "file_size",
    "observation_count",
    "feature_counts",
]

H5_INTROSPECT_COLUMNS = [
    "observation_count",
    "feature_counts",
    "feature_count_total",
    "unmapped_feature_types",
]

# Cell Ranger matrix/features/feature_type -> Lattice feature_counts enum.
CR_TO_LATTICE_FEATURE_TYPE = {
    "Gene Expression": "gene",
    "CRISPR Guide Capture": "guide capture",
    "Antibody Capture": "antibody capture",
    "Peaks": "peak",
}

H5_GENOME_COLUMN = "gene_counts_by_genome"

H5_METRICS_COLUMNS = [
    "metrics_cells",
    "metrics_cells_match",
]

TRANSIENT_ERROR_CODES = frozenset(
    {
        "Throttling",
        "ThrottlingException",
        "RequestThrottled",
        "SlowDown",
        "RequestTimeout",
        "RequestTimeoutException",
        "RequestTimeTooSkewed",
        "InternalError",
        "InternalServerError",
        "ServiceUnavailable",
        "ServiceUnavailableException",
        "500",
        "502",
        "503",
        "504",
    }
)

TRANSIENT_HTTP_STATUSES = frozenset({429, 500, 502, 503, 504})
