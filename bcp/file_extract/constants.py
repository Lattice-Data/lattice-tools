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

H5_BASE_COLUMNS = [
    "library",
    "sample",
    "s3_uri",
    "size_bytes",
    "crc64nvme_base64",
]

H5_INTROSPECT_COLUMNS = [
    "observation_count",
    "feature_counts",
    "feature_count_total",
]

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
