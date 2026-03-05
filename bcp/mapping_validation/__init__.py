from __future__ import annotations

from .constants import (
    ASSAYS_10X,
    ASSAYS_10X_NOVOGENE,
    ASSAYS_10X_PSOMAGEN,
    ASSAYS_BY_FAMILY,
    ASSAYS_SCALE,
    ASSAYS_SCI,
    CANONICAL_ASSAY,
    PROVIDER_BUCKET,
    PROVIDER_EXTRA_ASSAYS,
    PROVIDER_ORDER_RE,
    build_assay_regex,
    get_assays,
    get_order_pattern,
)
from .parsing import MappingRow, parse_mapping_file
from .sif_io import (
    _normalize_sif_groupid,
    load_sif_group_assays,
    load_sif_library_assays,
    load_sif_scale_group_assays,
    load_sif_scale_groupids,
)
from .uniqueness import validate_uniqueness
from .validators import (
    _is_run_metadata,
    compare_groupid_assays,
    find_unmatched_sif_paths_10x,
    validate_library_assay_consistency,
    validate_local_paths_scale_raw,
    validate_s3_10x_raw,
    validate_s3_local_consistency_scale,
    validate_s3_scale_raw,
    validate_sif_completeness_scale,
)
from .cli import main

__all__ = [
    # Per-family assay sets
    "ASSAYS_10X",
    "ASSAYS_SCALE",
    "ASSAYS_SCI",
    "ASSAYS_BY_FAMILY",
    "PROVIDER_EXTRA_ASSAYS",
    # Backward-compatible aliases
    "ASSAYS_10X_NOVOGENE",
    "ASSAYS_10X_PSOMAGEN",
    # Canonical spellings
    "CANONICAL_ASSAY",
    # Provider config
    "PROVIDER_ORDER_RE",
    "PROVIDER_BUCKET",
    # Helpers
    "get_assays",
    "build_assay_regex",
    "get_order_pattern",
    # Core types / parsing
    "MappingRow",
    "parse_mapping_file",
    # Uniqueness
    "validate_uniqueness",
    # Validators
    "compare_groupid_assays",
    "validate_s3_10x_raw",
    "validate_local_paths_scale_raw",
    "validate_s3_scale_raw",
    "validate_sif_completeness_scale",
    "validate_s3_local_consistency_scale",
    "validate_library_assay_consistency",
    "find_unmatched_sif_paths_10x",
    "_is_run_metadata",
    # SIF helpers
    "_normalize_sif_groupid",
    "load_sif_group_assays",
    "load_sif_scale_group_assays",
    "load_sif_scale_groupids",
    "load_sif_library_assays",
    # CLI
    "main",
]
