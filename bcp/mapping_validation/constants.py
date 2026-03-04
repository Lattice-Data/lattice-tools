from __future__ import annotations

import re

from typing import Final

# ---------------------------------------------------------------------------
# SOP-defined assay types per assay family
# ---------------------------------------------------------------------------

ASSAYS_10X_NOVOGENE: Final[set[str]] = {"GEX", "CRI", "ATAC"}
ASSAYS_10X_PSOMAGEN: Final[set[str]] = {"GEX", "CRI", "ATAC", "viral_ORF"}
ASSAYS_SCALE: Final[set[str]] = {"GEX", "hash_oligo", "GEX_hash_oligo"}
ASSAYS_SCI: Final[set[str]] = {"GEX", "hash_oligo", "GEX_hash_oligo"}

# Canonical (SOP) spellings – lower-cased key → canonical form
CANONICAL_ASSAY: Final[dict[str, str]] = {
    "gex": "GEX",
    "cri": "CRI",
    "atac": "ATAC",
    "hash_oligo": "hash_oligo",
    "gex_hash_oligo": "GEX_hash_oligo",
    "viral_orf": "viral_ORF",
}

# Regex alternation fragments used to anchor assay names inside S3 path
# regexes. Longest alternatives come first to avoid partial matches.
# Case-insensitive matching via (?i:...) so that ANY casing variant is
# captured (and then checked against the canonical SOP spelling).
_SCALE_ASSAY_RE: Final[str] = r"(?i:GEX_hash_oligo|hash_oligo|GEX)"
_10X_NOVOGENE_ASSAY_RE: Final[str] = r"(?i:ATAC|CRI|GEX)"
_10X_PSOMAGEN_ASSAY_RE: Final[str] = r"(?i:viral_ORF|ATAC|CRI|GEX)"

# Order number patterns (Novogene allows sub-order suffixes like -28-4)
_NOVOGENE_ORDER_RE: Final[str] = r"NVUS\d+-\d+(?:-\d+)*"
_PSOMAGEN_ORDER_RE: Final[str] = r"AN\d+"

# Run-level metadata filenames expected alongside sample data.
# The optional \d+_ prefix handles files placed at the order level
# with a RunID prefix (e.g. 439844_UploadCompleted.json).
_RUN_METADATA_RE: Final[re.Pattern[str]] = re.compile(
    r"^(?:\d+_)?(?:"
    r"LibraryInfo\.xml"
    r"|SequencingInfo\.json"
    r"|UploadCompleted\.json"
    r"|merged_trimmer-(?:failure_codes|stats)\.csv"
    r"|run_(?:SecondaryAnalysis|VariantCalling)\.txt"
    r")$"
)


__all__ = [
    "ASSAYS_10X_NOVOGENE",
    "ASSAYS_10X_PSOMAGEN",
    "ASSAYS_SCALE",
    "ASSAYS_SCI",
    "CANONICAL_ASSAY",
    "_SCALE_ASSAY_RE",
    "_10X_NOVOGENE_ASSAY_RE",
    "_10X_PSOMAGEN_ASSAY_RE",
    "_NOVOGENE_ORDER_RE",
    "_PSOMAGEN_ORDER_RE",
    "_RUN_METADATA_RE",
]

