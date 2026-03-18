from __future__ import annotations

import re

from typing import Final

# ---------------------------------------------------------------------------
# Assay types per assay family (provider-agnostic)
# ---------------------------------------------------------------------------

ASSAYS_10X: Final[set[str]] = {"GEX", "CRI", "ATAC"}
ASSAYS_SCALE: Final[set[str]] = {"GEX", "hash_oligo", "GEX_hash_oligo"}
ASSAYS_SCI: Final[set[str]] = {"GEX", "hash_oligo", "GEX_hash_oligo"}

ASSAYS_BY_FAMILY: Final[dict[str, set[str]]] = {
    "10x": ASSAYS_10X,
    "scale": ASSAYS_SCALE,
    "sci": ASSAYS_SCI,
}

# Provider-specific extra assays layered on top of the family base set.
# Only providers that add assays beyond the family base need entries here.
PROVIDER_EXTRA_ASSAYS: Final[dict[str, dict[str, set[str]]]] = {
    "10x": {
        "psomagen": {"viral_ORF"},
    },
}

# Canonical (SOP) spellings – lower-cased key → canonical form
CANONICAL_ASSAY: Final[dict[str, str]] = {
    "gex": "GEX",
    "cri": "CRI",
    "atac": "ATAC",
    "hash_oligo": "hash_oligo",
    "gex_hash_oligo": "GEX_hash_oligo",
    "viral_orf": "viral_ORF",
}

# ---------------------------------------------------------------------------
# Provider configuration
# ---------------------------------------------------------------------------

PROVIDER_ORDER_RE: Final[dict[str, str]] = {
    "novogene": r"NVUS\d+-\d+(?:-\d+)*",
    "psomagen": r"AN\d+",
}

PROVIDER_BUCKET: Final[dict[str, str]] = {
    "novogene": "czi-novogene",
    "psomagen": "czi-psomagen",
}

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


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def get_assays(assay_family: str, provider: str | None = None) -> set[str]:
    """Return valid assay types for an assay family, optionally extended by provider."""
    base = ASSAYS_BY_FAMILY.get(assay_family)
    if base is None:
        raise ValueError(
            f"Unknown assay family '{assay_family}', "
            f"expected one of {sorted(ASSAYS_BY_FAMILY)}"
        )
    assays = set(base)
    if provider:
        extras = PROVIDER_EXTRA_ASSAYS.get(assay_family, {}).get(
            provider.lower(), set()
        )
        assays |= extras
    return assays


def build_assay_regex(assays: set[str]) -> str:
    """Build a case-insensitive alternation regex from a set of assay names.

    Sorted longest-first to avoid partial matches (e.g. GEX_hash_oligo
    must be tried before hash_oligo or GEX).
    """
    sorted_assays = sorted(assays, key=len, reverse=True)
    return r"(?i:" + "|".join(sorted_assays) + r")"


def get_order_pattern(provider: str) -> str:
    """Return the order-number regex for a provider."""
    provider = provider.lower()
    pattern = PROVIDER_ORDER_RE.get(provider)
    if pattern is None:
        raise ValueError(
            f"Unknown provider '{provider}', expected one of {sorted(PROVIDER_ORDER_RE)}"
        )
    return pattern


# ---------------------------------------------------------------------------
# Backward-compatible aliases (derived from the new API)
# ---------------------------------------------------------------------------

ASSAYS_10X_NOVOGENE: Final[set[str]] = get_assays("10x", "novogene")
ASSAYS_10X_PSOMAGEN: Final[set[str]] = get_assays("10x", "psomagen")


__all__ = [
    # Per-family assay sets
    "ASSAYS_10X",
    "ASSAYS_SCALE",
    "ASSAYS_SCI",
    "ASSAYS_BY_FAMILY",
    "PROVIDER_EXTRA_ASSAYS",
    # Canonical spellings
    "CANONICAL_ASSAY",
    # Provider config
    "PROVIDER_ORDER_RE",
    "PROVIDER_BUCKET",
    # Metadata
    "_RUN_METADATA_RE",
    # Helpers
    "get_assays",
    "build_assay_regex",
    "get_order_pattern",
    # Backward-compatible aliases
    "ASSAYS_10X_NOVOGENE",
    "ASSAYS_10X_PSOMAGEN",
]
