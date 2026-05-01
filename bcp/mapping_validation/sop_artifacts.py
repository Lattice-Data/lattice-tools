"""
Single source of truth for SOP-required per-prefix raw file suffixes.

Each assay family (10x raw FASTQ, 10x_cram, sci, scale) has a fixed set of
artifacts the SOP requires per prefix.  Validators in this package classify
the suffix that follows the prefix against these tables so that:

- truncated S3 filenames (e.g. ``_codes.csv`` instead of
  ``_trimmer-failure_codes.csv``) are flagged as ``unexpected_sample_file``
- missing artifacts per prefix are flagged as ``missing_sample_artifacts``

Artifact keys (e.g. ``trimmer_failure_codes``, ``cram``, ``sample_fastq_r1``)
are stable identifiers used in error reports and inter-validator comparisons.
"""

from __future__ import annotations

import re
from typing import Final


# ---------------------------------------------------------------------------
# 10x_cram (per-sample CRAM bundle, per RunID metadata)
# ---------------------------------------------------------------------------

TENX_CRAM_SAMPLE_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    ".cram": "cram",
    ".csv": "csv",
    ".json": "json",
    "_extract_stats.h5": "extract_stats",
    "_SNVQ.metric": "snvq",
    "_FlowQ.metric": "flowq",
    "_trimmer-stats.csv": "trimmer_stats",
    "_trimmer-failure_codes.csv": "trimmer_failure_codes",
    "_trimmer-failure-codes.csv": "trimmer_failure_codes",
}

TENX_CRAM_REQUIRED_SAMPLE_ARTIFACTS: Final[set[str]] = {
    "cram",
    "csv",
    "json",
    "extract_stats",
    "snvq",
    "flowq",
    "trimmer_stats",
    "trimmer_failure_codes",
}

TENX_CRAM_CORE_SAMPLE_ARTIFACTS: Final[set[str]] = {"cram", "csv", "json"}


# ---------------------------------------------------------------------------
# 10x raw FASTQ (per-sample FASTQ bundle)
# ---------------------------------------------------------------------------
#
# SOP page 2-3 lists 21 expected files per ``{RunID}-{GroupID}_{Assay}-{UG-BC}``
# prefix.  STAR Solo logs and ``.scRNA.applicationQC.*`` are intentionally
# *optional* (kept in :data:`TENX_FASTQ_OPTIONAL_SUFFIX_TO_ARTIFACT` so they
# are still classifiable but never required).  The R1/R2 ``.fastq.gz``
# entries are tracked indirectly via :func:`validate_10x_raw_fastq_read_mates`
# (which handles ATAC's R3/I2 special case) so they do not appear in this
# table; we focus on the ancillary metadata files that share the same
# Illumina tail and are commonly truncated or omitted.

# Per main prefix (one set per `{RunID}-{GroupID}_{Assay}-{UG-BC}`).
TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    ".csv": "run_csv",
    ".json": "run_json",
    "_trimmer-stats.csv": "trimmer_stats",
    "_trimmer-failure_codes.csv": "trimmer_failure_codes",
    "_trimmer-failure-codes.csv": "trimmer_failure_codes",
    "_unmatched.cram": "unmatched_cram",
    "_unmatched.cram-metadata.json": "unmatched_cram_metadata",
    "_unmatched.csv": "unmatched_csv",
    "_unmatched.json": "unmatched_json",
}

TENX_FASTQ_REQUIRED_PER_PREFIX_ARTIFACTS: Final[set[str]] = {
    "run_csv",
    "run_json",
    "trimmer_stats",
    "trimmer_failure_codes",
    "unmatched_cram",
    "unmatched_csv",
    "unmatched_json",
}

# Per Illumina FASTQ tail (one set per `_S{N}_L{NNN}_R{X}_{NNN}` within a prefix).
TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    ".fastq.gz": "fastq",
    ".fastq.gz-metadata.json": "fastq_metadata",
    ".csv": "fastq_csv",
    ".json": "fastq_json",
    "_sample.fastq.gz": "sample_fastq",
    "_sample.fastq.gz-metadata.json": "sample_fastq_metadata",
}

TENX_FASTQ_REQUIRED_PER_TAIL_ARTIFACTS: Final[set[str]] = {
    "fastq",
    "fastq_csv",
    "fastq_json",
    "sample_fastq",
}

# Optional STAR Solo / application-QC files (recognised but never required).
TENX_FASTQ_OPTIONAL_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    ".scRNA.applicationQC.h5": "scrna_qc_h5",
    ".scRNA.applicationQC.html": "scrna_qc_html",
    "_Log.final.out": "star_log_final",
    "_Log.out": "star_log",
    "_Log.progress.out": "star_log_progress",
    "_ReadsPerGene.out.tab": "star_reads_per_gene",
    "_SJ.out.tab": "star_sj",
}


# ---------------------------------------------------------------------------
# sci raw (per-sample CRAM bundle)
# ---------------------------------------------------------------------------
#
# SOP page 4: ``{RunID}-{GroupID}_{Assay}-{UG-BC}`` family.  The SOP table
# spells the SNVQ metric file ``_trimmer-stats_SNVQ.metrix`` (with an "x"
# typo); we accept the canonical ``.metric`` spelling that producers actually
# use and treat the typo as a synonym.

SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    ".cram": "cram",
    ".cram-metadata.json": "cram_metadata",
    ".csv": "run_csv",
    ".json": "run_json",
    "_trimmer-failure_codes.csv": "trimmer_failure_codes",
    "_trimmer-failure-codes.csv": "trimmer_failure_codes",
    "_trimmer-stats.csv": "trimmer_stats",
    "_trimmer-stats_FlowQ.metric": "flowq",
    "_trimmer-stats_SNVQ.metric": "snvq",
    "_trimmer-stats_SNVQ.metrix": "snvq",
    "_FlowQ.metric": "flowq",
    "_SNVQ.metric": "snvq",
}

SCI_REQUIRED_PER_PREFIX_ARTIFACTS: Final[set[str]] = {
    "cram",
    "run_csv",
    "run_json",
    "trimmer_failure_codes",
    "trimmer_stats",
    "flowq",
    "snvq",
}


# ---------------------------------------------------------------------------
# Scale raw (two-tier: per-well UG_RT and per-UG aggregate)
# ---------------------------------------------------------------------------
#
# SOP page 5 splits Scale artifacts into two granularities:
#
# 1. Per ``{UG_RT}`` (with well code, e.g. ``QSR-8-SCALEPLEX_10A``):
#       .cram, .csv, .json
# 2. Per ``{UG}`` (no well code, e.g. ``QSR-8-SCALEPLEX``):
#       _trimmer-failure_codes.csv, _trimmer-stats.csv,
#       _unmatched.cram, _unmatched.csv, _unmatched.json
#
# The seahub regex captures ``ug_rt`` *without* the well code, so well code
# (separator + 1-2 digits + A-H + ext) is part of the suffix we classify here.

# Suffix shape ``[_-]<wellcode><ext>`` where wellcode = 1-2 digits + A-H.
SCALE_PER_WELL_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(
    r"^[_-](?P<wellcode>\d{1,2}[A-H])(?P<ext>\.(?:cram|csv|json))$"
)

# Optional per-well sidecar (Novogene convention, not SOP-required).
SCALE_PER_WELL_OPTIONAL_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(
    r"^[_-](?P<wellcode>\d{1,2}[A-H])\.cram-metadata\.json$"
)

SCALE_PER_WELL_EXT_TO_ARTIFACT: Final[dict[str, str]] = {
    ".cram": "cram",
    ".csv": "csv",
    ".json": "json",
}

SCALE_REQUIRED_PER_WELL_ARTIFACTS: Final[set[str]] = {"cram", "csv", "json"}

SCALE_PER_UG_SUFFIX_TO_ARTIFACT: Final[dict[str, str]] = {
    "_trimmer-failure_codes.csv": "trimmer_failure_codes",
    "_trimmer-failure-codes.csv": "trimmer_failure_codes",
    "_trimmer-stats.csv": "trimmer_stats",
    "_unmatched.cram": "unmatched_cram",
    "_unmatched.cram-metadata.json": "unmatched_cram_metadata",
    "_unmatched.csv": "unmatched_csv",
    "_unmatched.json": "unmatched_json",
}

SCALE_REQUIRED_PER_UG_ARTIFACTS: Final[set[str]] = {
    "trimmer_failure_codes",
    "trimmer_stats",
    "unmatched_cram",
    "unmatched_csv",
    "unmatched_json",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def classify_suffix(suffix: str, suffix_table: dict[str, str]) -> str | None:
    """Map a literal filename suffix to its canonical artifact key.

    Tries the longest suffixes first so e.g. ``_trimmer-failure_codes.csv``
    is matched before the bare ``.csv`` entry.  Returns ``None`` when no
    table entry matches.
    """
    for known_suffix, artifact in sorted(
        suffix_table.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if suffix == known_suffix:
            return artifact
    return None


def classify_basename(basename: str, suffix_table: dict[str, str]) -> str | None:
    """Find the longest table suffix that matches the *end* of ``basename``.

    Useful for classifying *local* basenames where the leading provider
    naming differs from the canonical S3 prefix.
    """
    for known_suffix, artifact in sorted(
        suffix_table.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if basename.endswith(known_suffix):
            return artifact
    return None


def expected_suffixes_message(suffix_table: dict[str, str]) -> str:
    """Render a stable, longest-first list of accepted suffixes for diagnostics."""
    return ", ".join(sorted(suffix_table, key=len, reverse=True))


def scale_classify_suffix(
    suffix: str,
) -> tuple[str, str, str | None] | None:
    """Classify a Scale suffix as either per-well or per-UG aggregate.

    Returns ``(scope, artifact, well_code)`` where:
    - ``scope`` is ``"per_well"`` or ``"per_ug"``
    - ``artifact`` is the canonical artifact key
    - ``well_code`` is the well position string (e.g. ``"10A"``) for per-well
      artifacts, ``None`` for per-UG aggregates

    Returns ``None`` when the suffix is unrecognised.
    """
    well_match = SCALE_PER_WELL_SUFFIX_RE.match(suffix)
    if well_match:
        ext = well_match.group("ext")
        artifact = SCALE_PER_WELL_EXT_TO_ARTIFACT[ext]
        return ("per_well", artifact, well_match.group("wellcode"))

    aggregate = classify_suffix(suffix, SCALE_PER_UG_SUFFIX_TO_ARTIFACT)
    if aggregate is not None:
        return ("per_ug", aggregate, None)

    return None


def is_scale_per_well_optional_suffix(suffix: str) -> bool:
    """True when ``suffix`` is a recognised optional per-well sidecar.

    These are files that Novogene uploads alongside per-well artifacts but
    that are not listed in the SOP (e.g. ``_10A.cram-metadata.json``).  We
    silently accept them rather than flagging every well's sidecar as an
    ``unexpected_sample_file``.
    """
    return bool(SCALE_PER_WELL_OPTIONAL_SUFFIX_RE.match(suffix))


__all__ = [
    # 10x_cram
    "TENX_CRAM_SAMPLE_SUFFIX_TO_ARTIFACT",
    "TENX_CRAM_REQUIRED_SAMPLE_ARTIFACTS",
    "TENX_CRAM_CORE_SAMPLE_ARTIFACTS",
    # 10x raw FASTQ
    "TENX_FASTQ_PER_PREFIX_SUFFIX_TO_ARTIFACT",
    "TENX_FASTQ_REQUIRED_PER_PREFIX_ARTIFACTS",
    "TENX_FASTQ_PER_TAIL_SUFFIX_TO_ARTIFACT",
    "TENX_FASTQ_REQUIRED_PER_TAIL_ARTIFACTS",
    "TENX_FASTQ_OPTIONAL_SUFFIX_TO_ARTIFACT",
    # sci
    "SCI_PER_PREFIX_SUFFIX_TO_ARTIFACT",
    "SCI_REQUIRED_PER_PREFIX_ARTIFACTS",
    # Scale
    "SCALE_PER_WELL_SUFFIX_RE",
    "SCALE_PER_WELL_EXT_TO_ARTIFACT",
    "SCALE_REQUIRED_PER_WELL_ARTIFACTS",
    "SCALE_PER_UG_SUFFIX_TO_ARTIFACT",
    "SCALE_REQUIRED_PER_UG_ARTIFACTS",
    # Helpers
    "classify_suffix",
    "classify_basename",
    "expected_suffixes_message",
    "scale_classify_suffix",
    "is_scale_per_well_optional_suffix",
]
