"""Protospacer set signature (guidesig).

Computes a deterministic string identity over the set of protospacer sequences
in a guide-template TSV or CSV. This is a string-identity check, not a
biological equivalence check: length, strand, and orientation are deliberately
not normalized. The input format must always be stated explicitly; it is never
inferred from the file name or contents.
"""

from __future__ import annotations

from .core import (
    DEFAULT_COLUMN,
    FILE_FORMATS,
    ID_COLUMN,
    GuideSigError,
    protospacer_set,
    protospacers_and_guide_ids,
    signature,
    signature_of,
)

__all__ = [
    "DEFAULT_COLUMN",
    "FILE_FORMATS",
    "ID_COLUMN",
    "GuideSigError",
    "protospacer_set",
    "protospacers_and_guide_ids",
    "signature",
    "signature_of",
]
