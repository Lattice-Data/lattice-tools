"""Protospacer set signature (guidesig).

Computes a deterministic string identity over the set of protospacer sequences
in a guide-template TSV. This is a string-identity check, not a biological
equivalence check: length, strand, and orientation are deliberately not
normalized.
"""

from __future__ import annotations

from .core import GuideSigError, protospacer_set, signature

__all__ = [
    "GuideSigError",
    "protospacer_set",
    "signature",
]
