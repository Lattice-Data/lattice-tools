from __future__ import annotations

# `classify_cas` is deliberately NOT re-exported. It lives in `cas_registry`, at the
# bcp/ top level, because CAS validation belongs to neither this package nor
# chebi_lookup -- and re-exporting it here would put `structure_check.classify_cas`
# back in the public surface of the package it was moved out of, recoupling the name
# the move decoupled. Import it from `cas_registry`.
from .cli import main
from .client import check_row, compare_structures, comparison_verdict_from_inchi
from .inchi import defined_side, defined_stereo, is_multi_component

__all__ = [
    "main",
    "check_row",
    "compare_structures",
    "comparison_verdict_from_inchi",
    "defined_side",
    "defined_stereo",
    "is_multi_component",
]
