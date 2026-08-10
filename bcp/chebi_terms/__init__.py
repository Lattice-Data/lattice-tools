from __future__ import annotations

from .cli import main
from .client import describe_chebi_id, verify_chebi_id

__all__ = ["main", "describe_chebi_id", "verify_chebi_id"]
