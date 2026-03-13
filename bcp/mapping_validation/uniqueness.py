from __future__ import annotations

from collections import defaultdict
from typing import Iterable, List

from .parsing import MappingRow


def validate_uniqueness(mappings: Iterable[MappingRow]) -> dict:
    """
    Check that each local path maps to exactly one S3 path and vice versa.

    Returns a dictionary with:
    - duplicate_locals: {local_path: [MappingRow, ...]}
    - duplicate_s3: {s3_path: [MappingRow, ...]}
    - total: total number of mappings inspected
    """
    local_to_rows: dict[str, List[MappingRow]] = defaultdict(list)
    s3_to_rows: dict[str, List[MappingRow]] = defaultdict(list)

    total = 0
    for row in mappings:
        total += 1
        local_to_rows[row.local_path].append(row)
        s3_to_rows[row.s3_path].append(row)

    duplicate_locals = {k: v for k, v in local_to_rows.items() if len(v) > 1}
    duplicate_s3 = {k: v for k, v in s3_to_rows.items() if len(v) > 1}

    return {
        "duplicate_locals": duplicate_locals,
        "duplicate_s3": duplicate_s3,
        "total": total,
    }


__all__ = ["validate_uniqueness"]
