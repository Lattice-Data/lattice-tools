"""96-well pairing for Scale samples.csv barcodes and sheet RT_index.

No I/O. Shared with future scale_cram.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Sequence

PLATE_ROWS = "ABCDEFGH"
PLATE_COLS = range(1, 13)
SCALEQUANT_PREFIX = "SCALEQUANT-"

# Column-wise: 1A, 1B, ... 1H, 2A, ... 12H.
COLUMN_WISE_WELLS = tuple(f"{col}{row}" for col in PLATE_COLS for row in PLATE_ROWS)
_WELL_INDEX = {well: i for i, well in enumerate(COLUMN_WISE_WELLS)}

_ROW_COL_RE = re.compile(r"^([A-H])(\d{1,2})$")
_COL_ROW_RE = re.compile(r"^(\d{1,2})([A-H])$")


class ScaleExtractError(Exception):
    """User-facing Scale extract failure."""


class WellOwnershipError(ScaleExtractError):
    """A well is claimed by more than one samples.csv row."""


@dataclass(frozen=True)
class ControlSample:
    sample: str
    barcodes: str


@dataclass(frozen=True)
class CorrelationResult:
    paired: tuple[str, ...]
    controls: tuple[ControlSample, ...]
    sample_names: dict[str, tuple[str, ...]]

    @property
    def paired_set(self) -> frozenset[str]:
        return frozenset(self.paired)

    @property
    def control_set(self) -> frozenset[str]:
        return frozenset(item.sample for item in self.controls)


def _require_col(col: int, *, raw: str) -> int:
    if col not in PLATE_COLS:
        raise ScaleExtractError(f"Invalid well {raw!r}: column must be 1-12")
    return col


def normalize_rt_index(value: str) -> str:
    """Strip SCALEQUANT- and flip A11 -> 11A. Already-flipped 11A is kept."""
    text = (value or "").strip()
    if text.upper().startswith(SCALEQUANT_PREFIX):
        text = text[len(SCALEQUANT_PREFIX) :]
    row_col = _ROW_COL_RE.fullmatch(text)
    if row_col:
        row = row_col.group(1)
        col = _require_col(int(row_col.group(2)), raw=value)
        return f"{col}{row}"
    col_row = _COL_ROW_RE.fullmatch(text)
    if col_row:
        col = _require_col(int(col_row.group(1)), raw=value)
        return f"{col}{col_row.group(2)}"
    raise ScaleExtractError(f"Invalid RT_index {value!r}")


def parse_rt_index_cell(value: str) -> tuple[str, ...]:
    """Normalize a comma-separated RT_index cell to column-row wells."""
    text = (value or "").strip()
    if not text:
        return ()
    wells: list[str] = []
    seen: set[str] = set()
    for token in text.split(","):
        part = token.strip()
        if not part:
            continue
        well = normalize_rt_index(part)
        if well not in seen:
            wells.append(well)
            seen.add(well)
    return tuple(wells)


def parse_well(token: str) -> str:
    """Parse a samples.csv well in column-row form (1A, 12H)."""
    text = (token or "").strip()
    match = _COL_ROW_RE.fullmatch(text)
    if not match:
        raise ScaleExtractError(
            f"Invalid barcodes well {token!r}: expected column-row form such as 1A"
        )
    col = _require_col(int(match.group(1)), raw=token)
    return f"{col}{match.group(2)}"


def expand_barcodes(barcodes: str) -> tuple[str, ...]:
    """Expand a barcodes cell into column-wise 96-well ids.

    Segments are split on ``;``. Each segment is a single well (``1A``,
    ``12H``) or a start-end range (``1A-2C``) walked in column-wise order.
    """
    text = (barcodes or "").strip()
    if not text:
        return ()
    wells: list[str] = []
    seen: set[str] = set()
    for segment in text.split(";"):
        part = segment.strip()
        if not part:
            continue
        if "-" in part:
            start_raw, end_raw = part.split("-", 1)
            start = parse_well(start_raw)
            end = parse_well(end_raw)
            start_i = _WELL_INDEX[start]
            end_i = _WELL_INDEX[end]
            if start_i > end_i:
                raise ScaleExtractError(
                    f"Invalid barcodes range {part!r}: start is after end"
                )
            for well in COLUMN_WISE_WELLS[start_i : end_i + 1]:
                if well not in seen:
                    wells.append(well)
                    seen.add(well)
        else:
            well = parse_well(part)
            if well not in seen:
                wells.append(well)
                seen.add(well)
    return tuple(wells)


def correlate_samples(
    sample_rows: Sequence[tuple[str, str]],
    sheet_wells: Iterable[str],
    sheet_names: Iterable[tuple[str, str]] = (),
) -> CorrelationResult:
    """Pair samples.csv rows to sheet wells; unpaired rows are controls.

    ``sample_rows`` is ``(sample, barcodes)``. ``sheet_wells`` are already
    normalized to column-row form (``11A``). ``sheet_names`` is
    ``(well, sample_name)`` from the sample template. A well claimed by two
    rows raises ``WellOwnershipError``.
    """
    sheet = set(sheet_wells)
    names_by_well: dict[str, list[str]] = {}
    for well, name in sheet_names:
        if not name:
            continue
        bucket = names_by_well.setdefault(well, [])
        if name not in bucket:
            bucket.append(name)

    owner: dict[str, str] = {}
    expanded: list[tuple[str, str, tuple[str, ...]]] = []
    for sample, barcodes in sample_rows:
        wells = expand_barcodes(barcodes)
        for well in wells:
            prior = owner.get(well)
            if prior is not None and prior != sample:
                raise WellOwnershipError(
                    f"Well {well} is claimed by both {prior!r} and {sample!r}"
                )
            owner[well] = sample
        expanded.append((sample, barcodes, wells))

    paired: list[str] = []
    controls: list[ControlSample] = []
    sample_names: dict[str, tuple[str, ...]] = {}
    for sample, barcodes, wells in expanded:
        seen: list[str] = []
        for well in wells:
            for name in names_by_well.get(well, []):
                if name not in seen:
                    seen.append(name)
        sample_names[sample] = tuple(seen)
        if any(well in sheet for well in wells):
            paired.append(sample)
        else:
            controls.append(ControlSample(sample=sample, barcodes=barcodes))
    return CorrelationResult(
        paired=tuple(paired),
        controls=tuple(controls),
        sample_names=sample_names,
    )
