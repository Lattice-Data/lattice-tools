"""Tests for file_extract.tsv_writer."""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from file_extract.tsv_writer import TsvWriter


def test_writes_header_and_rows(tmp_path: Path) -> None:
    out = tmp_path / "out.tsv"
    writer = TsvWriter(out, ["a", "b"])
    writer.append_row([1, 2])
    writer.append_row([3, 4])

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows == [["a", "b"], ["1", "2"], ["3", "4"]]


def test_concurrent_append(tmp_path: Path) -> None:
    out = tmp_path / "concurrent.tsv"
    writer = TsvWriter(out, ["col"])
    n = 50

    def write_row(i: int) -> None:
        writer.append_row([i])

    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(write_row, range(n)))

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert len(rows) == n + 1
    assert rows[0] == ["col"]
