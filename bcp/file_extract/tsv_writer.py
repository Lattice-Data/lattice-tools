from __future__ import annotations

import csv
import threading
from pathlib import Path
from typing import Sequence


class TsvWriter:
    """Thread-safe incremental TSV writer."""

    def __init__(self, path: str | Path, columns: Sequence[str]) -> None:
        self._path = Path(path)
        self._columns = list(columns)
        self._lock = threading.Lock()
        self._write_header()

    def _write_header(self) -> None:
        with self._path.open("w", newline="", encoding="utf-8") as fh:
            csv.writer(fh, delimiter="\t").writerow(self._columns)

    def append_row(self, row: Sequence[object]) -> None:
        with self._lock:
            with self._path.open("a", newline="", encoding="utf-8") as fh:
                csv.writer(fh, delimiter="\t").writerow(row)

    @property
    def path(self) -> Path:
        return self._path

    @property
    def columns(self) -> list[str]:
        return list(self._columns)
