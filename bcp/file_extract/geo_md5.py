"""GEO md5 CSV for processed files: one row per S3 object, named as GEO sees it.

Why a table in and a CSV out
----------------------------
GEO wants an MD5 beside every processed file it receives, and nothing upstream
holds one. DB2 files carry ``crc64nvme_base64`` only, and the other file_extract
modes fetch that same CRC from S3 object attributes without reading a byte of
content. So the md5 has to be computed from the bytes, which is what this module
does: a full streamed read of each object, the one expensive operation in the
package, kept behind its own subcommand so no header-only run pays for it.

The "file name" column is the harder half. Every Cell Ranger matrix is
``sample_filtered_feature_bc_matrix.h5``, so the S3 basename cannot name a GEO
file. The name GEO sees is the RawMatrixFile alias: db2-flattener writes it into
the GEO.csv "processed data file" column with the ``lab:`` namespace stripped,
and the upload step renames each download to it. That mapping lives in the
RawMatrixFile sheet, not in S3, so this module takes a table, any TSV or CSV with
an ``s3_uri`` column, plus optionally the column holding the name. Given
``--name-column aliases`` on the submission sheet, the two GEO files agree by
construction. Without it the basename is used, which is right for Scale h5ads
and wrong for Cell Ranger h5s, and the duplicate-name guard says so before a
byte is downloaded.

What stops the run, and why
---------------------------
Everything that would otherwise hand GEO a confident wrong answer fails before
any download: a missing or empty URI (a row is never skipped, because a GEO file
that cannot be fetched is exactly what this CSV exists to prevent), a name that
repeats (GEO rejects duplicate file names), two rows naming one object, a name
with a path separator (GEO takes bare file names), and an alias cell listing
more than one alias (one file, one name). A name whose extension differs from
the object's is only a warning: an alias without ``.h5`` is a legitimate
submitter choice, but it deserves a second look before upload.

A failed read leaves the checksum cell empty and is listed, the same contract as
every other mode; ``--strict`` turns that into exit 1.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from tqdm import tqdm

from .constants import GEO_MD5_COLUMNS, GEO_MD5_DEFAULT_WORKERS, GEO_MD5_URI_COLUMN
from .models import RunSummary
from .retry import retry_with_backoff
from .s3_utils import parse_s3_object_uri, s3_uri_for, stream_md5

# Suffixes that ride on another extension: ``.mtx.gz`` is one extension, not two.
_COMPRESSION_SUFFIXES = frozenset({".gz", ".bz2", ".xz", ".zst", ".zip"})

# The characters a file name can carry through FTP, a spreadsheet cell and a
# shell loop without quoting trouble. Anything else is worth a look, not a stop.
_SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")


class GeoMd5Error(Exception):
    """A table problem that stops the run before any object is read."""


@dataclass(frozen=True)
class GeoFile:
    """One processed file: where it is, and what GEO will call it."""

    row: int
    """1-based data row in the input table, for messages."""

    bucket: str
    key: str
    name: str

    @property
    def s3_uri(self) -> str:
        return s3_uri_for(self.bucket, self.key)

    @property
    def basename(self) -> str:
        return self.key.rsplit("/", 1)[-1]


def read_table(path: str | Path) -> list[dict[str, str]]:
    """Rows of a TSV or CSV as dicts; the delimiter follows the file suffix.

    ``utf-8-sig`` so a sheet exported from Excel, which prepends a byte-order
    mark, still has ``s3_uri`` as its first column name rather than
    ``\\ufeffs3_uri``.
    """
    table = Path(path)
    delimiter = "," if table.suffix.lower() == ".csv" else "\t"
    try:
        with table.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter=delimiter)
            if reader.fieldnames is None:
                raise GeoMd5Error(f"{table} is empty; expected a header row")
            named = [name for name in reader.fieldnames if name]
            repeated = sorted(n for n, count in Counter(named).items() if count > 1)
            if repeated:
                raise GeoMd5Error(
                    f"{table} repeats column names: "
                    + ", ".join(repr(name) for name in repeated)
                )
            # A row of nothing but delimiters, which a spreadsheet export often
            # leaves at the end, is not a data row. A row with any content is
            # kept, and plan_geo_files says what it lacks.
            return [dict(row) for row in reader if _has_content(row)]
    except UnicodeDecodeError as exc:
        raise GeoMd5Error(
            f"{table} is not UTF-8 text ({exc.reason}); export the sheet as a "
            "UTF-8 CSV or TSV"
        ) from exc


def _has_content(row: dict[str | None, Any]) -> bool:
    for value in row.values():
        # csv.DictReader parks cells beyond the header under None, as a list.
        values = value if isinstance(value, list) else [value]
        if any((item or "").strip() for item in values):
            return True
    return False


def geo_file_name(cell: str) -> str:
    """GEO file name from a name cell: ``name``, ``lab:name`` or ``["lab:name"]``.

    Mirrors db2-flattener's alias cleaning, everything after the first ``:``,
    so this CSV names a file exactly as the GEO.csv does.
    """
    value = (cell or "").strip()
    if value.startswith("["):
        try:
            items = json.loads(value)
        except json.JSONDecodeError as exc:
            raise GeoMd5Error(f"name cell {value!r} is not a JSON alias array") from exc
        if not isinstance(items, list) or not all(
            isinstance(item, str) for item in items
        ):
            raise GeoMd5Error(f"name cell {value!r} is not a list of alias strings")
        if len(items) != 1:
            raise GeoMd5Error(
                f"name cell {value!r} lists {len(items)} aliases; "
                "one GEO file name is needed"
            )
        value = items[0].strip()
    if ":" in value:
        value = value.split(":", 1)[1].strip()
    if not value:
        raise GeoMd5Error("name cell is empty")
    if set(value) == {"."}:
        raise GeoMd5Error(f"name cell {cell.strip()!r} is not a file name")
    return value


def duplicate_names_message(names: Sequence[str], name_column: str | None) -> str:
    listed = ", ".join(repr(name) for name in names)
    message = f"GEO file names repeat: {listed}. GEO rejects duplicate file names."
    if name_column is None:
        message += (
            " The S3 basename is not unique here; pass --name-column with the "
            "column holding each file's alias (e.g. the RawMatrixFile sheet's "
            "aliases)."
        )
    return message


def plan_geo_files(
    rows: Sequence[dict[str, str]],
    *,
    uri_column: str = GEO_MD5_URI_COLUMN,
    name_column: str | None = None,
) -> list[GeoFile]:
    """Resolve every row to a GeoFile, or raise before anything is downloaded."""
    if not rows:
        raise GeoMd5Error("table has no data rows")
    columns = [column for column in rows[0] if column]
    for column in (uri_column, name_column):
        if column is not None and column not in columns:
            raise GeoMd5Error(
                f"column {column!r} not in table; columns are: "
                + ", ".join(repr(name) for name in columns)
            )

    files: list[GeoFile] = []
    for row_number, row in enumerate(rows, start=1):
        uri = (row.get(uri_column) or "").strip()
        if not uri:
            raise GeoMd5Error(f"row {row_number}: empty {uri_column!r}")
        try:
            bucket, key = parse_s3_object_uri(uri)
        except ValueError as exc:
            raise GeoMd5Error(f"row {row_number}: {exc}") from exc
        if name_column is None:
            name = key.rsplit("/", 1)[-1]
        else:
            try:
                name = geo_file_name(row.get(name_column) or "")
            except GeoMd5Error as exc:
                raise GeoMd5Error(f"row {row_number}: {exc}") from exc
        if "/" in name or "\\" in name:
            raise GeoMd5Error(
                f"row {row_number}: GEO file name {name!r} contains a path separator"
            )
        files.append(GeoFile(row=row_number, bucket=bucket, key=key, name=name))

    uri_counts = Counter(file.s3_uri for file in files)
    repeated_uris = sorted(uri for uri, count in uri_counts.items() if count > 1)
    if repeated_uris:
        raise GeoMd5Error(
            "the same object is listed more than once: " + ", ".join(repeated_uris)
        )
    name_counts = Counter(file.name for file in files)
    repeated = sorted(name for name, count in name_counts.items() if count > 1)
    if repeated:
        raise GeoMd5Error(duplicate_names_message(repeated, name_column))
    return files


def _extension(filename: str) -> str:
    suffixes = Path(filename).suffixes
    if not suffixes:
        return ""
    if suffixes[-1].lower() in _COMPRESSION_SUFFIXES and len(suffixes) > 1:
        return "".join(suffixes[-2:])
    return suffixes[-1]


def extension_warning(file: GeoFile) -> str | None:
    """A warning when the GEO name's extension differs from the S3 object's."""
    name_ext = _extension(file.name)
    key_ext = _extension(file.basename)
    if name_ext == key_ext:
        return None
    return (
        f"row {file.row}: GEO file name {file.name!r} ends in "
        f"{name_ext or 'no extension'} but the S3 object ends in "
        f"{key_ext or 'no extension'}"
    )


def character_warning(file: GeoFile) -> str | None:
    """A warning when the GEO name holds whitespace or punctuation beyond ._-."""
    if _SAFE_NAME_RE.match(file.name):
        return None
    return (
        f"row {file.row}: GEO file name {file.name!r} has characters outside "
        "letters, digits, '.', '_' and '-'; the upload loop and the GEO "
        "template handle those names unevenly"
    )


def case_warnings(files: Sequence[GeoFile]) -> list[str]:
    """Warnings for names that differ only in case.

    GEO's FTP host keeps such names apart, so this is not a duplicate for
    GEO. But the renaming step runs on the submitter's machine, and the
    macOS and Windows default filesystems fold case: the second file would
    overwrite the first before the upload, leaving two names for one file.
    """
    groups: dict[str, list[GeoFile]] = {}
    for file in files:
        groups.setdefault(file.name.casefold(), []).append(file)
    warnings: list[str] = []
    for group in groups.values():
        if len(group) > 1:
            rows = ", ".join(str(file.row) for file in group)
            names = ", ".join(repr(file.name) for file in group)
            warnings.append(
                f"rows {rows}: GEO file names {names} differ only in case; a "
                "case-insensitive filesystem merges them during the renaming step"
            )
    return warnings


def name_warnings(files: Sequence[GeoFile]) -> list[str]:
    """Every advisory warning about the names, in table order."""
    warnings: list[str] = []
    for file in files:
        for warning in (extension_warning(file), character_warning(file)):
            if warning:
                warnings.append(warning)
    warnings.extend(case_warnings(files))
    return warnings


def default_geo_md5_output_name(table_path: str | Path) -> str:
    return f"{Path(table_path).stem}_geo_md5.csv"


def write_geo_md5_csv(path: str | Path, rows: Sequence[tuple[str, str]]) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(GEO_MD5_COLUMNS)
        writer.writerows(rows)


def extract_geo_md5(
    s3_client: Any,
    files: Sequence[GeoFile],
    output_path: str | Path,
    *,
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
) -> RunSummary:
    """Hash every file and write the GEO CSV, rows in table order."""
    summary = RunSummary(total=len(files))
    results: dict[int, tuple[str | None, int, str]] = {}
    max_workers = min(workers or GEO_MD5_DEFAULT_WORKERS, max(len(files), 1))

    # Header first, so an output path that cannot be written fails before the
    # first byte is streamed rather than after the last. The other modes open
    # their TSV before their pool for the same reason.
    write_geo_md5_csv(output_path, [])

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                retry_with_backoff,
                stream_md5,
                s3_client,
                file.bucket,
                file.key,
                retries=retries,
            ): file
            for file in files
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(files), desc="Hashing")
        try:
            for fut in iterator:
                file = futures[fut]
                value, err = fut.result()
                if err or value is None:
                    results[file.row] = (None, 0, err or "no result")
                else:
                    md5, nbytes = value
                    results[file.row] = (md5, nbytes, "")
        except BaseException:
            # Ctrl-C must not drain the queue. The executor's __exit__ waits
            # for every queued download to run to completion and the hashes
            # are then discarded; cancelling first means only the streams
            # already in flight finish, since a thread cannot be interrupted.
            executor.shutdown(wait=False, cancel_futures=True)
            raise

    rows: list[tuple[str, str]] = []
    for file in files:
        md5, nbytes, err = results[file.row]
        rows.append((file.name, md5 or ""))
        if err:
            summary.failures.append((file.s3_uri, err, ""))
        else:
            summary.md5_ok += 1
            summary.bytes_total += nbytes
    write_geo_md5_csv(output_path, rows)
    return summary
