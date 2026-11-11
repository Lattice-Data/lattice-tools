"""db2_flattener's SRA CSVs as chunker input, and the name contract.

The flattener writes two files per submission batch:

- ``sra_file_set{N}.csv``: the SRA file sheet, submittable as is. Each row is
  an SRA library; its ``filename``, ``filename2``, ... cells hold the names
  DB2 predicts for SRA (``sra_file_names``): a file's own name, or the names
  of its chunks when its set reaches the SRA limit.
- ``sra_sequence_fileset_set{N}.csv``: ``library_ID,sequence_file_set_uuid``,
  one row per SequenceFileSet whose reads that library submits.

``resolve`` fetches those sets and their read files from DB2 and writes the
same ``sets.json``/``files.json`` that ``plan`` already reads. The CSV only
says which sets; sizes, read counts, S3 URIs and slots come from DB2.

The contract: the names in the SRA sheet are what lands in the destination.
``resolve`` checks that the sheet equals DB2's ``sra_file_names`` for every
library, and ``plan`` checks that its own chunk names equal them too, so the
chunker either produces exactly the predicted names or refuses.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .plan import GB, Plan, PlanError
from .portal import READ_SLOTS

BATCH_LIMIT_BYTES = 5_000 * GB
FILESET_COLUMNS = ["library_ID", "sequence_file_set_uuid"]
LIBRARY_COLUMN = "library_ID"
FILENAME_COLUMN = re.compile(r"filename\d*")
UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")


def read_sra_sheet(path: Path) -> dict[str, list[str]]:
    """library_ID -> its filename cells, in column order, blanks dropped."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        if LIBRARY_COLUMN not in header:
            raise PlanError(f"{path}: no {LIBRARY_COLUMN} column")
        name_cols = [c for c in header if FILENAME_COLUMN.fullmatch(c)]
        if not name_cols:
            raise PlanError(f"{path}: no filename columns")
        rows: dict[str, list[str]] = {}
        problems = []
        for n, row in enumerate(reader, start=2):
            lib = (row.get(LIBRARY_COLUMN) or "").strip()
            if not lib:
                problems.append(f"{path}:{n}: empty {LIBRARY_COLUMN}")
                continue
            if lib in rows:
                problems.append(f"{path}:{n}: {LIBRARY_COLUMN} {lib} repeats")
                continue
            rows[lib] = [
                row[c].strip() for c in name_cols if row.get(c) and row[c].strip()
            ]
    if problems:
        raise PlanError("\n".join(problems))
    return rows


def read_fileset_sheet(path: Path) -> list[tuple[str, str]]:
    """(library_ID, set uuid) per row, in file order."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header != FILESET_COLUMNS:
            raise PlanError(
                f"{path}: header must be {','.join(FILESET_COLUMNS)}, got {header}"
            )
        rows = []
        problems = []
        seen: dict[str, int] = {}
        for n, row in enumerate(reader, start=2):
            if not any(cell.strip() for cell in row):
                continue
            if len(row) != 2:
                problems.append(f"{path}:{n}: expected 2 columns, got {len(row)}")
                continue
            lib, uuid = row[0].strip(), row[1].strip().lower()
            if not UUID.fullmatch(uuid):
                problems.append(f"{path}:{n}: not a uuid: {row[1]!r}")
                continue
            if uuid in seen:
                problems.append(
                    f"{path}:{n}: set {uuid} already listed on line {seen[uuid]}"
                )
                continue
            seen[uuid] = n
            rows.append((lib, uuid))
    if problems:
        raise PlanError("\n".join(problems))
    if not rows:
        raise PlanError(f"{path}: no rows")
    return rows


def uuid_of(at_id: str) -> str:
    m = UUID.search(at_id)
    if not m:
        raise PlanError(f"no uuid in {at_id!r}")
    return m.group(0)


def slot_refs(fset: dict) -> list[tuple[str, str]]:
    """(slot, file uuid) for every read slot the set fills."""
    out = []
    for slot in READ_SLOTS:
        ref = fset.get(slot)
        if not ref:
            continue
        if isinstance(ref, dict):
            ref = ref.get("@id") or ref.get("uuid") or ""
        out.append((slot, uuid_of(str(ref))))
    return out


@dataclass
class Resolved:
    sets: list[dict]
    files: list[dict]
    total_bytes: int
    n_libraries: int


def resolve(
    sheet: dict[str, list[str]],
    fileset_rows: list[tuple[str, str]],
    sets_by_uuid: dict[str, dict],
    files_by_uuid: dict[str, dict],
    batch_limit_bytes: int = BATCH_LIMIT_BYTES,
) -> Resolved:
    """Check the two sheets against DB2 and return the objects ``plan`` reads.

    Pure: the caller fetches. Every problem is collected into one PlanError.
    """
    problems: list[str] = []
    libs_with_sets = {lib for lib, _ in fileset_rows}
    for lib, _ in fileset_rows:
        if lib not in sheet:
            problems.append(f"{lib}: listed in the set sheet but not in the SRA sheet")
    for lib in sheet:
        if lib not in libs_with_sets:
            problems.append(f"{lib}: in the SRA sheet but has no set in the set sheet")

    sets: list[dict] = []
    files: dict[str, dict] = {}
    file_owner: dict[str, str] = {}
    predicted: dict[str, list[str]] = {lib: [] for lib in libs_with_sets}

    for lib, set_uuid in fileset_rows:
        fset = sets_by_uuid.get(set_uuid)
        if fset is None:
            problems.append(
                f"{lib}: set {set_uuid} not found in DB2 (missing or deleted)"
            )
            continue
        if fset.get("status") == "deleted":
            problems.append(f"{lib}: set {set_uuid} is deleted")
            continue
        refs = slot_refs(fset)
        if not refs:
            problems.append(f"{lib}: set {set_uuid} has no read slots")
            continue
        sets.append(fset)
        for slot, file_uuid in refs:
            where = f"{lib}: set {set_uuid} {slot} file {file_uuid}"
            if file_uuid in file_owner:
                problems.append(f"{where} is also in set {file_owner[file_uuid]}")
                continue
            file_owner[file_uuid] = set_uuid
            fobj = files_by_uuid.get(file_uuid)
            if fobj is None:
                problems.append(f"{where} not found in DB2 (missing or deleted)")
                continue
            if fobj.get("status") == "deleted":
                problems.append(f"{where} is deleted")
                continue
            names = fobj.get("sra_file_names")
            if not names:
                problems.append(
                    f"{where} has no sra_file_names (DB2 gives none when the set "
                    "would be refused, e.g. a file in two sets)"
                )
                continue
            files[file_uuid] = fobj
            predicted[lib].extend(names)

    for lib in libs_with_sets & set(sheet):
        want, got = Counter(predicted[lib]), Counter(sheet[lib])
        if want != got:
            missing = sorted((want - got).elements())
            extra = sorted((got - want).elements())
            detail = []
            if missing:
                detail.append(f"predicted but not in the sheet: {', '.join(missing)}")
            if extra:
                detail.append(f"in the sheet but not predicted: {', '.join(extra)}")
            problems.append(f"{lib}: names differ; " + "; ".join(detail))

    total = sum(int(f.get("file_size") or 0) for f in files.values())
    if total > batch_limit_bytes:
        problems.append(
            f"batch is {total / GB:.0f} GB, above the {batch_limit_bytes / GB:.0f} GB limit"
        )

    if problems:
        raise PlanError("\n".join(problems))
    return Resolved(
        sets=sets,
        files=list(files.values()),
        total_bytes=total,
        n_libraries=len(libs_with_sets),
    )


def name_contract(plan: Plan, files: list[dict]) -> None:
    """Refuse a plan whose names differ from DB2's ``sra_file_names``.

    Applies only when the file objects carry ``sra_file_names`` (DB2 exports
    and ``resolve`` output); then every file must carry it.
    """
    if not any("sra_file_names" in f for f in files):
        return
    by_uri = {f.get("s3_uri"): f for f in files}
    problems = []
    for g in plan.groups:
        for fp in g.files:
            actual = (
                [c.name for c in fp.chunks]
                if g.action == "split"
                else [fp.file.filename]
            )
            fobj = by_uri.get(fp.file.s3_uri, {})
            expected = fobj.get("sra_file_names")
            if not expected:
                problems.append(f"{fp.file.filename}: no sra_file_names in DB2")
            elif list(expected) != actual:
                problems.append(
                    f"{fp.file.filename}: plan names {actual} differ from DB2 "
                    f"sra_file_names {list(expected)}"
                )
    if problems:
        raise PlanError(
            "plan does not match DB2's predicted SRA names:\n" + "\n".join(problems)
        )
