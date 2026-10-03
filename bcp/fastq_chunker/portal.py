"""Build chunking groups from Lattice data portal JSON.

A ``sequence_file_set`` names its reads in the ``read1``/``read2``/``read3``/
``index1``/``index2`` slots, but embeds only their ``@id``; sizes, read counts and
S3 URIs live on the ``sequence_file`` objects, so both must be supplied.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from .plan import FastqFile, PlanError, stem_of

READ_SLOTS = ("read1", "read2", "read3", "index1", "index2")
FILE_COLLECTION = "/sequence_files/"


def load_objects(path: Path) -> list[dict]:
    """Accept one object, a JSON list, or a search export carrying ``@graph``."""
    data = json.loads(Path(path).read_text())
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and "@graph" in data:
        return list(data["@graph"])
    if isinstance(data, dict):
        return [data]
    raise PlanError(f"{path}: expected a JSON object or list")


def object_id(obj: dict, collection: str) -> str:
    if "@id" in obj:
        return obj["@id"]
    if "uuid" in obj:
        return f"{collection}{obj['uuid']}/"
    raise PlanError(f"object has neither @id nor uuid: {json.dumps(obj)[:120]}")


def ref_id(value: object) -> str:
    """A slot value is an ``@id`` string, a bare uuid, or an embedded object."""
    if isinstance(value, dict):
        return object_id(value, FILE_COLLECTION)
    if isinstance(value, str):
        if value.startswith("/"):
            return value
        return f"{FILE_COLLECTION}{value}/"
    raise PlanError(f"unrecognised file reference: {value!r}")


def set_label(fset: dict) -> str:
    return fset.get("summary") or (fset.get("aliases") or [object_id(fset, "")])[0]


def build_groups(sets: list[dict], files: list[dict]) -> list[list[FastqFile]]:
    """One list of :class:`FastqFile` per set that has at least one read slot.

    Every problem found is reported in a single :class:`PlanError` rather than
    one at a time, since a portal export tends to have several at once.
    """
    by_id = {object_id(f, FILE_COLLECTION): f for f in files}
    problems: list[str] = []
    owners: dict[str, list[str]] = defaultdict(list)
    groups: list[list[FastqFile]] = []

    for fset in sets:
        set_id = object_id(fset, "/sequence_file_sets/")
        members: list[FastqFile] = []
        for slot in READ_SLOTS:
            if slot not in fset:
                continue
            fid = ref_id(fset[slot])
            owners[fid].append(set_id)
            fobj = by_id.get(fid)
            if fobj is None:
                problems.append(f"{set_id} {slot}: file {fid} not supplied")
                continue
            member = _fastq_file(set_id, slot, fid, fobj, problems)
            if member is not None:
                members.append(member)
        if members:
            groups.append(members)

    for fid, set_ids in owners.items():
        if len(set_ids) > 1:
            problems.append(
                f"file {fid} belongs to more than one set: {', '.join(set_ids)}"
            )

    if problems:
        raise PlanError("\n".join(problems))
    return groups


def _fastq_file(
    set_id: str, slot: str, fid: str, fobj: dict, problems: list[str]
) -> FastqFile | None:
    where = f"{set_id} {slot} ({fid})"
    if fobj.get("file_format") != "fastq":
        problems.append(f"{where}: file_format is {fobj.get('file_format')!r}")
        return None
    missing = [k for k in ("s3_uri", "file_size", "read_count") if not fobj.get(k)]
    if missing:
        problems.append(f"{where}: missing {', '.join(missing)}")
        return None
    filename = fobj["s3_uri"].rsplit("/", 1)[-1]
    try:
        stem_of(filename)
    except PlanError as e:
        problems.append(f"{where}: {e}")
        return None
    return FastqFile(
        group=set_id,
        slot=slot,
        filename=filename,
        s3_uri=fobj["s3_uri"],
        size_bytes=int(fobj["file_size"]),
        reads=int(fobj["read_count"]),
    )


def labels(sets: list[dict]) -> dict[str, str]:
    return {object_id(s, "/sequence_file_sets/"): set_label(s) for s in sets}
