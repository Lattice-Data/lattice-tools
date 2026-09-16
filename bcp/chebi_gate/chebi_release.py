"""The ChEBI flat-file release: discovery, distillation to a queryable index, pinning.

"Is this substance already in ChEBI?" cannot be answered by the ChEBI REST API,
which resolves an identifier to a compound but offers no reverse lookup from a
structure. Only a bulk release does, so the gate reads a monthly release --
about 870 MB unpacked -- and distils it once into a small index it can query per
run. The release and the index are caches: neither is committed, both are hashed
into the run manifest so a verdict from three months ago can be explained.

Four things here were measured against release 255 rather than taken on trust,
because each contradicts something the port was handed.

**Paths are discovered, never hard-coded.** All four paths in the original plan
404 today; ``Flat_file_tab_delimited/`` is gone and ``chebiId_inchi.tsv`` no longer
exists in any form. Discovery reads the directory listing and matches *exact*
filenames. The prototype matched loose patterns and broke ties on filename
length, so its pick of ``structures.tsv`` over ``structure_registry.tsv`` was an
accident of 17 characters being fewer than 25.

**A missing table fails loudly.** The prototype printed "NOT FOUND
automatically" and exited 0, so a run that retrieved 0 of 6 tables looked like a
success to CI.

**There is no obsolete or merged structure match.** The handoff says liveness
comes from ``parent_id`` in ``compounds.tsv``. Measured: ``parent_id`` is
non-empty in 0 of 218,533 csv-parsed rows -- a naive line-split says 3, but those
are phantom rows produced by splitting the ``definition`` column on its embedded
newlines. ``status.tsv`` is 35 bytes holding CHECKED, OK and SUBMITTED, with no
deleted value, so the handoff's first half is right and its second half is not.
The real merge pointer is ``secondary_ids.tsv``, which the prototype never
downloaded, and the question resolves the other way: all 189,896 structure-bearing
compound ids are primary ids and none is a secondary, so a structure match can
never land on a merged entry. The prototype's "merged entry is not a duplicate"
branch was unreachable by construction, not merely dead. ``secondary_ids`` is
still loaded, for the different job of resolving an identifier somebody hands us.

**Names come from compounds.tsv, and ascii_name is preferred.** ``names.tsv``
holds synonyms, not the entry name. In ``compounds.tsv``, ``name`` carries markup
in 28,243 rows -- ``(<i>R</i>)-linalool`` -- against 1,349 in ``ascii_name``,
which is never blank.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

FLAT_FILES_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/"

# Exact filenames, with or without gzip. Loose patterns are how the prototype came
# to prefer one table over another by name length.
WANTED_TABLES = {
    "structures": "structures.tsv",
    "compounds": "compounds.tsv",
    "database_accession": "database_accession.tsv",
    "secondary_ids": "secondary_ids.tsv",
    "status": "status.tsv",
}

# Tables the gate can work without, reported rather than required.
OPTIONAL_TABLES = frozenset({"secondary_ids", "status"})

INDEX_FILES = {
    "inchikey": "inchikey_index.tsv",
    "compound": "compound_index.tsv",
    "cas": "cas_index.tsv",
    "secondary": "secondary_index.tsv",
}
INDEX_MANIFEST = "index_manifest.json"

INCHIKEY_COLUMNS = ("chebi_id", "inchikey")
# Names live in their own index rather than beside the structures, because an
# entry with no structure still has to be nameable: an identifier handed to us for
# resolution is often one we have no structure for, and reporting a bare
# "CHEBI:90" where "(-)-epicatechin" was available is a worse answer.
COMPOUND_COLUMNS = ("chebi_id", "stars", "ascii_name")
CAS_COLUMNS = ("cas", "chebi_id")
SECONDARY_COLUMNS = ("secondary_id", "primary_id")

# The accession type code for a CAS Registry Number. Not the human label: the
# column holds codes, and 40,203 of 422,942 accession rows carry this one.
CAS_ACCESSION_TYPE = "CAS"

_HREF = re.compile(r'href="([^"?/][^"]*)"', re.I)


class ReleaseError(Exception):
    """Raised when a release is incomplete or an index cannot be built."""


def _raise_field_size_limit() -> None:
    """structures.tsv embeds multi-line molfiles in a quoted column and is 755 MB.

    It must be parsed as real delimited text with the field-size limit raised, not
    by splitting lines on tabs. Splitting lines is also what produces the phantom
    rows in compounds.tsv, whose definition column contains newlines.
    """
    csv.field_size_limit(sys.maxsize)


@dataclass(frozen=True)
class Index:
    """A distilled release: enough to answer the gate's questions, small enough to hold."""

    by_inchikey: dict[str, tuple[str, ...]] = field(default_factory=dict)
    by_skeleton: dict[str, tuple[str, ...]] = field(default_factory=dict)
    by_cas: dict[str, str] = field(default_factory=dict)
    secondary_to_primary: dict[str, str] = field(default_factory=dict)
    labels: dict[str, str] = field(default_factory=dict)
    stars: dict[str, int] = field(default_factory=dict)
    manifest: dict = field(default_factory=dict)

    def exact(self, inchikey: str | None) -> tuple[str, ...]:
        """Every ChEBI id whose structure is exactly this InChIKey.

        A tuple, not a single id: 1,453 of 188,378 distinct keys are held by more
        than one ChEBI entry, so picking one silently hides the rest.
        """
        return self.by_inchikey.get(inchikey, ()) if inchikey else ()

    def skeleton(self, inchikey: str | None) -> tuple[str, ...]:
        """Every ChEBI id sharing this connectivity, stereo and charge ignored."""
        return self.by_skeleton.get(inchikey[:14], ()) if inchikey else ()

    def cas(self, cas: str | None) -> str | None:
        return self.by_cas.get(cas.strip()) if cas else None

    def resolve(self, chebi_id: str) -> str:
        """The primary id for a possibly-secondary ChEBI id.

        ChEBI's flat files carry only primary ids, so an identifier taken from an
        older document resolves to nothing unless it is mapped first. Handoff case
        5: one ChEBI identifier quoted in the original planning document does not
        exist, which is what this looks like.
        """
        numeric = chebi_id.split(":")[-1].strip()
        primary = self.secondary_to_primary.get(numeric)
        return f"CHEBI:{primary}" if primary else f"CHEBI:{numeric}"

    def label(self, chebi_id: str) -> str:
        return self.labels.get(chebi_id.split(":")[-1].strip(), "")

    def __len__(self) -> int:
        return len(self.by_inchikey)


EMPTY = Index()


def discover(listing: str, *, url: str = FLAT_FILES_URL) -> dict[str, str]:
    """Match the tables the gate needs in a directory listing.

    Raises when a required table is absent, with the listing in the message: an
    incomplete release has to stop a run, not weaken it silently.
    """
    names = set(_HREF.findall(listing))
    found: dict[str, str] = {}
    missing: list[str] = []
    for kind, filename in WANTED_TABLES.items():
        for candidate in (filename, f"{filename}.gz"):
            if candidate in names:
                found[kind] = candidate
                break
        else:
            missing.append(filename)

    required = [m for m in missing if _kind_of(m) not in OPTIONAL_TABLES]
    if required:
        raise ReleaseError(
            f"{url} is missing required table(s) {required}. "
            f"Listing held: {sorted(n for n in names if n.endswith(('.tsv', '.tsv.gz')))}"
        )
    for name in missing:
        log.warning("optional table %s absent from %s", name, url)
    return found


def _kind_of(filename: str) -> str:
    for kind, wanted in WANTED_TABLES.items():
        if wanted == filename:
            return kind
    return ""


def _open(path: Path):
    """Open a release table, transparently gunzipping."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open(newline="", encoding="utf-8", errors="replace")


def _reader(path: Path):
    """A csv reader over a release table, plus its lower-cased header."""
    _raise_field_size_limit()
    handle = _open(path)
    reader = csv.reader(handle, delimiter="\t", quotechar='"')
    try:
        header = [h.strip().lower() for h in next(reader)]
    except StopIteration as exc:
        handle.close()
        raise ReleaseError(f"{path.name} is empty") from exc
    return handle, reader, header


def _column(path: Path, header: list[str], *names: str) -> int:
    for name in names:
        if name in header:
            return header.index(name)
    raise ReleaseError(f"{path.name} has none of the columns {names}; header {header}")


def _resolve(directory: Path, filename: str) -> Path:
    """A release table on disk, gzipped or not."""
    for candidate in (directory / filename, directory / f"{filename}.gz"):
        if candidate.exists():
            return candidate
    raise ReleaseError(f"{filename} not found in {directory}")


def distil(
    release_dir: str | Path, index_dir: str | Path, *, generated: str = ""
) -> Index:
    """Build the queryable index from a release directory and write it out.

    Rows are written sorted, so distilling the same release twice produces
    byte-identical index files.
    """
    release_dir = Path(release_dir)
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)

    sources: dict[str, dict] = {}
    labels, stars = _read_compounds(release_dir, sources)
    structures = _read_structures(release_dir, sources)
    cas_rows = _read_cas(release_dir, sources)
    secondary = _read_secondary(release_dir, sources)

    _write(
        index_dir / INDEX_FILES["inchikey"],
        INCHIKEY_COLUMNS,
        # Deduplicated, not merely sorted. Nothing guarantees one structure row per
        # compound across releases -- `_read_structures` indexes every row with a
        # key and ignores `default_structure` -- and a repeated (id, key) pair makes
        # Index.exact return the same id twice, so EXT-02 names it twice in one
        # finding. Release 255 happens to have none; that is a measurement, not a
        # property of the format.
        sorted(set(structures)),
    )
    _write(
        index_dir / INDEX_FILES["compound"],
        COMPOUND_COLUMNS,
        sorted(
            (chebi_id, str(stars.get(chebi_id, "")), label)
            for chebi_id, label in labels.items()
        ),
    )
    _write(index_dir / INDEX_FILES["cas"], CAS_COLUMNS, sorted(cas_rows.items()))
    _write(
        index_dir / INDEX_FILES["secondary"],
        SECONDARY_COLUMNS,
        sorted(secondary.items()),
    )

    manifest = {
        "release_url": FLAT_FILES_URL,
        "release_dir": str(release_dir),
        "generated": generated,
        "sources": sources,
        "index": {
            name: _file_stats(index_dir / filename)
            for name, filename in INDEX_FILES.items()
        },
        "measured": {
            "structures_with_inchikey": len(structures),
            "distinct_cas": len(cas_rows),
            "secondary_ids": len(secondary),
            "compounds": len(labels),
        },
    }
    (index_dir / INDEX_MANIFEST).write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    log.info(
        "distilled %d structures, %d CAS, %d secondary ids into %s",
        len(structures),
        len(cas_rows),
        len(secondary),
        index_dir,
    )
    return load_index(index_dir)


def _read_compounds(
    release_dir: Path, sources: dict
) -> tuple[dict[str, str], dict[str, int]]:
    path = _resolve(release_dir, WANTED_TABLES["compounds"])
    handle, reader, header = _reader(path)
    labels: dict[str, str] = {}
    star_map: dict[str, int] = {}
    rows = 0
    try:
        i_id = _column(path, header, "id", "compound_id")
        # ascii_name first: `name` carries <i>/<sub> markup in 28,243 of 218,533
        # rows and ascii_name in 1,349, and ascii_name is never blank.
        i_name = _column(path, header, "ascii_name", "name")
        i_stars = header.index("stars") if "stars" in header else None
        for row in reader:
            if len(row) <= i_id:
                continue
            chebi_id = row[i_id].strip()
            if not chebi_id:
                continue
            rows += 1
            if len(row) > i_name:
                labels[chebi_id] = row[i_name].strip()
            if i_stars is not None and len(row) > i_stars:
                raw = row[i_stars].strip()
                if raw.isdigit():
                    star_map[chebi_id] = int(raw)
    finally:
        handle.close()
    sources["compounds"] = {**_file_stats(path), "rows": rows}
    return labels, star_map


def _read_structures(release_dir: Path, sources: dict) -> list[tuple[str, str]]:
    path = _resolve(release_dir, WANTED_TABLES["structures"])
    handle, reader, header = _reader(path)
    out: list[tuple[str, str]] = []
    rows = blank = 0
    try:
        i_cid = _column(path, header, "compound_id")
        i_key = _column(path, header, "standard_inchi_key")
        for row in reader:
            if len(row) <= max(i_cid, i_key):
                continue
            rows += 1
            key = row[i_key].strip()
            if not key:
                blank += 1
                continue
            out.append((row[i_cid].strip(), key))
    finally:
        handle.close()
    sources["structures"] = {**_file_stats(path), "rows": rows, "without_key": blank}
    return out


def _read_cas(release_dir: Path, sources: dict) -> dict[str, str]:
    path = _resolve(release_dir, WANTED_TABLES["database_accession"])
    handle, reader, header = _reader(path)
    out: dict[str, str] = {}
    rows = cas_rows = 0
    try:
        i_cid = _column(path, header, "compound_id")
        i_type = _column(path, header, "type")
        i_acc = _column(path, header, "accession_number")
        for row in reader:
            if len(row) <= max(i_cid, i_type, i_acc):
                continue
            rows += 1
            if row[i_type].strip().upper() != CAS_ACCESSION_TYPE:
                continue
            cas_rows += 1
            # First writer wins, so the index is stable under file order. 438 CAS
            # values map to more than one compound, which is why CAS is evidence
            # here and not an identity.
            out.setdefault(row[i_acc].strip(), row[i_cid].strip())
    finally:
        handle.close()
    sources["database_accession"] = {
        **_file_stats(path),
        "rows": rows,
        "cas_rows": cas_rows,
    }
    return out


def _read_secondary(release_dir: Path, sources: dict) -> dict[str, str]:
    try:
        path = _resolve(release_dir, WANTED_TABLES["secondary_ids"])
    except ReleaseError:
        log.warning(
            "secondary_ids table absent; an identifier given to us cannot be "
            "resolved to its primary"
        )
        return {}
    handle, reader, header = _reader(path)
    out: dict[str, str] = {}
    rows = 0
    try:
        i_primary = _column(path, header, "compound_id")
        i_secondary = _column(path, header, "secondary_id")
        for row in reader:
            if len(row) <= max(i_primary, i_secondary):
                continue
            rows += 1
            out[row[i_secondary].strip()] = row[i_primary].strip()
    finally:
        handle.close()
    sources["secondary_ids"] = {**_file_stats(path), "rows": rows}
    return out


def _write(path: Path, columns: tuple[str, ...], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def _file_stats(path: Path) -> dict:
    """Size and hash, so a run manifest can pin exactly what was read."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return {
        "name": path.name,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def load_index(index_dir: str | Path) -> Index:
    """Load a distilled index.

    ``by_cas``, ``stars`` and ``secondary_to_primary`` are built and hashed but no
    check consults them yet. They are kept deliberately: the CAS index is what a
    future EXT check needs to ask "is this CAS already in ChEBI" without a
    structure, and `secondary_ids.tsv` is the table whose absence made the
    prototype's merged-entry reasoning wrong. Distilling them now means the
    question is answerable against an index already in hand rather than after a
    re-download.
    """
    index_dir = Path(index_dir)
    manifest_path = index_dir / INDEX_MANIFEST
    if not manifest_path.exists():
        raise ReleaseError(
            f"no ChEBI index in {index_dir}; build one with distil() from a "
            f"release downloaded from {FLAT_FILES_URL}"
        )
    manifest = json.loads(manifest_path.read_text())
    _verify_index(index_dir, manifest)

    by_inchikey: dict[str, list[str]] = defaultdict(list)
    by_skeleton: dict[str, list[str]] = defaultdict(list)
    for row in _rows(index_dir / INDEX_FILES["inchikey"], INCHIKEY_COLUMNS):
        key = row["inchikey"]
        qualified = f"CHEBI:{row['chebi_id']}"
        by_inchikey[key].append(qualified)
        by_skeleton[key[:14]].append(qualified)

    labels: dict[str, str] = {}
    stars: dict[str, int] = {}
    for row in _rows(index_dir / INDEX_FILES["compound"], COMPOUND_COLUMNS):
        if row["ascii_name"]:
            labels[row["chebi_id"]] = row["ascii_name"]
        if row["stars"].isdigit():
            stars[row["chebi_id"]] = int(row["stars"])

    by_cas = {
        row["cas"]: f"CHEBI:{row['chebi_id']}"
        for row in _rows(index_dir / INDEX_FILES["cas"], CAS_COLUMNS)
    }
    secondary = {
        row["secondary_id"]: row["primary_id"]
        for row in _rows(index_dir / INDEX_FILES["secondary"], SECONDARY_COLUMNS)
    }

    return Index(
        by_inchikey={k: tuple(v) for k, v in by_inchikey.items()},
        by_skeleton={k: tuple(v) for k, v in by_skeleton.items()},
        by_cas=by_cas,
        secondary_to_primary=secondary,
        labels=labels,
        stars=stars,
        manifest=manifest,
    )


def _verify_index(index_dir: Path, manifest: dict) -> None:
    """Hash the index files being read, and refuse one that is not what it claims.

    ``manifest["index"]`` is written by :func:`distil` and was reported to the run
    manifest verbatim, so the hash in a run manifest was the hash of whatever was
    distilled -- not of the bytes the run actually queried. An index edited after
    distillation reported its original hash, and the run manifest, whose entire
    purpose is pinning what was judged against what, pinned a file that no longer
    existed.

    Refusing rather than re-reporting follows the gate's own rule for a malformed
    input: if the evidence is not the evidence the manifest names, every EXT-02 and
    EXT-03 verdict drawn from it is attributed to something that was not consulted.
    An index with no recorded stats -- one distilled before this existed -- is
    measured and filled in rather than refused, because nothing was claimed about
    it to contradict.
    """
    recorded = manifest.get("index")
    if not isinstance(recorded, dict):
        recorded = {}
        manifest["index"] = recorded
    for name, filename in INDEX_FILES.items():
        path = index_dir / filename
        if not path.exists():
            raise ReleaseError(f"index file missing: {path}")
        measured = _file_stats(path)
        claimed = recorded.get(name)
        if not claimed:
            recorded[name] = measured
            continue
        if claimed.get("sha256") != measured["sha256"]:
            raise ReleaseError(
                f"{path} does not match {INDEX_MANIFEST}: it records sha256 "
                f"{claimed.get('sha256')} and {claimed.get('bytes')} bytes, the "
                f"file is {measured['sha256']} and {measured['bytes']} bytes. "
                "Re-distil the release rather than editing an index in place."
            )


def _rows(path: Path, columns: tuple[str, ...]):
    if not path.exists():
        raise ReleaseError(f"index file missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header != columns:
            raise ReleaseError(f"{path.name} columns are {header}, expected {columns}")
        yield from reader
