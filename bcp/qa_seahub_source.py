"""
Well indexes for cross-bucket SeaHub trimming checks.

Builds a common ``(wafer, UG)`` index over each side of the comparison: the
untrimmed vendor delivery and the lab's trimmed upload.  These live in different
buckets with different layouts — note the group folder sits *before* ``raw`` on
the vendor side, so the two have different path depths::

    trimmed  s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/
                 {sublibrary}/{wafer}/{stem}.trim.*
    vendor   s3://czi-novogene/{project}/{order}/{ExperimentID}_{sublibrary}/raw/
                 {wafer}/{stem}.cram

Identity key
------------
Wells are matched on ``(wafer, UG)``.  Measured on a real delivery, each wafer
carries 96 CRAMs with 96 distinct UGs, so the pair is unique.  Matching on the
full filename would fail on every real defect seen in review (duplicated wafer
token, missing sublibrary type, vendor-versus-lab group formatting), whereas
wafer and UG survive all of them.  Barcode, sublibrary, and type are therefore
carried as verification fields, and :mod:`qa_seahub_recon` reports a mismatch in
them instead of letting them break the match.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse
import json
import tempfile

from qa_constants import SEAHUB_STEM_NO_TYPE_RE, SEAHUB_STEM_RE
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family

__all__ = [
    "IdentityKey",
    "SourceEntry",
    "TrimmedEntry",
    "index_trimmed_upload",
    "index_untrimmed_source",
    "load_source_read_counts",
    "parse_source_uri",
]

SOURCE_METADATA_MAX_WORKERS = 16

# Vendor wafer-level and delivery-level files that are not per-well artifacts.
_SOURCE_SKIP_SUFFIXES = (
    "_LibraryInfo.xml",
    "_UploadCompleted.json",
    "_merged_trimmer-stats.csv",
    "_merged_trimmer-failure_codes.csv",
    "_run_SecondaryAnalysis.txt",
    "_run_VariantCalling.txt",
    "file-manifest.json",
)

IdentityKey = tuple[str, str]


@dataclass
class SourceEntry:
    """One untrimmed vendor CRAM and its sidecar."""

    wafer: str
    ug: str
    barcode: str
    group: str
    assay: str | None
    cram_key: str
    metadata_key: str = ""
    size_bytes: int = 0
    read_count: int | None = None


@dataclass
class TrimmedEntry:
    """One well in the lab's trimmed upload."""

    wafer: str
    ug: str
    barcode: str
    sublibrary: str
    assay: str | None
    stem: str
    family: str
    raw_dir: str
    has_cram: bool = False
    read_count: int | None = None


def parse_source_uri(uri: str) -> tuple[str, str]:
    """Split an ``s3://bucket/prefix`` URI into bucket and listing prefix."""
    parsed = urlparse(uri.strip())
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError(f"Invalid untrimmed source path (expected s3://...): {uri!r}")
    prefix = parsed.path.strip("/")
    return parsed.netloc, f"{prefix}/" if prefix else ""


def _is_skippable_source_key(key: str) -> bool:
    name = key.split("/")[-1]
    if "unmatched" in name:
        return True
    return name.endswith(_SOURCE_SKIP_SUFFIXES)


def _parse_stem_fields(stem: str) -> dict[str, str | None] | None:
    """Parse a stem into wafer / group / assay / UG / barcode.

    Falls back to the type-less pattern so misnamed uploads still match.
    """
    match = SEAHUB_STEM_RE.match(stem)
    assay: str | None = None
    if match is not None:
        assay = match.group("assay")
    else:
        match = SEAHUB_STEM_NO_TYPE_RE.match(stem)
        if match is None:
            return None
    return {
        "wafer": match.group("wafer"),
        "group": match.group("group"),
        "assay": assay,
        "ug": match.group("ug"),
        "barcode": match.group("barcode"),
    }


def index_untrimmed_source(s3_client: Any, uri: str) -> dict[IdentityKey, SourceEntry]:
    """List the vendor delivery and index per-well CRAMs by ``(wafer, UG)``."""
    bucket, prefix = parse_source_uri(uri)
    paginator = s3_client.get_paginator("list_objects")
    index: dict[IdentityKey, SourceEntry] = {}
    metadata_keys: dict[IdentityKey, str] = {}

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if _is_skippable_source_key(key):
                continue
            name = key.split("/")[-1]
            is_metadata = name.endswith(".cram-metadata.json")
            if not is_metadata and not name.endswith(".cram"):
                continue
            stem = (
                name[: -len(".cram-metadata.json")]
                if is_metadata
                else name[: -len(".cram")]
            )
            fields = _parse_stem_fields(stem)
            if fields is None:
                continue
            identity = (str(fields["wafer"]), str(fields["ug"]))
            if is_metadata:
                metadata_keys[identity] = key
                continue
            index[identity] = SourceEntry(
                wafer=str(fields["wafer"]),
                ug=str(fields["ug"]),
                barcode=str(fields["barcode"]),
                group=str(fields["group"]),
                assay=fields["assay"],
                cram_key=key,
                size_bytes=int(obj.get("Size", 0) or 0),
            )

    for identity, metadata_key in metadata_keys.items():
        if identity in index:
            index[identity].metadata_key = metadata_key
    return index


def load_source_read_counts(
    s3_client: Any, bucket: str, index: dict[IdentityKey, SourceEntry]
) -> None:
    """Populate ``read_count`` on each entry from its vendor sidecar."""
    targets = [
        (identity, entry.metadata_key)
        for identity, entry in index.items()
        if entry.metadata_key
    ]
    if not targets:
        return

    def _fetch(identity: IdentityKey, key: str) -> tuple[IdentityKey, int | None]:
        with tempfile.NamedTemporaryFile(
            mode="w+b", delete=False, suffix=".json"
        ) as tf:
            local = tf.name
        try:
            s3_client.download_file(bucket, key, local)
            with open(local) as fh:
                payload = json.load(fh)
            value = payload.get("read_count")
            return identity, int(value) if value is not None else None
        finally:
            Path(local).unlink(missing_ok=True)

    workers = min(SOURCE_METADATA_MAX_WORKERS, len(targets))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_fetch, identity, key) for identity, key in targets]
        for future in as_completed(futures):
            identity, read_count = future.result()
            index[identity].read_count = read_count


def index_trimmed_upload(
    all_raw_files: list[str],
    read_metadata: dict[str, dict] | None = None,
    bucket: str = "",
) -> dict[IdentityKey, TrimmedEntry]:
    """Index a gathered SeaHub listing by ``(wafer, UG)``, both families."""
    index: dict[IdentityKey, TrimmedEntry] = {}
    for key in sorted(all_raw_files):
        parsed = seahub_stem_and_family(key.split("/")[-1])
        path_info = parse_seahub_raw_path(key)
        if parsed is None or path_info is None:
            continue
        stem, suffix, family = parsed
        fields = _parse_stem_fields(stem)
        if fields is None:
            continue
        identity = (str(fields["wafer"]), str(fields["ug"]))
        entry = index.get(identity)
        if entry is None:
            entry = TrimmedEntry(
                wafer=str(fields["wafer"]),
                ug=str(fields["ug"]),
                barcode=str(fields["barcode"]),
                sublibrary=path_info["sublibrary"],
                assay=fields["assay"],
                stem=stem,
                family=family,
                raw_dir="/".join(key.split("/")[:-1]),
            )
            index[identity] = entry
        if suffix in (".cram", ".trim.cram"):
            entry.has_cram = True
            if read_metadata:
                entry.read_count = _lookup_read_count(read_metadata, bucket, key)
    return index


def _lookup_read_count(
    read_metadata: dict[str, dict], bucket: str, cram_key: str
) -> int | None:
    """Find a CRAM's read_count, tolerating both metadata key conventions."""
    candidates = [
        f"s3://{bucket}/{cram_key}" if bucket else "",
        cram_key.split("/")[-1],
    ]
    for candidate in candidates:
        if candidate and candidate in read_metadata:
            value = read_metadata[candidate].get("read_count")
            return int(value) if value is not None else None
    return None
