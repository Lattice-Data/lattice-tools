"""
Well indexes for cross-bucket SeaHub trimming checks.

Builds a common ``(wafer, UG)`` index over each side of the comparison: the
untrimmed vendor delivery and the lab's trimmed upload.  These live in different
buckets with different layouts -- note the ExperimentID folder sits *before*
``raw`` on the vendor side, and the vendor path carries no sublibrary at all::

    trimmed  s3://czi-{lab}/{lastname}-{projectname}/{ExperimentID}/raw/
                 {sublibrary}/{wafer}/{stem}.trim.*
    vendor   s3://czi-novogene/{project}/{order}/{ExperimentID}/raw/
                 {wafer}/{stem}.cram

Measured against the real REF3 delivery (2597 objects under order
``NVUS2024101701-11``): six wafers, each carrying exactly 48 CRAMs with 48
distinct UGs, one sublibrary per wafer, every well typed ``GEX_hash_oligo``.

Identity key
------------
Wells are matched on ``(wafer, UG)``, which that delivery confirms is unique.
Matching on the full filename would fail on every real defect seen in review
(duplicated wafer token, missing sublibrary type, vendor-versus-lab group
formatting), whereas wafer and UG survive all of them.  Barcode, sublibrary and
type are therefore carried as verification fields, and :mod:`qa_seahub_recon`
reports a mismatch in them instead of letting them break the match.

Several sources
---------------
One experiment spans several vendor orders, and one order can hold several
experiments, so the untrimmed side is a *list* of prefixes.  REF3 makes this
concrete: order ``NVUS2024101701-11`` contains six of its seven sublibraries, and
``REF3_P05_1`` -- the only correctly-named one in the trimmed upload -- is not
there at all.  Listing that order alone would report all of its wells as
orphans, so :class:`SourceCoverage` exists to make incomplete input obvious
rather than letting it read as a completeness failure.

Because the sublibrary is absent from the vendor path, the authoritative
sublibrary name comes from the vendor *stem* with its trailing well token
stripped (``438514-REF3_P07_1_A3_GEX...`` gives ``REF3_P07_1``), not from any
path segment.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse
import json
import tempfile

from qa_constants import (
    SEAHUB_RAW_SEGMENT,
    SEAHUB_SEARCH_MAX_DEPTH,
    SEAHUB_SEARCH_MAX_LISTINGS,
    SEAHUB_STEM_NO_TYPE_RE,
    SEAHUB_STEM_RE,
    SEAHUB_UNKNOWN_ORDER_LABEL,
    SEAHUB_VENDOR_ORDER_RE,
    SEAHUB_WAFER_RE,
)
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family

__all__ = [
    "IdentityKey",
    "RawPrefix",
    "RawPrefixScan",
    "SourceCoverage",
    "SourceEntry",
    "TrimmedEntry",
    "UntrimmedSources",
    "WaferSearchPlan",
    "WaferSeeds",
    "derive_source_order",
    "discover_untrimmed_sources",
    "finding_row",
    "index_trimmed_upload",
    "index_untrimmed_sources",
    "load_source_read_counts",
    "normalize_search_roots",
    "normalize_source_uris",
    "parse_seahub_stem_fields",
    "parse_source_uri",
    "seahub_wafer_seeds",
    "source_order_by_wafer",
]

SOURCE_METADATA_MAX_WORKERS = 16

_CRAM_SUFFIX = ".cram"
_SIDECAR_SUFFIX = ".cram-metadata.json"

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
    # Provenance, so every finding can name the order the well came from.
    bucket: str = ""
    source_uri: str = ""
    source_order: str = ""

    @property
    def s3_uri(self) -> str:
        return f"s3://{self.bucket}/{self.cram_key}" if self.bucket else self.cram_key


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
    # No read_count here. The reconciliation compares the *vendor* count
    # (SourceEntry.read_count) with the trimmer's declared input totals, so a
    # trimmed-side count was filled in at one metadata lookup per CRAM and never
    # read. Dropping it is what lets index_trimmed_upload take neither
    # read_metadata nor bucket.
    # 0 means "unknown", not "empty": sizes are only collected in S3 mode.
    size_bytes: int = 0


@dataclass(frozen=True)
class WaferSeeds:
    """The wafers a vendor search will look for, and the tokens it refused."""

    wafers: tuple[str, ...] = ()
    # Tokens that reached the seed builder but are not wafer-shaped, so a
    # malformed upload cannot turn into a search for a garbage term.
    rejected: tuple[str, ...] = ()

    def __len__(self) -> int:
        return len(self.wafers)


@dataclass
class WaferSearchPlan:
    """What a wafer search looked for, where it found it, and what it missed.

    Carried on :class:`UntrimmedSources` so discovery has one return value and
    :func:`~qa_seahub_recon.reconcile_trimming` -- which reads only ``index``,
    ``coverage`` and ``findings`` -- needs no changes at all.
    """

    roots: tuple[str, ...] = ()
    seeds: WaferSeeds = field(default_factory=WaferSeeds)
    # wafer -> the raw/ prefix URIs it was found under. More than one is the
    # collision case, and is never silently collapsed to a choice.
    located: dict[str, tuple[str, ...]] = field(default_factory=dict)
    # (uri, experiment_segment) per delivery the descent saw and did not index,
    # because no uploaded wafer pointed at it. Free -- the walk enumerated them
    # anyway -- and the only signal that a delivery nothing points at exists.
    unseeded: tuple[tuple[str, str], ...] = ()
    listings: int = 0
    unreadable: tuple[tuple[str, str], ...] = ()
    budget_exhausted: bool = False
    depth_capped: tuple[str, ...] = ()

    @property
    def not_located(self) -> tuple[str, ...]:
        return tuple(w for w in self.seeds.wafers if w not in self.located)

    @property
    def collided(self) -> tuple[str, ...]:
        return tuple(sorted(w for w, at in self.located.items() if len(at) > 1))

    @property
    def complete(self) -> bool:
        """False when something stopped the search seeing everything it should."""
        return not (
            self.unreadable
            or self.budget_exhausted
            or self.depth_capped
            or self.not_located
        )

    def summary(self) -> str:
        """The one line an operator reads before opening any CSV."""
        parts = [
            f"located {len(self.located)} of {len(self.seeds.wafers)} seed wafer(s) "
            f"in {self.listings} listing(s)"
        ]
        if self.not_located:
            parts.append(f"not found: {', '.join(self.not_located)}")
        if self.collided:
            parts.append(f"{len(self.collided)} wafer(s) in more than one delivery")
        if self.unseeded:
            named = ", ".join(sorted({seg or "?" for _uri, seg in self.unseeded}))
            parts.append(
                f"{len(self.unseeded)} delivery(ies) under these roots were not "
                f"indexed because no uploaded wafer matched ({named})"
            )
        if self.seeds.rejected:
            parts.append(f"{len(self.seeds.rejected)} token(s) were not wafer-shaped")
        if self.unreadable:
            parts.append(f"{len(self.unreadable)} prefix(es) could not be listed")
        if self.budget_exhausted:
            parts.append("listing budget exhausted -- the search is incomplete")
        if self.depth_capped:
            parts.append(f"{len(self.depth_capped)} prefix(es) hit the depth cap")
        return "; ".join(parts)


@dataclass
class SourceCoverage:
    """What one untrimmed prefix contributed, so partial input is visible."""

    source_uri: str
    bucket: str
    prefix: str
    source_order: str = ""
    cram_keys: int = 0
    indexed: int = 0
    duplicate_losses: int = 0
    orders_seen: tuple[str, ...] = ()
    skipped_reason: str = ""
    # Filled in by the reconciliation, which is the only place that knows it.
    matched: int = 0
    # How this prefix came to be listed: "operator" for one that was handed in,
    # "seed" for one a wafer in the upload pointed at, "sibling" for a delivery
    # found beside a seeded one. A not_trimmed count that moves after a vendor
    # reorganisation is then explainable rather than mysterious.
    found_by: str = "operator"
    # The path segment before ``raw``, advisory, for naming a delivery in a report.
    experiment_segment: str = ""

    @property
    def unmatched(self) -> int:
        return max(0, self.indexed - self.matched)


@dataclass
class UntrimmedSources:
    """The merged vendor index plus per-prefix coverage and findings."""

    index: dict[IdentityKey, SourceEntry] = field(default_factory=dict)
    coverage: list[SourceCoverage] = field(default_factory=list)
    findings: list[dict[str, Any]] = field(default_factory=list)
    # None for the explicit-prefix path, so nothing that reads coverage or
    # findings has to know discovery exists.
    search: WaferSearchPlan | None = None

    def __len__(self) -> int:
        return len(self.index)

    @property
    def buckets(self) -> tuple[str, ...]:
        return tuple(sorted({e.bucket for e in self.index.values() if e.bucket}))


def finding_row(category: str, **fields: Any) -> dict[str, Any]:
    """Build one reconciliation row.

    Lives here rather than in :mod:`qa_seahub_recon` so the indexers can emit
    findings too; recon already imports this module, so the factory has to sit on
    this side to stay acyclic.
    """
    row: dict[str, Any] = {
        "category": category,
        "wafer": "",
        "ug": "",
        "sublibrary": "",
        "source_key": "",
        "trimmed_stem": "",
        "detail": "",
    }
    row.update(fields)
    return row


def parse_source_uri(uri: str) -> tuple[str, str]:
    """Split an ``s3://bucket/prefix`` URI into bucket and listing prefix."""
    parsed = urlparse(uri.strip())
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError(f"Invalid untrimmed source path (expected s3://...): {uri!r}")
    prefix = parsed.path.strip("/")
    return parsed.netloc, f"{prefix}/" if prefix else ""


def normalize_source_uris(uris: Any) -> list[str]:
    """Accept a bare string or a list, and validate each entry.

    A bucket-only URI is rejected rather than accepted: ``parse_source_uri``
    yields an empty prefix for it, which would paginate an entire bucket. A
    string is never comma-split -- an S3 prefix may legitimately contain a comma,
    and silently splitting one would list the wrong thing.
    """
    if uris is None:
        return []
    candidates = [uris] if isinstance(uris, str) else list(uris)

    normalized: list[str] = []
    for candidate in candidates:
        uri = str(candidate).strip()
        if not uri:
            continue
        _bucket, prefix = parse_source_uri(uri)
        if len([s for s in prefix.split("/") if s]) < 2:
            raise ValueError(
                f"Untrimmed source {uri!r} is too broad: expected at least "
                "s3://bucket/project/order, otherwise the whole bucket is listed"
            )
        if uri not in normalized:
            normalized.append(uri)
    return normalized


def normalize_search_roots(roots: Any) -> list[str]:
    """Validate the roots a wafer search may descend from; a sibling of the above.

    :func:`normalize_source_uris` is deliberately left alone rather than given a
    flag. Its "too broad" guard exists because the prefix it validates is handed
    straight to a flat object paginate, so a shallow one walks a whole bucket by
    accident -- and that is still true of every caller it has. A search root is a
    different kind of thing: it is descended with a delimiter, level by level, and
    a bucket root is the point of it. Two functions keep the two meanings apart;
    one function with a flag would have re-opened the accident for the old
    parameter.

    Accepted: any depth, including a bare bucket. ``s3://czi-novogene`` and
    ``s3://czi-novogene/`` normalize identically, so a trailing slash is never
    the difference between two roots. Rejected:

    * anything that is not ``s3://bucket[/prefix]`` -- as before, via
      :func:`parse_source_uri`.
    * a root containing a ``raw`` segment. That is a listing prefix, not a search
      root -- the shape the retired ``untrimmed_s3_paths`` took, so it is what an
      operator with muscle memory will paste. The descent looks for a
      child folder named ``raw``, so starting inside one finds nothing and would
      report zero located wafers, which is indistinguishable from a vendor bucket
      that genuinely holds none. Better to say so.

    A root inside another root in the same bucket is dropped, because descending
    the outer one already covers it. A string is never comma-split, for the same
    reason as above: an S3 prefix may contain a comma.
    """
    if roots is None:
        return []
    candidates = [roots] if isinstance(roots, str) else list(roots)

    parsed: list[tuple[str, str, str]] = []
    for candidate in candidates:
        uri = str(candidate).strip()
        if not uri:
            continue
        bucket, prefix = parse_source_uri(uri)
        if SEAHUB_RAW_SEGMENT in [s for s in prefix.split("/") if s]:
            raise ValueError(
                f"Search root {uri!r} is inside a {SEAHUB_RAW_SEGMENT!r} folder: "
                "give the project or bucket to search, not one delivery's listing "
                "prefix"
            )
        canonical = f"s3://{bucket}/{prefix}" if prefix else f"s3://{bucket}/"
        if canonical not in [c for c, _b, _p in parsed]:
            parsed.append((canonical, bucket, prefix))

    # Shortest prefix first, so a parent is always seen before its children.
    ordered = sorted(parsed, key=lambda item: (item[1], len(item[2]), item[2]))
    kept: list[tuple[str, str, str]] = []
    for canonical, bucket, prefix in ordered:
        if any(b == bucket and prefix.startswith(p) for _c, b, p in kept):
            continue
        kept.append((canonical, bucket, prefix))
    return [canonical for canonical, _b, _p in kept]


def derive_source_order(cram_key: str) -> str:
    """Read the vendor order out of a key, positionally.

    ``{project}/{order}/{ExperimentID}/raw/{wafer}/{file}`` -- so the order is
    the second segment. Returns ``""`` for any other shape, which degrades the
    coverage label without ever dropping the entry.
    """
    parts = cram_key.split("/")
    if len(parts) == 6 and parts[3] == "raw":
        return parts[1]
    return ""


def _order_label_from_prefix(prefix: str) -> str:
    """Best-effort order label for a prefix that produced no objects.

    Scans right to left for a segment shaped like a vendor order. Deliberately
    not ``segments[-1]``: for ``.../NVUS2024101701-11/REF3`` the last segment is
    the ExperimentID, and mislabelling coverage rows is worse than not labelling.
    """
    for segment in reversed([s for s in prefix.split("/") if s]):
        if SEAHUB_VENDOR_ORDER_RE.match(segment):
            return segment
    return SEAHUB_UNKNOWN_ORDER_LABEL


def _is_skippable_source_key(key: str) -> bool:
    name = key.split("/")[-1]
    if "unmatched" in name:
        return True
    return name.endswith(_SOURCE_SKIP_SUFFIXES)


def parse_seahub_stem_fields(stem: str) -> dict[str, str | None] | None:
    """Parse a stem into wafer / group / assay / UG / barcode.

    Falls back to the type-less pattern so misnamed uploads still match.

    Note that for a doubled-wafer stem the returned ``group`` still contains the
    second wafer: ``438514-438514-REF3_P07_1_A3-Z0305-...`` yields group
    ``438514-REF3_P07_1_A3``. Callers that care about the group must normalize
    the wafer first (:func:`qa_seahub_sop.normalize_doubled_wafer`).
    ``index_trimmed_upload`` is unaffected only because it takes the sublibrary
    from the folder.
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


def _dedupe_overlapping(
    uris: list[str],
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str, str]]]:
    """Drop prefixes already covered by a shorter prefix in the same bucket.

    Returns ``(survivors, skipped)`` as ``(uri, bucket, prefix)`` tuples, with
    ``skipped`` carrying the covering URI. Comparison is on the normalized
    prefix, which always ends in one slash, so ``.../REF3/`` cannot falsely
    cover ``.../REF3_P05_1/``.
    """
    parsed = [(uri, *parse_source_uri(uri)) for uri in uris]
    # Shortest prefix first so a parent is always seen before its children.
    ordered = sorted(parsed, key=lambda item: (item[1], len(item[2]), item[2]))

    survivors: list[tuple[str, str, str]] = []
    skipped: list[tuple[str, str, str, str]] = []
    for uri, bucket, prefix in ordered:
        covering = next(
            (
                kept_uri
                for kept_uri, kept_bucket, kept_prefix in survivors
                if kept_bucket == bucket and prefix.startswith(kept_prefix)
            ),
            None,
        )
        if covering is not None:
            skipped.append((uri, bucket, prefix, covering))
            continue
        survivors.append((uri, bucket, prefix))
    return survivors, skipped


@dataclass(frozen=True)
class RawPrefix:
    """One ``.../raw/`` prefix a descent located, and what was directly inside it.

    Both layouts are described by the same record, because a vendor delivery may
    use either and the descent cannot tell which until it looks: ``wafer_folders``
    is populated for the foldered layout, ``loose_objects`` for the flat one where
    the wafer appears only in the filename, and a real prefix can carry both --
    wafer folders plus wafer-level files beside them.
    """

    bucket: str
    prefix: str
    # The path segment immediately before ``raw``, read off the path so it
    # survives a key depth the six-segment positional readers cannot parse --
    # four of the six real vendor layouts. Advisory only: it names a delivery in
    # a report, and never decides whether one is indexed.
    experiment_segment: str = ""
    wafer_folders: tuple[str, ...] = ()
    loose_objects: tuple[str, ...] = ()

    @property
    def uri(self) -> str:
        return f"s3://{self.bucket}/{self.prefix}"


@dataclass
class RawPrefixScan:
    """Every ``raw/`` prefix a descent found, and how completely it looked."""

    prefixes: list[RawPrefix] = field(default_factory=list)
    listings: int = 0
    # (prefix, botocore error code) per listing that was refused, so a partly
    # readable bucket yields a partial answer plus a row rather than an exception.
    unreadable: tuple[tuple[str, str], ...] = ()
    budget_exhausted: bool = False
    # Roots whose descent hit SEAHUB_SEARCH_MAX_DEPTH without finding raw/.
    depth_capped: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        """False when anything stopped the walk seeing everything it should."""
        return not (self.unreadable or self.budget_exhausted or self.depth_capped)


def _descend_to_raw_prefixes(
    s3_client: Any,
    roots: Any,
    max_depth: int = SEAHUB_SEARCH_MAX_DEPTH,
    max_listings: int = SEAHUB_SEARCH_MAX_LISTINGS,
) -> RawPrefixScan:
    """Find the ``raw/`` prefixes under one or more search roots.

    Structural only: it locates deliveries and reports what is in them, and
    decides nothing about which are relevant. Filtering to the wafers an upload
    actually mentions is the caller's job, which is what lets the same walk
    answer both "where is my vendor data" and "what else is here that nothing
    points at".

    Every level is listed with ``Delimiter='/'`` through ``get_paginator``, which
    matters twice. The delimiter is what keeps this O(folders) instead of
    O(objects) -- the flat paginate ``index_untrimmed_sources`` uses would walk
    every object in the bucket, which is the whole reason a bucket root is
    refused there. And the paginator is what stops a level being silently
    truncated at 1000: ``qa_gather._list_folders`` issues a bare ``list_objects``
    and only appends a warning, so a project with more than 1000 orders would
    quietly walk a fraction of itself.

    Three things bound it, and each is reported rather than raised, because a
    partial answer an operator can see beats a dead cell:

    * ``max_depth`` levels below each root. A bucket root needs three to reach
      ``{project}/{order}/{ExperimentID}/raw``; the default tolerates one more.
    * ``max_listings`` calls in total.
    * a refused listing, recorded with its error code and stepped over, so one
      inaccessible project does not cost the rest of the bucket.
    """
    scan = RawPrefixScan()
    unreadable: list[tuple[str, str]] = []
    depth_capped: list[str] = []
    paginator = s3_client.get_paginator("list_objects")

    def _level(bucket: str, prefix: str) -> tuple[list[str], list[str]] | None:
        """One delimiter listing: (child folder prefixes, object keys), or None."""
        if scan.listings >= max_listings:
            scan.budget_exhausted = True
            return None
        scan.listings += 1
        folders: list[str] = []
        objects: list[str] = []
        try:
            for page in paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
                folders += [entry["Prefix"] for entry in page.get("CommonPrefixes", [])]
                objects += [entry["Key"] for entry in page.get("Contents", [])]
        except Exception as exc:  # noqa: BLE001 - botocore ClientError and kin
            unreadable.append((f"s3://{bucket}/{prefix}", _error_code(exc)))
            return None
        return folders, objects

    for root in normalize_search_roots(roots):
        bucket, root_prefix = parse_source_uri(root)
        # (prefix, depth), depth counted from the root.
        pending: list[tuple[str, int]] = [(root_prefix, 0)]
        while pending:
            prefix, depth = pending.pop(0)
            listed = _level(bucket, prefix)
            if listed is None:
                if scan.budget_exhausted:
                    break
                continue
            folders, _objects = listed
            for child in folders:
                if _segment_name(child) != SEAHUB_RAW_SEGMENT:
                    continue
                inside = _level(bucket, child)
                if inside is None:
                    if scan.budget_exhausted:
                        break
                    continue
                wafer_folders, loose = inside
                scan.prefixes.append(
                    RawPrefix(
                        bucket=bucket,
                        prefix=child,
                        experiment_segment=_segment_name(prefix),
                        wafer_folders=tuple(
                            _segment_name(folder) for folder in wafer_folders
                        ),
                        loose_objects=tuple(loose),
                    )
                )
            if scan.budget_exhausted:
                break
            deeper = [
                (child, depth + 1)
                for child in folders
                if _segment_name(child) != SEAHUB_RAW_SEGMENT
            ]
            if deeper and depth + 1 > max_depth:
                depth_capped.append(f"s3://{bucket}/{prefix}")
                continue
            pending += deeper
        if scan.budget_exhausted:
            break

    scan.unreadable = tuple(unreadable)
    scan.depth_capped = tuple(sorted(set(depth_capped)))
    return scan


def _segment_name(prefix: str) -> str:
    """Last path segment of a folder prefix, without its trailing slash."""
    return prefix.rstrip("/").split("/")[-1]


def _error_code(exc: Exception) -> str:
    """Botocore's error code when there is one, else the exception class name."""
    response = getattr(exc, "response", None)
    if isinstance(response, dict):
        code = response.get("Error", {}).get("Code")
        if code:
            return str(code)
    return type(exc).__name__


def _prefer_source_entry(
    a: SourceEntry, b: SourceEntry
) -> tuple[SourceEntry, SourceEntry]:
    """Choose between two vendor CRAMs claiming one ``(wafer, UG)``; ``(winner, loser)``.

    A named predicate rather than an inline sort key, because which copy of a
    re-delivered well QA compares against decides every downstream number for
    that well, and the rule needs somewhere to be read and argued with.

    Order of preference:

    1. A key whose layout exposes an order beats one that does not.  An
       unplaceable key -- wrong depth, no ``raw`` segment -- used to win
       outright, because ``derive_source_order`` returns ``""`` for it and the
       empty string sorts before every real order.
    2. The **newest** order.  This reverses the previous rule, which took the
       lexicographic minimum and so kept the *original* delivery and discarded
       the re-delivery -- backwards for the case duplicates mostly arise from,
       where the vendor re-sent a well because the first copy was wrong.
    3. The lower key, so the answer never depends on listing or caller order.

    Lexicographic comparison stands in for chronological here: a vendor order is
    ``[A-Z]{2,}\\d{6,}-\\d{2,}`` (SEAHUB_VENDOR_ORDER_RE), so a fixed-width date
    leads and orders from one vendor compare correctly.  Two orders whose date
    parts differ in width would not, which is why rule 3 exists rather than
    being unreachable.

    This is deliberately *not* the final word.  The evidence that actually
    settles which copy was trimmed is whether a candidate's ``read_count``
    matches one of the trimmer's declared input totals -- and neither is known
    here: ``load_source_read_counts`` runs after this function returns, and the
    trimmer's totals come from the trimmed side, which is indexed later still.
    Arbitrating on that evidence belongs after reconciliation has both.
    """
    if bool(a.source_order) != bool(b.source_order):
        return (a, b) if a.source_order else (b, a)
    if a.source_order != b.source_order:
        return (a, b) if a.source_order > b.source_order else (b, a)
    return (a, b) if a.cram_key < b.cram_key else (b, a)


def _index_prefix(
    paginator: Any,
    bucket: str,
    prefix: str,
    uri: str,
    sources: UntrimmedSources,
    coverage: SourceCoverage,
    sidecar_by_stem_key: dict[tuple[str, str], str],
    accept: Any = None,
) -> set[str]:
    """Index every per-well CRAM under one listed prefix into ``sources``.

    Lifted out of :func:`index_untrimmed_sources` unchanged so that the object
    loop -- the part wafer-seeded discovery needs to reuse against a prefix it
    located itself, rather than one an operator typed -- has a caller boundary.
    Everything around it stays where it was: the overlap pass, the per-prefix
    findings, and the coverage finalisation all remain in the caller.

    Mutates ``sources.index``, ``sources.findings``, ``coverage.cram_keys`` and
    ``sidecar_by_stem_key``; returns the orders it saw, which the caller needs to
    label the coverage row.

    ``accept(bucket, key, fields)`` filters CRAMs before they are indexed. The
    bucket is passed because two buckets can share a key path, so a filter keyed
    on the path alone would test a key against the wrong delivery. Discovery uses
    it to keep a located delivery to the wafers it is there for, so that walking
    a whole project does not drag every other experiment's wells into one
    upload's comparison -- the contamination the ExperimentID filter used to
    prevent, now prevented by the upload's own contents instead of by a name.
    """
    orders_seen: set[str] = set()

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if _is_skippable_source_key(key):
                continue
            name = key.split("/")[-1]
            is_sidecar = name.endswith(_SIDECAR_SUFFIX)
            if not is_sidecar and not name.endswith(_CRAM_SUFFIX):
                continue
            stem = (
                name[: -len(_SIDECAR_SUFFIX)]
                if is_sidecar
                else name[: -len(_CRAM_SUFFIX)]
            )
            fields = parse_seahub_stem_fields(stem)
            if fields is None:
                continue
            if is_sidecar:
                sidecar_by_stem_key[(bucket, key[: -len(_SIDECAR_SUFFIX)])] = key
                continue

            if accept is not None and not accept(bucket, key, fields):
                continue

            order = derive_source_order(key)
            orders_seen.add(order or SEAHUB_UNKNOWN_ORDER_LABEL)
            coverage.cram_keys += 1
            candidate = SourceEntry(
                wafer=str(fields["wafer"]),
                ug=str(fields["ug"]),
                barcode=str(fields["barcode"]),
                group=str(fields["group"]),
                assay=fields["assay"],
                cram_key=key,
                size_bytes=int(obj.get("Size", 0) or 0),
                bucket=bucket,
                source_uri=uri,
                source_order=order,
            )
            identity = (candidate.wafer, candidate.ug)
            incumbent = sources.index.get(identity)
            if incumbent is None:
                sources.index[identity] = candidate
                continue

            winner, loser = _prefer_source_entry(incumbent, candidate)
            sources.index[identity] = winner
            sources.findings.append(
                finding_row(
                    "duplicate_source_well",
                    wafer=candidate.wafer,
                    ug=candidate.ug,
                    source_key=loser.cram_key,
                    detail=(
                        f"well delivered by more than one source: kept "
                        f"{winner.s3_uri} (order {winner.source_order or '?'}), "
                        f"ignored {loser.s3_uri} "
                        f"(order {loser.source_order or '?'})"
                    ),
                )
            )

    return orders_seen


def index_untrimmed_sources(
    s3_client: Any, uris: Any, accept: Any = None
) -> UntrimmedSources:
    """List one or more vendor prefixes and index per-well CRAMs by ``(wafer, UG)``.

    An identity delivered by two prefixes is a re-delivery, not a merge: the
    winner is chosen by :func:`_prefer_source_entry` and the loser becomes a
    ``duplicate_source_well`` finding. Overwriting blindly, as the single-prefix
    version did, silently discarded one of them.

    ``accept(bucket, key, fields)`` filters CRAMs before indexing;
    :func:`discover_untrimmed_sources` uses it to hold a discovered delivery to
    the wafers it was located for.

    Callers hand in concrete listing prefixes, which is why
    :func:`normalize_source_uris` still refuses a shallow one: each surviving
    prefix is paginated *flat*, so a bucket or project root here would walk every
    object under it. Discovery satisfies that guard by construction -- it only
    ever emits a ``.../raw/`` prefix -- which makes the guard an invariant check
    on the descent rather than dead weight.
    """
    sources = UntrimmedSources()
    normalized = normalize_source_uris(uris)
    if not normalized:
        return sources

    survivors, skipped = _dedupe_overlapping(normalized)
    for uri, bucket, prefix, covering in skipped:
        sources.coverage.append(
            SourceCoverage(
                source_uri=uri,
                bucket=bucket,
                prefix=prefix,
                source_order=_order_label_from_prefix(prefix),
                skipped_reason=f"already covered by {covering}",
            )
        )
        sources.findings.append(
            finding_row(
                "overlapping_source_prefix",
                detail=(
                    f"untrimmed source {uri} is inside {covering}; listed once "
                    "under the broader prefix"
                ),
            )
        )

    paginator = s3_client.get_paginator("list_objects")
    # Keyed on (bucket, stem) so a sidecar can only ever attach to its own CRAM,
    # independent of listing order. The bucket belongs in the key because this
    # dict outlives the per-prefix loop: with only the stem, a CRAM in one bucket
    # is credited with a sidecar that exists at the same path in a *different*
    # one, so the entry claims metadata it does not have. load_source_read_counts
    # then fetches from the entry's own bucket, misses, and degrades to
    # metadata_unavailable -- which reads as "the vendor's sidecar has no
    # read_count" rather than "no sidecar was delivered here". Latent while one
    # bucket is listed; likely once prefixes are discovered rather than typed.
    sidecar_by_stem_key: dict[tuple[str, str], str] = {}

    for uri, bucket, prefix in survivors:
        coverage = SourceCoverage(source_uri=uri, bucket=bucket, prefix=prefix)
        orders_seen = _index_prefix(
            paginator,
            bucket,
            prefix,
            uri,
            sources,
            coverage,
            sidecar_by_stem_key,
            accept=accept,
        )

        coverage.orders_seen = tuple(sorted(orders_seen))
        coverage.source_order = (
            coverage.orders_seen[0]
            if len(coverage.orders_seen) == 1
            else _order_label_from_prefix(prefix)
        )
        sources.coverage.append(coverage)

    for identity, entry in sources.index.items():
        sidecar = sidecar_by_stem_key.get(
            (entry.bucket, entry.cram_key[: -len(_CRAM_SUFFIX)])
        )
        if sidecar:
            entry.metadata_key = sidecar

    for coverage in sources.coverage:
        if coverage.skipped_reason:
            continue
        coverage.indexed = sum(
            1 for e in sources.index.values() if e.source_uri == coverage.source_uri
        )
        coverage.duplicate_losses = coverage.cram_keys - coverage.indexed
        if coverage.cram_keys == 0:
            # Emitted as a row, not only printed: the category is documented as
            # one of the completeness CSV's, and a prefix that listed nothing is
            # the first thing to check when wells come back unaccounted for.
            sources.findings.append(
                finding_row(
                    "source_prefix_empty",
                    detail=(
                        f"untrimmed source {coverage.source_uri} yielded no vendor "
                        "CRAMs; check the prefix"
                    ),
                )
            )
    return sources


def _wafers_in_prefix(found: RawPrefix) -> tuple[set[str], set[str]]:
    """Wafers a located prefix holds, read both ways: (from folders, from names).

    Both, because the two layouts put the wafer in different places and a prefix
    may use either -- or, for a foldered delivery with wafer-level files beside
    its folders, both at once.
    """
    from_folders = {f for f in found.wafer_folders if SEAHUB_WAFER_RE.match(f)}
    from_names: set[str] = set()
    for key in found.loose_objects:
        name = key.split("/")[-1]
        if not name.endswith(_CRAM_SUFFIX) or _is_skippable_source_key(key):
            continue
        fields = parse_seahub_stem_fields(name[: -len(_CRAM_SUFFIX)])
        if fields is not None:
            from_names.add(str(fields["wafer"]))
    return from_folders, from_names


def discover_untrimmed_sources(
    s3_client: Any,
    roots: Any,
    seeds: WaferSeeds,
    expand_siblings: bool = True,
    max_depth: int = SEAHUB_SEARCH_MAX_DEPTH,
    max_listings: int = SEAHUB_SEARCH_MAX_LISTINGS,
) -> UntrimmedSources:
    """Index the vendor deliveries an upload's own wafers point at.

    Replaces naming the vendor orders. The identity of a delivery is taken from
    what is inside it -- a wafer, machine-assigned -- rather than from a path
    segment someone typed, which fails on every plausible vendor rename
    (``REF3-reupload``, ``reupload_REF3``, ``ref3``, a renumbered ``REF2``) and
    silently *succeeds* on some things it should not (``REF3_2`` matched
    ``REF3``; the real ``RNA3_098`` matches ``RNA3``).

    What replaces the ExperimentID filter as the guard against cross-experiment
    contamination is the seed set itself: a delivery is only indexed for the
    wafers the upload actually mentions, so walking a whole project cannot drag a
    foreign experiment's wells into this comparison. That is a positive assertion
    about contents rather than a string comparison, which is the point.

    ``expand_siblings`` indexes every wafer in a delivery that any seed pointed
    at, not just the seeded ones. It is on by default and it is what keeps the
    whole-wafer gap visible: a wafer the lab never trimmed cannot be a seed, so
    without expansion it would never be listed, never indexed, and produce
    neither a ``not_trimmed`` row nor a vendor-only ``DATA_GAP`` well -- and
    ``DATA_GAP`` is the only verdict the notebook writes to errors.txt. REF3
    delivers six wafers under one prefix, so one trimmed wafer drags in the other
    five. What expansion cannot reach is a delivery *no* uploaded wafer points at;
    those are reported as ``unseeded_vendor_delivery`` instead, which costs
    nothing because the descent enumerated them on its way past.
    """
    sources = UntrimmedSources()
    normalized_roots = normalize_search_roots(roots)
    plan = WaferSearchPlan(roots=tuple(normalized_roots), seeds=seeds)
    sources.search = plan
    if not normalized_roots or not seeds.wafers:
        return sources

    scan = _descend_to_raw_prefixes(
        s3_client, normalized_roots, max_depth=max_depth, max_listings=max_listings
    )
    plan.listings = scan.listings
    plan.unreadable = scan.unreadable
    plan.budget_exhausted = scan.budget_exhausted
    plan.depth_capped = scan.depth_capped

    wanted = set(seeds.wafers)
    located: dict[str, set[str]] = {}
    unseeded: list[tuple[str, str]] = []
    wanted_by_uri: dict[tuple[str, str], set[str]] = {}
    found_by_uri: dict[str, RawPrefix] = {}
    present_by_uri: dict[str, set[str]] = {}

    for found in sorted(scan.prefixes, key=lambda p: (p.bucket, p.prefix)):
        from_folders, from_names = _wafers_in_prefix(found)
        present = from_folders | from_names
        hits = present & wanted
        if not hits:
            # Present, enumerated, deliberately not indexed. Named so a delivery
            # nothing points at is a number rather than a silence.
            unseeded.append((found.uri, found.experiment_segment))
            continue
        for wafer in hits:
            located.setdefault(wafer, set()).add(found.uri)
        wanted_by_uri[(found.bucket, found.prefix)] = (
            present if expand_siblings else hits
        )
        found_by_uri[found.uri] = found
        present_by_uri[found.uri] = present

    if wanted_by_uri:
        # Delegated rather than reimplemented. Every prefix handed over is a
        # concrete `.../raw/`, four or five segments deep, so normalize_source_uris
        # accepts it -- and the "too broad" guard it still applies becomes an
        # invariant check on this walk: a bug in the descent cannot produce a
        # prefix shallow enough to flat-paginate a bucket without raising here.
        indexed = index_untrimmed_sources(
            s3_client,
            sorted(f"s3://{b}/{p}" for b, p in wanted_by_uri),
            accept=lambda bucket, key, fields: (
                str(fields["wafer"]) in _wanted_for(bucket, key, wanted_by_uri)
            ),
        )
        sources.index = indexed.index
        sources.coverage = indexed.coverage
        # source_prefix_empty says "check the prefix", which is the wrong advice
        # for a prefix nobody typed: the same fact under discovery is a located
        # delivery that yielded no wells, and wafer_folder_no_wells says so.
        sources.findings = [
            f for f in indexed.findings if f["category"] != "source_prefix_empty"
        ]
        for coverage in sources.coverage:
            found = found_by_uri.get(coverage.source_uri)
            if found is None:
                continue
            coverage.experiment_segment = found.experiment_segment
            indexed_wafers = {
                e.wafer
                for e in sources.index.values()
                if e.source_uri == coverage.source_uri
            }
            coverage.found_by = "sibling" if not indexed_wafers <= wanted else "seed"
            if coverage.cram_keys == 0:
                # Located by name, indexed nothing. The measured case is the
                # ScaleBio delivery whose 192 CRAMs carry no Z#### UG at all, so
                # every filename is refused by both stem patterns -- and the row
                # then reads cram_keys=0, exactly like a mistyped prefix. Present
                # and unreadable is a different problem from absent.
                hit = sorted(wanted & present_by_uri.get(coverage.source_uri, set()))
                coverage.skipped_reason = "no parseable wells"
                sources.findings.append(
                    finding_row(
                        "wafer_folder_no_wells",
                        wafer=", ".join(hit),
                        source_key=coverage.source_uri,
                        detail=(
                            f"{coverage.source_uri} was located for wafer(s) "
                            f"{', '.join(hit)} but yielded no vendor wells: no "
                            "object under it has a name that parses as "
                            "{wafer}-{sublibrary}_{type}-Z####-{barcode}.cram. "
                            "The delivery is present and unreadable, not absent"
                        ),
                    )
                )

    plan.located = {w: tuple(sorted(at)) for w, at in sorted(located.items()) if at}
    plan.unseeded = tuple(unseeded)
    _report_search(sources, plan)
    return sources


def _wanted_for(
    bucket: str, key: str, wanted_by_uri: dict[tuple[str, str], set[str]]
) -> set[str]:
    """The wafers a listed key's own delivery was located for.

    Keyed on ``(bucket, prefix)`` rather than on the prefix alone. Two search
    roots in different buckets can share a ``{project}/{order}/{exp}/raw/`` path
    -- a vendor delivery and a copy of it under another account -- and a
    prefix-only lookup would then test a key from one against the other's wanted
    set. Narrow, but the accept filter is the only thing keeping a foreign
    experiment out of this comparison, so it must not be approximate.
    """
    for (uri_bucket, prefix), wanted in wanted_by_uri.items():
        if uri_bucket == bucket and key.startswith(prefix):
            return wanted
    return set()


def _report_search(sources: UntrimmedSources, plan: WaferSearchPlan) -> None:
    """Turn the search accounting into findings rows.

    One row per fact, at the grain the fact has: a wafer that was not found is
    about that wafer, a delivery nothing pointed at is about that delivery, and
    an exhausted budget is about the whole search.
    """
    for wafer in plan.not_located:
        sources.findings.append(
            finding_row(
                "wafer_not_found",
                wafer=wafer,
                detail=(
                    f"no vendor delivery under {', '.join(plan.roots) or '(no root)'} "
                    f"holds wafer {wafer}; it cannot be reconciled, and if the "
                    "delivery exists elsewhere the search root does not cover it"
                ),
            )
        )
    for wafer in plan.collided:
        sources.findings.append(
            finding_row(
                "wafer_multiple_deliveries",
                wafer=wafer,
                detail=(
                    f"wafer {wafer} was found under {len(plan.located[wafer])} "
                    f"deliveries ({', '.join(plan.located[wafer])}); per-well "
                    "arbitration is reported as duplicate_source_well"
                ),
            )
        )
    for uri, segment in plan.unseeded:
        sources.findings.append(
            finding_row(
                "unseeded_vendor_delivery",
                detail=(
                    f"{uri} holds a vendor delivery"
                    + (f" for {segment!r}" if segment else "")
                    + " that no uploaded wafer points at, so it was not indexed; "
                    "if it belongs to this experiment its wells are missing from "
                    "the upload entirely"
                ),
            )
        )
    for token in plan.seeds.rejected:
        sources.findings.append(
            finding_row(
                "wafer_seed_rejected",
                detail=(
                    f"{token!r} came from the upload where a wafer was expected but "
                    "is not wafer-shaped, so it was not searched for"
                ),
            )
        )
    for uri, code in plan.unreadable:
        sources.findings.append(
            finding_row(
                "search_prefix_unreadable",
                source_key=uri,
                detail=f"listing {uri} was refused ({code}); the search skipped it",
            )
        )
    if plan.budget_exhausted:
        sources.findings.append(
            finding_row(
                "search_budget_exhausted",
                detail=(
                    f"the search stopped after {plan.listings} listings; it is "
                    "incomplete, so narrow the root to a project rather than a bucket"
                ),
            )
        )
    for uri in plan.depth_capped:
        sources.findings.append(
            finding_row(
                "search_depth_capped",
                source_key=uri,
                detail=(
                    f"stopped descending at {uri} without reaching a "
                    f"{SEAHUB_RAW_SEGMENT!r} folder; anything below it was not searched"
                ),
            )
        )


def load_source_read_counts(
    s3_client: Any,
    sources: UntrimmedSources | dict[IdentityKey, SourceEntry],
    bucket: str = "",
) -> None:
    """Populate ``read_count`` on each entry from its vendor sidecar.

    With several prefixes there is no single bucket, so each entry's own
    ``bucket`` is preferred and ``bucket`` is only a fallback for entries that
    predate it. An entry with neither is skipped rather than raising.

    One unreadable sidecar out of several hundred leaves that entry's
    ``read_count`` as None and lets the rest through. Raising instead took the
    whole notebook cell down over a single missing or malformed object, and the
    reconciliation already has a graceful answer for a well with no vendor count:
    it reports ``metadata_unavailable`` and carries on. A count that cannot be
    read and a count that was never there are the same fact to the report.
    """
    index = sources.index if isinstance(sources, UntrimmedSources) else sources
    targets = [
        (identity, entry.metadata_key, entry.bucket or bucket)
        for identity, entry in index.items()
        if entry.metadata_key and (entry.bucket or bucket)
    ]
    if not targets:
        return

    def _fetch(
        identity: IdentityKey, key: str, entry_bucket: str
    ) -> tuple[IdentityKey, int | None]:
        with tempfile.NamedTemporaryFile(
            mode="w+b", delete=False, suffix=".json"
        ) as tf:
            local = tf.name
        try:
            s3_client.download_file(entry_bucket, key, local)
            with open(local) as fh:
                payload = json.load(fh)
            value = payload.get("read_count")
            # int() on a non-numeric read_count fails the same way a missing
            # object does, and means the same thing: no usable count.
            return identity, int(value) if value is not None else None
        except Exception:  # noqa: BLE001 - any unreadable sidecar is "no count"
            return identity, None
        finally:
            Path(local).unlink(missing_ok=True)

    workers = min(SOURCE_METADATA_MAX_WORKERS, len(targets))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(_fetch, identity, key, entry_bucket)
            for identity, key, entry_bucket in targets
        ]
        for future in as_completed(futures):
            identity, read_count = future.result()
            index[identity].read_count = read_count


def source_order_by_wafer(
    index: dict[IdentityKey, SourceEntry],
) -> dict[str, str]:
    """Map each wafer to the order that delivered it.

    A wafer whose wells disagree gets every order joined by ``|`` -- deliberately
    ugly, because a wafer split across orders is worth noticing rather than
    silently picking one.
    """
    orders: dict[str, set[str]] = {}
    for entry in index.values():
        orders.setdefault(entry.wafer, set()).add(
            entry.source_order or SEAHUB_UNKNOWN_ORDER_LABEL
        )
    return {
        wafer: next(iter(seen)) if len(seen) == 1 else "|".join(sorted(seen))
        for wafer, seen in orders.items()
    }


def index_trimmed_upload(
    all_raw_files: list[str],
    sizes: dict[str, int] | None = None,
    findings: list[dict[str, Any]] | None = None,
) -> dict[IdentityKey, TrimmedEntry]:
    """Index a gathered SeaHub listing by ``(wafer, UG)``, both families."""
    index: dict[IdentityKey, TrimmedEntry] = {}
    duplicate_trimmed_seen: dict[IdentityKey, set[tuple[str, str]]] = {}
    for key in sorted(all_raw_files):
        parsed = seahub_stem_and_family(key.split("/")[-1])
        path_info = parse_seahub_raw_path(key)
        if parsed is None or path_info is None:
            continue
        stem, suffix, family = parsed
        fields = parse_seahub_stem_fields(stem)
        if fields is None:
            continue
        identity = (str(fields["wafer"]), str(fields["ug"]))
        raw_dir = "/".join(key.split("/")[:-1])
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
                raw_dir=raw_dir,
            )
            index[identity] = entry
        elif (path_info["sublibrary"], raw_dir, stem) != (
            entry.sublibrary,
            entry.raw_dir,
            entry.stem,
        ):
            # Not mere identity reuse -- that is normal, since the five artifacts
            # of a well all share it. This is one well appearing under two
            # different names or folders.
            if findings is not None:
                seen = duplicate_trimmed_seen.setdefault(identity, set())
                alt = (raw_dir, stem)
                if alt not in seen:
                    seen.add(alt)
                    findings.append(
                        finding_row(
                            "duplicate_trimmed_well",
                            wafer=entry.wafer,
                            ug=entry.ug,
                            sublibrary=entry.sublibrary,
                            trimmed_stem=entry.stem,
                            detail=(
                                f"well also present as {raw_dir}/{stem}; kept "
                                f"{entry.raw_dir}/{entry.stem}"
                            ),
                        )
                    )
        if suffix in (_CRAM_SUFFIX, ".trim.cram"):
            entry.has_cram = True
            if sizes:
                # Keep the largest when a well somehow carries several CRAMs.
                entry.size_bytes = max(entry.size_bytes, int(sizes.get(key, 0) or 0))
    return index


def seahub_wafer_seeds(
    trimmed_keys: list[str] | None = None,
    trimmed_index: dict[IdentityKey, TrimmedEntry] | None = None,
    discovered_wafers: set[str] | None = None,
) -> WaferSeeds:
    """Collect the wafers to look for on the vendor side, from the upload alone.

    Pure, and takes no S3 client: the three readings it unions are all already
    in hand by the time the notebook needs them, so locating the vendor
    deliveries costs no extra listing of the trimmed bucket.

    Why three readings rather than the identity index alone -- each rescues a
    case the others miss, verified against the real listings:

    * ``trimmed_index`` gives the wafer off the *filename*, so it is the only
      one that sees a well whose folder is not a wafer at all, or disagrees with
      its filename.  Both are real: ``wafer_mismatch`` is an SOP rule because
      folder and filename do diverge in practice.  Where they disagree, both
      tokens are searched for -- which of the two is right is not something this
      can decide, and the wrong one costs one lookup and reports as not found.
    * ``trimmed_keys`` gives the wafer off the *folder*, so it is the only one
      that sees a wafer whose filenames do not parse.  That is the ScaleBio
      delivery measured under ``RNA3_098``, whose 192 CRAMs carry no ``Z####``
      UG at all and therefore produce no identity, and it is any upload that
      misspells every artifact.
    * ``discovered_wafers`` comes from the folder walk, so it sees a wafer
      directory holding nothing this QA ingested -- an empty one included.  It
      is empty in manifest mode, which is why it cannot be the only reading.

    Neither path reading covers an object above or below wafer depth: the folder
    reading needs the SOP's six segments and the identity needs them too.  Such
    an object is reported as ``bad_path_depth`` and contributes no seed, which
    is deliberate -- there is no segment it could be read from.
    """
    tokens: list[str] = []
    if trimmed_index:
        tokens.extend(entry.wafer for entry in trimmed_index.values())
    for key in trimmed_keys or ():
        path_info = parse_seahub_raw_path(key)
        if path_info:
            tokens.append(path_info["wafer"])
    tokens.extend(discovered_wafers or ())

    wafers: set[str] = set()
    rejected: set[str] = set()
    for token in tokens:
        token = (token or "").strip()
        if not token:
            continue
        (wafers if SEAHUB_WAFER_RE.match(token) else rejected).add(token)
    return WaferSeeds(wafers=tuple(sorted(wafers)), rejected=tuple(sorted(rejected)))
