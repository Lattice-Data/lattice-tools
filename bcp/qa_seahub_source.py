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
    SEAHUB_STEM_NO_TYPE_RE,
    SEAHUB_STEM_RE,
    SEAHUB_UNKNOWN_ORDER_LABEL,
    SEAHUB_VENDOR_ORDER_RE,
    SEAHUB_WAFER_RE,
)
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family

__all__ = [
    "IdentityKey",
    "WaferSeeds",
    "SourceCoverage",
    "SourceEntry",
    "TrimmedEntry",
    "UntrimmedSources",
    "derive_source_experiment",
    "derive_source_order",
    "finding_row",
    "index_trimmed_upload",
    "index_untrimmed_sources",
    "load_source_read_counts",
    "normalize_source_uris",
    "parse_seahub_stem_fields",
    "parse_source_uri",
    "source_experiment_matches",
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

    @property
    def unmatched(self) -> int:
        return max(0, self.indexed - self.matched)


@dataclass
class UntrimmedSources:
    """The merged vendor index plus per-prefix coverage and findings."""

    index: dict[IdentityKey, SourceEntry] = field(default_factory=dict)
    coverage: list[SourceCoverage] = field(default_factory=list)
    findings: list[dict[str, Any]] = field(default_factory=list)

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


def derive_source_experiment(cram_key: str) -> str:
    """Read the experiment segment out of a key, positionally.

    ``{project}/{order}/{ExperimentID}/raw/{wafer}/{file}`` -- the segment before
    ``raw``, and the one an order-level prefix spans several of. Returns ``""``
    for any other shape; callers keep those entries rather than dropping a well
    they cannot place. Compare it with :func:`source_experiment_matches` rather
    than by equality, since the segment is not always the ExperimentID alone.
    """
    parts = cram_key.split("/")
    if len(parts) == 6 and parts[3] == "raw":
        return parts[2]
    return ""


def source_experiment_matches(segment: str, experiment_id: str) -> bool:
    """Does a pre-``raw`` segment belong to ``experiment_id``?

    The current vendor layout carries the ExperimentID alone (``REF3``); an older
    one appends the sublibrary (``REF5_P01``), which
    :func:`derive_source_experiment` above reads off the key, and a
    re-delivery appends something else again (``GENE7_reupload``, measured in
    order ``NVUS2024101701-11`` alongside ``REF3``). So the rule is any
    ``{ExperimentID}_...`` folder, not a sublibrary shape: a bare equality test
    would exclude every well of such a delivery and report the whole upload as
    orphaned. The underscore is required so ``REF50`` is not read as ``REF5``.
    """
    return segment == experiment_id or segment.startswith(f"{experiment_id}_")


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


def index_untrimmed_sources(
    s3_client: Any, uris: Any, experiment_id: str = ""
) -> UntrimmedSources:
    """List one or more vendor prefixes and index per-well CRAMs by ``(wafer, UG)``.

    An identity delivered by two prefixes is a re-delivery, not a merge: the
    winner is deterministic (lowest ``(source_order, cram_key)``) and the loser
    becomes a ``duplicate_source_well`` finding. Overwriting blindly, as the
    single-prefix version did, silently discarded one of them.

    ``experiment_id`` is the ExperimentID being QA'd. A prefix is documented as
    order-level, and one order holds several experiments, so without it the index
    also carries the *other* experiments' wells -- each of which then has no
    counterpart in the upload and surfaces as ``not_trimmed`` plus an ``UNKNOWN``
    well. Wells whose ExperimentID reads as a different one are excluded and
    reported once per foreign experiment. A key whose shape hides the
    ExperimentID is kept, since a spurious ``not_trimmed`` row is recoverable
    and a silently dropped vendor well is not.
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
    # Keyed on the stem so a sidecar can only ever attach to its own CRAM,
    # independent of listing order.
    sidecar_by_stem_key: dict[str, str] = {}

    for uri, bucket, prefix in survivors:
        coverage = SourceCoverage(source_uri=uri, bucket=bucket, prefix=prefix)
        orders_seen: set[str] = set()
        foreign_wells: dict[str, int] = {}
        unplaced_wells = 0

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
                    sidecar_by_stem_key[key[: -len(_SIDECAR_SUFFIX)]] = key
                    continue

                if experiment_id:
                    experiment = derive_source_experiment(key)
                    if not experiment:
                        unplaced_wells += 1
                    elif not source_experiment_matches(experiment, experiment_id):
                        foreign_wells[experiment] = foreign_wells.get(experiment, 0) + 1
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

                winner, loser = sorted(
                    (incumbent, candidate),
                    key=lambda e: (e.source_order, e.cram_key),
                )
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

        for experiment, count in sorted(foreign_wells.items()):
            sources.findings.append(
                finding_row(
                    "source_prefix_spans_experiments",
                    detail=(
                        f"untrimmed source {uri} also holds {count} well(s) of "
                        f"experiment {experiment!r}, which were excluded from the "
                        f"{experiment_id!r} comparison; narrow the prefix to "
                        f"{uri.rstrip('/')}/{experiment_id} to stop listing them"
                    ),
                )
            )
        if unplaced_wells:
            sources.findings.append(
                finding_row(
                    "source_experiment_unreadable",
                    detail=(
                        f"untrimmed source {uri} holds {unplaced_wells} CRAM(s) whose "
                        "key does not expose an ExperimentID (expected "
                        "{project}/{order}/{ExperimentID}/raw/{wafer}/{file}); kept "
                        f"in the {experiment_id!r} comparison rather than dropped"
                    ),
                )
            )

        coverage.orders_seen = tuple(sorted(orders_seen))
        coverage.source_order = (
            coverage.orders_seen[0]
            if len(coverage.orders_seen) == 1
            else _order_label_from_prefix(prefix)
        )
        sources.coverage.append(coverage)

    for identity, entry in sources.index.items():
        sidecar = sidecar_by_stem_key.get(entry.cram_key[: -len(_CRAM_SUFFIX)])
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


def _wafer_seeds(
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
