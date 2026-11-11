"""Extract Scale processed AnnData (.h5ad) metadata from a rundate prefix."""

from __future__ import annotations

import csv
import io
import json
import logging
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Callable, Sequence

from tqdm import tqdm

from .constants import MD5_COLUMN, SCALE_H5AD_COLUMNS
from .h5ad_introspect import count_h5ad_dims
from .metadata_sheet import fetch_sample_template, parse_sample_template
from .mtx_introspect import count_mtx_dims
from .models import NO_ERROR_MESSAGE, ListedObject, RunSummary, failure_detail
from .retry import retry_with_backoff
from .s3_utils import (
    fetch_crc64nvme,
    get_object_text,
    list_objects_with_size,
    parse_s3_uri,
    s3_uri_for,
    stream_md5,
)
from .scale_wells import (
    ScaleExtractError,
    correlate_samples,
    expand_barcodes,
    parse_well,
)
from .sheets import LabIdentity
from .tsv_writer import TsvWriter

log = logging.getLogger(__name__)

SAMPLES_CSV_NAME = "samples.csv"
SAMPLES_DIR = "samples"
SCALEPLEX_DIR = "scaleplex"
H5AD_SUFFIX = ".h5ad"
MERGED_ANNDATA_SUFFIX = "_anndata.h5ad"
MTX_BASENAMES = frozenset({"matrix.mtx.gz", "matrix.mtx"})
SCALEPLEX_MATRIX_DIR_RE = re.compile(r"\.QSR-\d+-SCALEPLEX\.filtered\.matrix$")
H5AD_FEATURE_TYPE = "gene"
SCALEPLEX_FEATURE_TYPE = "hash oligo"
# Thread count, and the S3 connection-pool size the CLI gives those threads.
DEFAULT_WORKERS = 16
CRAM_QSR_WELL_RE = re.compile(
    r"QSR-(\d+)(?:-SCALEPLEX)?[_-](\d{1,2}[A-H])\.cram$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ScaleCram:
    qsr: str
    well: str
    scaleplex: bool


@dataclass(frozen=True)
class FileResult:
    """What one listed file yielded; an error string is empty on success."""

    crc: str | None
    crc_error: str
    md5: str | None
    md5_error: str
    observation_count: int | None
    feature_counts: str | None
    introspect_error: str


@dataclass(frozen=True)
class RawCramSearch:
    """CRAMs found under ``--raw-subdirs``, and the entries that found none.

    ``empty_subdirs`` holds ``(subdir, searched_uris)`` so a mistyped or
    wrong-layout entry can be named in a warning -- alongside every
    prefix tried for it -- instead of quietly leaving ``derived_from``
    empty.
    """

    crams: list[tuple[str, str]]
    empty_subdirs: list[tuple[str, tuple[str, ...]]]


def format_samples_column(sample_names: Sequence[str], lab: str) -> str:
    """JSON list of sheet sample_name values prefixed with ``{lab}:``."""
    prefix = LabIdentity.parse(lab).name
    return json.dumps([f"{prefix}:{name}" for name in sample_names])


def format_feature_counts(feature_type: str, count: int) -> str:
    """Same list-of-dicts JSON shape as 10x h5 ``feature_counts``."""
    return json.dumps([{"feature_type": feature_type, "feature_count": count}])


def format_derived_from(labels: Sequence[str]) -> str:
    """JSON list of ``{lab}:{cram_filename}`` values attached to one file."""
    return json.dumps(list(labels))


def derived_from_label(lab: str, key: str) -> str:
    """``{lab}:{basename}`` for a CRAM object key."""
    prefix = LabIdentity.parse(lab).name
    return f"{prefix}:{key.rsplit('/', 1)[-1]}"


def is_deliverable_cram_key(key: str) -> bool:
    """True for ``*.cram`` objects that are not unmatched leftovers."""
    name = key.rsplit("/", 1)[-1]
    return name.lower().endswith(".cram") and "unmatched" not in name.lower()


def parse_scale_cram_name(filename: str) -> ScaleCram | None:
    """Read QSR#, well, and ScalePlex from a Scale CRAM basename."""
    name = filename.rsplit("/", 1)[-1]
    if not is_deliverable_cram_key(name):
        return None
    match = CRAM_QSR_WELL_RE.search(name)
    if not match:
        return None
    try:
        well = parse_well(match.group(2))
    except ScaleExtractError:
        return None
    return ScaleCram(
        qsr=str(int(match.group(1))),
        well=well,
        scaleplex="SCALEPLEX" in name.upper(),
    )


def derived_from_filenames(sample: str, parsed: ScaleCram) -> tuple[str, ...]:
    """TSV ``filename`` values a CRAM can attach to."""
    if parsed.scaleplex:
        dirname = f"{sample}.QSR-{parsed.qsr}-SCALEPLEX.filtered.matrix"
        return (f"{dirname}/matrix.mtx.gz", f"{dirname}/matrix.mtx")
    return (f"{sample}.QSR-{parsed.qsr}_anndata.h5ad",)


def well_to_sample_map(sample_rows: Sequence[tuple[str, str]]) -> dict[str, str]:
    """Map each expanded barcodes well to its samples.csv sample."""
    owner: dict[str, str] = {}
    for sample, barcodes in sample_rows:
        for well in expand_barcodes(barcodes):
            owner[well] = sample
    return owner


def resolve_raw_subdir(
    processed_bucket: str, processed_prefix: str, subdir: str
) -> tuple[str, str]:
    """Turn a ``--raw-subdirs`` value into ``(bucket, prefix/)``.

    A name such as ``426971`` is the folder under the ``raw/`` sibling of
    the processed rundate. An ``s3://`` URI is used as given (typically the
    group directory that contains ``raw/{numeric}/``).
    """
    item = (subdir or "").strip().rstrip("/")
    if not item:
        raise ScaleExtractError("--raw-subdirs must not contain empty values")
    if item.startswith("s3://"):
        location = parse_s3_uri(item + "/")
        return location.bucket, location.prefix
    parts = processed_prefix.rstrip("/").split("/")
    try:
        idx = parts.index("processed")
    except ValueError as exc:
        raise ScaleExtractError(
            f"processed prefix {processed_prefix!r} does not contain a "
            "processed/ segment; pass an s3:// URI in --raw-subdirs"
        ) from exc
    parent = "/".join(parts[:idx])
    prefix = f"{parent}/raw/{item}/" if parent else f"raw/{item}/"
    return processed_bucket, prefix


def raw_cram_search_prefixes(prefix: str) -> tuple[str, ...]:
    """Prefixes to try for the CRAMs under one ``--raw-subdirs`` value.

    A ``--raw-subdirs`` value is either the numeric run folder holding
    ``*.cram`` directly or a group directory whose ``raw/`` child holds
    the numeric folders. No segment name tells the two apart -- a group
    directory can itself sit under a ``raw/`` level -- so try the ``raw/``
    child first and fall back to the prefix itself. Whichever lists CRAMs
    wins, which makes the layout observed rather than inferred.

    The ``raw/`` child comes first because listing a group directory
    walks its whole subtree only to reject every CRAM under
    ``raw/{numeric}/``; the reverse order wastes one empty listing.

    A delivery holds its CRAMs at one of these levels, never both, so
    the first prefix that lists any CRAM has all of them. Searching the
    remaining prefixes as well would only re-walk a group subtree.

    Only a prefix whose final *segment* is ``raw`` is already at that
    level. A name that merely ends in ``raw``, such as ``ORD01_raw``,
    still needs its ``raw/`` child probed.

    An empty prefix is a bucket root. Probing ``raw/`` there would still
    accept ``raw/{numeric}/`` from any delivery in the bucket, so there is
    no candidate that both honours the request and stays inside one
    delivery -- hence the error. ``validate_raw_subdirs`` rejects a
    bucket-alone value before the CLI reaches here; this keeps a direct
    caller from searching a whole bucket instead.
    """
    trimmed = prefix.strip("/")
    if not trimmed:
        raise ScaleExtractError(
            "raw cram search needs a directory, not a bucket root: a "
            "bucket-wide search could match crams from another delivery"
        )
    normalized = f"{trimmed}/"
    if trimmed.rsplit("/", 1)[-1] == "raw":
        return (normalized,)
    return (f"{normalized}raw/", normalized)


def is_cram_in_raw_search(key: str, search_prefix: str) -> bool:
    """True for a deliverable CRAM in a numeric run folder.

    That folder is either a child of ``search_prefix`` or, when the
    prefix's own final segment is numeric, the prefix itself. The prefix
    need not contain a ``raw`` segment: since the search falls back to
    the bare ``--raw-subdirs`` prefix, it often does not.
    """
    if not is_deliverable_cram_key(key) or not key.startswith(search_prefix):
        return False
    rest = key[len(search_prefix) :]
    folder, sep, _name = rest.partition("/")
    if not sep:
        return search_prefix.rstrip("/").rsplit("/", 1)[-1].isdigit()
    return folder.isdigit()


def list_raw_crams(
    s3_client: Any,
    processed_bucket: str,
    processed_prefix: str,
    raw_subdirs: Sequence[str],
) -> RawCramSearch:
    """Find the CRAMs each ``--raw-subdirs`` value points at.

    A CRAM is read from a numeric run folder, which is either the prefix
    searched or a folder directly under it. Which prefix won is logged at
    debug level, so an unexpected layout is visible under ``-v`` even
    though the search itself succeeded.
    """
    found: list[tuple[str, str]] = []
    empty: list[tuple[str, tuple[str, ...]]] = []
    seen: set[tuple[str, str]] = set()
    for subdir in raw_subdirs:
        bucket, prefix = resolve_raw_subdir(processed_bucket, processed_prefix, subdir)
        candidates = raw_cram_search_prefixes(prefix)
        matched = 0
        for search in candidates:

            def _keep(key: str, _search: str = search) -> bool:
                return is_cram_in_raw_search(key, _search)

            for obj in list_objects_with_size(
                s3_client, bucket, search, predicate=_keep
            ):
                matched += 1
                ident = (bucket, obj.key)
                if ident in seen:
                    continue
                seen.add(ident)
                found.append(ident)
            if matched:
                log.debug(
                    "--raw-subdirs %r: %d crams under %s",
                    subdir,
                    matched,
                    s3_uri_for(bucket, search),
                )
                break
        if not matched:
            empty.append(
                (subdir, tuple(s3_uri_for(bucket, cand) for cand in candidates))
            )
    return RawCramSearch(crams=found, empty_subdirs=empty)


def derived_from_by_filename(
    crams: Sequence[tuple[str, str]],
    well_to_sample: dict[str, str],
    control_samples: set[str] | frozenset[str],
    lab: str,
) -> dict[str, list[str]]:
    """Group ``{lab}:{cram_filename}`` labels by the processed TSV filename."""
    grouped: dict[str, list[str]] = defaultdict(list)
    for _bucket, key in crams:
        parsed = parse_scale_cram_name(key)
        if parsed is None:
            continue
        sample = well_to_sample.get(parsed.well)
        if sample is None or sample in control_samples:
            continue
        label = derived_from_label(lab, key)
        for filename in derived_from_filenames(sample, parsed):
            if label not in grouped[filename]:
                grouped[filename].append(label)
    for filename in grouped:
        grouped[filename].sort()
    return grouped


def leftover_cram_uris(
    crams: Sequence[tuple[str, str]],
    derived_map: dict[str, list[str]],
    target_filenames: set[str] | frozenset[str],
    lab: str,
) -> list[str]:
    """S3 URIs of listed CRAMs that did not attach to an output row."""
    attached: set[str] = set()
    for filename in target_filenames:
        attached.update(derived_map.get(filename, ()))
    leftover: list[str] = []
    for bucket, key in crams:
        if derived_from_label(lab, key) not in attached:
            leftover.append(s3_uri_for(bucket, key))
    return leftover


def leftover_cram_warning(uri: str) -> str:
    return f"cram {uri} was not added to derived_from"


def empty_derived_from_warning(filename: str) -> str:
    return (
        f"derived_from is empty for {filename!r}; "
        "are any S3 raw cram directories missing from --raw-subdirs?"
    )


def empty_raw_subdir_warning(subdir: str, uris: Sequence[str]) -> str:
    """One ``--raw-subdirs`` entry found nothing to fill ``derived_from``.

    The advice has to hold for every prefix printed: a flat ``*.cram`` is
    read only when the prefix is itself the numeric run folder, so a
    ``raw/`` prefix never accepts one.
    """
    return (
        f"--raw-subdirs {subdir!r} matched no crams under {' or '.join(uris)}; "
        "*.cram must sit in a numeric run folder -- one of those prefixes "
        "itself, or a folder directly under one"
    )


def default_scale_h5ad_output_name(prefix: str) -> str:
    rundate = prefix.rstrip("/").rsplit("/", 1)[-1] or "output"
    return f"{rundate}_scale_h5ad_info.tsv"


def samples_csv_key(prefix: str) -> str:
    return prefix.rstrip("/") + f"/{SAMPLES_CSV_NAME}"


def samples_dir_prefix(prefix: str) -> str:
    return prefix.rstrip("/") + f"/{SAMPLES_DIR}/"


def scaleplex_dir_prefix(prefix: str) -> str:
    return prefix.rstrip("/") + f"/{SCALEPLEX_DIR}/"


def is_scale_h5ad_key(key: str, rundate_prefix: str) -> bool:
    """True for .h5ad files under {rundate}/samples/ whose name contains QSR."""
    base = samples_dir_prefix(rundate_prefix)
    if not key.startswith(base):
        return False
    name = key[len(base) :]
    return (
        bool(name) and "/" not in name and name.endswith(H5AD_SUFFIX) and "QSR" in name
    )


def is_scale_mtx_key(key: str, rundate_prefix: str) -> bool:
    """True for matrix.mtx(.gz) under a ScalePlex filtered.matrix directory."""
    base = scaleplex_dir_prefix(rundate_prefix)
    if not key.startswith(base):
        return False
    rest = key[len(base) :]
    parts = rest.split("/")
    if len(parts) != 2:
        return False
    dirname, basename = parts
    return basename in MTX_BASENAMES and bool(SCALEPLEX_MATRIX_DIR_RE.search(dirname))


def tsv_filename(key: str, rundate_prefix: str) -> str:
    """Basename for h5ad; directory-relative path for ScalePlex mtx."""
    if is_scale_mtx_key(key, rundate_prefix):
        return key[len(scaleplex_dir_prefix(rundate_prefix)) :]
    return key.rsplit("/", 1)[-1]


def sample_from_filename(filename: str) -> str:
    return filename.split(".", 1)[0]


def file_belongs_to_sample(filename: str, sample: str) -> bool:
    return filename == f"{sample}{MERGED_ANNDATA_SUFFIX}" or filename.startswith(
        f"{sample}."
    )


def is_control_file(filename: str, control_samples: set[str] | frozenset[str]) -> bool:
    return any(file_belongs_to_sample(filename, sample) for sample in control_samples)


def parse_samples_csv(text: str) -> list[tuple[str, str]]:
    reader = csv.DictReader(io.StringIO(text))
    fields = reader.fieldnames or []
    if "sample" not in fields or "barcodes" not in fields:
        raise ScaleExtractError("samples.csv must have sample and barcodes columns")
    rows: list[tuple[str, str]] = []
    for row in reader:
        sample = (row.get("sample") or "").strip()
        barcodes = (row.get("barcodes") or "").strip()
        if not sample:
            continue
        rows.append((sample, barcodes))
    return rows


def control_warning(sample: str, barcodes: str) -> str:
    return (
        f"control sample {sample!r} (barcodes {barcodes!r}) has no RT_index "
        "pairing; excluding from output"
    )


def _introspect_counts(s3_client: Any, bucket: str, key: str) -> tuple[int, str]:
    """n_obs and feature_counts JSON from one open of each listed file."""
    if key.endswith(H5AD_SUFFIX):
        n_obs, n_vars = count_h5ad_dims(bucket, key)
        return n_obs, format_feature_counts(H5AD_FEATURE_TYPE, n_vars)
    if key.rsplit("/", 1)[-1] in MTX_BASENAMES:
        n_obs, n_features = count_mtx_dims(s3_client, bucket, key)
        return n_obs, format_feature_counts(SCALEPLEX_FEATURE_TYPE, n_features)
    raise RuntimeError(f"no observation counter for {key}")


def _process_one(
    s3_client: Any, bucket: str, key: str, *, do_md5: bool, retries: int
) -> FileResult:
    crc, crc_err = retry_with_backoff(
        fetch_crc64nvme, s3_client, bucket, key, retries=retries
    )
    hashed, md5_err = None, None
    if do_md5:
        hashed, md5_err = retry_with_backoff(
            stream_md5, s3_client, bucket, key, retries=retries
        )
    intro, intro_err = retry_with_backoff(
        _introspect_counts, s3_client, bucket, key, retries=retries
    )
    obs, feature_counts = intro if intro is not None else (None, None)
    md5_error = ""
    if do_md5 and hashed is None:
        md5_error = md5_err or NO_ERROR_MESSAGE
    return FileResult(
        crc=crc,
        crc_error=crc_err or "",
        # A failed read leaves the cell empty, never a partial digest.
        md5=hashed[0] if hashed is not None else None,
        md5_error=md5_error,
        observation_count=obs,
        feature_counts=feature_counts,
        introspect_error=intro_err or "",
    )


def scale_h5ad_columns(*, do_md5: bool) -> list[str]:
    return [col for col in SCALE_H5AD_COLUMNS if do_md5 or col != MD5_COLUMN]


def extract_scale_h5ad(
    s3_client: Any,
    bucket: str,
    prefix: str,
    output_path: str,
    *,
    metadata_gid: str,
    metadata_experiment: str,
    lab: str,
    raw_subdirs: Sequence[str],
    do_md5: bool = True,
    workers: int | None = None,
    retries: int = 5,
    show_progress: bool = True,
    sheet_csv: str | None = None,
    fetch_sheet: Callable[[str], str] | None = None,
) -> RunSummary:
    """List non-control samples/*.h5ad and scaleplex mtx files, write the TSV.

    ``metadata_experiment`` selects sample template rows by
    ``experiment_name``. ``raw_subdirs`` are walked for ``*.cram`` files
    that fill ``derived_from``.
    """
    lab_name = LabIdentity.parse(lab).name

    csv_key = samples_csv_key(prefix)
    try:
        samples_text = get_object_text(s3_client, bucket, csv_key)
    except Exception as exc:
        raise ScaleExtractError(
            f"could not read {s3_uri_for(bucket, csv_key)}: {exc}"
        ) from exc
    sample_rows = parse_samples_csv(samples_text)

    if sheet_csv is None:
        loader = fetch_sheet or fetch_sample_template
        sheet_csv = loader(metadata_gid)
    template_rows = parse_sample_template(sheet_csv, experiment=metadata_experiment)
    sheet_wells = {row["well"] for row in template_rows}
    sheet_names = [(row["well"], row["sample_name"]) for row in template_rows]
    correlation = correlate_samples(sample_rows, sheet_wells, sheet_names)
    control_names = set(correlation.control_set)
    raw_search = list_raw_crams(s3_client, bucket, prefix, raw_subdirs)
    crams = raw_search.crams
    derived_map = derived_from_by_filename(
        crams,
        well_to_sample_map(sample_rows),
        control_names,
        lab_name,
    )

    summary = RunSummary()
    for control in correlation.controls:
        summary.warnings.append(control_warning(control.sample, control.barcodes))
    for subdir, searched_uris in raw_search.empty_subdirs:
        summary.empty_raw_subdirs.append(subdir)
        summary.warnings.append(empty_raw_subdir_warning(subdir, searched_uris))

    listed = list_objects_with_size(
        s3_client,
        bucket,
        samples_dir_prefix(prefix),
        predicate=lambda key: is_scale_h5ad_key(key, prefix),
    )
    listed.extend(
        list_objects_with_size(
            s3_client,
            bucket,
            scaleplex_dir_prefix(prefix),
            predicate=lambda key: is_scale_mtx_key(key, prefix),
        )
    )
    targets: list[ListedObject] = []
    for obj in listed:
        filename = tsv_filename(obj.key, prefix)
        if is_control_file(filename, control_names):
            continue
        targets.append(obj)

    target_names = {tsv_filename(obj.key, prefix) for obj in targets}
    for uri in leftover_cram_uris(crams, derived_map, target_names, lab_name):
        summary.warnings.append(leftover_cram_warning(uri))
    for filename in sorted(target_names):
        if not derived_map.get(filename):
            summary.warnings.append(empty_derived_from_warning(filename))

    summary.total = len(targets)
    if not targets:
        return summary

    writer = TsvWriter(output_path, scale_h5ad_columns(do_md5=do_md5))
    max_workers = min(workers or DEFAULT_WORKERS, len(targets))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_one,
                s3_client,
                bucket,
                obj.key,
                do_md5=do_md5,
                retries=retries,
            ): obj
            for obj in targets
        }
        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(targets), desc="Processing")

        # Keep TSV rows in listing order, not completion order.
        results: dict[str, FileResult] = {}
        try:
            for fut in iterator:
                obj = futures[fut]
                results[obj.key] = fut.result()
        except BaseException:
            # Ctrl-C must not drain the queue: the executor's __exit__ would
            # run every queued file, md5 streams included, to completion. Only
            # the files already in flight finish; a thread cannot be stopped.
            executor.shutdown(wait=False, cancel_futures=True)
            raise

    for obj in targets:
        filename = tsv_filename(obj.key, prefix)
        r = results[obj.key]
        sample = sample_from_filename(filename)
        row: list[object] = [
            filename,
            s3_uri_for(bucket, obj.key),
            r.crc if r.crc is not None else "",
        ]
        if do_md5:
            row.append(r.md5 if r.md5 is not None else "")
        row.extend(
            [
                sample,
                format_samples_column(
                    correlation.sample_names.get(sample, ()), lab_name
                ),
                obj.size_bytes,
                r.observation_count if r.observation_count is not None else "",
                r.feature_counts if r.feature_counts is not None else "",
                format_derived_from(derived_map.get(filename, ())),
            ]
        )
        writer.append_row(row)
        if not r.crc_error:
            summary.crc_ok += 1
        if do_md5 and not r.md5_error:
            summary.md5_ok += 1
        if not r.introspect_error:
            summary.enrichment_ok += 1
        if r.crc_error or r.introspect_error or r.md5_error:
            summary.failures.append(
                (
                    obj.key,
                    r.crc_error,
                    failure_detail(r.introspect_error, r.md5_error),
                )
            )

    return summary
