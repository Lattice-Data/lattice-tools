"""Check a guide metadata sheet's guide_ids against Cell Ranger h5 matrices.

Every ``guide_id`` in the sheet must appear among the ``CRISPR Guide Capture``
feature IDs of the h5 files found under an S3 order prefix. Guides present in
an h5 but absent from the sheet are reported as notes, not failures.

Two scopes decide which h5 files those are:

``multi`` (the default)
    The run-level ``raw_feature_bc_matrix.h5``, one per library, found under
    either ``outs/`` (Cell Ranger 10.x flat layout) or ``outs/multi/count/``
    (9.0.1 nested layout). Within a single ``cellranger multi`` run the feature
    reference is per-run, not per-sample, so every per-sample h5 carries an
    identical features table; reading the run-level file once per library is
    the whole guide_id check. A library with no such h5 fails, and a library
    with more than one fails as an ambiguity rather than being read.

``per-sample``
    Every ``per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix.h5``.
    Because the features table is identical across them, this is a file
    integrity sweep over an order's per-sample outputs, not a guide check: it
    catches missing, truncated or unreadable per-sample files. The sample
    folder checks and ``samples``/``read_sample_names`` filtering belong to
    this scope only.

S3 listing, retry and sample/library naming are reused from ``file_extract`` so
the two tools agree on which files an order holds and what they are called.
This module is imported only for ``--h5-uri`` runs, so a plain signature run
needs none of boto3, h5py, fsspec or s3fs.
"""

from __future__ import annotations

from collections.abc import Collection
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from file_extract.constants import DEFAULT_H5_TARGET_FILENAME
from file_extract.h5 import extract_library, extract_sample_name, h5_worker_ceiling
from file_extract.retry import retry_with_backoff
from file_extract.s3_utils import list_objects_with_size, parse_s3_uri, s3_uri_for

from .core import GuideSigError

SCOPE_MULTI = "multi"
SCOPE_PER_SAMPLE = "per-sample"
SCOPE_BOTH = "both"
SCOPES = (SCOPE_MULTI, SCOPE_PER_SAMPLE, SCOPE_BOTH)

#: Run-level matrix written by ``cellranger multi``, under ``outs/multi/count/``.
DEFAULT_MULTI_H5_FILENAME = "raw_feature_bc_matrix.h5"
#: Per-sample matrix, under ``outs/per_sample_outs/<sample>/count/``.
DEFAULT_PER_SAMPLE_H5_FILENAME = DEFAULT_H5_TARGET_FILENAME

# Where a run-level matrix sits depends on the Cell Ranger layout, both of
# which this repo supports (see qa_constants.CELLRANGER_FLAT_LAYOUT):
# 9.0.1 nests it under outs/multi/count/, while 10.0.0 and 10.1.0 write it
# flat at the top of outs/ with no multi/count/ subtree at all. The basename
# is the same either way, so both are accepted and the per-library cardinality
# rule decides the rest.
MULTI_COUNT_SEGMENT = "outs/multi/count/"
FLAT_OUTS_SEGMENT = "outs/"
RUN_LEVEL_LOCATIONS = f"{FLAT_OUTS_SEGMENT} and {MULTI_COUNT_SEGMENT}"
CRISPR_FEATURE_TYPE = "CRISPR Guide Capture"
EXAMPLE_LIMIT = 10

# Default thread/pool sizes, by scope; see default_workers_for. These are
# guidesig's own, so file_extract's shared ceiling keeps its tuned value for
# the high-fanout introspection it still does.
MULTI_SCOPE_WORKERS = 24
PER_SAMPLE_SCOPE_WORKERS = 8

# Bytes s3fs fetches per range request. Reading the features table takes only
# about eight seeks, but s3fs defaults to 50 MiB blocks, so on a multi-GB
# run-level matrix each of those seeks dragged down a whole 50 MiB block.
# Measured against one real 4.3 GiB raw_feature_bc_matrix.h5, reading the same
# features table at each block size:
#
#     block size   requests   fetched   elapsed
#       50 MiB            7   226.5 MiB   25.04s   <- s3fs default
#        8 MiB            8    56.0 MiB    2.77s
#        1 MiB            8     7.0 MiB    1.74s
#      256 KiB           10     2.3 MiB    1.95s
#
# Request count barely moves, so the cost was over-fetching rather than round
# trips. 256 KiB fetches less but is slower, being where the extra requests
# start to show; 1 MiB is the floor worth using.
FEATURES_BLOCK_SIZE = 1 * 2**20

# reader(bucket, key) -> CRISPR feature IDs in that h5.
H5Reader = Callable[[str, str], "set[str]"]


class H5CheckError(GuideSigError):
    """Raised when the h5 check cannot run at all (bad URI, listing failed)."""


def default_filename_for(scope: str) -> str:
    """The h5 basename a scope looks for when none is given."""
    if scope == SCOPE_MULTI:
        return DEFAULT_MULTI_H5_FILENAME
    return DEFAULT_PER_SAMPLE_H5_FILENAME


def default_workers_for(scope: str) -> int:
    """Threads, and S3 pool size, a scope uses when ``--workers`` is unset.

    Opening an h5 costs many range reads, so the ceiling is about how many
    files a scope opens at once rather than how fast one of them is. ``multi``
    opens one per library, which is few enough to widen; the per-sample sweep
    opens one per sample, which is the high-fanout case the shared
    ``file_extract`` ceiling of 8 was tuned for. ``both`` runs that sweep too,
    so it keeps the lower number.
    """
    if scope == SCOPE_MULTI:
        return MULTI_SCOPE_WORKERS
    return PER_SAMPLE_SCOPE_WORKERS


@dataclass(frozen=True)
class MultiLayout:
    """What one listing says about an order's run-level ``multi/count`` h5s.

    One library must hold exactly one. ``h5_keys`` are the libraries that do;
    ``libraries_without_h5`` and ``ambiguous_libraries`` are the two ways a
    library can fail instead.
    """

    h5_keys: list[str]
    libraries_without_h5: list[str]
    ambiguous_libraries: list[tuple[str, list[str]]]


@dataclass(frozen=True)
class PerSampleLayout:
    """What one listing says about an order's ``per_sample_outs`` h5s."""

    h5_keys: list[str]
    sample_dirs_without_h5: list[str]
    library_dirs_without_per_sample_outs: list[str]


@dataclass
class H5CheckReport:
    """Outcome of checking one sheet's guide_ids against an order's h5 files."""

    h5_total: int = 0
    h5_failed: int = 0
    failures: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failures


def check_dependencies() -> None:
    """Fail with an actionable message if an S3/h5 package is missing."""
    try:
        import boto3  # noqa: F401
        import fsspec  # noqa: F401
        import h5py  # noqa: F401
        import s3fs  # noqa: F401
    except ImportError as exc:
        raise H5CheckError(
            "--h5-uri needs boto3, h5py, fsspec and s3fs "
            f"(missing: {exc.name}); install them from requirements.txt"
        ) from exc


def _decode(values: Any) -> list[str]:
    return [
        v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v)
        for v in values
    ]


def guide_ids_from_h5(h5: Any) -> set[str]:
    """Return the CRISPR Guide Capture feature IDs from an open h5py File."""
    if "matrix" not in h5:
        raise RuntimeError(
            "no 'matrix' group; not a Cell Ranger feature-barcode matrix"
        )
    matrix = h5["matrix"]
    for dataset in ("features/id", "features/feature_type"):
        if dataset not in matrix:
            raise RuntimeError(f"no 'matrix/{dataset}' dataset")
    ids = _decode(matrix["features/id"][:])
    types = _decode(matrix["features/feature_type"][:])
    return {fid for fid, ftype in zip(ids, types) if ftype == CRISPR_FEATURE_TYPE}


def read_s3_h5_guide_ids(bucket: str, key: str, *, pool_size: int) -> set[str]:
    """Open an h5 over S3 and return its CRISPR feature IDs.

    Only the features table is read, never the count matrix, so the block size
    is set explicitly; see ``FEATURES_BLOCK_SIZE``. ``pool_size`` is passed the
    same way ``file_extract`` does, so s3fs reuses one filesystem.

    The block size has to go to ``fs.open``, not to ``fsspec.open``:
    ``fsspec.open`` forwards its keywords as *storage options* to the
    filesystem constructor, which takes ``default_block_size`` rather than
    ``block_size``, so a ``block_size`` passed there reaches the botocore
    session and raises ``TypeError``. ``fastq_chunker.s3io`` opens the same
    way for the same reason.
    """
    import fsspec
    import h5py

    fs = fsspec.filesystem("s3", config_kwargs={"max_pool_connections": pool_size})
    with fs.open(f"s3://{bucket}/{key}", "rb", block_size=FEATURES_BLOCK_SIZE) as fobj:
        with h5py.File(fobj, "r") as h5:
            return guide_ids_from_h5(h5)


def _library_dir(key: str) -> str:
    """The library folder a key sits under: the path before ``processed/``."""
    parts = key.split("/")
    if "processed" in parts:
        idx = parts.index("processed")
        if idx > 0:
            return "/".join(parts[:idx]) + "/"
    return ""


def _multi_group(key: str) -> str:
    """The library a run-level h5 counts towards.

    Normally the library folder. A key with no ``processed/`` segment still has
    to land in some group, so it falls back to whatever precedes ``outs/``;
    that keeps one h5 from being silently dropped.
    """
    library = _library_dir(key)
    if library:
        return library
    head, marker, _tail = key.partition(FLAT_OUTS_SEGMENT)
    return head if marker else key


def is_run_level_h5(key: str, filename: str) -> bool:
    """Whether ``key`` is a run-level ``filename`` in either Cell Ranger layout.

    Matched on path segments, not a string suffix, so ``per_sample_outs`` never
    counts as ``outs`` on a basename collision.
    """
    parts = key.split("/")
    if parts[-1] != filename:
        return False
    folders = parts[:-1]
    return folders[-1:] == ["outs"] or folders[-3:] == ["outs", "multi", "count"]


def find_multi_layout(keys: list[str], filename: str) -> MultiLayout:
    """Group the run-level ``filename`` keys by library, in either layout.

    Every library seen in the listing must hold exactly one. Libraries holding
    none, and libraries holding several (reprocessing left more than one run
    folder, so which one describes the delivery is undecidable), are returned
    for the caller to fail on rather than guessed at.
    """
    h5_keys = sorted(k for k in keys if is_run_level_h5(k, filename))

    library_dirs = {d for k in keys if (d := _library_dir(k))}
    by_library: dict[str, list[str]] = {}
    for key in h5_keys:
        by_library.setdefault(_multi_group(key), []).append(key)

    return MultiLayout(
        h5_keys=sorted(found[0] for found in by_library.values() if len(found) == 1),
        libraries_without_h5=sorted(library_dirs - set(by_library)),
        ambiguous_libraries=sorted(
            (library, sorted(found))
            for library, found in by_library.items()
            if len(found) > 1
        ),
    )


def find_per_sample_layout(keys: list[str], filename: str) -> PerSampleLayout:
    """Find per-sample h5 keys, and the folders that should hold one but do not.

    A sample folder is any ``.../per_sample_outs/<sample>/`` seen in the
    listing; one with no ``filename`` beneath it is reported, since it would
    otherwise be silently skipped. A library folder (the path before a
    ``processed`` segment, as in ``file_extract``) is reported when it holds no
    ``per_sample_outs/`` at all, which is a failed or half-uploaded run that no
    sample-folder check can see. Whether a library holds a run-level h5 is a
    separate question, answered by ``find_multi_layout``.
    """
    h5_keys = sorted(k for k in keys if k.rsplit("/", 1)[-1] == filename)

    sample_dirs: set[str] = set()
    library_dirs: set[str] = set()
    for key in keys:
        head, marker, tail = key.partition("per_sample_outs/")
        if marker:
            sample, slash, _rest = tail.partition("/")
            if sample and slash:
                sample_dirs.add(f"{head}per_sample_outs/{sample}/")
        if library := _library_dir(key):
            library_dirs.add(library)

    return PerSampleLayout(
        h5_keys=h5_keys,
        sample_dirs_without_h5=sorted(
            d for d in sample_dirs if not any(k.startswith(d) for k in h5_keys)
        ),
        library_dirs_without_per_sample_outs=sorted(
            lib
            for lib in library_dirs
            if not any(d.startswith(lib) for d in sample_dirs)
        ),
    )


def read_sample_names(path: str | Path) -> list[str]:
    """Read a sample list: one exact sample name per line, no header.

    Surrounding whitespace is stripped; blank lines and lines starting with
    ``#`` are skipped; repeated names are kept once, in first-seen order.

    Raises:
        H5CheckError: If the file cannot be read or is not UTF-8, a line holds
            a tab (a multi-column file, not one name per line), or no names
            remain.
    """
    file_path = Path(path)
    try:
        lines = file_path.read_text(encoding="utf-8-sig").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise H5CheckError(f"{file_path}: cannot read sample list: {exc}") from exc
    names: dict[str, None] = {}
    for number, line in enumerate(lines, start=1):
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        if "\t" in name:
            raise H5CheckError(
                f"{file_path}: line {number}: holds a tab; expected one sample "
                "name per line"
            )
        names.setdefault(name)
    if not names:
        raise H5CheckError(f"{file_path}: no sample names found")
    return list(names)


def _scope_name(key: str) -> str:
    """The name a sample list is matched against: sample, else library."""
    return extract_sample_name(key) or extract_library(key)


def _run_folder(key: str) -> str:
    """The run folder holding ``outs/``, as in ``.../cellranger/<run>/outs/``."""
    parts = key.split("/")
    if "outs" in parts:
        idx = parts.index("outs")
        if idx > 0:
            return parts[idx - 1]
    return ""


def _per_sample_label(key: str) -> str:
    """Name a per-sample h5 by its sample, falling back to library, then key."""
    sample = extract_sample_name(key)
    library = extract_library(key)
    if sample and library:
        return f"sample {sample} (library {library})"
    return f"sample {sample or library or key}"


def _multi_label(key: str) -> str:
    """Name a run-level h5 by library and run folder.

    The run folder is included because one library may hold several run
    folders; without it two h5s under one library read as the same thing.
    """
    library = extract_library(key)
    run = _run_folder(key)
    if library and run:
        return f"library {library} (run {run})"
    return f"library {library}" if library else key


def _examples(values: set[str]) -> str:
    return ", ".join(sorted(values)[:EXAMPLE_LIMIT])


def _make_s3_client(pool_size: int) -> Any:
    import boto3
    from botocore.config import Config

    return boto3.client("s3", config=Config(max_pool_connections=pool_size))


def _compare_h5_files(
    keys: list[str],
    *,
    bucket: str,
    guide_ids: set[str],
    reader: H5Reader,
    retries: int,
    pool_size: int,
    label: Callable[[str], str],
    report: H5CheckReport,
) -> None:
    """Read every key and record missing/extra guides against ``report``.

    Every h5 is read before anything is recorded, so one run shows every
    failing file rather than stopping at the first.
    """
    if not keys:
        return
    report.h5_total += len(keys)

    def fetch(key: str) -> tuple[set[str] | None, str | None]:
        return retry_with_backoff(reader, bucket, key, retries=retries)

    with ThreadPoolExecutor(max_workers=min(pool_size, len(keys))) as pool:
        results = list(pool.map(fetch, keys))

    for key, (h5_ids, err) in zip(keys, results):
        where = f"{label(key)} ({s3_uri_for(bucket, key)})"
        if err:
            report.h5_failed += 1
            report.failures.append(f"{where}: could not read h5: {err}")
            continue
        if not h5_ids:
            report.h5_failed += 1
            report.failures.append(
                f"{where}: no {CRISPR_FEATURE_TYPE!r} features in h5"
            )
            continue
        missing = guide_ids - h5_ids
        extra = h5_ids - guide_ids
        if missing:
            report.h5_failed += 1
            report.failures.append(
                f"{where}: {len(missing)} of {len(guide_ids)} guide_ids missing "
                f"from h5; examples: {_examples(missing)}"
            )
        if extra:
            report.notes.append(
                f"{where}: {len(extra)} guide(s) in h5 not in metadata; "
                f"examples: {_examples(extra)}"
            )


def _check_multi_scope(
    guide_ids: set[str],
    keys: list[str],
    *,
    bucket: str,
    h5_uri: str,
    filename: str,
    reader: H5Reader,
    retries: int,
    pool_size: int,
    samples: Collection[str] | None,
    samples_label: str,
    all_h5: int,
    report: H5CheckReport,
) -> None:
    """Check one run-level h5 per library.

    ``samples`` limits the check to those libraries. Filtering the keys before
    the layout is built does two jobs: only the named libraries are read, and
    the others stop being reported as missing a run-level h5. That second part
    is the point of the flag, since a library the sheet does not describe would
    otherwise fail either for lacking the sheet's guides or for lacking an h5.
    """
    if samples is not None:
        wanted = set(samples)
        keys = [k for k in keys if _scope_name(k) in wanted]
    layout = find_multi_layout(keys, filename)

    for library in layout.libraries_without_h5:
        report.failures.append(
            f"library folder {s3_uri_for(bucket, library)} has no run-level "
            f"{filename} under {RUN_LEVEL_LOCATIONS}"
        )
    for library, found in layout.ambiguous_libraries:
        locations = ", ".join(s3_uri_for(bucket, k) for k in found)
        report.failures.append(
            f"library folder {s3_uri_for(bucket, library)} has {len(found)} "
            f"run-level {filename} files; expected exactly one, so none were "
            f"read: {locations}"
        )
    if samples is not None:
        report.notes.append(
            f"checked {len(layout.h5_keys)} of {all_h5} run-level {filename} "
            f"file(s); libraries from {samples_label}"
        )
        found_names = {_scope_name(k) for k in layout.h5_keys}
        # A failure above already names a library whose folder was listed but
        # holds no usable h5, so this is kept for names with no folder at all.
        covered = {
            library.rstrip("/").rsplit("/", 1)[-1]
            for library in (
                layout.libraries_without_h5
                + [lib for lib, _found in layout.ambiguous_libraries]
            )
        }
        for name in samples:
            if name not in found_names and name not in covered:
                report.failures.append(
                    f"library {name!r} listed in {samples_label} has no "
                    f"run-level {filename} under {h5_uri}"
                )
    elif (
        not layout.h5_keys
        and not layout.libraries_without_h5
        and not layout.ambiguous_libraries
    ):
        report.failures.append(
            f"no run-level {filename} found under {h5_uri} "
            f"(looked under {RUN_LEVEL_LOCATIONS}); if this order has no "
            f"cellranger multi output, use --h5-scope {SCOPE_PER_SAMPLE}"
        )
        return

    _compare_h5_files(
        layout.h5_keys,
        bucket=bucket,
        guide_ids=guide_ids,
        reader=reader,
        retries=retries,
        pool_size=pool_size,
        label=_multi_label,
        report=report,
    )


def _check_per_sample_scope(
    guide_ids: set[str],
    keys: list[str],
    *,
    bucket: str,
    h5_uri: str,
    filename: str,
    reader: H5Reader,
    retries: int,
    pool_size: int,
    samples: Collection[str] | None,
    samples_label: str,
    all_h5: int,
    report: H5CheckReport,
) -> None:
    """Sweep every per-sample h5; a file integrity check, not a guide check."""
    if samples is not None:
        wanted = set(samples)
        keys = [k for k in keys if _scope_name(k) in wanted]
    layout = find_per_sample_layout(keys, filename)

    for folder in layout.library_dirs_without_per_sample_outs:
        report.failures.append(
            f"library folder {s3_uri_for(bucket, folder)} has no per_sample_outs/"
        )
    for folder in layout.sample_dirs_without_h5:
        report.failures.append(
            f"sample folder {s3_uri_for(bucket, folder)} has no {filename}"
        )
    if samples is not None:
        report.notes.append(
            f"checked {len(layout.h5_keys)} of {all_h5} {filename} file(s); "
            f"samples from {samples_label}"
        )
        found = {_scope_name(k) for k in layout.h5_keys}
        # A folder failure above already names a sample whose folder exists
        # but holds no h5, so the listed-sample message is kept for names with
        # no folder at all.
        covered = {
            folder.rstrip("/").rsplit("/", 1)[-1]
            for folder in (
                layout.sample_dirs_without_h5
                + layout.library_dirs_without_per_sample_outs
            )
        }
        for name in samples:
            if name not in found and name not in covered:
                report.failures.append(
                    f"sample {name!r} listed in {samples_label} has no {filename} "
                    f"under {h5_uri}"
                )
    if not layout.h5_keys:
        if samples is None:
            report.failures.append(f"no {filename} found under {h5_uri}")
        return

    _compare_h5_files(
        layout.h5_keys,
        bucket=bucket,
        guide_ids=guide_ids,
        reader=reader,
        retries=retries,
        pool_size=pool_size,
        label=_per_sample_label,
        report=report,
    )


def check_guide_ids_against_h5(
    guide_ids: set[str],
    h5_uri: str,
    *,
    scope: str = SCOPE_MULTI,
    filename: str | None = None,
    workers: int | None = None,
    retries: int = 5,
    samples: Collection[str] | None = None,
    samples_label: str = "the sample list",
    s3_client: Any = None,
    reader: H5Reader | None = None,
) -> H5CheckReport:
    """Check ``guide_ids`` against the h5 files under ``h5_uri``.

    One S3 listing serves every scope. ``scope`` selects which h5 files are
    read; see the module docstring. ``filename`` overrides that scope's default
    basename and is rejected under ``both``, where two different defaults are
    in play.

    ``samples`` limits every scope that runs to those exact names, so under
    ``both`` it narrows both passes. A name is matched against the
    ``per_sample_outs/<sample>/`` folder when there is one, else the library
    folder (the folder before ``processed/``); run-level keys have no
    ``per_sample_outs/`` segment, so under ``multi`` the names are library
    folders. A listed name with no h5 is a failure, so a typo is reported
    rather than quietly narrowing the check to nothing. Names not listed are
    excluded entirely, including from the missing-h5 reports: a library the
    sheet does not describe, such as one built from a different guide set,
    should drop out rather than fail. ``samples_label`` names the list in
    messages.

    ``s3_client`` and ``reader`` exist for tests; by default a boto3 client
    lists the prefix and each h5 is opened over S3.

    Raises:
        H5CheckError: If the scope is unknown, the URI is invalid, the prefix
            cannot be listed, or an argument does not apply to the scope.
    """
    if scope not in SCOPES:
        raise H5CheckError(
            f"unknown scope {scope!r}; expected one of {', '.join(SCOPES)}"
        )
    if filename is not None and scope == SCOPE_BOTH:
        raise H5CheckError(
            f"a filename override cannot apply to scope {SCOPE_BOTH!r}, which "
            "looks for two different h5 basenames"
        )
    try:
        location = parse_s3_uri(h5_uri)
    except ValueError as exc:
        raise H5CheckError(str(exc)) from exc

    pool_size = h5_worker_ceiling(
        do_introspect=True, workers=workers or default_workers_for(scope)
    )
    if s3_client is None or reader is None:
        check_dependencies()
    if s3_client is None:
        s3_client = _make_s3_client(pool_size)
    if reader is None:

        def reader(bucket: str, key: str) -> set[str]:
            return read_s3_h5_guide_ids(bucket, key, pool_size=pool_size)

    listed, list_err = retry_with_backoff(
        list_objects_with_size,
        s3_client,
        location.bucket,
        location.prefix,
        retries=retries,
    )
    if list_err:
        raise H5CheckError(f"could not list {h5_uri}: {list_err}")

    bucket = location.bucket
    keys = [obj.key for obj in listed or []]
    report = H5CheckReport()
    shared = {
        "bucket": bucket,
        "h5_uri": h5_uri,
        "reader": reader,
        "retries": retries,
        "pool_size": pool_size,
        "report": report,
    }

    if scope in (SCOPE_MULTI, SCOPE_BOTH):
        multi_filename = filename or DEFAULT_MULTI_H5_FILENAME
        _check_multi_scope(
            guide_ids,
            keys,
            filename=multi_filename,
            samples=samples,
            samples_label=samples_label,
            all_h5=sum(1 for k in keys if is_run_level_h5(k, multi_filename)),
            **shared,
        )
    if scope in (SCOPE_PER_SAMPLE, SCOPE_BOTH):
        per_sample_filename = filename or DEFAULT_PER_SAMPLE_H5_FILENAME
        _check_per_sample_scope(
            guide_ids,
            keys,
            filename=per_sample_filename,
            samples=samples,
            samples_label=samples_label,
            all_h5=sum(1 for k in keys if k.rsplit("/", 1)[-1] == per_sample_filename),
            **shared,
        )
    return report
