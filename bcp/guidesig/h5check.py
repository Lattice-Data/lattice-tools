"""Check a guide metadata sheet's guide_ids against Cell Ranger h5 matrices.

Every ``guide_id`` in the sheet must appear among the ``CRISPR Guide Capture``
feature IDs of each sample's h5 under an S3 order prefix. Guides present in an
h5 but absent from the sheet are reported as notes, not failures. A sample
list (``read_sample_names``) can limit the check to named samples.

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

DEFAULT_H5_FILENAME = DEFAULT_H5_TARGET_FILENAME
CRISPR_FEATURE_TYPE = "CRISPR Guide Capture"
EXAMPLE_LIMIT = 10

# reader(bucket, key) -> CRISPR feature IDs in that h5.
H5Reader = Callable[[str, str], "set[str]"]


class H5CheckError(GuideSigError):
    """Raised when the h5 check cannot run at all (bad URI, listing failed)."""


@dataclass(frozen=True)
class H5Layout:
    """What one listing of an order prefix says about its h5 files."""

    h5_keys: list[str]
    sample_dirs_without_h5: list[str]
    library_dirs_without_h5: list[str]


@dataclass
class H5CheckReport:
    """Outcome of checking one sheet's guide_ids against every h5 in an order."""

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

    Only the features table is read, never the count matrix. ``pool_size`` is
    passed the same way ``file_extract`` does, so s3fs reuses one filesystem.
    """
    import fsspec
    import h5py

    with fsspec.open(
        f"s3://{bucket}/{key}",
        "rb",
        config_kwargs={"max_pool_connections": pool_size},
    ) as fobj:
        with h5py.File(fobj, "r") as h5:
            return guide_ids_from_h5(h5)


def find_h5_layout(keys: list[str], filename: str) -> H5Layout:
    """Find the target h5 keys, and the folders that should hold one but do not.

    A sample folder is any ``.../per_sample_outs/<sample>/`` seen in the
    listing. A library folder is the path before a ``processed`` segment, as
    in ``file_extract``. Either one with no ``filename`` beneath it is
    reported, since such a folder would otherwise be silently skipped.
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
        parts = key.split("/")
        if "processed" in parts:
            idx = parts.index("processed")
            if idx > 0:
                library_dirs.add("/".join(parts[:idx]) + "/")

    def without_h5(dirs: set[str]) -> list[str]:
        return sorted(d for d in dirs if not any(k.startswith(d) for k in h5_keys))

    return H5Layout(
        h5_keys=h5_keys,
        sample_dirs_without_h5=without_h5(sample_dirs),
        library_dirs_without_h5=without_h5(library_dirs),
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


def _label(key: str) -> str:
    """Name an h5 by its sample, falling back to library, then the key."""
    sample = extract_sample_name(key)
    library = extract_library(key)
    if sample and library:
        return f"{sample} (library {library})"
    return sample or library or key


def _examples(values: set[str]) -> str:
    return ", ".join(sorted(values)[:EXAMPLE_LIMIT])


def _make_s3_client(pool_size: int) -> Any:
    import boto3
    from botocore.config import Config

    return boto3.client("s3", config=Config(max_pool_connections=pool_size))


def check_guide_ids_against_h5(
    guide_ids: set[str],
    h5_uri: str,
    *,
    filename: str = DEFAULT_H5_FILENAME,
    workers: int | None = None,
    retries: int = 5,
    samples: Collection[str] | None = None,
    samples_label: str = "the sample list",
    s3_client: Any = None,
    reader: H5Reader | None = None,
) -> H5CheckReport:
    """Check ``guide_ids`` against every ``filename`` under ``h5_uri``.

    Every h5 is checked before anything is reported, so one run shows every
    failing sample. ``samples``, when given, limits every part of the check
    (h5 reads and missing-folder reports) to those exact names, matched against
    the ``per_sample_outs/<sample>/`` folder (else the library folder); a
    listed name with no h5 is a failure. ``samples_label`` names the list in
    messages. ``s3_client`` and ``reader`` exist for tests; by default a boto3
    client lists the prefix and each h5 is opened over S3.

    Raises:
        H5CheckError: If the URI is invalid or the prefix cannot be listed.
    """
    try:
        location = parse_s3_uri(h5_uri)
    except ValueError as exc:
        raise H5CheckError(str(exc)) from exc

    pool_size = h5_worker_ceiling(do_introspect=True, workers=workers)
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
    all_h5 = sum(1 for k in keys if k.rsplit("/", 1)[-1] == filename)
    if samples is not None:
        wanted = set(samples)
        keys = [k for k in keys if _scope_name(k) in wanted]
    layout = find_h5_layout(keys, filename)
    report = H5CheckReport(h5_total=len(layout.h5_keys))

    for folder in layout.library_dirs_without_h5:
        report.failures.append(
            f"library folder {s3_uri_for(bucket, folder)} has no {filename}"
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
                layout.sample_dirs_without_h5 + layout.library_dirs_without_h5
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
        return report

    def fetch(key: str) -> tuple[set[str] | None, str | None]:
        return retry_with_backoff(reader, bucket, key, retries=retries)

    with ThreadPoolExecutor(max_workers=min(pool_size, len(layout.h5_keys))) as pool:
        results = list(pool.map(fetch, layout.h5_keys))

    for key, (h5_ids, err) in zip(layout.h5_keys, results):
        where = f"sample {_label(key)} ({s3_uri_for(bucket, key)})"
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

    return report
