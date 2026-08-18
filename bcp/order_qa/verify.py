"""Upload verification: what is in S3, when it last changed, and what moved.

Separate from the QA checks, and deliberately ahead of them. QA answers "is this
data any good"; the questions here are "did all of it arrive" and "has it stopped
arriving", and asking them second wastes a full gather on a half-uploaded order
and reports defects that are really just absences.

Everything here uses ListObjectsV2 and HeadObject only. Both authorize against
prefix-scoped ``s3:ListBucket`` / ``s3:GetObject``, which are confirmed available
on the vendor buckets; no bucket-level metadata call is ever made, because every
one of those is denied.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# A multipart ETag carries a ``-<partcount>`` suffix. For those the ETag is a
# digest of part digests, not the object's MD5, so it cannot be compared against
# a vendor-supplied MD5 -- see ``etag_is_md5``.
_MULTIPART_ETAG_RE = re.compile(r"-\d+$")

_S3_URI_RE = re.compile(r"s3://([^/\s,\t\"']+)/([^\s,\t\"']+)")

DEFAULT_QUIESCENCE_MINUTES = 15


@dataclass(frozen=True)
class S3Object:
    """One object from a listing."""

    key: str
    size_bytes: int
    last_modified: datetime
    etag: str
    storage_class: str = "STANDARD"

    @property
    def is_folder_marker(self) -> bool:
        """True for the zero-byte ``dir/`` keys the console and some tools create.

        These are not deliverables and must not be counted as files, or an order
        gains phantom "extra" objects that no vendor uploaded.
        """
        return self.key.endswith("/") and self.size_bytes == 0

    @property
    def etag_is_md5(self) -> bool:
        """True only when the ETag is genuinely this object's MD5.

        FASTQs and CRAMs arrive multipart, where it is not. Callers must gate any
        checksum claim on this, so that a report never says "checksum verified"
        about a comparison that could not have been made.
        """
        return not _MULTIPART_ETAG_RE.search(self.etag.strip('"'))

    @property
    def multipart_parts(self) -> int | None:
        """Part count parsed off a multipart ETag, or None if single-part.

        Free -- it comes from the listing -- and it is the cheap half of what
        ``GetObjectAttributes`` would tell us, so it is worth having even when
        that call turns out to be denied.
        """
        cleaned = self.etag.strip('"')
        match = _MULTIPART_ETAG_RE.search(cleaned)
        return int(match.group()[1:]) if match else None


@dataclass
class Listing:
    """The result of listing one order prefix."""

    bucket: str
    prefix: str
    objects: list[S3Object] = field(default_factory=list)
    complete: bool = True
    error: str | None = None

    @property
    def files(self) -> list[S3Object]:
        """Deliverables: every object that is not a folder marker."""
        return [o for o in self.objects if not o.is_folder_marker]

    @property
    def total_bytes(self) -> int:
        return sum(o.size_bytes for o in self.files)

    @property
    def newest(self) -> datetime | None:
        return max((o.last_modified for o in self.files), default=None)

    @property
    def keys(self) -> set[str]:
        return {o.key for o in self.files}

    def by_key(self) -> dict[str, S3Object]:
        return {o.key: o for o in self.files}


def list_order_objects(s3_client, bucket: str, prefix: str) -> Listing:
    """List one order prefix with ListObjectsV2, paginated.

    A listing that failed part way is returned with ``complete=False`` rather
    than raised. Everything downstream treats an incomplete listing as "no
    verdict": a partial listing looks exactly like a partial upload, and the one
    outcome we must never produce is a clean bill of health derived from half the
    objects.
    """
    listing = Listing(bucket=bucket, prefix=prefix)
    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            for obj in page.get("Contents", []):
                listing.objects.append(
                    S3Object(
                        key=obj["Key"],
                        size_bytes=obj["Size"],
                        last_modified=_as_utc(obj["LastModified"]),
                        etag=obj.get("ETag", ""),
                        storage_class=obj.get("StorageClass", "STANDARD"),
                    )
                )
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        listing.complete = False
        listing.error = f"{type(exc).__name__}: {exc}"
    return listing


def _as_utc(value: datetime) -> datetime:
    """Normalize to timezone-aware UTC.

    boto3 returns aware datetimes, but a stub or a replayed fixture may not, and
    subtracting a naive from an aware datetime raises -- which would turn the
    quiescence check into a crash on exactly the inputs used to test it.
    """
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


@dataclass(frozen=True)
class Quiescence:
    """Whether the upload has stopped changing."""

    checked: bool
    quiet: bool
    newest: datetime | None
    age_seconds: float | None
    threshold_seconds: int
    forced: bool = False

    @property
    def summary(self) -> str:
        if not self.checked:
            return "no objects to date-check"
        if self.forced:
            return f"quiescence forced past (newest {self.age_human} old)"
        if self.quiet:
            return f"quiet: newest object {self.age_human} old"
        return (
            f"NOT quiet: newest object {self.age_human} old, "
            f"under the {self.threshold_seconds // 60}m threshold"
        )

    @property
    def age_human(self) -> str:
        if self.age_seconds is None:
            return "unknown"
        minutes = self.age_seconds / 60
        if minutes < 60:
            return f"{minutes:.0f}m"
        return f"{minutes / 60:.1f}h"


def assess_quiescence(
    listing: Listing,
    *,
    threshold_minutes: int = DEFAULT_QUIESCENCE_MINUTES,
    force: bool = False,
    now: datetime | None = None,
) -> Quiescence:
    """Is the newest object older than the threshold?

    This is the whole distinction between "the vendor said the order is done" and
    "the order is actually done". Without it QA races the upload and reports
    missing files that are merely late.

    An empty prefix returns ``checked=False``: there is no timestamp to reason
    about. It is emphatically not quiet -- callers treat "nothing there" as a
    verification failure, not as a settled upload.
    """
    threshold_seconds = max(0, threshold_minutes) * 60
    newest = listing.newest
    if newest is None:
        return Quiescence(
            checked=False,
            quiet=False,
            newest=None,
            age_seconds=None,
            threshold_seconds=threshold_seconds,
            forced=force,
        )
    reference = _as_utc(now) if now is not None else datetime.now(timezone.utc)
    age = (reference - newest).total_seconds()
    return Quiescence(
        checked=True,
        # --force records that the gate was skipped; it does not fake a pass, so
        # the report can still say the upload may have been in flight.
        quiet=force or age >= threshold_seconds,
        newest=newest,
        age_seconds=age,
        threshold_seconds=threshold_seconds,
        forced=force,
    )


@dataclass
class ManifestComparison:
    """Listing vs. an expected-key manifest."""

    checked: bool
    expected: int = 0
    present: int = 0
    missing: list[str] = field(default_factory=list)
    extra: list[str] = field(default_factory=list)
    size_mismatches: list[tuple[str, int, int]] = field(default_factory=list)
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.checked and not self.missing and not self.size_mismatches


def read_manifest_keys(path: Path) -> dict[str, int | None]:
    """Read expected ``s3://bucket/key`` URIs, and sizes when the manifest has them.

    Accepts the same shape ``qa_mods`` manifest mode does -- CSV or TSV, with or
    without a header -- by scanning for ``s3://`` URIs anywhere in the file rather
    than by column index. A manifest is a hand-assembled artifact whose column
    order varies, and a key found in the wrong column is still a key; requiring a
    fixed index would report every file missing over a shifted column.

    A bare integer in another field on the same line is taken as that object's
    expected size, which is the common vendor layout.
    """
    expected: dict[str, int | None] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        # Every URI on the line, not just the first: a line holding two keys is
        # unusual but legal, and taking only the first silently drops the second
        # from the expected set -- which turns a missing file into one that looks
        # present, the exact failure this comparison exists to catch.
        matches = _S3_URI_RE.findall(line)
        if not matches:
            continue
        keys = [key for _bucket, key in matches]
        size: int | None = None
        # A size is attributed only when the line leaves no choice about which
        # number it is: exactly one object named, and exactly one bare integer
        # anywhere on the line. Taking the first integer instead read the lane out
        # of ``sampleA,1,s3://...,1048576`` and then reported every file in the
        # order as the wrong size. An unattributed size is checked by the
        # convention-driven checks anyway; a misattributed one invents failures.
        if len(keys) == 1:
            numbers = [
                token
                for field_text in re.split(r"[,\t]", line)
                if (token := field_text.strip().strip('"')).isdigit()
                and "s3://" not in field_text
            ]
            if len(numbers) == 1:
                size = int(numbers[0])
        for key in keys:
            expected[key] = size
    return expected


def compare_to_manifest(
    listing: Listing, manifest_path: Path | None
) -> ManifestComparison:
    """Compare the listing against a manifest of expected keys.

    With no manifest this reports ``checked=False`` and a reason, never a pass.
    There is no per-order manifest of expected contents in the vendor buckets --
    what an order should contain is derived from naming convention per assay by
    ``qa_checks.check_expected_raw_files`` -- so the absence of one is the normal
    case, and it has to be visible in the report rather than read as "verified".
    """
    if manifest_path is None:
        return ManifestComparison(
            checked=False,
            present=len(listing.files),
            reason=(
                "no --manifest given; expected contents come from the assay's "
                "naming convention instead (see the raw-file completeness check)"
            ),
        )
    try:
        expected = read_manifest_keys(manifest_path)
    except OSError as exc:
        return ManifestComparison(
            checked=False,
            present=len(listing.files),
            reason=f"manifest unreadable: {exc}",
        )
    if not expected:
        return ManifestComparison(
            checked=False,
            present=len(listing.files),
            reason=f"no s3:// URIs found in {manifest_path}",
        )

    by_key = listing.by_key()
    present_keys = set(by_key)
    expected_keys = set(expected)
    size_mismatches = [
        (key, size, by_key[key].size_bytes)
        for key, size in expected.items()
        if size is not None and key in by_key and by_key[key].size_bytes != size
    ]
    return ManifestComparison(
        checked=True,
        expected=len(expected_keys),
        present=len(present_keys),
        missing=sorted(expected_keys - present_keys),
        extra=sorted(present_keys - expected_keys),
        size_mismatches=sorted(size_mismatches),
    )
