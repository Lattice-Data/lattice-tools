"""What integrity checks are actually available, probed rather than assumed.

Three facts drive this module:

* ``ListObjectsV2``, ``ListBucket`` and ``HeadObject`` are confirmed on the vendor
  buckets. ``GetObjectAttributes`` is *not* -- it has never been exercised, and
  ``ListObjectVersions`` needs ``s3:ListBucketVersions``, which is a different
  permission from the ``s3:ListBucket`` we know we have.
* The objects are SSE-KMS encrypted with a key in a third account, so a call can
  fail on the key policy rather than the bucket policy, with a different error.
* FASTQs and CRAMs are multipart, where the ETag is not the object's MD5.

So every capability here is probed and reported, and a denied capability produces
``available=False`` with the reason attached -- never a silent fallback. A report
that says "checksums verified" when the checksum call was denied is worse than
one that says nothing, because someone will believe it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

# Error codes that mean "this identity may not do that", as opposed to a
# transient failure or a bad request. Kept distinct so the report can separate
# "denied" from "broken", which are different conversations with CZI Biohub.
_DENIAL_CODES = frozenset(
    {
        "AccessDenied",
        "AccessDeniedException",
        "AllAccessDisabled",
        "Forbidden",
        "KMS.AccessDeniedException",
        "UnauthorizedOperation",
    }
)

# Raised when the API exists but this object does not have what was asked for --
# a single-part object has no part 1 metadata to return. Not a permission
# problem, and it must not be reported as one.
_NOT_APPLICABLE_CODES = frozenset({"InvalidPartNumber", "InvalidRequest", "404"})


@dataclass(frozen=True)
class Capability:
    """One probed S3 capability."""

    name: str
    available: bool
    detail: str
    error_code: str = ""
    denied: bool = False

    @property
    def status(self) -> str:
        if self.available:
            return "available"
        return "denied" if self.denied else "unavailable"


@dataclass
class CapabilityReport:
    """Every capability probed for one bucket, with the key used to probe."""

    bucket: str
    probe_key: str = ""
    capabilities: list[Capability] = field(default_factory=list)
    reason: str = ""

    def get(self, name: str) -> Capability | None:
        return next((c for c in self.capabilities if c.name == name), None)

    def has(self, name: str) -> bool:
        cap = self.get(name)
        return bool(cap and cap.available)

    def as_dict(self) -> dict[str, object]:
        return {
            "bucket": self.bucket,
            "probe_key": self.probe_key,
            "reason": self.reason,
            "capabilities": {
                c.name: {
                    "status": c.status,
                    "detail": c.detail,
                    "error_code": c.error_code,
                }
                for c in self.capabilities
            },
        }


def _error_code(exc: Exception) -> str:
    """Pull the S3 error code out of a botocore ClientError, if that is what it is."""
    response = getattr(exc, "response", None)
    if isinstance(response, dict):
        code = (response.get("Error") or {}).get("Code")
        if code:
            return str(code)
    return type(exc).__name__


def _probe(name: str, fn) -> Capability:
    """Run one probe, classifying failure as denial vs. breakage vs. not-applicable."""
    try:
        detail = fn()
    except Exception as exc:  # noqa: BLE001 - classified and reported
        code = _error_code(exc)
        if code in _NOT_APPLICABLE_CODES:
            return Capability(
                name=name,
                available=False,
                detail=f"not applicable for the probe object: {exc}",
                error_code=code,
            )
        denied = code in _DENIAL_CODES
        return Capability(
            name=name,
            available=False,
            detail=f"{'denied' if denied else 'failed'}: {exc}",
            error_code=code,
            denied=denied,
        )
    return Capability(name=name, available=True, detail=detail)


def probe_capabilities(
    s3_client, bucket: str, *, probe_key: str = "", prefix: str = ""
) -> CapabilityReport:
    """Probe the integrity-related S3 calls against a real key.

    Called once per run, against one object, before any per-object work is
    planned -- so the plan reflects what is permitted instead of discovering it
    thousands of calls in.
    """
    report = CapabilityReport(bucket=bucket, probe_key=probe_key)
    if not probe_key:
        report.reason = "no object available to probe against"
        return report

    report.capabilities.append(
        _probe(
            "head_object",
            lambda: _describe_head(s3_client.head_object(Bucket=bucket, Key=probe_key)),
        )
    )
    report.capabilities.append(
        _probe(
            "head_object_part_number",
            lambda: _describe_parts(
                s3_client.head_object(Bucket=bucket, Key=probe_key, PartNumber=1)
            ),
        )
    )
    report.capabilities.append(
        _probe(
            "get_object_attributes",
            lambda: _describe_attributes(
                s3_client.get_object_attributes(
                    Bucket=bucket,
                    Key=probe_key,
                    ObjectAttributes=["ETag", "ObjectSize", "Checksum", "ObjectParts"],
                )
            ),
        )
    )
    if prefix:
        report.capabilities.append(
            _probe(
                "list_object_versions",
                lambda: _describe_versions(
                    s3_client.list_object_versions(
                        Bucket=bucket, Prefix=prefix, MaxKeys=1
                    )
                ),
            )
        )
    return report


def _describe_head(resp: dict) -> str:
    version = resp.get("VersionId") or "none"
    return (
        f"ContentLength={resp.get('ContentLength')} "
        f"VersionId={version} "
        f"SSE={resp.get('ServerSideEncryption', 'none')}"
    )


def _describe_parts(resp: dict) -> str:
    parts = resp.get("PartsCount")
    return f"PartsCount={parts} part1_bytes={resp.get('ContentLength')}"


def _describe_attributes(resp: dict) -> str:
    checksum = resp.get("Checksum") or {}
    algorithms = sorted(k for k in checksum if k.startswith("Checksum"))
    parts = (resp.get("ObjectParts") or {}).get("TotalPartsCount")
    return (
        f"ObjectSize={resp.get('ObjectSize')} "
        f"TotalPartsCount={parts} "
        f"checksums={','.join(algorithms) if algorithms else 'none stored'}"
    )


def _describe_versions(resp: dict) -> str:
    versions = resp.get("Versions") or []
    markers = resp.get("DeleteMarkers") or []
    return f"sample versions={len(versions)} delete_markers={len(markers)}"


@dataclass
class VersionState:
    """Current version identity per key, and any delete markers found.

    Version *identity* rather than version *count*: counting needs
    ``ListObjectVersions``, which may be denied, whereas ``HeadObject`` returns
    the current ``VersionId`` and is confirmed. A key whose VersionId differs from
    the previous run was re-uploaded -- which is the thing versioning otherwise
    hides, because a silent re-upload of the same key leaves the key count and the
    byte total untouched.
    """

    checked: bool = False
    method: str = ""
    version_by_key: dict[str, str] = field(default_factory=dict)
    delete_markers: list[str] = field(default_factory=list)
    reason: str = ""
    covered: int = 0
    total: int = 0

    @property
    def partial(self) -> bool:
        return self.checked and self.covered < self.total


def collect_versions_via_listing(s3_client, bucket: str, prefix: str) -> VersionState:
    """Version identities from ListObjectVersions: one paginated call set, all keys.

    Preferred when permitted because it is O(pages) rather than O(objects), and it
    is the only way to see delete markers -- a "deleted" object still exists as a
    marker over a live version, so an order can look complete while carrying a key
    nobody can read.
    """
    state = VersionState(method="list_object_versions")
    latest: dict[str, str] = {}
    try:
        paginator = s3_client.get_paginator("list_object_versions")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            for version in page.get("Versions", []):
                if version.get("IsLatest"):
                    latest[version["Key"]] = version.get("VersionId", "")
            for marker in page.get("DeleteMarkers", []):
                state.delete_markers.append(marker["Key"])
    except Exception as exc:  # noqa: BLE001 - degraded, not fatal
        state.reason = f"{_error_code(exc)}: {exc}"
        return state
    state.checked = True
    state.version_by_key = latest
    state.covered = len(latest)
    state.total = len(latest)
    return state


def collect_versions_via_head(
    s3_client, bucket: str, keys: list[str], *, limit: int
) -> VersionState:
    """Version identities from HeadObject, one call per key, bounded by ``limit``.

    The fallback for when ``ListObjectVersions`` is denied. It costs one request
    per object, so it is capped -- and the cap is recorded in ``covered``/``total``
    so a partial sweep is reported as partial. A bounded check reported as a
    complete one is how an order gets released on the strength of the first 200 of
    its 4000 files.

    Delete markers are invisible to this method: HeadObject on a delete-marked key
    404s, which is indistinguishable here from a key that was never uploaded.
    """
    state = VersionState(method="head_object", total=len(keys))
    if not keys:
        state.checked = True
        return state
    versions: dict[str, str] = {}
    failures = 0
    for key in keys[:limit]:
        try:
            resp = s3_client.head_object(Bucket=bucket, Key=key)
        except Exception as exc:  # noqa: BLE001 - counted, reported below
            failures += 1
            if failures == 1:
                state.reason = f"first HeadObject failure: {_error_code(exc)}: {exc}"
            continue
        versions[key] = resp.get("VersionId") or ""
    state.checked = bool(versions) or not keys
    state.version_by_key = versions
    state.covered = len(versions)
    if failures:
        state.reason = (
            f"{state.reason} ({failures} of {min(len(keys), limit)} HeadObject "
            "calls failed)"
        )
    return state


def collect_version_state(
    s3_client,
    bucket: str,
    prefix: str,
    keys: list[str],
    *,
    capabilities: CapabilityReport,
    head_limit: int,
) -> VersionState:
    """Best available version state, preferring the listing and saying which was used."""
    if capabilities.has("list_object_versions"):
        state = collect_versions_via_listing(s3_client, bucket, prefix)
        if state.checked:
            return state
        log.warning(
            "ListObjectVersions probed available but failed (%s); "
            "falling back to per-object HeadObject.",
            state.reason,
        )
    if not capabilities.has("head_object"):
        return VersionState(
            reason=(
                "neither ListObjectVersions nor HeadObject is available, so "
                "re-uploads of identical keys cannot be detected"
            )
        )
    return collect_versions_via_head(s3_client, bucket, keys, limit=head_limit)
