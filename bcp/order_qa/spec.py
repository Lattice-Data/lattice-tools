"""Resolving an order spec to a bucket and prefix.

Buckets and project prefixes are declared here, never discovered: this identity
is denied ``ListAllMyBuckets`` and every bucket-level metadata call, so the only
S3 operation available is a prefix-scoped listing of a path we already knew.

Order identifiers are deliberately matched loosely. Psomagen uses ``AN########``
and Novogene ``NVUS...-##``, but a vendor is free to invent a third shape and an
order QA cannot refuse to run on one it has not seen. So what is enumerated here
is the set of things that are *not* orders -- which is a closed, observed list --
and anything else under a watched project is treated as an order.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import urlparse

# Vendor buckets, and the projects under each that hold orders. Used to resolve a
# bare order ID; a fully-qualified spec needs no entry here and is not checked
# against it, so an order under a project nobody has added still runs.
WATCHED_PROJECTS: dict[str, tuple[str, ...]] = {
    "czi-psomagen": (
        "marson-crispra",
        "marson-macrophages-tregs-pilot",
        "marson-mapping-grns-perturb-seq",
    ),
    "czi-novogene": (),
    "czi-ultima": (),
}

# Children of a watched project that are not orders. Observed, not inferred.
# ``jsonlds`` and ``metadata-tsv`` are sibling data products; the ``test*`` names
# are vendor connectivity smoke tests that get refreshed periodically, so they
# will reappear after anyone deletes them.
NON_ORDER_NAMES: frozenset[str] = frozenset(
    {
        "jsonlds",
        "metadata-tsv",
        "test",
    }
)

# Smoke-test objects and folders share this prefix across both vendor buckets
# (``test.txt``, ``test/``, ``test_marson-crispra``).
_SMOKE_TEST_RE = re.compile(r"^test(?:[_.-]|$)", re.IGNORECASE)

# Recognised order shapes. Used only to *describe* a spec in the report, never to
# reject one -- see the module docstring.
_KNOWN_ORDER_SHAPES: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("psomagen_AN", re.compile(r"^AN\d{8}$")),
    ("novogene_NVUS", re.compile(r"^NVUS\d+-\d+$")),
)


class SpecError(ValueError):
    """The spec could not be resolved to a single order prefix."""


@dataclass(frozen=True)
class OrderSpec:
    """A resolved order. Identity is (bucket, project, order), never order alone.

    Two vendors can and do issue the same-looking identifier, and the same
    ExperimentID appears under both a lab bucket and the vendor bucket it was
    trimmed from, so an order ID on its own does not name a thing to QA.
    """

    bucket: str
    project: str
    order: str
    order_shape: str

    @property
    def prefix(self) -> str:
        """Listing prefix, trailing slash included so it cannot match a sibling.

        Without the slash, ``AN00028221`` also lists ``AN000282210``; with it the
        listing is exactly this order.
        """
        return f"{self.project}/{self.order}/"

    @property
    def s3_path(self) -> str:
        """``s3://`` form, which is what ``resolve_qa_run_context`` takes."""
        return f"s3://{self.bucket}/{self.project}/{self.order}"

    @property
    def key(self) -> str:
        """Stable identity for reports and re-run comparison."""
        return f"{self.bucket}/{self.project}/{self.order}"


def is_non_order_name(name: str) -> bool:
    """True if a project child is known not to be an order.

    Callers listing a project should log and skip these rather than raise: a
    watched project legitimately contains non-orders, and a new sibling data
    product appearing must not break QA for the orders beside it.
    """
    cleaned = name.strip().strip("/")
    if not cleaned:
        return True
    if cleaned.lower() in NON_ORDER_NAMES:
        return True
    return bool(_SMOKE_TEST_RE.match(cleaned))


def classify_order_shape(order: str) -> str:
    """Name the order-ID shape, or ``"unrecognized"``.

    Reported rather than enforced, so a new vendor format shows up as a note on a
    run that happened instead of an error on one that did not.
    """
    for label, pattern in _KNOWN_ORDER_SHAPES:
        if pattern.match(order):
            return label
    return "unrecognized"


def _split_spec(raw: str) -> list[str]:
    """Split a spec into path segments, accepting ``s3://`` or bare form."""
    text = raw.strip()
    if not text:
        raise SpecError("Empty order spec.")
    if "://" in text:
        parsed = urlparse(text)
        if parsed.scheme != "s3" or not parsed.netloc:
            raise SpecError(f"Only s3:// URIs are supported, got {text!r}.")
        return [parsed.netloc, *[p for p in parsed.path.split("/") if p]]
    return [p for p in text.split("/") if p]


def resolve_spec(raw: str) -> OrderSpec:
    """Resolve ``bucket/project/order``, ``s3://bucket/project/order``, or a bare
    order ID, to a single :class:`OrderSpec`.

    A bare ID is resolved by *name*, against ``WATCHED_PROJECTS`` -- not by
    listing S3 to see which project contains it. Searching would need a listing
    per candidate project and would still be ambiguous; more to the point, an ID
    that matches two watched projects is a question for the caller, not something
    to resolve by picking whichever listing answered first.
    """
    segments = _split_spec(raw)

    if len(segments) == 1:
        return _resolve_bare_order(segments[0])

    if len(segments) == 2:
        raise SpecError(
            f"Spec {raw!r} names a project but no order. Pass "
            "bucket/project/order, or a bare order ID."
        )

    # Deeper paths are accepted -- someone pasting a key out of the console gets
    # the order it belongs to. Project and order are taken from their fixed
    # positions in {bucket}/{project}/{order}/..., not off the tail: the tail of
    # a key is the filename, so reading backwards turned a pasted FASTQ path into
    # the bucket plus its last two path segments and QA'd a prefix nobody named.
    # The caller logs the resolved prefix, so truncation is visible.
    bucket, project, order = segments[0], segments[1], segments[2]

    if is_non_order_name(order):
        raise SpecError(
            f"{order!r} under {bucket}/{project} is not an order "
            "(sibling data product or vendor smoke test)."
        )
    return OrderSpec(
        bucket=bucket,
        project=project,
        order=order,
        order_shape=classify_order_shape(order),
    )


def _resolve_bare_order(order: str) -> OrderSpec:
    """Resolve a bare order ID against the declared watched projects."""
    if is_non_order_name(order):
        raise SpecError(f"{order!r} is not an order identifier.")

    shape = classify_order_shape(order)
    candidates = [
        (bucket, project)
        for bucket, projects in WATCHED_PROJECTS.items()
        for project in projects
    ]
    if not candidates:
        raise SpecError(
            f"No watched projects are configured, so bare order ID {order!r} "
            "cannot be resolved. Pass bucket/project/order."
        )
    if len(candidates) > 1:
        listing = ", ".join(f"{b}/{p}" for b, p in candidates)
        raise SpecError(
            f"Bare order ID {order!r} is ambiguous across {len(candidates)} "
            f"watched projects ({listing}). Pass bucket/project/order."
        )
    bucket, project = candidates[0]
    return OrderSpec(bucket=bucket, project=project, order=order, order_shape=shape)
