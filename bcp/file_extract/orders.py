"""Re-sequencing orders: parse them, find each file's GroupID, and route files.

An order is sometimes sequenced again under a new CRO order. ``fastq``, ``cram``
and ``h5`` can take those re-sequencing orders beside the main one, and each
re-sequenced file is routed by its GroupID, the first folder below its order:

- a GroupID the main order also has: the file joins the main order's outputs,
  next to the main order's own files, never replacing them;
- a GroupID new to the main order: the file goes into an output set of its own,
  named after the re-sequencing order it came from.

Routing by GroupID is only sound in the 10x layout, where that first folder
really is the GroupID::

    {order}/{GroupID}/raw/{file}              (fastq, cram)
    {order}/{GroupID}/processed/...           (h5)

Scale and sci deliveries put an experiment id there and the GroupID only in the
filename (``{order}/{experiment_id}/raw/{RunID}/{file}``), so every file is
checked against the layout and any other shape is refused rather than routed by
the wrong folder.

Pure by design: nothing in this module talks to S3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

from .models import ListedObject
from .s3_utils import parse_s3_uri, s3_uri_for

# Folder names that can sit right below an order but are never a GroupID. Seeing
# one means the prefix given was not the order folder itself.
STRUCTURAL_FOLDERS = frozenset({"raw", "processed", "per_sample_outs"})

# How many offending paths an error names before summarising the rest.
MAX_REPORTED_PATHS = 5


class OrderError(ValueError):
    """Raised when the orders given cannot be routed safely."""


@dataclass(frozen=True)
class OrderPrefix:
    """One S3 order folder and the CRO order identifier its rows report."""

    bucket: str
    prefix: str
    cro_order: str
    is_primary: bool = False

    @property
    def uri(self) -> str:
        return f"s3://{self.bucket}/{self.prefix}"

    @property
    def role(self) -> str:
        return "main" if self.is_primary else "re-sequencing"


def order_id_from_prefix(prefix: str) -> str:
    """The order identifier: the last folder of the prefix."""
    return prefix.rstrip("/").rsplit("/", 1)[-1]


def split_order_values(values: Sequence[str] | None) -> list[str]:
    """Split ``nargs="+"`` values that may also be comma-separated."""
    tokens: list[str] = []
    for value in values or ():
        tokens.extend(part.strip() for part in (value or "").split(","))
    return tokens


def _nested(a: str, b: str) -> bool:
    return a.startswith(b) or b.startswith(a)


def parse_order_prefixes(
    main_uri: str,
    values: Sequence[str] | None,
    *,
    main_cro_order: str | None = None,
    validate_order_id: Callable[[str], str] | None = None,
) -> list[OrderPrefix]:
    """The main order followed by each ``--resequencing-orders`` value.

    ``main_cro_order`` is the main order's identifier when the mode takes one
    (``--cro-order``); otherwise it is the main prefix's last folder. It is used
    as given, since the caller validates its own flag. A re-sequencing order's
    identifier is always its prefix's last folder, passed through
    ``validate_order_id`` when one is supplied.

    Everything is checked here, before any S3 call: a re-sequencing order must
    be a folder (not a whole bucket) in the main order's bucket, must not
    overlap the main order or another re-sequencing order, and must not repeat
    an identifier.
    """
    try:
        main = parse_s3_uri(main_uri)
    except ValueError as exc:
        raise OrderError(str(exc)) from exc
    if not main.prefix:
        raise OrderError(
            f"The main order {main_uri!r} names a whole bucket; with "
            "--resequencing-orders it must be the order folder, because each "
            "file's GroupID is the first folder below it"
        )
    main_id = (main_cro_order or "").strip() or order_id_from_prefix(main.prefix)
    orders = [
        OrderPrefix(
            bucket=main.bucket,
            prefix=main.prefix,
            cro_order=main_id,
            is_primary=True,
        )
    ]
    # The main order's folder name is reserved too: its info TSV is named
    # after it, so a re-sequencing order of that name would overwrite it.
    taken_ids = {main_id: main_uri, order_id_from_prefix(main.prefix): main_uri}

    tokens = split_order_values(values)
    if not tokens:
        raise OrderError("--resequencing-orders must name at least one order")
    for token in tokens:
        if not token:
            raise OrderError("--resequencing-orders must not contain empty values")
        try:
            location = parse_s3_uri(token)
        except ValueError as exc:
            raise OrderError(f"Invalid --resequencing-orders {token!r}: {exc}") from exc
        if not location.prefix:
            raise OrderError(
                f"Invalid --resequencing-orders {token!r}: name the order "
                "folder; a bucket alone would be searched whole and could "
                "match files from an unrelated delivery"
            )
        if location.bucket != main.bucket:
            raise OrderError(
                f"Re-sequencing order {token!r} is in bucket "
                f"{location.bucket!r}, but the main order is in "
                f"{main.bucket!r}; every order must be in the same bucket"
            )
        for existing in orders:
            if location.prefix == existing.prefix:
                raise OrderError(
                    f"{token!r} is given more than once, or is the main order"
                )
            if _nested(location.prefix, existing.prefix):
                raise OrderError(
                    f"Re-sequencing order {token!r} overlaps "
                    f"{existing.uri!r}; one sits inside the other, so their "
                    "files would be counted twice"
                )
        order_id = order_id_from_prefix(location.prefix)
        if validate_order_id is not None:
            order_id = validate_order_id(order_id)
        if order_id in taken_ids:
            raise OrderError(
                f"Re-sequencing order {token!r} has order ID {order_id!r}, "
                f"which {taken_ids[order_id]!r} already uses; output files "
                "are named by order ID, so they would overwrite each other"
            )
        taken_ids[order_id] = token
        orders.append(
            OrderPrefix(
                bucket=location.bucket,
                prefix=location.prefix,
                cro_order=order_id,
            )
        )
    return orders


def group_id_for(key: str, prefix: str) -> str:
    """The first folder of ``key`` below ``prefix``, or ``""`` if none."""
    if not key.startswith(prefix):
        raise ValueError(f"{key!r} is not under {prefix!r}")
    head, slash, _rest = key[len(prefix) :].partition("/")
    return head if slash else ""


@dataclass(frozen=True)
class Layout:
    """The 10x shape a mode's files must have, below the order folder."""

    name: str
    description: str
    matches: Callable[[list[str]], bool]


# fastq and cram: the file sits directly inside {GroupID}/raw/.
RAW_LAYOUT = Layout(
    name="raw",
    description="{GroupID}/raw/{file}",
    matches=lambda parts: len(parts) == 3 and parts[1] == "raw",
)

# h5: Cell Ranger output sits somewhere under {GroupID}/processed/.
PROCESSED_LAYOUT = Layout(
    name="processed",
    description="{GroupID}/processed/.../{file}",
    matches=lambda parts: len(parts) >= 3 and parts[1] == "processed",
)


def _paths_block(uris: Sequence[str]) -> str:
    shown = list(uris[:MAX_REPORTED_PATHS])
    more = len(uris) - len(shown)
    lines = "".join(f"\n  {uri}" for uri in shown)
    if more:
        lines += f"\n  ... and {more} more"
    return lines


def check_layout(
    listings: Sequence[tuple[OrderPrefix, Sequence[ListedObject]]],
    layout: Layout,
) -> None:
    """Refuse any file whose path does not have the 10x shape.

    Every file in every order is checked, so an order that mixes layouts is
    caught as well. Raises one OrderError naming each kind of problem found.
    """
    no_group: list[str] = []
    structural: list[str] = []
    other_shape: list[str] = []
    for order, objects in listings:
        for obj in objects:
            uri = s3_uri_for(order.bucket, obj.key)
            parts = obj.key[len(order.prefix) :].split("/")
            if len(parts) == 1:
                no_group.append(uri)
            elif parts[0] in STRUCTURAL_FOLDERS:
                structural.append(uri)
            elif not layout.matches(parts):
                other_shape.append(uri)

    problems: list[str] = []
    if no_group:
        problems.append(
            f"{len(no_group)} file(s) sit directly in an order folder, with no "
            "GroupID folder, so no GroupID can be determined; the order folder "
            f"should contain GroupID folders:{_paths_block(no_group)}"
        )
    if structural:
        problems.append(
            f"{len(structural)} file(s) would get a GroupID of 'raw', "
            "'processed' or 'per_sample_outs'; each S3 path must be the order "
            f"folder itself, not a folder inside it:{_paths_block(structural)}"
        )
    if other_shape:
        problems.append(
            f"{len(other_shape)} file(s) are not in the 10x layout "
            f"{{order}}/{layout.description}; --resequencing-orders supports "
            "only the 10x layout (Scale and sci deliveries are not "
            f"supported):{_paths_block(other_shape)}"
        )
    if problems:
        raise OrderError("\n".join(problems))


@dataclass(frozen=True)
class RoutedFile:
    """A listed file, the order it came from, and its GroupID."""

    order: OrderPrefix
    obj: ListedObject
    group_id: str


@dataclass
class OutputSet:
    """Files written to one set of outputs, named after ``owner``."""

    owner: OrderPrefix
    files: list[RoutedFile] = field(default_factory=list)


@dataclass
class OrderReport:
    """What routing did with one order's files, for the run summary."""

    order: OrderPrefix
    matched: int = 0
    # Re-sequencing orders only: GroupID -> file count.
    added_to_main: dict[str, int] = field(default_factory=dict)
    new_groups: dict[str, int] = field(default_factory=dict)
    # Main order only.
    group_ids: list[str] = field(default_factory=list)


@dataclass
class Routing:
    """Output sets, main first, and one report per order in the order given."""

    output_sets: list[OutputSet]
    reports: list[OrderReport]
    warnings: list[str] = field(default_factory=list)
    empty_orders: list[OrderPrefix] = field(default_factory=list)


def _normalized(group_id: str) -> str:
    return group_id.lower().replace("-", "_")


def route(
    listings: Sequence[tuple[OrderPrefix, Sequence[ListedObject]]],
) -> Routing:
    """Route every listed file into the main output set or its own order's.

    ``listings`` holds the main order first, then each re-sequencing order.
    The main order's GroupIDs come from the files it matched. An output set is
    only made for a re-sequencing order that brought new GroupIDs.
    """
    if not listings or not listings[0][0].is_primary:
        raise ValueError("route() needs the main order first")
    main_order, main_objects = listings[0]
    main_set = OutputSet(owner=main_order)
    main_report = OrderReport(order=main_order, matched=len(main_objects))
    for obj in main_objects:
        main_set.files.append(
            RoutedFile(main_order, obj, group_id_for(obj.key, main_order.prefix))
        )
    main_groups = {routed.group_id for routed in main_set.files}
    main_report.group_ids = sorted(main_groups)
    main_by_normalized = {_normalized(gid): gid for gid in main_groups}

    routing = Routing(output_sets=[main_set], reports=[main_report])
    if not main_objects:
        routing.warnings.append(
            f"The main order {main_order.cro_order} ({main_order.uri}) matched "
            "no files, so every re-sequenced file goes into its own order's "
            "output set"
        )

    for order, objects in listings[1:]:
        report = OrderReport(order=order, matched=len(objects))
        routing.reports.append(report)
        if not objects:
            routing.empty_orders.append(order)
            routing.warnings.append(
                f"Re-sequencing order {order.cro_order} ({order.uri}) matched no files"
            )
            continue
        own_set = OutputSet(owner=order)
        for obj in objects:
            routed = RoutedFile(order, obj, group_id_for(obj.key, order.prefix))
            if routed.group_id in main_groups:
                main_set.files.append(routed)
                tally = report.added_to_main
            else:
                own_set.files.append(routed)
                tally = report.new_groups
            tally[routed.group_id] = tally.get(routed.group_id, 0) + 1
        if own_set.files:
            routing.output_sets.append(own_set)
        if not report.added_to_main and main_groups:
            routing.warnings.append(
                f"None of re-sequencing order {order.cro_order}'s "
                f"{len(report.new_groups)} GroupID(s) are in the main order "
                f"{main_order.cro_order}, so all its files go into its own "
                "output set; check that this is the right order"
            )
        for gid in sorted(report.new_groups):
            near = main_by_normalized.get(_normalized(gid))
            if near is not None:
                routing.warnings.append(
                    f"GroupID {gid!r} in {order.cro_order} differs from "
                    f"{near!r} in {main_order.cro_order} only by '-'/'_' or "
                    "case; it is treated as a new GroupID"
                )
    return routing
