"""Pack chunked sets into SRA submission batches under a total-size limit.

Organisational only: nothing in the chunking path depends on batches. Every
chunk of a set lands in one batch, sets stay in plan order, and packing is
first-fit-in-order rather than optimal so batches are predictable.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .plan import GB, Plan, PlanError

BATCH_LIMIT_BYTES = 5_000 * GB
COLUMNS = ("batch", "group", "label", "name", "bytes", "uri")


@dataclass(frozen=True)
class Item:
    group: str
    label: str
    name: str
    bytes: int
    uri: str


@dataclass
class Batch:
    number: int
    items: list[Item]

    @property
    def bytes(self) -> int:
        return sum(i.bytes for i in self.items)


def actual_sizes(run_manifest: Path | None) -> dict[str, int]:
    if run_manifest is None:
        return {}
    sizes = {}
    with open(run_manifest, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["bytes"]:
                sizes[row["chunk"]] = int(row["bytes"])
    return sizes


def group_items(plan: Plan, sizes: dict[str, int]) -> list[list[Item]]:
    out = []
    for g in plan.groups:
        items = []
        for fp in g.files:
            if g.action == "split":
                for c in fp.chunks:
                    items.append(
                        Item(
                            g.group,
                            g.label,
                            c.name,
                            sizes.get(c.name, c.est_bytes),
                            fp.chunk_uri(plan.dst, c),
                        )
                    )
            else:
                items.append(
                    Item(
                        g.group,
                        g.label,
                        fp.file.filename,
                        fp.file.size_bytes,
                        fp.file.s3_uri,
                    )
                )
        out.append(items)
    return out


def build_batches(
    plan: Plan, run_manifest: Path | None = None, limit_bytes: int = BATCH_LIMIT_BYTES
) -> list[Batch]:
    batches: list[Batch] = []
    current = Batch(1, [])
    for items in group_items(plan, actual_sizes(run_manifest)):
        total = sum(i.bytes for i in items)
        if total > limit_bytes:
            raise PlanError(
                f"set {items[0].label} alone is {total / GB:.0f} GB, above the batch limit "
                f"{limit_bytes / GB:.0f} GB"
            )
        if current.items and current.bytes + total > limit_bytes:
            batches.append(current)
            current = Batch(current.number + 1, [])
        current.items.extend(items)
    if current.items:
        batches.append(current)
    return batches


def write_batches(path: Path, batches: list[Batch]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(COLUMNS)
        for b in batches:
            for i in b.items:
                w.writerow([b.number, i.group, i.label, i.name, i.bytes, i.uri])
