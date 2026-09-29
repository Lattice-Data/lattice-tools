"""Planning: chunk counts, reads per chunk, chunk names. Pure functions, no network."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

# Decimal gigabytes: SRA states 100 GB, and 100 GB < 100 GiB, so this is the
# conservative reading of the limit.
GB = 10**9
LIMIT_BYTES = 100 * GB
TARGET_BYTES = 80 * GB
ROUND_TO = 1_000_000
LINES_PER_READ = 4
MIN_SUFFIX_WIDTH = 3
EXTENSIONS = (".fastq.gz", ".fq.gz")


class PlanError(ValueError):
    pass


@dataclass(frozen=True)
class FastqFile:
    group: str
    slot: str
    filename: str
    s3_uri: str
    size_bytes: int
    reads: int


@dataclass(frozen=True)
class ChunkSpec:
    index: int
    name: str
    reads: int
    est_bytes: int


@dataclass
class FilePlan:
    file: FastqFile
    stem: str
    chunks: list[ChunkSpec]

    def chunk_uri(self, dst: str, chunk: ChunkSpec) -> str:
        return dst.rstrip("/") + "/" + chunk.name


@dataclass
class GroupPlan:
    group: str
    label: str
    action: str
    n_chunks: int
    reads_per_chunk: int
    lines_per_chunk: int
    suffix_width: int
    files: list[FilePlan]

    @property
    def reads(self) -> int:
        return self.files[0].file.reads

    def largest_est_bytes(self) -> int:
        if self.action != "split":
            return max(fp.file.size_bytes for fp in self.files)
        return max(c.est_bytes for fp in self.files for c in fp.chunks)


@dataclass
class Plan:
    created: str
    dst: str
    target_bytes: int
    limit_bytes: int
    round_to: int
    groups: list[GroupPlan]

    def split_files(self) -> list[tuple[GroupPlan, FilePlan]]:
        return [(g, fp) for g in self.groups if g.action == "split" for fp in g.files]


def stem_of(filename: str) -> str:
    for ext in EXTENSIONS:
        if filename.endswith(ext) and len(filename) > len(ext):
            return filename[: -len(ext)]
    raise PlanError(
        f"unsupported extension (want {' or '.join(EXTENSIONS)}): {filename}"
    )


def chunk_name(stem: str, index: int, width: int) -> str:
    return f"{stem}.part{index:0{width}d}.fastq.gz"


def suffix_width(n_chunks: int) -> int:
    return max(MIN_SUFFIX_WIDTH, len(str(n_chunks)))


def chunk_count(biggest_bytes: int, target_bytes: int, limit_bytes: int) -> int:
    """A file already below the limit is legal as-is and is never split."""
    if biggest_bytes < limit_bytes:
        return 1
    return math.ceil(biggest_bytes / target_bytes)


def reads_per_chunk(total_reads: int, n_chunks: int, round_to: int) -> int:
    """Reads in every chunk but the last, rounded up to a multiple of ``round_to``.

    Rounding up must not empty the last chunk; when it would, the exact value is
    used instead.
    """
    if n_chunks == 1:
        return total_reads
    if n_chunks > total_reads:
        raise PlanError(
            f"{total_reads} reads cannot fill {n_chunks} chunks; the file is too "
            "large for its read count"
        )
    exact = math.ceil(total_reads / n_chunks)
    rounded = math.ceil(exact / round_to) * round_to
    n = rounded if (n_chunks - 1) * rounded < total_reads else exact
    fits = (n_chunks - 1) * n < total_reads <= n_chunks * n
    if not fits or math.ceil(total_reads / n) != n_chunks:
        raise PlanError(
            f"{total_reads} reads cannot be spread over {n_chunks} chunks of {n}; "
            "the file is too large for its read count"
        )
    return n


def plan_group(
    files: list[FastqFile],
    target_bytes: int = TARGET_BYTES,
    limit_bytes: int = LIMIT_BYTES,
    round_to: int = ROUND_TO,
    label: str | None = None,
) -> GroupPlan:
    if not files:
        raise PlanError("empty group")
    if target_bytes >= limit_bytes:
        raise PlanError(f"target {target_bytes} must be below limit {limit_bytes}")
    groups = {f.group for f in files}
    if len(groups) != 1:
        raise PlanError(f"files from more than one group: {sorted(groups)}")
    group = files[0].group
    counts = {f.reads for f in files}
    if len(counts) != 1:
        detail = ", ".join(f"{f.filename}={f.reads}" for f in files)
        raise PlanError(f"read counts differ within group {group}: {detail}")
    total_reads = files[0].reads
    if total_reads <= 0:
        raise PlanError(f"group {group} has no reads")

    biggest = max(files, key=lambda f: f.size_bytes)
    n_chunks = chunk_count(biggest.size_bytes, target_bytes, limit_bytes)
    n = reads_per_chunk(total_reads, n_chunks, round_to)
    if n_chunks > 1 and biggest.size_bytes * n / total_reads >= limit_bytes:
        est_gb = biggest.size_bytes * n / total_reads / GB
        raise PlanError(
            f"group {group}: estimated chunk {est_gb:.1f} GB is not below the "
            f"{limit_bytes / GB:.0f} GB limit; lower the target"
        )

    width = suffix_width(n_chunks)
    action = "split" if n_chunks > 1 else "skip"
    file_plans = []
    for f in files:
        stem = stem_of(f.filename)
        chunks = []
        if action == "split":
            for k in range(1, n_chunks + 1):
                reads_k = n if k < n_chunks else total_reads - (n_chunks - 1) * n
                chunks.append(
                    ChunkSpec(
                        index=k,
                        name=chunk_name(stem, k, width),
                        reads=reads_k,
                        est_bytes=round(f.size_bytes * reads_k / total_reads),
                    )
                )
        file_plans.append(FilePlan(file=f, stem=stem, chunks=chunks))

    return GroupPlan(
        group=group,
        label=label or group,
        action=action,
        n_chunks=n_chunks,
        reads_per_chunk=n,
        lines_per_chunk=LINES_PER_READ * n,
        suffix_width=width,
        files=file_plans,
    )


def make_plan(
    groups: list[GroupPlan],
    dst: str,
    target_bytes: int,
    limit_bytes: int,
    round_to: int,
) -> Plan:
    owners: dict[str, list[str]] = defaultdict(list)
    for g in groups:
        for fp in g.files:
            owners[fp.stem].append(fp.file.s3_uri)
    shared = {stem: uris for stem, uris in owners.items() if len(uris) > 1}
    if shared:
        detail = "\n".join(
            f"  {stem}: {', '.join(uris)}" for stem, uris in sorted(shared.items())
        )
        raise PlanError(
            "chunk names would collide under one destination prefix; these files "
            f"share a stem:\n{detail}"
        )
    created = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return Plan(
        created=created,
        dst=dst.rstrip("/") + "/",
        target_bytes=target_bytes,
        limit_bytes=limit_bytes,
        round_to=round_to,
        groups=groups,
    )


def plan_to_dict(plan: Plan) -> dict:
    return asdict(plan)


def plan_from_dict(d: dict) -> Plan:
    groups = []
    for g in d["groups"]:
        files = []
        for fp in g["files"]:
            files.append(
                FilePlan(
                    file=FastqFile(**fp["file"]),
                    stem=fp["stem"],
                    chunks=[ChunkSpec(**c) for c in fp["chunks"]],
                )
            )
        groups.append(GroupPlan(**{**g, "files": files}))
    return Plan(**{**d, "groups": groups})


def write_plan(plan: Plan, path: Path) -> None:
    path.write_text(json.dumps(plan_to_dict(plan), indent=2) + "\n")


def load_plan(path: Path) -> Plan:
    return plan_from_dict(json.loads(Path(path).read_text()))
