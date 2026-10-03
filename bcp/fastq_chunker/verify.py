"""Verify chunks: quick (seconds per chunk) and full (a decompressing pass)."""

from __future__ import annotations

import concurrent.futures
import csv
import re
import shlex
import subprocess
import sys
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from . import s3io
from .plan import ChunkSpec, FilePlan, GroupPlan, Plan
from .run import child_env, chunk_listing

HEADER_BYTES = 64 * 1024
# Recompression at pigz -6 came out 15-22% smaller than the vendor stream on
# real data, so undershoot is expected and only a large one is worth a look;
# overshoot is what could threaten the limit.
DRIFT_OVER = 0.15
DRIFT_UNDER = 0.30
REPORT_COLUMNS = ("group", "filename", "chunk", "check", "status", "detail")


@dataclass
class Check:
    group: str
    filename: str
    chunk: str
    check: str
    status: str  # PASS, WARN, FAIL
    detail: str = ""

    def row(self) -> dict:
        return dict(
            zip(
                REPORT_COLUMNS,
                (
                    self.group,
                    self.filename,
                    self.chunk,
                    self.check,
                    self.status,
                    self.detail,
                ),
            )
        )


def read_id(header: str, id_regex: str | None = None) -> str:
    """The part of a FASTQ header that mates share.

    Default rule: up to the first whitespace, then a trailing ``/1`` or ``/2``
    removed. Illumina and Ultima both fit. ``id_regex`` (one capture group)
    replaces the rule for vendors that put the mate marker elsewhere.
    """
    header = header.rstrip("\r\n")
    if id_regex:
        m = re.search(id_regex, header)
        if not m:
            raise ValueError(f"id regex {id_regex!r} does not match {header!r}")
        return m.group(1)
    rid = header.split(None, 1)[0] if header.split() else header
    return re.sub(r"/[12]$", "", rid)


def first_record(fs, url: str) -> list[str]:
    """First four decompressed lines of a gzip object, reading as little as needed."""
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    text = b""
    offset = 0
    size = fs.size(url)
    with fs.open(url, "rb", block_size=HEADER_BYTES) as f:
        while text.count(b"\n") < 4 and offset < size:
            raw = f.read(HEADER_BYTES)
            if not raw:
                break
            offset += len(raw)
            text += d.decompress(raw)
    lines = text.decode("latin-1").split("\n")
    return lines[:4]


def header_sanity(lines: list[str]) -> str | None:
    if len(lines) < 4 or not all(lines[:4]):
        return "fewer than four lines"
    if not lines[0].startswith("@"):
        return f"line 1 does not start with @: {lines[0][:40]!r}"
    if not lines[2].startswith("+"):
        return f"line 3 does not start with +: {lines[2][:40]!r}"
    if len(lines[1]) != len(lines[3]):
        return "sequence and quality lengths differ"
    return None


def drift_detail(size: int, est_bytes: int) -> str | None:
    """A WARN detail when a chunk's size is far from its estimate, else None."""
    if not est_bytes:
        return None
    ratio = size / est_bytes - 1
    if ratio > DRIFT_OVER or ratio < -DRIFT_UNDER:
        return f"{size} is {ratio:+.0%} vs estimate {est_bytes}"
    return None


def quick(plan: Plan, fs, id_regex: str | None = None) -> list[Check]:
    checks: list[Check] = []
    for g in plan.groups:
        if g.action != "split":
            continue
        ids: dict[int, dict[str, str]] = defaultdict(dict)
        for fp in g.files:
            checks.extend(_quick_file(plan, fs, g, fp, ids, id_regex))
        for k in range(1, g.n_chunks + 1):
            seen = ids.get(k, {})
            if len(seen) < len(g.files):
                continue  # a missing chunk was already reported
            distinct = set(seen.values())
            if len(distinct) == 1:
                checks.append(
                    Check(
                        g.group,
                        "*",
                        f"part{k:0{g.suffix_width}d}",
                        "mate_alignment",
                        "PASS",
                        next(iter(distinct)),
                    )
                )
            else:
                detail = "; ".join(f"{f}={rid}" for f, rid in sorted(seen.items()))
                checks.append(
                    Check(
                        g.group,
                        "*",
                        f"part{k:0{g.suffix_width}d}",
                        "mate_alignment",
                        "FAIL",
                        detail,
                    )
                )
    return checks


def _quick_file(plan, fs, g: GroupPlan, fp: FilePlan, ids, id_regex) -> list[Check]:
    out = []
    have = chunk_listing(fs, plan.dst, fp)
    objects = {n for n in have if not n.endswith(".md5")}
    planned = {c.name for c in fp.chunks}
    extra = sorted(objects - planned)
    out.append(
        Check(
            g.group,
            fp.file.filename,
            "*",
            "chunk_count",
            "FAIL" if extra else "PASS",
            f"stray: {', '.join(extra)}" if extra else f"{len(planned)} chunks",
        )
    )
    for c in fp.chunks:
        if c.name not in have:
            out.append(
                Check(g.group, fp.file.filename, c.name, "exists", "FAIL", "missing")
            )
            continue
        size = have[c.name]
        out.append(
            Check(g.group, fp.file.filename, c.name, "exists", "PASS", str(size))
        )
        if size >= plan.limit_bytes:
            out.append(
                Check(
                    g.group,
                    fp.file.filename,
                    c.name,
                    "size_below_limit",
                    "FAIL",
                    str(size),
                )
            )
        else:
            out.append(
                Check(
                    g.group,
                    fp.file.filename,
                    c.name,
                    "size_below_limit",
                    "PASS",
                    str(size),
                )
            )
        drift = drift_detail(size, c.est_bytes)
        if drift:
            out.append(
                Check(g.group, fp.file.filename, c.name, "size_drift", "WARN", drift)
            )
        sidecar = c.name + ".md5"
        if sidecar not in have:
            out.append(
                Check(g.group, fp.file.filename, c.name, "sidecar", "FAIL", "missing")
            )
        else:
            _, name = s3io.read_sidecar(fs, plan.dst + c.name)
            ok = name == c.name
            out.append(
                Check(
                    g.group,
                    fp.file.filename,
                    c.name,
                    "sidecar",
                    "PASS" if ok else "FAIL",
                    "" if ok else f"names {name}",
                )
            )
        lines = first_record(fs, plan.dst + c.name)
        problem = header_sanity(lines)
        if problem:
            out.append(
                Check(g.group, fp.file.filename, c.name, "header", "FAIL", problem)
            )
            continue
        try:
            rid = read_id(lines[0], id_regex)
        except ValueError as e:
            out.append(
                Check(g.group, fp.file.filename, c.name, "header", "FAIL", str(e))
            )
            continue
        ids[c.index][fp.file.filename] = rid
        out.append(Check(g.group, fp.file.filename, c.name, "header", "PASS", rid))
    return out


# ---- full ---------------------------------------------------------------------


def count_stream(src, block: int) -> tuple[int, str, str]:
    """Count newlines; return (lines, first line, header of the last record).

    The last record's header is the fourth line from the end of a well-formed
    chunk; when the line count is not a multiple of four that choice is
    arbitrary, but the ``lines_mod_4`` check fails such a chunk anyway.
    """
    lines = 0
    first: bytes | None = None
    tail = b""
    last4: list[bytes] = []
    while True:
        b = src.read(block)
        if not b:
            break
        lines += b.count(b"\n")
        parts = (tail + b).split(b"\n")
        tail = parts.pop()
        if parts:
            if first is None:
                first = parts[0]
            last4 = (last4 + parts)[-4:]
    if first is None:
        first = tail
    last_header = last4[-4] if len(last4) == 4 else b""
    return lines, first.decode("latin-1"), last_header.decode("latin-1")


def count_main() -> int:
    lines, first, last = count_stream(sys.stdin.buffer, s3io.block_size())
    print(f"{lines}\t{first}\t{last}")
    return 0


@dataclass
class FullResult:
    lines: int
    first: str
    last: str
    md5: str


def full_one(pigz: str, url: str) -> FullResult:
    q = shlex.quote
    py = q(sys.executable)
    cmd = (
        "set -o pipefail\n"
        f"{py} -m fastq_chunker.s3io get --md5 {q(url)} | {q(pigz)} -dc | "
        f"{py} -m fastq_chunker.verify count"
    )
    proc = subprocess.run(
        ["bash", "-c", cmd], capture_output=True, text=True, env=child_env()
    )
    if proc.returncode != 0:
        raise RuntimeError(f"verify pipeline failed for {url}:\n{proc.stderr[-2000:]}")
    lines, first, last = proc.stdout.rstrip("\n").split("\t")
    md5 = ""
    for line in proc.stderr.splitlines():
        if line.startswith("md5\t"):
            md5 = line.split("\t", 1)[1].strip()
    return FullResult(int(lines), first, last, md5)


def full(
    plan: Plan, fs, pigz: str, workers: int = 4, id_regex: str | None = None
) -> list[Check]:
    jobs: list[tuple[GroupPlan, FilePlan, ChunkSpec]] = [
        (g, fp, c)
        for g in plan.groups
        if g.action == "split"
        for fp in g.files
        for c in fp.chunks
    ]
    results: dict[tuple[str, str], FullResult | Exception] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(full_one, pigz, plan.dst + c.name): (fp, c) for _, fp, c in jobs
        }
        for fut in concurrent.futures.as_completed(futs):
            fp, c = futs[fut]
            try:
                results[(fp.file.filename, c.name)] = fut.result()
            except Exception as e:
                results[(fp.file.filename, c.name)] = e

    checks: list[Check] = []
    for g in plan.groups:
        if g.action != "split":
            continue
        firsts: dict[int, dict[str, str]] = defaultdict(dict)
        lasts: dict[int, dict[str, str]] = defaultdict(dict)
        for fp in g.files:
            for c in fp.chunks:
                r = results[(fp.file.filename, c.name)]
                if isinstance(r, Exception):
                    checks.append(
                        Check(
                            g.group,
                            fp.file.filename,
                            c.name,
                            "pipeline",
                            "FAIL",
                            str(r)[:200],
                        )
                    )
                    continue
                ok = r.lines % 4 == 0
                checks.append(
                    Check(
                        g.group,
                        fp.file.filename,
                        c.name,
                        "lines_mod_4",
                        "PASS" if ok else "FAIL",
                        str(r.lines),
                    )
                )
                ok = r.lines // 4 == c.reads
                checks.append(
                    Check(
                        g.group,
                        fp.file.filename,
                        c.name,
                        "read_count",
                        "PASS" if ok else "FAIL",
                        f"{r.lines // 4} vs planned {c.reads}",
                    )
                )
                if fs.exists(s3io.sidecar_url(plan.dst + c.name)):
                    side, _ = s3io.read_sidecar(fs, plan.dst + c.name)
                    ok = side == r.md5
                    detail = r.md5 if ok else f"{r.md5} vs sidecar {side}"
                else:
                    ok = False
                    detail = f"{r.md5}; sidecar missing"
                checks.append(
                    Check(
                        g.group,
                        fp.file.filename,
                        c.name,
                        "md5",
                        "PASS" if ok else "FAIL",
                        detail,
                    )
                )
                try:
                    firsts[c.index][fp.file.filename] = read_id(r.first, id_regex)
                    lasts[c.index][fp.file.filename] = read_id(r.last, id_regex)
                except ValueError as e:
                    checks.append(
                        Check(
                            g.group, fp.file.filename, c.name, "header", "FAIL", str(e)
                        )
                    )
        for k in range(1, g.n_chunks + 1):
            label = f"part{k:0{g.suffix_width}d}"
            for name, table in (
                ("first_id_aligned", firsts),
                ("last_id_aligned", lasts),
            ):
                seen = table.get(k, {})
                if len(seen) < len(g.files):
                    continue
                distinct = set(seen.values())
                ok = len(distinct) == 1
                checks.append(
                    Check(
                        g.group,
                        "*",
                        label,
                        name,
                        "PASS" if ok else "FAIL",
                        next(iter(distinct))
                        if ok
                        else "; ".join(f"{f}={v}" for f, v in sorted(seen.items())),
                    )
                )
            if k < g.n_chunks:
                for fp in g.files:
                    last = lasts.get(k, {}).get(fp.file.filename)
                    nxt = firsts.get(k + 1, {}).get(fp.file.filename)
                    if last is None or nxt is None:
                        continue
                    ok = last != nxt
                    checks.append(
                        Check(
                            g.group,
                            fp.file.filename,
                            label,
                            "boundary_distinct",
                            "PASS" if ok else "FAIL",
                            "" if ok else f"{last} repeats at start of next chunk",
                        )
                    )
    return checks


def write_report(path: Path, checks: list[Check]) -> bool:
    failed = any(c.status == "FAIL" for c in checks)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=REPORT_COLUMNS, delimiter="\t")
        w.writeheader()
        for c in checks:
            w.writerow(c.row())
        w.writerow(
            dict(
                zip(
                    REPORT_COLUMNS,
                    ("*", "*", "*", "overall", "FAIL" if failed else "PASS", ""),
                )
            )
        )
    return not failed


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv == ["count"]:
        return count_main()
    print("usage: python -m fastq_chunker.verify count", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
