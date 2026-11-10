"""Build and launch one streaming split pipeline per file, several in parallel.

    s3io get --concurrency N SRC
    | pigz -dc                       (or rapidgzip -d -c -P T)
    | split -l LINES --numeric-suffixes=1 -a W
        --filter='s3io put DST/STEM.$FILE.fastq.gz --upload-concurrency M --compress pigz -c -p T -L'
        - part

Python touches bytes only at the two endpoints (ranged reads in, multipart
upload out); record counting is GNU split, inflation pigz or rapidgzip,
compression pigz. Sidecars are written here, after a file's pipeline exits 0.
"""

from __future__ import annotations

import concurrent.futures
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import fsspec

from . import s3io
from .plan import FilePlan, GroupPlan, Plan

BCP_DIR = Path(__file__).resolve().parents[1]
PROBE_NAME = ".fastq_chunker_probe"
MANIFEST_COLUMNS = (
    "filename",
    "chunk",
    "bytes",
    "md5",
    "status",
    "wall_seconds",
    "checked",
)
# what a copied file's check compared: its size, and S3's full-object
# CRC64NVME against DB2's crc64nvme_base64 when both exist
CHECKED_SIZE = "size"
CHECKED_CRC = "size+crc64nvme"


class RunError(Exception):
    pass


DECOMPRESSORS = ("pigz", "rapidgzip")


@dataclass
class Tools:
    pigz: str
    split: str
    python: str = sys.executable
    rapidgzip: str | None = None


@dataclass
class ChunkResult:
    filename: str
    chunk: str
    bytes: int | None
    md5: str | None
    status: str
    wall_seconds: float | None = None
    checked: str = ""
    detail: str = ""  # why a file failed; printed, not a manifest column

    def row(self) -> dict:
        return {
            "filename": self.filename,
            "chunk": self.chunk,
            "bytes": "" if self.bytes is None else self.bytes,
            "md5": self.md5 or "",
            "status": self.status,
            "wall_seconds": ""
            if self.wall_seconds is None
            else f"{self.wall_seconds:.1f}",
            "checked": self.checked,
        }


@dataclass
class RunOptions:
    workers: int = 4
    pigz_threads: int | None = None
    gzip_level: int = 6
    only: list[str] = field(default_factory=list)
    force: bool = False
    dry_run: bool = False
    copy_singletons: bool = False
    copy_workers: int = 4
    copy_parts: int = s3io.DEFAULT_COPY_PARTS
    log_dir: Path = Path("logs")
    manifest: Path = Path("run_manifest.tsv")
    decompressor: str = "pigz"
    read_concurrency: int = s3io.DEFAULT_CONCURRENCY
    upload_concurrency: int = s3io.DEFAULT_UPLOAD_CONCURRENCY

    def threads(self) -> int:
        if self.pigz_threads:
            return self.pigz_threads
        # one core per pipeline is taken by the single-threaded pigz -dc
        return max(2, (os.cpu_count() or 2) // self.workers - 1)


def find_split(explicit: str | None = None) -> str:
    for name in [explicit] if explicit else ["split", "gsplit"]:
        path = shutil.which(name)
        if not path:
            continue
        out = subprocess.run([path, "--version"], capture_output=True, text=True)
        if out.returncode == 0 and "GNU coreutils" in out.stdout:
            return path
    raise RunError(
        "GNU split (with --filter) not found; on macOS `brew install coreutils` "
        "provides it as gsplit"
    )


def find_pigz(explicit: str | None = None) -> str:
    path = shutil.which(explicit or "pigz")
    if not path:
        raise RunError("pigz not found on PATH")
    return path


def find_rapidgzip(explicit: str | None = None) -> str:
    path = shutil.which(explicit or "rapidgzip")
    if not path:
        raise RunError(
            "rapidgzip not found on PATH; `pip install rapidgzip` provides it"
        )
    return path


def find_tools(
    pigz: str | None = None,
    split: str | None = None,
    decompressor: str = "pigz",
    rapidgzip: str | None = None,
) -> Tools:
    if decompressor not in DECOMPRESSORS:
        raise RunError(
            f"unknown decompressor {decompressor!r}; use {' or '.join(DECOMPRESSORS)}"
        )
    return Tools(
        pigz=find_pigz(pigz),
        split=find_split(split),
        rapidgzip=find_rapidgzip(rapidgzip) if decompressor == "rapidgzip" else None,
    )


def decompress_command(tools: Tools, decompressor: str, threads: int) -> str:
    q = shlex.quote
    if decompressor == "rapidgzip":
        if tools.rapidgzip is None:
            raise RunError(
                "rapidgzip requested but not located; call find_tools with it"
            )
        return f"{q(tools.rapidgzip)} -d -c -P {threads}"
    return f"{q(tools.pigz)} -dc"


def child_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in (str(BCP_DIR), env.get("PYTHONPATH")) if p
    )
    return env


def pipeline_command(
    tools: Tools,
    group: GroupPlan,
    fp: FilePlan,
    dst: str,
    threads: int,
    level: int,
    decompressor: str = "pigz",
    read_concurrency: int = s3io.DEFAULT_CONCURRENCY,
    upload_concurrency: int = s3io.DEFAULT_UPLOAD_CONCURRENCY,
) -> str:
    q = shlex.quote
    py = q(tools.python)
    get = (
        f"{py} -m fastq_chunker.s3io get --concurrency {read_concurrency} "
        f"{q(fp.file.s3_uri)}"
    )
    inflate = decompress_command(tools, decompressor, threads)
    # $FILE is expanded by split's filter shell, not by ours: split sets it to the
    # output name, prefix included, so ``part`` + the numeric suffix
    # dst and stem are user and portal input: single-quote them so the filter
    # shell expands nothing but $FILE, whatever characters they contain
    put_url = q(f"{dst}{fp.stem}.") + '"$FILE"' + q(".fastq.gz")
    # put runs pigz itself: split's filter shell has no pipefail, so a shell pipe
    # there would turn a dead compressor into a truncated chunk with exit 0
    filt = (
        f"{py} -m fastq_chunker.s3io put {put_url} "
        f"--upload-concurrency {upload_concurrency} "
        f"--compress {q(tools.pigz)} -c -p {threads} -{level}"
    )
    split = (
        f"{q(tools.split)} -l {group.lines_per_chunk} --numeric-suffixes=1 "
        f"-a {group.suffix_width} --filter={q(filt)} - part"
    )
    return f"set -o pipefail\n{get} | {inflate} | {split}"


def chunk_listing(
    fs: fsspec.AbstractFileSystem, dst: str, fp: FilePlan
) -> dict[str, int]:
    """Existing ``<stem>.part*.fastq.gz`` objects and sidecars under dst, with sizes."""
    fs.invalidate_cache()
    if not fs.exists(dst.rstrip("/")):
        return {}
    # anchored, so a stem that is a prefix of another stem (S1 and S1.part2)
    # cannot pick up the other file's chunks and delete them as leftovers
    pattern = chunk_pattern(fp.stem)
    out = {}
    for info in fs.ls(dst.rstrip("/"), detail=True):
        name = s3io.basename(info["name"])
        if pattern.match(name):
            out[name] = info.get("size") or 0
    return out


def chunk_pattern(stem: str) -> re.Pattern[str]:
    return re.compile(rf"^{re.escape(stem)}\.part\d+\.fastq\.gz(\.md5)?$")


def chunks_complete(fs, plan: Plan, fp: FilePlan) -> bool:
    have = chunk_listing(fs, plan.dst, fp)
    want = {c.name for c in fp.chunks}
    objects = {n for n in have if not n.endswith(".md5")}
    if objects != want:
        return False
    for c in fp.chunks:
        if have[c.name] >= plan.limit_bytes or c.name + ".md5" not in have:
            return False
    return True


def delete_chunks(fs, plan: Plan, fp: FilePlan) -> None:
    for name in chunk_listing(fs, plan.dst, fp):
        fs.rm(plan.dst + name)


def preflight(
    fs, plan: Plan, opts: RunOptions, files: list[tuple[GroupPlan, FilePlan]]
) -> None:
    problems = []
    sources = [fp for _, fp in files]
    if opts.copy_singletons:
        sources += singleton_files(plan, opts)
    for fp in sources:
        uri = fp.file.s3_uri
        if not fs.exists(uri):
            problems.append(f"missing: {uri}")
            continue
        size = fs.size(uri)
        if size != fp.file.size_bytes:
            problems.append(
                f"size mismatch (plan {fp.file.size_bytes}, object {size}); the plan "
                f"is stale: {uri}"
            )
    if problems:
        raise RunError("\n".join(problems))
    probe = plan.dst + PROBE_NAME
    try:
        fs.pipe(probe, b"ok")
        fs.rm(probe)
    except Exception as e:
        raise RunError(f"cannot write under {plan.dst}: {e}") from e


def parse_put_lines(stdout: str) -> list[tuple[str, int, str]]:
    out = []
    for line in stdout.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 3 and parts[1].isdigit():
            out.append((parts[0], int(parts[1]), parts[2]))
    return out


def run_file(
    tools: Tools, plan: Plan, group: GroupPlan, fp: FilePlan, opts: RunOptions, fs
) -> list[ChunkResult]:
    cmd = pipeline_command(
        tools,
        group,
        fp,
        plan.dst,
        opts.threads(),
        opts.gzip_level,
        opts.decompressor,
        opts.read_concurrency,
        opts.upload_concurrency,
    )
    opts.log_dir.mkdir(parents=True, exist_ok=True)
    log_path = opts.log_dir / f"{fp.stem}.log"
    t0 = time.monotonic()
    with open(log_path, "a") as log:
        log.write(f"# {time.strftime('%Y-%m-%dT%H:%M:%S')} start\n{cmd}\n")
        log.flush()
        proc = subprocess.run(
            ["bash", "-c", cmd],
            stdout=subprocess.PIPE,
            stderr=log,
            env=child_env(),
            text=True,
        )
        log.write(proc.stdout)
        log.write(f"# exit {proc.returncode} after {time.monotonic() - t0:.1f}s\n")
    wall = time.monotonic() - t0
    if proc.returncode != 0:
        delete_chunks(fs, plan, fp)
        return [
            ChunkResult(fp.file.filename, c.name, None, None, "failed", wall)
            for c in fp.chunks
        ]
    uploaded = {name: (n, md5) for name, n, md5 in parse_put_lines(proc.stdout)}
    results = []
    for c in fp.chunks:
        if c.name in uploaded:
            n, md5 = uploaded[c.name]
            results.append(ChunkResult(fp.file.filename, c.name, n, md5, "ok", wall))
        else:
            results.append(
                ChunkResult(fp.file.filename, c.name, None, None, "missing", wall)
            )
    extra = set(uploaded) - {c.name for c in fp.chunks}
    for name in sorted(extra):
        n, md5 = uploaded[name]
        results.append(ChunkResult(fp.file.filename, name, n, md5, "unplanned", wall))
    if any(r.status != "ok" for r in results):
        delete_chunks(fs, plan, fp)
        return results
    # a sidecar is written only once the whole file succeeded, so resume can
    # trust "chunk plus sidecar" to mean a chunk from a completed pipeline
    for r in results:
        s3io.write_sidecar(fs, plan.dst + r.chunk, r.md5)
    return results


def skipped_results(fs, plan: Plan, fp: FilePlan) -> list[ChunkResult]:
    have = chunk_listing(fs, plan.dst, fp)
    out = []
    for c in fp.chunks:
        md5, _ = s3io.read_sidecar(fs, plan.dst + c.name)
        out.append(ChunkResult(fp.file.filename, c.name, have[c.name], md5, "skipped"))
    return out


def write_manifest(path: Path, results: list[ChunkResult]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in results:
            w.writerow(r.row())


def singleton_files(plan: Plan, opts: RunOptions) -> list[FilePlan]:
    """Files that need no split, filtered by --only: copied as they are."""
    return [
        fp
        for g in plan.groups
        if g.action == "skip"
        for fp in g.files
        if not opts.only or fp.file.filename in opts.only
    ]


def check_copy(fs, target: str, fp: FilePlan, s3_client) -> tuple[str | None, str]:
    """(problem or None, what was compared) for a copied file at ``target``.

    The size must match the plan. When the destination is S3, DB2 recorded a
    CRC64NVME, and S3 holds a full-object CRC64NVME for the copy, the two must
    be equal; otherwise only the size was compared, and ``checked`` says so.
    """
    fs.invalidate_cache()
    if not fs.exists(target):
        return "not found after copy", CHECKED_SIZE
    size = fs.size(target)
    if size != fp.file.size_bytes:
        return f"size {size}, plan {fp.file.size_bytes}", CHECKED_SIZE
    if s3_client is None or not fp.file.crc64nvme:
        return None, CHECKED_SIZE
    crc, kind = s3io.head_crc64nvme(s3_client, target)
    if not crc or kind != "FULL_OBJECT":
        return None, CHECKED_SIZE
    if crc != fp.file.crc64nvme:
        return f"CRC64NVME {crc}, DB2 {fp.file.crc64nvme}", CHECKED_CRC
    return None, CHECKED_CRC


def copy_one(fs, plan: Plan, fp: FilePlan, opts: RunOptions, s3_client) -> ChunkResult:
    name = fp.file.filename
    target = plan.dst + name
    t0 = time.monotonic()

    def result(status: str, checked: str, detail: str = "") -> ChunkResult:
        return ChunkResult(
            name,
            name,
            fp.file.size_bytes if status in ("ok", "skipped") else None,
            None,
            status,
            time.monotonic() - t0,
            checked,
            detail,
        )

    try:
        if not opts.force and fs.exists(target):
            problem, checked = check_copy(fs, target, fp, s3_client)
            if problem is None:
                return result("skipped", checked)
        if s3_client is not None and fp.file.s3_uri.startswith("s3://"):
            s3io.managed_copy(
                s3_client, fp.file.s3_uri, target, fp.file.size_bytes, opts.copy_parts
            )
        else:
            fs.copy(fp.file.s3_uri, target)
        problem, checked = check_copy(fs, target, fp, s3_client)
    except Exception as e:  # one bad copy must not stop the others
        problem, checked = f"{type(e).__name__}: {e}", ""
    if problem is None:
        return result("ok", checked)
    try:
        if fs.exists(target):
            fs.rm(target)
    except Exception:
        pass
    return result("failed", checked, problem)


def copy_singletons(fs, plan: Plan, opts: RunOptions, out) -> list[ChunkResult]:
    """Server-side copy of every unsplit file into dst, several at once.

    ``copy_workers`` files at a time, each as ``copy_parts`` parts in flight,
    over one client whose pool holds a connection for every part. A file
    already in dst that passes the check is skipped; one that fails it is
    copied again. A failed copy is removed and reported; the rest go on.
    """
    todo = singleton_files(plan, opts)
    s3_client = (
        s3io.make_s3_client(opts.copy_workers * opts.copy_parts)
        if s3io.is_s3(fs)
        else None
    )
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=opts.copy_workers) as pool:
        futures = [pool.submit(copy_one, fs, plan, fp, opts, s3_client) for fp in todo]
        for fut in concurrent.futures.as_completed(futures):
            r = fut.result()
            results.append(r)
            if r.status == "failed":
                out.write(f"FAILED copy {r.filename}: {r.detail}\n")
            elif r.status == "skipped":
                out.write(f"present, skipping: {r.filename} ({r.checked})\n")
            else:
                out.write(f"copied {r.filename} ({r.checked})\n")
    return results


def run_plan(plan: Plan, opts: RunOptions, tools: Tools | None = None, out=None) -> int:
    out = out or sys.stdout
    if opts.copy_singletons and min(opts.copy_workers, opts.copy_parts) < 1:
        raise RunError("--copy-workers and --copy-parts must be at least 1")
    files = plan.split_files()
    if opts.only:
        known = {fp.file.filename for _, fp in files} | {
            fp.file.filename for g in plan.groups for fp in g.files
        }
        unknown = set(opts.only) - known
        if unknown:
            raise RunError(f"--only names not in plan: {', '.join(sorted(unknown))}")
        files = [(g, fp) for g, fp in files if fp.file.filename in opts.only]
    # largest first, so the long pole starts before the small index reads
    files.sort(key=lambda gf: gf[1].file.size_bytes, reverse=True)

    if files:
        # a plan with nothing to split only copies, and needs no pigz or split
        tools = tools or find_tools(decompressor=opts.decompressor)
    if opts.dry_run:
        for g, fp in files:
            out.write(
                pipeline_command(
                    tools,
                    g,
                    fp,
                    plan.dst,
                    opts.threads(),
                    opts.gzip_level,
                    opts.decompressor,
                    opts.read_concurrency,
                    opts.upload_concurrency,
                )
                + "\n\n"
            )
        if opts.copy_singletons:
            for fp in singleton_files(plan, opts):
                out.write(f"copy {fp.file.s3_uri} {plan.dst}{fp.file.filename}\n")
        return 0

    fs = s3io.make_fs(plan.dst)
    preflight(fs, plan, opts, files)

    results: list[ChunkResult] = []
    failed: list[str] = []
    if opts.copy_singletons:
        copies = copy_singletons(fs, plan, opts, out)
        results.extend(copies)
        failed.extend(r.filename for r in copies if r.status == "failed")

    todo = []
    for g, fp in files:
        if not opts.force and fp.chunks and chunks_complete(fs, plan, fp):
            out.write(f"complete, skipping: {fp.file.filename}\n")
            results.extend(skipped_results(fs, plan, fp))
            continue
        delete_chunks(fs, plan, fp)
        todo.append((g, fp))

    with concurrent.futures.ThreadPoolExecutor(max_workers=opts.workers) as pool:
        futures = {
            pool.submit(run_file, tools, plan, g, fp, opts, fs): fp for g, fp in todo
        }
        for fut in concurrent.futures.as_completed(futures):
            fp = futures[fut]
            res = fut.result()
            results.extend(res)
            bad = [r for r in res if r.status != "ok"]
            if bad:
                failed.append(fp.file.filename)
                out.write(
                    f"FAILED {fp.file.filename}: {bad[0].status} "
                    f"(see {opts.log_dir / (fp.stem + '.log')})\n"
                )
            else:
                out.write(
                    f"done {fp.file.filename}: {len(res)} chunks in {res[0].wall_seconds:.0f}s\n"
                )

    results.sort(key=lambda r: (r.filename, r.chunk))
    write_manifest(opts.manifest, results)
    out.write(f"run manifest written to {opts.manifest}\n")
    if failed:
        only = " ".join(shlex.quote(f) for f in sorted(failed))
        out.write(f"{len(failed)} file(s) failed; rerun with: --only {only}\n")
        return 1
    return 0
