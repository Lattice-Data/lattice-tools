"""Command line entry point: plan | run | verify | batch | test-data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import batch as batch_mod
from . import portal, s3io, testdata, verify
from .plan import (
    GB,
    LIMIT_BYTES,
    ROUND_TO,
    TARGET_BYTES,
    Plan,
    PlanError,
    load_plan,
    make_plan,
    plan_group,
    write_plan,
)
from .run import RunError, RunOptions, find_tools, run_plan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastq_chunker",
        description=(
            "Split large gzipped FASTQs in S3 into record-aligned chunks below "
            "the SRA per-file limit, writing the chunks to a second S3 prefix."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser(
        "plan",
        help="compute chunk counts and names from portal JSON; touches no network",
    )
    p.add_argument(
        "--file-sets",
        required=True,
        type=Path,
        help="sequence_file_set JSON: one object, a list, or a search @graph",
    )
    p.add_argument(
        "--files",
        required=True,
        type=Path,
        help="sequence_file JSON for every file the sets reference",
    )
    p.add_argument(
        "--dst", required=True, help="destination prefix, s3://bucket/prefix/"
    )
    p.add_argument("--out", type=Path, default=Path("plan.json"))
    size = p.add_mutually_exclusive_group()
    size.add_argument("--target-gb", type=float, default=TARGET_BYTES / GB)
    size.add_argument("--target-bytes", type=int, help="exact target, for tests")
    limit = p.add_mutually_exclusive_group()
    limit.add_argument("--limit-gb", type=float, default=LIMIT_BYTES / GB)
    limit.add_argument("--limit-bytes", type=int, help="exact limit, for tests")
    p.add_argument("--round-to", type=int, default=ROUND_TO)
    p.add_argument(
        "--names-only",
        action="store_true",
        help="print only the expected chunk names, one per line",
    )
    p.set_defaults(func=cmd_plan)

    r = sub.add_parser("run", help="split every planned file, several in parallel")
    r.add_argument("--plan", required=True, type=Path)
    r.add_argument("--workers", type=int, default=4, help="pipelines at once")
    r.add_argument(
        "--pigz-threads",
        type=int,
        help="compression threads per pipeline (default: cores/workers - 1)",
    )
    r.add_argument("--gzip-level", type=int, default=6, choices=range(1, 10))
    r.add_argument("--only", nargs="+", default=[], metavar="FILENAME")
    r.add_argument("--force", action="store_true", help="re-split complete files too")
    r.add_argument("--dry-run", action="store_true", help="print the bash per file")
    r.add_argument(
        "--copy-singletons",
        action="store_true",
        help="server-side copy files that need no splitting into --dst as well",
    )
    r.add_argument("--log-dir", type=Path, default=Path("logs"))
    r.add_argument("--manifest", type=Path, default=Path("run_manifest.tsv"))
    r.add_argument("--pigz", help="path to pigz (default: from PATH)")
    r.add_argument("--split", help="path to GNU split (default: split, then gsplit)")
    r.set_defaults(func=cmd_run)

    v = sub.add_parser("verify", help="check the chunks under the plan's destination")
    v.add_argument("--plan", required=True, type=Path)
    v.add_argument("--level", choices=("quick", "full"), default="quick")
    v.add_argument("--workers", type=int, default=4)
    v.add_argument("--out", type=Path, default=Path("verify_report.tsv"))
    v.add_argument(
        "--id-regex",
        help="regex with one capture group that extracts the shared read ID from a header",
    )
    v.add_argument("--pigz", help="path to pigz (full level only)")
    v.add_argument("--split", help="path to GNU split (unused, accepted for symmetry)")
    v.set_defaults(func=cmd_verify)

    b = sub.add_parser(
        "batch", help="pack sets into submission batches under a size limit"
    )
    b.add_argument("--plan", required=True, type=Path)
    b.add_argument(
        "--run-manifest",
        type=Path,
        help="run_manifest.tsv for actual chunk sizes (default: plan estimates)",
    )
    b.add_argument(
        "--batch-limit-gb", type=float, default=batch_mod.BATCH_LIMIT_BYTES / GB
    )
    b.add_argument("--out", type=Path, default=Path("batches.tsv"))
    b.set_defaults(func=cmd_batch)

    t = sub.add_parser(
        "test-data", help="write a small synthetic quadruple plus portal JSON"
    )
    t.add_argument("--reads", type=int, default=20_000)
    t.add_argument(
        "--out-prefix", required=True, help="where the FASTQs go (any fsspec URL)"
    )
    t.add_argument(
        "--meta-dir",
        type=Path,
        default=Path("."),
        help="where sets.json and files.json go",
    )
    t.add_argument("--seed", type=int, default=353)
    t.set_defaults(func=cmd_test_data)
    return parser


def cmd_plan(args: argparse.Namespace) -> int:
    target = args.target_bytes or int(args.target_gb * GB)
    limit = args.limit_bytes or int(args.limit_gb * GB)
    sets = portal.load_objects(args.file_sets)
    files = portal.load_objects(args.files)
    label_of = portal.labels(sets)
    groups = [
        plan_group(
            members,
            target_bytes=target,
            limit_bytes=limit,
            round_to=args.round_to,
            label=label_of.get(members[0].group),
        )
        for members in portal.build_groups(sets, files)
    ]
    plan = make_plan(groups, args.dst, target, limit, args.round_to)
    write_plan(plan, args.out)
    if args.names_only:
        for g in plan.groups:
            for fp in g.files:
                for c in fp.chunks:
                    print(c.name)
        return 0
    print_plan_summary(plan, args.out, sys.stdout)
    return 0


def print_plan_summary(plan: Plan, out_path: Path, out) -> None:
    n_split = sum(1 for g in plan.groups if g.action == "split")
    out.write(
        f"{len(plan.groups)} groups, {n_split} to split; target "
        f"{plan.target_bytes / GB:.0f} GB, limit {plan.limit_bytes / GB:.0f} GB; "
        f"plan written to {out_path}\n"
    )
    for g in plan.groups:
        largest = g.largest_est_bytes()
        if g.action == "split":
            margin = 100 * (1 - largest / plan.limit_bytes)
            out.write(
                f"  {g.label}: {g.n_chunks} chunks x {g.reads_per_chunk:,} reads, "
                f"largest chunk ~{largest / GB:.1f} GB ({margin:.0f}% below limit)\n"
            )
            for fp in g.files:
                out.write(f"    {fp.file.filename}\n")
                for c in fp.chunks:
                    out.write(f"      {c.name}  ~{c.est_bytes / GB:.1f} GB\n")
        else:
            out.write(
                f"  {g.label}: skip, largest file {largest / GB:.1f} GB is below the limit\n"
            )


def cmd_run(args: argparse.Namespace) -> int:
    plan = load_plan(args.plan)
    opts = RunOptions(
        workers=args.workers,
        pigz_threads=args.pigz_threads,
        gzip_level=args.gzip_level,
        only=args.only,
        force=args.force,
        dry_run=args.dry_run,
        copy_singletons=args.copy_singletons,
        log_dir=args.log_dir,
        manifest=args.manifest,
    )
    return run_plan(plan, opts, tools=find_tools(args.pigz, args.split))


def cmd_verify(args: argparse.Namespace) -> int:
    plan = load_plan(args.plan)
    fs = s3io.make_fs(plan.dst)
    checks = verify.quick(plan, fs, args.id_regex)
    if args.level == "full":
        tools = find_tools(args.pigz, args.split)
        checks.extend(verify.full(plan, fs, tools, args.workers, args.id_regex))
    ok = verify.write_report(args.out, checks)
    fails = [c for c in checks if c.status == "FAIL"]
    warns = [c for c in checks if c.status == "WARN"]
    for c in fails[:50]:
        print(f"FAIL {c.filename} {c.chunk} {c.check}: {c.detail}")
    for c in warns[:50]:
        print(f"WARN {c.filename} {c.chunk} {c.check}: {c.detail}")
    print(
        f"{'PASS' if ok else 'FAIL'}: {len(checks)} checks, {len(fails)} failed, "
        f"{len(warns)} warnings; report written to {args.out}"
    )
    return 0 if ok else 1


def cmd_batch(args: argparse.Namespace) -> int:
    plan = load_plan(args.plan)
    batches = batch_mod.build_batches(
        plan, args.run_manifest, int(args.batch_limit_gb * GB)
    )
    batch_mod.write_batches(args.out, batches)
    for b in batches:
        groups = len({i.group for i in b.items})
        print(
            f"batch {b.number}: {groups} sets, {len(b.items)} files, {b.bytes / GB:.0f} GB"
        )
    print(f"batches written to {args.out}")
    return 0


def cmd_test_data(args: argparse.Namespace) -> int:
    sets_path, files_path = testdata.generate(
        args.reads, args.out_prefix, args.meta_dir, seed=args.seed
    )
    print(f"wrote {sets_path} and {files_path}; FASTQs under {args.out_prefix}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (PlanError, RunError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
