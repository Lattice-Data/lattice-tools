"""``check-submission``: is the destination exactly what the SRA sheet names?

The last gate before an SRA upload. Reads the destination itself rather than
trusting the run manifest, and changes nothing. Every one of these must hold:

- the plan is for this sheet: the names the plan produces (a chunk's name, or
  an unsplit file's own name) are exactly the sheet's ``filename`` cells;
- the destination holds exactly those objects, nothing missing and nothing
  extra; the only other objects allowed are the ``.md5`` sidecars of chunks,
  and the empty marker the S3 console leaves at the prefix itself;
- every unsplit file has its planned size and, when DB2 has a CRC64NVME
  for it, a full-object CRC64NVME in S3 equal to it (the ``run`` copy check);
- every chunk is below the limit, has its sidecar, and passed ``verify
  --level full`` (read count and md5), with no FAIL for it, its file's part,
  or its part across mates (the alignment checks);
- a verify report, when given, has no FAIL at all and ends in an overall PASS,
  and records for every chunk the size and ETag of the object it read, which
  must still be the object in the destination;
- the batch, in actual bytes, is within the batch limit.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from . import s3io
from .plan import GB, Plan
from .run import PROBE_NAME, check_copy

COLUMNS = ("name", "check", "status", "detail")
FULL_CHECKS = ("read_count", "md5")


@dataclass
class Row:
    name: str
    check: str
    status: str  # PASS or FAIL
    detail: str = ""


def expected_names(plan: Plan) -> dict[str, tuple]:
    """name -> (kind, group, file plan, chunk or None) for everything the plan produces."""
    out = {}
    for g in plan.groups:
        for fp in g.files:
            if g.action == "split":
                for c in fp.chunks:
                    out[c.name] = ("chunk", g, fp, c)
            else:
                out[fp.file.filename] = ("file", g, fp, None)
    return out


@dataclass
class VerifyReport:
    passed: dict[str, set[str]]  # chunk -> full-level checks that PASSed
    failed: set[str]  # chunks a FAIL row names, by itself or by its part
    fails: list[tuple[str, str, str]]  # every FAIL row: (filename, chunk, check)
    overall_pass: bool  # the report ends in an overall PASS
    objects: dict[str, str]  # chunk -> object_identity verify read


def chunk_labels(plan: Plan) -> dict[tuple[str, str, str], list[str]]:
    """(group, filename, label) -> chunks a verify row with that key is about.

    ``verify --level full`` writes the mate checks (first_id_aligned,
    last_id_aligned) with filename ``*`` and the part label (``part001``) in
    the chunk column, and boundary_distinct with the file and the part label.
    """
    out: dict[tuple[str, str, str], list[str]] = {}
    for g in plan.groups:
        if g.action != "split":
            continue
        for fp in g.files:
            for c in fp.chunks:
                label = f"part{c.index:0{g.suffix_width}d}"
                out.setdefault((g.group, fp.file.filename, label), []).append(c.name)
                out.setdefault((g.group, "*", label), []).append(c.name)
    return out


def read_verify_report(report: Path | None, plan: Plan) -> VerifyReport:
    out = VerifyReport({}, set(), [], False, {})
    if report is None:
        return out
    if not report.exists():
        # e.g. the gate run before verify finished: a FAIL, not a traceback
        out.fails.append(("*", "*", f"{report} not found"))
        return out
    labels = chunk_labels(plan)
    with open(report, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["check"] == "overall":
                out.overall_pass = r["status"] == "PASS"
                continue
            if r["status"] == "FAIL":
                out.fails.append((r["filename"], r["chunk"], r["check"]))
                out.failed.add(r["chunk"])
                out.failed.update(
                    labels.get((r["group"], r["filename"], r["chunk"]), [])
                )
            elif r["status"] == "PASS" and r["check"] in FULL_CHECKS:
                out.passed.setdefault(r["chunk"], set()).add(r["check"])
            elif r["status"] == "PASS" and r["check"] == "object":
                out.objects[r["chunk"]] = r["detail"]
    return out


def check(
    plan: Plan,
    sheet: dict[str, list[str]],
    fs,
    verify_report: Path | None = None,
    batch_limit_bytes: int = 5_000 * GB,
    s3_client=None,
) -> list[Row]:
    rows: list[Row] = []
    want = expected_names(plan)

    sheet_names = Counter(n for names in sheet.values() for n in names)
    for name, k in sorted(sheet_names.items()):
        if k > 1:
            rows.append(
                Row(name, "sheet_unique", "FAIL", f"named {k} times in the sheet")
            )
    for name in sorted(set(sheet_names) - set(want)):
        rows.append(
            Row(
                name,
                "plan_matches_sheet",
                "FAIL",
                "in the sheet, not produced by the plan",
            )
        )
    for name in sorted(set(want) - set(sheet_names)):
        rows.append(
            Row(
                name,
                "plan_matches_sheet",
                "FAIL",
                "produced by the plan, not in the sheet",
            )
        )

    fs.invalidate_cache()
    root = plan.dst.rstrip("/")
    own = fs._strip_protocol(root).rstrip("/")
    listing: dict[str, dict] = {}  # name -> the listing's entry for it
    if fs.exists(root):
        for info in fs.ls(root, detail=True):
            if info["name"].rstrip("/") == own and not info.get("size"):
                # the empty marker the S3 console's "Create folder" leaves at
                # the prefix itself: not an object anyone uploads
                continue
            if info.get("type") == "directory":
                listing[s3io.basename(info["name"]) + "/"] = info
                continue
            listing[s3io.basename(info["name"])] = info
    chunk_names = {n for n, v in want.items() if v[0] == "chunk"}
    allowed_extra = {n + ".md5" for n in chunk_names}
    for name in sorted(listing):
        if name in want or name in allowed_extra:
            continue
        detail = "leftover probe" if name == PROBE_NAME else "not named by the sheet"
        rows.append(Row(name, "nothing_extra", "FAIL", detail))

    report = read_verify_report(verify_report, plan)
    passed, failed = report.passed, report.failed
    if verify_report is not None:
        # any FAIL fails the batch, whatever it names: a row the mapping to
        # chunks above misses must still stop the upload
        problems = [f"{f} {c} {k}" for f, c, k in report.fails]
        if not report.overall_pass:
            problems.append("no overall PASS line")
        rows.append(
            Row(
                "(verify report)",
                "verify_clean",
                "FAIL" if problems else "PASS",
                "; ".join(problems[:5])
                + (f" (+{len(problems) - 5} more)" if len(problems) > 5 else ""),
            )
        )
    total = 0
    for name in sorted(set(want) & set(sheet_names)):
        kind, g, fp, c = want[name]
        if name not in listing:
            rows.append(Row(name, "exists", "FAIL", "missing from the destination"))
            continue
        size = listing[name].get("size") or 0
        total += size
        if kind == "file":
            problem, checked = check_copy(fs, plan.dst + name, fp, s3_client)
            rows.append(
                Row(name, "copy", "FAIL" if problem else "PASS", problem or checked)
            )
            continue
        if size >= plan.limit_bytes:
            rows.append(Row(name, "size_below_limit", "FAIL", f"{size / GB:.1f} GB"))
        if name + ".md5" not in listing:
            rows.append(Row(name, "sidecar", "FAIL", "missing"))
        if verify_report is None:
            rows.append(Row(name, "verified_full", "FAIL", "no --verify-report given"))
        elif name in failed:
            rows.append(Row(name, "verified_full", "FAIL", "verify report has a FAIL"))
        elif passed.get(name, set()) != set(FULL_CHECKS):
            rows.append(
                Row(name, "verified_full", "FAIL", "not in a full-level verify report")
            )
        else:
            rows.append(Row(name, "verified_full", "PASS", f"{size} bytes"))
        if verify_report is not None:
            # the report must be about this object, not one it replaced
            now = s3io.object_identity(listing[name])
            seen = report.objects.get(name)
            if seen is None:
                problem = "the verify report records no object for it"
            elif seen != now:
                problem = f"verify read {seen}, the destination holds {now}"
            else:
                problem = None
            rows.append(
                Row(
                    name,
                    "same_as_verified",
                    "FAIL" if problem else "PASS",
                    problem or now,
                )
            )

    status = "FAIL" if total > batch_limit_bytes else "PASS"
    rows.append(
        Row(
            "(batch)",
            "batch_limit",
            status,
            f"{total / GB:.1f} GB of {batch_limit_bytes / GB:.0f} GB",
        )
    )
    return rows


def write_report(path: Path, rows: list[Row]) -> bool:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r.name, r.check, r.status, r.detail])
    return not any(r.status == "FAIL" for r in rows)
