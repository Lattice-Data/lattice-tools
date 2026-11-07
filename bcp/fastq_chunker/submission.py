"""``check-submission``: is the destination exactly what the SRA sheet names?

The last gate before an SRA upload. Reads the destination itself rather than
trusting the run manifest, and changes nothing. Every one of these must hold:

- the plan is for this sheet: the names the plan produces (a chunk's name, or
  an unsplit file's own name) are exactly the sheet's ``filename`` cells;
- the destination holds exactly those objects, nothing missing and nothing
  extra; the only other objects allowed are the ``.md5`` sidecars of chunks;
- every unsplit file has its planned size, and its DB2 CRC64NVME where S3
  holds a full-object checksum to compare it with (the ``run`` copy check);
- every chunk is below the limit, has its sidecar, and passed ``verify
  --level full`` (read count and md5), with no FAIL anywhere for it;
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


def verify_passes(report: Path | None) -> tuple[dict[str, set[str]], set[str]]:
    """(chunk -> full-level checks that PASSed, chunks with any FAIL)."""
    passed: dict[str, set[str]] = {}
    failed: set[str] = set()
    if report is None:
        return passed, failed
    with open(report, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] == "FAIL":
                failed.add(r["chunk"])
            elif r["status"] == "PASS" and r["check"] in FULL_CHECKS:
                passed.setdefault(r["chunk"], set()).add(r["check"])
    return passed, failed


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
    listing = {}
    if fs.exists(root):
        for info in fs.ls(root, detail=True):
            if info.get("type") == "directory":
                listing[s3io.basename(info["name"]) + "/"] = 0
                continue
            listing[s3io.basename(info["name"])] = info.get("size") or 0
    chunk_names = {n for n, v in want.items() if v[0] == "chunk"}
    allowed_extra = {n + ".md5" for n in chunk_names}
    for name in sorted(listing):
        if name in want or name in allowed_extra:
            continue
        detail = "leftover probe" if name == PROBE_NAME else "not named by the sheet"
        rows.append(Row(name, "nothing_extra", "FAIL", detail))

    passed, failed = verify_passes(verify_report)
    total = 0
    for name in sorted(set(want) & set(sheet_names)):
        kind, g, fp, c = want[name]
        if name not in listing:
            rows.append(Row(name, "exists", "FAIL", "missing from the destination"))
            continue
        size = listing[name]
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
