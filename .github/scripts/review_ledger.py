#!/usr/bin/env python3
"""Carry the PR review's findings from one round to the next.

claude-pr-review.yml runs one review per push. Before this script each run
started with no memory of the one before and deleted its comment, so a finding
could be raised, fixed, and raised again in new words; a finding the author
had declined came back every round; and nothing ever checked that a fix had
closed the finding it was written for. The ledger is that memory: a JSON record
of every finding ever raised on the PR with its current status, embedded in the
review comment the way the usage report embeds its figures - inside an HTML
comment, so it renders as nothing and travels with the comment.

Two subcommands, both invoked by the workflow:

    review_ledger.py extract <previous-comment-body> --out previous-findings.json
        Recover the ledger from the last completed round's comment and print a
        one-line JSON summary ({"found": bool, "round": int, "head_sha": str,
        ...}) for the workflow to turn into step outputs. Exits 0 whether or not
        a ledger is there: no ledger means the round starts from scratch, not
        that it fails.

    review_ledger.py embed --review review.md --findings findings.json \\
        --previous previous-findings.json --round N --head-sha SHA \\
        --marker '<!-- claude-pr-review -->' --out review-body.md
        Validate the ledger the reviewer wrote, carry forward every previous
        finding it left out, and assemble the comment: marker, ledger, review,
        footer. A missing or malformed findings.json degrades to re-embedding
        the previous ledger with a visible note, so a round that produced a
        readable review but a broken ledger still posts and loses no state.

Everything the reviewer can get wrong about the ledger is handled here rather
than trusted. The reviewer is a model with a diff to read and eighty turns to
read it in, and "every previous finding appears exactly once with a valid
status" is the kind of invariant to enforce in code, not to ask for nicely.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

LEDGER_PREFIX = "<!-- claude-review-ledger: "
LEDGER_SUFFIX = " -->"
SCHEMA = 1

SEVERITIES = ("blocking", "should-fix", "nit", "question", "pre-existing")
STATUSES = ("open", "resolved", "declined", "withdrawn")
ID_RE = re.compile(r"^F([1-9][0-9]*)$")

DECLINE_HINT = (
    "To decline a finding, reply on the PR with its ID, for example "
    "`F3: by design, <reason>`. A declined finding is not raised again."
)


def warn(message: str) -> None:
    """A GitHub Actions warning annotation.

    Goes to stderr because `extract` owns stdout for its JSON summary. The
    runner scans both streams for workflow commands.
    """
    print(f"::warning::{message}", file=sys.stderr)


def id_number(finding_id: str) -> int:
    return int(ID_RE.match(finding_id).group(1))


def plural(count: int, noun: str) -> str:
    return f"{count} {noun}" + ("" if count == 1 else "s")


def escape_for_comment(blob: str) -> str:
    """Make a JSON string safe inside an HTML comment.

    `--` is not allowed inside an HTML comment and `-->` would end it early.
    `\\u002d` is a hyphen to any JSON parser, so the escaped text parses back
    to the original with no unescape step: extract() just calls json.loads.
    """
    while "--" in blob:
        blob = blob.replace("--", "-\\u002d")
    return blob


def ledger_comment(ledger: dict) -> str:
    blob = json.dumps(ledger, sort_keys=True, separators=(",", ":"))
    return f"{LEDGER_PREFIX}{escape_for_comment(blob)}{LEDGER_SUFFIX}"


def find_ledger(body: str) -> tuple[dict | None, str]:
    """The ledger embedded in a review comment, or (None, why not)."""
    start = body.find(LEDGER_PREFIX)
    if start < 0:
        return None, "the previous review carries no ledger"
    start += len(LEDGER_PREFIX)
    # The JSON cannot contain the suffix: every `--` in it was escaped.
    end = body.find(LEDGER_SUFFIX, start)
    if end < 0:
        return None, "the ledger comment is unterminated"
    try:
        ledger = json.loads(body[start:end])
    except json.JSONDecodeError as exc:
        return None, f"the ledger is not valid JSON: {exc}"
    if not isinstance(ledger, dict) or not isinstance(ledger.get("findings"), list):
        return None, "the ledger has no findings list"
    return ledger, ""


# ---------------------------------------------------------------------------
# Validation. One entry at a time, so one malformed finding costs that finding
# and not the round's whole ledger.
# ---------------------------------------------------------------------------


def _word(value) -> str:
    """Normalise an enum-like field: case and underscores are not disagreements."""
    if not isinstance(value, str):
        return ""
    return value.strip().lower().replace("_", "-")


def _line(value) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 1 else None
    if isinstance(value, str) and value.strip().isdigit():
        return _line(int(value.strip()))
    return None


def _text(value) -> str | None:
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def clean_finding(raw, position: int) -> tuple[dict | None, str]:
    """A finding in canonical shape, or (None, what was wrong with it)."""
    if not isinstance(raw, dict):
        return None, f"entry {position} is not an object"
    finding_id = str(raw.get("id", "")).strip()
    if not ID_RE.match(finding_id):
        return None, f"entry {position} has no valid id (got {raw.get('id')!r})"
    severity = _word(raw.get("severity"))
    if severity not in SEVERITIES:
        return (
            None,
            f"{finding_id}: severity {raw.get('severity')!r} is not one of {SEVERITIES}",
        )
    status = _word(raw.get("status"))
    if status not in STATUSES:
        return (
            None,
            f"{finding_id}: status {raw.get('status')!r} is not one of {STATUSES}",
        )
    path = _text(raw.get("path"))
    if path is None:
        return None, f"{finding_id}: path is missing"
    title = _text(raw.get("title"))
    if title is None:
        return None, f"{finding_id}: title is missing"
    first_round = raw.get("first_round")
    if (
        isinstance(first_round, bool)
        or not isinstance(first_round, int)
        or first_round < 1
    ):
        first_round = None
    return {
        "id": finding_id,
        "severity": severity,
        "status": status,
        "path": path,
        "line": _line(raw.get("line")),
        "title": title,
        "note": _text(raw.get("note")),
        "first_round": first_round,
    }, ""


def clean_all(entries: list) -> tuple[list[dict], list[str]]:
    findings, problems = [], []
    for position, raw in enumerate(entries, start=1):
        finding, problem = clean_finding(raw, position)
        if finding is None:
            problems.append(problem)
        else:
            findings.append(finding)
    return findings, problems


def load_findings(path: str | None) -> list | None:
    """The raw list from a findings file, or None if there is nothing usable.

    Accepts either a bare list or a whole ledger object, since the reviewer
    writes the former and extract() writes the latter.
    """
    if not path:
        return None
    file = Path(path)
    if not file.is_file():
        return None
    try:
        data = json.loads(file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        warn(f"{path} could not be read as JSON: {exc}")
        return None
    if isinstance(data, dict):
        data = data.get("findings")
    if not isinstance(data, list):
        warn(f"{path} holds no findings list")
        return None
    return data


# ---------------------------------------------------------------------------
# Merge. Previous findings keep their order and identity; the reviewer's entry
# for one replaces it; anything the reviewer did not mention is carried forward
# unchanged; new findings must take fresh ids above every old one.
# ---------------------------------------------------------------------------


def merge(
    previous: list[dict], reported: list[dict], round_no: int
) -> tuple[list[dict], list[str], list[str]]:
    """Return (findings, carried_forward_ids, problems)."""
    problems: list[str] = []
    merged = {f["id"]: f for f in previous}
    highest = max((id_number(f["id"]) for f in previous), default=0)
    seen: set[str] = set()
    new: list[dict] = []
    for finding in reported:
        finding_id = finding["id"]
        if finding_id in seen:
            problems.append(
                f"{finding_id} appears twice in findings.json; kept the first"
            )
            continue
        seen.add(finding_id)
        if finding_id in merged:
            previous_entry = merged[finding_id]
            # The round the finding was first raised is history, and the
            # reviewer is not the record of it. A note or line the reviewer
            # left blank is not a retraction of the one already recorded.
            finding["first_round"] = previous_entry["first_round"]
            for key in ("note", "line"):
                if finding[key] is None:
                    finding[key] = previous_entry[key]
            merged[finding_id] = finding
        elif id_number(finding_id) <= highest:
            # An id below the high-water mark that is not in the ledger was
            # either dropped as malformed in an earlier round or invented.
            # Reusing it would make the review text ambiguous against older
            # comments, so it is refused rather than renumbered.
            problems.append(
                f"{finding_id} is new but not above the previous highest id "
                f"F{highest}; dropped"
            )
        else:
            finding["first_round"] = round_no
            new.append(finding)
    carried = [fid for fid in merged if fid not in seen]
    new.sort(key=lambda f: id_number(f["id"]))
    findings = list(merged.values()) + new
    for finding in findings:
        if finding["first_round"] is None:
            finding["first_round"] = round_no
    return findings, carried, problems


# ---------------------------------------------------------------------------
# Rendering.
# ---------------------------------------------------------------------------


def footer_lines(
    round_no: int,
    head_sha: str,
    findings: list[dict],
    carried: list[str],
    notes: list[str],
) -> list[str]:
    open_count = sum(1 for f in findings if f["status"] == "open")
    settled = len(findings) - open_count
    summary = (
        f"Round {round_no}, reviewed at {head_sha[:7]}: "
        f"{plural(open_count, 'finding')} open, {settled} settled."
    )
    lines = [f"<sub>{summary} {DECLINE_HINT}</sub>"]
    if carried:
        lines.append(
            "<sub>Carried forward without a verdict this round: "
            f"{', '.join(carried)}.</sub>"
        )
    lines.extend(f"<sub>{note}</sub>" for note in notes)
    return lines


def assemble(marker: str, ledger: dict, review: str, footer: list[str]) -> str:
    # The marker stays the first line: it is how the next run recognises this
    # comment as a review to supersede. The ledger comes second so extract()
    # finds it before any text that might mention the prefix.
    lines = [marker, ledger_comment(ledger), "", review.rstrip("\n"), ""]
    lines.extend(footer)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Commands.
# ---------------------------------------------------------------------------


def cmd_extract(args: argparse.Namespace) -> int:
    body_file = Path(args.body)
    body = body_file.read_text(encoding="utf-8") if body_file.is_file() else ""
    if body.strip():
        ledger, reason = find_ledger(body)
    else:
        ledger, reason = None, "there is no previous review"

    summary = {
        "found": False,
        "reason": reason or None,
        "round": None,
        "head_sha": None,
        "findings": 0,
        "open": 0,
    }
    if ledger is not None:
        findings, problems = clean_all(ledger["findings"])
        for problem in problems:
            warn(f"previous ledger: {problem}")
        round_no = ledger.get("round")
        head_sha = ledger.get("head_sha")
        cleaned = {
            "schema": SCHEMA,
            "round": round_no if isinstance(round_no, int) else None,
            "head_sha": head_sha if isinstance(head_sha, str) and head_sha else None,
            "findings": findings,
        }
        Path(args.out).write_text(
            json.dumps(cleaned, indent=2) + "\n", encoding="utf-8"
        )
        summary.update(
            found=True,
            reason=None,
            round=cleaned["round"],
            head_sha=cleaned["head_sha"],
            findings=len(findings),
            open=sum(1 for f in findings if f["status"] == "open"),
        )
    print(json.dumps(summary))
    return 0


def cmd_embed(args: argparse.Namespace) -> int:
    review_file = Path(args.review)
    review = review_file.read_text(encoding="utf-8") if review_file.is_file() else ""
    if not review.strip():
        print(
            f"::error::{args.review} is missing or empty; nothing to post",
            file=sys.stderr,
        )
        return 1

    previous_raw = load_findings(args.previous)
    previous, problems = clean_all(previous_raw or [])
    for problem in problems:
        warn(f"previous ledger: {problem}")

    notes: list[str] = []
    reported_raw = load_findings(args.findings)
    if reported_raw is None:
        reported: list[dict] = []
        if previous:
            notes.append(
                "The reviewer did not record a findings ledger this round; "
                "the previous round's ledger was carried forward unchanged."
            )
        else:
            notes.append("The reviewer did not record a findings ledger this round.")
    else:
        reported, problems = clean_all(reported_raw)
        for problem in problems:
            warn(f"findings.json: {problem}")
        if problems:
            count = len(problems)
            notes.append(
                f"{count} malformed ledger {'entry was' if count == 1 else 'entries were'} "
                "dropped this round; see the job log."
            )

    findings, carried, problems = merge(previous, reported, args.round)
    for problem in problems:
        warn(f"findings.json: {problem}")
    if carried and reported_raw is not None:
        warn(f"carried forward without a verdict: {', '.join(carried)}")
    # When the whole ledger was missing, the note above already says so; a
    # second line listing every id would only repeat it.
    shown_carried = carried if reported_raw is not None else []

    ledger = {
        "schema": SCHEMA,
        "round": args.round,
        "head_sha": args.head_sha,
        "findings": findings,
    }
    body = assemble(
        args.marker,
        ledger,
        review,
        footer_lines(args.round, args.head_sha, findings, shown_carried, notes),
    )
    Path(args.out).write_text(body, encoding="utf-8")
    open_count = sum(1 for f in findings if f["status"] == "open")
    print(
        f"round {args.round}: {len(findings)} findings in the ledger, "
        f"{open_count} open, {len(carried)} carried forward, "
        f"{len(findings) - len(previous)} new"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    extract = sub.add_parser(
        "extract", help="recover the ledger from a review comment body"
    )
    extract.add_argument("body", help="file holding the previous review comment's body")
    extract.add_argument(
        "--out", required=True, help="where to write the recovered ledger"
    )

    embed = sub.add_parser("embed", help="build the comment body for this round")
    embed.add_argument("--review", required=True, help="the review markdown")
    embed.add_argument("--findings", required=True, help="the reviewer's findings.json")
    embed.add_argument("--previous", help="the previous round's ledger, if any")
    embed.add_argument("--round", type=int, required=True)
    embed.add_argument("--head-sha", required=True)
    embed.add_argument(
        "--marker", required=True, help="first line of every review comment"
    )
    embed.add_argument("--out", required=True, help="where to write the comment body")

    args = parser.parse_args(argv)
    if args.command == "extract":
        return cmd_extract(args)
    return cmd_embed(args)


if __name__ == "__main__":
    sys.exit(main())
