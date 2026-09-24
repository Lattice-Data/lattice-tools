"""The review ledger: what one round of the PR review hands to the next.

``.github/scripts/review_ledger.py`` is not a package, so it is loaded here by
path. It has to hold two invariants the workflow relies on: a ledger that went
into a comment comes back out unchanged, and every previous finding survives a
round exactly once whatever the reviewer wrote or failed to write.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / ".github" / "scripts" / "review_ledger.py"
)
MARKER = "<!-- claude-pr-review -->"
SHA = "0123456789abcdef0123456789abcdef01234567"


def _load():
    spec = importlib.util.spec_from_file_location("review_ledger", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ledger = _load()


def finding(fid, status="open", severity="should-fix", **extra):
    base = {
        "id": fid,
        "severity": severity,
        "status": status,
        "path": "bcp/foo.py",
        "line": 42,
        "title": f"finding {fid}",
        "note": None,
        "first_round": 1,
    }
    base.update(extra)
    return base


def embed(
    tmp_path,
    review="## Review\n\nfine.\n",
    findings=None,
    previous=None,
    round_no=2,
    fetched_at=None,
):
    """Run ``embed`` with files laid out as the workflow lays them out."""
    (tmp_path / "review.md").write_text(review)
    argv = [
        "embed",
        "--review",
        str(tmp_path / "review.md"),
        "--findings",
        str(tmp_path / "findings.json"),
        "--round",
        str(round_no),
        "--head-sha",
        SHA,
        "--marker",
        MARKER,
        "--out",
        str(tmp_path / "body.md"),
    ]
    if findings is not None:
        (tmp_path / "findings.json").write_text(json.dumps(findings))
    if previous is not None:
        (tmp_path / "previous.json").write_text(json.dumps(previous))
        argv += ["--previous", str(tmp_path / "previous.json")]
    if fetched_at is not None:
        argv += ["--comments-fetched-at", fetched_at]
    code = ledger.main(argv)
    body = (tmp_path / "body.md").read_text() if (tmp_path / "body.md").exists() else ""
    return code, body


def extract(tmp_path, body):
    (tmp_path / "comment.md").write_text(body)
    out = tmp_path / "recovered.json"
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        code = ledger.main(["extract", str(tmp_path / "comment.md"), "--out", str(out)])
    summary = json.loads(stdout.getvalue())
    recovered = json.loads(out.read_text()) if out.exists() else None
    return code, summary, recovered


def _ledger_in(body):
    comment = body.splitlines()[1]
    return json.loads(comment[len(ledger.LEDGER_PREFIX) : -len(ledger.LEDGER_SUFFIX)])


# --- shape of the posted comment -------------------------------------------


def test_marker_is_first_then_ledger_then_review_then_footer(tmp_path):
    code, body = embed(tmp_path, findings=[finding("F1")], round_no=1)
    assert code == 0
    lines = body.splitlines()
    assert lines[0] == MARKER, "the superseding run matches on startswith(marker)"
    assert lines[1].startswith(ledger.LEDGER_PREFIX)
    assert lines[1].endswith(ledger.LEDGER_SUFFIX)
    assert "## Review" in body
    assert body.index("## Review") < body.index("<sub>Round 1, reviewed at 0123456")
    assert "1 finding open, 0 settled." in body
    assert "F3: by design" in body, "the footer tells the author how to decline"
    # A clean round carries no degradation notes.
    assert "malformed" not in body
    assert "could not be recorded" not in body
    assert "Carried forward" not in body


def test_empty_review_is_refused_with_its_own_exit_code(tmp_path):
    """Exit 3 is what tells the workflow not to fall back to a bare post. It is
    not 2, which Python itself uses for a usage error or an unopenable script,
    and which must fall back."""
    code, body = embed(tmp_path, review="  \n\t\n", findings=[])
    assert code == ledger.EXIT_EMPTY_REVIEW == 3
    assert body == ""


# --- round trip ---------------------------------------------------------------


def test_embed_then_extract_returns_the_same_findings(tmp_path):
    reported = [
        finding("F1", title="`--alias-namespace` ignored -- see cli.py"),
        finding(
            "F2", severity="nit", status="declined", note="author: by design", line=None
        ),
    ]
    code, body = embed(
        tmp_path, findings=reported, round_no=1, fetched_at="2026-09-24T20:00:00Z"
    )
    assert code == 0

    code, summary, recovered = extract(tmp_path, body)
    assert code == 0
    assert summary["found"] is True
    assert summary["round"] == 1
    assert summary["head_sha"] == SHA
    assert summary["comments_fetched_at"] == "2026-09-24T20:00:00Z"
    assert summary["open"] == 1
    assert recovered["findings"] == reported
    assert recovered["comments_fetched_at"] == "2026-09-24T20:00:00Z"


def test_ledger_comment_holds_no_double_hyphen(tmp_path):
    """`--` ends or corrupts an HTML comment, and titles quote CLI flags."""
    reported = [finding("F1", title="--flag --> other -- and ---")]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    comment = body.splitlines()[1]
    inner = comment[len(ledger.LEDGER_PREFIX) : -len(ledger.LEDGER_SUFFIX)]
    assert "--" not in inner
    assert json.loads(inner)["findings"][0]["title"] == reported[0]["title"]


def test_ledger_is_read_only_from_the_second_line(tmp_path):
    """A ledger quoted in the review text (say, a review of this script) must
    not pass for the real one when the comment has none, and must lose to the
    real one when it has."""
    decoy = ledger.ledger_comment(
        {"schema": 1, "round": 99, "head_sha": "b" * 40, "findings": [finding("F9")]}
    )
    bare = f"{MARKER}\n\n## Review\n\nLine two normally reads `{decoy}`.\n"
    _, summary, recovered = extract(tmp_path, bare)
    assert summary["found"] is False
    assert recovered is None

    real = ledger.ledger_comment(
        {"schema": 1, "round": 1, "head_sha": "a" * 40, "findings": [finding("F1")]}
    )
    both = f"{MARKER}\n{real}\n\n## Review\n\n{decoy}\n"
    _, summary, recovered = extract(tmp_path, both)
    assert summary["round"] == 1
    assert [f["id"] for f in recovered["findings"]] == ["F1"]


# --- extract without a usable ledger -----------------------------------------


def test_extract_reports_not_found_for_a_pre_ledger_review(tmp_path):
    body = f"{MARKER}\n\n## Review: old style (#1)\n\n**Verdict:** fine.\n"
    code, summary, recovered = extract(tmp_path, body)
    assert code == 0
    assert summary["found"] is False
    assert recovered is None


def test_extract_reports_not_found_for_an_empty_body(tmp_path):
    code, summary, recovered = extract(tmp_path, "")
    assert (code, summary["found"], recovered) == (0, False, None)


def test_extract_tolerates_a_corrupt_ledger(tmp_path):
    body = (
        f"{MARKER}\n{ledger.LEDGER_PREFIX}{{not json{ledger.LEDGER_SUFFIX}\n\nreview\n"
    )
    code, summary, recovered = extract(tmp_path, body)
    assert code == 0
    assert summary["found"] is False
    assert "not valid JSON" in summary["reason"]
    assert recovered is None


def test_ledger_without_a_findings_list_is_not_found(tmp_path):
    body = f'{MARKER}\n{ledger.LEDGER_PREFIX}{{"schema":1}}{ledger.LEDGER_SUFFIX}\n'
    _, summary, recovered = extract(tmp_path, body)
    assert summary["found"] is False
    assert "no findings list" in summary["reason"]
    assert recovered is None


def test_extract_drops_malformed_entries_and_keeps_the_rest(tmp_path):
    bad = {
        "id": "F2",
        "severity": "urgent",
        "status": "open",
        "path": "p",
        "title": "t",
    }
    body = f"{MARKER}\n" + ledger.ledger_comment(
        {
            "schema": 1,
            "round": 1,
            "head_sha": "abc1234",
            "findings": [finding("F1"), bad],
        }
    )
    _, summary, recovered = extract(tmp_path, body)
    assert summary["findings"] == 1
    assert [f["id"] for f in recovered["findings"]] == ["F1"]


def test_extract_drops_a_doctored_head_sha_and_timestamp(tmp_path):
    """Both values are spliced into shell by the workflow, so only a strict
    shape gets through; anything else reads as absent."""
    body = f"{MARKER}\n" + ledger.ledger_comment(
        {
            "schema": 1,
            "round": 2,
            "head_sha": "abc\ndef; rm -rf /",
            "comments_fetched_at": "yesterday",
            "findings": [],
        }
    )
    _, summary, recovered = extract(tmp_path, body)
    assert summary["found"] is True
    assert summary["round"] == 2
    assert summary["head_sha"] is None
    assert summary["comments_fetched_at"] is None
    assert recovered["head_sha"] is None


# --- merging the reviewer's ledger with the previous one ----------------------


def test_unmentioned_open_findings_are_named_and_settled_ones_carried_quietly(
    tmp_path,
):
    previous = [
        finding("F1"),
        finding("F2", status="resolved", note="guarded"),
        finding("F3"),
    ]
    reported = [finding("F1", status="resolved", note="fixed at 44"), finding("F4")]
    _, body = embed(tmp_path, findings=reported, previous=previous, round_no=2)
    got = {f["id"]: f for f in _ledger_in(body)["findings"]}
    assert got["F1"]["status"] == "resolved"
    assert got["F2"] == previous[1], "settled and not mentioned: exactly as it was"
    assert got["F3"] == previous[2], "open and not mentioned: carried, but named"
    assert got["F4"]["first_round"] == 2
    assert "Carried forward without a verdict this round: F3." in body


def test_previous_first_round_wins_over_what_the_reviewer_wrote(tmp_path):
    previous = [finding("F1", first_round=1)]
    reported = [finding("F1", first_round=7)]
    _, body = embed(tmp_path, findings=reported, previous=previous, round_no=3)
    assert _ledger_in(body)["findings"][0]["first_round"] == 1


def test_a_previous_entry_without_first_round_is_backfilled(tmp_path):
    previous = [finding("F1", first_round=None)]
    _, body = embed(tmp_path, findings=[], previous=previous, round_no=3)
    assert _ledger_in(body)["findings"][0]["first_round"] == 3


def test_a_blank_note_or_line_keeps_the_previous_one(tmp_path):
    """Re-listing a declined finding without repeating the author's reason must
    not erase the reason; a line the reviewer did not re-locate stays put."""
    previous = [finding("F1", status="declined", note="author: by design", line=42)]
    reported = [finding("F1", status="declined", note=None, line=None)]
    _, body = embed(tmp_path, findings=reported, previous=previous, round_no=2)
    got = _ledger_in(body)["findings"][0]
    assert (got["note"], got["line"]) == ("author: by design", 42)


def test_new_ids_must_sit_above_every_previous_id(tmp_path):
    """F2 is not in the ledger but is below F4: reusing it would collide with
    older comments that may have shown an F2, so it is refused, not renumbered,
    and the reader is told the review text names an ID the ledger lacks."""
    previous = [finding("F1"), finding("F4")]
    reported = [finding("F1"), finding("F4"), finding("F2"), finding("F5")]
    _, body = embed(tmp_path, findings=reported, previous=previous, round_no=2)
    assert [f["id"] for f in _ledger_in(body)["findings"]] == ["F1", "F4", "F5"]
    assert "1 finding in this review could not be recorded in the ledger" in body


def test_duplicate_ids_keep_the_first(tmp_path):
    reported = [finding("F1", title="first"), finding("F1", title="second")]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    assert [f["title"] for f in _ledger_in(body)["findings"]] == ["first"]
    assert "could not be recorded in the ledger" in body


def test_malformed_entries_are_dropped_with_a_visible_note(tmp_path):
    reported = [finding("F1"), {"id": "F2", "severity": "should-fix"}]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    assert [f["id"] for f in _ledger_in(body)["findings"]] == ["F1"]
    assert "1 malformed ledger entry was dropped this round" in body

    reported = [finding("F1"), {"id": "F2"}, {"id": "F3", "status": "open"}]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    assert "2 malformed ledger entries were dropped this round" in body


def test_missing_path_or_title_keeps_the_entry_with_a_placeholder(tmp_path):
    """Dropping the entry would orphan its ID in the review text and make the
    ID unusable next round; a placeholder loses less."""
    reported = [finding("F1", path=None), finding("F2", title="   ")]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    got = {f["id"]: f for f in _ledger_in(body)["findings"]}
    assert got["F1"]["path"] == ledger.MISSING_PATH
    assert got["F2"]["title"] == ledger.MISSING_TITLE
    assert "malformed" not in body


def test_enum_spelling_slips_are_normalised(tmp_path):
    reported = [finding("F1", severity="Should_Fix", status="OPEN", line="42")]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    got = _ledger_in(body)["findings"][0]
    assert (got["severity"], got["status"], got["line"]) == ("should-fix", "open", 42)


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("should fix", "should-fix"),
        ("[Q]", "question"),
        ("Q", "question"),
        ("PRE_EXISTING", "pre-existing"),
        ("Outside this PR", "pre-existing"),
        ("blocker", "blocking"),
    ],
)
def test_severity_synonyms_a_reviewer_plausibly_writes(raw, expected):
    cleaned, _ = ledger.clean_finding(finding("F1", severity=raw), 1)
    assert cleaned["severity"] == expected


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("fixed", "resolved"),
        ("WONTFIX", "declined"),
        ("won't fix", "declined"),
        ("by design", "declined"),
        ("retracted", "withdrawn"),
        ("reopened", "open"),
    ],
)
def test_status_synonyms_a_reviewer_plausibly_writes(raw, expected):
    cleaned, _ = ledger.clean_finding(finding("F1", status=raw), 1)
    assert cleaned["status"] == expected


def test_footer_counts_pre_existing_bugs_separately(tmp_path):
    reported = [
        finding("F1"),
        finding("F2", severity="pre-existing"),
        finding("F3", status="resolved"),
    ]
    _, body = embed(tmp_path, findings=reported, round_no=1)
    assert (
        "1 finding open, 1 settled, 1 pre-existing bug noted outside this PR." in body
    )


def test_missing_findings_file_reembeds_the_previous_ledger(tmp_path):
    previous = [finding("F1"), finding("F2", status="declined")]
    _, body = embed(tmp_path, findings=None, previous=previous, round_no=2)
    assert _ledger_in(body)["findings"] == previous
    assert "did not record a findings ledger" in body
    assert "Carried forward without a verdict" not in body, "one note, not two"
    assert "1 finding open, 1 settled." in body


def test_missing_findings_file_on_round_one_still_posts(tmp_path):
    _, body = embed(tmp_path, findings=None, round_no=1)
    assert body.startswith(MARKER)
    assert _ledger_in(body)["findings"] == []


def test_a_byte_order_mark_in_findings_json_is_tolerated(tmp_path):
    (tmp_path / "findings.json").write_text(
        json.dumps([finding("F1")]), encoding="utf-8-sig"
    )
    (tmp_path / "review.md").write_text("## Review\n")
    code = ledger.main(
        [
            "embed",
            "--review",
            str(tmp_path / "review.md"),
            "--findings",
            str(tmp_path / "findings.json"),
            "--round",
            "1",
            "--head-sha",
            SHA,
            "--marker",
            MARKER,
            "--out",
            str(tmp_path / "body.md"),
        ]
    )
    body = (tmp_path / "body.md").read_text()
    assert code == 0
    assert [f["id"] for f in _ledger_in(body)["findings"]] == ["F1"]
    assert "did not record" not in body


def test_oversized_review_text_is_truncated_but_ledger_and_footer_survive(tmp_path):
    review = "## Review\n\n" + ("x" * 70_000)
    _, body = embed(tmp_path, review=review, findings=[finding("F1")], round_no=1)
    assert len(body) <= 65_536
    assert body.startswith(MARKER)
    assert _ledger_in(body)["findings"][0]["id"] == "F1"
    assert ledger.TRUNCATION_NOTICE.strip() in body
    assert "<sub>Round 1, reviewed at" in body


def test_a_bare_list_and_a_whole_ledger_object_are_both_accepted_as_previous(tmp_path):
    whole = {"schema": 1, "round": 1, "head_sha": "abc", "findings": [finding("F1")]}
    (tmp_path / "previous.json").write_text(json.dumps(whole))
    assert ledger.load_findings(str(tmp_path / "previous.json")) == [finding("F1")]
    (tmp_path / "previous.json").write_text(json.dumps([finding("F1")]))
    assert ledger.load_findings(str(tmp_path / "previous.json")) == [finding("F1")]


@pytest.mark.parametrize("bad_id", ["F0", "F01", "3", "f3", "F-3", ""])
def test_ids_are_f_plus_a_positive_integer(bad_id):
    cleaned, problems = ledger.clean_finding(finding(bad_id), 1)
    assert cleaned is None
    assert "valid id" in problems[0]
